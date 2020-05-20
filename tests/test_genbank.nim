import os
import strutils
import strtabs
import unittest

import bio/genbank


suite "Test SequenceRecord operation":
  setup:
    discard
    #var dna = Sequence(chain: "ACGTGGGGT", class: scDna)
    #var dnaRecord = SequenceRecord(record: dna, name: "Sample")

  test "Load a Genbank file as a Sequence Record":
    var gb: seq[SequenceRecord] = load(getAppDir() / "test_files/noref.gb")

    check gb[0].name == "Homo sapiens dynein, cytoplasmic, " &
                        "light intermediate polypeptide 2 (DNCLI2), mRNA."
    check len(gb[0]) == 1622
    check gb[0].record.chain.startsWith("GGCAAGATGGCG")
    check gb[0].record.chain.endsWith("AAAAAAAAAAAA")

  test "Load a Genbank file and get the Features included":
    var gb: seq[SequenceRecord] = load(getAppDir() / "test_files/noref.gb")

    check len(gb[0].features) == 3
    check gb[0].features[0].key == "source"
    check gb[0].features[1].key == "gene"
    check gb[0].features[2].key == "CDS"

    check len(gb[0].features[0].qualifiers) == 3
    check len(gb[0].features[1].qualifiers) == 3
    check len(gb[0].features[2].qualifiers) == 7

    check gb[0].features[2].qualifiers["db_xref"] == "LocusID:1783,GI:5453634"
