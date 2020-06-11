import os
import posix_utils
import strformat
import strutils
import unittest

import bio/sequences
import bio/fasta


suite "Test SequenceRecord operation":
  setup:
    var dna = Sequence(chain: "ACGTGGGGT", class: scDna)
    var dnaRecord = SequenceRecord(record: dna, name: "Sample")
    var tmpOutput: (string, File) = mkstemp("tests/file_")

  teardown:
    tmpOutput[1].close()
    removeFile(tmpOutput[0])

  test "Write records as a FASTA through filename":
    dnaRecord.dumpTo(tmpOutput[0], "fasta")

    check existsFile(tmpOutput[0])

    let fastaFile = tmpOutput[0].readLines(2)

    check fastaFile == @[">Sample", "ACGTGGGGT"]

  test "Write records as a FASTA through handler":
    dnaRecord.dumpTo(tmpOutput[1], "fasta")

    let fastaFile = tmpOutput[0].readLines(2)

    check fastaFile == @[">Sample", "ACGTGGGGT"]

  test "Automatically wrap the long lines to 60 chars":
    dna.chain = repeat(dna.chain, 8)
    dnaRecord.dumpTo(tmpOutput[1], "fasta")

    let fastaFile = tmpOutput[0].readLines(3)

    check fastaFile.len == 3
    check fastaFile[1].len == 60
    check fastaFile[1] == dna.chain[0 .. 59]

  test "Load records from FASTA through filename":
    let records: seq[SequenceRecord] = load(getAppDir() /
                                            "test_files/regular.fas")

    check records.len == 2
    check records[0].name == "First sequence"
    check records[1].name == "Second sequence"
    check records[0].record.chain.len == 120
    check records[1].record.chain.len == 120
    check records[0].record.class == scDna

  test "Load record from FASTA through iterator + filename":
    var records: seq[SequenceRecord]

    for s in sequences(getAppDir() / "test_files/regular.fas"):
      records.add s

    check records.len == 2
    check records[0].name == "First sequence"
    check records[1].name == "Second sequence"
    check records[0].record.chain.len == 120
    check records[1].record.chain.len == 120
    check records[0].record.class == scDna

  test "Load records from FASTA: special cases":
    test "The empty file: returns empty seq":
      var records: seq[SequenceRecord]

      for s in sequences(getAppDir() / "test_files/empty.fas"):
        records.add s

      check records.len == 0

    test "The lonely sequence: returns a seq with only one sequence":
      var records: seq[SequenceRecord]
      for s in sequences(getAppDir() / "test_files/one_sequence.fas"):
        records.add s

      check records.len == 1

    test "The lonely sequence but with extra blank line at end":
      var records: seq[SequenceRecord]
      for s in sequences(getAppDir() / "test_files/one_sequence_extraline.fas"):
        records.add s

      check records.len == 1

  test "Write a bunch of records as a FASTA through filehandler":
    let multiRecords: seq[SequenceRecord] = @[dnaRecord, dnaRecord]

    multiRecords.dumpTo(tmpOutput[1], "fasta")

    check existsFile(tmpOutput[0])

    var fastaFile: seq[string]
    for line in tmpOutput[0].lines:
      fastaFile.add line

    check fastaFile == @[">Sample", "ACGTGGGGT", ">Sample", "ACGTGGGGT"]
