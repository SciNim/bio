import os
import posix_utils
import strformat
import strutils
import unittest

import bio/sequences


suite "Test Sequence operation":
  setup:
    var dna = Sequence(chain: "ACGTGGGGT", class: "DNA")
    let dnaT: Dna = initDna("ACGTGGGGT")
    let rnaT: Rna = initRna("ACGUGGGGU")
    let proteinT: Protein = initProtein("TWG")

  test "Echo for the objects":
    check $dnaT == "DNA: ACGTGGGGT"
    check $rnaT == "RNA: ACGUGGGGU"
    check $proteinT == "Protein: TWG"
    let genericSequence: Sequence = Sequence(chain: "GENERIC")
    check $genericSequence == "Sequence: GENERIC"
    check genericSequence.class == ""

  test "Echo truncates very long sequences":
    let s: Sequence = Sequence(chain: 'A'.repeat(200))
    check $s == &"Sequence: {'A'.repeat(60)}…"

  test "DNA objects construction":
    check dnaT of Dna
    check dnaT of Sequence

    check dnaT.class == "DNA"

  test "RNA objects construction":
    check rnaT of Rna
    check rnaT of Sequence

    check rnaT.class == "RNA"

  test "Protein objects construction":
    check proteinT of Protein
    check proteinT of Sequence

    check proteinT.class == "Protein"

  test "DNA complement":
    check dnaT.complement ?= initDna("TGCACCCCA")

  test "DNA reverse complement":
    check dnaT.reverseComplement ?= initDna("ACCCCACGT")

  test "DNA translation":
    ## Remember: DNA to protein
    check dnaT.translate ?= proteinT

  test "DNA translations with indefinitions":
    # Check the indefinitions
    let oddDna: Dna = initDna("ACGTGGGGTT")
    check oddDna.translate ?= initProtein("TWGX")

    let delDna: Dna = initDna("ACGT-GGGT")
    check delDna.translate ?= initProtein("TXG")

  test "DNA transcription":
    ## Remember: DNA to RNA
    check dnaT.transcript ?= rnaT

  test "RNA back transcribe":
    ## Remember: RNA to DNA
    check rnaT.backTranscribe ?= dnaT

  test "RNA complement":
    check rnaT.complement ?= initRna("UGCACCCCA")

  test "RNA reverse complement":
    check rnaT.reverseComplement ?= initRna("ACCCCACGU")

  test "RNA translation":
    check rnaT.translate ?= proteinT

  test "Sequence class guessed":
    check guess(dnaT.chain) ?= initDna("ACGTGGGGT")
    check guess(rnaT.chain) ?= initRna("ACGUGGGGU")
    check guess("FSYWLSCPIK") ?= initProtein("FSYWLSCPIK")

    check guess(proteinT.chain) ?= Sequence(chain: proteinT.chain, class: "")

  test "Sequence len":
    check dnaT.len == 9
    check proteinT.len == 3

suite "Test more complex sequence operations":
  setup:
    let dnaShifted: Dna = initDna("ACGT--GGGGT")
    let rnaShifted: Rna = initRna("ACGU--GGGGU")
    let dnaGapped: Dna = initDna("ACG---GGGGT")
    let dnaStops: Dna = initDna("ACGTAAGGGGT")
    let dnaComplement: Dna = initDna("ACCCC--ACGT")
    let proteinShifted: Protein = initProtein("TXGX")
    let proteinGapped: Protein = initProtein("T-GX")
    let proteinStops: Protein = initProtein("T*GX")

  test "DNA translation with gaps":
    check dnaShifted.translate ?= proteinShifted

  test "DNA translation with stop codons":
    check dnaStops.translate ?= proteinStops

  test "DNA translation with gaps in sync":
    check dnaGapped.translate ?= proteinGapped

  test "DNA complement with gaps":
    check dnaShifted.reverseComplement ?= dnaComplement

  test "RNA complement with gaps":
    check rnaShifted.reverseComplement ?=
      dnaShifted.transcript.reverseComplement

suite "Test sequenceRecord operations":
  setup:
    let sequence: Dna = initDna("ACTGGTGGA")
    let sr: SequenceRecord = SequenceRecord(name: "SR1",
                                            record: sequence)

  test "SequenceRecord length can be checked":
    check sr.len == sequence.len
    check sr.len == 9

  test "SequenceRecord humanized string":
    check $sr == "SR1 [DNA: ACTGGTGGA]"
