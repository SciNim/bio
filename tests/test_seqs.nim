import strformat
import strutils
import unittest

import bio/sequences


suite "Test Sequence operation":
  setup:
    var dna = Sequence(chain: "ACGTGGGGT", class: scDna)
    let dnaT: Sequence = initDna("ACGTGGGGT")
    let rnaT: Sequence = initRna("ACGUGGGGU")
    let proteinT: Sequence = initProtein("TWG")

  test "Echo for the objects":
    check $dnaT == "DNA: ACGTGGGGT"
    check $rnaT == "RNA: ACGUGGGGU"
    check $proteinT == "Protein: TWG"
    let genericSequence: Sequence = Sequence(chain: "GENERIC")
    check $genericSequence == "Sequence: GENERIC"
    check genericSequence.class == scSequence

  test "Echo truncates very long sequences":
    let s: Sequence = Sequence(chain: 'A'.repeat(200))
    check $s == &"Sequence: {'A'.repeat(60)}…"

  test "DNA objects construction":
    check dnaT of Sequence
    check dnaT.class == scDna

  test "RNA objects construction":
    check rnaT of Sequence
    check rnaT.class == scRna

  test "Protein objects construction":
    check proteinT of Sequence
    check proteinT.class == scProtein

  test "DNA complement":
    check dnaT.complement ?= initDna("TGCACCCCA")

  test "DNA reverse complement":
    check dnaT.reverseComplement ?= initDna("ACCCCACGT")

  test "DNA translation":
    ## Remember: DNA to protein
    check dnaT.translate ?= proteinT

  test "DNA translations with indefinitions":
    # Check the indefinitions
    let oddDna: Sequence = initDna("ACGTGGGGTT")
    check oddDna.translate ?= initProtein("TWGX")

    let delDna: Sequence = initDna("ACGT-GGGT")
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

    check guess(proteinT.chain) ?= Sequence(chain: proteinT.chain, class: scSequence)

  test "A slightly more complex guessing (dashes)":
    check guess("---------------------ACGTGGGGT").class == scDna
    check guess("--------------------FSYWLSCPIK").class == scProtein

    # Indefinitions are resolved
    check guess("NNNNNNNNNNNNNNNNNNNNNNNNNNNNACTGGTGG").class == scDna
    check guess("XXXXXXXXXXXXXXXXXXXXXXXXXXFSYWLSCPIK").class == scProtein

  test "Operations over sequences guessed as DNA":
    let myDna = guess("ACGTGGGGT")

    check myDna.complement ?= dnaT.complement
    check myDna.transcript ?= initRna("ACGUGGGGU")
    check myDna.translate ?= initProtein("TWG")

  test "Operations over sequences guessed as RNA":
    let myRna = guess("ACGUGGGGU")

    echo myRna.complement
    check myRna.complement ?= initRna("UGCACCCCA")
    check myRna.backTranscribe ?= initDna("ACGTGGGGT")

  test "Sequence len":
    check dnaT.len == 9
    check proteinT.len == 3

suite "Test more complex sequence operations":
  setup:
    let dnaShifted: Sequence = initDna("ACGT--GGGGT")
    let rnaShifted: Sequence = initRna("ACGU--GGGGU")
    let dnaGapped: Sequence = initDna("ACG---GGGGT")
    let dnaStops: Sequence = initDna("ACGTAAGGGGT")
    let dnaComplement: Sequence = initDna("ACCCC--ACGT")
    let proteinShifted: Sequence = initProtein("TXGX")
    let proteinGapped: Sequence = initProtein("T-GX")
    let proteinStops: Sequence = initProtein("T*GX")

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
    let sequence: Sequence = initDna("ACTGGTGGA")
    let sr: SequenceRecord = SequenceRecord(name: "SR1",
                                            record: sequence)

  test "SequenceRecord length can be checked":
    check sr.len == sequence.len
    check sr.len == 9

  test "SequenceRecord humanized string":
    check $sr == "SR1 [DNA: ACTGGTGGA]"

  test "SequenceRecord point getting":
    check sr[0] == 'A'
    check sr[^1] == 'A'

  test "SequenceRecord slicing":
    let slice = sr[1..^2]
    check slice.record.chain == "CTGGTGG"
    check slice.record.class == sequence.class
    check slice.name == sr.name & " (1 .. ^2)"
