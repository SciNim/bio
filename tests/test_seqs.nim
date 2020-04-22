import os
import posix_utils
import strformat
import strutils
import unittest

import bio


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


suite "Test SequenceRecord operation":
  setup:
    var dna = Sequence(chain: "ACGTGGGGT", class: "DNA")
    var dnaRecord = SequenceRecord(record: dna, name: "Sample")
    var tmpOutput: (string, File) = mkstemp("tests/file_")

  teardown:
    tmpOutput[1].close()
    removeFile(tmpOutput[0])

  test "Write records as a FASTA through filename":
    dnaRecord.write(tmpOutput[0], "fasta")

    check existsFile(tmpOutput[0])

    let fastaFile = tmpOutput[0].readLines(2)

    check fastaFile == @[">Sample", "ACGTGGGGT"]

  test "Write records as a FASTA through handler":
    dnaRecord.write(tmpOutput[1], "fasta")

    let fastaFile = tmpOutput[0].readLines(2)

    check fastaFile == @[">Sample", "ACGTGGGGT"]

  test "Automatically wrap the long lines to 60 chars":
    dna.chain = repeat(dna.chain, 8)
    dnaRecord.write(tmpOutput[1], "fasta")

    let fastaFile = tmpOutput[0].readLines(3)

    check fastaFile.len == 3
    check fastaFile[1].len == 60
    check fastaFile[1] == dna.chain[0 .. 59]
