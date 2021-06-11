import math
import sugar
import strformat
import strtabs
import strutils
import unittest

import bio/sequences


suite "Test Sequence operation":
  setup:
    var dna = Sequence(chain: "ACGTGGGGT", class: scDna)
    let dnaT: Sequence = newDna("ACGTGGGGT")
    let rnaT: Sequence = newRna("ACGUGGGGU")
    let proteinT: Sequence = newProtein("TWG")

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

  test "Remove any non printable character from the sequence":
    ## When you load a DNA string from a file with newlines, they shouldn't be
    ## carried to the final sequence.
    let dnaString: string = "ACTGCG\nGCAT"

    check newDna(dnaString).chain == "ACTGCGGCAT"

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
    check dnaT.complement ?= newDna("TGCACCCCA")

  test "DNA reverse complement":
    check dnaT.reverseComplement ?= newDna("ACCCCACGT")

  test "DNA translation":
    ## Remember: DNA to protein
    check dnaT.translate ?= proteinT

  test "DNA translations with indefinitions":
    # Check the indefinitions
    let oddDna: Sequence = newDna("ACGTGGGGTT")
    check oddDna.translate ?= newProtein("TWGX")

    let delDna: Sequence = newDna("ACGT-GGGT")
    check delDna.translate ?= newProtein("TXG")

  test "DNA transcription":
    ## Remember: DNA to RNA
    check dnaT.transcript ?= rnaT

  test "RNA back transcribe":
    ## Remember: RNA to DNA
    check rnaT.backTranscript ?= dnaT

  test "RNA complement":
    check rnaT.complement ?= newRna("UGCACCCCA")

  test "RNA reverse complement":
    check rnaT.reverseComplement ?= newRna("ACCCCACGU")

  test "RNA translation":
    check rnaT.translate ?= proteinT

  test "Sequence class guessed":
    check guess(dnaT.chain) ?= newDna("ACGTGGGGT")
    check guess(rnaT.chain) ?= newRna("ACGUGGGGU")
    check guess("FSYWLSCPIK") ?= newProtein("FSYWLSCPIK")

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
    check myDna.transcript ?= newRna("ACGUGGGGU")
    check myDna.translate ?= newProtein("TWG")

  test "Operations over sequences guessed as RNA":
    let myRna = guess("ACGUGGGGU")

    check myRna.complement ?= newRna("UGCACCCCA")
    check myRna.backTranscript ?= newDna("ACGTGGGGT")

  test "Sequence len":
    check dnaT.len == 9
    check proteinT.len == 3

  test "Sequence point getting":
    check dna[2] == 'G'
    check dna[^1] == 'T'

  test "Sequence slice getting":
    check dna[1 .. 4].chain == "CGTG"

  test "Sequence nucleotide iteration":
    var newChain: string

    for nucl in dna:
      newChain.add nucl

    check newChain == "ACGTGGGGT"

  test "Sequence nucleotide iteration by pairs":
    var newChain: string
    var count: seq[int]

    for i, position in dna:
      newChain.add position
      count.add i

    check count == @[0, 1, 2, 3, 4, 5, 6, 7, 8]
    doAssert newChain == dna.chain

  test "Sequence position replacement":
    dna[0] = 'T'
    check dna.chain == "TCGTGGGGT"

    dna[^1] = 'A'
    check dna.chain == "TCGTGGGGA"

  test "Sequence position asignment using HSlices":
    dna[2 .. ^3] = "AAAA"
    check dna.chain == "ACAAAAGT"

    dna[2 .. 5] = "T"
    check dna.chain == "ACTGT"

  test "Calculate sequence GC content":
    check dna.gc == 6 / 9
    check dnaT.gc == 6 / 9
    check rnaT.gc == 6 / 9

    # Lets do a slice
    check dna[2 .. ^3].gc == 4 / 5

    # Lets do a zero GC
    check newDna("AAAAAATTTT").gc == 0

    # Lets do a zero length sequence
    check classify(newDna("").gc) == fcNan

    # Lets try a upper/lower case mixed
    check newDna("AcGTGGCT").gc == 5 / 8

    # Lets try a protein
    expect SequenceClassError:
      discard proteinT.gc

  test "Calculate sequence GC content at last position of codon":
    let dnaTruncated = newDna("ACGACCTAACC")
    check dna.gc3 == 2 / 3
    check rnaT.gc3 == 2 / 3
    check dnaTruncated.gc3 == 2 / 3

    # Lets do a slice
    check dna[2 .. ^2].gc3 == 1

    # Lets do a zero GC
    check newDna("AAAAAATTTT").gc3 == 0

    # Lets do a zero length sequence
    check classify(newDna("").gc3) == fcNan

    # Lets try a protein
    expect SequenceClassError:
      discard proteinT.gc3

suite "Test more complex sequence operations":
  setup:
    let dnaShifted: Sequence = newDna("ACGT--GGGGT")
    let rnaShifted: Sequence = newRna("ACGU--GGGGU")
    let dnaGapped: Sequence = newDna("ACG---GGGGT")
    let dnaStops: Sequence = newDna("ACGTAAGGGGT")
    let dnaComplement: Sequence = newDna("ACCCC--ACGT")
    let proteinShifted: Sequence = newProtein("TXGX")
    let proteinGapped: Sequence = newProtein("T-GX")
    let proteinStops: Sequence = newProtein("T*GX")

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
    let sequence: Sequence = newDna("ACTGGTGGA")
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

  test "SequenceRecord nucleotide iteration":
    var newChain: string

    for nucl in sr:
      newChain.add nucl

    check newChain == "ACTGGTGGA"

  test "SequenceRecord nucleotide iteration by pairs":
    var newChain: string
    var count: seq[int]

    for i, nucl in sr:
      count.add i
      newChain.add nucl

    check count == @[0, 1, 2, 3, 4, 5, 6, 7, 8]
    check newChain == "ACTGGTGGA"

  test "SequenceRecord position replacement":
    sr[0] = 'T'
    doAssert sr.record.chain == "TCTGGTGGA"

    sr[^1] = 'T'
    doAssert sr.record.chain == "TCTGGTGGT"

  test "SequenceRecord slice assignment":
    sr[0 .. 4] = "AAAA"
    doAssert sr.record.chain == "AAAATGGA"

  test "SequenceRecord getting a codon":
    # From plain sequences (no gaps)
    let expectedA: Sequence = newDna("ACT")
    let expectedB: Sequence = newDna("GGA")
    check sequence.codon(0) ?= expectedA
    check sequence.codon(1) ?= expectedA
    check sequence.codon(2) ?= expectedA
    check sequence.codon(6) ?= expectedB

    let gapped: Sequence = newDna("A-CTGGTG-GA")
    check gapped.codon(2) ?= expectedA
    check gapped.codon(6) ?= expectedB

  test "Zipping two SequenceRecords":
    let sequence2: Sequence = newDna("GTGCATCAGT")
    let sr2: SequenceRecord = SequenceRecord(name: "S2", record: sequence2)

    let pairs = collect(newSeq):
      for s in zip(sr, sr2):
        s

    for i, _ in sequence.pairs:
      check pairs[i] == (sequence[i], sequence2[i])

suite "Test SequenceRecord Features":
  setup:
    var dna = Sequence(chain: "ACGTGGGGT", class: scDna)

  test "Add Features to SequenceRecord":
    var qualifiers = {"gene": "DNCLI2"}.newStringTable
    let feature = Feature(key: "CDS", location: "7..1485",
                          qualifiers: qualifiers)
    var dnaRec = SequenceRecord(name: "DNA sample",
                                record: dna,
                                features: @[feature])

    check dnaRec.features[0].key == "CDS"

suite "Sequence Catchable Errors":
  setup:
    let proteinT: Sequence = newProtein("TWG")
    let dnaT: Sequence = newDna("ACACTCAGCAT")
    let rnaT: Sequence = newRna("CACAGUGUGCA")

  test "Protein Translate":
    expect SequenceClassError:
      discard proteinT.translate

  test "DNA or Protein Backtranscription":
    expect SequenceClassError:
      discard dnaT.backTranscript
    expect SequenceClassError:
      discard proteinT.backTranscript

  test "RNA or Protein Transcription":
    expect SequenceClassError:
      discard rnaT.transcript
    expect SequenceClassError:
      discard proteinT.transcript

  test "Protein reverseComplement":
    expect SequenceClassError:
      discard proteinT.reverseComplement

  test "Protein complement":
    expect SequenceClassError:
      discard proteinT.complement
