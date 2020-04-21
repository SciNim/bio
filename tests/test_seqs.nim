import os
import posix_utils
import strutils
import unittest

import bio


suite "Test Seq operation":
  setup:
    var dna = Seq(chain: "ACGTGGGGT", class: "DNA")
    let dnaT: Dna = initDna("ACGTGGGGT")

  test "DNA objects construction":
    check dnaT of Dna
    check dnaT of Seq

    check dnaT.class == "DNA"

  test "DNA complement":
    check dnaT.complement == initDna("TGCACCCCA")

  test "DNA reverse complement":
    check dnaT.reverseComplement == initDna("ACCCCACGT")

  #test "DNA translation":
  #  fail "TDB"
  #
  #test "DNA transcription":
  #  fail "TDB"
  #

suite "Test SeqRecord operation":
  setup:
    var dna = Seq(chain: "ACGTGGGGT", class: "DNA")
    var dnaRecord = SeqRecord(record: dna, name: "Sample")
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
