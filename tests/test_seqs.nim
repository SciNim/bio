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
