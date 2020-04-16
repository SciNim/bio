import unittest

import bio


suite "Test Seq creation":
  setup:
    var dna = Seq(chain: "ACGTGGGGT", class: "DNA")
    let dnaT: Dna = initDna("ACGTGGGGT")
  test "DNA objects construction":
    check dnaT of Dna
    check dnaT of Seq

    check dnaT.class == "DNA"
