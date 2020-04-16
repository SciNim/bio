type
  Seq* = ref object of RootObj
    chain*, class*: string

  Dna* = ref object of Seq

proc initDna*(chain: string): Dna =
  Dna(chain: chain, class: "DNA")
