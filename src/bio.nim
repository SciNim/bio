import unicode
import tables
import bio/data

type
  Seq* = ref object of RootObj
    chain*, class*: string

  Dna* = ref object of Seq

proc initDna*(chain: string): Dna =
  ## Initializes a new Dna object
  Dna(chain: chain, class: "DNA")

proc `==`*(a, b: Dna): bool =
  a.chain == b.chain

proc complement*(dna: Dna): Dna =
  result = initDna("")
  for base in dna.chain:
    result.chain.add(dnaAmbiguousComplement[base])

proc reverseComplement*(dna: Dna): Dna =
  result = dna.complement()
  result.chain = reversed(result.chain)
