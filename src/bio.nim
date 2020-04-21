import unicode
import tables
import bio/data

type
  Sequence* = ref object of RootObj
    chain*, class*: string

  SequenceRecord* = ref object of RootObj
    name*: string
    record*: Sequence

  Dna* = ref object of Sequence

proc initDna*(chain: string): Dna =
  ## Initializes a new Dna object
  Dna(chain: chain, class: "DNA")

proc `==`*(a, b: Dna): bool =
  a.chain == b.chain

proc write*(record: SequenceRecord, fHandler: File, kind: string) =
  # TODO kind is a string as in "fasta", to support different formats
  const wrapSize: int = 60
  fHandler.write(">", record.name)

  for i, base in record.record.chain.pairs:
    if i mod wrapSize == 0:
      fHandler.write("\n")
    fHandler.write(base)
  fHandler.write("\n")
  fHandler.flushFile

proc write*(record: SequenceRecord, fName, kind: string) =
  let fHandler: File = open(fName, fmWrite)
  defer: fHandler.close()

  write(record, fHandler, kind)

proc complement*(dna: Dna): Dna =
  result = initDna("")
  for base in dna.chain:
    result.chain.add(dnaAmbiguousComplement[base])

proc reverseComplement*(dna: Dna): Dna =
  result = dna.complement()
  result.chain = reversed(result.chain)
