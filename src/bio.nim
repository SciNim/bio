import unicode
import tables
import bio/data

type
  Seq* = ref object of RootObj
    chain*, class*: string

  SeqRecord* = ref object of RootObj
    name*: string
    record*: Seq

  Dna* = ref object of Seq

proc initDna*(chain: string): Dna =
  ## Initializes a new Dna object
  Dna(chain: chain, class: "DNA")

proc `==`*(a, b: Dna): bool =
  a.chain == b.chain

proc write*(record: SeqRecord, fHandler: File, kind: string) =
  #kind is a string as in "fasta", to support different formats
  const wrapSize: int = 60
  fHandler.write(">", record.name, "\n")
  var currentLen: int
  for i in 0 .. (record.record.chain.len div wrapSize):
    currentLen = min(wrapSize, record.record.chain[(i * wrapSize) .. ^1].len) - 1
    fHandler.write(
      record.record.chain[(i * wrapSize) .. ((i * wrapSize) + currentLen)])
    fHandler.write("\n")
  fHandler.flushFile

proc write*(record: SeqRecord, fName, kind: string) =
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
