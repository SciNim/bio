import strformat
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
  Rna* = ref object of Sequence
  Protein* = ref object of Sequence

proc initDna*(chain: string): Dna =
  ## Initializes a new Dna object
  runnableExamples:
    var dna: Dna = initDna("TGCACCCCA")
  Dna(chain: chain, class: "DNA")

proc initRna*(chain: string): Rna =
  ## Initializes a new Rna object
  runnableExamples:
    var rna: Rna = initDna("ACGUGGGGU")
  Rna(chain: chain, class: "RNA")

proc initProtein*(chain: string): Protein =
  ## Initializes a new Protein object
  runnableExamples:
    var protein: Protein = initProtein("TWG")
  Protein(chain: chain, class: "Protein")

proc `$`*(s: Sequence): string =
  if s.class == "": &"Sequence: {s.chain}" else: &"{s.class}: {s.chain}"

proc `?=`*(a, b: Sequence): bool =
  (a.chain == b.chain) and (a.class == b.class)

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

proc transcript*(dna: Dna): Rna =
  result = initRna("")
  for base in dna.chain:
    if base == 'T':
      result.chain.add('U')
    else:
      result.chain.add(base)

proc translate*(dna: Dna): Protein =
  result = initProtein("")
  var codon: string
  for base in dna.chain:
    codon.add base
    if codon.len == 3:
      result.chain.add codonTable.getOrDefault(codon, 'X')
      codon = ""
  if codon.len > 0:
    result.chain.add 'X'

proc complement*(rna: Rna): Rna =
  result = initRna("")
  for base in rna.chain:
    result.chain.add(rnaAmbiguousComplement[base])

proc reverseComplement*(rna: Rna): Rna =
  result = rna.complement()
  result.chain = reversed(result.chain)

proc backTranscribe*(rna: Rna): Dna =
  result = initDna("")
  for base in rna.chain:
    if base == 'U':
      result.chain.add('T')
    else:
      result.chain.add(base)

proc translate*(rna: Rna): Protein =
  translate(rna.backTranscribe)
