## :Author: |author|
## :Version: |libversion|
import sequtils
import strformat
import strtabs
import strutils
import unicode
import tables

import data
export data

type
  Feature* = ref object of RootObj
    key*: string
    location*: string
    qualifiers*: StringTableRef

  MetaKind* = enum
    mkInt,
    mkFloat,
    mkSeqInt8,
    mkString,
    mkTableString

  MetaObj* = object
    case kind*: MetaKind
    of mkInt: metaInt*: int
    of mkFloat: metaFloat*: float
    of mkSeqInt8: metaSeqInt8*: seq[int8]
    of mkString: metaString*: string
    of mkTableString: metaTableString*: StringTableObj

  MetaRef* = ref MetaObj

  SequenceClass* = enum
    ## The class of the sequence, being `scSequence` a generic one.
    scSequence = "Sequence",
    scDna = "DNA",
    scRna = "RNA",
    scProtein = "Protein"

  Sequence* = ref object of RootObj
    ## A generic object to store the chain of the `Sequence<#Sequence>`_
    ## (either nucleotides or aminoacids).
    chain*: string
    class*: SequenceClass

  SequenceRecord* = ref SequenceRecordObj

  SequenceRecordObj* = object of Sequence
    ## An intermediate construct to expand `Sequence<#Sequence>`_ with other
    ## data:
    ##
    ## * `name` is an arbitrary identifier given to the `Sequence`.
    ## * `Features<#Feature>`_ are the *FEATURES* fields of a GenBank record.
    ## * `Meta<#MetaObj>`_ are other `metadata`. To get access to them you need
    ##    to import the std module `tables<https://nim-lang.org/docs/tables.html>`_.
    ##
    name*: string
    features*: seq[Feature]
    meta*: Table[string, MetaObj]

  SequenceClassError* = object of ValueError

proc stripSeq*(s: string): string =
  ## Remove `Whitespace<https://nim-lang.org/docs/strutils.html#Whitespace>`_
  ## from string (targeting newlines and the like).
  for c in s:
    if likely(c notin Whitespace):
      result.add c

proc newDna*(chain: string): Sequence =
  ## Initializes a new `Sequence<#Sequence>`_ object, autoadding the class
  ## "DNA".
  runnableExamples:
    let dna: Sequence = newDna("TGCACCCCA")
    doAssert dna.class == scDna
  Sequence(chain: stripSeq(chain), class: scDna)

proc newRna*(chain: string): Sequence =
  ## Initializes a new `Sequence<#Sequence>`_ object, autoadding the class
  ## "RNA".
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")
    doAssert rna.class == scRna
  Sequence(chain: stripSeq(chain), class: scRna)

proc newProtein*(chain: string): Sequence =
  ## Initializes a new `Sequence<#Sequence>`_ object, autoadding the class
  ## "Protein".
  runnableExamples:
    let protein: Sequence = newProtein("TWG")
    doAssert protein.class == scProtein
  Sequence(chain: stripSeq(chain), class: scProtein)

proc truncate(chain: string): string =
  ## Truncate a `Sequence<#Sequence>`_ chain to less than 60 chars if longer
  ##
  let limit: int = 59
  result = chain[0 .. min(chain.len - 1, limit)]
  if chain.len > limit:
    result.add "???"

proc `$`*(s: Sequence): string =
  ## `Sequence<#Sequence>`_ representation is limited to "`Class` + 60 chars"
  ## of sequence.
  ##
  ## If you need the whole sequence, access `Sequence.chain<#Sequence>`_ directly
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")
    doAssert $rna == "RNA: ACGUGGGGU"

  runnableExamples:
    import strutils

    let rna: Sequence = newRna("ACGUGGGGU".repeat(10))
    doAssert $rna ==
      "RNA: ACGUGGGGUACGUGGGGUACGUGGGGUACGUGGGGUACGUGGGGUACGUGGGGUACGUGG???"

  &"{s.class}: {truncate(s.chain)}"

proc `$`*(sr: SequenceRecord): string =
  ## Same as the `$<#$,Sequence>`_ for `Sequence<#Sequence>`_, but adding the
  ## `name` of the `SequenceRecord<#SequenceRecord>`_.
  runnableExamples:
    let sequenceRecord = SequenceRecord(
      name: "MyRna", chain: "ACGUGGGGU", class: scRna)
    doAssert $sequenceRecord == "MyRna [RNA: ACGUGGGGU]"

  &"{sr.name} [{sr.class}: {truncate(sr.chain)}]"

proc `$`*(bi: BackwardsIndex): string = &"^{bi.int}"

proc `==`*(a, b: Sequence): bool =
  ## Compare two `Sequence<#Sequence>`_. `true` if both `class` and `chain` are
  ## the same.
  runnableExamples:
    doAssert newDna("AAACGGG") == Sequence(chain: "AAACGGG", class: scDna)
    doAssert (newDna("AAACGGG") == newRna("AAACGGG")) == false
  (a.chain == b.chain) and (a.class == b.class)

proc `[]`*(s: Sequence, i: int|BackwardsIndex): char =
  ## Gets a position of a `Sequence<#Sequence>`_ as if it was a string.
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")

    doAssert rna[0] == 'A'
    doAssert rna[^1] == 'U'
  s.chain[i]

proc `[]`*(s: Sequence, hs: HSlice): Sequence =
  ## Gets a Slice position of a `Sequence<#Sequence>`_
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")

    doAssert rna[1 .. 3] == newRna("CGU")
  Sequence(chain: s.chain[hs], class: s.class)

proc `[]`*(sr: SequenceRecord, s: HSlice): SequenceRecord =
  ## Gets a slice of a `SequenceRecord<#SequenceRecord>`_, returning a new
  ## SequenceRecord that keeps the same Features and Meta, but changes the name
  ## to show the slice.
  runnableExamples:
    let sequenceRecord = SequenceRecord(
      name: "MyRna", chain: "ACGUGGGGU", class: scRna)

    doAssert sequenceRecord[1 .. ^2].chain == "CGUGGGG"
    doAssert sequenceRecord[1 .. ^2].name == "MyRna (1 .. ^2)"
  SequenceRecord(
    chain: sr.chain[s],
    class: sr.class,
    features: sr.features,
    meta: sr.meta,
    name: &"{sr.name} ({s.a} .. {s.b})")

proc `[]=`*(s: Sequence, i: int|BackwardsIndex, val: char) =
  ## Replaces a position in the chain of a `Sequence<#Sequence>`_.
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")

    rna[0] = 'U'
    doAssert rna.chain == "UCGUGGGGU"

  s.chain[i] = val

proc `[]=`*[T, U](s: Sequence, x: HSlice[T, U], b: string) =
  ## Replaces a stretch in the chain of a `Sequence<#Sequence>`_.
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")

    rna[1 .. 3] = "AT"
    doAssert rna.chain == "AATGGGGU"

    rna[^2 .. ^1] = "TT"
    doAssert rna.chain == "AATGGGTT"
  s.chain[x] = b

iterator items*(s: Sequence): char =
  ## Iterates a `Sequence<#Sequence>`_ yielding `chars` from the chain.
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")

    var newSeq: string
    for position in rna:
      newSeq.add position

    doAssert newSeq == rna.chain

  for c in s.chain:
    yield c

iterator pairs*(s: Sequence): tuple[a: int, b: char] =
  ## Iterates a `Sequence<#Sequence>`_ yielding the position and the char in a
  ## tuple.
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")

    var count: seq[int]
    var newSeq: string
    for i, position in rna:
      count.add i
      newSeq.add position

    doAssert newSeq == "ACGUGGGGU"
    doAssert count == @[0, 1, 2, 3, 4, 5, 6, 7, 8]

  var i: int = 0
  for c in s:
    yield (a: i, b: c)
    inc i

proc zip*(sr1, sr2: Sequence): seq[(char, char)] =
  ## Zips_ two `Sequences<#Sequence>`_ in tuples of bases.
  ##
  ##  .. _Zips: https://nim-lang.org/docs/sequtils.html#zip%2C%2C
  runnableExamples:
    let srRna = SequenceRecord(name: "MyRna", chain: "ACGUGGGGU", class: scRna)
    let srDna = SequenceRecord(name: "MyDna", chain: "ACGTGGGGT", class: scDna)

    var pairs: seq[tuple[a, b: char]]
    for pair in zip(srDna, srRna):
      pairs.add pair

    doAssert pairs[0] == ('A', 'A')
    doAssert pairs[^1] == ('T', 'U')

  zip(sr1.chain, sr2.chain)

proc len*(s: Sequence): int =
  ## Get the length of a `Sequence<#Sequence>`_ chain.
  runnableExamples:
    doAssert len(newDna("AAACGGG")) == len("AAACGGG")
    doAssert len(newDna("AAACGGG")) == 7
  len(s.chain)

proc guess*(s: string): SequenceClass =
  ## Guesses what is a chain of characters (DNA, RNA or Protein), inferring
  ## from its first bases.
  ##
  ## 'N' and 'X' are not considered for the guess, even if 'N' is a valid
  ## aminoacid.
  ##
  ## If the class cannot be inferred, a generic Sequence (class
  ## `scSequence<#SequenceClass>`_) is returned.
  runnableExamples:
    doAssert guess("ATCGGCATCG") == scDna
    doAssert guess("AUCGGCAUCG") == scRna
    doAssert guess("FSYWLSCPIK") == scProtein
  if s.len < 5: return scSequence  # Sequence too short to guess

  var q: string

  for c in s:
    if q.len == 80:
      break
    if c in {'-', 'X', 'N'}:  # Add, after tests, N and X
      continue
    q.add c

  let limit: int = int(float(len(q)) * 0.9)

  if countIt(q, it in dnaLetters) >= limit:
    return scDna
  elif countIt(q, it in rnaLetters) >= limit:
    return scRna
  elif countIt(q, it in proteinLetters) >= limit:
    return scProtein
  return scSequence

proc complement*(s: Sequence): Sequence =
  ## Return a new `Dna or Rna<#Sequence>`_ with the complement sequence.
  runnableExamples:
    doAssert newDna("TGCACCCCA").complement.chain == "ACGTGGGGT"
  result = Sequence()
  result.class = s.class

  case s.class
  of scDna:
    for base in s.chain:
      result.chain.add(dnaAmbiguousComplement[base])
  of scRna:
    for base in s.chain:
      result.chain.add(rnaAmbiguousComplement[base])
  else:
    raise newException(SequenceClassError,
                       &"Operation available only for {scDna} or {scRna}.")

proc reverseComplement*(s: Sequence): Sequence =
  ## Return a new `Dna or Rna<#Sequence>`_ with the reverse complement sequence.
  runnableExamples:
    doAssert newDna("ACGTGGGGT").reverseComplement.chain == "ACCCCACGT"

  result = s.complement()

  case s.class
  of scDna, scRna:
    result.chain = reversed(result.chain)
  else:
    raise newException(SequenceClassError,
                       &"Operation available only for {scDna} or {scRna}.")

proc transcript*(s: Sequence): Sequence =
  ## Return a new `Rna<#Sequence>`_ with the transcribed sequence.
  runnableExamples:
    doAssert newDna("ACGTGGGGT").transcript.chain == "ACGUGGGGU"

  result = Sequence()
  case s.class
  of scDna:
    result.class = scRna
    for base in s.chain:
      if base == 'T':
        result.chain.add('U')
      else:
        result.chain.add(base)
  else:
    raise newException(SequenceClassError,
                       &"Operation available only for {scDna}.")

proc backTranscript*(s: Sequence): Sequence =
  ## Return a new `Dna<#Sequence>`_ with the sequence transcribed from RNA.
  runnableExamples:
    doAssert newRna("ACGUGGGGU").backTranscript.chain == "ACGTGGGGT"

  result = Sequence()
  case s.class
  of scRna:
    result.class = scDna
    for base in s.chain:
      if base == 'U':
        result.chain.add('T')
      else:
        result.chain.add(base)
  else:
    raise newException(SequenceClassError,
                       &"Operation available only for {scRna}.")

proc translate*(s: Sequence): Sequence =
  ## Return a new `Protein<#Sequence>`_ with the translated sequence.
  ##
  ## In-frame deletions are translated as `-`, out-of-frame deletions as `X`,
  ## stop codons as `*`.
  ##
  runnableExamples:
    doAssert newDna("ACGTGGGGT").translate.chain == "TWG"
    doAssert newDna("ACGT--GGGGTA").translate.chain == "TXGV"  # Out of frame
    doAssert newDna("ACG---GGGGT").translate.chain == "T-GX"  # In frame
    doAssert newDna("ACGTAAGGGGT").translate.chain == "T*GX"

  result = Sequence()
  var codon: string
  case s.class
  of scDna:
    result.class = scProtein
    for base in s.chain:
      codon.add base
      if codon.len == 3:
        result.chain.add codonTable.getOrDefault(codon, 'X')
        codon = ""
    if codon.len > 0:
      result.chain.add 'X'
  of scRna:
    result = translate(s.backTranscript)
  else:
    raise newException(SequenceClassError,
                       &"Operation available only for {scDna} or {scRna}.")

proc codon*(s: Sequence, position: int): Sequence =
  ## Return a new `Sequence<#Sequence>`_ which is a triplete that includes
  ## the position.
  ##
  ## `position` parameter doesn't account for gaps, meaning that in the
  ## sequence `---------A---T` the `A` is in position 0, and the `T` is in 1.
  ##
  ## All the operations are zero-indexed: if you target the last base of the
  ## first codon, the position is **2**.
  runnableExamples:
    doAssert newDna("ACGTGGGGT").codon(6).chain == "GGT"
    doAssert newRna("A--CGUGGGGU").codon(1).chain == "ACG"

  runnableExamples:
    doAssert newDna("ACGTGGGGT").codon(0).chain == "ACG"
    doAssert newDna("ACGTGGGGT").codon(1).chain == "ACG"
    doAssert newDna("ACGTGGGGT").codon(2).chain == "ACG"
    doAssert newDna("ACGTGGGGT").codon(3).chain == "TGG"

  result = Sequence()

  case s.class
  of scDna, scRna:
    result.class = s.class
    var idx: int
    var codon = newSeqOfCap[char](3)
    for base in s:
      if base in Letters:
        codon.add base
        inc idx
      if idx > position and codon.len == 3:
        result.chain = join(codon)
        break
      if codon.len == 3:
        codon.delete(0..2)
  else:
    raise newException(SequenceClassError,
                       &"Operation available only for {scDna} or {scRna}.")

func gc_any(s: Sequence, stride: int = 1): float =
  ## Return the GC content for the bases at `stride` positions.
  ##
  if s.class == scProtein:
    raise newException(SequenceClassError,
                       &"Operation available only for {scDna} or {scRna}.")

  let stop: int = stride - 1
  var cgs, length: int

  for i, base in s.pairs:
    if i mod stride == stop:
      if base in {'C', 'c', 'G', 'g'}:
        inc cgs
      inc length
  cgs / length

proc gc*(s: Sequence): float =
  ## Return the GC content for a given sequence, only for Dna and Rna.
  ##
  runnableExamples:
    doAssert newDna("AcGTGGCT").gc == 5 / 8
    doAssert newRna("ACGuGCAU").gc == 4 / 8

  return s.gc_any(1)

proc gc3*(s: Sequence): float =
  ## Return the GC3 content for a given sequence, only for Dna and Rna.
  ##
  ## GC3 is the content of GC for the third codon positions.
  ##
  runnableExamples:
    doAssert newDna("AcGTGGCT").gc3 == 2 / 2
    doAssert newRna("ACGuGCAUA").gc3 == 2 / 3

  s.gc_any(3)
