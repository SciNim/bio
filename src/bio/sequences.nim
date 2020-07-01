## :Author: |author|
## :Version: |libversion|
import sequtils
import strformat
import strtabs
import strutils
import unicode
import tables

import data

type
  Feature* = ref object of RootObj
    key*: string
    location*: string
    qualifiers*: StringTableRef

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

  SequenceRecord* = ref object of RootObj
    ## An intermediate construct to hold a `Sequence<#Sequence>`_ while naming
    ## it.
    name*: string
    record*: Sequence
    features*: seq[Feature]

  SequenceClassError* = object of ValueError

proc newDna*(chain: string): Sequence =
  ## Initializes a new `Sequence<#Sequence>`_ object, autoadding the class
  ## "DNA".
  runnableExamples:
    let dna: Sequence = newDna("TGCACCCCA")
    doAssert dna.class == scDna
  Sequence(chain: chain, class: scDna)

proc newRna*(chain: string): Sequence =
  ## Initializes a new `Sequence<#Sequence>`_ object, autoadding the class
  ## "RNA".
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")
    doAssert rna.class == scRna
  Sequence(chain: chain, class: scRna)

proc newProtein*(chain: string): Sequence =
  ## Initializes a new `Sequence<#Sequence>`_ object, autoadding the class
  ## "Protein".
  runnableExamples:
    let protein: Sequence = newProtein("TWG")
    doAssert protein.class == scProtein
  Sequence(chain: chain, class: scProtein)

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
      "RNA: ACGUGGGGUACGUGGGGUACGUGGGGUACGUGGGGUACGUGGGGUACGUGGGGUACGUGG…"

  let limit: int = 59
  var e: string
  let reprChain: string = s.chain[.. min(s.chain.len - 1, limit)]
  if s.chain.len > limit:
    e = "…"  # This cannot be a char because is over 1 byte
  &"{s.class}: {reprChain}{e}"

proc `$`*(sr: SequenceRecord): string =
  ## Same as the `$<#$,Sequence>`_ for `Sequence<#Sequence>`_, but adding the
  ## `name` of the `SequenceRecord<#SequenceRecord>`_.
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")
    let sequenceRecord = SequenceRecord(name: "MyRna", record: rna)
    doAssert $sequenceRecord == "MyRna [RNA: ACGUGGGGU]"

  &"{sr.name} [{sr.record}]"

proc `$`*(bi: BackwardsIndex): string = &"^{bi.int}"

proc `?=`*(a, b: Sequence): bool =
  ## Compare two `Sequence<#Sequence>`_. `true` if both `class` and `chain` are
  ## the same.
  runnableExamples:
    doAssert newDna("AAACGGG") ?= Sequence(chain: "AAACGGG", class: scDna)
    doAssert (newDna("AAACGGG") ?= newRna("AAACGGG")) == false
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

    doAssert rna[1 .. 4] == newRna("CGU", scRna)
  Sequence(chain: s.chain[hs], class: s.class)

proc `[]`*(sr: SequenceRecord, i: int|BackwardsIndex): char =
  ## Gets a position of a `SequenceRecord<#SequenceRecord>`_ as if it was a
  ## string.
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")
    let sequenceRecord = SequenceRecord(name: "MyRna", record: rna)

    doAssert sequenceRecord[0] == 'A'
    doAssert sequenceRecord[^1] == 'U'
  sr.record[i]

proc `[]`*(sr: SequenceRecord, s: HSlice): SequenceRecord =
  ## Gets a slice of a `SequenceRecord<#SequenceRecord>`_, returning a new
  ## SequenceRecord.
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")
    let sequenceRecord = SequenceRecord(name: "MyRna", record: rna)

    doAssert sequenceRecord[1 .. ^2].record.chain == "CGUGGGG"
  result = SequenceRecord()
  let newSeq = Sequence(chain: sr.record.chain[s], class: sr.record.class)
  result.record = newSeq
  result.name = &"{sr.name} ({s.a} .. {s.b})"

  return result

proc `[]=`*(s: Sequence, i: int|BackwardsIndex, val: char) =
  ## Replaces a position in the chain of a `Sequence<#Sequence>`_.
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")

    rna[0] = 'U'
    doAssert rna.chain == "UCGUGGGGU"

  s.chain[i] = val

proc `[]=`*(sr: SequenceRecord, i: int|BackwardsIndex, val: char) =
  ## Replaces a position in the chain of a `SequenceRecord<#SequenceRecord>`_.
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")
    let sequenceRecord = SequenceRecord(name: "MyRna", record: rna)

    sequenceRecord[0] = 'U'
    doAssert sequenceRecord.record.chain == "UCGUGGGGU"

  sr.record[i] = val

proc `[]=`*[T, U](s: Sequence, x: HSlice[T, U], b: string) =
  ## Replaces a stretch in the chain of a `Sequence<#Sequence>`_.
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")

    rna[1 .. 3] = "AT"
    doAssert rna.chain == "AATGGGGU"

    rna[^2 .. ^1] = "TT"
    doAssert rna.chain == "AATGGGTT"
  s.chain[x] = b

proc `[]=`*[T, U](s: SequenceRecord, x: HSlice[T, U], b: string) =
  ## Replaces a stretch in the chain of a `Sequence<#Sequence>`_.
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")
    let sequenceRecord = SequenceRecord(name: "MyRna", record: rna)

    sequenceRecord[1 .. 3] = "AT"
    doAssert rna.chain == "AATGGGGU"

    rna[^2 .. ^1] = "TT"
    doAssert rna.chain == "AATGGGTT"
  s.record.chain[x] = b

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

iterator items*(sr: SequenceRecord): char =
  ## Iterates a `SequenceRecord<#SequenceRecord>`_ yielding `chars` from the
  ## chain.
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")
    let srRna = SequenceRecord(name: "MyRna", record: rna)

    var newSeq: string
    for position in srRna:
      newSeq.add position

    doAssert newSeq == rna.chain

  for c in sr.record:
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

iterator pairs*(sr: SequenceRecord): tuple[a: int, b: char] =
  ## Iterates a `SequenceRecord<#SequenceRecord>`_ yielding the position and
  ## the char in a tuple.
  runnableExamples:
    let rna: Sequence = newRna("ACGUGGGGU")
    let srRna = SequenceRecord(name: "MyRna", record: rna)

    var count: seq[int]
    var newSeq: string
    for i, position in srRna:
      count.add i
      newSeq.add position

    doAssert newSeq == "ACGUGGGGU"
    doAssert count == @[0, 1, 2, 3, 4, 5, 6, 7, 8]

  var i: int = 0
  for c in sr.record:
    yield (a: i, b: c)
    inc i

proc len*(s: Sequence): int =
  ## Get the length of a `Sequence<#Sequence>`_ chain.
  runnableExamples:
    doAssert len(newDna("AAACGGG")) == len("AAACGGG")
    doAssert len(newDna("AAACGGG")) == 7
  len(s.chain)

proc len*(sr: SequenceRecord): int =
  ## Get the length of a `SequenceRecord<#SequenceRecord>`_ chain.
  runnableExamples:
    let s = newDna("ACTGGTGGA")
    let sequenceRecord = SequenceRecord(name: "SR1", record: s)
    doAssert len("ACTGGTGGA") == len(sequenceRecord)
  len(sr.record)

proc guess*(s: string): Sequence =
  ## Guesses what is a `Sequence<#Sequence>`_ (DNA, RNA or Protein), inferring
  ## from its first bases.
  ##
  ## 'N' and 'X' are not considered for the guess, even if 'N' is a valid
  ## aminoacid.
  ##
  ## If the class cannot be inferred, a generic Sequence (class
  ## `scSequence<#SequenceClass>`_) is returned.
  runnableExamples:
    doAssert guess("ATCGGCATCG") ?= newDna("ATCGGCATCG")
    doAssert guess("AUCGGCAUCG") ?= newRna("AUCGGCAUCG")
    doAssert guess("FSYWLSCPIK") ?= newProtein("FSYWLSCPIK")
  result = Sequence(chain: s, class: scSequence)
  if s.len < 5: return  # Sequence too short to guess

  var q: string

  for c in s:
    if q.len == 80:
      break
    if c in {'-', 'X', 'N'}:  # Add, after tests, N and X
      continue
    q.add c

  let limit: int = int(float(len(q)) * 0.9)

  if countIt(q, it in dnaLetters) >= limit:
    result.class = scDna
  elif countIt(q, it in rnaLetters) >= limit:
    result.class = scRna
  elif countIt(q, it in proteinLetters) >= limit:
    result.class = scProtein

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

    doAssert newDna("ACGT--GGGGT").translate.chain == "TXGX"

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
