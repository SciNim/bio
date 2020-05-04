## :Author: |author|
## :Version: |libversion|
import sequtils
import strformat
import unicode
import tables

import data


type
  Sequence* = ref object of RootObj
    ## A generic object to store the chain of the Sequence (either nucleotides
    ## or aminoacids).
    chain*, class*: string

  SequenceRecord* = ref object of RootObj
    ## An intermediate construct to hold a Sequence while naming it.
    name*: string
    record*: Sequence

  Dna* = ref object of Sequence
  Rna* = ref object of Sequence
  Protein* = ref object of Sequence

proc initDna*(chain: string): Dna =
  ## Initializes a new Dna object
  runnableExamples:
    let dna: Dna = initDna("TGCACCCCA")
  Dna(chain: chain, class: "DNA")

proc initRna*(chain: string): Rna =
  ## Initializes a new Rna object
  runnableExamples:
    let rna: Rna = initRna("ACGUGGGGU")
  Rna(chain: chain, class: "RNA")

proc initProtein*(chain: string): Protein =
  ## Initializes a new Protein object
  runnableExamples:
    let protein: Protein = initProtein("TWG")
  Protein(chain: chain, class: "Protein")

proc `$`*(s: Sequence): string =
  ## Sequence representation is limited to "`Class` + 60 chars" of sequence.
  ##
  ## If you need the whole sequence, access `Sequence.chain` directly
  runnableExamples:
    let rna: Rna = initRna("ACGUGGGGU")
    doAssert $rna == "RNA: ACGUGGGGU"

  runnableExamples:
    import strutils

    let rna: Rna = initRna("ACGUGGGGU".repeat(10))
    doAssert $rna ==
      "RNA: ACGUGGGGUACGUGGGGUACGUGGGGUACGUGGGGUACGUGGGGUACGUGGGGUACGUGG…"

  let limit: int = 59
  var e: string
  let reprChain: string = s.chain[.. min(s.chain.len - 1, limit)]
  if s.chain.len > limit:
    e = "…"  # This cannot be a char because is over 1 byte
  if s.class == "":
    &"Sequence: {reprChain}{e}"
  else:
    &"{s.class}: {reprChain}{e}"

proc `?=`*(a, b: Sequence): bool =
  ## Compare two sequences. `true` if both `class` and `chain` are the same.
  runnableExamples:
    doAssert initDna("AAACGGG") ?= Sequence(chain: "AAACGGG", class: "DNA")
    doAssert (initDna("AAACGGG") ?= initRna("AAACGGG")) == false
  (a.chain == b.chain) and (a.class == b.class)

proc len*(s: Sequence): int =
  ## Get the length of a Sequence chain
  runnableExamples:
    doAssert len(initDna("AAACGGG")) ?= len("AAACGGG")
  len(s.chain)

proc guess*(s: string): Sequence =
  ## Guesses what is a sequence: DNA, RNA or Protein, inferring from its first
  ## bases
  runnableExamples:
    doAssert guess("ATCGGCATCG") ?= initDna("ATCGGCATCG")
    doAssert guess("AUCGGCAUCG") ?= initRna("AUCGGCAUCG")
    doAssert guess("FSYWLSCPIK") ?= initProtein("FSYWLSCPIK")
  if s.len < 5: return Sequence(chain: s) # Sequence too short to guess

  let q = s[0 .. min(80, s.len - 1)]  # Run at most with 80 positions
  let limit: int = int(float(len(q)) * 0.9)

  if countIt(q, it in dnaLetters) >= limit:
    result = initDna(s)
  elif countIt(q, it in rnaLetters) >= limit:
    result = initRna(s)
  elif countIt(q, it in proteinLetters) >= limit:
    result = initProtein(s)
  else:
    result = Sequence(chain: s)

proc complement*(dna: Dna): Dna =
  ## Return a new `Dna object<#Dna>`_ with the complement sequence.
  runnableExamples:
    doAssert initDna("TGCACCCCA").complement.chain == "ACGTGGGGT"

  result = initDna("")
  for base in dna.chain:
    result.chain.add(dnaAmbiguousComplement[base])

proc reverseComplement*(dna: Dna): Dna =
  ## Return a new `Dna object<#Dna>`_ with the reverse complement sequence.
  runnableExamples:
    doAssert initDna("ACGTGGGGT").reverseComplement.chain == "ACCCCACGT"

  result = dna.complement()
  result.chain = reversed(result.chain)

proc transcript*(dna: Dna): Rna =
  ## Return a new `Rna object<#Rna>`_ with the transcribed sequence.
  runnableExamples:
    doAssert initDna("ACGTGGGGT").transcript.chain == "ACGUGGGGU"
  result = initRna("")
  for base in dna.chain:
    if base == 'T':
      result.chain.add('U')
    else:
      result.chain.add(base)

proc translate*(dna: Dna): Protein =
  ## Return a new `Protein object<#Protein>`_ with the translated sequence.
  ##
  ## In-frame deletions are translated as `-`, out-of-frame deletions as `X`,
  ## stop codons as `*`.
  ##
  runnableExamples:
    doAssert initDna("ACGTGGGGT").translate.chain == "TWG"

    doAssert initDna("ACGT--GGGGT").translate.chain == "TXGX"

    doAssert initDna("ACGTAAGGGGT").translate.chain == "T*GX"

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
  ## Return a new `Rna object<#Rna>`_ with the complement sequence.
  runnableExamples:
    doAssert initRna("ACGUGGGGU").complement.chain == "UGCACCCCA"
  result = initRna("")
  for base in rna.chain:
    result.chain.add(rnaAmbiguousComplement[base])

proc reverseComplement*(rna: Rna): Rna =
  ## Return a new `Rna object<#Rna>`_ with the reverse complement sequence.
  runnableExamples:
    doAssert initRna("ACGUGGGGU").reverseComplement.chain == "ACCCCACGU"
  result = rna.complement()
  result.chain = reversed(result.chain)

proc backTranscribe*(rna: Rna): Dna =
  ## Return a new `Dna object<#Dna>`_ with the sequence transcribed from RNA.
  runnableExamples:
    doAssert initRna("ACGUGGGGU").backTranscribe.chain == "ACGTGGGGT"
  result = initDna("")
  for base in rna.chain:
    if base == 'U':
      result.chain.add('T')
    else:
      result.chain.add(base)

proc translate*(rna: Rna): Protein =
  ## Return a new `Protein object<#Protein>`_ with the translated sequence.
  runnableExamples:
    doAssert initRna("ACGUGGGGU").translate.chain == "TWG"
  translate(rna.backTranscribe)
