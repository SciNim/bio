## :Author: |author|
## :Version: |libversion|
import strformat
import strutils
import tables

import sequences
export sequences


type
 Index* = ref object of RootObj
    ## The index is a `table` with the name of the sequence as a string, except
    ## the initial '>', and its position in the file. The filename which the
    ## Index refers to is stored in `source`.
    source*: string
    table*: TableRef[string, int64]

iterator sequences*(fName: string, kind: string="fasta"):
  SequenceRecord {.inline.} =
  ## Iterate through all the `Sequences<sequences.html#Sequence>`_ in a given
  ## filename, yielding `SequenceRecords`.
  ##
  ## .. code-block::
  ##
  ##   import bio/fasta
  ##
  ##   for sequence in sequences("path/to/file.fas"):
  ##     doAssert(sequence of SequenceRecord)

  let fileIn: File = open(fName)
  defer: fileIn.close

  var name, sequence: string
  for line in fileIn.lines:
    if (line.startsWith('>') and sequence.len > 0) or fileIn.endOfFile:
      if fileIn.endOfFile:
        sequence.add line
      yield SequenceRecord(name: name, record: guess(sequence.toUpperAscii))

    if line.startsWith('>'):
      sequence = ""
      name = line[1 .. ^1]
    else:
      sequence.add line

proc newIndex*(fName: string): Index =
  ## Build an `Index<#Index>`_ for a given fasta file.
  ##
  ## You should build an index if you plan to access a few sequences by name
  ## randomly and repeatedly. If you'll need all the sequences and specially
  ## if you only need it once, use `sequences
  ## iterator<#sequences.i,string,string>`_
  ##
  ## This index doesn't attempt to compete with powerful indexing like
  ## `faidx<https://www.htslib.org/doc/samtools-faidx.html>`_ . E.g. for FASTAs
  ## like the human genome sequence, with a few huge sequences, it's almost
  ## useless.
  ##
  ## .. code-block::
  ##
  ##   import bio/fasta
  ##
  ##   let index: Index = newIndex("path/to/file.fas")

  var table = newTable[string, int64]()

  let fileIn: File = open(fName)
  defer: fileIn.close

  const newLine = 1

  for line in fileIn.lines:
    if line.startsWith('>'):
      table[line[1 .. ^1]] = fileIn.getFilePos - len(line) - newLine

  Index(source: fName, table: table)

proc `[]`*(index: Index, fileIn: File, name: string): SequenceRecord {.inline.} =
  ## Returns a `SequenceRecords<sequences.html#SequenceRecord>`_ from an
  ## `Index<#Index>`_ using an already opened file handler.
  ##
  ## Ignores the file path provided with `Index.source` field, and doesn't
  ## check the validity of `fileIn` handler (e.g. if it's in `fmRead` mode).
  ##
  ## .. code-block::
  ##
  ##   import bio/fasta
  ##
  ##   let index: Index = newIndex("path/to/file.fas")
  ##
  ##   block:
  ##     let fastaFile = open("path/to/file.fas")
  ##     defer: fastaFile.close
  ##     echo index[fastaFile, "My Record"]
  ##     echo index[fastaFile, "Other Record"]
  fileIn.setFilePos(index.table[name])
  let name = fileIn.readLine[1 .. ^1]
  var sequence: string
  for line in fileIn.lines:
    if line.startsWith('>') or fileIn.endOfFile:
      if fileIn.endOfFile:
        sequence.add line
      return SequenceRecord(name: name,
                            record: guess(sequence.toUpperAscii))
    sequence.add line


proc `[]`*(index: Index, name: string): SequenceRecord {.inline.} =
  ## Returns a `SequenceRecords<sequences.html#SequenceRecord>`_ from an
  ## `Index<#Index>`_.
  ##
  ## Doesn't check if the file provided at `Index.source` exists.
  ##
  ## .. code-block::
  ##
  ##   import bio/fasta
  ##
  ##   let index: Index = newIndex("path/to/file.fas")
  ##   echo index["My Record"]

  let fileIn: File = open(index.source)
  defer: fileIn.close

  return index[fileIn, name]

proc `$`*(index: Index): string {.inline.} =
  ## Returns a brief description of the `Index<#Index>`_.

  &"Index for {index.source}, lenght: {len(index.table)}"

proc load*(fName: string, kind: string="fasta"): seq[SequenceRecord] =
  ## Load a `seq` of `SequenceRecords<sequences.html#SequenceRecord>`_ from a
  ## filename.
  ##
  ## *Careful*: this loads **all** the file at once. If you need an interator
  ## over the file, use the `sequences iterator<#sequences.i,string,string>`_
  ##
  ## .. code-block::
  ##
  ##    import bio/fasta
  ##
  ##    let mySeqs = load("path/to/file.fas")

  for sequence in sequences(fName):
    result.add sequence

proc dumpTo*(record: SequenceRecord, fHandler: File, kind: string="fasta") =
  ## Write a `SequenceRecords<sequences.html#SequenceRecord>`_ to `fHandler`,
  ## wrapping the sequence by 60 positions.
  ## The name of the SequenceRecord remains untouched.
  ##
  ## TBD: `kind` is a string as in "fasta", to support different formats.
  ## Right now only FASTA files are supported.
  ##
  ## To write a FASTA file `myOutput.fasta` with the contents:
  ##
  ## .. code-block::
  ##   `import bio/fasta`
  ##
  ##   let mySeq = newDna("TGCACCCCA")
  ##   let myRec = SequenceRecord(name: "My DNA sequence", record: mySeq)
  ##
  ##   block:
  ##     let fastaOut = open("myOutput.fasta", fmWrite)
  ##     defer: fastaOut.close
  ##     myRec.dumpTo(fastaOut)
  ##
  const wrapSize: int = 60
  fHandler.write(">", record.name)

  for i, base in record.record.chain.pairs:
    if i mod wrapSize == 0:
      fHandler.write('\n')
    fHandler.write(base)
  fHandler.write('\n')
  fHandler.flushFile

proc dumpTo*(records: seq[SequenceRecord], fHandler: File, kind: string="fasta") =
  ## A shortcut to avoid the explicit cycle to write a `seq` of
  ## `SequenceRecords<sequences.html#SequenceRecord>`_
  ##
  ## .. code-block::
  ##   import bio/fasta
  ##
  ##   let mySeqA = newDna("TGCACCCCA")
  ##   let mySeqB = newDna("GTGAGAGTG")
  ##   let myRecA = SequenceRecord(name: "My DNA sequence", record: mySeqA)
  ##   let myRecB = SequenceRecord(name: "My DNA sequence", record: mySeqB)
  ##
  ##   block:
  ##     let fastaOut = open("myOutput.fasta", fmWrite)
  ##     defer: fastaOut.close
  ##     @[myRecA, myRecB].dumpTo(fastaOut)
  ##
  for sr in records:
    sr.dumpTo(fHandler, kind)

proc dumpTo*(record: SequenceRecord, fName: string, kind: string="fasta") =
  ## Same as `write-through-handler proc<#dumpTo,SequenceRecord,File,string>`_
  ## but you only need to point out the name of the file.
  ##
  ## If the file exists, it will be silently overwritten.
  ##
  ## .. code-block::
  ##   import bio/fasta
  ##
  ##   let mySeq = newDna("TGCACCCCA")
  ##   let myRec = SequenceRecord(name: "My DNA sequence", record: mySeq)
  ##
  ##   myRec.dumpTo("myOutput.fasta")
  ##
  let fHandler: File = open(fName, fmWrite)
  defer: fHandler.close()

  record.dumpTo(fHandler, kind)

proc dumpTo*(records: seq[SequenceRecord], fName: string, kind: string="fasta") =
  ## Same as `write-through-handler proc<#dumpTo,seq[SequenceRecord],File,string>`_
  ## but you only need to point out the name of the file.
  ##
  ## If the file exists, it will be silently overwritten.
  ##
  ## .. code-block::
  ##   import bio/fasta
  ##
  ##   let mySeqA = newDna("TGCACCCCA")
  ##   let mySeqB = newDna("GTGAGAGTG")
  ##   let myRecA = SequenceRecord(name: "My DNA sequence", record: mySeqA)
  ##   let myRecB = SequenceRecord(name: "My DNA sequence", record: mySeqB)
  ##
  ##   myRecs = @[myRecA, myRecB]
  ##   myRecs.dumpTo("myOutput.fasta")
  ##
  let fHandler: File = open(fName, fmWrite)
  defer: fHandler.close()

  for record in records:
    record.dumpTo(fHandler, kind)
