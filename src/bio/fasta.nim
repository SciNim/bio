## :Author: |author|
## :Version: |libversion|
import strutils

import sequences
export sequences


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
      name = line[1..^1]
    else:
      sequence.add line

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
