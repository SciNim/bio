## :Author: |author|
## :Version: |libversion|
##
## If your code works with different formats, you can use this module to access
## different file formats
##
import streams

import fasta, fastq
import sequences
export sequences

type
  FileType* = enum
    ftFasta = "fasta"
    ftFastq = "fastq"

proc dumpTo*(record: SequenceRecord, strm: Stream, kind: FileType) =
  ## Writes the `SequenceRecords<sequences.html#SequenceRecord>`_ to the given
  ## stream, following the format requested in `kind`.
  ##
  ## .. code-block::
  ##
  ##   import streams
  ##   import bio / io
  ##
  ##   let mySeq = newDna("TGCACCCCA")
  ##   let myRec = SequenceRecord(name: "My DNA sequence", record: mySeq)
  ##
  ##   block:
  ##     let strm = newFileStream("path/to/file.fas")
  ##     defer: strm.close
  ##
  ##     myRec.dumpTo(strm, ftFasta)
  ##
  case kind
  of ftFasta:
    fasta.dumpTo(record, strm)
  of ftFastq:
    fastq.dumpTo(record, strm)
  else:
    echo "Not implemented"

proc dumpTo*(records: seq[SequenceRecord], strm: Stream, kind: FileType) =
  ## A shortcut to avoid the explicit cycle to write a `seq` of
  ## `SequenceRecords<sequences.html#SequenceRecord>`_
  ##
  ## .. code-block::
  ##
  ##   import bio / io
  ##
  ##   let mySeqA = newDna("TGCACCCCA")
  ##   let mySeqB = newDna("GTGAGAGTG")
  ##   let myRecA = SequenceRecord(name: "My DNA sequence", record: mySeqA)
  ##   let myRecB = SequenceRecord(name: "My DNA sequence", record: mySeqB)
  ##
  ##   block:
  ##     let strm = newFileStream("myOutput.fasta", fmWrite)
  ##     defer: strm.close
  ##
  ##     @[myRecA, myRecB].dumpTo(strm, ftFasta)
  ##
  for record in records:
    dumpTo(record, strm, kind)

iterator sequences*(strm: Stream, kind: FileType, platform: PlatformName = pnNone):
    SequenceRecord {.inline.} =
  ## Iterate through all the `Sequences<sequences.html#Sequence>`_ in a given
  ## stream, yielding `SequenceRecords<sequences.html#SequenceRecord>`_.
  ##
  ## .. code-block::
  ##
  ##   import bio / io
  ##
  ##   var strm = newFileStream("path/to/file.fas")
  ##
  ##   for sequence in sequences(strm, ftFasta):
  ##     doAssert(sequence of SequenceRecord)
  ##
  ## .. code-block::
  ##
  ##   import bio / io
  ##
  ##   var strm = newFileStream("path/to/file.fasq")
  ##
  ##   for sequence in sequences(strm, ftFastq, pnIllumina):
  ##     doAssert(sequence of SequenceRecord)
  ##     echo sequence.meta.getOrDefault("quality")
  ##
  case kind
  of ftFasta:
    for sr in fasta.sequences(strm):
      yield sr
  of ftFastq:
    for sr in fastq.sequences(strm, platform):
      yield sr
  else:
    echo "Not implemented"

iterator sequences*(fName: string, kind: FileType, platform: PlatformName = pnNone):
    SequenceRecord {.inline.} =
  ## Iterate through all the `Sequences<sequences.html#Sequence>`_ in a given
  ## filename, yielding `SequenceRecords<sequences.html#SequenceRecord>`_.
  ##
  ## .. code-block::
  ##
  ##   import bio / fasta
  ##
  ##   for sequence in sequences("path/to/file.fas"):
  ##     doAssert(sequence of SequenceRecord)
  ##
  let strm = newFileStream(fName)
  case kind
  of ftFastq, ftFasta:
    for sr in sequences(strm, kind, platform):
      yield sr
  else:
    echo "Not implemented"
