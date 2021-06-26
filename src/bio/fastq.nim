# Straight from the Wikipedia:
# A FASTQ file normally uses four lines per sequence.
#
#  Line 1 begins with a '@' character and is followed by a sequence identifier
#    and an optional description (like a FASTA title line).
#  Line 2 is the raw sequence letters.
#  Line 3 begins with a '+' character and is *optionally* followed by the same
#    sequence identifier (and any description) again.
#  Line 4 encodes the quality values for the sequence in Line 2, and must
#    contain the same number of symbols as letters in the sequence.
#
# Quality:
#   ASCII character by adding 64 to the Phred value
#   A Phred score of a base is: Qphred = −10*log(e) , e is the probability of
#     the base being wrong.
#
# Illumina
#   @HWUSI-EAS100R:6:73:941:1973#0/1
#    Instrument:Flowcell lane:tile:X:Y#IdxMultiplex/paired
#
#   @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
#    Instrument:RunId:FlowcellId:Flowcell lane:Tile:X:Y Pair:Filtered:Control:IdxSequence
#
# NCBI
#
#   @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=72
#   ID IlluminaTag length
#
## :Author: |author|
## :Version: |libversion|
import streams
import strutils
import tables

import sequences
export sequences


type
  FileType* = enum
    ftFastq = "fastq"

iterator sequences*(strm: Stream, kind: FileType=ftFastq): SequenceRecord {.inline.} =
  ## Iterate through all the `Sequences<sequences.html#Sequence>`_ in a given
  ## stream, yielding `SequenceRecords<sequences.html#SequenceRecord>`_.
  ##
  ## .. code-block::
  ##
  ##   import streams
  ##   import bio/fastq
  ##
  ##   var strm = newFileStream("path/to/file.fasq")
  ##
  ##   for sequence in sequences(strm):
  ##     doAssert(sequence of SequenceRecord)
  ##
  ## Reading Streams allows you to consume compressed files directly, e.g.:
  ##
  ## .. code-block::
  ##
  ##   import zip/gzipfiles # Requires https://github.com/nim-lang/zip
  ##
  ##   import streams
  ##   import bio/fasta
  ##
  ##   var strm = newGzFileStream("path/to/file.fasq.gz")
  ##
  ##   for sequence in sequences(strm):
  ##     doAssert(sequence of SequenceRecord)
  ##
  var name, sequence, quality, line: string
  var firstPlusFound: bool

  while strm.readLine(line):
    if (line.startsWith('@') and sequence.len > 0 and quality.len > 0) or strm.atEnd:
      if strm.atEnd:
        quality.add line
      let meta = {"quality": Meta(kind: mkString, metaString: quality)}.toTable
      yield SequenceRecord(name: name,
                           record: guess(sequence.toUpperAscii),
                           meta: meta)
      name = ""
    if line.startsWith('@'): # Can be Name or Quality
      if name.isEmptyOrWhitespace:
        name = line[1 .. ^1]
        sequence = ""
        quality = ""
        firstPlusFound = false
        continue
      if quality.isEmptyOrWhitespace:
        quality = line.strip()
        # At this point, SequenceRecord is ready, yield at top
        continue
    elif line.startsWith('+'): # Can be marker or Quality
      if firstPlusFound: # is Quality
        quality = line.strip()
      else:
        firstPlusFound = true
      continue
    else: # Either Sequence or Quality lines
      if sequence.isEmptyOrWhitespace:
        sequence = line.strip()
      else:
        quality = line.strip()

iterator sequences*(fName: string, kind: FileType=ftFastq): SequenceRecord {.inline.} =
  ## Iterate through all the `Sequences<sequences.html#Sequence>`_ in a given
  ## fileame, yielding `SequenceRecords<sequences.html#SequenceRecord>`_.
  ##
  ## .. code-block::
  ##
  ##   import bio/fastq
  ##
  ##   for sequence in sequences("path/to/file.fasq"):
  ##     doAssert(sequence of SequenceRecord)
  ##
  let strm = newFileStream(fName)

  for sr in strm.sequences():
    yield sr
