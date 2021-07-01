# Straight from the Wikipedia: https://en.wikipedia.org/wiki/FASTQ_format
#
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
# The logic has to change:
#  - Find the first @
#  - Read Sequence until +
#  - Read Quality until len(quality) == len(sequence)
#  - Repeat (check that a new @ exists)
#
## :Author: |author|
## :Version: |libversion|
import sets
import streams
import strscans
import strutils
import tables

import io
import sequences
export sequences


type
  TagName* = enum
    tnNone,
    tnIllumina,
    tnIlluminaOld,
    tnNcbiSra

template scan(charSet: typed) =
  while result + start < input.len and input[result + start] in charSet:
    value.add input[result + start]
    result.inc

proc yn(input: string, value: var string, start: int): int =
  scan({'Y', 'N'})

proc bases(input: string, value: var string, start: int): int =
  scan(Digits + {'A' .. 'Z'})

proc alnum(input: string, value: var string, start: int): int =
  # As defined in http://maq.sourceforge.net/fastq.shtml
  scan(Letters + Digits + {'_', '.', '-'})

proc parseTag*(tag: string, tagName: TagName = tnNone): Table[string, MetaObj] =
  ## Parses the sequence identifier into its parts
  ##
  ## TODO example here
  var instrumentId, runId, flowCellId, index, filtered, empty: string
  var flowCellLane, tileNumber, xTile, yTile, pairing, readNumber: int
  case tagName
  of tnIllumina:
    if scanf(tag, "${alnum}:${alnum}:${alnum}:$i:$i:$i:$i $i:${yn}:$i:${alnum}",
             instrumentId, runId, flowCellId, flowCellLane, tileNumber, xTile,
             yTile, readNumber, filtered, pairing, index):
      result["instrumentId"] = MetaObj(kind: mkString, metaString: instrumentId)
      result["runId"] = MetaObj(kind: mkString, metaString: runId)
      result["flowCellId"] = MetaObj(kind: mkString, metaString: flowCellId)
      result["flowCellLane"] = MetaObj(kind: mkInt, metaInt: flowCellLane)
      result["tileNumber"] = MetaObj(kind: mkInt, metaInt: tileNumber)
      result["xTile"] = MetaObj(kind: mkInt, metaInt: xTile)
      result["yTile"] = MetaObj(kind: mkInt, metaInt: yTile)

      result["readNumber"] = MetaObj(kind: mkInt, metaInt: readNumber)
      result["isFiltered"] = MetaObj(kind: mkString, metaString: filtered)
      result["control"] = MetaObj(kind: mkInt, metaInt: pairing)
      result["barcode"] = MetaObj(kind: mkString, metaString: index)
  of tnIlluminaOld:
    if scanf(tag, "${alnum}:$i:$i:$i:$i#${bases}/$i", instrumentId, flowCellLane,
             tileNumber, xTile, yTile, index, pairing):
      result["instrumentId"] = MetaObj(kind: mkString, metaString: instrumentId)
      result["flowCellLane"] = MetaObj(kind: mkInt, metaInt: flowCellLane)
      result["tileNumber"] = MetaObj(kind: mkInt, metaInt: tileNumber)
      result["xTile"] = MetaObj(kind: mkInt, metaInt: xTile)
      result["yTile"] = MetaObj(kind: mkInt, metaInt: yTile)
      result["index"] = MetaObj(kind: mkString, metaString: index)
      result["pairing"] = MetaObj(kind: mkInt, metaInt: pairing)
  of tnNcbiSra:
    if scanf(tag, "${alnum} ${alnum}:$i:$i:$i:$i", empty, instrumentId, flowCellLane,
             tileNumber, xTile, yTile, index, pairing):
      result["instrumentId"] = MetaObj(kind: mkString, metaString: instrumentId)
      result["flowCellLane"] = MetaObj(kind: mkInt, metaInt: flowCellLane)
      result["tileNumber"] = MetaObj(kind: mkInt, metaInt: tileNumber)
      result["xTile"] = MetaObj(kind: mkInt, metaInt: xTile)
      result["yTile"] = MetaObj(kind: mkInt, metaInt: yTile)
  else:
    return

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
  var loadQuality: bool

  while strm.readLine(line):
    if (line.startsWith('@') and
        len(sequence) > 0 and len(quality) == len(sequence)) or strm.atEnd:
      if strm.atEnd:
        quality.add line
      let meta = {"quality": MetaObj(kind: mkString, metaString: quality)}.toTable
      yield SequenceRecord(name: name,
                           record: guess(sequence.toUpperAscii),
                           meta: meta)
      loadQuality = false
      sequence = ""

    if loadQuality:  # '+' was detected, we are in quality lines
      if len(quality) < len(sequence):
        quality.add line
    else:  # Is either name, sequence or the '+' separator line
      if line.startsWith('@'):  # Is the name line
        name = line[1 .. ^1]
      elif line.startsWith('+'):  # The '+' separator
        loadQuality = true
        quality = newStringOfCap(len(sequence))
      else:  # Is sequence
        sequence.add line

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
