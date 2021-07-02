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
#   Paper on qualities: https://doi.org/10.1093/nar/gkp1137
#
# Illumina
#   @HWUSI-EAS100R:6:73:941:1973#0/1
#    Instrument:FlowcellLane:tile:X:Y#IdxMultiplex/paired
#
#   @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
#    Instrument:RunId:FlowcellId:FlowcellLane:Tile:X:Y Pair:Filtered:Control:IdxSequence
#
# NCBI
#
#   @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=72
#   ID IlluminaTag length
#
## :Author: |author|
## :Version: |libversion|
import math
import sets
import streams
import strscans
import strutils
import sugar
import tables

import io
import sequences
export sequences


type
  PlatformName* = enum
    pnNone,
    pnIllumina,
    pnIlluminaOld,
    pnNcbiSra

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

func getPes*(scores: seq[int8], platform: PlatformName = pnNone): seq[float] =
  ## Calculates the Probability of Error (Pe) of each base from the letter seq
  ##
  ## Mainly for internal use, can be called with a sequence of Phred scores:
  ##
  runnableExamples:
    import std / [math, sugar]
    let scores = @[33'i8, 40, 50, 60, 80, 126]

    let errors = collect(newSeqOfCap(scores.len)):
      for s in getPes(scores):
        round(s, 5)

    doAssert errors == @[1.0, 0.19953, 0.01995, 0.002, 2e-05, 0.0]

  runnableExamples:
    # In 2004 Solexa introduced their own scoring system. Illumina older than
    # 1.3 uses this format, specified with the tag 'pnIlluminaOld'.
    #
    # More about this wonderful event at https://doi.org/10.1093/nar/gkp1137
    #
    import std / [math, sugar]
    let scores = @[59'i8, 64, 70, 80, 90, 126]

    let errors = collect(newSeqOfCap(scores.len)):
      for s in getPes(scores, pnIlluminaOld):
        round(s, 5)

    doAssert errors == @[0.75975, 0.5, 0.20076, 0.0245, 0.00251, 0.0]

  runnableExamples:
    # From Illumina 1.3 or newer, the minimum value for a Qscore is 64.
    #
    import std / [math, sugar]
    let scores = @[64'i8, 80, 90, 126]

    let errors = collect(newSeqOfCap(scores.len)):
      for s in getPes(scores, pnIllumina):
        round(s, 5)

    doAssert errors == @[1.0, 0.02512, 0.00251, 0.0]

  let offset = case platform
    of pnIlluminaOld, pnIllumina: 64
    else: 33

  collect(newSeqOfCap(len(scores))):
    for i in scores:
      case platform
      of pnIlluminaOld:
        1 / (pow(10, (i - offset).float / 10.0) + 1)
      else:
        pow(10, -(i - offset).float / 10.0)

func parseQuality*(quality: string): seq[int8] =
  ## Parses a quality string into a seq of signed integers of size 8.
  ##
  ## Doesn't take into account the differences between platforms. Refer to <TBD>
  ##
  collect(newSeqOfCap(len(quality))):
    for c in quality:
      int8(ord(c))

func parseTag*(tag: string, platformName: PlatformName = pnNone): Table[string, MetaObj] =
  ## Parses the tag identifier into its parts.
  ##
  ## There are some docs explaining how this identifier work, and none of them
  ## seems to agree. This procedure attempts the following:
  ##
  ## - For Illumina > 1.4
  ##
  ##   `@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG`
  ##   `Instrument:Run:FlowCell:Lane:Tile:TileX:TileY Number:Filter:Control:Barcode`
  ##
  ## - For Illumina < 1.4
  ##
  ##   `@HWUSI-EAS100R:6:73:941:1973#0/1`
  ##   `Instrument:FlowCell:Tile:TileX:TileY#Index/Pairing`
  ##
  ## - For Ncbi Sra
  ##
  ##   `@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=72`
  ##   `Id(discarded) Instrument:FlowCell:Tile:TileX:TileY (discarded)`
  ##
  ## If your tag doesn't match any of the above, you get an empty
  ## `Table<https://nim-lang.org/docs/tables.html>`_. This is the default if no
  ## `tagName<#tagName>`_ is specified.
  ##
  runnableExamples:
    import tables

    let meta = parseTag("HWUSI-EAS100R:6:73:941:1973#0/1", pnIlluminaOld)

    doAssert meta["instrumentId"].metaString == "HWUSI-EAS100R"
    doAssert meta["pairing"].metaInt == 1

  runnableExamples:
    import tables

    let meta = parseTag("EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG",
                        pnIllumina)

    doAssert meta["barcode"].metaString == "ATCACG"

  var instrumentId, runId, flowCellId, index, filtered, empty: string
  var flowCellLane, tileNumber, xTile, yTile, pairing, readNumber: int
  case platformName
  of pnIllumina:
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
  of pnIlluminaOld:
    if scanf(tag, "${alnum}:$i:$i:$i:$i#${bases}/$i", instrumentId, flowCellLane,
             tileNumber, xTile, yTile, index, pairing):
      result["instrumentId"] = MetaObj(kind: mkString, metaString: instrumentId)
      result["flowCellLane"] = MetaObj(kind: mkInt, metaInt: flowCellLane)
      result["tileNumber"] = MetaObj(kind: mkInt, metaInt: tileNumber)
      result["xTile"] = MetaObj(kind: mkInt, metaInt: xTile)
      result["yTile"] = MetaObj(kind: mkInt, metaInt: yTile)
      result["index"] = MetaObj(kind: mkString, metaString: index)
      result["pairing"] = MetaObj(kind: mkInt, metaInt: pairing)
  of pnNcbiSra:
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
  ## The quality line is stored at `meta<sequences.html#MetaObj>`_ property.
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
  ##   import bio/fastq
  ##
  ##   var strm = newGzFileStream("path/to/file.fasq.gz")
  ##
  ##   for sequence in sequences(strm):
  ##     doAssert(sequence of SequenceRecord)
  ##
  ## Accessing the quality of the `SequenceRecord` is done through `tables`:
  ##
  ## .. code-block::
  ##
  ##   import streams
  ##   import tables
  ##   import bio/fastq
  ##
  ##   var strm = newFileStream("path/to/file.fasq")
  ##
  ##   for sequence in sequences(strm):
  ##     echo sequence.meta.getOrDefault("quality")
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
