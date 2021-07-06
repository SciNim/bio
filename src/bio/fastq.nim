## :Author: |author|
## :Version: |libversion|
##
## A FASTQ file `normally`_ uses four lines per sequence, but multilines for
## sequence and quality "lines" are allowed:
##
## * Line 1 begins with a `@` character and is followed by a sequence identifier
##   and an optional description (like a FASTA title line).
## * Line 2 is the raw sequence of letters.
## * Line 3 begins with a `+` character and is *optionally* followed by the same
##   sequence identifier (and any description) again. Some parsers enforce that
##   the optional text here match the contents of line 1, but this one ignores
##   the line entirely.
## * Line 4 encodes the quality values for the sequence in Line 2, and must
##   contain the same number of symbols as letters in the sequence.
##
## .. _normally: https://en.wikipedia.org/wiki/FASTQ_format

import math
import sequtils
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

func getErrors*(scores: seq[int8], platform: PlatformName = pnNone): seq[float] =
  ## Calculates the Probability of Error (Pe) of each base from the sequence of
  ## Phred scores encoded as Ascii `chars`_ (*i.e. not shifted*).
  ##
  ## The default score encoding is a char between 33 and 126 (printable chars),
  ## with 33 the worst base quality (the base is as good as choosen randomly)
  ## and 126 the best. This encoding is referred to as `Sanger`.
  ##
  ## In 2004 Solexa introduced their own scoring system. They used chars
  ## between 59 and 126 and a different formula. Illumina older than 1.3 uses
  ## this format, specified here with the tag `pnIlluminaOld`.
  ##
  ## From Illumina 1.3 to 1.7, the minimum value for a Qscore is 64, and the
  ## formula to calculate the Phred changed again. Use `pnIllumina` for those
  ## sequences.
  ##
  ## Both scales are nearly equivalent for Phred scores above 10 (probability
  ## of a base being wrong < 0.1), as shown here:
  ##
  ## .. image:: PhredVsError.svg
  ##
  ## From Illumina 1.8, they returned to the default scores. Use `pnNone` or
  ## nothing for those sequences.
  ##
  ## To recap, a `D` in a quality string can mean: a high quality base (Q35) in
  ## the default Sanger scale, a low quality (Q4) in Illumina pre-1.8, and a
  ## again a high quality (Q35) in Illumina 1.8 or newer.
  ##
  ## More about this wonderful events at https://doi.org/10.1093/nar/gkp1137
  ## and https://en.wikipedia.org/wiki/FASTQ_format#Variations
  ##
  ## .. _chars: https://nim-lang.org/docs/system.html#ord%2CT
  ##
  runnableExamples:
    import std / [math, sugar]
    let scores = @[33'i8, 40, 50, 60, 80, 126]

    let errors = collect(newSeqOfCap(scores.len)):
      for s in getErrors(scores):
        round(s, 5)

    doAssert errors == @[1.0, 0.19953, 0.01995, 0.002, 2e-05, 0.0]

  runnableExamples:
    import std / [math, sugar]
    let scores = @[59'i8, 64, 70, 80, 90, 126]

    let errors = collect(newSeqOfCap(scores.len)):
      for s in getErrors(scores, pnIlluminaOld):
        round(s, 5)

    doAssert errors == @[0.75975, 0.5, 0.20076, 0.0245, 0.00251, 0.0]

  runnableExamples:
    import std / [math, sugar]
    let scores = @[64'i8, 80, 90, 126]

    let errors = collect(newSeqOfCap(scores.len)):
      for s in getErrors(scores, pnIllumina):
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
        clamp(pow(10, -(i - offset).float / 10.0), 0.0, 1.0)

func getQstring*(values: seq[float], platform: PlatformName = pnNone): string =
  ## Returns the Qstring for a given sequence of errors expressed as floats
  ## between 0.0 and 1.0
  ##
  ## Refer to `getErrors`_ for more info about differences between platforms.
  ##
  ## .. _getErrors: #getErrors,seq[int8],PlatformName
  ##
  runnableExamples:
    let errors = @[1e-14, 0.5, 0.75, 1]

    doAssert getQstring(errors) == "~$\"!"
    doAssert getQstring(errors, pnIlluminaOld) == "~@;;"
    doAssert getQstring(errors, pnIllumina) == "~CA@"

  let offset = case platform
    of pnIlluminaOld, pnIllumina: 64
    else: 33

  let s = collect(newSeqOfCap(values.len)):
    for value in values:
      case platform
      of pnIlluminaOld:
        char(clamp(int(round(-10 * log10(value / (1 - value)))),
                   -5, 62) + offset)
      of pnIllumina:
        char(clamp(int(round(-10 * log10(value))), 0, 62) + offset)
      else:
        char(clamp(int(round(-10 * log10(value))), 0, 93) + offset)

  join(s)

func parseQuality*(quality: string): seq[int8] =
  ## Parses a quality string into a seq of signed integers of size 8.
  ##
  ## Doesn't take into account the differences between platforms. Refer to
  ## `getErrors`_ for more info.
  ##
  ## .. _getErrors: #getErrors,seq[int8],PlatformName
  ##
  runnableExamples:
    let qualityString = "!'A~"

    doAssert parseQuality(qualityString) == @[33'i8, 39, 65, 126]

  collect(newSeqOfCap(len(quality))):
    for c in quality:
      int8(ord(c))

func qString*(quality: string, src, tgt: PlatformName): string =
  ## Transforms a quality string from the `src` format to the `tgt` format.
  ##
  ## Due to scale constrains, the transformation might cause loss of
  ## information.
  ##
  runnableExamples:
    let sanger = "ABCDEFGH"

    # Shifts about 30 chars in the Ascii table
    doAssert qString(sanger, pnNone, pnIllumina) == "`abcdefg"

  let rawQuality = parseQuality(quality)

  getQString(getErrors(rawQuality, src), tgt)

func parseTag*(tag: string, platformName: PlatformName = pnNone): Table[string, MetaObj] =
  ## Parses the identifier into its parts.
  ##
  ## There are some docs explaining how this identifier work, and none of them
  ## seems to agree. This procedure attempts the following:
  ##
  ## - For Illumina < 1.4
  ##
  ##   `@HWUSI-EAS100R:6:73:941:1973#0/1`
  ##   `Instrument:FlowCell:Tile:TileX:TileY#Index/Pairing`
  ##
  ## - For Illumina > 1.4
  ##
  ##   `@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG`
  ##   `Instrument:Run:FlowCell:Lane:Tile:TileX:TileY Number:Filter:Control:Barcode`
  ##
  ## - For Ncbi Sra
  ##
  ##   `@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=72`
  ##   `Id(discarded) Instrument:FlowCell:Tile:TileX:TileY (discarded)`
  ##
  ## If the identifier doesn't match any of the above, you get an empty
  ## `Table<https://nim-lang.org/docs/tables.html>`_. This is the default if no
  ## `platformName<#platformName>`_ is specified.
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
  ##  Note: when compiled with `-d:danger`, some checks implemented as
  ##  assertions do not run.
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
  ## Accessing the quality of the `SequenceRecord` is done through `tables`_:
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
  ## .. _tables: https://nim-lang.org/docs/tables.html
  ##
  var name, sequence, quality, line: string
  var loadQuality: bool

  template addQual(line: typed): untyped =
    assert line.allIt(ord(it) > 32 and ord(it) < 127), "Invalid quality character"
    quality.add line

  while strm.readLine(line):
    if (line.startsWith('@') and
        len(sequence) > 0 and len(quality) == len(sequence)) or strm.atEnd:
      if strm.atEnd:
        # Once reached the EOF, check the last record is correct. File can be
        # truncated.
        addQual line
        assert len(quality) == len(sequence)
        assert not sequence.isEmptyOrWhitespace
      let meta = {"quality": MetaObj(kind: mkString, metaString: quality)}.toTable
      yield SequenceRecord(name: name,
                           record: guess(sequence.toUpperAscii),
                           meta: meta)
      loadQuality = false
      name = ""
      sequence = ""
      quality = ""

    if loadQuality:  # '+' was detected, we are in quality lines.
      if len(quality) < len(sequence):
        addQual line
      else:
        raise newException(AssertionDefect, "Extra quality data found.")
    else:  # Is either name, sequence or the '+' separator line
      if line.startsWith('@'):  # Is the name line
        assert name.isEmptyOrWhitespace  # Or we have a double @ line
        name = line[1 .. ^1]
      elif line.startsWith('+'):  # The '+' separator
        loadQuality = true
        quality = newStringOfCap(len(sequence))
      else:  # Is sequence
        sequence.add line

proc dumpTo*(record: SequenceRecord, strm: Stream, kind: FileType=ftFastq) =
  ## Write a `SequenceRecord<sequences.html#SequenceRecord>`_ to `fHandler`,
  ## wrapping the sequence by 60 positions.
  ##
  ## The identifier used to save is the `name` of the SequenceRecord, and not
  ## constructed from the `Metadata<sequences.html#MetaObj>`.
  ##
  ## The quality line comes from the `Metadata` value for `quality`.
  ##
  ## .. code-block::
  ##   `import bio/fastq`
  ##
  ##   let mySeq = newDna("TGCACCCCA")
  ##   var myRec = SequenceRecord(name: "My DNA Read", record: mySeq)
  ##   let quality = MetaObj(kind: mkString, metaString: "ABCDEFGHI")
  ##   myRec.meta["quality"] = quality
  ##
  ##   block:
  ##     let fastqOut = newFileStream("somefile.fastq", fmWrite)
  ##     defer: fastqOut.close
  ##     myRec.dumpTo(fastqOut)
  ##

  const wrapSize: int = 60

  strm.write("@", record.name)
  for i, base in record.record.chain.pairs:
    if i mod wrapSize == 0:
      strm.write('\n')
    strm.write(base)
  strm.write('\n')
  strm.write("+")
  for i, qVal in record.meta.getOrDefault("quality").metaString.pairs:
    if i mod wrapSize == 0:
      strm.write('\n')
    strm.write(qVal)
  strm.write('\n')
  strm.flush()

proc dumpTo*(records: seq[SequenceRecord], strm: Stream, kind: FileType=ftFastq) =
  ## A shortcut to avoid the explicit cycle to write a `seq` of
  ## `SequenceRecords<sequences.html#SequenceRecord>`_
  ##
  ## .. code-block::
  ##   `import bio/fastq`
  ##
  ##   let mySeq = newDna("TGCACCCCA")
  ##   var myRecA = SequenceRecord(name: "My DNA Read1", record: mySeq)
  ##   var myRecB = SequenceRecord(name: "My DNA Read2", record: mySeq)
  ##   let quality = MetaObj(kind: mkString, metaString: "ABCDEFGHI")
  ##   myRecA.meta["quality"] = quality
  ##   myRecB.meta["quality"] = quality
  ##
  ##   let myRecs = @[myRecA, myRecB]
  ##
  ##   block:
  ##     let fastqOut = newFileStream("somefile.fastq", fmWrite)
  ##     defer: fastqOut.close
  ##     myRecs.dumpTo(fastqOut)
  ##
  for sr in records:
    sr.dumpTo(strm, kind)

proc dumpTo*(record: SequenceRecord, fName: string, kind: FileType=ftFasta) =
  ## Same as `write-through-handler proc<#dumpTo,SequenceRecord,Stream,string>`_
  ## but you only need to point out the name of the file.
  ##
  ## If the file exists, it will be silently overwritten.
  ##
  ## .. code-block::
  ##   `import bio/fastq`
  ##
  ##   let mySeq = newDna("TGCACCCCA")
  ##   var myRec = SequenceRecord(name: "My DNA Read", record: mySeq)
  ##   let quality = MetaObj(kind: mkString, metaString: "ABCDEFGHI")
  ##   myRec.meta["quality"] = quality
  ##
  ##   block:
  ##     let fastqOut = newFileStream("somefile.fastq", fmWrite)
  ##     defer: fastqOut.close
  ##     myRec.dumpTo("somefile.fastq")
  ##
  let strm: Stream = newFileStream(fName, fmWrite)
  defer: strm.close()

  record.dumpTo(strm, kind)

proc dumpTo*(records: seq[SequenceRecord], fName: string, kind: FileType=ftFasta) =
  ## Same as `write-through-handler proc<#dumpTo,SequenceRecord,Stream,string>`_
  ## but you only need to point out the name of the file.
  ##
  ## If the file exists, it will be silently overwritten.
  ##
  ## .. code-block::
  ##   import bio / fastq
  ##
  ##   let mySeq = newDna("TGCACCCCA")
  ##   var myRecA = SequenceRecord(name: "My DNA Read1", record: mySeq)
  ##   var myRecB = SequenceRecord(name: "My DNA Read2", record: mySeq)
  ##   let quality = MetaObj(kind: mkString, metaString: "ABCDEFGHI")
  ##   myRecA.meta["quality"] = quality
  ##   myRecB.meta["quality"] = quality
  ##
  ##   let myRecs = @[myRecA, myRecB]
  ##
  ##   myRecs.dumpTo("somefile.fastq")
  ##
  let strm: Stream = newFileStream(fName, fmWrite)
  defer: strm.close()

  for sr in records:
    sr.dumpTo(strm, kind)

iterator sequences*(fName: string, kind: FileType=ftFastq): SequenceRecord {.inline.} =
  ## Iterate through all the `Sequences<sequences.html#Sequence>`_ in a given
  ## filename, yielding `SequenceRecords<sequences.html#SequenceRecord>`_.
  ##
  ## .. code-block::
  ##
  ##   import bio/fastq
  ##
  ##   for sequence in sequences("path/to/file.fasq"):
  ##     doAssert(sequence of SequenceRecord)
  ##
  let strm = newFileStream(fName)

  for sr in strm.sequences(ftFastq):
    yield sr

iterator sequences*(strm: Stream, platform: PlatformName): SequenceRecord {.inline.} =
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
  ##   for sequence in sequences(strm, pnIllumina):
  ##     doAssert(sequence of SequenceRecord)
  ##
  for sr in strm.sequences(ftFastq):
    case platform
    of pnIlluminaOld:
      assert allIt(sr.meta["quality"].metaString, ord(it) >= 59)
    of pnIllumina:
      assert allIt(sr.meta["quality"].metaString, ord(it) >= 64)
    else:
      discard
    yield sr

iterator sequences*(fName: string, platform: PlatformName): SequenceRecord {.inline.} =
  ## Iterate through all the `Sequences<sequences.html#Sequence>`_ in a given
  ## filename, yielding `SequenceRecords<sequences.html#SequenceRecord>`_,
  ## *enforcing* the quality line to comply with the Platform requested.
  ##
  ## .. code-block::
  ##
  ##   import bio/fastq
  ##
  ##   for sequence in sequences("path/to/file.fasq", pnIlluminaOld):
  ##     doAssert(sequence of SequenceRecord)
  ##
  let strm = newFileStream(fName)

  for sr in strm.sequences(platform):
    yield sr
