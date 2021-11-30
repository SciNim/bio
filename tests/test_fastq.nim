import algorithm
import os
import posix_utils
import streams
import strformat
import strtabs
import strutils
import tables
import unittest

import bio / [fastq, io]


proc `$`(m: MetaObj): string =
  case m.kind
  of mkInt:
    return $m.kind & ": " & $m.metaInt
  of mkFloat:
    return $m.kind & ": " & $m.metaFloat
  of mkSeqInt8:
    return $m.kind & ": " & $m.metaSeqInt8
  of mkString:
    return $m.kind & ": " & $m.metaString
  of mkTableString:
    return $m.kind & ": " & $m.metaTableString

proc `==`(x, y: Table[string, MetaObj]): bool =
  if len(x) != len(y): return false
  for k, v in x.pairs:
    if v.kind != y[k].kind:
      return false
    if $v != $y[k]:
      return false
  return true

suite "Operations with FASTQ files":
  setup:
    let fastqF: string = currentSourcePath().parentDir / "test_files" /
      "one_sequence.fsq"
    let fastqF2: string = currentSourcePath().parentDir / "test_files" /
      "two_sequence.fsq"
    let fastqFM: string = currentSourcePath().parentDir / "test_files" /
      "two_multiline.fsq"
      # one_illumina_new.fsq, one_illumina_old.fsq, two_nbcis.fsq,
      # two_sequence.fsq

    let emptyTable = initTable[string, MetaObj]()

  test "Read records from a stream of one record":
    let strm = newFileStream(fastqF)
    var records: seq[SequenceRecord]

    let metaQuality = "!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65"

    var expMeta = MetaObj(kind: mkString, metaString: metaQuality)

    for s in sequences(strm):
      records.add(s)

    check records.len == 1
    check records[0].name == "SEQ_ID"
    check records[0].chain.len == 60
    check records[0].class == scDna
    check records[0].meta["quality"].kind == mkString
    check records[0].meta["quality"].metaString == expMeta.metaString

  test "Read records from a stream of two records":
    let strm = newFileStream(fastqF2)
    var records: seq[SequenceRecord]

    let metaQuality = "@''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCFFF"

    var expMeta = MetaObj(kind: mkString, metaString: metaQuality)

    for s in sequences(strm):
      records.add(s)

    check records.len == 2
    check records[0].name == "SEQ_ID1"
    check records[0].chain.len == 60
    check records[0].class == scDna
    check records[0].meta["quality"].kind == mkString
    check records[0].meta["quality"].metaString.len == 60
    check records[1].name == "SEQ_ID2"
    check records[1].chain.len == 60
    check records[1].class == scDna
    check records[1].meta["quality"].kind == mkString
    check records[1].meta["quality"].metaString == expMeta.metaString

  test "Read an empty file":
    let strm = newFileStream(currentSourcePath().parentDir / "test_files" / "empty.fas")
    var records: seq[SequenceRecord]

    for s in sequences(strm):
      records.add s

    check records.len == 0

  test "Load records from a filename":
    var records: seq[SequenceRecord]
    for s in sequences(fastqF2):
      records.add s

    check records.len == 2

  test "Read records from a stream of two multiline records":
    let strm = newFileStream(fastqFM)
    var records: seq[SequenceRecord]

    for s in sequences(strm):
      records.add(s)

    check records.len == 2
    check records[0].name == "SEQ_ID1"
    check records[0].chain.len == 120
    check records[0].meta["quality"].metaString.len == 120
    check records[1].name == "SEQ_ID2"
    check records[1].chain.len == 120
    check records[1].meta["quality"].metaString.len == 120

  test "Load the Illumina / Solexa new tags into meta data":
    let newTag = "EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG"
    let expected = {
      "instrumentId": MetaObj(kind: mkString, metaString: "EAS139"),
      "runId": MetaObj(kind: mkString, metaString: "136"),
      "flowCellId": MetaObj(kind: mkString, metaString: "FC706VJ"),
      "flowCellLane": MetaObj(kind: mkInt, metaInt: 2),
      "tileNumber": MetaObj(kind: mkInt, metaInt: 2104),
      "xTile": MetaObj(kind: mkInt, metaInt: 15343),
      "yTile": MetaObj(kind: mkInt, metaInt: 197393),
      #
      "readNumber": MetaObj(kind: mkInt, metaInt: 1),
      "isFiltered": MetaObj(kind: mkString, metaString: "Y"),
      "control": MetaObj(kind: mkInt, metaInt: 18),
      "barcode": MetaObj(kind: mkString, metaString: "ATCACG")}.toTable
    check parseTag(newTag, pnIllumina) == expected

    # Assert also that this is the default mode
    check parseTag(newTag) == emptyTable

  test "Save a sequence of SequenceRecord to a file":
    var tmpOutput: (string, File) = mkstemp("tests/file_")
    defer:
      tmpOutput[1].close()
      removeFile(tmpOutput[0])

    var strm = newFileStream(tmpOutput[0], fmWrite)

    var records: seq[SequenceRecord]
    for s in sequences(fastqF2):
      records.add s

    records.dumpTo(strm)

    let dataSaved = tmpOutput[0].readLines(5)

    check dataSaved[0] == "@SEQ_ID1"
    check dataSaved[1][.. 20] == "GATTTGGGGTTCAAAGCAGTA"
    check dataSaved[2] == "+"
    check dataSaved[3][.. 20] == "+''*((((***+))%%%++)("

  test "Save a single SequenceRecord to a file":
    var tmpOutput: (string, File) = mkstemp("tests/file_")
    defer:
      tmpOutput[1].close()
      removeFile(tmpOutput[0])

    var strm = newFileStream(tmpOutput[0], fmWrite)

    var records: seq[SequenceRecord]
    for s in sequences(fastqF2):
      records.add s
      break

    records[0].dumpTo(strm)

    let dataSaved = tmpOutput[0].readLines(4)

    check dataSaved[0] == "@SEQ_ID1"
    check dataSaved[1][.. 20] == "GATTTGGGGTTCAAAGCAGTA"
    check dataSaved[2] == "+"
    check dataSaved[3][.. 20] == "+''*((((***+))%%%++)("

  test "Save a sequence of SequenceRecord to a file (using name)":
    var tmpOutput: (string, File) = mkstemp("tests/file_")
    defer:
      tmpOutput[1].close()
      removeFile(tmpOutput[0])

    var records: seq[SequenceRecord]
    for s in sequences(fastqF2):
      records.add s

    records.dumpTo(tmpOutput[0])

    let dataSaved = tmpOutput[0].readLines(5)

    check dataSaved[0] == "@SEQ_ID1"
    check dataSaved[1][.. 20] == "GATTTGGGGTTCAAAGCAGTA"
    check dataSaved[2] == "+"
    check dataSaved[3][.. 20] == "+''*((((***+))%%%++)("

  test "Save a single SequenceRecord to a file (using name)":
    var tmpOutput: (string, File) = mkstemp("tests/file_")
    defer:
      tmpOutput[1].close()
      removeFile(tmpOutput[0])

    var records: seq[SequenceRecord]
    for s in sequences(fastqF2):
      records.add s
      break

    records[0].dumpTo(tmpOutput[0])

    let dataSaved = tmpOutput[0].readLines(4)

    check dataSaved[0] == "@SEQ_ID1"
    check dataSaved[1][.. 20] == "GATTTGGGGTTCAAAGCAGTA"
    check dataSaved[2] == "+"
    check dataSaved[3][.. 20] == "+''*((((***+))%%%++)("

  test "Empty barcode":
    let newTag = "EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:1"
    let expected = {
      "instrumentId": MetaObj(kind: mkString, metaString: "EAS139"),
      "runId": MetaObj(kind: mkString, metaString: "136"),
      "flowCellId": MetaObj(kind: mkString, metaString: "FC706VJ"),
      "flowCellLane": MetaObj(kind: mkInt, metaInt: 2),
      "tileNumber": MetaObj(kind: mkInt, metaInt: 2104),
      "xTile": MetaObj(kind: mkInt, metaInt: 15343),
      "yTile": MetaObj(kind: mkInt, metaInt: 197393),
      #
      "readNumber": MetaObj(kind: mkInt, metaInt: 1),
      "isFiltered": MetaObj(kind: mkString, metaString: "Y"),
      "control": MetaObj(kind: mkInt, metaInt: 18),
      "barcode": MetaObj(kind: mkString, metaString: "1")}.toTable
    check parseTag(newTag, pnIllumina) == expected

  test "Load the Illumina / Solexa old tags into meta data":
    let oldTag = "HWUSI-EAS100R:6:73:941:1973#ATCACG/1"

    let expected = {
      "instrumentId": MetaObj(kind: mkString, metaString: "HWUSI-EAS100R"),
      "flowCellLane": MetaObj(kind: mkInt, metaInt: 6),
      "tileNumber": MetaObj(kind: mkInt, metaInt: 73),
      "xTile": MetaObj(kind: mkInt, metaInt: 941),
      "yTile": MetaObj(kind: mkInt, metaInt: 1973),
      "index": MetaObj(kind: mkString, metaString: "ATCACG"),
      "pairing": MetaObj(kind: mkInt, metaInt: 1)}.toTable

    check parseTag(oldTag, pnIlluminaOld) == expected

  test "Load the NCBI SRA tags into meta data":
    let sraTag = "SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36"

    let expected = {
      "instrumentId": MetaObj(kind: mkString, metaString: "071112_SLXA-EAS1_s_7"),
      "flowCellLane": MetaObj(kind: mkInt, metaInt: 5),
      "tileNumber": MetaObj(kind: mkInt, metaInt: 1),
      "xTile": MetaObj(kind: mkInt, metaInt: 817),
      "yTile": MetaObj(kind: mkInt, metaInt: 345)}.toTable

    check parseTag(sraTag, pnNcbiSra) == expected

  test "Quality code to int8 values":
    let quality = "!''*((((***+))%%%++)(%%%%~"

    let qualities = parseQuality(quality)
    check qualities[0] == 33'i8
    check qualities[^1] == 126'i8

  test "Get the scores from a quality string":
    let quality = "!&~"

    let pes = getErrors(parseQuality(quality))  # The default is plain Phred
    check pes[0] == 1 # The lowest quality, !
    check pes[1] > 0.31 and pes[1] < 0.32 # Very low quality of Phred 5, &
    check pes[^1] < 1e-9  # The highest quality, ~

  test "Get the scores from a quality string, Solexa/Old Illumina":
    let quality = ";@~"

    let pes = getErrors(parseQuality(quality), pnIlluminaOld)
    check pes[0] > 0.75 and pes[0] < 0.76  # The lowest quality, ;
    check pes[1] == 0.5 # Very low quality of Phred 5, @
    check pes[^1] < 1e-6  # The highest quality, ~

  test "Get the scores from a quality string, New Illumina":
    let quality = "@E~"

    let pes = getErrors(parseQuality(quality), pnIllumina)
    check pes[0] == 1 # The lowest quality, @
    check pes[1] > 0.31 and pes[1] < 0.32 # Very low quality of Phred 5, E
    check pes[^1] < 1e-6  # The highest quality, ~

  test "Don't allow errors greater than 1.0":
    check getErrors(@[0'i8]) == @[1.0]
    check getErrors(@[0'i8], pnIllumina) == @[1.0]

  test "Get the Qstring from a seq of errors":
    let values = @[1e-14, 0.5, 0.75, 1]

    check getQstring(values) == "~$\"!" # The default is plain Phred
    check getQstring(values, pnIlluminaOld) == "~@;;"
    check getQstring(values, pnIllumina) == "~CA@"

  test "Translate a Qstring from Default mode to Illumina Old, and back":
    let defaultString = "!\"#$%&'()~"
    let oldIlluminaString = ";;>@BCEFG~"

    check qString(defaultString, pnNone, pnIlluminaOld) == oldIlluminaString
    check qString(oldIlluminaString, pnIlluminaOld, pnNone) == "\"\"#$%&'()_"

  test "Translate a Qstring from Default mode to Illumina New, and back":
    let defaultString = "!\"#$%&'()~"
    let newIlluminaString = "@ABCDEFGH~"

    check qString(defaultString, pnNone, pnIllumina) == newIlluminaString
    check qString(newIlluminaString, pnIllumina, pnNone) == "!\"#$%&'()_"

  test "Reversing a Sequence Record also reverses Quality string":
    let strm = newFileStream(fastqF)
    var records: seq[SequenceRecord]

    for s in sequences(strm):
      records.add(s)

    let revCompRecord = reverseComplement(records[0])
    check revCompRecord.meta["quality"].metaString ==
      join(reversed(records[0].meta["quality"].metaString))

suite "File set test":
  ## Files provided in
  ## P.J.A. Cock, C.J. Fields, N. Goto, M.L. Heuer and P.M. Rice (2009).
  ## The Sanger FASTQ file format for sequences with quality scores, and the
  ## Solexa/Illumina FASTQ variants.
  ##
  ## Ignored the "error_diff_ids", as the + line content is optional
  ##
  test "Expect to be correct":
    for k, f in walkDir(currentSourcePath().parentDir / "test_files" / "fastq" /
                        "right"):

      var records: seq[SequenceRecord]
      if k == pcFile:
        for s in sequences(f):
          records.add s

  test "Expected to be wrong":
    let wrongDir = currentSourcePath().parentDir / "test_files" / "fastq" /
      "wrong"
    test "Double Quality or Id line":
      for fName in ["qual", "seq"]:
        let f = wrongDir / &"error_double_{fName}.fastq"

        expect AssertionDefect:
          for s in sequences(f):
            discard

    test "Quality data != Sequence line":
      for fName in ["long", "short", "no"]:
        let f = wrongDir / &"error_{fName}_qual.fastq"

        expect AssertionDefect:
          for s in sequences(f):
            discard

    test "Non-printable chars in quality line":
      for fName in ["del", "escape", "space", "tab", "unit_sep", "vtab"]:
        let f = wrongDir / &"error_qual_{fName}.fastq"

        expect AssertionDefect:
          for s in sequences(f):
            discard

      for fName in ["spaces", "tabs"]:
        let f = wrongDir / &"error_{fName}.fastq"

        expect AssertionDefect:
          for s in sequences(f):
            discard

    test "Truncated files i.e. after downloading":
      for fName in ["at_plus", "at_qual", "at_seq", "in_plus", "in_qual",
                    "in_seq", "in_title"]:
        let f = wrongDir / &"error_trunc_{fName}.fastq"

        expect AssertionDefect:
          for s in sequences(f):
            discard

  test "Expected to be correct in some conditions":
    let srcDir = currentSourcePath().parentDir / "test_files" / "fastq" /
      "right"

    test "Full range Sanger read as Sanger, not Illumina":
      let f = srcDir / "sanger_full_range_original_sanger.fastq"

      # Ensure the file can be read under Sanger conditions
      for s in sequences(f, pnNone):
        discard

      expect AssertionDefect:
        for s in sequences(f, pnIlluminaOld):
          discard

      expect AssertionDefect:
        for s in sequences(f, pnIllumina):
          discard

    test "Full range IlluminaOld read as Sanger or IlluminaOld, not Illumina":
      let f = srcDir / "solexa_full_range_original_solexa.fastq"

      # Ensure the file can be read under Sanger conditions
      for s in sequences(f, pnNone):
        discard

      # Ensure the file can be read under IlluminaOld conditions
      for s in sequences(f, pnIlluminaOld):
        discard

      expect AssertionDefect:
        for s in sequences(f, pnIllumina):
          discard

    test "Full range Illumina read as Sanger, IlluminaOld and Illumina":
      let f = srcDir / "illumina_full_range_original_illumina.fastq"

      # Ensure the file can be read under Sanger conditions
      for s in sequences(f, pnNone):
        discard

      # Ensure the file can be read under IlluminaOld conditions
      for s in sequences(f, pnIlluminaOld):
        discard

      # Ensure the file can be read under IlluminaNew conditions
      for s in sequences(f, pnIllumina):
        discard
