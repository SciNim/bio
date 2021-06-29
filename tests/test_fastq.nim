import os
import streams
import strtabs
import tables
import unittest

import bio / fastq


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
    check records[0].record.chain.len == 60
    check records[0].record.class == scDna
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
    check records[0].record.chain.len == 60
    check records[0].record.class == scDna
    check records[0].meta["quality"].kind == mkString
    check records[0].meta["quality"].metaString.len == 60
    check records[1].name == "SEQ_ID2"
    check records[1].record.chain.len == 60
    check records[1].record.class == scDna
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
    check records[0].record.chain.len == 120
    check records[0].meta["quality"].metaString.len == 120
    check records[1].name == "SEQ_ID2"
    check records[1].record.chain.len == 120
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
    check parseTag(newTag, tnIllumina) == expected

    # Assert also that this is the default mode
    check parseTag(newTag) == emptyTable

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
    check parseTag(newTag, tnIllumina) == expected

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

    check parseTag(oldTag, tnIlluminaOld) == expected

  test "Load the NCBI SRA tags into meta data":
    let sraTag = "SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36"

    let expected = {
      "instrumentId": MetaObj(kind: mkString, metaString: "071112_SLXA-EAS1_s_7"),
      "flowCellLane": MetaObj(kind: mkInt, metaInt: 5),
      "tileNumber": MetaObj(kind: mkInt, metaInt: 1),
      "xTile": MetaObj(kind: mkInt, metaInt: 817),
      "yTile": MetaObj(kind: mkInt, metaInt: 345)}.toTable

    check parseTag(sraTag, tnNcbiSra) == expected

  test "Operations on quality code":
    check false
