import os
import streams
import tables
import unittest

import bio / fastq


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

  test "Read records from a stream of one record":
    let strm = newFileStream(fastqF)
    var records: seq[SequenceRecord]

    let metaQuality = "!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65"

    var expMeta = Meta(kind: mkString, metaString: metaQuality)

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

    var expMeta = Meta(kind: mkString, metaString: metaQuality)

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
    let oldTag = "HWUSI-EAS100R:6:73:941:1973#0/1"
    check parseTag(oldTag, tnIllumina) == true

    # Assert also that this is the default mode
    check parseTag(oldTag) ==

  test "Load the Illumina / Solexa old tags into meta data":
    let newTag = "EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG"
    check parseTag(newTag, tnIlluminaOld)

  test "Load the NCBI SRA tags into meta data":
    let sraTag = "SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36"
    check parseTag(sraTag, tnNcbiSra)

  test "Operations on quality code":
    check false
