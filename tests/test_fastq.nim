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

    let metaQuality = "!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCFFF"

    var expMeta = Meta(kind: mkString, metaString: metaQuality)

    for s in sequences(strm):
      records.add(s)

    check records.len == 2
    check records[0].name == "SEQ_ID"
    check records[0].record.chain.len == 60
    check records[0].record.class == scDna
    check records[0].meta["quality"].kind == mkString
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

  test "Load the Illumina / Solexa tags into meta data":
    check false

  test "Operations on quality code":
    check false
