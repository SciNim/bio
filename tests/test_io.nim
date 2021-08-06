import os
import posix_utils
import streams
import unittest

import bio / [fastq, io]

suite "IO Operations with FASTQ files":
  setup:
    let fastqFM: string = currentSourcePath().parentDir / "test_files" /
      "two_multiline.fsq"
    var tmpOutput: (string, File) = mkstemp("tests/file_")
    let fasta: string = currentSourcePath().parentDir / "test_files" /
      "regular.fas"

  teardown:
    defer:
      tmpOutput[1].close
      removeFile(tmpOutput[0])

  test "Read Fasta from a filename":
    var records: seq[SequenceRecord]
    for s in sequences(fasta, ftFasta):
      records.add s

    check len(records) == 2
    check records[0].name == "First sequence"

  test "Read Fasta from a stream":
    let strm = newFileStream(fasta)

    var records: seq[SequenceRecord]
    for s in sequences(strm, ftFasta):
      records.add s

    check len(records) == 2
    check records[0].name == "First sequence"

  test "Write to a stream using the io interface":
    var strm = newFileStream(tmpOutput[0], fmWrite)

    var records: seq[SequenceRecord]
    for s in sequences(fastqFM, ftFastq):
      records.add s

    records.dumpTo(strm, ftFastQ)

    let dataSaved = tmpOutput[0].readLines(5)

    check dataSaved[0] == "@SEQ_ID1"
    check dataSaved[1][.. 20] == "GATTTGGGGTTCAAAGCAGTA"
    check dataSaved[2][.. 20] == "GATTTGGGGTTCAAAGCAGTA"
    check dataSaved[3] == "+"
    check dataSaved[4][.. 20] == "+''*((((***+))%%%++)("
