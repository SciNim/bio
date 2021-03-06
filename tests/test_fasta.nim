import std / [os, posix_utils, streams, strformat, strutils, tables, unittest]

import zip / gzipfiles

import bio / [fasta, io]


suite "Test SequenceRecord operation":
  setup:
    #var dna = Sequence(chain: "ACGTGGGGT", class: scDna)
    var dnaRecord = SequenceRecord(
      chain: "ACGTGGGGT", class: scDna, name: "Sample")
    var tmpOutput: (string, File) = mkstemp("tests/file_")

  teardown:
    tmpOutput[1].close
    removeFile(tmpOutput[0])

  test "Write records as a FASTA through filename":
    dnaRecord.dumpTo(tmpOutput[0])

    check fileExists(tmpOutput[0])

    let fastaFile = tmpOutput[0].readLines(2)

    check fastaFile == @[">Sample", "ACGTGGGGT"]

  test "Write records as a FASTA through stream":
    let output = newFileStream(tmpOutput[1])
    dnaRecord.dumpTo(output)

    let fastaFile = tmpOutput[0].readLines(2)

    check fastaFile == @[">Sample", "ACGTGGGGT"]

  test "Automatically wrap the long lines to 60 chars":
    dnaRecord.chain = repeat(dnaRecord.chain, 8)
    let output = newFileStream(tmpOutput[1])
    dnaRecord.dumpTo(output)

    let fastaFile = tmpOutput[0].readLines(3)

    check fastaFile.len == 3
    check fastaFile[1].len == 60
    check fastaFile[1] == dnaRecord.chain[0 .. 59]

  test "Load records from FASTA through filename":
    let records: seq[SequenceRecord] = load(currentSourcePath().parentDir /
                                            "test_files" / "regular.fas")

    check records.len == 2
    check records[0].name == "First sequence"
    check records[1].name == "Second sequence"
    check records[0].chain.len == 120
    check records[1].chain.len == 120
    check records[0].class == scDna

  test "Load record from FASTA through iterator + filename":
    var records: seq[SequenceRecord]

    for s in sequences(currentSourcePath().parentDir / "test_files" / "regular.fas"):
      records.add s

    check records.len == 2
    check records[0].name == "First sequence"
    check records[1].name == "Second sequence"
    check records[0].chain.len == 120
    check records[1].chain.len == 120
    check records[0].class == scDna

  test "Load records from FASTA: special cases":
    test "The empty file: returns empty seq":
      var records: seq[SequenceRecord]

      for s in sequences(currentSourcePath().parentDir / "test_files" / "empty.fas"):
        records.add s

      check records.len == 0

    test "The lonely sequence: returns a seq with only one sequence":
      var records: seq[SequenceRecord]
      for s in sequences(currentSourcePath().parentDir / "test_files" / "one_sequence.fas"):
        records.add s

      check records.len == 1

    test "The lonely sequence but with extra blank line at end":
      var records: seq[SequenceRecord]
      for s in sequences(currentSourcePath().parentDir /
                         "test_files" / "one_sequence_extraline.fas"):
        records.add s

      check records.len == 1

  test "Load record slurping from a file":
    # File opening is not available when compiling, but you can slurp the file
    # into a seq[string] and then load that into a seq[SequenceRecords]
    const seqLines = slurp(currentSourcePath().parentDir / "test_files" /
                           "regular.fas")

    var records: seq[SequenceRecord]

    for r in sequences(seqLines.splitLines):
      records.add r

    check records.len == 2
    check records[0].name == "First sequence"
    check records[1].name == "Second sequence"
    check records[0].chain.len == 120
    check records[1].chain.len == 120
    check records[0].class == scDna

  test "Write a bunch of records as a FASTA through filehandler":
    let multiRecords: seq[SequenceRecord] = @[dnaRecord, dnaRecord]

    let output = newFileStream(tmpOutput[1])

    multiRecords.dumpTo(output)

    check fileExists(tmpOutput[0])

    var fastaFile: seq[string]
    for line in tmpOutput[0].lines:
      fastaFile.add line

    check fastaFile == @[">Sample", "ACGTGGGGT", ">Sample", "ACGTGGGGT"]

  test "Read records from a stream":
    let strm = newFileStream(currentSourcePath().parentDir / "test_files" / "regular.fas")
    var records: seq[SequenceRecord]

    for s in sequences(strm):
      records.add(s)

    check records.len == 2
    check records[0].name == "First sequence"
    check records[1].name == "Second sequence"
    check records[0].chain.len == 120
    check records[1].chain.len == 120
    check records[0].class == scDna

  test "Read records from a compressed file through streams":
    let strm = newGzFileStream(currentSourcePath().parentDir / "test_files" / "regular.fas.gz")
    var records: seq[SequenceRecord]

    for s in sequences(strm):
      records.add(s)

    check records.len == 2
    check records[0].name == "First sequence"
    check records[1].name == "Second sequence"
    check records[0].chain.len == 120
    check records[1].chain.len == 120
    check records[0].class == scDna

  test "Yielded record can be modified":
    for s in sequences(currentSourcePath().parentDir / "test_files" /
                       "regular.fas"):
      s.name = "New Name"

suite "Operations with FASTA files":
  setup:
    let fastaF: string = currentSourcePath().parentDir / "test_files" / "regular.fas"

  test "Index building":
    let tbl = {"First sequence": 0'i64, "Second sequence": 138'i64}.newTable
    let index: Index = Index(table: tbl)

    check newIndex(fastaF).table == tbl
    check newIndex(fastaF).source == fastaF

  test "Built index is pointing correctly":
    let f: File = open(fastaF)
    defer: f.close

    let index: Index = newIndex(fastaF)

    setFilePos(f, index.table["Second sequence"])
    check f.readChar == '>'

  test "Access the Indexes by sequence name":
    let index: Index = newIndex(fastaF)
    let expectedRec = SequenceRecord(
      name: "Second sequence",
      chain: "GGGGCATGCATCGACATACGCATCAGCAGACGACTACGACTCAGACTACGACTCAGCGGG" &
             "TTTGCATGCATCGACATACGCATCAGCAGACGACTACGACTCAGACTACGACTCAGCTTT",
      class: scDna)

    check index["Second sequence"].chain == expectedRec.chain
    check index["Second sequence"].class == expectedRec.class

  test "Use the index with a file handler":
    let index: Index = newIndex(fastaF)
    let expectedRec = SequenceRecord(
      name: "Second sequence",
      chain: "GGGGCATGCATCGACATACGCATCAGCAGACGACTACGACTCAGACTACGACTCAGCGGG" &
             "TTTGCATGCATCGACATACGCATCAGCAGACGACTACGACTCAGACTACGACTCAGCTTT",
      class: scDna)
    let fastaHandler = open(fastaF)
    defer: fastaHandler.close

    check index[fastaHandler, "Second sequence"].chain == expectedRec.chain
    check index[fastaHandler, "Second sequence"].class == expectedRec.class

  test "Convenience 'items' procedure to deal with the Index.table":
    let index: Index = newIndex(fastaF)

    var seqNames: seq[string]

    for s in index:
      seqNames.add s

    check seqNames == @["First sequence", "Second sequence"]

  test "Give the index an useful string":
    let index: Index = newIndex(fastaF)

    check $index == &"Index for {fastaF}, length: 2"
