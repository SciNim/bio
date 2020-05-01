import strutils
import ../bio


proc load*(fName: string, kind: string="fasta"): seq[SequenceRecord] =
  ## Load a sequence of files from a filename.
  ##
  ## If you need an interator over a file, use `sequences iterator<#sequences>`_
  runnableExamples:
    let mySeqs: SequenceRecord = load("test_files/regular_fasta.fas")

  let fileIn: File = open(fName)
  defer: fileIn.close

  let name = fileIn.readLine
  var sequence = Sequence()
  var seqRecord = SequenceRecord(name: name.strip()[1..^1], record: sequence)
  for line in fileIn.lines:
    if line[0] == '>':
      seqRecord.record.class = guess(seqRecord.record.chain)
      result.add seqRecord
      sequence = Sequence()
      seqRecord = SequenceRecord(name: line.strip()[1..^1], record: sequence)
    else:
      seqRecord.record.chain.add line.strip().toUpper()
  result.add seqRecord

proc write*(record: SequenceRecord, fHandler: File, kind: string="fasta") =
  ## Write a SequenceRecord to fHandler, wrapping the sequence by 60 positions.
  ## The name of the SequenceRecord remains untouched.
  ##
  ## TBD: `kind` is a string as in "fasta", to support different formats.
  ## Right now only FASTA files are supported.
  ##
  ## To write a FASTA file `myOutput.fasta` with the contents:
  ##
  ## .. code-block::
  ##   >My DNA sequence
  ##   TGCACCCCA
  ##
  ## use the following:
  runnableExamples:
    let fastaOut = open("myOutput.fasta", fmWrite)
    let mySeq = initDna("TGCACCCCA")
    let myRec = SequenceRecord(name: "My DNA sequence", record: mySeq)

    myRec.write(fastaOut)

  const wrapSize: int = 60
  fHandler.write(">", record.name)

  for i, base in record.record.chain.pairs:
    if i mod wrapSize == 0:
      fHandler.write("\n")
    fHandler.write(base)
  fHandler.write("\n")
  fHandler.flushFile

proc write*(record: SequenceRecord, fName: string, kind: string="fasta") =
  ## Same as `write-through-handler proc<#write,SequenceRecord,string,string>`_
  ## but you only need to point out the name of the file.
  ##
  runnableExamples:
    let mySeq = initDna("TGCACCCCA")
    let myRec = SequenceRecord(name: "My DNA sequence", record: mySeq)

    myRec.write("myOutput.fasta")

  let fHandler: File = open(fName, fmWrite)
  defer: fHandler.close()

  write(record, fHandler, kind)

iterator sequences*(fName: string, kind: string="fasta"):
  SequenceRecord {.inline.} =
  ## Iterate through all the sequences in a given filename
  runnableExamples:
    for sequence in sequences("test_files/regular_fasta.fas"):
      doAssert(sequence of SequenceRecord)

  let fileIn: File = open(fName)
  defer: fileIn.close

  let name = fileIn.readLine
  var sequence = Sequence()
  var seqRecord = SequenceRecord(name: name.strip()[1..^1], record: sequence)
  for line in fileIn.lines:
    if line[0] == '>':
      seqRecord.record.class = guess(seqRecord.record.chain)
      yield seqRecord
      sequence = Sequence()
      seqRecord = SequenceRecord(name: line.strip()[1..^1], record: sequence)
    else:
      seqRecord.record.chain.add line.strip().toUpper()
  yield seqRecord
