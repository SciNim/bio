import strutils
import ../bio


iterator sequences*(fName: string, kind: string="fasta"):
  SequenceRecord {.inline.} =
  ## Iterate through all the sequences in a given filename, yielding
  ## `SequenceRecords`.
  ##
  ## .. code-block::
  ##
  ##   for sequence in sequences("path/to/file.fas"):
  ##     doAssert(sequence of SequenceRecord)

  let fileIn: File = open(fName)
  defer: fileIn.close

  var sequence: Sequence = Sequence()
  var seqRecord: SequenceRecord
  for line in fileIn.lines:
    if (line[0] == '>' and sequence.chain.len > 0) or fileIn.endOfFile:
      sequence.class = guess(sequence.chain)
      yield seqRecord

    if line[0] == '>':
      sequence = Sequence()
      seqRecord = SequenceRecord(name: line.strip()[1..^1], record: sequence)
    else:
      sequence.chain.add line.strip().toUpper

proc load*(fName: string, kind: string="fasta"): seq[SequenceRecord] =
  ## Load a sequence of files from a filename.
  ##
  ## If you need an interator over a file, use the `sequences iterator<#sequences.i,string,string>`_
  ##
  ## .. code-block::
  ##
  ##  let mySeqs = load("path/to/file.fas")

  for sequence in sequences(fName):
    result.add sequence

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
    import ../bio

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
    import ../bio

    let mySeq = initDna("TGCACCCCA")
    let myRec = SequenceRecord(name: "My DNA sequence", record: mySeq)

    myRec.write("myOutput.fasta")

  let fHandler: File = open(fName, fmWrite)
  defer: fHandler.close()

  write(record, fHandler, kind)
