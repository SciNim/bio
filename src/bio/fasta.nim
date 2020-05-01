import ../bio

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

proc write*(record: SequenceRecord, fName, kind: string="fasta") =
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
