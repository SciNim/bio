=================
Bio: the cookbook
=================

:Author: Xabier Bello
:Version: |libversion|
:Nim Version: |nimversion|

.. contents::


Building, storing and retrieving an Index
=========================================

When dealing with Fasta files you might found your code retrieving random files
from the same file, maybe in chained code. It can be useful to have the Fasta
`indexed<fasta.html#Index>`_ indexed once and then just retrieve the index in
the subsequent accessess.

To create and **Store** an Index:

.. code-block::

  import marshal, streams
  import bio/fasta

  let idx: Index = newIndex("YourFasta.fst")
  let idxFile = newFileStream("YourFastaIndex.idx", fmWrite)

  idxFile.store(idx)
  idxFile.close

To **Load** the index. The loading of an index is quite fast. The Index include
the name of the source file, so if the file changes between loads it becomes
obsolete.

.. code-block::

  import marshal, streams
  import bio/fasta

  let idxFile = newFileStream("YourFastaIndex.idx", fmWrite)
  let idx = to[Index](idxFile.readAll)
  idxFile.close
