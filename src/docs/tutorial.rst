=================
Bio: the tutorial
=================

:Author: Xabier Bello
:Version: |libversion|
:Nim Version: |nimversion|

.. contents::


Creating sequences
==================


.. code-block:: nim

    import bio/sequences

    let mySeq = Sequence(chain: "ACGTGTGAAC")


.. code-block:: nim

    import bio/sequences

    let myDNA = guess("ACGTGTGAAC")

    echo myDNA.complement

API: `sequences <sequences.html>`_.

Working with files
===================

TODO

* `fasta <fasta.html>`_ is a set of utilities to work with Fasta files.
