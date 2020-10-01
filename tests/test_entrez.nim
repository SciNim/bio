import httpclient
import unittest

import bio / [entrez, sequences]


suite "Test getting fasta sequences from Entrez":
  test "Creating the Entrez object":
    let etrz = newEntrez(apiKey="dummy_api_key_string")
    check etrz.apiKey == "dummy_api_key_string"

  test "Get a single object through an Entrez object":
    let etrz = newEntrez(apiKey="dummy_api_key_string")

    let seq = etrz.fetchId(db="nuccore", id="MN908947")

    check etrz.client.urlRequest ==
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" &
      "db=nuccore&api_key=dummy_api_key_string&id=MN908947&rettype=fasta"

    check seq.name == "MN908947.3"
    check seq.record.chain == "ATTAAAGGTTTATACCTTCCCAGGTAACAAACC"
    check seq.record.class == scDna

  test "Get a 404 when the Id is not found":
    let etrz = newEntrez(apiKey="dummy_api_key_string")

    expect HttpRequestError:
      discard etrz.fetchId(db="nuccore", id="NONE")
