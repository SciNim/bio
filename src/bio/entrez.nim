## :Author: |author|
## :Version: |libversion|
##
## ..
##
##   Entrez documentation https://www.ncbi.nlm.nih.gov/books/NBK25497/
##
##   The Entrez Programming Utilities (E-utilities) are a set of nine
##   server-side programs that provide a stable interface into the Entrez query
##   and database system at the National Center for Biotechnology Information
##   (NCBI).
##
## This module provides a shallow interface with some of the Entrez tools.
##
## If you skip the Entrez documentation linked above, you should at least know
## that if your client does more than 3 request per second you will receive
## an error or get IP-banned. This module encourages you to get a
## `registered account <http://www.ncbi.nlm.nih.gov/account/>`_ that will allow
## your code to get up to 10 reqs/s.
##
## Once registered, you can view/create an
## `API key in settings<https://www.ncbi.nlm.nih.gov/account/settings/>`_, and
## use that key in your code.
##
## **Important**
##
## Any code that imports this module should get compiled with the `-d:ssl flag
## <https://nim-lang.org/docs/nimc.html#compiler-usage-compile-time-symbols>`_.
## It can be automated by creating a new file besides your code with the same
## filename but adding a `.cfg` extension, containing the line `--define:ssl`.
##
import httpclient
import strformat
import strutils
import uri

import fasta


type
  Entrez* = ref EntrezObj
    ## This object is used to store your Entrez API key and a httpClient_.
    ##
    ## .. _httpClient: https://nim-lang.org/docs/httpclient.html#HttpClient
    ##
  EntrezObj = object
    apiKey*: string
    client*: HttpClient

const entrezUrl: string = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

proc newEntrez*(apiKey: string): Entrez =
  ## Returns a new Entrez object to query the NCBI-Entrez database.
  ##
  ## The API key is a 36 character length string, that you should use to create
  ## this object.
  ##
  ## .. code-block::
  ##
  ##   import entrez
  ##
  ##   let entrezClient = newEntrez("thisismyapikey1337")
  ##
  let client = newHttpClient()
  result = Entrez(apiKey: apiKey,
                  client: client)

proc fetchId*(e: Entrez, db, id: string): SequenceRecord =
  ## Get the SequenceRecord for a given `id` in the given `db`.
  ##
  ## The values allowed to `db` can be checked
  ## `here<https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch>`_.
  ##
  ## .. code-block::
  ##
  ##   import entrez
  ##
  ##   let entrezClient = newEntrez("thisismyapikey1337")
  ##
  ##   let mycobacterium = entrezClient.getId("nuccore", "AP018036")
  ##
  let params = encodeQuery(
    {"db": db,
     "api_key": e.apiKey,
     "id": id,
     "rettype": "fasta"})
  let url = entrezUrl & "efetch.fcgi?" & params

  let resp = e.client.get(url)

  if resp.code == Http200:
    for sr in sequences(resp.body.splitLines):
      return sr
  else:
    raise newException(HttpRequestError, &"Error found when fetching {url}")
