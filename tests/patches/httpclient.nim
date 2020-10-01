include "system/inclrtl"
import httpcore, strutils
export httpcore

type
  Response* = ref object
    status*: string
    headers*: HttpHeaders
    body: string

  HttpClient* = ref object
    urlRequest*: string

  HttpRequestError* = object of IOError

proc newHttpClient*: HttpClient =
  new result

proc body*(response: Response): string =
  return response.body

proc code*(response: Response): HttpCode =
  return response.status[0 .. 2].parseInt.HttpCode

proc get*(client: HttpClient, url: string): Response =
  new result
  result.headers = newHttpHeaders()
  client.urlRequest = url

  if "MN908947" in url:
    result.status = "200"
    result.body = ">MN908947.3\nATTAAAGGTTTATACCTTCCCAGGTAACAAACC\n"
  else:
    result.status = "404"
    result.body = ""
