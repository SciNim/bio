import tables

const dnaLettersAmbiguous: string = "GATCRYWSMKHBVDN"
const dnaLetters: string = "ACGT"

#   B == 5-bromouridine
#   D == 5,6-dihydrouridine
#   S == thiouridine
#   W == wyosine
const dnaLettersExtended: string = "GATCBDSW"

const dnaValuesAmbiguous: Table[char, string] =
  {'A': "A", 'C': "C", 'G': "G", 'T': "T",
   'M': "AC", 'R': "AG", 'W': "AT", 'S': "CG", 'Y': "CT", 'K': "GT",
   'V': "ACG", 'H': "ACT", 'D': "AGT", 'B': "CGT",
   'X': "GATC", 'N': "GATC"}.toTable

const dnaAmbiguousComplement*: Table[char, char] =
  {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
   'M': 'K', 'R': 'Y', 'W': 'W', 'S': 'S',
   'Y': 'R', 'K': 'M', 'V': 'B', 'H': 'D',
   'D': 'H', 'B': 'V', 'X': 'X', 'N': 'N'}.toTable
