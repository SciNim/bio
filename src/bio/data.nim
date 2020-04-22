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
   'D': 'H', 'B': 'V', 'X': 'X', 'N': 'N', '-': '-'}.toTable

const codonTable*: Table[string, char] =
  {"TTT": 'F', "TTC": 'F', "TTA": 'L', "TTG": 'L',
   "TCT": 'S', "TCC": 'S', "TCA": 'S', "TCG": 'S',
   "TAT": 'Y', "TAC": 'Y', "TGT": 'C', "TGC": 'C',
   "TGG": 'W', "CTT": 'L', "CTC": 'L', "CTA": 'L',
   "CTG": 'L', "CCT": 'P', "CCC": 'P', "CCA": 'P',
   "CCG": 'P', "CAT": 'H', "CAC": 'H', "CAA": 'Q',
   "CAG": 'Q', "CGT": 'R', "CGC": 'R', "CGA": 'R',
   "CGG": 'R', "ATT": 'I', "ATC": 'I', "ATA": 'I',
   "ATG": 'M', "ACT": 'T', "ACC": 'T', "ACA": 'T',
   "ACG": 'T', "AAT": 'N', "AAC": 'N', "AAA": 'K',
   "AAG": 'K', "AGT": 'S', "AGC": 'S', "AGA": 'R',
   "AGG": 'R', "GTT": 'V', "GTC": 'V', "GTA": 'V',
   "GTG": 'V', "GCT": 'A', "GCC": 'A', "GCA": 'A',
   "GCG": 'A', "GAT": 'D', "GAC": 'D', "GAA": 'E',
   "GAG": 'E', "GGT": 'G', "GGC": 'G', "GGA": 'G',
   "GGG": 'G', "---": '-', "TAA": '*', "TAG": '*', "TGA": '*'}.toTable

const rnaAmbiguousComplement*: Table[char, char] =
  {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A',
   'M': 'K', 'R': 'Y', 'W': 'W', 'S': 'S',
   'Y': 'R', 'K': 'M', 'V': 'B', 'H': 'D',
   'D': 'H', 'B': 'V', 'X': 'X', 'N': 'N', '-': '-'}.toTable
