# importar arquivo fasta DNA
with open('DNA.fasta', 'r') as DNA:
  DNAseq = []
  for line in DNA:
    if line.startswith('>'):
      continue
    else:
      DNAseq.append(line)

# importar arquivo fasta introns
with open('introns.fasta', 'r') as introns:
  intronsSeq = []
  for line in introns:
    if line.startswith('>'):
      continue
    else:
      line = line.replace('\n','')
      intronsSeq.append(line)

# removendo introns do DNA
tempList = []
i = 0
while i < len(intronsSeq):
    if len(tempList) == 0:
        DNAprocessed = DNAseq[0].replace(intronsSeq[i], '')
        tempList.append(DNAprocessed)
        i += 1
    else:
        DNAprocessed = tempList[len(tempList)-1].replace(intronsSeq[i], '')
        i += 1

# transcription
mRNA = DNAprocessed.replace('T', 'U')

# translation
RNAcodonTable = {
    'UUU':'F', 'UUC':'F', 'UUA':'L', 'UUG':'L',
    'UCU':'S', 'UCC':'S', 'UCA':'S', 'UCG':'S',
    'UAU':'Y', 'UAC':'Y', 'UAA':'STOP', 'UAG':'STOP',
    'UGU':'C', 'UGC':'C', 'UGA':'STOP', 'UGG':'W',
    'CUU':'L', 'CUC':'L', 'CUA':'L', 'CUG':'L',
    'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'CAU':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
    'CGU':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
    'AUU':'I', 'AUC':'I', 'AUA':'I', 'AUG':'M',
    'ACU':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'AAU':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
    'AGU':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
    'GUU':'V', 'GUC':'V', 'GUA':'V', 'GUG':'V',
    'GCU':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'GAU':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
    'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
    }

protein = ''
for codonPosition in range(0, len(mRNA), 3):
    codon = mRNA[codonPosition: codonPosition + 3]
    if RNAcodonTable[codon] == 'STOP':
      break
    protein += RNAcodonTable[codon]

print(protein)