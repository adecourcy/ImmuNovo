def getMass(peptide):
  masses = { 'G': 57.02147,
'A': 71.03712,
'S': 87.03203,
'P': 97.05277,
'V': 99.06842,
'T': 101.04768,
'L': 113.08407,
'I': 113.08407,
'N': 114.04293,
'D': 115.02695,
'Q': 128.05858,
'K': 128.09497,
'E': 129.04260,
'M': 131.04049,
'H': 137.05891,
'B': 147.03549,
'F': 147.06842,
'R': 156.10112,
'C': 160.03065,
'Y': 163.06333,
'W': 186.07932 }

  mass = 0

  for amino in peptide:
    #print(peptide)
    mass += masses[amino]
  return str(mass)

if __name__ == "__main__":

  with open('C501.newMatrix.modification.output.tsv') as f:
    resultsFile = f.read()

  resultsFile = resultsFile.split("\n")
  del resultsFile[0]
  if resultsFile[-1] == "":
    del resultsFile[-1]

  output = open("resultsAsDecoy.csv", 'w')

  for line in resultsFile:
    if "No Peptide Found" in line:
      continue
    splitLine = line.split(",")
    peptide = splitLine[2].replace("M+15.995", "B")[::-1]
    output.write(str(len(peptide)) + " " + peptide + " " + getMass(peptide) + "\n")
