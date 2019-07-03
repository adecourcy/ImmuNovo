import sys

def getTitles(spectrumFile):
  titles = []
  for line in spectrumFile:
    line = line.strip()
    if "TITLE" in line:
      title = line.split("=")[1].split(":")[0].split()[0]
      titles.append(title)
  return titles

def createNewTitle(peptideFile, appendTitle, currentPepLine):
  newLine = ""
  while True:
    if peptideFile[currentPepLine] == "\n":
      currentPepLine += 1
      return currentPepLine, newLine
    else:
      line = peptideFile[currentPepLine].strip()
      newLine += line.ljust(20) + appendTitle + "\n"
      currentPepLine += 1

def checkNewPeptide(peptideFile, spectrum, newFile):
  currentSpecLine = 0
  currentPepLine = 0
  currentLine = ""
  while currentPepLine < len(peptideFile):
    line = peptideFile[currentPepLine].strip()
    if "No Peptide Found" in line:
      print(currentSpecLine)
      currentPepLine += 1
      currentSpecLine += 1
    elif line != "":
      currentPepLine, pepLine = createNewTitle(peptideFile,
                                               spectrum[currentSpecLine],
                                               currentPepLine)
      print("\n--------")
      print(currentSpecLine)
      print("--------\n")
      currentSpecLine += 1
      newFile.write(pepLine)
    else:
      currentPepLine += 1
      

if __name__ == "__main__":

  try:
    peptideFile = open(sys.argv[1], "r")
  except FileNotFoundError:
    print("Invalid file path")
    exit()

  peptideLines = peptideFile.readlines()

  try:
    spectrumFile = open(sys.argv[2], "r")
  except FileNotFoundError:
    print("Invalid file path")
    exit()

  spectrum = getTitles(spectrumFile)

  outputFile = open("converted", "w")

  checkNewPeptide(peptideLines, spectrum, outputFile) 

