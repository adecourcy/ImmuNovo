from main import parseInput
from typing import *
from random import shuffle
from decimal import Decimal
from fileParser import parseAminoAcidScores
import sys, os


def scramblePositionalMatrix(aminoScoresList: List[Dict[str, List[Decimal]]]) \
                                              -> List[Dict[str, List[Decimal]]]:

  aminos = []
  for key in aminoScoresList[0]:
    aminos.append(key)

  shuffle(aminos)

  newAminoScoresList = []

  for table in aminoScoresList:
    tempTable = []
    tempDict = {}
    for key in table:
      tempTable.append(table[key])
    #for scores in tempTable:
     # shuffle(scores)
    shuffle(tempTable)
    for key, index in zip(table, range(0, len(tempTable))):
      tempDict[key] = tempTable[index]
    newAminoScoresList.append(tempDict)

  return newAminoScoresList


def writeAminoScoresList(aminoScoresList, directory):
  if not os.path.exists(directory):
    os.makedirs(directory)
  falseScoresFile = open(directory + "falseScores", "w")

  for table in aminoScoresList:
    falseScoresFile.write("START\n")
    for key in table:
      scores = [str(x) for x in table[key]]
      falseScoresFile.write("{} {}\n".format(key, " ".join(scores)))
    falseScoresFile.write("END\n\n")

  falseScoresFile.close()


def getCurrentSpectrum(file, fileLine):
  if fileLine >= len(file):
    return ''
  splitLine = file[fileLine].split(",")
  return splitLine[1]


def getNextBatch(file, fileLine):

  currentSpectrum = getCurrentSpectrum(file, fileLine)

  results = []

  while(True):
    if fileLine >= len(file):
      break
    if currentSpectrum == getCurrentSpectrum(file, fileLine):
      results.append(file[fileLine])
      fileLine += 1
    else:
      break

  return results, fileLine

def getBestScoringPeptide(originalBatch, falseBatch):

  originalBatch = [x.split(",") for x in originalBatch]
  falseBatch = [x.split(",") for x in falseBatch]

  originalBatch.sort(key=lambda x: float(x[5]))
  originalBatch.reverse()
  falseBatch.sort(key=lambda x: float(x[5]))
  falseBatch.reverse()

  return originalBatch[0][2], falseBatch[0][2]


def comparePeptides(originalPeptide, falsePeptide, match, noMatch):
  originalPeptide = [x for x in originalPeptide]
  falsePeptide = [x for x in falsePeptide]
  originalPeptide.sort()
  falsePeptide.sort()

  if originalPeptide == falsePeptide:
    match += 1
  else:
    noMatch += 1

  return match, noMatch


def analyzeResults(originalOutput, falseOutput):
  with open(originalOutput) as f:
    original = f.read()
  with open(falseOutput) as f:
    false = f.read()

  original = original.split("\n")
  false = false.split("\n")
  del original[0]
  if original[-1] == '':
    del original[-1]
  del false[0]
  if false[-1] == '':
    del false[-1]

  originalFileLine = 0
  falseFileLine = 0

  match = 0
  noMatch = 0

  while(True):
    originalBatch, originalFileLine = getNextBatch(original, originalFileLine)
    falseBatch, falseFileLine = getNextBatch(false, falseFileLine)
    if originalBatch == []:
      break
    if "No Peptide Found" in originalBatch[0]:
      continue

    originalPeptide, falsePeptide = \
    	  getBestScoringPeptide(originalBatch, falseBatch)

    match, noMatch = \
  	  comparePeptides(originalPeptide, falsePeptide, match, noMatch)

  print("Of {} peptides, {}% (or {}) matched".format((match + noMatch),
                                                    int((match / (match + noMatch)) * 100),
                                                    match))

if __name__ == "__main__":
  defaultParameters, massFile, aaMasses, aaScoresDir = \
                                          parseInput(sys.argv)

  args = sys.argv[:4]
  args[0] = "main.py"

  if aaScoresDir[-1] != "/":
    aaScoresDir += "/"
  aminoScoresList = \
  	parseAminoAcidScores(aaScoresDir + os.listdir(aaScoresDir)[0])

  falseAminoScoresList = scramblePositionalMatrix(aminoScoresList)

  falseAminoScoresDirectory = "./falseTable/"
  writeAminoScoresList(falseAminoScoresList, falseAminoScoresDirectory)

  originalArgs = " ".join(args)
  originalOutput = defaultParameters["OUT"]
  for key in defaultParameters:
    originalArgs += " {}={}".format(key, defaultParameters[key])

  falseOutput = originalOutput+"falseOutput.csv"
  defaultParameters["OUT"] = falseOutput
  falseTableArgs = args[:4]
  falseTableArgs[3] = falseAminoScoresDirectory
  falseTableArgs = " ".join(falseTableArgs)
  for key in defaultParameters:
    falseTableArgs += " {}={}".format(key, defaultParameters[key])

  os.system("python3 {}".format(originalArgs))
  os.system("python3 {}".format(falseTableArgs))

  #analyzeResults(originalOutput, falseOutput)
