from sys import argv
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
from decimal import Decimal
from math import log


def getScore(splitLine):
  globalScore = float(splitLine[5])
  aminoScore = (float(splitLine[9]) / getPepLength(splitLine[2]))
  score = globalScore * aminoScore
  return score


def batchSpectrums(file):
  newBatch = {}
  for line in file:
    splitLine = line.split(",")
    if splitLine[1] not in newBatch:
      newBatch[splitLine[1]] = []
    if "No Peptide Found" in line:
      continue
    newBatch[splitLine[1]].append((splitLine[2], getScore(splitLine)))

  return newBatch

def getPepLength(peptide):
  return len(peptide.replace("M+15.995", "M"))


def findLowestFDR(trueScores, falseScores, cutoffFDR):
  trueScores.sort(reverse=True)
  falseScores.sort(reverse=True)
  cutoffList = []

  FDRList = []

  for numTrue, trueValue in zip(range(1, len(trueScores)+1), trueScores):
    for numFalse, falseValue in zip(range(len(falseScores)), falseScores):
      if falseValue < trueValue:
        FDRList.append((numFalse, numTrue, trueValue, numFalse/numTrue))
        if (numFalse / numTrue) < cutoffFDR:
          cutoffList.append((numTrue, numFalse, trueValue, (numFalse / numTrue)))
        break

  FDRList.sort(key=lambda x: x[3])
  
  for elm in FDRList:
    print("Num False: " + str(elm[0]))
    print("Num True: " + str(elm[1]))
    print("True Value: " + str(elm[2]))
    print("FDR: " + str(elm[3]))
    print()

  cutoffList.sort(key=lambda x: x[0], reverse=True)
  for val in cutoffList:
    print(str(val[0]) + ", " + str(val[1]) + ", " + str(val[2]) + ", " + str(val[3]))

def plotCombination(combinedTrue, combinedFalse, cutoffFDR):

  print("True total: " + str(len(combinedTrue)))
  print("False total: " + str(len(combinedFalse)))
  print()

  #findLowestFDR(combinedTrue, combinedFalse, cutoffFDR)

  combinedFinal = []
  for x in combinedTrue:
    if x < 0.0155:
      combinedFinal.append(x)

  plt.hist(combinedFinal, 
           bins=50,
           histtype='step',
           color="blue",
           label="DeNovo Output")

  plt.hist(combinedFalse,
           bins=50,
           histtype='step',
           color="red",
           label="Random Output")

  blue_patch = mpatches.Patch(color="blue", label="De Novo Score Frequency")
  red_patch = mpatches.Patch(color="red", label="Decoy Score Frequency")
  plt.legend(handles=[blue_patch, red_patch])
  plt.xlabel("Global Score")
  plt.ylabel("Frequency")
  #plt.title("Global Score Plot")
  #plt.show()
  plt.savefig("./Results/decoyRAndRComp.png", dpi=300)

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

  originalBatch = batchSpectrums(original)

  uniquePeptidesOriginal = {}
  uniquePeptidesDecoy = {}
  originalScores = []
  decoyScores = []


  for spectrum in originalBatch:
    if originalBatch[spectrum] == []:
      continue
    originalBatch[spectrum].sort(key= lambda x: x[1], reverse=True)
    if originalBatch[spectrum][0][0] not in uniquePeptidesOriginal:
      uniquePeptidesOriginal[originalBatch[spectrum][0][0]] = 1
    else:
      uniquePeptidesOriginal[originalBatch[spectrum][0][0]] += 1
    originalScores.append(originalBatch[spectrum][0][1])


  for spectrum in false:

    if "No Peptide Found" in spectrum:
      continue

    spectrumList = spectrum.split(",")
    if spectrumList[1] not in uniquePeptidesDecoy:
      uniquePeptidesDecoy[spectrumList[1]] = 1
    else:
      uniquePeptidesDecoy[spectrumList[1]] += 1
    decoyScores.append(float(spectrumList[2]))

  print("Number of unique spectrums: " + str(len(originalBatch)))
  print("Number of unique DeNovo Peptides: " \
          + str(len(uniquePeptidesOriginal)))
  print("Number of unique Decoy Peptides: " + str(len(uniquePeptidesDecoy)))

  plotCombination(originalScores, decoyScores, 0.01)


if __name__ == "__main__":
  analyzeResults(argv[1], argv[2])