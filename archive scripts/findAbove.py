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
  return (score, float(globalScore), float(aminoScore))


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

def analyzeResults(originalOutput):
  with open(originalOutput) as f:
    original = f.read()

  original = original.split("\n")
  del original[0]
  if original[-1] == '':
    del original[-1]

  originalBatch = batchSpectrums(original)

  uniquePeptidesOriginal = {}
  uniquePeptidesDecoy = {}
  originalScores = []
  decoyScores = []

  aboveGlobal = 0
  abovePrior = 0
  aboveBoth = 0

  total = 0


  for spectrum in originalBatch:
    if originalBatch[spectrum] == []:
      continue
    originalBatch[spectrum].sort(key= lambda x: x[1][0], reverse=True)
    if originalBatch[spectrum][0][1][0] >= 0.0058:
      total += 1
      if originalBatch[spectrum][0][1][1] >= 0.1 and originalBatch[spectrum][0][1][2] >= 0.1:
        aboveBoth += 1
        aboveGlobal += 1
        abovePrior += 1
      elif originalBatch[spectrum][0][1][1] >= 0.1:
        aboveGlobal += 1
      elif originalBatch[spectrum][0][1][2] >= 0.1:
        abovePrior += 1

  print("Global above 0.1: " + str(aboveGlobal))
  print("Prior above 0.1: " + str(abovePrior))
  print("Both above 0.1: " + str(aboveBoth))
  print(total)

  sujunData = {}
  sujunFile = open("C501_3Cat.tsv", 'r')
  for line in sujunFile:
    line = line.split()
    if line[0] == "3":
      if line[1] not in sujunData:
        sujunData[line[1].strip()] = line[2].strip()
      else:
        print("Duplicate Spectrum: " + line[1])

  print()
  for spectrum in originalBatch:
    if originalBatch[spectrum] == []:
      continue
    if originalBatch[spectrum][0][1][0] < 0.0058:
      continue
    if spectrum not in sujunData:
      print("Sujun file missing spectrum: " + spectrum)

  print()
  for spectrum in sujunData:
    if spectrum not in originalBatch:
      print("Original file missing spectrum: " + spectrum)

  print()
  for spectrum in sujunData:
    if spectrum in originalBatch:
      if sujunData[spectrum] != originalBatch[spectrum][0][0]:
        print("For " + spectrum + " Sujun file: " + sujunData[spectrum] + " Original file: " + originalBatch[spectrum][0][0])


if __name__ == "__main__":
  analyzeResults("C501.newMatrix.modification.output.tsv")