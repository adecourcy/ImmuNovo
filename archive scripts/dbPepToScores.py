from sys import argv
from fileParser import *
from processSpectrum import *
from processResults import *
from main import *
from decimal import Decimal
import random
import gc

def massToleranceMaxDiff(calculatedMass):

  minExp = int(Decimal(calculatedMass) / (Decimal("0.000035") + Decimal(1)))
  maxExp = int(Decimal(calculatedMass) / (Decimal(1) - Decimal("0.000035")))

  return (minExp, maxExp)


def getDecoys(mass, pepFile):
  possibles = []

  minDiff, maxDiff = massToleranceMaxDiff(mass)

  for mass in range(minDiff, maxDiff+1):
    if mass in pepFile:
      possibles += pepFile[mass]

  if len(possibles) <= 10:
    return possibles

  else:
    return random.sample(possibles, 10)

def calculatedAminoScore(peptide, aminoScoreList):
  aminoScores = 0
  aminoDict = aminoScoreList[len(peptide) - 9]
  for position, amino in zip(range(len(peptide)), peptide):
    aminoScores += aminoDict[amino][position]
  return (float(aminoScores) / len(peptide))


def fileToDict(pepFile):
  pepDict = {}
  allValid = True

  total = 0

  for pepLine in pepFile:
    if pepLine == "":
      break
    length, peptide, pepMass = pepLine.split()
    if len(peptide) < 9 or len(peptide) > 12:
      continue
    total += 1
    for amino in peptide:
      if amino not in ["G", "A", "S", "P", "V",
                       "T", "L", "I", "N", "D",
                       "Q", "K", "E", "M", "H",
                       "F", "R", "C", "Y", "W", "B"]:
        allValid = False
        break
    if allValid == False:
      allValid = True
      continue
    pepMass = int(float(pepMass) * 10**4)
    if pepMass in pepDict:
      pepDict[pepMass].append(peptide)
    else:
      pepDict[pepMass] = [peptide]

  return pepDict

if __name__ == "__main__":
  if len(argv) < 5:
    print("peptideFile spectrumFile aminoMassFile aminoScoreFile")
    exit()

  with open(argv[1]) as f:
    pepFile = f.read()
  pepFile = pepFile.split("\n")

  pepFile = fileToDict(pepFile)
  gc.collect()
  print("got here")

  parsedMassFile, aminoMasses = parseFiles(argv[2], argv[3])

  output = open("decoyPeps.csv", 'w')

  print("got here")

  iterations = 0

  for trial in parsedMassFile:
    iterations += 1
    pepMass, charge, spectrumMasses, spectrumScores, title = trial

    aminoScoresList = parseAminoAcidScores(argv[4])
    aminoScoresList2 = parseAminoAcidScores(argv[4])

    spectrumMasses, spectrumMassesDouble, \
    spectrumScores, spectrumScoresDouble = \
                              processSpectrum(spectrumMasses,
                                              spectrumScores,
                                              pepMass,
                                              H2OMASS,
                                              PROTONMASS,
                                              charge,
                                              2)

    aminoMassesAdjusted, aminoScoresList2, \
    spectrumMasses, spectrumMassesDouble, \
    spectrumScores, spectrumScoresDouble, \
    protonMassAdjusted, H2OMassAdjusted = \
              adjustAllForPrecision(aminoMasses,
                                    aminoScoresList2,
                                    spectrumMasses,
                                    spectrumMassesDouble,
                                    spectrumScores,
                                    spectrumScoresDouble,
                                    PROTONMASS,
                                    H2OMASS,
                                    4)

    experimentalSpectrum, experimentalScores = \
          mergeMassesScores(spectrumMasses, spectrumMassesDouble,
                            spectrumScores, spectrumScoresDouble)

    peptides = getDecoys(spectrumMasses[-1], pepFile)
    scoredPeptides = []
    globalScore = 0
    aminoScore = 0
    
    for peptide in peptides:
      globalScore = calculateGlobalScore(aminoMassesAdjusted,
                                         experimentalSpectrum,
                                         experimentalScores,
                                         0.5,
                                         peptide,
                                         protonMassAdjusted,
                                         H2OMassAdjusted,
                                         int((NH3MASS * (10 ** 4))),
                                         35)

      aminoScore = (calculatedAminoScore(peptide[::-1], aminoScoresList))
      scoredPeptides.append((peptide, (globalScore * aminoScore)))

    if scoredPeptides == []:
      output.write(title + "," + "No Peptide Found\n")
    else:
      scoredPeptides.sort(key=lambda x: x[1], reverse=True)
      output.write(title + "," + str(scoredPeptides[0][0]) + "," + str(scoredPeptides[0][1]) + "\n")
