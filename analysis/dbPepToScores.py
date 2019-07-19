import sys
import os
# cheap hack until I can figure out how to do this properly in
# python 3.6
sys.path.append(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])

import random
from decimal import Decimal

import backend.Structures.pssmIO as PSSMIO
import backend.Structures.acidMassTableIO as AcidMassTableIO
import backend.Structures.spectrumIO as SpectrumIO
import backend.Structures.pssm as PSSM
import backend.Structures.acidMassTable as AcidMassTable
import backend.Structures.spectrum as Spectrum
import backend.PreProcessing.acidConversion as AcidConversion
import backend.PreProcessing.spectrumConversion as SpectrumConversion
import backend.PostProcessing.processResults as ProcessResults
from backend.constants import *
from backend.userInput import *


def massToleranceMaxDiff(calculatedMass, massTolerance):

  minExp = int(Decimal(calculatedMass) / (Decimal(massTolerance) + Decimal(1)))
  maxExp = int(Decimal(calculatedMass) / (Decimal(1) - Decimal(massTolerance)))

  return (minExp, maxExp)


def getDecoys(mass, peptideDict, massTolerance, maxDecoys):
  possibles = []

  minDiff, maxDiff = massToleranceMaxDiff(mass, massTolerance)

  for mass in range(minDiff, maxDiff+1):
    if mass in peptideDict:
      possibles += peptideDict[mass]

  if len(possibles) <= maxDecoys:
    return possibles

  else:
    return random.sample(possibles, maxDecoys)


def calculatedAminoScore(peptide, aminoScoreList):
  aminoScores = 0
  aminoDict = aminoScoreList[len(peptide) - 9]
  for position, amino in zip(range(len(peptide)), peptide):
    aminoScores += aminoDict[amino][position]
  return (float(aminoScores) / len(peptide))


def convertPeptide(peptide,
                   acidConversion,
                   conversionOrder):
  for conversion in conversionOrder:
    newPep = peptide.replace(conversion,
                             acidConversion[conversion])
  return newPep


def getPepMass(peptide,
               acidMassTable):
  mass = 0
  for acid in peptide:
    mass += AcidMassTable.getMass(acidMassTable, acid)
  return mass


def fileToDict(pepFile,
               precision,
               minP,
               maxP,
               acidConversion,
               conversionOrder,
               acidMassTable,
               reverse):

  pepDict = {}

  for pepLine in pepFile:
    if pepLine == "":
      break
    peptide = pepLine.strip()
    if 'U' in peptide:
      continue
    if len(peptide) < minP or len(peptide) > maxP:
      continue
    peptide = convertPeptide(peptide,
                             acidConversion,
                             conversionOrder)
    if reverse == True:
      peptide = peptide[::-1]
    if peptide == '':
      continue
    pepMass = getPepMass(peptide, acidMassTable)
    if pepMass in pepDict:
      pepDict[pepMass].append(peptide)
    else:
      pepDict[pepMass] = [peptide]

  return pepDict


def getPSSMScore(peptide,
                 pssmTitle,
                 allPSSM):

  if len(peptide) in PSSM.getIgnoredLengths(pssmTitle, allPSSM):
    return 0

  score = 0
  matrix = PSSM.getMatrixOfLength(pssmTitle, len(peptide), allPSSM)
  for position, acid in enumerate(peptide):
    score += matrix[acid][position]

  return score


def programUsageOutput():
  print("\nUsage: dbPepToScores.py SPECDIR ACIDMASSFILE PSSMDIR DCYDBDIR REVERSE\n\n")

  padding = '{:<20}{}'
  print("Positional arguments:")
  print(str.format(padding, 'SPECDIR',
        'A directory containing 1 or more spectrum files'))
  print(str.format(padding, 'ACIDMASSFILE',
        'A file containing mass spectrometry data'))
  print(str.format(padding, 'PSSMDIR',
        'A directory containing a positional scoring matrix'))
  print(str.format(padding, 'DCYDBDIR',
        'A directory containing a databases of decoy peptides'))
  print(str.format(padding, 'REVERSE',
        'Should decoy peptides be reversed?'))
  print('{:<20}{}'.format('','"t" or "true" or "f" or "false" (not case sensitive)'))
  print(str.format(padding, 'MAX DECOYS',
        'What is the maximum number of decoy peptides we consider for a spectrum?'))

  printOptionalArguments()
  exit()


def checkInputLength(args):
  if len(args) == 1:
    programUsageOutput()

  elif len(args) < 7:
    print("\nThis script takes at least 7 arguments.\n")
    print("Call this script with no arguments for usage details\n")
    print("Exiting...\n")



if __name__ == '__main__':

  fileLocation = 'analysis'
  fileName = 'decoyScores.csv'
  checkInputLength(sys.argv)

  spectrumDirectory, acidMassFile, \
    pssmDirectory, decoyPeptideDirectory, \
    reverse, maxDecoys = sys.argv[1:7]
  maxDecoys = int(maxDecoys)
  defaultParameters = parseParameterInput(sys.argv[7:])

  MMTString = str(defaultParameters['MMT'])
  massTolerance = '0.' + ('0' * (6 - len(MMTString))) + MMTString

  reverse = reverse.lower()
  if reverse == 'true' or reverse == 't':
    reverse = True
  elif reverse == 'false' or reverse == 'f':
    reverse = False
  else:
    print("unknown reverse parameter '{}'".format(reverse))

  decoyNames = os.listdir(decoyPeptideDirectory)
  decoyFiles = \
      [os.path.join(decoyPeptideDirectory, name) for name in decoyNames]
  tryFiles(decoyFiles)

  acidMassTable = \
      AcidMassTable.adjustForPrecision(
          AcidMassTableIO.getAminoMasses(acidMassFile),
          defaultParameters['PREC'])

  allPSSM = \
      PSSM.adjustForPrecision(
          PSSMIO.getAllPSSM(pssmDirectory,
                            defaultParameters['minP'],
                            defaultParameters['maxP']),
          defaultParameters['PREC'])

  conversionTable = \
    AcidConversion.createAcidConversionTable([acid for acid in acidMassTable])
  conversionOrder = []
  for conversion in conversionTable:
    conversionOrder.append(conversion)
  conversionOrder.sort(key=lambda x: len(x), reverse=True)

  if conversionTable == {}:
    print("Currently this program only accepts up to 23 amino acids total\n")
    print("The mass file you have provided "
          "has more than 23 amino acids defined\n")
    print("Other files may also have too "
          "many masses but haven't been checked\n"
          "by the program at this time\n")
    exit()

  # quick and dirty way to make sure all our amino masses match for now
  sanityCheck(allPSSM, acidMassTable)

  acidMassTable, allPSSM = \
      AcidConversion.convertAcidModifications(acidMassTable,
                                              allPSSM,
                                              conversionTable)

  H2OMassAdjusted = int(H2OMASS * (10**defaultParameters['PREC']))
  NH3MassAdjusted = int(NH3MASS * (10**defaultParameters['PREC']))
  protonMassAdjusted = int(PROTONMASS * (10**defaultParameters['PREC']))

  allDecoyPeptides = {}

  for decoyFile in decoyFiles:
    with open(decoyFile, 'r') as f:
      decoys = f.read()
    decoys = decoys.split('\n')

    decoys = fileToDict(decoys,
                        defaultParameters['PREC'],
                        defaultParameters['minP'],
                        defaultParameters['maxP'],
                        conversionTable,
                        conversionOrder,
                        acidMassTable,
                        reverse)
    

    for spectrumFileName in os.listdir(spectrumDirectory):

      if spectrumFileName not in allDecoyPeptides:
        allDecoyPeptides[spectrumFileName] = {}

      spectrumFile = os.path.join(spectrumDirectory, spectrumFileName)

      for spectrum in SpectrumIO.getSpectrums(spectrumFile):
        spectrumTitle = Spectrum.getTitle(spectrum)

        if spectrumTitle not in allDecoyPeptides[spectrumFileName]:
          allDecoyPeptides[spectrumFileName][spectrumTitle] = []

        spectrumMasses, spectrumMassesDouble, \
        spectrumIntensities, spectrumIntensitiesDouble = \
                                                         \
          SpectrumConversion.adjustSpectrumPrecision(
              *SpectrumConversion.processSpectrum(
                            Spectrum.getMasses(spectrum),
                            Spectrum.getIntensities(spectrum),
                            Spectrum.getPrecursorMass(spectrum),
                            H2OMASS,
                            PROTONMASS,
                            Spectrum.getCharge(spectrum),
                            defaultParameters['COMP']),
              defaultParameters['PREC'])

        peptides = \
            getDecoys(spectrumMasses[-1], decoys, massTolerance, maxDecoys)

        for peptide in peptides:

          experimentalSpectrum = \
              spectrumMasses + spectrumMassesDouble
          experimentalIntensities = \
              spectrumIntensities + spectrumIntensitiesDouble

          globalScore = \
              ProcessResults.calculateGlobalScore(acidMassTable,
                                                  experimentalSpectrum,
                                                  experimentalIntensities,
                                                  peptide,
                                                  protonMassAdjusted,
                                                  H2OMassAdjusted,
                                                  NH3MassAdjusted,
                                                  defaultParameters['MMT'])

          allDecoyPeptides[spectrumFileName][spectrumTitle].append((peptide,
                                                                   globalScore))


  with open(os.path.join(fileLocation, fileName), 'w') as f:

    f.write('{},{},{},{},{},{},{},{},{}\n'.format(PEPTIDE_GLOBAL_BEST, SCORE_GLOBAL,
                                  PEPTIDE_PSSM_BEST, SCORE_PSSM,
                                  PEPTIDE_COMBINED_BEST, SCORE_COMBINED,
                                  TITLE_PSSM, TITLE_SPECTRUM,
                                  FILE_SPECTRUM))

    for pssmTitle in allPSSM:
      ignoredLengths = PSSM.getIgnoredLengths(pssmTitle, allPSSM)
      bestAminos = []
      bestCombined = []
      bestGlobal = []
      spectrums = []

      for spectrumFileName in allDecoyPeptides:

        for spectrumTitle in allDecoyPeptides[spectrumFileName]:
          addedScores = []

          if allDecoyPeptides[spectrumFileName][spectrumTitle] == []:
            continue

          for peptide in allDecoyPeptides[spectrumFileName][spectrumTitle]:
            positionalScore = \
              getPSSMScore(peptide[0], pssmTitle, allPSSM) / \
              (10 ** defaultParameters['PREC'] * len(peptide[0]))

            addedScores.append((peptide[0],
                                peptide[1],
                                positionalScore,
                                positionalScore * peptide[1]))

          addedScores.sort(key=lambda x: x[1], reverse=True)
          bestGlobal.append((addedScores[0][0],
                              addedScores[0][1]))
          addedScores.sort(key=lambda x: x[2], reverse=True)
          bestAminos.append((addedScores[0][0],
                              addedScores[0][2]))
          addedScores.sort(key=lambda x: x[3], reverse=True)
          bestCombined.append((addedScores[0][0],
                                addedScores[0][3]))
          spectrums.append((spectrumTitle, spectrumFileName))
    
      for gBest, aBest, cBest, spec in zip(bestGlobal,
                                           bestAminos,
                                           bestCombined,
                                           spectrums):
        f.write('{},{},{},{},{},{},{},{},{}\n'.format(gBest[0], gBest[1],
                                                      aBest[0], aBest[1],
                                                      cBest[0], cBest[1],
                                                      pssmTitle, spec[0],
                                                      spec[1]))