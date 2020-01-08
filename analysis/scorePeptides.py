#!/usr/bin/env python3
# Give a set of peptides, PSSMs, and associated Spectra, score peptides

# cheap hack until I can figure out how to do this properly in
# python 3.6
import sys, os
sys.path.append(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])

import argparse
import pandas as pd
import numpy as np

from copy import deepcopy

import backend.userInput as userInput
import backend.Structures.spectrum as Spectrum
import backend.Structures.spectrumIO as SpectrumIO
import backend.PreProcessing.spectrumConversion as SpectrumConversion
import backend.Structures.pssm as PSSM
import backend.PostProcessing.processResults as processResults
import backend.PostProcessing.processResultsOld as processResultsOld
import backend.PreProcessing.acidConversion as AcidConversion
import backend.PreProcessing.misc as Misc
from backend.constants import *


def getUniquePeptides(peptideList):
  uniquePeptides = {}
  for pep in peptideList:
    uniquePeptides[pep] = 0
  return uniquePeptides

def scorePeptides(uniquePeptideScoreDict, conversionDict, scoreFunc):
  for peptide in uniquePeptideScoreDict:
    uniquePeptideScoreDict[peptide] = scoreFunc(conversionDict[peptide])
  return uniquePeptideScoreDict

def pssmScore(peptideString, allPSSM, title):
  # Assumes converted Peptide String
  matrix = PSSM.getMatrixOfLength(title, len(peptideString), allPSSM)
  score = 0
  for i in range(len(peptideString)):
    score += PSSM.getAcidProbabilities(matrix, peptideString[i])[i]
  score /= len(peptideString)
  return score

def pssmScoreForPeptides(uniquePeptideScoreDict, conversionDict, allPSSM, title):
  scoringFunction = \
      lambda peptideString: pssmScore(peptideString, allPSSM, title)
  return scorePeptides(uniquePeptideScoreDict, conversionDict, scoringFunction)

def allPeptideAllPssmScores(uniquePeptideScoreDict, conversionDict, allPSSM):
  pssmScores = []

  for title in allPSSM:
    pssmScores.append((title, pssmScoreForPeptides(deepcopy(uniquePeptideScoreDict),
                                                   conversionDict,
                                                   allPSSM,
                                                   title)))
  
  return pssmScores

def getMaximalPssmScoreForPeptide(pssmScores, peptide):
  # Output from allPeptideAllPssmScores
  maxScore = -1
  for entry in pssmScores:
    currentScore = entry[1][peptide]
    if currentScore > maxScore:
      maxScore = currentScore
      maxTitle = entry[0]
  return (maxTitle, maxScore)

def getMaximalPssmScores(uniquePeptideScoreDict, conversionDict, allPSSM):

  pssmPeptideScoreDict = deepcopy(uniquePeptideScoreDict)

  pssmScores = \
      allPeptideAllPssmScores(uniquePeptideScoreDict, conversionDict, allPSSM)
  
  for peptide in pssmPeptideScoreDict:
    pssmPeptideScoreDict[peptide] = \
        getMaximalPssmScoreForPeptide(pssmScores, peptide)
  
  return pssmPeptideScoreDict

def spectrumScore(acidMassTable,
                  experimentalSpectrum,
                  experimentalScores,
                  protonMassModified,
                  H2OMassModified,
                  NH3MassModified,
                  maxMassTolerance,
                  peptide):


  globalScore  = processResults.calculateGlobalScore(acidMassTable,
                                                     experimentalSpectrum,
                                                     experimentalScores,
                                                     peptide,
                                                     protonMassModified,
                                                     H2OMassModified,
                                                     NH3MassModified,
                                                     maxMassTolerance)

  return globalScore

def spectrumScoreOld(acidMassTable,
                     experimentalSpectrum,
                     experimentalScores,
                     protonMassModified,
                     H2OMassModified,
                     NH3MassModified,
                     maxMassTolerance,
                     peptide):


  globalScore = processResultsOld.calculateGlobalScore(acidMassTable,
                                                       experimentalSpectrum,
                                                       experimentalScores,
                                                       0.5,
                                                       peptide,
                                                       protonMassModified,
                                                       H2OMassModified,
                                                       NH3MassModified,
                                                       maxMassTolerance)

  return globalScore

def allPeptideSpectrumScores(acidMassTable,
                             experimentalSpectrum,
                             experimentalScores,
                             protonMassModified,
                             H2OMassModified,
                             NH3MassModified,
                             maxMassTolerance,
                             uniquePeptideScoreDict,
                             conversionDict):

  scoreFunction = lambda peptide: spectrumScore(acidMassTable,
                                                experimentalSpectrum,
                                                experimentalScores,
                                                protonMassModified,
                                                H2OMassModified,
                                                NH3MassModified,
                                                maxMassTolerance,
                                                peptide)
  
  return scorePeptides(deepcopy(uniquePeptideScoreDict), conversionDict, scoreFunction)

def allPeptideSpectrumScoresOld(acidMassTable,
                                experimentalSpectrum,
                                experimentalScores,
                                protonMassModified,
                                H2OMassModified,
                                NH3MassModified,
                                maxMassTolerance,
                                uniquePeptideScoreDict,
                                conversionDict):

  scoreFunction = lambda peptide: spectrumScoreOld(acidMassTable,
                                                   experimentalSpectrum,
                                                   experimentalScores,
                                                   protonMassModified,
                                                   H2OMassModified,
                                                   NH3MassModified,
                                                   maxMassTolerance,
                                                   peptide)
  
  return scorePeptides(deepcopy(uniquePeptideScoreDict), conversionDict, scoreFunction)

def spectrumVariableSetup(acidMassFile,
                          pssmDir,
                          minPepLength,
                          maxPepLength,
                          precision):

  allPSSM, acidMassTable, conversionTable = \
                              Misc.getAminoVariables(acidMassFile,
                                                     precision,
                                                     pssmDir,
                                                     minPepLength,
                                                     maxPepLength)
  H2OMassAdjusted = int(H2OMASS * (10**precision))
  NH3MassAdjusted = int(NH3MASS * (10**precision))
  protonMassAdjusted = int(PROTONMASS * (10**precision))

  return allPSSM, acidMassFile, conversionTable, \
         H2OMASS, NH3MASS, protonMassAdjusted

def peptideSpectrumDict(peptideDF):
  specToPeptideDict = {}
  spectrumList = list(peptideDF[TITLE_SPECTRUM])
  peptideList = list(peptideDF[PEPTIDE])

  for spectrum, peptide in zip(spectrumList, peptideList):
    if spectrum not in specToPeptideDict:
      specToPeptideDict[spectrum] = []
    specToPeptideDict[spectrum].append(peptide)
  
  for spectrum in specToPeptideDict:
    specToPeptideDict[spectrum] = getUniquePeptides(specToPeptideDict[spectrum])
  
  return specToPeptideDict

###############################################################################

def scorePeptideSpectrumOverlap(peptideDF,
                                 conversionDict,
                                 acidMassTable,
                                 conversionTable,
                                 H2OMassAdjusted,
                                 NH3MassAdjusted,
                                 protonMassAdjusted,
                                 spectrumDirectory,
                                 scoringFunction,
                                 precision=4,
                                 minPepLength=9,
                                 maxPepLength=12,
                                 maxMassTolerance=35,
                                 compression=2):

  specToPeptideDict = peptideSpectrumDict(peptideDF)

  spectrumTitles = set()

  for spectrumFileName in os.listdir(spectrumDirectory):
    spectrumFile = os.path.join(spectrumDirectory, spectrumFileName)
    for spectrum in SpectrumIO.getSpectrums(spectrumFile):
      spectrumTitle = Spectrum.getTitle(spectrum)

      if spectrumTitle in spectrumTitles:
        continue # some of these files seem to have duplicates
      spectrumTitles.add(spectrumTitle)
      if spectrumTitle not in specToPeptideDict:
        continue
      spectrumPeptides = specToPeptideDict[spectrumTitle]

      spectrumMasses, spectrumMassesDouble, \
      spectrumIntensities, spectrumIntensitiesDouble = \
          SpectrumConversion.adjustSpectrumPrecision(
              *SpectrumConversion.processSpectrum(Spectrum.getMasses(spectrum),
                                                  Spectrum.getIntensities(spectrum),
                                                  Spectrum.getPrecursorMass(spectrum),
                                                  H2OMASS,
                                                  PROTONMASS,
                                                  Spectrum.getCharge(spectrum),
                                                  compression),
              precision)

      experimentalSpectrum = spectrumMasses + spectrumMassesDouble
      experimentalIntensities = spectrumIntensities + spectrumIntensitiesDouble

      # Scoring function passed in, should be one of ours
      specToPeptideDict[spectrumTitle] = \
          scoringFunction(acidMassTable,
                          experimentalSpectrum,
                          experimentalIntensities,
                          protonMassAdjusted,
                          H2OMassAdjusted,
                          NH3MassAdjusted,
                          maxMassTolerance,
                          specToPeptideDict[spectrumTitle],
                          conversionDict)
  
  peptideDF[SCORE_GLOBAL] = \
    peptideDF.apply(lambda x: specToPeptideDict[x[TITLE_SPECTRUM]][x[PEPTIDE]], axis=1)

  return peptideDF

def scorePeptidesPssm(peptideDF, conversionDict, allPSSM):
  uniquePeptideDict = \
      getUniquePeptides(list(peptideDF[PEPTIDE]))
  scoredPeptideDict = \
      getMaximalPssmScores(uniquePeptideDict, conversionDict, allPSSM)

  peptideDF[TITLE_PSSM] = \
    peptideDF.apply(lambda x: scoredPeptideDict[x[PEPTIDE]][0], axis=1)
  peptideDF[SCORE_PSSM] = \
    peptideDF.apply(lambda x: scoredPeptideDict[x[PEPTIDE]][1], axis=1)
  
  return peptideDF


def getPeptideScores(peptideDF,
                     peptideConversionDict,
                     acidMassTable,
                     conversionTable,
                     H2OMassAdjusted,
                     NH3MassAdjusted,
                     protonMassAdjusted,
                     spectrumDirectory,
                     scoringFunction,
                     allPSSM,
                     precision=4,
                     minPepLength=9,
                     maxPepLength=12,
                     maxMassTolerance=35,
                     compression=2):

  if scoringFunction == CALCULATION_NN:
    pass
  elif scoringFunction == CALCULATION_OVERLAP_NEW:
    peptideDF = scorePeptideSpectrumOverlap(peptideDF,
                                            peptideConversionDict,
                                            acidMassTable,
                                            conversionTable,
                                            H2OMassAdjusted,
                                            NH3MassAdjusted,
                                            protonMassAdjusted,
                                            spectrumDirectory,
                                            allPeptideSpectrumScores,
                                            precision,
                                            minPepLength,
                                            maxPepLength,
                                            maxMassTolerance,
                                            compression)

  elif scoringFunction == CALCULATION_OVERLAP_OLD:
    peptideDF = scorePeptideSpectrumOverlap(peptideDF,
                                            peptideConversionDict,
                                            acidMassTable,
                                            conversionTable,
                                            H2OMassAdjusted,
                                            NH3MassAdjusted,
                                            protonMassAdjusted,
                                            spectrumDirectory,
                                            allPeptideSpectrumScoresOld,
                                            precision,
                                            minPepLength,
                                            maxPepLength,
                                            maxMassTolerance,
                                            compression)

  peptideDF = scorePeptidesPssm(peptideDF, peptideConversionDict, allPSSM)

  peptideDF[SCORE_COMBINED] = peptideDF[SCORE_GLOBAL] * peptideDF[SCORE_PSSM]
  
  return peptideDF
