#!/usr/bin/env python3

import sys
import os
# cheap hack until I can figure out how to do this properly in
# python 3.6
sys.path.append(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])

import random
import pandas as pd

import backend.Structures.acidMassTableIO as AcidMassTableIO
import backend.Structures.spectrumIO as SpectrumIO
import backend.Structures.acidMassTable as AcidMassTable
import backend.Structures.spectrum as Spectrum
import backend.PreProcessing.acidConversion as AcidConversion
import backend.PreProcessing.spectrumConversion as SpectrumConversion
import backend.PostProcessing.processResults as ProcessResults
import backend.PostProcessing.processResultsOld as ProcessResultsOld
from backend.constants import *
from backend.userInput import *

LENGTH = 'LENGTH'
MASS = 'MASS'
DECOYS = 'DECOYS'
PEPTIDE_CONVERTED = 'PEPTIDE_CONVERTED'


# Read decoy file*
# convert peptides (can be skipped?) [in AcidConversion]*
# Filter by length*
# Filter by unknown peptides*
# Calculate mass*
# separate by mass (dictionary)*

# read spectrum file*
# read each spectrum*
# select all peptides within mass tolerance of spectrum*
# randomly select peptides within mass tolerance of spectrum*
# Add peptides + spectrum title to dataframe*
# deconvert peptides*
# return results


# def getAllDecoyPeptides(decoyPeptideDirectory):
#   # Expects files to include masses, peptide, and peptide length with headers
#   # as stated in the constants

#   # Input: decoy database directory
#   # Output: dataframe
#   allData = []
#   for fileName in os.listdir(decoyPeptideDirectory):
#     currentFile = os.path.join(decoyPeptideDirectory, fileName)
#     allData.append(pd.read_csv(currentFile))
#   return pd.concat(allData)

# def filterByLength(df, minLength, maxLength):
#   # Input: Dataframe with column of peptides
#   # Output: Dataframe with column of peptides
#   df = df[df.apply(lambda x: x[LENGTH] <= maxLength, axis=1)]
#   df = df[df.apply(lambda x: x[LENGTH] >= minLength, axis=1)]
#   return df

# def validPeptide(peptide, validAcids):
#   for acid in peptide:
#     if acid not in validAcids:
#       return False
#   return True

# def removeUnknownPeptides(df, acidMasstable):
#   validAcids = AcidMassTable.getAcids(acidMasstable)
#   return df[df.apply(lambda x: validPeptide(x[PEPTIDE], validAcids), axis=1)]

# def calculateMasses(df, acidMassTable):
#   # Input: Dataframe with column of peptides
#   # Output: Dataframe with column of peptides and masses
#   df[MASS] = df.apply(lambda x: AcidMassTable.peptideMass(x[PEPTIDE]), axis=1)
#   return df

# def separateByMass(df, precision):
#   # Input: Dataframe with masses, peptides
#   # Output: Dictionary with key masses (int, adjusted for precision), list peptides
#   dataDict = {}
#   df[MASS] = df.apply(lambda x: int(round(x[MASS] * (10**precision))), axis=1)
#   for elm in df.groupby(MASS):
#     dataDict[mass] = list(elm[1][PEPTIDE])
#   return dataDict

def extractSpectrumInformation(spectrumFileDirectory, precision):
  # Return a dataframe with spectrum title and mass columns
  allInfo = []
  def extractFileInformation(spectrumFile):
    precursorMasses = []
    titles = []
    for spectrum in SpectrumIO.getSpectrums(spectrumFile):
      charge = Spectrum.getCharge(spectrum)
      precursorMasses.append(int(round((10 ** precision) * \
          ((float(Spectrum.getPrecursorMass(spectrum) * charge) - (PROTONMASS * (charge - 1)))))))
      titles.append(Spectrum.getTitle(spectrum))
    return pd.DataFrame({TITLE_SPECTRUM: titles,
                         PRECURSOR_MASS: precursorMasses})
  
  for fileName in os.listdir(spectrumFileDirectory):
    allInfo.append(extractFileInformation(os.path.join(spectrumFileDirectory, fileName)))

  return pd.concat(allInfo)

def massToleranceMaxDiff(calculatedMass, massTolerance):

  minExp = \
    round(calculatedMass - ((massTolerance * calculatedMass) / 1000000))
  maxExp = \
    round(calculatedMass + ((massTolerance * calculatedMass) / 1000000))

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

def peptidesForSpectrum(spectrumData, peptideDict, massTolerance, maxDecoys):
  # Input: Dataframe with Spectrum, masses
  # Output: Dataframe with Spectrum, decoy peptides
  spectrumTitles = []
  decoyPeptides = []
  for spec, mass in zip(spectrumData[TITLE_SPECTRUM], spectrumData[PRECURSOR_MASS]):
    newDecoyPeptides = getDecoys(mass, peptideDict, massTolerance, maxDecoys)
    for i in range(len(newDecoyPeptides)):
      spectrumTitles.append(spec)
    decoyPeptides += newDecoyPeptides
  df = pd.DataFrame({TITLE_SPECTRUM: spectrumTitles, PEPTIDE: decoyPeptides})
  return df

def getPeptideDict(decoyPeptideDirectory,
                   precision,
                   minP,
                   maxP,
                   reverse):

  pepDict = {}

  for fileName in os.listdir(decoyPeptideDirectory):
    currentFile = os.path.join(decoyPeptideDirectory, fileName)
    pepDict = fileToDict(currentFile, precision, minP, maxP, pepDict, reverse)
  
  return pepDict

def fileToDict(peptideFileName,
               precision,
               minP,
               maxP,
               pepDict,
               reverse):
  
  with open(peptideFileName, 'r') as f:
    pepFile = f.read()
  pepFile = pepFile.split('\n')


  for pepLine in pepFile[1:]:
    if pepLine == "":
      break
    length, peptide, mass = pepLine.strip().split(',')
    mass = int(round((10 ** precision) * float(mass)))
    length = int(length)

    if reverse:
      peptide = peptide[::-1]

    if 'U' in peptide or 'X' in peptide:
      continue
    if length < minP or length > maxP:
      continue

    if mass in pepDict:
      pepDict[mass].append(peptide)
    else:
      pepDict[mass] = [peptide]

  return pepDict

def selectDecoyPeptides(decoyPeptideDirectory,
                        spectrumFileDirectory,
                        acidMassTable,
                        conversionTable,
                        massTolerance=35,
                        minPeptideLength=9,
                        maxPeptideLength=12,
                        maxDecoys=10,
                        precision=4,
                        reverse=False):
  
  spectrumData = extractSpectrumInformation(spectrumFileDirectory, precision)

  peptideDict = getPeptideDict(decoyPeptideDirectory,
                               precision,
                               minPeptideLength,
                               maxPeptideLength,
                               reverse)

  decoyDataframe = \
      peptidesForSpectrum(spectrumData, peptideDict, massTolerance, maxDecoys)
    
  return decoyDataframe


if __name__ == "__main__":

  acidMassTable = \
      AcidMassTable.adjustForPrecision(
          AcidMassTableIO.getAminoMasses(sys.argv[3]),
          4)
  conversionTable = \
      AcidConversion.createAcidConversionTable([acid for acid in acidMassTable])

  selectDecoyPeptides(sys.argv[1],
                      sys.argv[2],
                      sys.argv[3],
                      conversionTable)