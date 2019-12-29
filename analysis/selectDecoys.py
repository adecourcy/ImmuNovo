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
from backend.constants import *
from backend.userInput import *

MASS = 'MASS'
DECOYS = 'DECOYS'


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


def getAllDecoyPeptides(decoyPeptideDirectory):
  # Input: decoy database directory
  # Output: dataframe
  allData = []
  for fileName in os.listdir(decoyPeptideDirectory):
    currentFile = os.path.join(decoyPeptideDirectory, fileName)
    with open(currentFile, 'r') as f:
      allData += f.read().strip().split('\n')
  return pd.DataFrame({PEPTIDE: allData})


def filterByLength(df, minLength, maxLength):
  # Input: Dataframe with column of peptides
  # Output: Dataframe with column of peptides
  df = df[df.apply(lambda x: len(x[PEPTIDE]) <= maxLength, axis=1)]
  df = df[df.apply(lambda x: len(x[PEPTIDE]) >= minLength, axis=1)]
  return df


def validPeptide(peptide, validAcids):
  for acid in peptide:
    if acid not in validAcids:
      return False
  return True


def removeUnknownPeptides(df, acidMasstable):
  validAcids = AcidMassTable.getAcids(acidMasstable)
  return df[df.apply(lambda x: validPeptide(x[PEPTIDE], validAcids), axis=1)]


def calculateMasses(df, acidMassTable):
  # Input: Dataframe with column of peptides
  # Output: Dataframe with column of peptides and masses
  return df.apply(lambda x: AcidMassTable.peptideMass(x[PEPTIDE]), axis=1)


def separateByMass(df, precision):
  # Input: Dataframe with masses, peptides
  # Output: Dictionary with key masses (int, adjusted for precision), list peptides
  dataDict = {}
  df[MASS] = df.apply(lambda x: int(round(x[MASS] * (10**precision))), axis=1)
  for elm in df.groupby(MASS):
    dataDict[mass] = list(elm[1][PEPTIDE])
  return dataDict

def massToleranceMaxDiff(calculatedMass, massTolerance):

  minExp = int(round(calculatedMass / (massTolerance + 1)))
  maxExp = int(round(calculatedMass / (1 - massTolerance)))

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


def extractSpectrumInformation(spectrumFileDirectory):
  # Return a dataframe with spectrum title and mass columns
  allInfo = []
  def extractFileInformation(spectrumFile):
    precursorMasses = []
    titles = []
    for spectrum in SpectrumIO.getSpectrums(spectrumFile):
      precursorMasses.append(Spectrum.getPrecursorMass(spectrum))
      titles.append(Spectrum.getTitle(spectrum))
    return pd.DataFrame({TITLE_SPECTRUM: titles,
                         PRECURSOR_MASS: precursorMasses})
  
  for fileName in os.listdir(spectrumFileDirectory):
    allInfo.append(extractFileInformation(os.path.join(spectrumFileDirectory, fileName)))

  return pd.concat(allInfo)


def peptidesForSpectrum(df, peptideDict, massTolerance, maxDecoys):
  # Input: Dataframe with Spectrum, masses
  # Output: Dataframe with Spectrum, decoy peptides
  spectrumTitles = []
  decoyPeptides = []
  for spec, mass in zip(df[TITLE_SPECTRUM], df[MASS]):
    newDecoyPeptides = getDecoys(mass, peptideDict, massTolerance, maxDecoys)
    for i in range(len(newDecoyPeptides)):
      spectrumTitles.append(spec)
    decoyPeptides += newDecoyPeptides
  df = pd.DataFrame({TITLE_SPECTRUM: spectrumTitles, DECOYS: decoyPeptides})
  return df


def selectDecoyPeptides(decoyPeptideDirectory,
                        spectrumFileDirectory,
                        acidMassTable,
                        massTolerance,
                        minPeptideLength,
                        maxPeptideLength,
                        maxDecoys,
                        conversionTable,
                        precision,):
  pass