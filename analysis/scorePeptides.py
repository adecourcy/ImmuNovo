#!/usr/bin/env python3
# Give a set of peptides, PSSMs, and associated Spectra, score peptides

# cheap hack until I can figure out how to do this properly in
# python 3.6
import sys, os
sys.path.append(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])

import argparse
import pandas as pd
import numpy as np

import backend.userInput as userInput
import backend.Structures.spectrum as Spectrum
import backend.Structures.spectrumIO as SpectrumIO
import backend.PreProcessing.spectrumConversion as SpectrumConversion
import backend.Structures.pssm as PSSM
#import backend.PostProcessing.processResults as processResults
import backend.PostProcessing.processResultsOld as processResults
from main import getAminoVariables
from backend.constants import *


def parseArguments():
  abspath = os.path.abspath(__file__)
  dname = os.path.dirname(abspath)

  parser = argparse.ArgumentParser()
  parser.add_argument('-o', '--Output-File',
                      dest='output_file',
                      default=os.path.join(dname, 'scored_peptides'))
  parser.add_argument('peptide_file',
                        help='A csv of peptides to be scored')
  parser.add_argument('spec_dir',
                        help='A directory of spectra to use in scoring')
  parser.add_argument('pssm_dir',
                        help='A directory of PSSM files for scoring')
  parser.add_argument('acid_mass_file',
      help='A file containing mass spectrometry data')
  parser = userInput.parseOptionalArguments(parser)

  arguments = parser.parse_args()
  arguments.peptide_file = os.path.abspath(arguments.peptide_file)
  arguments.spec_dir = os.path.abspath(arguments.spec_dir)
  arguments.pssm_dir = os.path.abspath(arguments.pssm_dir)

  if not arguments.output_file.endswith('.csv'):
    arguments.output_file += '.csv'
  
  return arguments


def deConvertPeptideString(peptideString, acidConversion):

  reverseConversionTable = {acidConversion[x]: x for x in acidConversion}

  convertedPeptide = ''
  for acid in peptideString:
    convertedPeptide += reverseConversionTable[acid]
  
  return convertedPeptide


def convertPeptidString(peptideString, acidConversion):
  acidConversion = [(x, acidConversion[x]) for x in acidConversion]
  acidConversion.sort(key=lambda x: len(x[0]), reverse=True)
  for item in acidConversion:
    peptideString = peptideString.replace(item[0], item[1])

  return peptideString


def pssmScore(peptideString, allPSSM, title):
  matrix = PSSM.getMatrixOfLength(title, len(peptideString), allPSSM)
  score = 0
  for i in range(len(peptideString)):
    score += PSSM.getAcidProbabilities(matrix, peptideString[i])[i]
  score /= len(peptideString)
  return score


def globalScore(acidMassTable,
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


def scorePeptides(peptideDF,
                  acidMassFile,
                  pssmDir,
                  spectrumDirectory,
                  precision,
                  minPepLength,
                  maxPepLength,
                  maxMassTolerance,
                  compression):

  peptideDF = peptideDF[peptideDF[PEPTIDE] != NO_PEP]
  pepBySpectrum = {}

  allPSSM, acidMassTable, conversionTable = \
                              getAminoVariables(acidMassFile,
                                                precision,
                                                pssmDir,
                                                minPepLength,
                                                maxPepLength)
  peptideDF[PEPTIDE] = \
    peptideDF.apply(lambda row: convertPeptidString(row[PEPTIDE],
                    conversionTable),
                    axis=1)
  peptideDF['pepFilter'] = \
    peptideDF.apply(lambda row: False if len(row[PEPTIDE]) < minPepLength or \
                                         len(row[PEPTIDE]) > maxPepLength else True,
                    axis=1)
  peptideDF = peptideDF[peptideDF['pepFilter'] == True]
  peptideDF = peptideDF.drop(columns=['pepFilter'])

  for item in peptideDF.groupby(by=TITLE_SPECTRUM):
    pepBySpectrum[item[0]] = item[1]

  H2OMassAdjusted = int(H2OMASS * (10**precision))
  NH3MassAdjusted = int(NH3MASS * (10**precision))
  protonMassAdjusted = int(PROTONMASS * (10**precision))

  processedPeptides = []

  for spectrumFileName in os.listdir(spectrumDirectory):
    spectrumFile = os.path.join(spectrumDirectory, spectrumFileName)
    for spectrum in SpectrumIO.getSpectrums(spectrumFile):
      spectrumTitle = Spectrum.getTitle(spectrum)
      # I don't remember why I was doing this
      spectrumTitle = spectrumTitle.replace(",", '').split()[0] 
      if spectrumTitle not in pepBySpectrum:
        continue
      spectrumPeptides = pepBySpectrum[spectrumTitle]

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
      # applyGlobalScore = lambda x: globalScore(acidMassTable,
      #                                          experimentalSpectrum,
      #                                          experimentalIntensities,
      #                                          protonMassAdjusted,
      #                                          H2OMassAdjusted,
      #                                          NH3MassAdjusted,
      #                                          maxMassTolerance,
      #                                          x)

      applyGlobalScore = lambda x: globalScore(acidMassTable,
                                               experimentalSpectrum,
                                               experimentalIntensities,
                                               0.5,
                                               protonMassAdjusted,
                                               H2OMassAdjusted,
                                               NH3MassAdjusted,
                                               maxMassTolerance,
                                               x)

      spectrumPeptides[SCORE_GLOBAL] = \
          spectrumPeptides.apply(lambda row: applyGlobalScore(row[PEPTIDE]), axis=1)

      pssmTitles = [title for title in allPSSM]
      spectrumPeptides = \
        (spectrumPeptides.loc[spectrumPeptides.index.repeat(len(pssmTitles))]
                         # place holder column names
                         .assign(pt = np.tile(pssmTitles, len(spectrumPeptides)))
                         .assign(ps = lambda x: [pssmScore(peptideString, allPSSM, title) for peptideString, title in zip(x[PEPTIDE], np.tile(pssmTitles, len(spectrumPeptides)))])
                         )
      spectrumPeptides = \
        spectrumPeptides.rename(columns = {'pt': TITLE_PSSM, 'ps': SCORE_PSSM})
      spectrumPeptides[SCORE_PSSM] /= 10**precision
      spectrumPeptides[SCORE_COMBINED] = \
            spectrumPeptides[SCORE_GLOBAL] * spectrumPeptides[SCORE_PSSM]
      
      processedPeptides.append(spectrumPeptides)
  
  peptideDF = pd.concat(processedPeptides)
  peptideDF[PEPTIDE] = \
    peptideDF.apply(lambda row: deConvertPeptideString(row[PEPTIDE],
                                                       conversionTable),
                    axis=1)

  return peptideDF


if __name__ == '__main__':
  """
  pass a CSV with headers for spectra and peptides, headers should conform to
  the constants file standards

  Pass directory of spectra files, and a directory of PSSMs
  """
  arguments = parseArguments()

  peptideDF = scorePeptides(pd.read_csv(arguments.peptide_file),
                            arguments.acid_mass_file,
                            arguments.pssm_dir,
                            arguments.spec_dir,
                            arguments.prec,
                            arguments.minP,
                            arguments.maxP,
                            arguments.mmt,
                            arguments.comp)

  peptideDF.to_csv(arguments.output_file, index=False)
