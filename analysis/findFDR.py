#!/usr/bin/env python3

# Find the scores associated with a given FDR
# Right now, this file is hard-coded to use the combined peptide score
#
# Outputs a histogram with FDR scores marked with lines
#
# Take intersection of all spectrums in which both a decoy and results were found?
# Or use a score of 0 if no decoy was found?

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import sys, os
import argparse

# cheap hack until I can figure out how to do this properly in
# python 3.6
sys.path.append(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])

from backend.constants import *


RESULT_IMMUNO = 'IMMUNOSCORES'
RESULTS_DECOY = 'DECOYSCORES'
RESULTS_DELTA = 'SCOREDELTA'



def parseArguments():
  parser = argparse.ArgumentParser()
  parser.add_argument('-f', '--FDR',
                        dest='fdr',
                        default = 0.05,
                        help='The target FDR (default=0.05)', type=float)
  parser.add_argument('-pn', '--Plot-Name',
                        dest='plotName',
                        default='fdrPlot',
                        help='Output name of our FDR plot (default="fdrPlot")')
  parser.add_argument('-st', '--Score-Type',
                        dest='scoreType',
                        default=SCORE_COMBINED,
                        help='Score type to compare (defaults to combined score)')
  parser.add_argument('results_directory',
                        help='A directory of results files from ImmuNovo')
  parser.add_argument('decoy_file',
                        help='A decoy scores file from dbPepToScores')

  arguments = parser.parse_args()
  
  if not (0 < arguments.fdr < 1):
    print(str.format('{} is not a valid FDR. Should be between (0, 1)', arguments.fdr))
    exit()
  
  return arguments


def getResultsScores(dataFrame, scoreType):
  df_filtered = \
      dataFrame[dataFrame[PEPTIDE] != NO_PEP][[scoreType, TITLE_SPECTRUM]]

  df_grouped = \
      df_filtered.groupby(TITLE_SPECTRUM)[scoreType].transform(max) == \
          df_filtered[scoreType]

  return df_filtered[df_grouped][[scoreType, TITLE_SPECTRUM]]


def getDecoyScores(dataFrame, scoreType):
  df_filtered = \
      dataFrame[[scoreType, TITLE_SPECTRUM]]

  df_grouped = \
      df_filtered.groupby(TITLE_SPECTRUM)[scoreType].transform(max) == \
          df_filtered[scoreType]

  return df_filtered[df_grouped][[scoreType, TITLE_SPECTRUM]]


def getScores(resultsDF, decoyDF, scoreType, deltaThreshold=0):
  results = getResultsScores(resultsDF, scoreType)
  decoy = getDecoyScores(decoyDF, scoreType)

  resultsFiltered = \
      results[results[TITLE_SPECTRUM].isin(decoy[TITLE_SPECTRUM])]
  decoyFiltered = \
      decoy[decoy[TITLE_SPECTRUM].isin(resultsFiltered[TITLE_SPECTRUM])]

  merged = \
    pd.merge(resultsFiltered.rename(columns={scoreType: RESULT_IMMUNO}),
              decoyFiltered.rename(columns={scoreType: RESULTS_DECOY}))
  
  thresholdCalc = \
    lambda row: True if row[RESULT_IMMUNO] - row[RESULTS_DECOY] > deltaThreshold \
                     else False

  merged[RESULTS_DELTA] = \
      merged.apply(thresholdCalc, axis=1)

  return merged


def findFDR(dataFrame, FDR, fdrDelta=0.005, precision=3):
  
  def reduceFDR(fdrScores, fdrDelta, precision):
    fdrScores = [round(x, precision) for x in fdrScores]

    for i in range(len(fdrScores) - 1, 0, -1):
      if fdrScores[i] - fdrScores[i-1] <= fdrDelta:
        del fdrScores[i]
    
    return fdrScores
  

  immuNovoScores = list(mergedScores[RESULT_IMMUNO])
  deltas = [1 if x == True else 0 for x in list(mergedScores[RESULTS_DELTA])]
  forward = sum(deltas)
  decoys = len(deltas) - forward

  print(forward)
  print(decoys)

  fdrScores = []
  scoreCheck = list(zip(immuNovoScores, deltas))
  scoreCheck.sort(key=lambda x: x[0])

  for entry in scoreCheck:
    print(decoys)
    print(forward)
    print(decoys / forward)
    if decoys / forward <= FDR:
      fdrScores.append(entry[0])
    print(entry)
    input()
    forward -= entry[1]
    decoys -= (entry[1] - 1)
  
  return reduceFDR(fdrScores, fdrDelta, precision)


def plotResults(mergedScores,
                fileName,
                fdrScores,
                FDR):
  
  plt.hist(mergedScores[RESULT_IMMUNO],
            bins=50,
            histtype='step',
            color='blue',
            label='ImmuNovo')

  plt.hist(mergedScores[RESULTS_DECOY],
            bins=50,
            histtype='step',
            color='red',
            label='Decoys')

  blue_patch = mpatches.Patch(color="blue", label="ImmuNovo Scores")
  red_patch = mpatches.Patch(color="red", label="Decoy Scores")

  plt.legend(handles=[blue_patch, red_patch])
  plt.title('ImmuNovo Results vs Decoy')
  plt.xlabel(scoreType)
  plt.ylabel("Frequency")

  fdrScoresString = '\n'.join([str(round(x,3)) for x in fdrScores])

  plt.text(0.62,
           0.75,
           'FDR Scores at {}:\n{}'.format(FDR, fdrScoresString),
            va='top',
            ha='left',
            fontsize=12,
            transform=plt.gcf().transFigure)

  for fdrScore in fdrScores:
    plt.axvline(x=fdrScore, color='g', linestyle='dashed', linewidth=1)

  plt.savefig(fileName + '.png', dpi=300)
  #plt.show()


if __name__ == '__main__':

  ## Add PSSM detection to this?
  arguments = parseArguments()

  scoreType = arguments.scoreType
  plotName = arguments.plotName
  fdr = arguments.fdr

  resultsDirectory = arguments.results_directory
  combinedResults = []
  for fileName in os.listdir(resultsDirectory):
    combinedResults.append(pd.read_csv(os.path.join(resultsDirectory, fileName)))
  resultsDF = pd.concat(combinedResults)
  decoyDF = pd.read_csv(arguments.decoy_file)

  mergedScores = getScores(resultsDF, decoyDF, scoreType)
  
  fdrScores = findFDR(mergedScores, fdr)
  
  plotResults(mergedScores,
              plotName,
              fdrScores,
              fdr)
