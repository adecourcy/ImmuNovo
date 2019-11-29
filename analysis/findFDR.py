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
import math

from copy import deepcopy

# cheap hack until I can figure out how to do this properly in
# python 3.6
sys.path.append(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])

from backend.constants import *


def parseArguments():
  abspath = os.path.abspath(__file__)
  dname = os.path.dirname(abspath)

  parser = argparse.ArgumentParser()
  parser.add_argument('-f', '--FDR',
                        dest='fdr',
                        default = 0.01,
                        help='The target FDR', type=float)
  # parser.add_argument('-pn', '--Plot-Name',
  #                       dest='plotName',
  #                       default='fdrPlot',
  #                       help='Output name of our FDR plot (default="fdrPlot")')
  parser.add_argument('-o', '--Output-File',
                      dest='output_file',
                      default=os.path.join(dname, 'added_fdrs.csv'))
  parser.add_argument('-st', '--Score-Type',
                        dest='scoreType',
                        default=SCORE_COMBINED,
                        help='Score type to compare (defaults to combined score)')
  parser.add_argument('results_directory',
                        help='A directory of results files from ImmuNovo')
  parser.add_argument('decoy_file',
                        help='A decoy scores file from dbPepToScores')

  arguments = parser.parse_args()
  arguments.results_directory = os.path.abspath(arguments.results_directory)
  arguments.decoy_file = os.path.abspath(arguments.decoy_file)

  if not arguments.output_file.endswith('.csv'):
    arguments.output_file += '.csv'
  
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
    lambda row: True if row[RESULT_IMMUNO] - row[RESULTS_DECOY] >= deltaThreshold \
                     else False

  merged[RESULTS_DELTA] = \
      merged.apply(thresholdCalc, axis=1)

  return merged


def dynamicFDR(maxFDR, scoreList, calculatedFDRs, increment=0.01, scoreIndex=0, fdrIndex=1):
  # Given a sets of scores vs FDR calculations, use dynamic programming to find
  # an optimal, monotonically decreasing step function to match scores vs FDRs

  def getMinScoreIndex(calcs):
    minScore = math.inf
    minIndex = -1
    for i in range(len(calcs)):
      if calcs[i] < minScore:
        minScore = calcs[i]
        minIndex = i
    return minScore, minIndex
  

  def createThresholdList(theoreticalFDRs):
    thresholds = {}
    for tFDR in theoreticalFDRs:
      thresholds[tFDR] = 0
    thresholdList = {}
    for tFDR in theoreticalFDRs:
      thresholdList[tFDR] = deepcopy(thresholds)
    return thresholdList

  maxFDR = round(maxFDR, 2)
  calculatedFDRs = [round(x, 2) for x in calculatedFDRs]
  theoreticalFDRs = [round((x * increment), 2) for x in range(1, round((maxFDR+increment)/increment))]
  # print(calculatedFDRs)
  # input()
  # print(theoreticalFDRs)
  # input()

  # (pathIndex, calculation)
  prevCalc = [abs(calculatedFDRs[0] - x) for x in theoreticalFDRs]

  # We don't need to track paths, only score to FDR thresholds
  prevThresholdList = createThresholdList(theoreticalFDRs)

  for cIndex in range(1, len(calculatedFDRs)):
    currentCalc = []
    cFDR = calculatedFDRs[cIndex]
    prevScore = scoreList[cIndex - 1]
    #print(prevScore)
    newThresholdList = {}
    # print(prevCalc)
    # input()
    for tIndex in range(len(theoreticalFDRs)):
      minPrevScore, minPrevIndex = getMinScoreIndex(prevCalc[:tIndex+1])
      tFDR = theoreticalFDRs[tIndex]
      prevFDR = theoreticalFDRs[minPrevIndex]
      # print(minPrevScore)
      # print(minPrevIndex)
      # print(prevFDR)
      # input()

      currentCalc.append(abs(cFDR - tFDR) + minPrevScore)
      newThresholdList[tFDR] = deepcopy(prevThresholdList[prevFDR])
      newThresholdList[tFDR][prevFDR] = prevScore
    
    prevThresholdList = deepcopy(newThresholdList)
    prevCalc = deepcopy(currentCalc)
  
  finalThresholds = prevThresholdList[maxFDR]
  if finalThresholds[maxFDR] == 0:
    finalThresholds[maxFDR] = scoreList[-1]
  
  thresholdList = []
  for fdr in finalThresholds:
    if finalThresholds[fdr] == 0:
      continue
    else:
      thresholdList.append((finalThresholds[fdr], fdr))
  
  thresholdList.sort(key=lambda x: x[1])

  return thresholdList


def findFDR(dataFrame, FDR, fdrDelta=0.005, precision=3):

  immuNovoScores = list(mergedScores[RESULT_IMMUNO])
  deltas = [1 if x == True else 0 for x in list(mergedScores[RESULTS_DELTA])]
  forward = sum(deltas) # How many scores >= decoy scores
  decoys = len(deltas) - forward # How many scores < decoy scores

  fdrList = []
  scoreCheck = list(zip(immuNovoScores, deltas))
  scoreCheck.sort(key=lambda x: x[0])
  immuNovoScores.sort()

  for entry in scoreCheck:
    fdrList.append(decoys / forward)
    forward -= entry[1]
    decoys -= (1 - entry[1])

  
  immuNovoScores, fdrList = zip(*sorted(zip(immuNovoScores, fdrList), reverse=True))
  maxFDR = max(fdrList)

  return dynamicFDR(maxFDR, immuNovoScores, fdrList)


def addFDR(dataFrame, fdrCutoffs, scoreType, increment=0.01):
  # fdrCutoffs = (scores, fdr)

  def findCutoff(score, fdrCutoffs, increment):
    for elm in fdrCutoffs:
      if score > elm[0]:
        return round(elm[1] + increment, 2)

  fdrCutoffs.sort(key=lambda x: x[0], reverse=True)
  
  dataFrame[FDR] = \
    dataFrame.apply(lambda row: findCutoff(row[scoreType], fdrCutoffs, increment), axis=1)
  
  return dataFrame


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
  #plotName = arguments.plotName
  fdr = arguments.fdr

  resultsDirectory = arguments.results_directory
  combinedResults = []
  for fileName in os.listdir(resultsDirectory):
    combinedResults.append(pd.read_csv(os.path.join(resultsDirectory, fileName)))
  resultsDF = pd.concat(combinedResults)
  decoyDF = pd.read_csv(arguments.decoy_file)

  mergedScores = getScores(resultsDF, decoyDF, scoreType)
  
  fdrCutoffs = findFDR(mergedScores, fdr)

  finalResults = addFDR(resultsDF, fdrCutoffs, scoreType)
  finalResults[finalResults[PEPTIDE] != NO_PEP].to_csv(arguments.output_file)
  
  # plotResults(mergedScores,
  #             plotName,
  #             fdrScores,
  #             fdr)
