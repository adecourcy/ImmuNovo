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


def getBestScoringEntries(dataFrame, scoreType):
  df_filtered = \
      dataFrame[[scoreType, TITLE_SPECTRUM]]

  # Create a dataframe of bools which indicate whether each entry is the max
  # score for that spectrum
  df_grouped = \
      df_filtered.groupby(TITLE_SPECTRUM)[scoreType].transform(max) == \
          df_filtered[scoreType]

  # return only the maximum scoring entry for each spectrum
  return df_filtered[df_grouped][[scoreType, TITLE_SPECTRUM]]

def hitOrMissDataframe(resultsDF, decoyDF, scoreType):

  resultsDF = \
      resultsDF[resultsDF[TITLE_SPECTRUM].isin(decoyDF[TITLE_SPECTRUM])]
  decoyDF = \
      decoyDF[decoyDF[TITLE_SPECTRUM].isin(resultsDF[TITLE_SPECTRUM])]

  merged = \
    pd.merge(resultsDF.rename(columns={scoreType: RESULT_IMMUNO})[[TITLE_SPECTRUM, RESULT_IMMUNO]],
              decoyDF.rename(columns={scoreType: RESULTS_DECOY})[[TITLE_SPECTRUM, RESULTS_DECOY]])
  
  thresholdCalc = \
    lambda row: True if row[RESULT_IMMUNO] - row[RESULTS_DECOY] >= 0 \
                     else False

  if len(merged) != 0:
    merged[RESULTS_DELTA] = \
        merged.apply(thresholdCalc, axis=1)

  return merged

def dynamicFDR(maxFDR, scoreList, calculatedFDRs, increment=0.01, scoreIndex=0, fdrIndex=1):
  # Given a sets of scores vs FDR calculations, use dynamic programming to find
  # an optimal, monotonically decreasing step function to match scores vs FDRs

  # Assumes scoreList and calculatedFDRs are ordered by greatest to least based
  # on scoreList entries

  # Return (score, fdr) sorted lowest to highest by fdr

  # def getMinScoreIndex(calcs):
  #   minScore = math.inf
  #   minIndex = -1
  #   for i in range(len(calcs)):
  #     if calcs[i] < minScore:
  #       minScore = calcs[i]
  #       minIndex = i
  #   return minScore, minIndex
  

  def createThresholdList(theoreticalFDRs):
    # Create a matrix of optimal threshold scores for each FDR
    thresholds = {}
    for tFDR in theoreticalFDRs:
      thresholds[tFDR] = 0
    thresholdList = {}
    for tFDR in theoreticalFDRs:
      thresholdList[tFDR] = deepcopy(thresholds)
    return thresholdList

  maxFDR = round(maxFDR, 2)
  # Rounded version of the actual, calculated FDRs
  calculatedFDRs = [round(x, 2) for x in calculatedFDRs]
  # List going from the min to max FDR by "increment"
  theoreticalFDRs = [round((x * increment), 2) for x in range(1, round((maxFDR+increment)/increment))]


  # We don't need to track paths, only score to FDR thresholds
  prevThresholdList = createThresholdList(theoreticalFDRs)

  # Dinstance between the calculate FDR, and each potential FDR label
  # for the first entry. Used to start off the dynamic programming
  prevCalc = [abs(calculatedFDRs[0] - x) for x in theoreticalFDRs]

  scoreIndexTracker = [0 for i in range(len(theoreticalFDRs))]

  for tIndex in range(len(theoreticalFDRs)):
    currentScore = abs(calculatedFDRs[0] - theoreticalFDRs[tIndex])
    if tIndex == 0:
      scoreIndexTracker[0] = (currentScore, 0)
    elif currentScore < oldScore:
      scoreIndexTracker[tIndex] = (currentScore, tIndex)
    else:
      scoreIndexTracker[tIndex] = (oldScore, oldIndex)
    oldScore, oldIndex = scoreIndexTracker[tIndex]


  for cIndex in range(1, len(calculatedFDRs)):
    cFDR = calculatedFDRs[cIndex] #### !!!!! Is this right? Should be cIndex-1? !!!!! 
    prevScore = scoreList[cIndex - 1]
    newThresholdList = {}
    for tIndex in range(len(theoreticalFDRs)):
      #minPrevScore, minPrevIndex = getMinScoreIndex(prevCalc[:tIndex+1])
      minPrevScore, minPrevIndex = scoreIndexTracker[tIndex]

      tFDR = theoreticalFDRs[tIndex]
      prevFDR = theoreticalFDRs[minPrevIndex]

      currentScore = abs(cFDR - tFDR) + minPrevScore
      if tIndex == 0:
        scoreIndexTracker[0] = (currentScore, 0)
      elif currentScore < oldScore:
        scoreIndexTracker[tIndex] = (currentScore, tIndex)
      else:
        scoreIndexTracker[tIndex] = (oldScore, oldIndex)
      oldScore, oldIndex = scoreIndexTracker[tIndex]

      newThresholdList[tFDR] = deepcopy(prevThresholdList[prevFDR])
      newThresholdList[tFDR][prevFDR] = prevScore
    
    prevThresholdList = newThresholdList
  

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

def findFDR(mergedScores, precision=3):

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
        return round(elm[1], 2)
        # pretty sure this was a bug
        # return round(elm[1] + increment, 2)

  # sort (score, fdr) from highest to lowest by score
  fdrCutoffs.sort(key=lambda x: x[0], reverse=True)
  
  dataFrame[FDR] = \
    dataFrame.apply(lambda row: findCutoff(row[scoreType], fdrCutoffs, increment), axis=1)
  
  return dataFrame

def addFDRToDataframe(targetDF, decoyDF, scoreType, precision=3, increment=0.01):

  hitOrMissDF = hitOrMissDataframe(targetDF, decoyDF, scoreType)
  if len(hitOrMissDF) == 0:
    return hitOrMissDF
  fdrCutoffs = findFDR(hitOrMissDF, precision)
  fdrDF = addFDR(targetDF, fdrCutoffs, scoreType, increment)

  return fdrDF


def findFDROld(targetScores, decoyScores):

  def nextDecoyIndex(currentDecoyIndex, decoyScores, targetScore):
    if currentDecoyIndex == -1:
      return -1
    while decoyScores[currentDecoyIndex] > targetScore:
      currentDecoyIndex += 1
      if currentDecoyIndex >= len(decoyScores):
        return -1
    return currentDecoyIndex

  targetScores.sort(reverse=True)
  decoyScores.sort(reverse=True)
  fdrList = []

  targetIndex = 0
  decoyIndex = 0

  numTargetPeptides = len(targetScores)
  numDecoyPeptides = len(decoyScores)

  while targetIndex < len(targetScores):
    decoyIndex = \
        nextDecoyIndex(decoyIndex, decoyScores, targetScores[targetIndex])
    if decoyIndex == -1:
      fdrList.append(0.0)
    else:
      fdrList.append((numDecoyPeptides - decoyIndex) / \
            ((numDecoyPeptides - decoyIndex) + (numTargetPeptides - targetIndex)))
    targetIndex += 1
  
  return dynamicFDR(max(fdrList), targetScores, fdrList)


def addFDRToDataframeOld(targetDF, decoyDF, scoreType, precision=3, increment=0.01):
  
  fdrCutoffs = findFDROld(list(targetDF[scoreType]), list(decoyDF[scoreType]))
  return addFDR(targetDF, fdrCutoffs, scoreType, increment)
  



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
