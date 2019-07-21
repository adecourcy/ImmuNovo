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
RESULTS_RATIO = 'SCORERATIO'



def parseArguments():
  parser = argparse.ArgumentParser()
  parser.add_argument('-f', '--FDR',
                        default = 0.05,
                        help='The target FDR (default=0.05)', type=float)
  parser.add_argument('results_directory',
                        help='A directory of results files from ImmuNovo')
  parser.add_argument('decoy_file',
                        help='A decoy scores file from dbPepToScores')

  arguments = parser.parse_args()
  
  if not (0 < arguments.FDR < 1):
    print(str.format('{} is not a valid FDR. Should be between (0, 1)', arguments.FDR))
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


def getScores(resultsDF, decoyDF, scoreType):
  results = getResultsScores(resultsDF, scoreType)
  decoy = getDecoyScores(decoyDF, scoreType)
  resultsFiltered = \
      results[results[TITLE_SPECTRUM].isin(decoy[TITLE_SPECTRUM])]
  decoyFiltered = \
      decoy[decoy[TITLE_SPECTRUM].isin(resultsFiltered[TITLE_SPECTRUM])]
      
  merged = \
    pd.merge(resultsFiltered.rename(columns={scoreType: RESULT_IMMUNO}),
              decoyFiltered.rename(columns={scoreType: RESULTS_DECOY}))

  merged[RESULTS_RATIO] = \
      merged.apply(lambda row: row[RESULT_IMMUNO] / row[RESULTS_DECOY], axis=1)

  return (resultsFiltered[scoreType],
          decoyFiltered[scoreType],
          merged)


def findFDR(resultScores, decoyScores, FDR):

  def findNextIndex(scoreList, currentScore):
    if currentScore > scoreList[-1]:
      return len(scoreList)
    for i in range(len(scoreList)):
      if scoreList[i] >= currentScore:
        return i
  
  def findNextScore(scoresList, currentScore, currentIndex):
    if currentScore >= scoresList[-1]:
      return currentScore, len(scoresList)
    for i in range(currentIndex, len(scoresList)):
      if scoresList[i] > currentScore:
        return scoresList[i], i
  
  def convertFDRScores(fdrIndex, resultScores):
    fdrScores = []

    # If we have multiple valid FDR scores in a row, take the lowest
    haveFound = False
    for i in range(len(fdrIndex)):
      if fdrIndex[i] == False:
        haveFound = False
        continue
      elif fdrIndex[i] == True and haveFound == True:
        continue
      else:
        haveFound = True
        fdrScores.append(resultScores[i])
    
    for i in range(len(fdrScores)-1, 0, -1):
      if (fdrScores[i] - fdrScores[i-1]) < 0.03:
        del fdrScores[i]

    return fdrScores


  sortedResults = sorted([round(x, 3) for x in list(resultScores)])
  sortedDecoy = sorted([round(x,3) for x in list(decoyScores)])

  currentScore = sortedResults[0]
  fdrIndex = [False for x in range(len(resultScores))]

  resultsIndex = 0
  decoyIndex = 0

  while resultsIndex < len(sortedResults):
    if sortedDecoy[decoyIndex] < currentScore:
      decoyIndex = findNextIndex(sortedDecoy, currentScore)
      if decoyIndex >= len(sortedDecoy):
        # Everything following has an FDR less than our desired rate
        for i in range(resultsIndex, len(fdrIndex)):
          fdrIndex[i] = True
        break

    resultSum = len(sortedResults[resultsIndex:])
    decoySum = len(sortedDecoy[decoyIndex:])
    currentFDR = round(decoySum / (resultSum + decoySum), 2)
    
    if currentFDR <= FDR:
      fdrIndex[resultsIndex] = True
    
    currentScore, resultsIndex = \
        findNextScore(sortedResults, currentScore, resultsIndex)

  return convertFDRScores(fdrIndex, sortedResults), sortedResults, sortedDecoy


def plotResults(resultScores,
                decoyScores,
                sortedResults,
                sortedDecoys,
                scoreType,
                fileName,
                fdrScores,
                FDR):
  
  plt.hist(resultScores,
            bins=50,
            histtype='step',
            color='blue',
            label='ImmuNovo')

  plt.hist(decoyScores,
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

  scoreType = SCORE_COMBINED
  plotName = 'FDRPlot'
  ## Add PSSM detection to this
  arguments = parseArguments()

  resultsDirectory = arguments.results_directory
  combinedResults = []
  for fileName in os.listdir(resultsDirectory):
    combinedResults.append(pd.read_csv(os.path.join(resultsDirectory, fileName)))
  resultsDF = pd.concat(combinedResults)
  decoyDF = pd.read_csv(arguments.decoy_file)

  resultScores, decoyScores, mergedScores = \
          getScores(resultsDF, decoyDF, scoreType)
  
  plt.hist(mergedScores[RESULTS_RATIO])
  plt.show()
  exit()
  
  fdrScores, sortedResults, sortedDecoys = \
          findFDR(resultScores, decoyScores, arguments.FDR)
  
  plotResults(resultScores,
              decoyScores,
              sortedResults,
              sortedDecoys,
              scoreType,
              plotName,
              fdrScores,
              arguments.FDR)
