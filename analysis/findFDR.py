# Find the scores associated with a given FDR
# Right now, this file is hard-coded to use the combined peptide score
#
# Outputs a histogram with FDR scores marked with lines
#
# Take intersection of all spectrums in which both a decoy and results were found?
# Or use a score of 0 if no decoy was found?

import matplotlib.pyplot as plt
import pandas as pd
import sys, os
import argparse

# cheap hack until I can figure out how to do this properly in
# python 3.6
sys.path.append(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])

from backend.constants import *



def parseArguments():
  parser = argparse.ArgumentParser()
  parser.add_argument('-f', '--FDR',
                        default = 0.05,
                        help='The target FDR (default=0.05)', type=float)
  parser.add_argument('results_file',
                        help='The results file from ImmuNovo')
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

  return df_filtered[df_grouped][scoreType]


def getDecoyScores(dataFrame, scoreType):
  return dataFrame[scoreType]


def getScores(resultsDF, decoyDF, scoreType):
  return (getResultsScores(resultsDF, scoreType),
          getDecoyScores(decoyDF, scoreType))


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

    return fdrScores


  sortedResults = sorted(list(resultScores))
  sortedDecoy = sorted(list(decoyScores))

  currentScore = sortedResults[0]
  fdrIndex = [False for x in range(len(resultScores))]

  resultsIndex = 0
  decoyIndex = 0

  while resultsIndex < len(sortedResults):
    if sortedDecoy[decoyIndex] < currentScore:
      decoyIndex = findNextIndex(sortedDecoy, currentScore)
      if decoyIndex == len(sortedDecoy):
        # Everything following has an FDR less than our desired rate
        for i in range(resultsIndex, len(fdrIndex)):
          fdrIndex[i] = True
        break

    resultSum = len(sortedResults[resultsIndex:])
    currentFDR = resultSum / (resultSum + len(sortedDecoy[decoyIndex:]))
    
    if currentFDR <= FDR:
      fdrIndex[resultsIndex] = True
    
    currentScore, resultsIndex = \
        findNextScore(sortedResults, currentScore, resultsIndex)

  return convertFDRScores(fdrIndex, sortedResults), sortedResults, sortedDecoy


def plotResults(resultScores,
                decoyScores,
                sortedResults,
                sortedDecoys,
                fdrScores,
                FDR):
  
  plt.hist(resultScores, histtype='step', color='blue')
  plt.hist(decoyScores, histtype='step', color='red')
  plt.show()


if __name__ == '__main__':
  ## Add PSSM detection to this
  arguments = parseArguments()

  resultsDF = pd.read_csv(arguments.results_file)
  decoyDF = pd.read_csv(arguments.decoy_file)

  resultScores, decoyScores = getScores(resultsDF,
                                          decoyDF,
                                          SCORE_COMBINED)
  
  fdrScores, sortedResults, sortedDecoys = \
          findFDR(resultScores, decoyScores, arguments.FDR)
  print(fdrScores)
  
  plotResults(resultScores,
              decoyScores,
              sortedResults,
              sortedDecoys,
              fdrScores,
              arguments.FDR)
