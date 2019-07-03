from sys import argv
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
from decimal import Decimal

def getCurrentSpectrum(file, fileLine):
  if fileLine >= len(file):
    return ''
  splitLine = file[fileLine].split(",")
  return splitLine[1]


def getNextBatch(file, fileLine):

  currentSpectrum = getCurrentSpectrum(file, fileLine)

  results = []

  while(True):
    if fileLine >= len(file):
      break
    if currentSpectrum == getCurrentSpectrum(file, fileLine):
      results.append(file[fileLine])
      fileLine += 1
    else:
      break

  return results, fileLine


def checkMatch(oPep, fPep):
  if oPep == fPep:
    return True
  else:
    if len(oPep) != len(fPep):
      return False
    else:
      # Check for I/L differences
      #print(oPep)
      #print(fPep)
      newPep = ""
      for amino, i in zip(oPep, range(len(oPep))):
        if amino == "I":
          newPep = newPep + "L"
        else:
          newPep = newPep + amino
      oPep = newPep
      newPep = ""
      for amino, i in zip(fPep, range(len(fPep))):
        if amino == "I":
          newPep = newPep + "L"
        else:
          newPep = newPep + amino
      fPep = newPep
      misMatches = 0
      misIndex = []
      if oPep == fPep:
        #print("I/L!")
        return True

      #check for swaps

      for oA, fA, i in zip(oPep, fPep, range(len(oPep))):
        if oA != fA:
          misMatches += 1
          misIndex.append(i)
      if misMatches == 2:
        if oPep[misIndex[0]] == fPep[misIndex[1]] and \
           oPep[misIndex[1]] == fPep[misIndex[0]]:
          #print("swap!")
          return True
        else:
          return False

      if misMatches == 4:
        if oPep[misIndex[0]] == fPep[misIndex[1]] and \
           oPep[misIndex[1]] == fPep[misIndex[0]] and \
           oPep[misIndex[2]] == fPep[misIndex[3]] and \
           oPep[misIndex[3]] == fPep[misIndex[2]]:
          #print("swap!")
          return True
        else:
          return False

      else:
        return False

def getPepLength(peptide):
  return len(peptide.replace("M+15.995", "M"))

# def getBestScores(falseBatch, originalBatch,
#                   originalGlobalScore, falseGlobalScore,
#                   originalAminoScore, falseAminoScore,
#                   matchGlobalScore, matchAminoScore,
#                   misGlobal, misAmino,
#                   totalGlobal, totalAmino,
#                   uniquePeptides, aminoScoreCutoff):
  
#   for i in range(len(originalBatch)-1, -1, -1):
#     if "No Peptide Found" in originalBatch[i]:
#       del originalBatch[i]
  
#   for i in range(len(falseBatch)-1, -1, -1):
#     if "No Peptide Found" in falseBatch[i]:
#       del falseBatch[i]

#   originalBatch = [x.split(",") for x in originalBatch]
#   falseBatch = [x.split(",") for x in falseBatch]

#   #print("original batch: \n" + str(originalBatch))
#   #print("original global score:\n" + str(originalGlobalScore))
#   #print("original amino score:\n" + str(originalAminoScore))
#   #print("false batch: \n" + str(falseBatch))
#   #print("false global score:\n" + str(falseGlobalScore))
#   #print("false amino score:\n" + str(falseAminoScore))
#   #input()

#   if originalBatch != []:
#     originalBatch.sort(key=lambda x: float(x[5]))
#     originalGlobalPeptide = originalBatch[-1][2]

#   if falseBatch != []:
#     falseBatch.sort(key=lambda x: float(x[5]))
#     falseGlobalPeptide = falseBatch[-1][2]

#   if originalBatch != []:
#     if falseBatch == []:
#       matchGlobalScore.append(float(originalBatch[-1][5]))
#       originalGlobalScore.append(float(originalBatch[-1][5]))

#       ###############################################################
#       # Modified code to sort by global score and compare amino score
#       ###############################################################
#       matchAminoScore.append(float(originalBatch[-1][-1]) /
#                              getPepLength(originalGlobalPeptide))
#       originalAminoScore.append(float(originalBatch[-1][-1]) /
#                                 getPepLength(originalGlobalPeptide))

#       if originalGlobalPeptide not in uniquePeptides and \
#          originalAminoScore[-1] >= aminoScoreCutoff:
#         uniquePeptides.append(originalGlobalPeptide)
      
#       totalAmino += 1
#       ###############################################################
#       # Modified code to sort by global score and compare amino score
#       ###############################################################

#       totalGlobal += 1
#     elif checkMatch(originalGlobalPeptide, falseGlobalPeptide):


#     else:
#       #print(originalGlobalPeptide + "  " + falseGlobalPeptide)
#       originalGlobalScore.append(float(originalBatch[-1][5]))
#       falseGlobalScore.append(float(falseBatch[-1][5]))


#       ###############################################################
#       # Modified code to sort by global score and compare amino score
#       ###############################################################
#       originalAminoScore.append(float(originalBatch[-1][-1]))
#       falseAminoScore.append(float(falseBatch[-1][-1]))

#       if originalGlobalPeptide not in uniquePeptides and \
#          originalAminoScore[-1] >= aminoScoreCutoff:
#         uniquePeptides.append(originalGlobalPeptide)
        
#       totalAmino += 1
#       misAmino += 1
#       ###############################################################
#       # Modified code to sort by global score and compare amino score
#       ###############################################################

#       totalGlobal +=1
#       misGlobal += 1
  
#   elif falseBatch != []:
#     misGlobal += 1
#     totalGlobal += 1
#     falseGlobalScore.append(float(falseBatch[-1][5]))

#     ###############################################################
#     # Modified code to sort by global score and compare amino score
#     ###############################################################
#     falseAminoScore.append(float(falseBatch[-1][-1]))
#     misAmino += 1
#     totalAmino += 1
#     ###############################################################
#     # Modified code to sort by global score and compare amino score
#     ###############################################################
  
        
#   return (misGlobal, misAmino, totalGlobal, totalAmino)

def getBestCombineScores(trueBatch, falseBatch,
                         combinedTrue, combinedFalse,
                         pepMatchNumber):

  for i in range(len(trueBatch)-1, -1, -1):
    if "No Peptide Found" in trueBatch[i]:
      del trueBatch[i]


  for i in range(len(falseBatch)-1, -1, -1):
    if "No Peptide Found" in falseBatch[i]:
      del falseBatch[i]


  if trueBatch == [] and falseBatch == []:
    return (combinedTrue,
            combinedFalse,
            pepMatchNumber)


  if trueBatch != []:

    trueBatch = [x.split(",") for x in trueBatch]

    trueBatch.sort(key=lambda x: float(x[5]), reverse=True)
    truePeptide = trueBatch[0][2]
    trueScore = float(trueBatch[0][5]) * \
                    (float(trueBatch[0][9]) / \
                     float(getPepLength(truePeptide)))
    combinedTrue.append(trueScore)


  if falseBatch != []:

    falseBatch = [x.split(",") for x in falseBatch]

    falseBatch.sort(key=lambda x: float(x[5]), reverse=True)
    falsePeptide = falseBatch[0][2]
    falseScore = float(falseBatch[0][5]) * \
                    (float(falseBatch[0][9]) / \
                     float(getPepLength(falsePeptide)))

    if trueBatch != [] and checkMatch(truePeptide, falsePeptide):
      pepMatchNumber += 1
      if len(falseBatch) > 1:
        falsePeptide = falseBatch[1][2]
        falseScore = float(falseBatch[1][5]) * \
                        (float(falseBatch[1][9]) / \
                         float(getPepLength(falsePeptide)))
      else:
        falseScore = 0

    if falseScore != 0:
      combinedFalse.append(falseScore)


  return(combinedTrue, combinedFalse, pepMatchNumber)


def getBestCombineScoresRandDecoy(trueBatch, falseBatch,
                                  combinedTrue, combinedFalse,
                                  totalTrue, totalFalse,
                                  pepSpecAboveFDR):

  for i in range(len(trueBatch)-1, -1, -1):
    if "No Peptide Found" in trueBatch[i]:
      del trueBatch[i]


  for i in range(len(falseBatch)-1, -1, -1):
    if "No Peptide Found" in falseBatch[i]:
      del falseBatch[i]


  if trueBatch == [] and falseBatch == []:
    return (combinedTrue,
            combinedFalse,
            totalTrue,
            totalFalse)


  if trueBatch != []:

    trueBatch = [x.split(",") for x in trueBatch]

    trueBatch.sort(key=lambda x: float(x[5]), reverse=True)
    truePeptide = trueBatch[0][2]
    trueScore = float(trueBatch[0][5]) * \
                    (float(trueBatch[0][9]) / \
                     float(getPepLength(truePeptide)))

    #if trueScore >= 0.0046:
      #print(trueBatch[0][1] + ',' + truePeptide + '\n')
    combinedTrue.append(trueScore)
    totalTrue += 1


  if falseBatch != []:

    falseBatch = [x.split(",") for x in falseBatch]
    falseScore = float(falseBatch[0][2])

    if falseScore != 0:
      totalFalse += 1
      combinedFalse.append(falseScore)


  return (combinedTrue, combinedFalse, totalTrue, totalFalse)


def getMatchedScores(trueBatch, falseBatch,
                     matchedTrue, matchedFalse,
                     pepMatchNumber):

  for i in range(len(trueBatch)-1, -1, -1):
    if "No Peptide Found" in trueBatch[i]:
      del trueBatch[i]


  for i in range(len(falseBatch)-1, -1, -1):
    if "No Peptide Found" in falseBatch[i]:
      del falseBatch[i]


  if trueBatch == [] and falseBatch == []:
    return (matchedTrue,
            matchedFalse,
            pepMatchNumber)


  if trueBatch != []:

    trueBatch = [x.split(",") for x in trueBatch]

    trueBatch.sort(key=lambda x: float(x[5]), reverse=True)
    truePeptide = trueBatch[0][2]
    trueScore = float(trueBatch[0][5]) * \
                    (float(trueBatch[0][9]) / \
                     float(getPepLength(truePeptide)))
  else:
    trueScore = 0

  if falseBatch != []:

    falseBatch = [x.split(",") for x in falseBatch]

    falseBatch.sort(key=lambda x: float(x[5]), reverse=True)
    falsePeptide = falseBatch[0][2]

    # if trueScore == 0:
    #   trueScore = float(falseBatch[0][5]) * \
    #                 (float(falseBatch[0][9]) / \
    #                  float(getPepLength(falsePeptide)))

    if trueBatch != [] and checkMatch(truePeptide, falsePeptide):
      pepMatchNumber += 1
      matchedTrue.append(trueScore)

    else:
      if trueBatch != []:
        print(str(truePeptide) + " - " + str(falsePeptide))
        print()
      matchedFalse.append(trueScore)


  return(matchedTrue, matchedFalse, pepMatchNumber)




def plotBestGlobalScores(originalGlobalScore,
                         falseGlobalScore,
                         matchGlobalScore,
                         matchPct):

  plt.hist(originalGlobalScore,
           bins=50,
           histtype='step',
           color="blue",
           label="True Output")

  plt.hist(falseGlobalScore,
           bins=50, histtype='step',
           color="red",
           label="False Output")

  plt.hist(matchGlobalScore,
           bins=50,
           histtype='step',
           color="green",
           label="True/False Match")

  blue_patch = mpatches.Patch(color="blue", label="True Output")
  red_patch = mpatches.Patch(color="red", label="False Output")
  green_patch = mpatches.Patch(color="green", label="True/False Match")
  plt.legend(handles=[blue_patch, red_patch, green_patch])
  plt.xlabel("Global Score")
  plt.ylabel("Frequency")
  plt.title("Global Score Comparison\n(" + matchPct + "% match between true and false Peptides)")
  plt.show()


def plotBestAminoScores(originalAminoScore,
                        falseAminoScore,
                        matchAminoScore, 
                        matchPct):

  plt.hist(originalAminoScore, 
           bins=50,
           histtype='step',
           color="blue",
           label="True Output")

  plt.hist(falseAminoScore,
           bins=50,
           histtype='step',
           color="red",
           label="False Output")

  plt.hist(matchAminoScore,
           bins=50,
           histtype='step',
           color="green",
           label="True/False Match")

  blue_patch = mpatches.Patch(color="blue", label="True Output")
  red_patch = mpatches.Patch(color="red", label="False Output")
  green_patch = mpatches.Patch(color="green", label="True/False Match")
  plt.legend(handles=[blue_patch, red_patch, green_patch])
  plt.xlabel("Global Score")
  plt.ylabel("Frequency")
  plt.title("Amino Score Comparison\n(" + matchPct + "% match between true and false Peptides)")
  plt.show()


def plotCorrelation(originalAminoScore, falseAminoScore,
                    originalGlobalScore, falseGlobalScore):
  plt.xlabel("Amino Score")
  plt.ylabel("Global Score")
  plt.title("Global/Amino Score scatter plot")
  plt.plot(originalAminoScore, originalGlobalScore, "bo", markersize=2, marker="x")
  plt.plot(falseAminoScore, falseGlobalScore, "ro", markersize=2, marker="x")
  plt.show()

def findLowestFDR(trueScores, falseScores, cutoffFDR):
  trueScores.sort(reverse=True)
  falseScores.sort(reverse=True)
  cutoffList = []

  FDRList = []

  for numTrue, trueValue in zip(range(1, len(trueScores)+1), trueScores):
    for numFalse, falseValue in zip(range(len(falseScores)), falseScores):
      if falseValue < trueValue:
        FDRList.append((numFalse, numTrue, trueValue, numFalse/numTrue))
        if (numFalse / numTrue) < cutoffFDR:
          cutoffList.append((numTrue, numFalse, trueValue, (numFalse / numTrue)))
        break

  FDRList.sort(key=lambda x: x[3])
  
  for elm in FDRList:
    print("Num False: " + str(elm[0]))
    print("Num True: " + str(elm[1]))
    print("True Value: " + str(elm[2]))
    print("FDR: " + str(elm[3]))
    print()

  cutoffList.sort(key=lambda x: x[0], reverse=True)
  for val in cutoffList:
    print(str(val[0]) + ", " + str(val[1]) + ", " + str(val[2]) + ", " + str(val[3]))

def plotCombination(combinedTrue, combinedFalse, matched, cutoffFDR):

  print("True total: " + str(len(combinedTrue)))
  print("False total: " + str(len(combinedFalse)))
  print("Matched: " + str(matched))
  print()

  findLowestFDR(combinedTrue, combinedFalse, cutoffFDR)

  plt.hist(combinedTrue, 
           bins=50,
           histtype='step',
           color="blue",
           label="DeNovo Output")

  plt.hist(combinedFalse,
           bins=50,
           histtype='step',
           color="red",
           label="Random Output")

  blue_patch = mpatches.Patch(color="blue", label="DeNovo Output")
  red_patch = mpatches.Patch(color="red", label="Random Output")
  plt.legend(handles=[blue_patch, red_patch])
  plt.xlabel("Global Score")
  plt.ylabel("Frequency")
  plt.title("Global Score Plot")
  plt.show()

def analyzeResults(originalOutput, falseOutput):
  with open(originalOutput) as f:
    original = f.read()
  with open(falseOutput) as f:
    false = f.read()

  original = original.split("\n")
  false = false.split("\n")
  del original[0]
  if original[-1] == '':
    del original[-1]
  del false[0]
  if false[-1] == '':
    del false[-1]

  originalFileLine = 0
  falseFileLine = 0

  misGlobal = 0
  misAmino = 0
  totalGlobal = 0
  totalAmino = 0

  originalGlobalScore = []
  originalAminoScore = []

  falseGlobalScore = []
  falseAminoScore = []

  matchGlobalScore = []
  matchAminoScore = []

  combinedTrue = []
  combinedFalse = []
  matched = 0

  matchedTrue = []
  matchedFalse = []
  matchedMatched = 0

  newList = []

  totalTrue = 0
  totalFalse = 0

  uniquePeptides = set()
  pepSpecAboveFDR = ''
  specNums = {}

  iterations = 0

  while(True):
    originalBatch, originalFileLine = getNextBatch(original, originalFileLine)
    falseBatch, falseFileLine = getNextBatch(false, falseFileLine)
    if originalBatch == []:
      break

    specLine = originalBatch[0]
    if "No Peptide Found" not in specLine:
      specLine = specLine.split(',')
      specNums[specLine[1]] = 1

    iterations += 1

    # combinedTrue, combinedFalse, matched = \
    #     getBestCombineScores(originalBatch, falseBatch,
    #                          combinedTrue, combinedFalse,
    #                          matched)

    combinedTrue, combinedFalse, totalTrue, totalFalse = \
        getBestCombineScoresRandDecoy(originalBatch, falseBatch,
                                      combinedTrue, combinedFalse,
                                      totalTrue, totalFalse, uniquePeptides)

    # matchedTrue, matchedFalse, matchedMatched = \
    #     getMatchedScores(originalBatch, falseBatch,
    #                      matchedTrue, matchedFalse,
    #                      matchedMatched)

  #   (misGlobal, misAmino, totalGlobal, totalAmino) = \
  #     getBestScores(falseBatch,
  #                   originalBatch,
  #                   originalGlobalScore,
  #                   falseGlobalScore,
  #                   originalAminoScore,
  #                   falseAminoScore,
  #                   matchGlobalScore,
  #                   matchAminoScore,
  #                   misGlobal,
  #                   misAmino,
  #                   totalGlobal,
  #                   totalAmino,
  #                   newList,
  #                   0.15)
  
  # print(len(newList))
    
  # plotBestGlobalScores(originalGlobalScore,
  #                      falseGlobalScore,
  #                      matchGlobalScore,
  #                      str(int( ((totalGlobal - misGlobal) / totalGlobal) * 100)))
  # plotBestAminoScores(originalAminoScore,
  #                     falseAminoScore,
  #                     matchAminoScore,
  #                     str(int( ((totalAmino - misAmino) / totalAmino) * 100)))

  # plotCorrelation(originalAminoScore, falseAminoScore,
  #                 originalGlobalScore, falseGlobalScore)

  print(iterations)
  exit()

  plotCombination(combinedTrue, combinedFalse, matched, 0.01)
  #plotCombination(matchedTrue, matchedFalse, matchedMatched, 0.01)

  #print("Global score matches: " + str((totalGlobal - misGlobal) / totalGlobal) + "%")
  #print("Amino score matches: " + str((totalAmino - misAmino) / totalAmino) + "%")


if __name__ == "__main__":
  analyzeResults(argv[1], argv[2])