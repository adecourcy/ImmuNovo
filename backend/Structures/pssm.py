import copy

def getMatrixOfLength(title, length, allPSSM):
  return copy.deepcopy(allPSSM[title][1][length])

def getAcidProbabilities(matrix, acid):
  return copy.deepcopy(matrix[acid])

def getIgnoredLengths(title, allPSSM):
  return copy.deepcopy(allPSSM[title][0])

def _getDirection(title, allPSSM):
  return allPSSM[title][2]

def setDirection(forward, allPSSM):
  # newPSSM = {}
  # for title in allPSSM:
  #   if _getDirection(title, allPSSM) == forward:
  #     newPSSM[title] = copy.deepcopy(allPSSM[title])
  #   allLengths = {}
  #   for length in allPSSM[title][1]:
  #     matrix = getMatrixOfLength(title, length, allPSSM)
  #     for acid in matrix:
  #       matrix[acid].reverse()
  #     allLengths[length] = matrix
  #   newPSSM[title] = (allPSSM[0], allLengths, forward)
  
  def _reverseList(newMatrix, acidProbabilities, acid):
    newMatrix[acid] = acidProbabilities
    newMatrix[acid].reverse()
    return newMatrix

  for title in allPSSM:
    if _getDirection(title, allPSSM) == forward:
      return allPSSM
    else:
      break

  return _adjustPSSM(allPSSM, _reverseList)


def _adjustPSSM(allPSSM, adjustFunc):

  adjustedAllPSSM = {}

  for title in allPSSM:
    newPSSM = {}

    for length in allPSSM[title][1]:
      newMatrix = {}
      oldMatrix = getMatrixOfLength(title, length, allPSSM)

      for acid in oldMatrix:
        newMatrix = adjustFunc(newMatrix, getAcidProbabilities(oldMatrix, acid), acid)

      newPSSM[length] = newMatrix

    adjustedAllPSSM[title] = (allPSSM[title][0], newPSSM, allPSSM[title][2])

  return adjustedAllPSSM

def adjustForPrecision(allPSSM, precision):

  def _precAdjust(matrix, probabilities, acid):
    matrix[acid] = [int(10**precision * x) for x in probabilities]
    return matrix

  return _adjustPSSM(allPSSM, _precAdjust)

def convertAcids(allPSSM, acidConversion):

  def _acidConvert(matrix, probabilities, acid):
    matrix[acidConversion[acid]] = probabilities
    return matrix

  return _adjustPSSM(allPSSM, _acidConvert)