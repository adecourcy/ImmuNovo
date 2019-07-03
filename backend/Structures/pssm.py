def getMatrixOfLength(title, length, allPSSM):
  return allPSSM[title][1][length]

def getAcidProbabilities(matrix, acid):
  return matrix[acid]

def getIgnoredLengths(title, allPSSM):
  return allPSSM[title][0]

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

    adjustedAllPSSM[title] = (allPSSM[title][0], newPSSM)

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