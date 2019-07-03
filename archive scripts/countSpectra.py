def getCurrentSpectrum(file, fileLine):
  if fileLine >= len(file):
    return ''
  splitLine = file[fileLine].split(",")
  return splitLine[1]


def getNextBatch(file, fileLine, spectraDict, duplicateCount):

  currentSpectrum = getCurrentSpectrum(file, fileLine)
  if currentSpectrum not in spectraDict:
    spectraDict[currentSpectrum] = 0
  else:
    print(currentSpectrum)
    duplicateCount += 1


  results = []

  while(True):
    if fileLine >= len(file):
      break
    if currentSpectrum == getCurrentSpectrum(file, fileLine):
      results.append(file[fileLine])
      fileLine += 1
    else:
      break

  return results, fileLine, spectraDict, duplicateCount


if __name__ == "__main__":
  with open("C501.newMatrix.modification.output.tsv") as f:
    original = f.read()

  original = original.split("\n")
  del original[0]
  if original[-1] == '':
    del original[-1]

  fileLine = 0
  duplicateCount = 0
  spectraDict = {}
  results = []

  while(True):
    results, fileLine, spectraDict, duplicateCount = \
          getNextBatch(original, fileLine, spectraDict, duplicateCount)
    if results == []:
      break

  print(duplicateCount)


