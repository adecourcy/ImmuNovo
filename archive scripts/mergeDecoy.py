if __name__ == "__main__":
  with open("decoyPeps501FromResultsMerged.csv") as f:
    newDecoys = f.read()
  newDecoys = newDecoys.split("\n")
  with open("decoyPeps501RandMerged.csv") as f:
    oldDecoys = f.read()
  oldDecoys = oldDecoys.split("\n")

  output = open("decoyPeps501RandAndResults.csv", 'w')

  oldDict = {}
  newDict = {}

  for lineOld, lineNew in zip(oldDecoys, newDecoys):
    splitLineOld = lineOld.split(",")
    splitLineNew = lineNew.split(",")

    oldDict[splitLineOld[0]] = lineOld
    newDict[splitLineNew[0]] = lineNew


  for spectrum in oldDict:
    splitLineOld = oldDict[spectrum].split(",")
    splitLineNew = newDict[spectrum].split(",")

    if splitLineNew[1] == "No Peptide Found":
      output.write(oldDict[spectrum] + '\n')
      continue
    elif splitLineOld[1] == "No Peptide Found":
      output.write(newDict[spectrum] + '\n')
      continue

    if float(splitLineOld[2]) > float(splitLineNew[2]):
      output.write(oldDict[spectrum] + '\n')
    else:
      output.write(newDict[spectrum] + '\n')
