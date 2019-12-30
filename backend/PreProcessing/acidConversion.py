import backend.Structures.acidMassTable as acidMassTable
import backend.Structures.pssm as pssm

def createAcidConversionTable(acidMasses):
  standardAA = ["G", "A", "S", "P", "V",
                "T", "L", "I", "N", "D",
                "Q", "K", "E", "M", "H",
                "F", "R", "C", "Y", "W"]
  
  modifiedAA = ["B", "J", "O"]

  conversionList = {}

  if len(acidMasses) > 23:
    return {}
  
  else:
    modIndex = 0
    for acid in acidMasses:
      if acid not in standardAA:
        conversionList[acid] = modifiedAA[modIndex]
        modIndex += 1
      else:
        conversionList[acid] = acid
  
  return conversionList

def convertAcidModifications(acidMasses, allPSSM, conversionTable):

  newPSSM = pssm.convertAcids(allPSSM, conversionTable)
  newAcidMasses = acidMassTable.convertAcids(acidMasses, conversionTable)

  return newAcidMasses, newPSSM


def convertPeptideString(peptideString, acidConversion):
  allConversions = []
  for key in acidConversion:
    if key != acidConversion[key]:
      allConversions.append(key)
  allConversions.sort(key=lambda x: len(x), reverse=True)
  print(allConversions)

  for conversion in allConversions:
    peptideString = peptideString.replace(conversion, acidConversion[conversion])

  return peptideString


def deConvertPeptideString(peptideString, acidConversion):

  reverseConversionTable = {acidConversion[x]: x for x in acidConversion}

  convertedPeptide = ''
  for acid in peptideString:
    convertedPeptide += reverseConversionTable[acid]
  
  return convertedPeptide