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