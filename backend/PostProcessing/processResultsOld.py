from typing import *
from math import sqrt

# cheap hack until I can figure out how to do this properly in
# python 3.6
sys.path.append(os.path.split(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])[0])

import backend.PostProcessing.processResults as ProcessResults

class Node:
  def __init__(self,
               peptideString: str,
               mass: float,
               precursorError: float,
               adjustedScore: float,
               totalScore: float,
               massScore: float,
               aminoScore: float,
               globalScore):
    self.peptideString = peptideString
    self.mass = mass
    self.precursorError = precursorError
    self.adjustedScore = adjustedScore
    self.totalScore = totalScore
    self.massScore = massScore
    self.aminoScore = aminoScore
    self.globalScore = globalScore

  def __eq__(self, other) -> bool:
    if type(other) is type(self):
      if self.globalScore == other.globalScore:
        return self.adjustedScore == other.adjustedScore
      else:
        return False

  def __gt__(self, other) -> bool:
    if type(other) is type(self):
      if self.globalScore == other.globalScore:
        return self.adjustedScore > other.adjustedScore
      else:
        return self.globalScore > other.globalScore

  def __lt__(self, other) -> bool:
    if type(other) is type(self):
      if self.globalScore == other.globalScore:
        return self.adjustedScore < other.adjustedScore
      else:
        return self.globalScore < other.globalScore

def deConvertPeptideString(peptideString, acidConversion):
  standardAA = ["G", "A", "S", "P", "V",
                "T", "L", "I", "N", "D",
                "Q", "K", "E", "M", "H",
                "F", "R", "C", "Y", "W"]

  convertedPeptide = ''
  for acid in peptideString:
    if acid in standardAA:
      convertedPeptide += acid
    else:
      for tup in acidConversion:
        if acid == tup[1]:
          convertedPeptide += tup[0]
          break
  
  return convertedPeptide



def processResults(resultsFile: str,
                   aminoMassesModified: Dict[str, int],
                   experimentalSpectrum: List[int],
                   experimentalScores: List[int],
                   bPenalty: float,
                   protonMassModified: int,
                   H2OMassModified: int,
                   NH3MassModified: int,
                   maxMassTolerance: int,
                   precision: int,
                   acidConversion: List[Tuple[str, str]]) -> List[Node]:

  decPrec = float(10 ** precision)
  results = open(resultsFile, "r")
  nodes = []

  while True:
    peptideString = results.readline().strip()
    if peptideString == "":
      break
    mass = results.readline()
    precursorMass = results.readline()
    adjustedScore = results.readline()
    totalScore = results.readline()
    massScore = results.readline()
    aminoScore = results.readline()

    globalScore  = calculateGlobalScore(aminoMassesModified,
                                        experimentalSpectrum,
                                        experimentalScores,
                                        bPenalty,
                                        peptideString[::-1],
                                        protonMassModified,
                                        H2OMassModified,
                                        NH3MassModified,
                                        maxMassTolerance)
      
    nodes.append(Node(deConvertPeptideString(peptideString[::-1], acidConversion),
                      float(mass) / decPrec,
                      float(precursorMass),
                      float(adjustedScore) / decPrec,
                      float(totalScore) / decPrec,
                      float(massScore) / decPrec,
                      float(aminoScore) / decPrec,
                      globalScore))

  nodes.sort()

  return nodes[::-1]


def calculateGlobalScore(aminoMassesModified: Dict[str, int],
                         experimentalSpectrum: List[int],
                         experimentalScores: List[int],
                         bPenalty: float,
                         peptideString: str,
                         protonMassModified: int,
                         H2OMassModified: int,
                         NH3MassModified: int,
                         maxMassTolerance: int) -> float:

  theoreticalSpectrum, yEnd, bEnd = \
      ProcessResults.generateTheoreticalSpectrum(aminoMassesModified,
                                                 peptideString,
                                                 protonMassModified,
                                                 H2OMassModified,
                                                 NH3MassModified,
                                                 peptideString)

  theoreticalVector = createTheoreticalVector(theoreticalSpectrum,
                                              experimentalSpectrum,
                                              maxMassTolerance)
  
  theoreticalMassToleranceDifferences = []
  for i in range(len(experimentalSpectrum)):
    if theoreticalVector[i] == 0:
      theoreticalMassToleranceDifferences.append(maxMassTolerance)
    else:
      massToleranceResult = \
        massTolerance(theoreticalVector[i], experimentalSpectrum[i])
      theoreticalMassToleranceDifferences.append(massToleranceResult)
  

  sumSpectrumScores = sum(experimentalScores)
  if sumSpectrumScores == 0:
    sumSpectrumScores = 1
  spectrumScoresNormalized = \
    [x / sumSpectrumScores for x in experimentalScores]


  euclideanDistance = 0


  for mt, weight in zip(theoreticalMassToleranceDifferences,
                        spectrumScoresNormalized):
    euclideanDistance += mt * weight


  return (1 - (euclideanDistance / maxMassTolerance))


def createTheoreticalVector(theoreticalSpectrum: List[int],
                            experimentalSpectrum: List[int],
                            maxMassTolerance: int) -> List[int]:

  theoreticalSpectrum.sort()
  experimentalSpectrum.sort()
  
  theoreticalVector = []

  i = 0
  k = 0
  while (i < len(theoreticalSpectrum)) and (k < len(experimentalSpectrum)):

    if massTolerance(theoreticalSpectrum[i], experimentalSpectrum[k]) \
        < maxMassTolerance:

      theoreticalVector.append(theoreticalSpectrum[i])
      i += 1
      k += 1

    elif theoreticalSpectrum[i] > experimentalSpectrum[k]:
      k += 1
      theoreticalVector.append(0)

    else:
      i += 1

  while len(theoreticalVector) < len(experimentalSpectrum):
    theoreticalVector.append(0)

  return theoreticalVector


def massTolerance(calculatedMass: int,
                  experimentalMass: int) -> float:
  return abs((((calculatedMass - experimentalMass)/experimentalMass) \
              * 1000000))


# def generateTheoreticalSpectrum(aminoMassesModified: Dict[str, int],
#                                 peptideString: str,
#                                 protonMassModified: int,
#                                 H2OMassModified: int,
#                                 NH3MassModified: int,
#                                 peptideCharge: int) -> List[int]:

#   bEnd = generateSpectrumList(aminoMassesModified,
#                               peptideString,
#                               protonMassModified,
#                               H2OMassModified,
#                               NH3MassModified,
#                               peptideCharge,
#                               False)

#   yEnd = generateSpectrumList(aminoMassesModified,
#                               peptideString[::-1],
#                               protonMassModified,
#                               H2OMassModified,
#                               NH3MassModified,
#                               peptideCharge,
#                               True)

#   spectrum = yEnd + bEnd
#   final = protonMassModified
#   for i in range(0, len(peptideString)):
#     final += aminoMassesModified[peptideString[i]]

#   spectrum.append(final)


#   return (spectrum, yEnd, bEnd)


# def generateSpectrumList(aminoMassesModified: Dict[str, int],
#                          peptideString: str,
#                          protonMassModified: int,
#                          H2OMassModified: int,
#                          NH3MassModified: int,
#                          peptideCharge: int,
#                          yEnd: bool) -> List[int]:

#   if yEnd:
#     additiveMass = H2OMassModified + protonMassModified
#   else:
#     additiveMass = protonMassModified

#   spectrum = [aminoMassesModified[peptideString[0]] + additiveMass]

#   for i in range(1, len(peptideString)-1):
#     spectrum.append(spectrum[-1] + aminoMassesModified[peptideString[i]])

#   spectrumMinusNH3 = [x - NH3MassModified for x in spectrum]
#   spectrumMinusH2O = [x - H2OMassModified for x in spectrum]

#   spectrum += spectrumMinusNH3
#   spectrum += spectrumMinusH2O

#   return spectrum
