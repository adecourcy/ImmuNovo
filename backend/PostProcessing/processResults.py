from typing import *
from decimal import Decimal
from math import sqrt

class Node:
  def __init__(self,
               peptideString: str,
               mass: Decimal,
               precursorError: Decimal,
               adjustedScore: Decimal,
               unadjustedScore: Decimal,
               massScore: Decimal,
               aminoScore: Decimal,
               globalScore,
               combinedScore):
    self.peptideString = peptideString
    self.mass = mass
    self.precursorError = precursorError
    self.adjustedScore = adjustedScore
    self.unadjustedScore = unadjustedScore
    self.massScore = massScore
    self.aminoScore = aminoScore
    self.globalScore = globalScore
    self.combinedScore = combinedScore

  def __eq__(self, other) -> bool:
    if type(other) is type(self):
      if self.combinedScore == other.combinedScore:
        return self.combinedScore == other.combinedScore
      else:
        return False

  def __gt__(self, other) -> bool:
    if type(other) is type(self):
      if self.combinedScore == other.combinedScore:
        return self.combinedScore > other.combinedScore
      else:
        return self.combinedScore > other.combinedScore

  def __lt__(self, other) -> bool:
    if type(other) is type(self):
      if self.combinedScore == other.combinedScore:
        return self.combinedScore < other.combinedScore
      else:
        return self.combinedScore < other.combinedScore

def deConvertPeptideString(peptideString, acidConversion):

  reverseConversionTable = {acidConversion[x]: x for x in acidConversion}

  convertedPeptide = ''
  for acid in peptideString:
    convertedPeptide += reverseConversionTable[acid]
  
  return convertedPeptide



def processResults(resultsFile: str,
                   acidMassTable: Dict[str, int],
                   experimentalSpectrum: List[int],
                   experimentalScores: List[int],
                   protonMassModified: int,
                   H2OMassModified: int,
                   NH3MassModified: int,
                   maxMassTolerance: int,
                   precision: int,
                   acidConversion: List[Tuple[str, str]]) -> List[Node]:

  decPrec = Decimal(10 ** precision)
  results = open(resultsFile, "r")
  nodes = []

  while True:
    peptideString = results.readline().strip()
    if peptideString == "":
      break

    mass = Decimal(results.readline()) / decPrec
    precursorMass = results.readline()
    adjustedScore = Decimal(results.readline()) / decPrec
    unadjustedScore = Decimal(results.readline()) / decPrec
    massScore = Decimal(results.readline()) / decPrec
    aminoScore = (Decimal(results.readline()) / decPrec) / len(peptideString)
    globalScore  = calculateGlobalScore(acidMassTable,
                                        experimentalSpectrum,
                                        experimentalScores,
                                        peptideString[::-1],
                                        protonMassModified,
                                        H2OMassModified,
                                        NH3MassModified,
                                        maxMassTolerance)

    combinedScore = aminoScore * Decimal(globalScore)
      
    nodes.append(Node(deConvertPeptideString(peptideString[::-1], acidConversion),
                      mass,
                      precursorMass,
                      adjustedScore,
                      unadjustedScore,
                      massScore,
                      aminoScore,
                      globalScore,
                      combinedScore))

  nodes.sort()

  return nodes[::-1]


# Taken directly from Lei's code
def cosineSimilarity(spectrum_lhs: List[Tuple[int, int]],
                     spectrum_rhs: List[Tuple[int, int]],
                     maxMassTolerance) -> float:

  similarity = 0
  magnitude_lhs = 0
  magnitude_rhs = 0

  i_lhs = 0
  i_rhs = 0

  while (i_lhs < len(spectrum_lhs) and i_rhs < len(spectrum_rhs)):

    mz_lhs = spectrum_lhs[i_lhs][0]
    mz_rhs = spectrum_rhs[i_rhs][0]

    intensity_lhs = spectrum_lhs[i_lhs][1]
    intensity_rhs = spectrum_rhs[i_rhs][1]

    if massTolerance(mz_lhs, mz_rhs) < maxMassTolerance:
      similarity += intensity_lhs * intensity_rhs
      magnitude_lhs += intensity_lhs * intensity_lhs
      magnitude_rhs += intensity_rhs * intensity_rhs
      i_lhs += 1
      i_rhs += 1
    elif (mz_lhs < mz_rhs):
      magnitude_lhs += intensity_lhs * intensity_lhs
      i_lhs += 1
    else:
      magnitude_rhs += intensity_rhs * intensity_rhs
      i_rhs += 1
  
  while (i_lhs < len(spectrum_lhs)):
    magnitude_lhs += spectrum_lhs[i_lhs][1] ** 2
    i_lhs += 1
  
  while (i_rhs < len(spectrum_rhs)):
    magnitude_rhs += spectrum_rhs[i_rhs][1] ** 2
    i_rhs += 1
  
  similarity /= (sqrt(magnitude_rhs) * sqrt(magnitude_lhs))

  return round(similarity, 4)


def calculateGlobalScore(acidMassTable: Dict[str, int],
                         experimentalSpectrum: List[int],
                         experimentalIntensities: List[int],
                         peptideString: str,
                         protonMassModified: int,
                         H2OMassModified: int,
                         NH3MassModified: int,
                         maxMassTolerance: int) -> float:

  theoreticalSpectrum, yEnd, bEnd = \
      generateTheoreticalSpectrum(acidMassTable,
                                  peptideString,
                                  protonMassModified,
                                  H2OMassModified,
                                  NH3MassModified,
                                  peptideString)

  experimentalVector, experimentalIntenseVector = \
        createExperimentalVector(theoreticalSpectrum,
                                 experimentalSpectrum,
                                 experimentalIntensities,
                                 maxMassTolerance)
  
  theoreticalMassToleranceDifferences = []
  for i in range(len(theoreticalSpectrum)):
    if experimentalVector[i] == 0:
      theoreticalMassToleranceDifferences.append(maxMassTolerance)
    else:
      massToleranceResult = \
        massTolerance(theoreticalSpectrum[i], experimentalVector[i])
      theoreticalMassToleranceDifferences.append(massToleranceResult)
  

  sumIntensities = sum(experimentalIntenseVector)
  if sumIntensities == 0:
    return 0
  spectrumIntensitiesNormalized = \
    [x / sumIntensities for x in experimentalIntenseVector]


  euclideanDistance = 0

  for mt, weight in zip(theoreticalMassToleranceDifferences,
                        spectrumIntensitiesNormalized):
    euclideanDistance += (1 - mt/maxMassTolerance) * weight

  return euclideanDistance


def createExperimentalVector(theoreticalSpectrum: List[int],
                             experimentalSpectrum: List[int],
                             experimentalIntensities: List[int],
                             maxMassTolerance: int) -> List[int]:

  theoreticalSpectrum.sort()

  expSpecIntense = []
  for x,y in zip(experimentalSpectrum, experimentalIntensities):
    expSpecIntense.append((x, y))

  expSpecIntense.sort(key=lambda x: x[0])
  
  experimentalVector = []
  experimentalIntenseVector = []

  i = 0
  k = 0
  while (i < len(theoreticalSpectrum)) and (k < len(expSpecIntense)):

    if massTolerance(theoreticalSpectrum[i], expSpecIntense[k][0]) \
        < maxMassTolerance:

      experimentalVector.append(expSpecIntense[k][0])
      experimentalIntenseVector.append(expSpecIntense[k][1])
      i += 1
      k += 1

    elif expSpecIntense[k][0] > theoreticalSpectrum[i]:
      i += 1
      experimentalVector.append(0)
      experimentalIntenseVector.append(0)

    else:
      k += 1

  while len(experimentalVector) < len(theoreticalSpectrum):
    experimentalVector.append(0)
    experimentalIntenseVector.append(0)

  return experimentalVector, experimentalIntenseVector


def massTolerance(calculatedMass: int,
                  experimentalMass: int) -> float:
  return abs((((calculatedMass - experimentalMass)/experimentalMass) \
              * 1000000))


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


def generateTheoreticalSpectrum(acidMassTable: Dict[str, int],
                                peptideString: str,
                                protonMassModified: int,
                                H2OMassModified: int,
                                NH3MassModified: int,
                                peptideCharge: int) -> List[int]:

  bEnd = generateSpectrumList(acidMassTable,
                              peptideString,
                              protonMassModified,
                              H2OMassModified,
                              NH3MassModified,
                              peptideCharge,
                              False)

  yEnd = generateSpectrumList(acidMassTable,
                              peptideString[::-1],
                              protonMassModified,
                              H2OMassModified,
                              NH3MassModified,
                              peptideCharge,
                              True)

  spectrum = yEnd + bEnd
  final = protonMassModified
  for i in range(0, len(peptideString)):
    final += acidMassTable[peptideString[i]]

  spectrum.append(final)


  return (spectrum, yEnd, bEnd)


def generateSpectrumList(aminoMassesModified: Dict[str, int],
                         peptideString: str,
                         protonMassModified: int,
                         H2OMassModified: int,
                         NH3MassModified: int,
                         peptideCharge: int,
                         yEnd: bool) -> List[int]:

  if yEnd:
    additiveMass = H2OMassModified + protonMassModified
  else:
    additiveMass = protonMassModified

  spectrum = [aminoMassesModified[peptideString[0]] + additiveMass]

  for i in range(1, len(peptideString)-1):
    spectrum.append(spectrum[-1] + aminoMassesModified[peptideString[i]])

  spectrumMinusNH3 = [x - NH3MassModified for x in spectrum]
  spectrumMinusH2O = [x - H2OMassModified for x in spectrum]

  spectrum += spectrumMinusNH3
  spectrum += spectrumMinusH2O

  return spectrum
