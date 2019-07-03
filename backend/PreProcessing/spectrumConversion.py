from decimal import Decimal
from math import log
from typing import *

def processSpectrum(spectrumMasses: List[Decimal],
                    spectrumScores: List[Decimal],
                    pepMass: Decimal,
                    H2OMass: Decimal,
                    protonMass: Decimal,
                    charge: int,
                    compressionRate: int) -> Tuple[List[Decimal], List[Decimal]]:

  spectrumMasses, spectrumScores = _filterSpectrum(spectrumMasses,
                                                  spectrumScores)
  
  spectrumMasses, spectrumScores = _eliminatePepMass(spectrumMasses,
                                                    spectrumScores,
                                                    pepMass,
                                                    H2OMass,
                                                    charge)


  if compressionRate !=0:
    spectrumScores = _compressScores(spectrumScores, compressionRate)


  spectrumMassesDouble = \
      _createDoubleChargeMassSpectrum(spectrumMasses,
                                     protonMass)

  spectrumMasses, spectrumMassesDouble, \
  spectrumScores, spectrumScoresDouble = \
      _addFinalMass(spectrumMasses,
                   spectrumMassesDouble,
                   spectrumScores,
                   pepMass,
                   charge,
                   protonMass,
                   H2OMass)

  return (spectrumMasses, spectrumMassesDouble,
          spectrumScores, spectrumScoresDouble)


def _filterSpectrum(spectrumMasses: List[Decimal],
                   spectrumScores: List[Decimal]) -> Tuple[List[Decimal], List[Decimal]]:

  for i in range(len(spectrumMasses)-1, -1, -1):
    if spectrumScores[i] == 0:
      del spectrumScores[i]
      del spectrumMasses[i]

  return (spectrumMasses, spectrumScores)


def _compressScores(spectrumScores: List[Decimal],
                   compressionRate: int) -> List[Decimal]:

  for i in range(0, len(spectrumScores)):
    spectrumScores[i] = Decimal(log(float(spectrumScores[i]), compressionRate))

  return spectrumScores


def _eliminatePepMass(spectrumMasses: List[Decimal],
                     spectrumScores: List[Decimal],
                     pepMass: Decimal,
                     H2OMass: Decimal,
                     charge: int) -> Tuple[List[Decimal], List[Decimal]]:

  increment = int(10 / charge)
  if charge == 2:
    increment = 5
  if charge == 3:
    increment = 3

  comparisonPepMass = int(pepMass * 10)
  eliminations = []

  for i in range(0, len(spectrumMasses)):
    comparisonSpectrumMass = int(spectrumMasses[i] * 10)
    if comparisonPepMass == comparisonSpectrumMass:

      if (comparisonPepMass - increment) == \
            int(spectrumMasses[i-1] * increment):
        eliminations.append(i-1)

      eliminations.append(i)
      comparisonPepMass += increment

      for k in range(i+1, len(spectrumMasses)):
        comparisonSpectrumMass = int(spectrumMasses[k] * 10)
        if comparisonPepMass == comparisonSpectrumMass:
          eliminations.append(k)
          comparisonPepMass += increment
        elif comparisonSpectrumMass < comparisonPepMass:
          eliminations.append(k)
        else:
          break

      break

  eliminations.reverse()

  for i in eliminations:
    del spectrumMasses[i]
    del spectrumScores[i]

  eliminations = []
  comparisonPepMass = int((pepMass - (H2OMass / charge)) * 10)

  for i in range(0, len(spectrumMasses)):
    comparisonSpectrumMass = int(spectrumMasses[i] * 10)
    if comparisonPepMass == comparisonSpectrumMass:

      eliminations.append(i)
      comparisonPepMass += increment

      for k in range(i+1, len(spectrumMasses)):
        comparisonSpectrumMass = int(spectrumMasses[k] * 10)
        if comparisonPepMass == comparisonSpectrumMass:
          eliminations.append(k)
          comparisonPepMass += increment
        elif comparisonSpectrumMass < comparisonPepMass:
          eliminations.append(k)
        else:
          break

      break

  eliminations.reverse()

  for i in eliminations:
    del spectrumMasses[i]
    del spectrumScores[i]

  return (spectrumMasses, spectrumScores)


def _createDoubleChargeMassSpectrum(spectrumMasses: List[Decimal],
                                   protonMass: Decimal) -> List[Decimal]:

  doubleChargedSpectrum = [2*mass - protonMass for mass in spectrumMasses]

  return doubleChargedSpectrum


def _addFinalMass(spectrumMasses: List[Decimal],
                 spectrumMassesDouble: List[Decimal],
                 spectrumScores: List[Decimal],
                 pepMass: Decimal,
                 charge: int,
                 protonMass: Decimal,
                 H2OMass: Decimal) -> Tuple[List[Decimal], List[Decimal]]:

  spectrumScoresDouble = [x for x in spectrumScores]
  pepMass = (pepMass * charge) - (protonMass * (charge - 1))

  splitIndex = 0
  for i in range(0, len(spectrumMasses)):
    if spectrumMasses[i] > pepMass:
      splitIndex = i
      break

  if splitIndex != 0:
    spectrumMasses = spectrumMasses[:splitIndex]
    spectrumScores = spectrumScores[:splitIndex]
    
  spectrumMasses.append(pepMass)
  spectrumScores.append(Decimal(1))

  splitIndex = 0
  for i in range(0, len(spectrumMassesDouble)):
    if spectrumMassesDouble[i] > pepMass:
      splitIndex = i
      break

  if splitIndex != 0:
    spectrumMassesDouble = spectrumMassesDouble[:splitIndex]
    spectrumScoresDouble = spectrumScoresDouble[:splitIndex]

  spectrumMassesDouble.append(pepMass)
  spectrumScoresDouble.append(1)

  return (spectrumMasses, spectrumMassesDouble, spectrumScores, spectrumScoresDouble)


def adjustSpectrumPrecision(spectrumMasses: List[Decimal],
                            spectrumMassesDouble: List[Decimal],
                            spectrumScores: List[Decimal],
                            spectrumScoresDouble: List[Decimal],
                            precision):

  spectrumMasses = [int(x * (10**precision)) for x in spectrumMasses]
  spectrumMassesDouble = [int(x * (10**precision)) for x in spectrumMassesDouble]
  spectrumScores = [int(x * (10**precision)) for x in spectrumScores]
  spectrumScoresDouble = [int(x * (10**precision)) for x in spectrumScoresDouble]

  return spectrumMasses, spectrumMassesDouble, \
          spectrumScores, spectrumScoresDouble

def mergeMassesScores(spec, specDouble, score, scoreDouble):
  specTup = []
  specDoubleTup = []
  for mass, index in zip(specDouble, range(len(specDouble))):
    if mass not in spec:
      specDoubleTup.append((mass, scoreDouble[index]))

  for m, s in zip(spec, score):
    specTup.append((m, s))

  merged = specDoubleTup + specTup

  merged.sort(key=lambda x: x[0])

  eSpec = [x[0] for x in merged]
  eScore = [x[1] for x in merged]

  return eSpec, eScore