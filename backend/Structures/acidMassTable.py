def getMass(masses, acid):
  return masses[acid]

def getAcids(table):
  acids = []
  for key in table:
    acids.append(key)
  return acids

def adjustForPrecision(masses, precision):
  newMasses = {}

  for acid in masses:
    newMasses[acid] = int(masses[acid] * 10**precision)

  return newMasses

def convertAcids(masses, acidConversion):
  newMasses = {}

  for acid in masses:
    newMasses[acidConversion[acid]] = int(masses[acid])

  return newMasses


def peptideMass(peptide,
                acidMassTable):
  mass = 0
  for acid in peptide:
    mass += getMass(acidMassTable, acid)
  return mass