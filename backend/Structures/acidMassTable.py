def getMass(masses, acid):
  return masses[acid]

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