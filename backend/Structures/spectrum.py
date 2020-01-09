def getPrecursorMass(spectrum):
  return spectrum['pepMass']

def getCharge(spectrum):
  return spectrum['charge']

def getTitle(spectrum):
  return spectrum['title']

def getMasses(spectrum):
  return spectrum['masses']

def getIntensities(spectrum):
  return spectrum['intensities']

# def adjustForPrecision(spectrum, precision):

#   newSpectrum = {'title': spectrum['title'],
#                  'pepMass': int(spectrum['pepMass'] * precision),
#                  'charge': spectrum['charge'],
#                  'masses': [int(precision * x) for x in spectrum['masses']],
#                  'intensities': [int(precision * x) for x in spectrum['intensities']]}

#   return newSpectrum