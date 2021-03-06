
def getSpectrums(file):


  try:
    spectrumFile = open(file, 'r')
  except FileNotFoundError:
    print("Unable to open spectrum file: {}".format(file))
    exit()
  spectrumFile.close()

  lastLocation = 0

  while True:
    spectrumFile = open(file, 'r')
    spectrumFile.seek(lastLocation)
    nextSpectrum = _parseSpectrumFile(spectrumFile)
    lastLocation = spectrumFile.tell()
    spectrumFile.close()
    if nextSpectrum == '':
      break
    yield nextSpectrum


def _parseSpectrumFile(spectrumFile):

  spectrum = {'title': '',
              'pepMass': '',
              'charge': 2,
              'masses': [],
              'intensities': []}

  nextLine = spectrumFile.readline()
  while nextLine.isspace():
    nextLine = spectrumFile.readline()
  if nextLine == '':
    return ''

  while "END IONS" not in nextLine:
    if '###' in nextLine:
      nextLine = spectrumFile.readline()
      continue
    if 'PEPMASS=' in nextLine:
      spectrum['pepMass'] = float(nextLine.split('=')[1].split()[0])
    elif 'CHARGE=' in nextLine:
      spectrum['charge'] = int(nextLine.split('=')[1][0])
    elif 'TITLE=' in nextLine:
      spectrum['title'] = nextLine.split('=')[1].strip().replace(",", '').split()[0]
    elif nextLine[0].isdigit():
      mass, intensity = nextLine.split()
      spectrum['masses'].append(float(mass))
      spectrum['intensities'].append(float(intensity))
    nextLine = spectrumFile.readline()

  return spectrum