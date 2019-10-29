from decimal import Decimal
from os import listdir, path

def getAllPSSM(directory, minPepLength, maxPepLength, acids):
  allPSSM = {}
  missing = {}

  for file in listdir(directory):
    missing = {}
    pssm = _parsePSSM(path.join(directory, file))

    missing = _assessPSSM(missing,
                          pssm,
                          maxPepLength,
                          minPepLength)

    pssm = _adjustPSSM(pssm,
                       minPepLength,
                       maxPepLength,
                       acids)

    if missing == [] or _promptUser(missing, file) == 'uniform':
      allPSSM[file] = ([], pssm)
    else:
      allPSSM[file] = (missing, pssm)

  allPSSM['uniform'] = ([], _adjustPSSM({}, minPepLength, maxPepLength, acids))
  print(allPSSM['uniform'])

  return allPSSM


def _parsePSSM(file):
  # A pssm will be a dictionary of dictionaries of lists. Explained thusly:
  #
  # pssm[length] will return a pssm for all amino acids of length 'length'
  # pssm[length][acid] will return a list
  # of probabilities for acid 'acid' of length 'length'

  try:
    pssmFile = open(file, 'r')
  except FileNotFoundError:
    print('unable to open pssm file: {}'.format(file))
    exit()

  pssm = {}


  for line in pssmFile:
    if "START" in line:
      length = int(line.split()[1])
      pssm[length] = {}

    else:
      lineList = line.split()
      acid = lineList[0]
      # We're just always assuming a y-end peptide for now
      scores = lineList[1:]
      scores.reverse()

      if len(scores) != length:
        print(str.format("Table length for '{}' in file {} should be {} but is {}",
                         acid,
                         file,
                         length,
                         len(scores)))
        exit()

      pssm[length][acid] = [Decimal(x) for x in scores]

  return pssm


def _assessPSSM(missing,
                pssm,
                maxPepLength,
                minPepLength):

  missing = []

  for i in range(minPepLength, maxPepLength+1):
    if i not in pssm:
      missing.append(i)

  return missing


def _promptUser(missing, title):
  return 'ignore'

  print(str.format("For {} there is no PSSM information for peptides of length {}",
                   title,
                   " ".join([str(x) for x in missing])))

  while True:
    print("Would you like to:")
    response = input(str.format("(e)xit, (u)se a uniform PSSM, or "
                                "(d)iscard peptides of this length for {}? ",
                                title))
    if response.lower() == 'e':
      exit()
    elif response.lower() == 'u':
      return 'uniform'
    elif response.lower() == 'd':
      return 'ignore'
    else:
      print("Please choose a valid response")


def _adjustPSSM(pssm, minPepLength, maxPepLength, aminoAcids):
  
  newPSSM = {}

  for i in range(minPepLength, maxPepLength+1):
    if i in pssm:
      newPSSM[i] = pssm[i]
    else:
      newPSSM[i] = _createUniformMatrix(pssm, i, aminoAcids)

  return newPSSM

def _createUniformMatrix(pssm, length, aminoAcids):
  newMatrix = {}

  numAcids = len(aminoAcids)
  for acid in aminoAcids:
    newMatrix[acid] = []
    for i in range(0, length):
      newMatrix[acid].append(1/numAcids)

  return newMatrix