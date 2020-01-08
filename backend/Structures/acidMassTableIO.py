
def getAminoMasses(file):
  try:
    aminoMassFile = open(file, 'r')
  except FileNotFoundError:
    print("Unable to open amino acide mass file: {}".format(file))
    exit()

  aminos = {}

  for line in aminoMassFile:
    amino, mass = line.split()
    aminos[amino] = float(mass)
  return aminos