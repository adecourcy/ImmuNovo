#!/usr/bin/env python3
from random import shuffle
from decimal import Decimal
import sys, os
from copy import deepcopy

def decoyMatrixString(matrix):
  aminos = []
  probs = []

  for entry in matrix:
    entry = entry.split()
    aminos.append(entry[0])
    probs.append([Decimal(x) for x in entry[1:]])

  columnMatrix = []
  for col in range(len(probs[0])):
    colList = []
    for row in range(len(probs)):
      colList.append(probs[row][col])
    columnMatrix.append(colList)

  oldCol = deepcopy(columnMatrix)
  while columnMatrix == oldCol:
   shuffle(columnMatrix)

  rowMatrix = []
  for col in range(len(columnMatrix[0])):
    colList = []
    for row in range(len(columnMatrix)):
      colList.append(columnMatrix[row][col])
    rowMatrix.append(colList)

  outputString = ''
  for amino, prob in zip(aminos, rowMatrix):
    outputString += "{} {}\n".format(amino, ' '.join([str(x) for x in prob]))
  return outputString


def createDecoyPSSM(pssmFile, pssmFileDecoy):
  matrix = []

  for line in pssmFile:
    if 'START' in line and matrix == []:
      pssmFileDecoy.write(line)
    elif 'START' in line:
      pssmFileDecoy.write(decoyMatrixString(matrix))
      matrix = []
      pssmFileDecoy.write(line)
    elif line == '':
      break
    else:
      matrix.append(line.strip())


if __name__ == "__main__":

  for file in sys.argv[1:]:
    try:
      pssmFile = open(file, 'r')
    except FileNotFoundError:
      print('cannot open {}, skipping...'.format(file))
      continue

    try:
      pssmFileDecoy = open(file + ".decoy", 'w')
    except FileNotFoundError:
      print('cannot open {} to write decoy file, skipping...'.format(file + '.decoy'))

    createDecoyPSSM(pssmFile, pssmFileDecoy)

