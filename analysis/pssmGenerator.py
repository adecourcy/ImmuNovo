#!/usr/bin/env python3
import sys


def parseFile(file, fileName, minPeptides):
  allPeps = {}

  encountered = set()
  unknownAminos = set()

  for line in file:
    line = line.strip()
    line = line.split(',')
    peptide = line[2].strip('"')
    peptide = peptide.split()
    peptide = peptide[0]

    if peptide.isalpha() and peptide.isupper():

      if len(peptide) not in allPeps:
        allPeps[len(peptide)] = {peptide}
      else:
        allPeps[len(peptide)].add(peptide)

  outputName = fileName + ".pssm"
  reportName = fileName + ".report"

  output = open(outputName, "w")

  allLengths = []
  for length in allPeps:
    allLengths.append(length)
  allLengths.sort()

  for length in allLengths:
    if len(allPeps[length]) < minPeptides:
      continue
    else:
      output.write("START {}\n".format(length))
      output.write(createTable(list(allPeps[length]), unknownAminos))

  reportName = open(reportName, "w")
  reportName.write("Unique peptides found: " + str(len(encountered)) + "\n\n")
  for length in allLengths:
    reportName.write("Peptides length {}: {}\n".format(length, len(allPeps[length])))

  reportName.write("Peptides with unknown amino acids:\n")
  for pep in unknownAminos:
      reportName.write(pep + "\n")


def createTable(peptides, unknownAminos):

  table = []
  acids = ["G",
           "A",
           "S",
           "P",
           "V",
           "T",
           "C",
           "C+57.021",
           "L",
           "I",
           "N",
           "D",
           "Q",
           "K",
           "E",
           "M",
           "M+15.995",
           "H",
           "F",
           "R",
           "Y",
           "W"]


  for i in range(0, len(peptides[0])):
    column = {}

    for acid in acids:
      column[acid] = 1
 
    table.append(column)


  for peptide in peptides:
    for i in range(0, len(peptide)):
      if peptide[i] not in table[i]:
        unknownAminos.add(peptide)
        continue
      else:
        table[i][peptide[i]] += 1
        if peptide[i] == "M":
            table[i]["M+15.995"] += 1
        if peptide[i] == "C":
            table[i]["C+57.021"] += 1


  for column in table:
    total = 0

    for key in column:
      total += column[key]

    for key in column:
      column[key] = column[key] / total

  output = ""

  for acid in acids:
    output += acid + " "

    for column in table:
      output += str(column[acid]) + " "

    output = output.strip()
    output += "\n"



  return output





if __name__ == "__main__":

  start = 1
  minPeptides = 100
  if sys.argv[1].isdigit():
    minPeptides = int(sys.argv[1])
    start = 2

  for i in range(start, len(sys.argv)):
    try:
      file = open(sys.argv[i], "r")
    except FileNotFoundError:
      print("Could not find file " + sys.argv[i])
      continue

    parseFile(file, sys.argv[i], minPeptides)