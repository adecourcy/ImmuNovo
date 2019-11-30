"""
This script finds all spectra with peptides scoring above a given FDR cutoff,
picks the top peptide for each of those spectra, and finds the unique peptides
from that subset. It separates the unique peptides by length and returns a
dictionary with lengths as keys (ints), and a general iteratable as contents.

This script changes all "I"s to "L"s and simply strips all AA modifications 
(i.e., M+15995 simply becomes M). For future expansion, it returns a, currently
unused, variable from the AA modification stripping process in case more
sophisticate behavior is to be implemented.

For input, this script requires a DeNovo output file, and an FDR cutoff
"""

import sys, os
sys.path.append(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])

from backend.constants import *


def removeModifications(peptideDF):
  removalDict = {} # placeholder for future functionality
  #print(peptideDF)
  removal = lambda row: ''.join([x for x in row[PEPTIDE] if x.isalpha()])
  
  peptideDF[PEPTIDE] = peptideDF.apply(removal, axis=1)

  return peptideDF, removalDict


def getTopPeptides(spectralGroups):
  topPeptides = set()
  for elm in spectralGroups:
    df = elm[1]
    df = df[df[SCORE_COMBINED] == df[SCORE_COMBINED].max()]
    # If equal, pick greastest PSSM Score
    df = [df[SCORE_PSSM] == df[SCORE_PSSM].max()]
    # If all equal, just pick the first one
    topPeptides.add(df[PEPTIDE][0])
  return topPeptides


def separateByLength(topPeptides):
  lengthDict = {}
  for pep in topPeptides:
    if len(pep) not in lengthDict:
      lengthDict[len(pep)] = []
    lengthDict[len(pep)].append(pep)
  
  return lengthDict


def getPeptideDict(peptideDF, fdrCutoff, hardFDR=True):
  peptideDF = \
    peptideDF[peptideDF[FDR] <= fdrCutoff]
  if not hardFDR and len(peptideDF) == 0:
    while len(peptideDF) == 0:
      fdrCutoff += 0.01
      peptideDF = \
        peptideDF[peptideDF[FDR] <= fdrCutoff]

  peptideDF[PEPTIDE] = \
      peptideDF.apply(lambda row: row[PEPTIDE].replace('I', 'L'), axis=1)
  peptideDF, removalDict = removeModifications(peptideDF)

  spectralGroups = peptideDF.groupby(by=TITLE_SPECTRUM)
  topPeptides = getTopPeptides(spectralGroups)

  return separateByLength(topPeptides), fdrCutoff