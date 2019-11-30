"""
Run a series of programs in the analysis directory and create a results report.

This program takes a decoy databases, the results from database search, and
the results of ImmuNovo, filters all ImmuNovo and database search results below
a certain FDR threshold, and compares the results. It creates a directory with
a text report and plots of the results. At the moment, the following are
included in the report:

All PSSMs specified in the search process
The distribution of each AA in the peptides found by ImmuNovo (per length)
The distribution of each AA in the peptides found by database search
The number of unique peptides above the threshold found by each method
A plot of length distributions of peptides foud by ImmuNovo
A plot of the number of exact matches, matches with <=2AA differences, and 
more than 2AA differences found by each method
"""

import argparse
import subprocess
import shutil
import pandas as pd
import matplotlib.pyplot as plt

import sys, os
# cheap hack until I can figure out how to do this properly in
# python 3.6
sys.path.append(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])

import analysis.findFDR as findFDR
import analysis.scorePeptides as scorePeptides
import analysis.findUniquePeptides as findUniquePeptides
from backend.constants import *

def parseArguments():

  def checkExists(directory):
    if not os.path.exists(directory):
      print("Directory: {} does not exist".format(directory))
      sys.exit()
  
  def createDir(directory):
    if not os.path.exists(directory):
      os.mkdir(directory)

  abspath = os.path.abspath(__file__)
  dname = os.path.dirname(abspath)
  
  parser = argparse.ArgumentParser()
  parser.add_argument('immunovo_results_dir',
                      help='A directory containing results of the ImmuNovo Program')
  parser.add_argument('database_results_dir',
                      help='A directory containing results of database search')
  parser.add_argument('spec_dir',
                      help='A directory containing the spectrum files that were searched')
  parser.add_argument('pssm_dir',
                      help='A directory containing pssm files searched')
  parser.add_argument('acid_mass_file',
                      help='A file containing amino acid mass data')
  parser.add_argument('decoy_dir',
                      help='A directory containing decoy peptides for FDR calculation')
  parser.add_argument('plt_title',
                      help='A title to put on all the plots')


  parser.add_argument('-f', '--FDR',
                      dest='fdr',
                      default = 0.01,
                      help='The target FDR', type=float)
  parser.add_argument('-minP', '--Minimum-Peptide-Length',
                      dest='minP',
                      type=int,
                      default=9,
                      help='The minimum length of a peptide in our search')
  parser.add_argument('-maxP', '--Maximum-Peptide-Length',
                      dest='maxP',
                      type=int,
                      default=12,
                      help='The maximum length of a peptide in our search')
  parser.add_argument('-o', '--Output-Directory',
                      dest='output_dir',
                      default=os.path.join(dname, './reports'))
  parser.add_argument('-db', '--database-search-program',
                        dest='database',
                        default=MSGF,
                        help='The database search program used (only MSGF)')
  parser.add_argument('-d', '--max-decoys',
                        dest='decoys',
                        default=10,
                        help='The number of decoy peptides to be considered in calculating FDR')
  parser.add_argument('-rev', '--reversed_decoys',
                        dest='rev',
                        default='t',
                        help='Are the decoy peptides in the database reversed? ("t" or "f")')
  parser.add_argument('-st', '--Score-Type',
                        dest='scoreType',
                        default=SCORE_COMBINED,
                        help='Score type to compare (defaults to combined score)')
  parser.add_argument('-prec', '--Precision',
      dest='prec',
      type=int,
      default=4,
      help='Decimial precision of our search')


  arguments = parser.parse_args()
  arguments.immunovo_results_dir = os.path.abspath(arguments.immunovo_results_dir)
  arguments.database_results_dir = os.path.abspath(arguments.database_results_dir)
  arguments.spec_dir = os.path.abspath(arguments.spec_dir)
  arguments.pssm_dir = os.path.abspath(arguments.pssm_dir)
  arguments.acid_mass_file = os.path.abspath(arguments.acid_mass_file)
  arguments.decoy_dir = os.path.abspath(arguments.decoy_dir)
  arguments.output_dir = os.path.abspath(arguments.output_dir)

  checkExists(arguments.immunovo_results_dir)
  checkExists(arguments.database_results_dir)
  checkExists(arguments.spec_dir)
  checkExists(arguments.pssm_dir)
  checkExists(arguments.acid_mass_file)
  checkExists(arguments.decoy_dir)

  createDir(arguments.output_dir)
  
  if not (0 < arguments.fdr < 1):
    print(str.format('{} is not a valid FDR. Should be between (0, 1)', arguments.fdr))
    sys.exit()
  
  return arguments


def runPepToScores(arguments, dname):
  
  subprocess.run([os.path.join(dname, 'dbPepToScores.py'),
                  str(arguments.prec),
                  arguments.acid_mass_file,
                  arguments.pssm_dir,
                  arguments.decoy_dir,
                  arguments.rev])
  
  # Keep the decoy search results in case we want to re-run anything later
  os.rename('decoyScores.csv', os.path.join(arguments.output_dir, 'decoyScores.csv'))

  return os.path.join(arguments.output_dir, 'decoyScores.csv')


def runFDR(arguments, decoyPeptides):
  combinedResults = []
  for fileName in os.listdir(arguments.immunovo_results_dir):
    combinedResults.append(pd.read_csv(os.path.join(arguments.immunovo_results_dir, fileName)))
  resultsDF = pd.concat(combinedResults)
  decoyDF = pd.read_csv(decoyPeptides)

  mergedScores = findFDR.getScores(resultsDF, decoyDF, arguments.scoreType)
  
  fdrCutoffs = findFDR.findFDR(mergedScores, arguments.fdr)

  finalResults = findFDR.addFDR(resultsDF, fdrCutoffs, arguments.scoreType)
  finalResults[finalResults[PEPTIDE] != NO_PEP].to_csv(arguments.output_file)

  return fdrCutoffs, finalResults


def processDatabaseData(arguments, databaseDir, fdrCutoffs, scoreType, databaseType):
  # databaseType only MSGF for now

  combinedResults = []
  for fileName in os.listdir(databaseDir):
    combinedResults.append(pd.read_csv(os.path.join(databaseDir, fileName), sep='\t'))
  dbDF = pd.concat(combinedResults)

# def scorePeptides(peptideDF,
#                 acidMassFile,
#                 pssmDir,
#                 spectrumDirectory,
#                 precision,
#                 minPepLength,
#                 maxPepLength,
#                 maxMassTolerance,
#                 compression)

  dbDF = dbDF[['Title', 'Peptide']]
  dbDF = scorePeptides.scorePeptides(dbDF,
                                    arguments.acid_mass_file,
                                    arguments.pssm_dir,
                                    arguments.spec_dir,
                                    arguments.prec,
                                    arguments.minP,
                                    arguments.maxP,
                                    arguments.mmt,
                                    arguments.comp)

  return findFDR.addFDR(dbDF, fdrCutoffs, scoreType)


def compareResults(immuNovoDict, databaseDict):

  def hammingDistance(pep1, pep2, maxDist=2):
    distance = 0
    for aa1, aa2 in zip(pep1, pep2):
      if aa1 != aa2:
        distance += 1
      if distance > maxDist:
        return distance

  numIdentical = 0
  num2AA = 0
  total = sum([len(immuNovoDict[length]) for length in immuNovoDict])

  for length in immuNovoDict:
    for iPep in immuNovoDict[length]:
      if iPep in databaseDict[length]:
        numIdentical += 1
      else:
        for dPep in databaseDict[length]:
          if hammingDistance(iPep, dPep) <= 2:
            num2AA += 1
            break
      
  return numIdentical, num2AA


def copyPSSM(pssmDir, ouputDir):
  for f in os.listdir(pssmDir):
    shutil.copy(f, ouputDir)


def aminoAcidDistribution(peptideDict):
  def analyzeLength(peptides, length):
    acidDistribution = {}
    for pep in peptides:
      for i in range(length):
        if pep[i] not in acidDistribution:
          acidDistribution[i] = [0 for x in range(length)]
        acidDistribution[pep[i]][i] += 1

  distribution = {}
  for length in peptideDict:
    distribution[length] = analyzeLength(peptideDict[length])
  return distribution


def acidDistributionToString(distribution):
  distString = []
  for aa in distribution:
    distString.append(aa + ' '.join([str(x) for x  in distribution]))
  return '\n'.join(distString)


def pepCountToString(immuNovoDict, databaseDict, numMatch, num2AA):
  outString = 'Exact Matches: {}\n'.format(numMatch)
  outString += '2AA Difference: {}\n'.format(num2AA)
  outString += 'No Match: {}\n'.format(sum([len(immuNovoDict[x] for x in immuNovoDict)]) - (num2AA + numMatch))
  outString += 'ImmuNovo Lengths\n'
  outString += ' '.join(['{}: {}'.format(x, len(immuNovoDict[x])) for x in immuNovoDict])
  outString += 'Database Lengths\n'
  outString += ' '.join(['{}: {}'.format(x, len(databaseDict[x])) for x in databaseDict])


def plotMatches(immuNovoDict, numMatch, num2AA, outputDir, plotTitle):
  totalPeptides = sum([len(immuNovoDict[x] for x in immuNovoDict)])

  matchString = 'Identical'
  less2AAString = 'Different with\nno more than 2AA'
  more2AAString = 'Different with\nmore than 2AA'

  plt.bar([1,2,3],
          [numMatch, num2AA, totalPeptides - (numMatch + num2AA)],
          width=0.4)
  plt.xticks([1,2,3],[matchString, less2AAString, more2AAString])
  plt.ylabel('count')
  plt.title(plotTitle)

  plt.savefig(os.path.join(outputDir, 'matches.png', dpi=600))


def plotLengths(immuNovoDict, outputDir, plotTitle):
  lengthStrings = []
  lengths = []

  for length in immuNovoDict:
    lengths.append(len(immuNovoDict[length]))
    lengthStrings.append('length {}'.format(length))
  
  plt.bar([i+1 for i in range(len(lengths))], lengths, width=0.4)
  plt.xticks([i+1 for i in range(len(lengths))], [lengthStrings])
  plt.ylabel('count')
  plt.title(plotTitle)

  plt.savefig(os.path.join(outputDir, 'lengths.png', dpi=600))


if __name__ == '__main__':
  arguments = parseArguments()

  abspath = os.path.abspath(__file__)
  dname = os.path.dirname(abspath)

  decoyPeptides = runPepToScores(arguments)

  fdrCutoffs, fdrImmuNovo = runFDR(arguments, decoyPeptides)
  fdrDatabase = \
      processDatabaseData(arguments.database_results_dir,
                          fdrCutoffs,
                          SCORE_COMBINED,
                          MSGF)

  immuNovoDict = findUniquePeptides.getPeptideDict(fdrImmuNovo, arguments.fdr)
  databaseDict = findUniquePeptides.getPeptideDict(fdrDatabase, arguments.fdr)

  numIdentical, num2AA = compareResults(immuNovoDict, databaseDict)

  # Create report
  copyPSSM(arguments.pssm_dir, arguments.output_dir)
  plotMatches(immuNovoDict,
              numIdentical,
              num2AA,
              arguments.output_dir,
              arguments.plt_title)
  plotLengths(immuNovoDict, arguments.outputDir, arguments.plt_title)

  with open(os.path.join(arguments.output_dir, 'report.txt'), 'w') as f:
    f.write(pepCountToString(immuNovoDict, databaseDict, numIdentical, num2AA))
    f.write('\n\n')
    f.write('ImmuNovo Lengths\n')
    for length in immuNovoDict:
      f.write('Length {}\n'.format(length))
      f.write(acidDistributionToString(immuNovoDict[length]))

    f.write('Database Lengths\n')
    for length in databaseDict:
      f.write('Length {}\n'.format(length))
      f.write(acidDistributionToString(databaseDict[length]))
  