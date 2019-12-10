"""
***** Part of this code was written specifically for the Darwin server     *****
***** Part of this code relies on a fixed amino acid list                  *****

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
import math
import statistics
import subprocess
import string
import random

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

  abspath = os.path.abspath(__file__)
  dname = os.path.dirname(abspath)

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
  parser.add_argument('dataset_name',
                      help='The name of the current dataset being processed')
  parser.add_argument('plt_title',
                      help='A title to put on all the length distribution plots')


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
  parser.add_argument('-mmt', '--Maximum-Mass-Tolerance',
                      dest='mmt',
                      type=int,
                      default=35,
                      help='The maximum deviance between a peptide mass and spectrum mass '
                      'before we consider the masses to not be matched (PPM)')
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
  parser.add_argument('-comp', '--Spectrum-Compression',
                      dest='comp',
                      type=int,
                      default=2,
                      help='The spectrum intensity log compression level')

  parser.add_argument('-tsl', '--TSL_Location',
                      dest='tsl',
                      type=str,
                      default=os.path.join(dname, './tsl/cgi-bin/tsl'),
                      help='Path to the tsl binary')
  parser.add_argument('--update', action='store_true')


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
  
  os.system(' '.join(['python3',
                      os.path.join(dname, 'dbPepToScores.py'),
                      arguments.spec_dir,
                      arguments.acid_mass_file,
                      arguments.pssm_dir,
                      arguments.decoy_dir,
                      arguments.rev,
                      arguments.decoys]))

  # subprocess.run([os.path.join(dname, 'dbPepToScores.py'),
  #                 str(arguments.prec),
  #                 arguments.acid_mass_file,
  #                 arguments.pssm_dir,
  #                 arguments.decoy_dir,
  #                 arguments.rev,
  #                 arguments.decoys])
  
  # Keep the decoy search results in case we want to re-run anything later
  os.rename(os.path.join(dname, 'decoyScores.csv'), os.path.join(arguments.output_dir, 'decoyScores.csv'))

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

  return fdrCutoffs, finalResults


def processDatabaseData(databaseDir, fdrCutoffs, scoreType, databaseType):
  # databaseType only MSGF for now

  combinedResults = []
  for fileName in os.listdir(databaseDir):
    combinedResults.append(pd.read_csv(os.path.join(databaseDir, fileName), sep='\t'))
  dbDF = pd.concat(combinedResults)

  dbDF = dbDF[['Title', 'Peptide']]
  dbDF = dbDF.rename(columns = {'Title': TITLE_SPECTRUM, 'Peptide': PEPTIDE})
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
      # if distance > maxDist:
      #   return distance
    return distance
  
  overlapPeptides = set()

  numIdentical = 0
  num2AA = 0

  averageSimilarity = 0
  numSimilarity = 0

  for length in immuNovoDict:
    for iPep in immuNovoDict[length]:
      if iPep in databaseDict[length]:
        overlapPeptides.add(iPep)
        numIdentical += 1
      else:
        bestHamming = math.inf
        for dPep in databaseDict[length]:
          currentHamming = hammingDistance(iPep, dPep)
          if currentHamming <= 2:
            overlapPeptides.add(iPep)
            num2AA += 1
            break
          elif bestHamming > currentHamming:
            bestHamming = currentHamming
        if currentHamming > 2:
          numSimilarity += 1
          averageSimilarity += (bestHamming / length)

  return numIdentical, num2AA, (averageSimilarity / numSimilarity), overlapPeptides


def copyPSSM(pssmDir, ouputDir):
  for f in os.listdir(pssmDir):
    shutil.copy(os.path.join(pssmDir, f), ouputDir)


def aminoAcidDistribution(peptideDict):
  def analyzeLength(peptides, length):
    acidDistribution = {}
    for pep in peptides:
      for i in range(length):
        if pep[i] not in acidDistribution:
          acidDistribution[pep[i]] = [0 for x in range(length)]
        acidDistribution[pep[i]][i] += 1
    return acidDistribution

  distribution = {}
  for length in peptideDict:
    distribution[length] = analyzeLength(peptideDict[length], length)
  return distribution


def acidDistributionToString(distribution):
  distString = []
  for aa in distribution:
    distString.append(aa + ': ' + ' '.join([str(x) for x  in distribution[aa]]))
  return '\n'.join(distString)


def pepCountToString(immuNovoDict, databaseDict, numMatch, num2AA, similarity):
  immuNovoLengths = [len(immuNovoDict[x]) for x in immuNovoDict]
  totalImmuNovo = sum(immuNovoLengths)
  immuNovoPercentages = [round((x / totalImmuNovo) * 10**2, 2) for x in immuNovoLengths]
  databaseLengths = [len(databaseDict[x]) for x in databaseDict]
  totalDatabase = sum(databaseLengths)
  databasePercentages = [round((x / totalDatabase) * 10**2, 2) for x in databaseLengths]

  outString = 'Exact Matches: {}\n'.format(numMatch)
  outString += '2AA Difference: {}\n'.format(num2AA)
  outString += 'No Match: {}\n'.format(sum([len(immuNovoDict[x]) for x in immuNovoDict]) - (num2AA + numMatch))
  outString += 'Remaining Similarity: {}\n'.format(similarity)
  outString += 'ImmuNovo Lengths\n'
  outString += ' '.join(['{}: {} ({}%)'.format(x, len(immuNovoDict[x]), y) for x, y in zip(immuNovoDict, immuNovoPercentages)])
  outString += '\n'
  outString += "Total Peptides: {}\n".format(totalImmuNovo)
  outString += '\n\n'
  outString += 'Database Lengths\n'
  outString += ' '.join(['{}: {} ({}%)'.format(x, len(databaseDict[x]), y) for x, y in zip(databaseDict, databasePercentages)])
  outString += '\n'
  outString += "Total Peptides: {}\n".format(totalDatabase)
  outString += '\n\n'

  return outString


def plotMatches(immuNovoDict, numMatch, num2AA, outputDir, plotTitle):
  totalPeptides = sum([len(immuNovoDict[x]) for x in immuNovoDict])

  matchString = 'Identical'
  less2AAString = 'Different with\nno more than 2AA'
  more2AAString = 'Different with\nmore than 2AA'

  plt.bar([1,2,3],
          [numMatch, num2AA, totalPeptides - (numMatch + num2AA)],
          width=0.4)
  plt.xticks([1,2,3],[matchString, less2AAString, more2AAString])
  plt.ylabel('count')
  plt.title(plotTitle)

  plt.savefig(os.path.join(outputDir, 'matches.png'), dpi=600)
  plt.clf()


def plotLengths(immuNovoDict, outputDir, plotTitle):
  lengthStrings = []
  lengths = []

  for length in immuNovoDict:
    lengths.append(len(immuNovoDict[length]))
    lengthStrings.append('length {}'.format(length))
  
  plt.bar([i+1 for i in range(len(lengths))], lengths, width=0.4)
  plt.xticks([i+1 for i in range(len(lengths))], lengthStrings)
  plt.ylabel('count')
  plt.title(plotTitle)

  plt.savefig(os.path.join(outputDir, 'lengths.png'), dpi=600)
  plt.clf()


def getScores(immuNovoDF, databaseDF, overlappingPeptides, fdrCutoff):
  def getStats(df, fdrCutoff):
    # The following filters should already be applied, just making sure
    df = df[df[PEPTIDE] != NO_PEP]
    df = df[df[FDR] <= fdrCutoff]
    meanSpec = statistics.mean(list(df[SCORE_GLOBAL]))
    meanPSSM = statistics.mean(list(df[SCORE_PSSM]))

    return meanSpec, meanPSSM

  separateImmuNovo = \
      immuNovoDF[~immuNovoDF[PEPTIDE].isin(overlappingPeptides)]
  separateDatabase = \
      databaseDF[~databaseDF[PEPTIDE].isin(overlappingPeptides)]
  overlapping = \
      pd.concat([immuNovoDF[immuNovoDF[PEPTIDE].isin(overlappingPeptides)],
                 databaseDF[databaseDF[PEPTIDE].isin(overlappingPeptides)]])
  
  iSpec, iPSSM = getStats(separateImmuNovo, fdrCutoff)
  dSpec, dPSSM = getStats(separateDatabase, fdrCutoff)
  oSpec, oPSSM = getStats(overlapping, fdrCutoff)

  return iSpec, iPSSM, dSpec, dPSSM, oSpec, oPSSM


# def getTopPeptides(spectralGroups):
#   topPeptides = set()
#   for elm in spectralGroups:
#     df = elm[1]
#     df = df[df[SCORE_COMBINED] == df[SCORE_COMBINED].max()]
#     # If equal, pick greastest PSSM Score
#     df = df[df[SCORE_PSSM] == df[SCORE_PSSM].max()]
#     # If all equal, just pick the first one
#     topPeptides.add(list(df[PEPTIDE])[0])
#   return topPeptides


def groupByPSSM(df, fdrCutoff):

  df = df[df[FDR] <= fdrCutoff]

  topPeptides = []
  for elm in df.groupby(TITLE_SPECTRUM):
    curDF = elm[1]
    curDF = curDF[curDF[SCORE_COMBINED] == curDF[SCORE_COMBINED].max()]
    # If equal, pick the greatest PSSM Score
    curDF = curDF[curDF[SCORE_PSSM] == curDF[SCORE_PSSM].max()]
    # If all equal, just pick the first one
    topPeptides.append(curDF)
  df = pd.concat(topPeptides)

  return df.groupby(TITLE_PSSM)


def getPSSMDistribution(groupedDF):
  pssmDistribution = {}

  totalFinds = 0
  for elm in groupedDF:
    totalFinds += len(elm[1])
    pssmDistribution[elm[0]] = len(elm[1])
  for pssm in pssmDistribution:
    pssmDistribution[pssm] /= totalFinds
  return pssmDistribution


def getPSSMPeptides(groupedDF):
  def removeModifications(peptide):
    return ''.join([x for x in peptide if x in string.ascii_letters])
  
  def groupByLength(peptideList):
    groupedPeptides = {}
    for peptide in peptideList:
      if len(peptide) not in groupedPeptides:
        groupedPeptides[len(peptide)] = []
      groupedPeptides[len(peptide)].append(peptide)
    return groupedPeptides
  
  def processDF(df):
    return groupByLength([removeModifications(x) for x in list(df[PEPTIDE])])

  # Assume all modifications are non-alpha characters
  peptideByPSSM = {}
  for item in groupedDF:
    peptideByPSSM[item[0]] = processDF(item[1])
  
  return peptideByPSSM


def generateDecoyPeptides(length, number):
  # Fixed Amino Acid List for now
  aaList = ['G', 'A', 'S', 'P',
            'V', 'T', 'L', 
            'N', 'D', 'Q', 'K',
            'E', 'M', 'H', 'F',
            'R', 'C', 'Y', 'W']
  
  decoys = set()
  while len(decoys) < number:
     decoys.add(''.join([random.choice(aaList) for i in range(length)]))
  return decoys


def printGraphicTSL(groupedDF, dataSetName, tslLocation, outputDirectory):
  peptideByPSSM = getPSSMPeptides(groupedDF)
  tmpPeptideFile = os.path.join(outputDirectory, 'tmpPep.txt')
  tmpDecoyFile = os.path.join(outputDirectory, 'tmpDecoy.txt')
  for pssmName in peptideByPSSM:
    for length in peptideByPSSM[pssmName]:
      peptideList = peptideByPSSM[pssmName][length]
      pssmName = pssmName.replace('.csv.pssm', '')
      decoys = generateDecoyPeptides(len(peptideList[0]), len(peptideList))
      outputName = \
          os.path.join(outputDirectory, '{}_{}_{}'.format(dataSetName,
                                                          pssmName,
                                                          length))
      with open(tmpPeptideFile, 'w') as f:
        f.write('\n'.join(peptideList))
      with open(tmpDecoyFile, 'w') as f:
        f.write('\n'.join(decoys))
      subprocess.run("ruby {} -P {} -N {} -K A -O {} -R 600 -x -y -I {}".format(os.path.abspath(tslLocation),
                                                                                os.path.abspath(tmpPeptideFile),
                                                                                os.path.abspath(tmpDecoyFile),
                                                                                os.path.abspath(outputName + '.png'),
                                                                                pssmName).split())
  os.remove(tmpPeptideFile)
  os.remove(tmpDecoyFile)


if __name__ == '__main__':
  arguments = parseArguments()

  abspath = os.path.abspath(__file__)
  dname = os.path.dirname(abspath)

  if not arguments.update:
    decoyPeptides = runPepToScores(arguments, dname)
  decoyPeptides = os.path.join(arguments.output_dir, 'decoyScores.csv')

  fdrCutoffs, fdrImmuNovo = runFDR(arguments, decoyPeptides)
  fdrImmuNovo.to_csv(os.path.join(arguments.output_dir, 'processedImmunovo.csv'), index=0)

  if arguments.update:
    fdrDatabase = pd.read_csv(os.path.join(arguments.output_dir, 'processedDatabase.csv'))
  else: 
    fdrDatabase = \
        processDatabaseData(arguments.database_results_dir,
                            fdrCutoffs,
                            SCORE_COMBINED,
                            MSGF)
    
    # This is taking a while, so in case something crashes
    fdrDatabase.to_csv(os.path.join(arguments.output_dir, 'processedDatabase.csv'), index=0)

  
  immuNovoDict, fdrCutoff = \
    findUniquePeptides.getPeptideDict(fdrImmuNovo, arguments.fdr, False)
  databaseDict, fdrCutoff = \
    findUniquePeptides.getPeptideDict(fdrDatabase, fdrCutoff)

  numIdentical, num2AA, similarity, overlapPeptides = \
                        compareResults(immuNovoDict, databaseDict)
  
  groupedDF = groupByPSSM(fdrImmuNovo, fdrCutoff)
  
  # Create tsl images. Written for Darwin Server, specifically
  os.system("module load ruby/2.1.0")
  oldCWD = os.getcwd()
  os.chdir(os.path.dirname(arguments.tsl))
  printGraphicTSL(groupedDF, arguments.dataset_name, arguments.tsl, arguments.output_dir)
  os.chdir(oldCWD)


  # Create report
  copyPSSM(arguments.pssm_dir, arguments.output_dir)
  plotMatches(immuNovoDict,
              numIdentical,
              num2AA,
              arguments.output_dir,
              arguments.plt_title)
  plotLengths(immuNovoDict, arguments.output_dir, arguments.plt_title)

  # Output summary statistics
  with open(os.path.join(arguments.output_dir, 'report.txt'), 'w') as f:
    f.write('FDR cutoff used: {}\n\n'.format(fdrCutoff))
    iSpec, iPSSM, dSpec, dPSSM, oSpec, oPSSM = \
          getScores(fdrImmuNovo, fdrDatabase, overlapPeptides, fdrCutoff)
    f.write('ImmuNovo Spectrum Score: {}\nImmuNovo PSSM Score: {}\n'.format(iSpec, iPSSM))
    f.write('Database Spectrum Score: {}\nDatabase PSSM Score: {}\n'.format(dSpec, dPSSM))
    f.write('Overlapping Spectrum Score: {}\nOverlapping PSSM Score: {}\n'.format(oSpec, oPSSM))
    f.write(pepCountToString(immuNovoDict, databaseDict, numIdentical, num2AA, similarity))
    f.write('\n\n')

    denovoPSSM = getPSSMDistribution(groupedDF)
    f.write('PSSM Distribution\n')
    f.write('\n'.join(['{}: {}%'.format(pssm, round(denovoPSSM[pssm], 2)) for pssm in denovoPSSM]))
    f.write('\n\n')

    f.write('ImmuNovo Lengths\n')
    immuNovoDistribution = aminoAcidDistribution(immuNovoDict)
    databaseDistribution = aminoAcidDistribution(databaseDict)
    for length in immuNovoDistribution:
      f.write('Length {}\n'.format(length))
      f.write(acidDistributionToString(immuNovoDistribution[length]))
      f.write('\n')
    f.write('\n\n')
    

    f.write('Database Lengths\n')
    for length in databaseDistribution:
      f.write('Length {}\n'.format(length))
      f.write(acidDistributionToString(databaseDistribution[length]))
      f.write('\n')
  