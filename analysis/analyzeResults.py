"""
***** Part of this code was written specifically for the Darwin server     *****
***** Part of this code relies on a fixed amino acid list                  *****
***** Modifications have been removed from reported peptides               *****

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

# read Denovo Results, remove duplicates
# read database Results, remove duplicates
# read decoy peptides, remove duplicates



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

import analysis.findFDR as FindFDR
import analysis.scorePeptides as ScorePeptides
import analysis.selectDecoys as SelectDecoys
import analysis.findUniquePeptides as FindUniquePeptides
from main import getAminoVariables
from backend.constants import *

QVALUE = 'QValue'
SOURCE = 'SOURCE'
SOURCE_DENOVO = 'SOURCE_DENOVO'
SOURCE_DECOY = 'SOURCE_DECOY'
SOURCE_DATABASE = 'SOURCE_DATABASE'

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
  parser.add_argument('spec_dir',
                      help='A directory containing the spectrum files that were searched')
  parser.add_argument('pssm_dir',
                      help='A directory containing pssm files searched')

  parser.add_argument('--update', action='store_true')


  #####
  #
  # Add optional arguments here
  #
  #####


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




def copyPSSM(pssmDir, ouputDir):
  for f in os.listdir(pssmDir):
    shutil.copy(os.path.join(pssmDir, f), ouputDir)

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

def getPssmDistribution(pssmByLengthDict):
  # take a dictionary: dict[pssmTitle][length] = peptides
  pssmDistribution = {}

  totalFinds = 0
  for pssmTitle in pssmByLengthDict:
    pssmDistribution[pssmTitle] = 0
    for length in pssmByLengthDict[pssmTitle]:
      totalFinds += len(pssmByLengthDict[pssmTitle][length])
      pssmDistribution[pssmTitle] += len(pssmByLengthDict[pssmTitle][length])
  for pssm in pssmDistribution:
    pssmDistribution[pssm] = \
        round((pssmDistribution[pssm] / totalFinds) * 10**2, 2)
  return pssmDistribution

def getPssmLengthDistribution(pssmByLengthDict):
  # Get a distribution of peptides by pssm and length
  def getLengthDistribution(lengthDict):
    totalFinds = 0
    lengthDistribution = {}
    for length in lengthDict:
      totalFinds += lengthDict[length]
    for length in lengthDict:
      lengthDistribution[length] = \
          round((lengthDict[length] / totalFinds) * 10**2, 2)
    return lengthDistribution
  pssmLengthDistribution = {}
  for pssmTitle in pssmByLengthDict:
    pssmLengthDistribution[pssmTitle] = \
          getLengthDistribution(pssmByLengthDict[pssmTitle])
  return pssmLengthDistribution

def getLengthDistribution(lengthDict, uniquePeptides):
  lengthDistribution = {}
  for length in lengthDict:
    lengthDistribution[length] = lengthDict[length] / len(uniquePeptides)
  return lengthDistribution

def getLengthCountDict(peptideDF):
  lengthCount = {}
  for item in peptideDF.groupby(LENGTH):
    lengthCount[item[0]] = len(item[1])
  return lengthCount

def getLengthDistributionString(lengthCounts, lengthPercents):
  return ' '.join(['{}: {} ({}%)'.format(length, lengthCounts[length], lengthPercents[length]) for length in lengthCounts])

def getPssmDistributionString(pssmDistribution):
  return '\n'.join(['{}: {}%'.format(pssmTitle, pssmDistribution[pssmTitle]) for pssmTitle in pssmDistribution])

def getPssmLengthDistributionString(pssmLengthDict, pssmLengthDistribution):
  def distributionString(lengthCounts, lengthPercents):
    return ' '.join(['{}: {} ({}%)'.format(length, lengthCounts[length], lengthPercents[length]) for length in lengthCounts])
  return '\n'.joint(['{} -- {}'.format(pssmTitle, distributionString(pssmLengthDict[pssmTitle], pssmLengthDistribution[pssmTitle])) for pssmTitle in pssmLengthDict])



def getSpectrumHits(peptideDF):
  return set(peptideDF[TITLE_SPECTRUM].drop_duplicates())

def uniquePeptides(peptideDF):
  return set(peptideDF[PEPTIDE].drop_duplicates())

def generateTSLPeptides(length, number):
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

def printGraphicTSL(pssmDistributionDict, dataSetName, tslLocation, outputDirectory):

  def stripModifications(peptide):
    # We're just going to strip the peptide modifications and call it good
    # enough for this
    return [x for x in peptide if x in string.ascii_letters]

  tmpPeptideFile = os.path.join(outputDirectory, 'tmpPep.txt')
  tmpDecoyFile = os.path.join(outputDirectory, 'tmpDecoy.txt')
  # Let's report the number of peptides for each pssm
  for pssmName in pssmDistributionDict:
    for length in pssmDistributionDict[pssmName]:
      peptideList = [stripModifications(x) for x in pssmDistributionDict[pssmName][length]]
      decoys = generateTSLPeptides(len(peptideList[0]), len(peptideList))
      outputName = \
          os.path.join(outputDirectory, '{}_{}_{}'.format(dataSetName,
                                                          pssmName.replace('.csv.pssm', ''),
                                                          length))
      with open(tmpPeptideFile, 'w') as f:
        f.write('\n'.join(peptideList))
      with open(tmpDecoyFile, 'w') as f:
        f.write('\n'.join(decoys))
      subprocess.run("ruby {} -P {} -N {} -K A -O {} -R 600 -x -y -I {}".format(os.path.abspath(tslLocation),
                                                                                os.path.abspath(tmpPeptideFile),
                                                                                os.path.abspath(tmpDecoyFile),
                                                                                os.path.abspath(outputName + '.png'),
                                                                                pssmName.replace('.csv.pssm', '')).split())
  os.remove(tmpPeptideFile)
  os.remove(tmpDecoyFile)

def getPSSMDistributionDict(peptideDF):
  # Return a dictionary: dict[pssmTitel][peptideLength] = set(peptides)
  pssmDistributionDict = {}
  for pssmItem in peptideDF.groupby(TITLE_PSSM):
    pssmDistributionDict[pssmItem[0]] = {}
    for lengthItem in pssmItem[1].groupby(LENGTH):
      pssmDistributionDict[pssmItem[0]][lengthItem[0]] = \
                                      set(lengthItem[1][PEPTIDE])
  return pssmDistributionDict

def convertPeptide(peptideString, allConversions, acidConversion):
  # peptideString, list AAs we actualy want to convert, conversion dictionary
  for conversion in allConversions:
    peptideString = \
        peptideString.replace(conversion, acidConversion[conversion])
    return peptideString

def convertAllPeptides(peptideList, acidConversion):

  allConversions = []
  for key in acidConversion:
    if key != acidConversion[key]:
      allConversions.append(key)
  allConversions.sort(key=lambda x: len(x), reverse=True)
  
  convertedPeptides = []
  for peptide in peptideList:
    convertedPeptides.append(convertPeptide(peptide, allConversions, acidConversion))

  return convertedPeptides

def createConversionDict(peptideList, conversionTable):
  convertedPeptides = convertAllPeptides(peptideList, conversionTable)
  conversionDict = {}
  for pep, convertedPep in zip(peptideList, convertedPeptides):
    conversionDict[pep] = convertedPep
  return conversionDict

def createPeptideByLengthDict(peptideDF):
  # Return Dictionary with length keys and peptide set entries
  peptideDict = {}
  for item in peptideDF.groupby(LENGTH):
    peptideDict[item[0]] = set(item[1][PEPTIDE])
  return peptideDict

def getOverlap(denovoPeptides, databasePeptides):
  return denovoPeptides.intersection(databasePeptides)

def closestFDR(resultsDF, fdrCutoff, increment=0.01):
  while len(resultsDF[resultsDF[FDR] >= fdrCutoff]) == 0:
    fdrCutoff += increment
  return resultsDF[resultsDF[FDR] >= fdrCutoff], fdrCutoff

def separateDataFrames(df, qValue):
  dfDict = {}
  for item in df.groupby(SOURCE):
    dfDict[item[0]] = item[1]
  if not qValue:
    return dfDict[SOURCE_DENOVO].drop(SOURCE, axis=1), \
           dfDict[SOURCE_DECOY].drop(SOURCE, axis=1)
  return dfDict[SOURCE_DENOVO].drop(SOURCE, axis=1), \
         dfDict[SOURCE_DECOY].drop(SOURCE, axis=1), \
         dfDict[SOURCE_DATABASE].drop(SOURCE, axis=1)

def mergeDataFrames(denovoDF, decoyDF, databaseDF, qValue):
  denovoDF[SOURCE] = [SOURCE_DENOVO for i in range(len(denovoDF))]
  decoyDF[SOURCE] = [SOURCE_DECOY for i in range(len(decoyDF))]
  combinedDF = pd.concat(denovoDF, decoyDF)
  if not qValue:
    databaseDF[SOURCE] = [SOURCE_DATABASE for i in range(databaseDF)]
    combinedDF = pd.concat(combinedDF, databaseDF)
  return combinedDF

def importDenovoData(denovoDir):
  dataframeList = []
  for fileName in os.listdir(denovoDir):
    fileLocation = os.path.join(denovoDir, fileName)
    dataframeList.append(pd.read_csv(fileLocation))
  return pd.concat(dataframeList)

def importDatabaseData(databaseDir, databaseType, minPepLength, maxPepLength):
  # databaseType only MSGF for now

  combinedResults = []
  for fileName in os.listdir(databaseDir):
    combinedResults.append(pd.read_csv(os.path.join(databaseDir, fileName), sep='\t'))
  dbDF = pd.concat(combinedResults)

  if QVALUE in dbDF:
    dbDF = dbDF[['Title', 'Peptide', QVALUE]]
  else:
    dbDF = dbDF[['Title', 'Peptide']] # Add QValue
  dbDF = dbDF.rename(columns = {'Title': TITLE_SPECTRUM, 'Peptide': PEPTIDE, QVALUE: FDR})

  dbDF['length'] = \
    dbDF.apply(lambda x: len([y for y in x[PEPTIDE] if y in string.ascii_letters]), axis=1)
  dbDF = dbDF[dbDF['length'] >= minPepLength]
  dbDF = dbDF[dbDF['length'] <= maxPepLength]

  return dbDF.drop('length', axis=1)

def resultsFilter(peptideDF):
  # Assume column headers are standard for peptides and spectrum titles
  # Assume lengths are filtered
  # Assume not using QValue
  peptideDF = peptideDF[[PEPTIDE, TITLE_SPECTRUM]]
  peptideDF = peptideDF[peptideDF[PEPTIDE] != NO_PEP]
  peptideDF[PEPTIDE] = peptideDF.apply(lambda x: x[PEPTIDE].replace('I', 'L'))
  peptideDF = peptideDF.drop_duplicates()

  return peptideDF

def filterTopPeptides(peptideDF, scoreType):
  # Take dataframe filtered by FDR, with scores, and return a dataframe
  # with only the top scoring peptide for each spectrum

  return peptideDF[peptideDF.groupby([TITLE_SPECTRUM])[scoreType].transform(max) == peptideDF[scoreType]]

def addPeptideLength(peptideDF, conversionDict):
  peptideDF[LENGTH] = \
      peptideDF.apply(lambda x: len(conversionDict(x[PEPTIDE])), axis=1)
  return peptideDF

def dataframeSetup(denovoResultsDirectory,
                   databaseResultsDirectory,
                   decoyPeptideDirectory,
                   spectrumFileDirectory,
                   acidMassTable,
                   acidConversionTable,
                   massTolerance=35,
                   minPeptideLength=9,
                   maxPeptideLength=12,
                   maxDecoys=10,
                   precision=4,
                   reverse=False,
                   databaseType=MSGF):
  
  denovoDF = resultsFilter(importDenovoData(denovoResultsDirectory))

  databaseDF = resultsFilter(importDatabaseData(databaseResultsDirectory,
                                                databaseType,
                                                minPeptideLength,
                                                maxPeptideLength))

  decoyDF = resultsFilter(SelectDecoys(decoyPeptideDirectory,
                                       spectrumFileDirectory,
                                       acidMassTable,
                                       acidConversionTable,
                                       massTolerance,
                                       minPeptideLength,
                                       maxPeptideLength,
                                       maxDecoys,
                                       precision,
                                       reverse))

  return denovoDF, decoyDF, databaseDF

def spectrumVariableSetup(acidMassFile,
                          pssmDir,
                          minPepLength,
                          maxPepLength,
                          precision):

  allPSSM, acidMassTable, conversionTable = \
                              getAminoVariables(acidMassFile,
                                                precision,
                                                pssmDir,
                                                minPepLength,
                                                maxPepLength)
  H2OMassAdjusted = int(H2OMASS * (10**precision))
  NH3MassAdjusted = int(NH3MASS * (10**precision))
  protonMassAdjusted = int(PROTONMASS * (10**precision))

  return allPSSM, acidMassTable, conversionTable, \
         H2OMASS, NH3MASS, protonMassAdjusted



def getAnalysis(denovoResultsDirectory,
                databaseResultsDirectory,
                decoyPeptideDirectory,
                spectrumFileDirectory,
                qValue,
                fdrCutoff,
                scoringFunctionChoice,
                scoreComparisionType,
                pssmDir,
                acidMassFile,
                outputDirectory,
                dataSetName,
                tslLocation,
                increment=0.01,
                massTolerance=35,
                minPeptideLength=9,
                maxPeptideLength=12,
                maxDecoys=10,
                precision=4,
                compression=2,
                reverse=False,
                databaseType=MSGF):

  allPSSM, acidMassTable, acidConversionTable, H2OMASSAdjusted, \
    NH3MASSAdjusted, protonMassAdjusted = spectrumVariableSetup(acidMassFile,
                                                                  pssmDir,
                                                                  minPeptideLength,
                                                                  maxPeptideLength,
                                                                  precision)
  
  denovoDF, decoyDF, databaseDF = dataframeSetup(denovoResultsDirectory,
                                                 databaseResultsDirectory,
                                                 decoyPeptideDirectory,
                                                 spectrumFileDirectory,
                                                 acidMassTable,
                                                 acidConversionTable,
                                                 massTolerance,
                                                 minPeptideLength,
                                                 maxPeptideLength,
                                                 maxDecoys,
                                                 precision,
                                                 reverse,
                                                 databaseType)
  
  mergedDF = mergeDataFrames(denovoDF, decoyDF, databaseDF, qValue)

  peptideConversionDict = \
      createConversionDict(list(mergedDF[PEPTIDE]), acidConversionTable)
  
  denovoDF = addPeptideLength(denovoDF, peptideConversionDict)
  databaseDF = addPeptideLength(denovoDF, peptideConversionDict)
  databaseDF = addPeptideLength(denovoDF, peptideConversionDict)
  mergedDF = addPeptideLength(denovoDF, peptideConversionDict)

  mergedDF = ScorePeptides.getPeptideScores(mergedDF,
                                            peptideConversionDict,
                                            acidMassTable,
                                            acidConversionTable,
                                            H2OMASSAdjusted,
                                            NH3MASSAdjusted,
                                            protonMassAdjusted,
                                            spectrumFileDirectory,
                                            scoringFunctionChoice,
                                            precision,
                                            minPeptideLength,
                                            maxPeptideLength,
                                            massTolerance,
                                            compression)
  
  if qValue:
    denovoDF, decoyDF, databaseDF = separateDataFrames(mergedDF, qValue)
  else:
    denovoDF, decoyDF = separateDataFrames(mergedDF, qValue)
  
  # No FDR yet, but we're going to add it to this variable
  fdrDenovoDF = filterTopPeptides(denovoDF, scoreComparisionType)
  fdrDecoyDF = filterTopPeptides(decoyDF, scoreComparisionType)

  if not qValue:
    fdrDatabaseDF = filterTopPeptides(databaseDF, scoreComparisionType)
  else:
    fdrDatabaseDF = databaseDF
  
  denovoDF = FindFDR.addFDRToDataframe(denovoDF,
                                       fdrDecoyDF,
                                       scoreComparisionType,
                                       precision,
                                       increment)

  if not qValue:
    fdrDatabaseDF = FindFDR.addFDRToDataframe(fdrDatabaseDF,
                                              fdrDecoyDF,
                                              scoreComparisionType,
                                              precision,
                                              increment)
  

  fdrDenovoDF, fdrCutoff = closestFDR(fdrDenovoDF, fdrCutoff, increment)
  fdrDatabaseDF = fdrDatabaseDF[fdrDatabaseDF[FDR] >= fdrCutoff]

  ############# Do Analysis ####################

  # Using our original Dataframes
  outputFileName = os.path.join(outputDirectory, 'report.txt')

  with open(outputFileName, 'w') as f:
    uniquePeptidesDenovo = uniquePeptides(fdrDenovoDF)
    uniquePeptidesDatabase = uniquePeptides(fdrDatabaseDF)
    f.write('FDR used: {}\n'.format(fdrCutoff))
    f.write('\n')
    f.write('Denovo Spectrum Matches: {}\n'.format(getSpectrumHits(denovoDF)))
    f.write('Database Spectrum Matches: {}\n'.format(getSpectrumHits(databaseDF)))
    f.write('\n')
    f.write('Denovo Unique Peptides Found: {}\n'.format(uniquePeptidesDenovo))
    f.write('Database Unique Peptides Found: {}\n'.format(uniquePeptidesDatabase))
    f.write('Overlap: {}\n'.format(getOverlap(fdrDatabaseDF, fdrDatabaseDF)))
    f.write('\n')

    f.write('Length Distribution:\n')
    lengthCountDenovo = getLengthCountDict(fdrDenovoDF)
    lengthCountDatabase = getLengthCountDict(fdrDatabaseDF)
    lengthDistributionDenovo = \
      getLengthDistribution(lengthCountDenovo, uniquePeptidesDenovo)
    lengthDistributionDatabase = \
      getLengthDistribution(lengthCountDatabase, uniquePeptidesDatabase)
    f.write('Denovo Length Distribution:\n')
    f.write(getLengthDistributionString(lengthCountDenovo, lengthDistributionDenovo))
    f.write('\n')
    f.write('Database Length Distribution:\n')
    f.write(getLengthDistributionString(lengthCountDatabase, lengthDistributionDatabase))
    f.write('\n\n')


    pssmDistributionDict = getPSSMDistributionDict(fdrDenovoDF)
    f.write('PSSM Distribution:\n')
    f.write(getPssmDistributionString(getPssmDistribution(pssmDistributionDict)))
    f.write('\n\n')

    f.write('PSSM By Length Distribution:\n')
    pssmLengthDistributionDict = getPssmLengthDistribution(pssmDistributionDict)
    f.write(getPssmLengthDistributionString(pssmDistributionDict, pssmLengthDistributionDict))

  printGraphicTSL(getPSSMDistributionDict(fdrDenovoDF),
                  dataSetName,
                  tslLocation,
                  outputDirectory)






# if __name__ == '__main__':
#   arguments = parseArguments()

#   abspath = os.path.abspath(__file__)
#   dname = os.path.dirname(abspath)

#   if not arguments.update:
#     decoyPeptides = runPepToScores(arguments, dname)
#   decoyPeptides = os.path.join(arguments.output_dir, 'decoyScores.csv')

#   fdrCutoffs, fdrImmuNovo = runFDR(arguments, decoyPeptides)
#   fdrImmuNovo.to_csv(os.path.join(arguments.output_dir, 'processedImmunovo.csv'), index=0)

#   immuNovoDict, fdrCutoff = \
#     findUniquePeptides.getPeptideDict(fdrImmuNovo, arguments.fdr, False)

#   groupedDF = groupByPSSM(fdrImmuNovo, fdrCutoff)

#   # Create tsl images. Written for Darwin Server, specifically
#   os.system("module load ruby/2.1.0")
#   oldCWD = os.getcwd()
#   os.chdir(os.path.dirname(arguments.tsl))
#   pssmPeptides = printGraphicTSL(groupedDF, arguments.dataset_name,
#                                  arguments.tsl, arguments.output_dir)
#   os.chdir(oldCWD)

#   ##### Done With ImmuNovo Only analysis #####
#   if not arguments.immunovo_only:

#     if arguments.update:
#       fdrDatabase = pd.read_csv(os.path.join(arguments.output_dir, 'processedDatabase.csv'))
#     else: 
#       fdrDatabase = \
#           processDatabaseData(arguments.database_results_dir,
#                               fdrCutoffs,
#                               SCORE_COMBINED,
#                               MSGF)
      
#       # This is taking a while, so in case something crashes
#       fdrDatabase.to_csv(os.path.join(arguments.output_dir, 'processedDatabase.csv'), index=0)

#     databaseDict, fdrCutoff = \
#       findUniquePeptides.getPeptideDict(fdrDatabase, fdrCutoff)

#     numIdentical, num2AA, similarity, overlapPeptides = \
#                           compareResults(immuNovoDict, databaseDict)

#   #### Start creating report ####

#   for pssm in pssmPeptides:
#     for length in pssmPeptides[pssm]:
#       with open(os.path.join(arguments.output_dir, '{}_{}_peptides.txt'.format(pssm.replace('.csv.pssm', ''), length)), 'w') as f:
#         f.write('\n'.join(pssmPeptides[pssm][length]))

#   if not arguments.immunovo_only:
#     copyPSSM(arguments.pssm_dir, arguments.output_dir)
#     plotMatches(immuNovoDict,
#                 numIdentical,
#                 num2AA,
#                 arguments.output_dir,
#                 arguments.plt_title)
#     plotLengths(immuNovoDict, arguments.output_dir, arguments.plt_title)

#   # Output summary statistics
#   with open(os.path.join(arguments.output_dir, 'report.txt'), 'w') as f:
#     f.write('FDR cutoff used: {}\n\n'.format(fdrCutoff))

#     if not arguments.immunovo_only:
#       iSpec, iPSSM, dSpec, dPSSM, oSpec, oPSSM = \
#             getScores(fdrImmuNovo, fdrDatabase, overlapPeptides, fdrCutoff)
#       f.write('ImmuNovo Spectrum Score: {}\nImmuNovo PSSM Score: {}\n'.format(iSpec, iPSSM))
#       f.write('Database Spectrum Score: {}\nDatabase PSSM Score: {}\n'.format(dSpec, dPSSM))
#       f.write('Overlapping Spectrum Score: {}\nOverlapping PSSM Score: {}\n'.format(oSpec, oPSSM))
#       f.write(pepCountToString(immuNovoDict, databaseDict, numIdentical, num2AA, similarity))
#       f.write('\n\n')
#     else:
#       iSpec, iPSSM, dSpec, dPSSM, oSpec, oPSSM = \
#             getScores(fdrImmuNovo, fdrImmuNovo, set(), fdrCutoff)
#       immuNovoLengths = [len(immuNovoDict[x]) for x in immuNovoDict]
#       totalImmuNovo = sum(immuNovoLengths)
#       immuNovoPercentages = [round((x / totalImmuNovo) * 10**2, 2) for x in immuNovoLengths]
#       f.write(' '.join(['{}: {} ({}%)'.format(x, len(immuNovoDict[x]), y) for x, y in zip(immuNovoDict, immuNovoPercentages)]))
#       f.write('\n')
#       f.write('ImmuNovo Spectrum Score: {}\nImmuNovo PSSM Score: {}\n'.format(iSpec, iPSSM))
#       f.write('\n\n')

#     denovoPSSM = getPSSMDistribution(groupedDF)
#     f.write('PSSM Distribution\n')
#     f.write('\n'.join(['{}: {}%'.format(pssm, round(denovoPSSM[pssm], 2)) for pssm in denovoPSSM]))
#     f.write('\n')

#     for pssm in pssmPeptides:
#       f.write(pssm.replace('.csv.pssm', '') + ' -- ')
#       for length in pssmPeptides[pssm]:
#         f.write('{}: {} '.format(length, len(pssmPeptides[pssm][length])))
#       f.write('\n')
#     f.write('\n')
        

#     f.write('ImmuNovo Lengths\n')
#     immuNovoDistribution = aminoAcidDistribution(immuNovoDict)
#     for length in immuNovoDistribution:
#       f.write('Length {}\n'.format(length))
#       f.write(acidDistributionToString(immuNovoDistribution[length]))
#       f.write('\n')
#     f.write('\n\n')
    

#     if not arguments.immunovo_only:
#       f.write('Database Lengths\n')
#       databaseDistribution = aminoAcidDistribution(databaseDict)
#       for length in databaseDistribution:
#         f.write('Length {}\n'.format(length))
#         f.write(acidDistributionToString(databaseDistribution[length]))
#         f.write('\n')
  