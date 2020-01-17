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
import backend.Structures.pssm as PSSM
import backend.userInput as UserInput
import backend.PreProcessing.misc as Misc
from backend.constants import *
from copy import deepcopy

QVALUE = 'QValue'
SOURCE = 'SOURCE'
SOURCE_DENOVO = 'SOURCE_DENOVO'
SOURCE_DECOY = 'SOURCE_DECOY'
SOURCE_DATABASE = 'SOURCE_DATABASE'

def parseArguments():

# denovoResultsDirectory    -- required (add to analyzeResults)

  abspath = os.path.abspath(__file__)
  dname = os.path.dirname(abspath)

  def checkExists(directory):
    if not os.path.exists(directory):
      print("{} does not exist".format(directory))
      sys.exit()
  
  def createDir(directory):
    if not os.path.exists(directory):
      os.mkdir(directory)

  parser = argparse.ArgumentParser()

  parser.add_argument('immunovo_results',
                      help='A directory or file containing results of the ImmuNovo Program')
  parser = UserInput.parseDefaultArguments(parser)
  parser.add_argument('-df', '--decoy-file',
      dest='decoy_file',
      default='',
      help='File of decoy peptides for a spectrum file')
  parser = UserInput.parseOptionalArguments(parser, os.getcwd())
  
  parser.add_argument('--update',
                       action='store_true',
                       help='Used for diagnostic purposes')
  
  arguments = UserInput.parseArgumentsSetup(os.getcwd(), parser)

  arguments.immunovo_results = \
        os.path.abspath(arguments.immunovo_results)
  checkExists(arguments.immunovo_results)
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
      totalFinds += len(lengthDict[length])
    for length in lengthDict:
      lengthDistribution[length] = \
          round((len(lengthDict[length]) / totalFinds) * 10**2, 2)
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
  return ' '.join(['{}: {} ({}%)'.format(length, lengthCounts[length], round(lengthPercents[length], 2)) for length in lengthCounts])

def getPssmDistributionString(pssmDistribution):
  return '\n'.join(['{}: {}%'.format(pssmTitle, pssmDistribution[pssmTitle]) for pssmTitle in pssmDistribution])

def getPssmLengthDistributionString(pssmLengthDict, pssmLengthDistribution):
  def distributionString(lengthCounts, lengthPercents):
    return ' '.join(['{}: {} ({}%)'.format(length, len(lengthCounts[length]), lengthPercents[length]) for length in lengthCounts])
  return '\n'.join(['{} -- {}'.format(pssmTitle, distributionString(pssmLengthDict[pssmTitle], pssmLengthDistribution[pssmTitle])) for pssmTitle in pssmLengthDict])



def getSpectrumHits(peptideDF):
  return set(peptideDF[TITLE_SPECTRUM])

def uniquePeptides(peptideDF):
  return set(peptideDF[PEPTIDE])

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
    return ''.join([x for x in peptide if x in string.ascii_letters])

  oldCWD = os.getcwd()
  os.chdir(os.path.dirname(tslLocation))

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
  os.chdir(oldCWD)

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
  for peptide in set(peptideList):
    convertedPeptides.append((peptide, convertPeptide(peptide, allConversions, acidConversion)))

  return convertedPeptides

def createConversionDict(peptideList, conversionTable):
  convertedPeptides = convertAllPeptides(peptideList, conversionTable)
  conversionDict = {}
  for pepInfo in convertedPeptides:
    conversionDict[pepInfo[0]] = pepInfo[1]
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
  if qValue:
    return dfDict[SOURCE_DENOVO].drop(SOURCE, axis=1), \
           dfDict[SOURCE_DECOY].drop(SOURCE, axis=1)
  return dfDict[SOURCE_DENOVO].drop(SOURCE, axis=1), \
         dfDict[SOURCE_DECOY].drop(SOURCE, axis=1), \
         dfDict[SOURCE_DATABASE].drop(SOURCE, axis=1)

def mergeDataFrames(denovoDF, decoyDF, databaseDF, qValue):
  denovoDF[SOURCE] = [SOURCE_DENOVO for i in range(len(denovoDF))]
  decoyDF[SOURCE] = [SOURCE_DECOY for i in range(len(decoyDF))]
  combinedDF = pd.concat([denovoDF, decoyDF])
  if not qValue:
    if type(databaseDF) != type(''):
      databaseDF[SOURCE] = [SOURCE_DATABASE for i in range(len(databaseDF))]
      combinedDF = pd.concat([combinedDF, databaseDF[[PEPTIDE, TITLE_SPECTRUM, SOURCE]]])
  return combinedDF

def importDenovoData(denovoDir):
  if os.path.isdir(denovoDir):
    dataframeList = []
    for fileName in os.listdir(denovoDir):
      fileLocation = os.path.join(denovoDir, fileName)
      dataframeList.append(pd.read_csv(fileLocation))
    return pd.concat(dataframeList)
  elif os.path.isfile(denovoDir):
    return pd.read_csv(denovoDir)

def importDatabaseData(databaseDir, databaseType, minPepLength, maxPepLength):
  # databaseType only MSGF for now

  if os.path.isdir(databaseDir):
    combinedResults = []
    for fileName in os.listdir(databaseDir):
      combinedResults.append(pd.read_csv(os.path.join(databaseDir, fileName), sep='\t'))
    dbDF = pd.concat(combinedResults)
  elif os.path.isfile(databaseDir):
    dbDF = pd.read_csv(databaseDir)

  if QVALUE in dbDF:
    dbDF = dbDF[['Title', 'Peptide', QVALUE]]
    dbDF = dbDF.rename(columns = {'Title': TITLE_SPECTRUM, 'Peptide': PEPTIDE, QVALUE: FDR})
  else:
    dbDF = dbDF[['Title', 'Peptide']] # Add QValue
    dbDF = dbDF.rename(columns = {'Title': TITLE_SPECTRUM, 'Peptide': PEPTIDE})

  dbDF['length'] = \
    dbDF.apply(lambda x: len([y for y in x[PEPTIDE] if y in string.ascii_letters]), axis=1)
  dbDF = dbDF[dbDF['length'] >= minPepLength]
  dbDF = dbDF[dbDF['length'] <= maxPepLength]

  return dbDF.drop('length', axis=1)

def resultsFilter(peptideDF):
  # Assume column headers are standard for peptides and spectrum titles
  # Assume lengths are filtered
  if FDR in peptideDF:
    peptideDF = peptideDF[[PEPTIDE, TITLE_SPECTRUM, FDR]]
  else:
    peptideDF = peptideDF[[PEPTIDE, TITLE_SPECTRUM]]
  peptideDF = peptideDF[peptideDF[PEPTIDE] != NO_PEP]
  peptideDF[PEPTIDE] = peptideDF.apply(lambda x: x[PEPTIDE].replace('I', 'L'), axis=1)
  peptideDF = peptideDF.drop_duplicates()

  return peptideDF

def filterTopPeptides(peptideDF, scoreType):
  # Take dataframe filtered by FDR, with scores, and return a dataframe
  # with only the top scoring peptide for each spectrum

  return peptideDF.loc[[entry[1][scoreType].idxmax() for entry in peptideDF.groupby(TITLE_SPECTRUM)]]
  

def addPeptideLength(peptideDF, conversionDict):
  peptideDF[LENGTH] = \
      peptideDF.apply(lambda x: len(conversionDict[x[PEPTIDE]]), axis=1)
  return peptideDF

def dataframeSetup(denovoResultsDirectory,
                   databaseResultsDirectory,
                   decoyPeptideDirectory,
                   spectrumFileDirectory,
                   acidMassTable,
                   acidConversionTable,
                   fdrOnly,
                   update,
                   outputDirectory,
                   massTolerance=35,
                   minPeptideLength=9,
                   maxPeptideLength=12,
                   maxDecoys=10,
                   precision=4,
                   reverse=False,
                   databaseType=MSGF,
                   decoyFile=''):
  
  if not update:
    denovoDF = resultsFilter(importDenovoData(denovoResultsDirectory))
    denovoDF.to_csv(os.path.join(outputDirectory, 'denovo_peptides.csv'))
  else:
    denovoDF = pd.read_csv(os.path.join(outputDirectory, 'denovo_peptides.csv'))

  if len(denovoDF) == 0:
    sys.exit()

  if fdrOnly or databaseResultsDirectory == '':
    databaseDF = ''
  else:
    databaseDF = resultsFilter(importDatabaseData(databaseResultsDirectory,
                                                  databaseType,
                                                  minPeptideLength,
                                                  maxPeptideLength))

  if (not update) and (decoyFile == ''):
    decoyDF = resultsFilter(SelectDecoys.selectDecoyPeptides(decoyPeptideDirectory,
                                                            spectrumFileDirectory,
                                                            acidMassTable,
                                                            acidConversionTable,
                                                            massTolerance,
                                                            minPeptideLength,
                                                            maxPeptideLength,
                                                            maxDecoys,
                                                            precision,
                                                            reverse))
    decoyDF.to_csv(os.path.join(outputDirectory, 'decoy_peptides.csv'))
  elif update:
    decoyDF = pd.read_csv(os.path.join(outputDirectory, 'decoy_peptides.csv'))
  else:
    decoyDF = pd.read_csv(decoyFile)

  return denovoDF, decoyDF, databaseDF

def spectrumVariableSetup(acidMassFile,
                          pssmDir,
                          minPepLength,
                          maxPepLength,
                          precision):

  allPSSM, acidMassTable, conversionTable = \
                              Misc.getAminoVariables(acidMassFile,
                                                     precision,
                                                     pssmDir,
                                                     minPepLength,
                                                     maxPepLength,
                                                     False)

  allPSSM = PSSM.setForward(allPSSM)
  
  H2OMassAdjusted = int(H2OMASS * (10**precision))
  NH3MassAdjusted = int(NH3MASS * (10**precision))
  protonMassAdjusted = int(PROTONMASS * (10**precision))

  return allPSSM, acidMassTable, conversionTable, \
         H2OMassAdjusted, NH3MassAdjusted, protonMassAdjusted


def addFDR(fdrTargetDF,
           fdrDecoyDF,
           scoreComparisionType,
           precision,
           increment,
           fdrCalculationType):

  if fdrCalculationType == FDR_TDA:
    return FindFDR.addFDRToDataframe(fdrTargetDF,
                                     fdrDecoyDF,
                                     scoreComparisionType,
                                     precision,
                                     increment)
  elif fdrCalculationType == FDR_DISTRIBUTION:
    return FindFDR.addFDRToDataframeOld(fdrTargetDF,
                                        fdrDecoyDF,
                                        scoreComparisionType,
                                        precision,
                                        increment)


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
                fdrCalculationType,
                fdrOnly,
                update,
                increment=0.01,
                massTolerance=35,
                minPeptideLength=9,
                maxPeptideLength=12,
                maxDecoys=10,
                precision=4,
                compression=2,
                reverse=False,
                databaseType=MSGF,
                decoyFile=''):

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
                                                 fdrOnly,
                                                 update,
                                                 outputDirectory,
                                                 massTolerance,
                                                 minPeptideLength,
                                                 maxPeptideLength,
                                                 maxDecoys,
                                                 precision,
                                                 reverse,
                                                 databaseType,
                                                 decoyFile)
  
  if not update:
    mergedDF = mergeDataFrames(denovoDF, decoyDF, databaseDF, qValue)

    if type(databaseDF) != type(''):
      peptideConversionDict = \
          createConversionDict(list(mergedDF[PEPTIDE]) + list(databaseDF[PEPTIDE]), acidConversionTable)
    else:
      peptideConversionDict = \
          createConversionDict(list(mergedDF[PEPTIDE]), acidConversionTable)
    
    with open(os.path.join(outputDirectory, 'conversionDictionary'), 'w') as f:
      f.write(str(peptideConversionDict))
  # Also need to write decoy selection, otherwised this won't match
  else:
    with open(os.path.join(outputDirectory, 'conversionDictionary')) as f:
      peptideConversionDict = eval(f.read())
  
  if qValue:
    databaseDF = addPeptideLength(databaseDF, peptideConversionDict)

  if not update:

    mergedDF = addPeptideLength(mergedDF, peptideConversionDict)
    mergedDF = mergedDF[mergedDF[LENGTH] >= minPeptideLength]
    mergedDF = mergedDF[mergedDF[LENGTH] <= maxPeptideLength]

    mergedDF = ScorePeptides.getPeptideScores(mergedDF,
                                              peptideConversionDict,
                                              acidMassTable,
                                              acidConversionTable,
                                              H2OMASSAdjusted,
                                              NH3MASSAdjusted,
                                              protonMassAdjusted,
                                              spectrumFileDirectory,
                                              scoringFunctionChoice,
                                              allPSSM,
                                              precision,
                                              minPeptideLength,
                                              maxPeptideLength,
                                              massTolerance,
                                              compression)
    
    mergedDF.to_csv(os.path.join(outputDirectory, 'scoredPeptides'))
  else:
    mergedDF = pd.read_csv(os.path.join(outputDirectory, 'scoredPeptides'))
  
  if not qValue and type(databaseDF) != type(''):
    denovoDF, decoyDF, databaseDF = separateDataFrames(mergedDF, qValue)
  else:
    denovoDF, decoyDF = separateDataFrames(mergedDF, qValue)
  
  # No FDR yet, but we're going to add it to this variable
  fdrDenovoDF = filterTopPeptides(denovoDF, scoreComparisionType)
  fdrDecoyDF = filterTopPeptides(decoyDF, scoreComparisionType)

  if type(databaseDF) != type(''):
    if not qValue:
      fdrDatabaseDF = filterTopPeptides(databaseDF, scoreComparisionType)
      fdrDatabaseDF.reset_index(drop=True, inplace=True)
    else:
      fdrDatabaseDF = databaseDF
  fdrDenovoDF.reset_index(drop=True, inplace=True)
  fdrDecoyDF.reset_index(drop=True, inplace=True)
  
  fdrDenovoDF = addFDR(fdrDenovoDF,
                       fdrDecoyDF,
                       scoreComparisionType,
                       precision,
                       increment,
                       fdrCalculationType)

  if not qValue and type(databaseDF) != type(''):
    fdrDatabaseDF = addFDR(fdrDatabaseDF,
                           fdrDecoyDF,
                           scoreComparisionType,
                           precision,
                           increment,
                           fdrCalculationType)
  
  if len(fdrDenovoDF) == 0:
    sys.exit()
  

  fdrDenovoDF, fdrCutoff = closestFDR(fdrDenovoDF, fdrCutoff, increment)
  if type(databaseDF) != type('') and len(fdrDatabaseDF) > 0:
    fdrDatabaseDF = fdrDatabaseDF[fdrDatabaseDF[FDR] >= fdrCutoff]

  ############# Do Analysis ####################

  outputFileName = os.path.join(outputDirectory, 'report.txt')

  with open(outputFileName, 'w') as f:
    uniquePeptidesDenovo = uniquePeptides(fdrDenovoDF)
    if type(databaseDF) != type('') and len(fdrDatabaseDF) > 0:
      uniquePeptidesDatabase = uniquePeptides(fdrDatabaseDF)
    f.write('FDR used: {}\n'.format(fdrCutoff))
    f.write('\n')
    f.write('Denovo Spectrum Matches: {}\n'.format(len(getSpectrumHits(fdrDenovoDF))))
    if type(databaseDF) != type('') and len(fdrDatabaseDF) > 0:
      f.write('Database Spectrum Matches: {}\n'.format(len(getSpectrumHits(fdrDatabaseDF))))
    f.write('\n')
    f.write('Denovo Unique Peptides Found: {}\n'.format(len(uniquePeptidesDenovo)))
    if type(databaseDF) != type('') and len(fdrDatabaseDF) > 0:
      f.write('Database Unique Peptides Found: {}\n'.format(len(uniquePeptidesDatabase)))
      f.write('Overlap: {}\n'.format(len(getOverlap(uniquePeptidesDenovo, \
                                                    uniquePeptidesDatabase))))
    f.write('\n')

    f.write('Length Distribution:\n')
    lengthCountDenovo = getLengthCountDict(fdrDenovoDF)
    lengthDistributionDenovo = \
      getLengthDistribution(lengthCountDenovo, uniquePeptidesDenovo)
    f.write('Denovo Length Distribution:\n')
    f.write(getLengthDistributionString(lengthCountDenovo, lengthDistributionDenovo))
    f.write('\n')

    if type(databaseDF) != type('') and len(fdrDatabaseDF) > 0:
      lengthCountDatabase = getLengthCountDict(fdrDatabaseDF)
      lengthDistributionDatabase = \
        getLengthDistribution(lengthCountDatabase, uniquePeptidesDatabase)
      f.write('Database Length Distribution:\n')
      f.write(getLengthDistributionString(lengthCountDatabase, lengthDistributionDatabase))
      f.write('\n')
    f.write('\n')


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


if __name__ == '__main__':

  arguments = parseArguments()

  getAnalysis(arguments.immunovo_results,
              arguments.database_results_dir,
              arguments.decoy_dir,
              arguments.spec_dir,
              arguments.qValue,
              arguments.fdr,
              arguments.spectrumComparison,
              arguments.scoreType,
              arguments.pssm_dir,
              arguments.acid_mass_file,
              arguments.output_dir,
              arguments.dataset_name,
              arguments.tsl,
              arguments.fdrCalculation,
              arguments.fdrOnly,
              arguments.update,
              arguments.increment,
              arguments.mmt,
              arguments.minP,
              arguments.maxP,
              arguments.decoys,
              arguments.prec,
              arguments.comp,
              arguments.reverse,
              arguments.database,
              arguments.decoy_file)