import argparse
import subprocess
import sys, os
# cheap hack until I can figure out how to do this properly in
# python 3.6
sys.path.append(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])

import pandas as pd
import backend.Structures.spectrumIO as SpectrumIO
import backend.Structures.spectrum as Spectrum
import backend.constants as Constants
import backend.PostProcessing.processResults as ProcessResults



def writeNNFormattedFile(nnResultsDF, tmpFileName):
  with open(tmpFileName, 'w') as f:
    f.write(nnResultsDF.to_csv(sep='\t'))


def writeFinal(mergedResults, finalName):
  with open(finalName, 'w') as f:
    f.write(mergedResults.to_csv())


def callNN(inputFileName, outputFileName, predfullLocation):
  # /l/glibc-2.17/lib/ld-2.17.so --library-path /l/glibc-2.17/lib:/usr/lib64:/lib64
  # /l/python3.6.6/bin/python predfull.py --input single.tsv --output single.mgf
  subprocess.run(['/l/glibc-2.17/lib/ld-2.17.so',
                  '--library-path',
                  '/l/glibc-2.17/lib:/usr/lib64:/lib64',
                  '/l/python3.6.6/bin/python',
                  os.path.join(predfullLocation, 'predfull.py'),
                  '--input',
                  inputFileName,
                  '--output',
                  outputFileName])


def getNNFormattedDF(reducedResultsDF, reducedSpectrums, collisionType):
  nnDF = pd.merge(reducedResultsDF,
                  reducedSpectrums,
                  how='left',
                  on=Constants.TITLE_SPECTRUM)
  nnDF['Type'] = [collisionType for i in range(len(nnDF))]
  nnDF['NCE'] = [25 for i in range(len(nnDF))]
  # Just in case the names don't match
  nnDF = nnDF.rename({Constants.PEPTIDE : 'Peptide'})
  nnDF = nnDF.drop(columns='index')
  
  return nnDF


def mergeCosineScores(resultsDF, reducedResultsDF, cosineScores, collisionType):
  if collisionType == 'HCD':
    nnScoreHeader = Constants.SCORE_HCD_NN
    nnCombinedHeader = Constants.SCORE_HCD_COMBINED
  elif collisionType == 'ETD':
    nnScoreHeader = Constants.SCORE_ETD_NN
    nnCombinedHeader = Constants.SCORE_ETD_COMBINED
  
  reducedResultsDF[nnScoreHeader] = cosineScores

  merged = pd.merge(resultsDF,
                    reducedResultsDF.drop(columns=[Constants.PEPTIDE, Constants.CHARGE]),
                    how='left',
                    on='index')
  
  merged[nnCombinedHeader] = \
      merged[nnScoreHeader] * merged[Constants.SCORE_PSSM]

  return merged


def getCosineScores(originalMGF, nnMGF, binSize, maxMassTolerance):
  allScores = []

  origSpectrumGenerator = SpectrumIO.getSpectrums(originalMGF)
  
  for nnSpectrum in SpectrumIO.getSpectrums(nnMGF):
    origSpectrum = next(origSpectrumGenerator)
    while Spectrum.getTitle(nnSpectrum) != Spectrum.getTitle(origSpectrum):
      origSpectrum = next(origSpectrumGenerator)
    cosineScore = \
      ProcessResults.cosineSimilarity(generateSpectralVector(origSpectrum, maxMassTolerance),
                                      generateSpectralVector(nnSpectrum, maxMassTolerance),
                                      binSize,
                                      'difference')
    allScores.append(cosineScore)
  return allScores


def generateSpectralVector(spectrum, maxMassTolerance):
  spectralVector = \
    [(float(x), float(y)) for x, y in zip(Spectrum.getMasses(spectrum),
                                          Spectrum.getIntensities(spectrum))]
  processedSpectralVector = \
    ProcessResults.normalize(ProcessResults.removeAdjacentPeaks(ProcessResults.keepTopKPeaks(spectralVector, 100),
                                  maxMassTolerance))
  
  return processedSpectralVector


def getSpectrumCharges(fileName):
  charges = []
  titles = []

  for spectrum in SpectrumIO.getSpectrums(fileName):
    charges.append(Spectrum.getCharge(spectrum))
    titles.append(Spectrum.getTitle(spectrum))
  
  return pd.DataFrame({Constants.TITLE_SPECTRUM: titles,
                       Constants.CHARGE: charges})


def reduceSpectrums(reducedResults, spectrums):
  results = \
    reducedResults[Constants.TITLE_SPECTRUM].drop_duplicates()

  return spectrums[spectrums[Constants.TITLE_SPECTRUM].isin(results)]


def parseResults(fileName):
  results = pd.read_csv(fileName)
  results['index'] = results.index
  return results


def reduceResults(results):
  results = results[[Constants.TITLE_SPECTRUM, Constants.PEPTIDE, 'index']]
  results = results[results[Constants.PEPTIDE] != Constants.NO_PEP]
  
  return results


def parseArguments():
  parser = argparse.ArgumentParser()
  parser.add_argument('results_file',
                      help='The file of results for scores to be added')
  parser.add_argument('mgf_file',
                      help='The MGF file with the spectrums for the results file')
  parser.add_argument('predfull_location',
                       help='The directory containing the predfull.py file')
  parser.add_argument('fragmentation_type',
                      help='The fragmentation type to run for the NN (HCD or ETD)',
                      choices=['HCD', 'ETD', 'BOTH'])
  parser.add_argument('-s', '--suffix',
                      dest='suffix',
                      help='The suffix of the output results file (defaults to -hcd, -etd, or -nn for the "both" option)')
  parser.add_argument('-d', '--Comparison-Distance',
                      default=0.1,
                      dest='dist',
                      help='The maximum comparison distance for the cosine correlation (daltons)')
  parser.add_argument('-mt', '--Mass-Tolerance',
                      default=20,
                      dest='mt',
                      help='The maximum mass tolerance for 2 m/z peaks to be considered the same (ppm)')
  
  return parser.parse_args()
  

if __name__ == '__main__':

  tmpNNInput = 'tmp.tsv'
  tmpNNOutput = 'tmp.mgf'
  
  arguments = parseArguments()

  results = parseResults(arguments.results_file)
  reducedResults = reduceResults(results)

  spectrumCharges = \
    reduceSpectrums(reducedResults, getSpectrumCharges(arguments.mgf_file))
  
  if arguments.fragmentation_type == 'BOTH':
    fragTypes = ['HCD', 'ETD']
    if arguments.suffix == '':
      arguments.suffix = '-nn'
  else:
    fragTypes = [arguments.fragmentation_type]
    if arguments.suffix == '':
      arguments.suffix = '-' + arguments.suffix
  
  suffix = arguments.suffix
  
  for fragmentationType in fragTypes:
    fragmentationType = arguments.fragmentation_type
    nnFormattedFile = \
      getNNFormattedDF(reducedResults, spectrumCharges, fragmentationType)
    writeNNFormattedFile(nnFormattedFile, tmpNNInput)

    callNN(tmpNNInput, tmpNNOutput, arguments.predfull_location)

    cosineScores = \
        getCosineScores(arguments.mgf_file,
                        tmpNNOutput,
                        arguments.dist,
                        arguments.mt)
    
    mergedDF = mergeCosineScores(results,
                                reducedResults,
                                cosineScores,
                                fragmentationType)
  
  finalFileName = \
    arguments.results_file.replace('.csv', arguments.suffix + '.csv')
  writeFinal(mergedDF, finalFileName)

