#!/usr/bin/env python3
import sys, os

import backend.Structures.pssmIO as PSSMIO
import backend.Structures.acidMassTableIO as AcidMassTableIO
import backend.Structures.spectrumIO as SpectrumIO
import backend.Structures.pssm as PSSM
import backend.Structures.acidMassTable as AcidMassTable
import backend.Structures.spectrum as Spectrum
import backend.PreProcessing.acidConversion as AcidConversion
import backend.PreProcessing.spectrumConversion as SpectrumConversion
import backend.PostProcessing.processResults as ProcessResults
import backend.PreProcessing.misc as Misc
import analysis.analyzeResults as AnalyzeResults
from backend.constants import *
from backend.userInput import *

def writeToFile(acidTable,
                pssm,
                spectrumMasses,
                spectrumMassesDouble,
                spectrumScores,
                spectrumScoresDouble,
                aminoAcids,
                specFileName,
                specDoubleFileName,
                acidFileName,
                massFileName,
                scoreFileName):

  specFile = open(specFileName, "w")
  specDoubleFile = open(specDoubleFileName, "w")
  acidFile = open(acidFileName, "w")
  massFile = open(massFileName, "w")
  scoreFile = open(scoreFileName, "w")

  writeLine = ""
  specFile.write(str(len(spectrumMasses)) + "\n")
  for mass, score in zip(spectrumMasses, spectrumScores):
    writeLine += str(mass) + " " + str(score) + "\n"
  specFile.write(writeLine.strip())

  writeLine = ""
  specDoubleFile.write(str(len(spectrumMassesDouble)) + "\n")
  for mass, score in zip(spectrumMassesDouble, spectrumScoresDouble):
    writeLine += str(mass) + " " + str(score) + "\n"
  specDoubleFile.write(writeLine.strip())

  writeLine = ""
  acidFile.write(str(len(acidTable)) + "\n")
  for acid in aminoAcids:
    if acid in acidTable:
      writeLine += (acid + "\n")
  acidFile.write(writeLine.strip())

  writeLine = ""
  for acid in aminoAcids:
    if acid in acidTable:
      writeLine += str(acidTable[acid]) + "\n"
  massFile.write(writeLine.strip())

  writeLine = ""
  for length in pssm:
    for acid in aminoAcids:
      if acid in pssm[length]:
        for prob in pssm[length][acid]:
          writeLine += str(prob) + " "
        writeLine = writeLine.strip() + "\n"
  scoreFile.write(writeLine.strip())

  specFile.close()
  specDoubleFile.close()
  acidFile.close()
  massFile.close()
  scoreFile.close()


def programUsageOutput():
  print("\nUsage: main.py SPECDIR ACIDMASSFILE PSSMDIR [OPT ARGS]\n\n")

  padding = '{:<20}{}'
  print("Positional arguments:")
  print(str.format(padding, 'SPECDIR',
        'A directory containing 1 or more spectrum files'))
  print(str.format(padding, 'ACIDMASSFILE',
        'A file containing mass spectrometry data'))
  print(str.format(padding, 'PSSMDIR',
        'A directory containing a positional scoring matrix'))

  printOptionalArguments()
  exit()


def improperNumArgumentsOutput():
  print("\nThis script takes a minimum of 3 arguments,")
  print("with up to 13 additional, optional arguments")
  print("Call this script with no arguments for usage details")
  print("Exiting...\n")
  exit()

def clean(fileList):

  for f in fileList:
    try:
      os.remove(f)
    except:
      continue


def checkInputLength(args):
  if len(args) == 1:
    programUsageOutput()

  if len(args) < 4 or len(args) > 17:
    improperNumArgumentsOutput()


if __name__ == '__main__':

  programDir = os.path.dirname(os.path.abspath(__file__))

  cPath = os.path.join(os.path.join(programDir, "scoring"))
  specFile = os.path.join(cPath, "specFile")
  specDoubleFile = os.path.join(cPath, "specDoubleFile")
  acidFile = os.path.join(cPath, "acidFile")
  massFile = os.path.join(cPath, "massFile")
  scoreFile = os.path.join(cPath, "scoreFile")
  resultsFile = os.path.join(cPath, "results")
  aminoAcids = ["G", "A", "S", "P", "V",
                "T", "L", "I", "N", "D",
                "Q", "K", "E", "M", "H",
                "F", "R", "C", "Y", "W",
                "B", "J", "O"]

  fileList = [specFile, specDoubleFile, acidFile,
              massFile, scoreFile, resultsFile]

  try:
    clean(fileList)
  except:
    pass

  arguments = parseArguments(os.path.dirname(os.path.abspath(os.path.realpath(__file__))))
  allPSSM, acidMassTable, conversionTable = \
                              Misc.getAminoVariables(arguments.acid_mass_file,
                                                     arguments.prec,
                                                     arguments.pssm_dir,
                                                     arguments.minP,
                                                     arguments.maxP,
                                                     arguments.bEnd)


  H2OMassAdjusted = int(H2OMASS * (10**arguments.prec))
  NH3MassAdjusted = int(NH3MASS * (10**arguments.prec))
  protonMassAdjusted = int(PROTONMASS * (10**arguments.prec))


  for spectrumFileName in os.listdir(arguments.spec_dir):
    spectrumFile = os.path.join(arguments.spec_dir, spectrumFileName)

    try:
      outputFile = \
          open(os.path.join(arguments.output_dir, spectrumFileName + arguments.os), 'w')
    except:
      print(str.format("Error opening {} to write output, skpping {} spectrum file",
                        spectrumFileName + arguments.os,
                        spectrumFileName))
      continue

    outputFile.write("Filename,{},{},{},"
                     "Precursor Mass,{}\n".format(TITLE_SPECTRUM,
                                                  TITLE_PSSM,
                                                  PEPTIDE,
                                                  PRECURSOR_ERROR))

    for spectrum in SpectrumIO.getSpectrums(spectrumFile):
      spectrumMasses, spectrumMassesDouble, \
      spectrumIntensities, spectrumIntensitiesDouble = \
          SpectrumConversion.adjustSpectrumPrecision(
              *SpectrumConversion.processSpectrum(Spectrum.getMasses(spectrum),
                                                  Spectrum.getIntensities(spectrum),
                                                  Spectrum.getPrecursorMass(spectrum),
                                                  H2OMASS,
                                                  PROTONMASS,
                                                  Spectrum.getCharge(spectrum),
                                                  arguments.comp),
              arguments.prec)


      for pssmTitle in allPSSM:

        for acid in acidMassTable:
          if acid not in aminoAcids:
            aminoAcids.append(acid)

        writeToFile(acidMassTable,
                    allPSSM[pssmTitle][1],
                    spectrumMasses,
                    spectrumMassesDouble,
                    spectrumIntensities,
                    spectrumIntensitiesDouble,
                    aminoAcids,
                    specFile,
                    specDoubleFile,
                    acidFile,
                    massFile,
                    scoreFile)

        oldCWD = os.getcwd()
        os.chdir(programDir)

        cOutput = os.system("./scoring/findPep {} {} {} {} {} \
                                     {} {} {} {} {} \
                                     {} {} {} {} {} \
                                     {} {} {} {} {} \
                                     {} {} {}".format(specFile,
                                                      specDoubleFile,
                                                      acidFile,
                                                      massFile,
                                                      scoreFile,
                                                      resultsFile,
                                                      arguments.imc,
                                                      arguments.amc,
                                                      Spectrum.getCharge(spectrum) - 1,
                                                      arguments.prec,
                                                      float(arguments.bPen),
                                                      arguments.mmt,
                                                      arguments.amt,
                                                      arguments.mmt,
                                                      arguments.tolp,
                                                      spectrumMasses[-1],
                                                      spectrumMasses[0],
                                                      arguments.minP,
                                                      arguments.maxP,
                                                      H2OMassAdjusted,
                                                      protonMassAdjusted,
                                                      1,
                                                      arguments.bin))

        if cOutput != 0:
          clean(fileList)
          exit()
      
        os.chdir(oldCWD)

        experimentalSpectrum = spectrumMasses + spectrumMassesDouble
        experimentalIntensities = spectrumIntensities + spectrumIntensitiesDouble
        results = ProcessResults.processResults(resultsFile,
                                                acidMassTable,
                                                experimentalSpectrum,
                                                experimentalIntensities,
                                                protonMassAdjusted,
                                                H2OMassAdjusted,
                                                int((NH3MASS * (10 ** arguments.prec))),
                                                arguments.mmt,
                                                arguments.prec,
                                                conversionTable)


        spectrumTitle = Spectrum.getTitle(spectrum)
        if len(results) == 0:
          outputString = spectrumFileName + "," + spectrumTitle + \
                         ',' + pssmTitle + ",{}".format(NO_PEP) + "\n"
          outputFile.write(outputString)

        else:
          ignoredLengths = PSSM.getIgnoredLengths(pssmTitle, allPSSM)

          for node in results:
            if len(node.peptideString) in ignoredLengths:
              continue

            outputString = \
              str.format('{},{},{},{},{},{}',
                         spectrumFileName,
                         spectrumTitle,
                         pssmTitle,
                         node.peptideString,
                         Spectrum.getPrecursorMass(spectrum),
                         node.precursorError)

            outputFile.write(outputString)

        clean(fileList)

    outputFile.close()
  
  analyzeOutputDir = os.path.join(arguments.output_dir, 'reports')
  if not os.path.exists(analyzeOutputDir):
    os.makedirs(analyzeOutputDir)
  AnalyzeResults.getAnalysis(os.path.join(arguments.output_dir, spectrumFileName + arguments.os),
                             arguments.database_results_dir,
                             arguments.decoy_dir,
                             arguments.spec_dir,
                             arguments.qValue,
                             arguments.fdr,
                             arguments.spectrumComparison,
                             arguments.scoreType,
                             arguments.pssm_dir,
                             arguments.acid_mass_file,
                             analyzeOutputDir,
                             arguments.dataset_name,
                             arguments.tsl,
                             arguments.fdrCalculation,
                             arguments.fdrOnly,
                             False,
                             arguments.increment,
                             arguments.mmt,
                             arguments.minP,
                             arguments.maxP,
                             arguments.decoys,
                             arguments.prec,
                             arguments.comp,
                             arguments.reverse,
                             arguments.database)