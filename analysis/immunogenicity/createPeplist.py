import pandas as pd
import string
import random

import subprocess
import os

# Let's keep it simple from know and just append the output of MixMHCpred
# to the outputfile

# This is really rough, and basically assumes that all of the output files
# are in the same, hard-coded directory

# Assume the MixMHCpred script is in our directory


def tmpGetAllele(pssmTitle: str) -> str:
  """
    A temporary function to extract the allele label from the PSSM title in the
    output file (assuming the allele is in the PSSM title)
    
    This function assumes the following format:
      "HLA-[AlleleName] [etc].pssm"
    
    With an example AlleleName "A*55:02"
  """

  unformattedAllele = pssmTitle.split()[0].split(sep='-')[1]
  convertedAlleleTitle = \
                unformattedAllele[0] + \
                unformattedAllele[2:4] + \
                unformattedAllele[5:7]
  
  return convertedAlleleTitle


def replaceAminoAcidModifications(peptide):
  return ''.join([x for x in peptide if x.isalpha()])


def extractPeptides(df):  
  # Find all peptides, and write them to a temporary file
  allPeps = df['Peptide'].to_csv(index=False).strip().split('\n')
  return '\n'.join([replaceAminoAcidModifications(x) for x in allPeps])


def extractAlleles(df):
  # Find all unique PSSM titles, extract allele entries from them, and create
  # an argument string for MixMHCpred
  allAlleles = \
    [tmpGetAllele(x) for x in \
         df['PSSM Name'].drop_duplicates().to_csv(index=False).strip().split('\n')]
  return ','.join(allAlleles)


def setUpData(dataFrame, temporaryFileName):
  peptides = extractPeptides(dataFrame)
  alleles = extractAlleles(dataFrame)
  with open(temporaryFileName, 'w') as f:
    f.write(peptides)
  return alleles


def runMixMHC(alleles, inputFileName, outFileName):
  subprocess.run([os.path.join('./' 'MixMHCpred', 'MixMHCpred'),
                        '-i', inputFileName,
                        '-o', outFileName,
                        '-a', alleles])


def createImmunovoDataFrame(directoryName):
  dataFrameList = []
  for file in os.listdir(directoryName):
    dataFrameList.append(pd.read_csv(os.path.join(directoryName, file)))
  df = pd.concat(dataFrameList)
  return df[df.Peptide != "No Peptide Found"]


def filterTitle(fileData):
  commentIndex = 0
  while '#' in fileData[commentIndex]:
    commentIndex += 1
  return fileData[commentIndex:]


def reformatMixData(fileData):
  def reformatLine(line):
    return ','.join(line.split())
  return [reformatLine(x) for x in fileData]


def createMixDataFrame(fileName):
  with open(fileName, 'r') as f:
    fileData = filterTitle(f.read().strip().split('\n'))
  with open(fileName, 'w') as f:
    f.write('\n'.join(reformatMixData(fileData)))
  return pd.read_csv(fileName)


def getUniqueFileName(numChars=10):
  # Just try random file names until a non-existant one is found
  name = ''
  for i in range(numChars):
    name += random.choice(string.ascii_letters)
  if os.path.isfile(name):
    return getUniqueFileName(numChars)
  else:
    return name



if __name__ == "__main__":
  inputFileName = getUniqueFileName()
  with open(inputFileName, 'w') as f:
    f.write('')
  outputFileName = getUniqueFileName()
  with open(outputFileName, 'w') as f:
    f.write('')

  ImmunovoDataFrame = createImmunovoDataFrame('outputFiles')

  alleles = setUpData(ImmunovoDataFrame, inputFileName)
  runMixMHC(alleles, inputFileName, outputFileName)
  mixedDataFrame = createMixDataFrame(outputFileName)

  mergedDF = pd.concat([ImmunovoDataFrame, mixedDataFrame])

  os.remove(inputFileName)
  os.remove(outputFileName)

  with open('data_out.csv', 'w') as f:
    f.write(mergedDF.to_csv(index=False))


  














  