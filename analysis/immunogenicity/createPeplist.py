import pandas as pd
import tempfile

# Let's keep it simple from know and just append the output of MixMHCpred
# to the outputfile


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


def extractPeptides(df, tmpFileName):
  # Find all unique peptides, and write them to a temporary file
  df[df.Peptide != "No Peptide Found"].Peptide.drop_duplicates().to_csv(tmpFileName, index=False)


def extractAlleles(df):
  # Find all unique PSSM titles, extract allele entries from them, and create
  # an argument string for MixMHCpred
  allAlleles = \
    [tmpGetAllele(x) for x in \
         df['PSSM Name'].drop_duplicates().to_csv(index=False).split()]
  return allAlleles














  