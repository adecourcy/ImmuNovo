#!/usr/bin/env python3
# Give a set of peptides, PSSMs, and associated Spectra, score peptides

# cheap hack until I can figure out how to do this properly in
# python 3.6
import sys, os
sys.path.append(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])

import argparse

import backend.userInput as userInput
from backend.constants import *


def parseArguments():
  parser = argparse.ArgumentParser()
  parser.add_argument('-o', '--Output-File',
                      dest='output_file',
                      default='scored_peptides.csv')
  parser.add_argument('peptide_file',
                        help='A csv of peptides to be scored')
  parser.add_argument('spectra_directory',
                        help='A directory of spectra to use in scoring')
  parser.add_argument('PSSM_directory',
                        help='A directory of PSSM files for scoring')
  parser.add_argument('acid_mass_file',
      help='A file containing mass spectrometry data')
  parser = userInput.parseOptionalArguments(parser)

  arguments = parser.parse_args()
  arguments.peptide_file = os.path.abspath(arguments.peptide_file)
  arguments.spectra_directory = os.path.abspath(arguments.spectra_directory)
  arguments.PSSM_directory = os.path.abspath(arguments.PSSM_directory)

  if not arguments.output_file.endswith('.csv'):
    arguments.output_file += '.csv'
  
  return arguments


if __name__ == '__main__':
  """
  pass a CSV with headers for spectra and peptides, headers should conform to
  the constants file standards

  Pass directory of spectra files, and a directory of PSSMs
  """
  parseArguments()
  pass