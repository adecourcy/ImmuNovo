from decimal import Decimal

import argparse

# Keep the legacy, dictionary-based argument parser for dependent files

def tryFiles(files):
  cannotOpen = []
  for file in files:
    try:
      f = open(file, 'r')
    except:
      cannotOpen.append(file)

  if cannotOpen != []:
    print("Cannot open the following files.\n"
           "Please check the names and try again:")
    for file in cannotOpen:
      print(file)
    exit()

def printOptionalArguments():
  padding = '{:<20}{}'
  print("\nAll optional arguments must be passed as ARG=VALUE, in any order")
  print("Optional arguments:")
  print(str.format(padding, 'minP',
        'The minimum peptide length'))
  print(str.format(padding, 'maxP',
        'The maximum peptide length'))
  print(str.format(padding, 'PREC',
        'Decimal precision'))
  print(str.format(padding, 'OUT',
        'File output suffix'))
  print(str.format(padding, 'MMT',
        'The mass tolerance limit (maximum mass tolerance)'))
  print(str.format(padding, 'AMT',
        'The average (ideal) mass tolerance limit'))
  print('\nDon\'t touch these unless you know what you\'re doing:')

  print(str.format(padding, 'IMC',
        'Number of consectutive miscleavages before a peptide is discarded'))
  print(str.format(padding, 'TMC',
        'Number of total miscleavages before a peptide is discarded'))
  print(str.format(padding, 'BIN',
        'The maximum "bin" size of a spectrum mass'))
  print(str.format(padding, 'BPEN',
        'The score reduction if a b-ion, but not a y-ion, is matched'))
  print(str.format(padding, 'COMP',
        'The spectrum intensity log compression level'))
  print(str.format(padding, 'TOLP',
        'The maximum mass tolerance penalty'))



def parseArguments():
  parser = argparse.ArgumentParser()
  parser = parseDefaultArguments(parser)
  parser = parseOptionalArguments(parser)
  return parser.parse_args()


def parseDefaultArguments(parser):

  parser.add_argument('spec_dir',
      help='A directory containing 1 or more spectrum files')
  parser.add_argument('acid_mass_file',
      help='A file containing amino acid mass data')
  parser.add_argument('pssm_dir',
      help='A directory containing a positional scoring matrix')
  return parser


def parseOptionalArguments(parser):
  
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
  parser.add_argument('-prec', '--Precision',
      dest='prec',
      type=int,
      default=4,
      help='Decimial precision of our search')
  parser.add_argument('-os', '--Output-File-Suffix',
      dest='os',
      default='.out.csv',
      help='Suffix of the output file')
  parser.add_argument('-mmt', '--Maximum-Mass-Tolerance',
      dest='mmt',
      type=int,
      default=35,
      help='The maximum deviance between a peptide mass and spectrum mass '
      'before we consider the masses to not be matched (PPM)')
  parser.add_argument('-amt', '--Average-Mass-Tolerance',
      dest='amt',
      type=int,
      default=10,
      help='The deviance between a peptide and spectrum mass we would expect '
      'to see given instrument accuracy (PMM)')
  
  parser.add_argument('-imc', '--Intermediate-Miscleavage',
      dest='imc',
      type=int,
      default=1,
      help='Number of consectutive miscleavages until a peptide is discarded')
  parser.add_argument('-amc', '--Absolute-Miscleavage',
      dest='amc',
      type=int,
      default=3,
      help='Number of total miscleavages before a peptide is discarded')
  parser.add_argument('-bin', '--Number-of-Bins',
      dest='bin',
      type=int,
      default=5,
      help='The maximum "bin" size of a spectrum mass')
  parser.add_argument('-bPen', '--B-ion-Penalty',
      dest='bPen',
      type=float,
      default=0.5,
      help='The score reduction if a b-ion, but not a y-ion, is matched')
  parser.add_argument('-comp', '--Spectrum-Compression',
      dest='comp',
      type=int,
      default=2,
      help='The spectrum intensity log compression level')
  parser.add_argument('-tolp', '--Total-Mass-Tolerance-Penalty',
      dest='tolp',
      type=float,
      default=0.5,
      help='The maximum mass tolerance penalty')
  parser.add_argument('-db', '--debug',
      dest='db',
      type=bool,
      default=False,
      help='Run the program in debug mode')
      
  return parser

def parseParameterInput(args):

  defaultParameters = \
      {"IMC": 1,
       "TMC": 3,
       "MMT": 35, # maximum mass tolerance rate
       "AMT": 10, # average mass tolerance rate
       "minP": 9,
       "maxP": 12,
       "BIN": 5, # bin size
       "BPEN": Decimal(0.5), # b-ion penalty
       "COMP": 2, # log compression rate
       "TOLP": 0.5, # mass tolerance penalty
       "DEBUG": 0,
       "PREC": 4,
       "OUT": ".out.csv"}

  badArgs = []
  badValues = []

  for arg in args:
    try:
      argKey, argValue = arg.split("=")
    except ValueError:
      badArgs.append(arg)
      continue

    if argKey in defaultParameters:
      if argKey == "BPEN" or argKey == "TOLP" or argKey == 'COMP':

        try:
          converted = float(argValue)
        except ValueError:
          badValues.append((argKey, argValue))
          continue

        defaultParameters[argKey] = Decimal(converted)

      elif argKey == "OUT":
        defaultParameters[argKey] = argValue

      else:

        try:
          converted = int(argValue)
        except ValueError:
          badValues.append((argKey, argValue))
          continue

        defaultParameters[argKey] = converted

    else:
      badArgs.append(argKey)

  if len(badArgs) != 0 or len(badValues) != 0:

    if len(badArgs) != 0:
      print("Bad optional arguments given:")
      for b in badArgs:
        print(b)
      print("Call the script with no arguments for a list of argument options")

    if len(badValues) != 0:
      print("Bad values given for the following arguments: ")
      for b in badValues:
        print(b[0] + ": " + b[1])
      print("All argument values should be able to be converted to")
      print("integers, or float in the case of the mass tolerance argument")

    print("Exiting...\n")
    exit()

  minP = defaultParameters["minP"]
  maxP = defaultParameters["maxP"]

  if minP > maxP:
    print("The minimum entered peptide length is " + str(minP))
    print("The maximum entered peptide length is " + str(maxP))
    print("The mimimum peptide length must be less than the maximum length.")
    print("Switching minimum and maximum values")
    tmp = defaultParameters['minP']
    defaultParameters['minP'] = defaultParameters['maxP']
    defaultParameters['maxP'] = tmp

  return defaultParameters

def sanityCheck(allPSSM, acidMassTable):

  def notMatching(pssmName):
    print("Amino acids in mass table do not match amino acids in PSSM {}".format(pssmName))

  fail = False
  acids = [x for x in acidMassTable]

  for pssm in allPSSM:
    for length in allPSSM[pssm][1]:

      for acid in acids:
        if acid not in allPSSM[pssm][1][length]:
          fail = True
          notMatching(pssm)

      for acid in allPSSM[pssm][1][length]:
        if acid not in acids:
          fail = True
          notMatching(pssm)

  if fail:
    exit()