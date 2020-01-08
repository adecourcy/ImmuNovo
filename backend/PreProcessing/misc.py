import sys
# cheap hack until I can figure out how to do this properly in
# python 3.6
sys.path.append(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])


import backend.Structures.pssmIO as PSSMIO
import backend.Structures.acidMassTableIO as AcidMassTableIO
import backend.Structures.acidMassTable as AcidMassTable
import backend.Structures.pssm as PSSM
import backend.PreProcessing.acidConversion as AcidConversion
import backend.userInput as UserInput


def getAminoVariables(acidMassFile,
                      precision,
                      pssmDir,
                      minPepLength,
                      maxPepLength,
                      bEnd):

  acidMassTable = \
      AcidMassTable.adjustForPrecision(
          AcidMassTableIO.getAminoMasses(acidMassFile),
          precision)
  
  aminoAcids = AcidMassTable.getAcids(acidMassTable)

  allPSSM = \
      PSSM.adjustForPrecision(
          PSSMIO.getAllPSSM(pssmDir,
                            minPepLength,
                            maxPepLength,
                            aminoAcids,
                            bEnd),
          precision)

  conversionTable = \
      AcidConversion.createAcidConversionTable([acid for acid in acidMassTable])

  if conversionTable == {}:
    print("Currently this program only accepts up to 23 amino acids total\n")
    print("The mass file you have provided has more than 23 amino acids defined\n")
    print("Other files may also have too many masses but haven't been checked\n"
          "by the program at this time\n")
    exit()

  # quick and dirty way to make sure all our amino masses match for now
  UserInput.sanityCheck(allPSSM, acidMassTable)

  acidMassTable, allPSSM = \
      AcidConversion.convertAcidModifications(acidMassTable,
                                              allPSSM,
                                              conversionTable)
  
  return allPSSM, acidMassTable, conversionTable