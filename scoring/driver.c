#include <stdio.h>
#include <stdlib.h>
#include "includes.h"
#include "inputFunctions.h"

int main(int argc, char *argv[])
{
  int i = 1;

  FILE* specFile = fopen(argv[i], "r"); i++;
  FILE* specDoubleFile = fopen(argv[i], "r"); i++;
  FILE* acidFile = fopen(argv[i], "r"); i++;
  FILE* massFile = fopen(argv[i], "r"); i++;
  FILE* scoreFile = fopen(argv[i], "r"); i++;
  FILE* resultsFile = fopen(argv[i], "w"); i++;
  int maxInstanceMiscleavage = (int) strtol(argv[i], NULL, 10); i++;
  int maxTotalMiscleavage = (int) strtol(argv[i], NULL, 10); i++;
  int maxCharge = (int) strtol(argv[i], NULL, 10); i++;
  int precision = (int) strtol(argv[i], NULL, 10); i++;
  double BPenalty = strtod(argv[i], NULL); i++;
  float massTolerance = strtof(argv[i], NULL); i++;
  int averageMassTolerance = (int) strtol(argv[i], NULL, 10); i++;
  int maximumMassTolerance = (int) strtol(argv[i], NULL, 10); i++;
  float lowestPenalty = strtof(argv[i], NULL); i++;
  long peptideMass = strtol(argv[i], NULL, 10); i++;
  long lowestSpecMass = strtol(argv[i], NULL, 10); i++;
  int minPeptideLength = (int) strtol(argv[i], NULL, 10); i++;
  int maxPeptideLength = (int) strtol(argv[i], NULL, 10); i++;
  long H2OMass = strtol(argv[i], NULL, 10); i++;
  long protonMass = strtol(argv[i], NULL, 10); i++;
  short yFirst = (short) strtol(argv[i], NULL, 10); i++;
  int binSize = (int) strtol(argv[i], NULL, 10); i++;

  SpectrumHash specHash;
  SpectrumHash specHashDouble;
  AminoAcids aminoAcids;
  AminoMasses aminoMasses;
  AminoScores aminoScoresArray[(maxPeptideLength - minPeptideLength) + 1];

  if (specFile == NULL) {
    printf("\n%s\n", argv[1]);
    printf("Invalid file for specFile\n");
    exit(0);
  }

  if (specDoubleFile == NULL) {
    printf("Invalid file for specDoubleFile\n");
    exit(0);
  }

  if (acidFile == NULL) {
    printf("Invalid file for acidFile\n");
    exit(0);
  }

  if (massFile == NULL) {
    printf("Invalid file for massFile\n");
    exit(0);
  }

  if (scoreFile == NULL) {
    printf("Invalid file for scoreFile\n");
    exit(0);
  }

  specHash = fileToSpectrumHash(specFile, precision);
  specHashDouble = fileToSpectrumHash(specDoubleFile, precision);

  aminoAcids = fileToAminoAcids(acidFile);
  aminoMasses = fileToAminoMasses(massFile, aminoAcids);

  AminoScores aminoScores;
  for (int i = (minPeptideLength - 1); i < maxPeptideLength; i++) {
    aminoScores = fileToAminoScores(scoreFile, aminoAcids, i+1);
    aminoScoresArray[i - (minPeptideLength-1)] = aminoScores;
  }

  fclose(specFile);
  fclose(specDoubleFile);
  fclose(acidFile);
  fclose(massFile);
  fclose(scoreFile);


  ViablePeptides viablePeptides = createViablePeptides(peptideMass,
                                                       specHash.size,
                                                       specHashDouble.size,
                                                       binSize,
                                                       precision);

  /* INPUT DIAGNOSTICS */

  /*

  printf("Amino Masses\n");

  for (int i = 0; i < aminoMasses.size; i++) {
    printf("%ld\n", aminoMasses.aminoAcid[i]);
  }
  printf("************************\n\n");

  printf("Amino Acids\n");
  for (int i = 0; i < aminoAcids.size; i++) {
    printf("%c\n", aminoAcids.array[i]);
  }
  printf("************************\n\n");

  int test = 0;
  scanf("%d", &test);
  */

  Results results = findBestPeptides(&viablePeptides,
                                     aminoMasses,
                                     aminoScoresArray,
                                     aminoAcids,
                                     specHash,
                                     specHashDouble,
                                     maxInstanceMiscleavage,
                                     maxTotalMiscleavage,
                                     maxCharge,
                                     precision,
                                     BPenalty,
                                     massTolerance,
                                     averageMassTolerance,
                                     maximumMassTolerance,
                                     lowestPenalty,
                                     peptideMass,
                                     lowestSpecMass,
                                     minPeptideLength,
                                     maxPeptideLength,
                                     H2OMass,
                                     protonMass,
                                     yFirst);

  LinkedList *refLink;

  for (int i = 0; i < results.singleResults.size; i++) {
    if (results.singleResults.bin[i].size != 0) {
      refLink = results.singleResults.bin[i].next;
      while (refLink != NULL) {
        fprintf(resultsFile, "%s\n", refLink->node.peptideString);
        fprintf(resultsFile, "%ld\n", refLink->node.mass);
        fprintf(resultsFile, "%f\n",
            massToleranceCalculation(refLink->node.mass, peptideMass));
        fprintf(resultsFile, "%lld\n",
          (long long) ((double) refLink->node.adjustedScore /
                       (double) refLink->node.pepLength));
        fprintf(resultsFile, "%lld\n", refLink->node.totalScore);
        fprintf(resultsFile, "%lld\n", refLink->node.massScore);
        fprintf(resultsFile, "%ld\n", refLink->node.aminoScore);
        refLink = refLink->next;
      }
    }
  }


  for (int i = 0; i < results.doubleResults.size; i++) {
    if (results.doubleResults.bin[i].size != 0) {
      refLink = results.doubleResults.bin[i].next;
      while (refLink != NULL) {
        fprintf(resultsFile, "%s\n", refLink->node.peptideString);
        fprintf(resultsFile, "%ld\n", refLink->node.mass);
        fprintf(resultsFile, "%f\n",
            massToleranceCalculation(refLink->node.mass, peptideMass));
        fprintf(resultsFile, "%lld\n",
          (long long) ((double) refLink->node.adjustedScore /
                       (double) refLink->node.pepLength));
        fprintf(resultsFile, "%lld\n", refLink->node.totalScore);
        fprintf(resultsFile, "%lld\n", refLink->node.massScore);
        fprintf(resultsFile, "%ld\n", refLink->node.aminoScore);
        refLink = refLink->next;
      }
    }
  }

  fclose(stdout);

  return 0;
}
