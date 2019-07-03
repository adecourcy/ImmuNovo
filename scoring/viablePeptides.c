#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "includes.h"

void printViablePeptides(ViablePeptides viablePeptides,
                         int debug)
{
  printf("Bin Size = %d\n", viablePeptides.binSize);
  printf("Spectrum Bins:\n");
  for (int i = 0; i < viablePeptides.spectrumBins.size; i++) {
    printf("********** Bin %d **********\n", i);
    printList(viablePeptides.spectrumBins.bin[i].next);
    printf("\n\n");
  }
  printf("Double Charge Bins:\n");
  for (int i = 0; i < viablePeptides.spectrumDoubleBins.size; i++) {
    printf("********** Bin %d **********\n", i);
    printList(viablePeptides.spectrumDoubleBins.bin[i].next);
    printf("\n\n");
  }
  if (debug != 1) {
    printf("Miscleavage Bins:\n");
    for (int i = 0; i < viablePeptides.miscleavageBins.size; i++) {
      printf("********** Bin %d **********\n", i);
      printList(viablePeptides.miscleavageBins.bin[i].next);
      printf("\n\n");
    }
  }
}

SpectrumBin createMiscleavageBins(long largestSpectrum,
                                  int precision)
{
  int arraySize = (int) largestSpectrum / pow(10, precision) + 1;
  Head *bin = (Head*) GC_MALLOC(arraySize * sizeof(Head));
  //Head *bin = (Head*) malloc(arraySize * sizeof(Head));
  for (int i = 0; i < arraySize; i++) {
    bin[i] = newHead();
  }
  SpectrumBin spectrumBin;
  spectrumBin.size = arraySize;
  spectrumBin.bin = bin;
  return spectrumBin;
}

SpectrumBin createSpectrumBins(int spectrumLength)
{
  SpectrumBin spectrumBin;
  spectrumBin.bin = (Head*) GC_MALLOC(spectrumLength * sizeof(Head));
  //spectrumBin.bin = (Head*) malloc(spectrumLength * sizeof(Head));
  spectrumBin.size = spectrumLength;
  for (int i=0; i < spectrumLength; i++) {
    spectrumBin.bin[i] = newHead();
  }

  return spectrumBin;
}

ViablePeptides createViablePeptides(long largestSpectrum,
                                    int singleLength,
                                    int doubleLength,
                                    int binSize,
                                    int precision)
{
  ViablePeptides viablePeptides;
  viablePeptides.binSize = binSize;
  viablePeptides.spectrumBins = createSpectrumBins(singleLength);
  viablePeptides.spectrumDoubleBins = createSpectrumBins(doubleLength);
  viablePeptides.miscleavageBins = createMiscleavageBins(largestSpectrum,
                                                         precision);
  return viablePeptides;
}


void destroySpectrumBin(SpectrumBin spectrumBin)
{
  for (int i = 0; i < spectrumBin.size; i++) {
    destroyLinkedList(&spectrumBin.bin[i]);
  }
}


void destroyViablePeptides(ViablePeptides viablePeptides)
{
  destroySpectrumBin(viablePeptides.spectrumBins);
  destroySpectrumBin(viablePeptides.spectrumDoubleBins);
  destroySpectrumBin(viablePeptides.miscleavageBins);
}


int getMiscleavageBinNumber(Node node, int precision)
{
  return (int) node.mass / pow(10, precision);
}

void addToMiscleavageBin(ViablePeptides* viablePeptides,
                         Node node,
                         int precision)
{
  int binNumber = getMiscleavageBinNumber(node, precision);

  if (binNumber < viablePeptides->miscleavageBins.size) {
    addLink(&viablePeptides->miscleavageBins.bin[binNumber], node);
  }
}


long long calculateAdjustedScore(Node node,
                                 int averageMassTolerance,
                                 int maximumMassTolerance,
                                 float lowestPenalty)
{
  float adjustment =
        ((((lowestPenalty - 1.0) /
         ((float) ((node.pepLength * maximumMassTolerance) -
                   (node.pepLength * averageMassTolerance)))) *
                    node.massTolerance) + 1.0);

  return (long long) ((float) node.totalScore * adjustment);
}

void addToSpectrumBinHelper(Head *head,
                            RejectedDatabase *rejectedDatabase,
                            Node node,
                            int charge,
                            int binNumber,
                            int binSize)
{

  addLink(head, node);

  if (head->size > binSize) {
    Node rejectedNode = removeLowest(head);

    addToRejectedDatabase(rejectedDatabase,
                          rejectedNode,
                          binNumber,
                          charge);
  }

}


int trySpectrumBin(SpectrumBin spectrumBin,
                   RejectedDatabase *rejectedDatabase,
                   SpectrumHash spectrumHash,
                   Node node,
                   long mass,
                   double BPenalty,
                   float massTolerance,
                   int averageMassTolerance,
                   int maximumMassTolerance,
                   float lowestPenalty,
                   int charge,
                   int binSize,
                   short bCompliment,
                   int precision,
                   short doubleCharge)
{
  float massToleranceResult;
  long long bAdjustedScore;
  SpectrumNode lookupNode = lookupMass(mass,
                                       precision,
                                       spectrumHash,
                                       massTolerance,
                                       &massToleranceResult);

  if (lookupNode.mass == 0) {
    return 0;
  } else {
    if (doubleCharge == 1) {
      node.singleCharge = 0;
    }
    if (bCompliment == 1) {
    bAdjustedScore =
          (long long) (BPenalty * (long double) lookupNode.score);
    } else {
      bAdjustedScore = lookupNode.score;
    }
    node.totalScore += bAdjustedScore;
    node.massTolerance += massToleranceResult;
    node.massScore += bAdjustedScore;

    long long adjustedScore =
        calculateAdjustedScore(node,
                               averageMassTolerance,
                               maximumMassTolerance,
                               lowestPenalty);

    node.currentMiscleavage = 0;
    node.adjustedScore = adjustedScore;

    addToSpectrumBinHelper(&spectrumBin.bin[lookupNode.positionNumber],
                           rejectedDatabase,
                           node,
                           charge,
                           lookupNode.positionNumber,
                           binSize);

    return 1;
  }

}


int addToSpectrumBin(ViablePeptides* viablePeptides,
                     RejectedDatabase* rejectedDatabase,
                     SpectrumHash spectrumHash,
                     SpectrumHash spectrumHashDouble,
                     Node node,
                     double BPenalty,
                     int precision,
                     float massTolerance,
                     int averageMassTolerance,
                     int maximumMassTolerance,
                     float lowestPenalty,
                     long peptideMass,
                     int charge,
                     long protonMass,
                     long lowestSpecMass)
{
  int result;
  long mass = node.mass;
  short bCompliment = 0;

  for (int i = 0; i < 2; i++) {

    if (node.singleCharge == 1) {

      result = trySpectrumBin(viablePeptides->spectrumBins,
                              rejectedDatabase,
                              spectrumHash,
                              node,
                              mass,
                              BPenalty,
                              massTolerance,
                              averageMassTolerance,
                              maximumMassTolerance,
                              lowestPenalty,
                              charge,
                              viablePeptides->binSize,
                              bCompliment,
                              precision,
                              0);

      if (result == 1) {

        return result;
      }
    }

    if (node.doubleCharge == 1) {

      result = trySpectrumBin(viablePeptides->spectrumDoubleBins,
                              rejectedDatabase,
                              spectrumHashDouble,
                              node,
                              mass,
                              BPenalty,
                              massTolerance,
                              averageMassTolerance,
                              maximumMassTolerance,
                              lowestPenalty,
                              charge,
                              viablePeptides->binSize,
                              bCompliment,
                              precision,
                              1);

      if (result == 1) {
        return result;
      }
    }

    if (BPenalty != 0.0) {
      mass = findBCompliment(peptideMass, node.mass, protonMass);
      bCompliment = 1;
    } else {
      break;
    }

    if (mass < lowestSpecMass) {
      return 0;
    }

  }

  return 0;

}

long findBCompliment(long peptideMass, long nodeMass, long protonMass)
{
  return (peptideMass - nodeMass) + protonMass;
}



void addToViablePeptides(Node node,
                         ViablePeptides* viablePeptides,
                         RejectedDatabase* rejectedDatabase,
                         SpectrumHash spectrumHash,
                         SpectrumHash spectrumHashDouble,
                         int maxInstanceMiscleavage,
                         int maxTotalMiscleavage,
                         int maxCharge,
                         double BPenalty,
                         int precision,
                         float massTolerance,
                         int averageMassTolerance,
                         int maximumMassTolerance,
                         float lowestPenalty,
                         long peptideMass,
                         int charge,
                         long protonMass,
                         long lowestSpecMass)
{
  short addedCharge = 0;

  if (node.charge == 0) {
    node.charge = 1;
    node.mass += protonMass;
    addedCharge = 1;
  }

  int result = addToSpectrumBin(viablePeptides,
                                rejectedDatabase,
                                spectrumHash,
                                spectrumHashDouble,
                                node,
                                BPenalty,
                                precision,
                                massTolerance,
                                averageMassTolerance,
                                maximumMassTolerance,
                                lowestPenalty,
                                peptideMass,
                                charge,
                                protonMass,
                                lowestSpecMass);

  if ((result == 0) &&
      (node.currentMiscleavage < maxInstanceMiscleavage) &&
      (node.totalMiscleavage < maxTotalMiscleavage)) {

    if (node.peptideString[node.pepLength - 1] == 'P') {
      node.miscleavageSwap = 1;
    }

    node.currentMiscleavage += 1;
    node.totalMiscleavage += 1;
    node.massTolerance += maximumMassTolerance;
    
    if (addedCharge == 1) {
      node.charge = 0;
      node.mass -= protonMass;
    }

    addToMiscleavageBin(viablePeptides, node, precision);
  }
}


ViablePeptides getNextGeneration(ViablePeptides* viablePeptides,
                                 RejectedDatabase* rejectedDatabase,
                                 AminoMasses aminoMasses,
                                 AminoScores aminoScores,
                                 AminoAcids aminoAcids,
                                 SpectrumHash spectrumHash,
                                 SpectrumHash spectrumHashDouble,
                                 int maxInstanceMiscleavage,
                                 int maxTotalMiscleavage,
                                 int maxCharge,
                                 int precision,
                                 double BPenalty,
                                 float massTolerance,
                                 int averageMassTolerance,
                                 int maximumMassTolerance,
                                 float lowestPenalty,
                                 long peptideMass,
                                 long protonMass,
                                 long lowestSpecMass)
{
  LinkedList *link;
  Node *nodeList;

  long miscleavageBinSize =
      (viablePeptides->miscleavageBins.size - 1) * pow(10, precision);

  ViablePeptides newViablePeptides =
      createViablePeptides(miscleavageBinSize,
                           viablePeptides->spectrumBins.size,
                           viablePeptides->spectrumDoubleBins.size,
                           viablePeptides->binSize,
                           precision);

  for (int i = 0; i < viablePeptides->spectrumBins.size; i++) {
    link = viablePeptides->spectrumBins.bin[i].next;
    while (link != NULL) {
      nodeList = branchPeptide(link->node,
                               aminoScores,
                               aminoMasses,
                               aminoAcids);

      for (int k = 0; k < aminoAcids.size; k++) {
        addToViablePeptides(nodeList[k],
                            &newViablePeptides,
                            rejectedDatabase,
                            spectrumHash,
                            spectrumHashDouble,
                            maxInstanceMiscleavage,
                            maxTotalMiscleavage,
                            maxCharge,
                            BPenalty,
                            precision,
                            massTolerance,
                            averageMassTolerance,
                            maximumMassTolerance,
                            lowestPenalty,
                            peptideMass,
                            nodeList[k].charge,
                            protonMass,
                            lowestSpecMass);
      }
      //free(nodeList);

      link = link->next;
    }
  }

  for (int i = 0; i < viablePeptides->spectrumDoubleBins.size; i++) {
    link = viablePeptides->spectrumDoubleBins.bin[i].next;
    while (link != NULL) {
      nodeList = branchPeptide(link->node,
                               aminoScores,
                               aminoMasses,
                               aminoAcids);

      for (int k = 0; k < aminoAcids.size; k++) {
        addToViablePeptides(nodeList[k],
                            &newViablePeptides,
                            rejectedDatabase,
                            spectrumHash,
                            spectrumHashDouble,
                            maxInstanceMiscleavage,
                            maxTotalMiscleavage,
                            maxCharge,
                            BPenalty,
                            precision,
                            massTolerance,
                            averageMassTolerance,
                            maximumMassTolerance,
                            lowestPenalty,
                            peptideMass,
                            nodeList[k].charge,
                            protonMass,
                            lowestSpecMass);
      }
      //free(nodeList);

      link = link->next;
    }
  }

  for (int i = 0; i < viablePeptides->miscleavageBins.size; i++) {
    link = viablePeptides->miscleavageBins.bin[i].next;
    while (link != NULL) {
      nodeList = branchPeptide(link->node,
                               aminoScores,
                               aminoMasses,
                               aminoAcids);

      for (int k = 0; k < aminoAcids.size; k++) {
        addToViablePeptides(nodeList[k],
                            &newViablePeptides,
                            rejectedDatabase,
                            spectrumHash,
                            spectrumHashDouble,
                            maxInstanceMiscleavage,
                            maxTotalMiscleavage,
                            maxCharge,
                            BPenalty,
                            precision,
                            massTolerance,
                            averageMassTolerance,
                            maximumMassTolerance,
                            lowestPenalty,
                            peptideMass,
                            nodeList[k].charge,
                            protonMass,
                            lowestSpecMass);
      }
      //free(nodeList);

      link = link->next;
    }
  }

  destroyViablePeptides(*viablePeptides);
  return newViablePeptides;

}
