#include <stdlib.h>
#include "includes.h"

Node rescoreNode(Node node,
                 AminoScores aminoScores,
                 int averageMassTolerance,
                 int maximumMassTolerance,
                 float lowestPenalty)
{
  node.aminoScore = 0;

  for (int i = 0; i < node.pepLength; i++) {
    node.aminoScore += getAminoScore(aminoScores,
                                     node.peptideString[i],
                                     i);
  }

  node.totalScore = node.aminoScore + node.massScore;
  node.adjustedScore = calculateAdjustedScore(node,
                                              averageMassTolerance,
                                              maximumMassTolerance,
                                              lowestPenalty);

  return node;

}


void addNodesToViablePeptides(ViablePeptides *viablePeptides,
                              RejectedDatabase *rejectedDatabase,
                              RejectedLinkedList *nodes)
{

  RejectedLinkedList *listHead = nodes;
  Head *binHead;
  Node node;
  int binNumber;
  int charge;

  while (nodes != NULL) {
    node = copyNode(nodes->node);
    //nodes->node = createNode(0);
    binNumber = nodes->binNumber;
    charge = nodes->charge;

    if (node.singleCharge == 1) {
      binHead = &viablePeptides->spectrumBins.bin[binNumber];
    } else {
      binHead = &viablePeptides->spectrumDoubleBins.bin[binNumber];
    }

    addToSpectrumBinHelper(binHead,
                           rejectedDatabase,
                           node,
                           charge,
                           binNumber,
                           viablePeptides->binSize);

    nodes = nodes->next;

  }

  destroyRejectedList(listHead);
}


void iterateRejectedPeptides(RejectedDatabase* rejectedDatabase,
                             ViablePeptides* viablePeptides,
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
                             long lowestSpecMass,
                             int previousPeptideLength)
{
  RejectedLinkedList *rejectedPeptides;
  RejectedLinkedList *head;
  for (int i = 0; i < previousPeptideLength; i++) {

    head = rejectedPeptides = rejectedDatabase->array[i];
    rejectedDatabase->array[i] = NULL;


    if (rejectedPeptides != NULL) {

      while (rejectedPeptides != NULL) {
        rejectedPeptides->node = rescoreNode(rejectedPeptides->node,
                                             aminoScores,
                                             averageMassTolerance,
                                             maximumMassTolerance,
                                             lowestPenalty);

        rejectedPeptides = rejectedPeptides->next;
      }

      rejectedPeptides = head;

      /*
      #include <stdio.h>
      printViablePeptides(*viablePeptides, 1);
      printf("#####################################\n");
      printRejectedDatabase(*rejectedDatabase);
      printf("#####################################\n");
      while (rejectedPeptides != NULL) {
        printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
        printf("node:\n");
        printNode(rejectedPeptides->node);
        printf("\n");
        printf("bin: %d\n", rejectedPeptides->binNumber);
        printf("charge %d\n", rejectedPeptides->charge);

        rejectedPeptides = rejectedPeptides->next;
      }
      rejectedPeptides = head;
      */
      /*
      #include <stdio.h>
      if (previousPeptideLength == 10 && i == 8) {
        //printViablePeptides(*viablePeptides, 0);
        //printRejectedDatabase(*rejectedDatabase);
        printRejectedList(rejectedPeptides);
        exit(0);
      }
      */

      addNodesToViablePeptides(viablePeptides,
                               rejectedDatabase,
                               rejectedPeptides);

      /*
      #include <stdio.h>
      printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
      printViablePeptides(*viablePeptides, 1);
      */
      
      *viablePeptides = getNextGeneration(viablePeptides,
                                          rejectedDatabase,
                                          aminoMasses,
                                          aminoScores,
                                          aminoAcids,
                                          spectrumHash,
                                          spectrumHashDouble,
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
                                          protonMass,
                                          lowestSpecMass);
      
      /*
      printf("###############################################\n");
      printViablePeptides(*viablePeptides, 1);
      */
      
    }
  }

  head = rejectedPeptides =
      rejectedDatabase->array[previousPeptideLength];
  rejectedDatabase->array[previousPeptideLength] = NULL;
  if (rejectedPeptides != NULL) {
    rejectedPeptides->node = rescoreNode(rejectedPeptides->node,
                                         aminoScores,
                                         averageMassTolerance,
                                         maximumMassTolerance,
                                         lowestPenalty);

    rejectedPeptides = rejectedPeptides->next;
  }

  rejectedPeptides = head;

  addNodesToViablePeptides(viablePeptides,
                           rejectedDatabase,
                           rejectedPeptides);

}


void addAllToRejectedDatabase(ViablePeptides* viablePeptides,
                              RejectedDatabase* rejectedDatabase)
{
  Head head;
  LinkedList *link;
  Node node;

  for (int i = 0; i < viablePeptides->spectrumBins.size; i++) {
    head = viablePeptides->spectrumBins.bin[i];
    if (head.size != 0) {
      link = head.next;
      while (link != NULL) {
        node = copyNode(link->node);
        addToRejectedDatabase(rejectedDatabase,
                              node,
                              i,
                              node.charge);
        link = link->next;
      }
    }
    destroyLinkedList(&viablePeptides->spectrumBins.bin[i]);
  }

  for (int i = 0; i < viablePeptides->spectrumDoubleBins.size; i++) {
    head = viablePeptides->spectrumDoubleBins.bin[i];
    if (head.size != 0) {
      link = head.next;
      while (link != NULL) {
        node = copyNode(link->node);
        addToRejectedDatabase(rejectedDatabase,
                              node,
                              i,
                              node.charge);
        link = link->next;
      }
    }
    destroyLinkedList(&viablePeptides->spectrumDoubleBins.bin[i]);
  }

}

void mergeMiscleavages(ViablePeptides* viablePeptides,
                       RejectedDatabase* rejectedDatabase,
                       SpectrumBin oldMiscleavages,
                       AminoScores aminoScores,
                       SpectrumHash spectrumHash,
                       SpectrumHash spectrumHashDouble,
                       int maxInstanceMiscleavage,
                       int maxTotalMiscleavage,
                       int maxCharge,
                       int BPenalty,
                       int precision,
                       float massTolerance,
                       int averageMassTolerance,
                       int maximumMassTolerance,
                       float lowestPenalty,
                       long peptideMass,
                       long protonMass,
                       long lowestSpecMass)
{
  Node node;
  LinkedList *link;
  Head head;

  for (int i = 0; i < viablePeptides->miscleavageBins.size; i++) {
    head = viablePeptides->miscleavageBins.bin[i];
    if (head.size != 0) {
      link = head.next;
      while (link != NULL) {
        node = rescoreNode(copyNode(link->node),
                           aminoScores,
                           averageMassTolerance,
                           maximumMassTolerance,
                           lowestPenalty);

        addToViablePeptides(node,
                            viablePeptides,
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
                            node.charge,
                            protonMass,
                            lowestSpecMass);

        link = link->next;
      }
    }
  }

  destroySpectrumBin(oldMiscleavages);

}