#include <stdlib.h>
#include "includes.h"

/*
  The variable massTolerance appears many times throughout this program.
  It is a mistake that needs to be removed at some point. It's essentially
  a floating point version of maximumMassTolerance, and denotes the maximum
  acceptable mass tolerance. I don't remember why I put it in, but it's
  propogated through the entire program and will take some time to remove.
*/

Results createResults(int size)
{
  Results results;
  results.singleResults = createSpectrumBins(size);
  results.doubleResults = createSpectrumBins(size);

  return results;
}

Results findBestPeptides(ViablePeptides* viablePeptides,
                         AminoMasses aminoMasses,
                         AminoScores* aminoScoresArray,
                         AminoAcids aminoAcids,
                         SpectrumHash spectrumHash,
                         SpectrumHash spectrumHashDouble,
                         int maxInstanceMiscleavage,
                         int maxTotalMiscleavage,
                         int maxCharge,
                         int precision,
                         double BPenalty,
                         float massTolerance, // see comment at top of file
                         int averageMassTolerance,
                         int maximumMassTolerance,
                         float lowestPenalty,
                         long peptideMass,
                         long lowestSpecMass,
                         int minPeptideLength,
                         int maxPeptideLength,
                         long H20Mass,
                         long protonMass,
                         short yFirst)
{
  SpectrumBin oldMiscleavages;
  Results results = createResults((maxPeptideLength - minPeptideLength) + 1);
  int resultsIndex = 0;
  LinkedList *link;
  RejectedDatabase rejectedDatabase = createRejectedDatabase(maxPeptideLength);

  int aminoScoresIndex = 0;
  AminoScores aminoScores = aminoScoresArray[aminoScoresIndex];

  for (int currentPeptideLength = 0;
       currentPeptideLength < maxPeptideLength;
       currentPeptideLength++) {

    /*
    #include <stdio.h>
    printf("%d\n", currentPeptideLength);
    printf("########## ViablePeptides %d ##########\n", currentPeptideLength);
    printViablePeptides(*viablePeptides, 1);
    printf("##########################################################\n");
    */
    

    if (currentPeptideLength == 0) {

      Node seedNode = createNode(maxPeptideLength);
      if (yFirst == 1) {
        seedNode.mass = H20Mass;
      }
      addLink(&viablePeptides->miscleavageBins.bin[0], seedNode);


    } else if (currentPeptideLength >= minPeptideLength) {

      aminoScoresIndex++;
      aminoScores = aminoScoresArray[aminoScoresIndex];

      link =
        viablePeptides->spectrumBins.bin[viablePeptides->spectrumBins.size - 1].next;

      while (link != NULL) {
        addLink(&results.singleResults.bin[resultsIndex], copyNode(link->node));
        link = link->next;
      }

      if (viablePeptides->spectrumDoubleBins.size > 0) {
        link = 
          viablePeptides->spectrumDoubleBins.bin[viablePeptides->spectrumDoubleBins.size - 1].next;

        while (link != NULL) {
          addLink(&results.doubleResults.bin[resultsIndex], copyNode(link->node));
          link = link->next;
        }
      }

      resultsIndex++;

      addAllToRejectedDatabase(viablePeptides, &rejectedDatabase);
      oldMiscleavages = viablePeptides->miscleavageBins;
      viablePeptides->miscleavageBins =
            createMiscleavageBins(peptideMass, precision);
      iterateRejectedPeptides(&rejectedDatabase,
                              viablePeptides,
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
                              lowestSpecMass,
                              currentPeptideLength - 1);

      mergeMiscleavages(viablePeptides,
                        &rejectedDatabase,
                        oldMiscleavages,
                        aminoScores,
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
                        protonMass,
                        lowestSpecMass);

    }

    *viablePeptides = getNextGeneration(viablePeptides,
                                        &rejectedDatabase,
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

  }

  link =
    viablePeptides->spectrumBins.bin[viablePeptides->spectrumBins.size - 1].next;

  while (link != NULL) {
    addLink(&results.singleResults.bin[resultsIndex], copyNode(link->node));
    link = link->next;
  }

  if (viablePeptides->spectrumDoubleBins.size > 0) {
    link = 
      viablePeptides->spectrumDoubleBins.bin[viablePeptides->spectrumDoubleBins.size - 1].next;

    while (link != NULL) {
      addLink(&results.doubleResults.bin[resultsIndex], copyNode(link->node));
      link = link->next;
    }
  }

  //exit(0);
  return results;

}