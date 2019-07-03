#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "includes.h"

Node createNode(int maxPeptideLength)
{
  Node node;
  node.peptideString = (char*) GC_MALLOC(maxPeptideLength * sizeof(char));
  //node.peptideString = (char*) malloc(maxPeptideLength * sizeof(char));

  node.currentMiscleavage = 0;
  node.totalMiscleavage = 0;
  node.massTolerance = 0;
  node.mass = 0;
  node.massScore = 0;
  node.aminoScore = 0;
  node.totalScore = 0;
  node.adjustedScore = 0;
  node.pepLength = 0;
  node.charge = 0;
  node.singleCharge = 1;
  node.doubleCharge = 0;
  node.miscleavageSwap = 0;
  node.precursorError = 0;
  node.maxPeptideLength = maxPeptideLength;

  return node;

}

void destroyNode(Node node)
{
  //free(node.peptideString);
}

Node copyNode(Node node)
{
  Node newNode = createNode(node.maxPeptideLength);

  newNode.currentMiscleavage = node.currentMiscleavage;
  newNode.totalMiscleavage = node.totalMiscleavage;
  newNode.massTolerance = node.massTolerance;
  newNode.mass = node.mass;
  newNode.massScore = node.massScore;
  newNode.aminoScore = node.aminoScore;
  newNode.totalScore = node.totalScore;
  newNode.adjustedScore = node.adjustedScore;
  strcpy(newNode.peptideString, node.peptideString);
  newNode.pepLength = node.pepLength;
  newNode.charge = node.charge;
  newNode.singleCharge = node.singleCharge;
  newNode.doubleCharge = node.doubleCharge;
  newNode.miscleavageSwap = node.miscleavageSwap;
  newNode.precursorError = node.precursorError;
  newNode.maxPeptideLength = node.maxPeptideLength;

  return newNode;

}

void addAmino(Node node, char* amino)
{
  strcat(node.peptideString, amino);
}

void printNode(Node node)
{
  printf("currentMiscleavage=%d\n", node.currentMiscleavage);
  printf("totalMiscleavage=%d\n", node.totalMiscleavage);
  printf("massTolerance=%f\n", node.massTolerance);
  printf("mass=%ld\n", node.mass);
  printf("massScore=%lld\n", node.massScore);
  printf("aminoScore=%ld\n", node.aminoScore);
  printf("totalScore=%lld\n", node.totalScore);
  printf("adjustedScore=%lld\n", node.adjustedScore);
  printf("peptideString=%s\n", node.peptideString);
  printf("charge=%d\n", node.charge);
  if (node.singleCharge == 0) {
    printf("singleCharge=%s\n", "False");
  }
  else {
    printf("singleCharge=%s\n", "True");
  }
  if (node.doubleCharge == 0) {
    printf("doubleCharge=%s\n", "False");
  }
  else {
    printf("doubleCharge=%s\n", "True");
  }
  printf("\n");

}



Node* branchPeptide(Node node,
                    AminoScores aminoScores,
                    AminoMasses aminoMasses,
                    AminoAcids aminoAcids)
{
  Node* nodeList = (Node*) GC_MALLOC(aminoAcids.size * sizeof(Node));
  //Node* nodeList = (Node*) malloc(aminoAcids.size * sizeof(Node));

  for (int i = 0; i < aminoAcids.size; i++) {
    Node newNode = copyNode(node);
    newNode.pepLength += 1;

    newNode.peptideString[newNode.pepLength - 1] = aminoAcids.array[i];
    newNode.peptideString[newNode.pepLength] = '\0';
    newNode.mass += aminoMasses.aminoAcid[i];
    long long aminoScore = getAminoScore(aminoScores,
                                         aminoAcids.array[i],
                                         newNode.pepLength - 1);
    newNode.aminoScore += aminoScore;
    newNode.totalScore += aminoScore;
    if (aminoAcids.array[i] == 'K') {
      newNode.doubleCharge = 1;
    }

    if (newNode.miscleavageSwap == 1) {
      char tmp1, tmp2;
      tmp1 = newNode.peptideString[newNode.pepLength - 1];
      tmp2 = newNode.peptideString[newNode.pepLength - 2];
      newNode.peptideString[newNode.pepLength - 1] = tmp2;
      newNode.peptideString[newNode.pepLength - 2] = tmp1;
      newNode.miscleavageSwap = 0;
    }

    nodeList[i] = newNode;
  }

  return nodeList;

}