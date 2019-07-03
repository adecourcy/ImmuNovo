#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "includes.h"

SpectrumHash newSpectrumTable(int size)
{
  SpectrumHash table;
  table.array = (SpectrumNode*) GC_MALLOC(size * sizeof(SpectrumNode));
  //table.array = (SpectrumNode*) malloc(size * sizeof(SpectrumNode));
  table.size = size;

  SpectrumNode node;
  node.mass = 0;

  for (int i = 0; i < size; i++) {
    table.array[i] = node;
  }

  return table;
}

void addSpectrumNode(long mass,
                     long long score,
                     int positionNumber,
                     SpectrumHash *table,
                     int precision)
{
  SpectrumNode* node = (SpectrumNode*) GC_MALLOC(sizeof(SpectrumNode));
  //SpectrumNode* node = (SpectrumNode*) malloc(sizeof(SpectrumNode));
  node->mass = mass;
  node->score = score;
  node->positionNumber = positionNumber;
  node->next = NULL;

  int nodeHash = getSpectrumHash(mass, table->size, precision);
  if (table->array[nodeHash].mass == 0) {
    table->array[nodeHash] = *node;
  } else {
    SpectrumNode *headNode = &table->array[nodeHash];
    while (1) {
      if (headNode->next == NULL) {
        headNode->next = node;
        break;
      } else {
        headNode = headNode->next;
      }

    }
  }
}

int getSpectrumHash(long mass, int size, int precision)
{
  mass = (mass / pow(10, precision));
  return (mass % size);
}

void printSpectrumTable(SpectrumHash* table)
{
  for (int i = 0; i < table->size; i++) {
    if (table->array[i].mass == 0) {
      printf("%d:\n\n", i);
    } else {
      printf("%d: ", i);
      SpectrumNode node = table->array[i];
      while (1) {
        printf("[%ld, %lld, %d], ", node.mass, node.score, node.positionNumber);
        if (node.next == NULL) {
          printf("\n");
          break;
        }
        node = *node.next;
      }
    }
  }
}

void destroySpectrumNode(SpectrumNode* node)
{/*
  if (node->next != NULL) {
    destroySpectrumNode(node->next);
    free(node->next);
  }
  */
}

void destroySpectrumHashTable(SpectrumHash* table)
{/*
  for (int i = 0; i < table->size; i++) {
    destroySpectrumNode(&table->array[i]);
  }
  free(table->array);
  */
}

float massToleranceCalculation(long calculatedMass,
                               long experimentalMass)
{
  return ((((double) (calculatedMass - experimentalMass)) *
              1000000.0) /
              ((double) experimentalMass));
}

/*
  Try to find a spectrum mass given a calculated mass. Return a node if
  a spectrum mass is found. Otherwise return a dummy node with a mass of 0
*/
SpectrumNode lookupMass(long mass,
                        int precision,
                        SpectrumHash table,
                        float massTolerance,
                        float *massToleranceResult)
{
  int hashNumber = getSpectrumHash(mass, table.size, precision);
  SpectrumNode *node = &table.array[hashNumber];

  if (node->mass == 0) {
    return *node;
  } else {
    while (1) {
      *massToleranceResult = massToleranceCalculation(mass, node->mass);
      if (*massToleranceResult < 0) {
        *massToleranceResult *= -1.0;
      }
      if (*massToleranceResult <= massTolerance) {
        return *node;
      } else if(node->next == NULL) {
        SpectrumNode dummyNode;
        dummyNode.mass = 0;
        return dummyNode;
      } else {
        node = node->next;
      }
    }
  }
}