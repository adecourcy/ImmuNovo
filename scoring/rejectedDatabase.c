#include <stdio.h>
#include <stdlib.h>
#include "includes.h"

RejectedDatabase createRejectedDatabase(int maxPeptideLength)
{
  RejectedDatabase database;
  database.maxPeptideLength = maxPeptideLength;
  database.array =
    (RejectedLinkedList**) GC_MALLOC(maxPeptideLength * sizeof(RejectedLinkedList*));
    //(RejectedLinkedList**) malloc(maxPeptideLength * sizeof(RejectedLinkedList*));

  for (int i = 0; i < maxPeptideLength; i++) {
    database.array[i] = NULL;
  }

  return database;
}

void printRejectedNode(Node node)
{
  printf("%s, ", node.peptideString);
}

void printRejectedList(RejectedLinkedList *list)
{
  RejectedLinkedList *next =
      (RejectedLinkedList*) GC_MALLOC(sizeof(RejectedLinkedList));
      //(RejectedLinkedList*) malloc(sizeof(RejectedLinkedList));
  next = list;

  while (next != NULL) {
    printRejectedNode(next->node);
    next = next->next;
  }
}

void printRejectedDatabase(RejectedDatabase database)
{
  for (int i = 0; i < database.maxPeptideLength; i++) {
    printf("Bin Number %d:\n\n", (i+1));
    printRejectedList(database.array[i]);
    printf("\n\n");
  }
}

void destroyRejectedList(RejectedLinkedList *link)
{/*
  if (link != NULL) {
    destroyRejectedList(link->next);
    destroyNode(link->node);
    free(link);
  }
  */
}

void destroyRejectedDatabase(RejectedDatabase database)
{/*
  for (int i = 0; i < database.maxPeptideLength; i++) {
    destroyRejectedList(database.array[i]);
    database.array[i] = NULL;
  }
  free(database.array);
  */
}

void addToRejectedDatabase(RejectedDatabase* database,
                           Node node,
                           int binNumber,
                           int charge)
{
  if (node.pepLength <= database->maxPeptideLength) {
  
    RejectedLinkedList *link =
          (RejectedLinkedList*) GC_MALLOC(sizeof(RejectedLinkedList));
          //(RejectedLinkedList*) malloc(sizeof(RejectedLinkedList));
    link->next = database->array[node.pepLength - 1];
    link->node = node;
    link->binNumber = binNumber;
    link->charge = charge;

    database->array[node.pepLength - 1] = link;

    } 
}