#include <stdlib.h>
#include "includes.h"

void addLink(Head *head, Node node)
{

  LinkedList *link = (LinkedList*) GC_MALLOC(sizeof(LinkedList));
  //LinkedList *link = (LinkedList*) malloc(sizeof(LinkedList));
  link->next = NULL;
  link->node = node;

  LinkedList* current = head->next;
  LinkedList* prev = NULL;

  while (1) {

    if (current == NULL) {
      head->size++;
      if (prev == NULL) {
        head->next = link;
      } else {
        prev->next = link;
      }
      break;

    } else if (node.mass > current->node.mass) {
      head->size++;
      link->next = current;
      if (prev == NULL) {
        head->next = link;
      } else {
        prev->next = link;
      }
      break;

    } else if (node.mass < current->node.mass) {
      prev = current;
      current = current->next;
      continue;

    } else {
      if (node.adjustedScore >= current->node.adjustedScore) {
        link->next = current->next;
        if (prev == NULL) {
          head->next = link;
        } else {
          prev->next = link;
        }
        destroyLink(current);
        break;
      } else {
        destroyLink(link);
        break;
      }

    }
  }
}

Node removeLowest(Head *head)
{

  LinkedList* current = NULL;
  LinkedList* lowestPrev = NULL;
  LinkedList* prev = NULL;
  LinkedList* lowest = current;
  long long lowestScore;

  head->size--;
  current = head->next;
  lowestScore = current->node.adjustedScore;
  lowest = current;

  while (current != NULL) {
    if (current->node.adjustedScore < lowestScore) {
      lowest = current;
      lowestPrev = prev;
      lowestScore = current->node.adjustedScore;
    }
    prev = current;
    current = current->next;

  }

  if (lowestPrev == NULL) {
      head->next = head->next->next;
  } else if (lowest->next == NULL) {
      lowestPrev->next = NULL;
  } else {
      lowestPrev->next = lowest->next;
  }
  
  Node lowestNode = lowest->node;
  //free(lowest);
  return lowestNode;

}

void printList(LinkedList *link)
{
  while (link != NULL) {

    printNode(link->node);
    link = link->next;

  }
}

void destroyLink(LinkedList *link)
{ /*
  destroyNode(link->node);
  free(link);
  */
}

void destroyLinkedListHelper(LinkedList *link)
{/*
  if (link != NULL) {
    destroyLinkedListHelper(link->next);
    destroyLink(link);
  }
  */
}

void destroyLinkedList(Head *head)
{
  destroyLinkedListHelper(head->next);
  head->next = NULL;
  head->size = 0;
}

Head newHead()
{
  Head head;
  head.next = NULL;
  head.size = 0;
  return head;
}