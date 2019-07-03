typedef struct linkedList_tag {
  struct linkedList_tag *next;
  struct node_tag node;
} LinkedList;

typedef struct head_tag {
  int size;
  struct linkedList_tag *next;
} Head;

void addLink(Head *head, Node node);
Node removeLowest(Head *head);
void printList(LinkedList *list);
void destroyLinkedListHelper(LinkedList *link);
void destroyLinkedList(Head *list);
void destroyLink(LinkedList *list);
Head newHead();