typedef struct rejectLinkedList_tag {
  struct rejectLinkedList_tag *next;
  struct node_tag node;
  int binNumber;
  int charge;
} RejectedLinkedList;

typedef struct rejected_tag
{
  int maxPeptideLength;
  struct rejectLinkedList_tag **array;
  
} RejectedDatabase;

RejectedDatabase createRejectedDatabase(int maxPeptideLength);
void printRejectedNode(Node node);
void printRejectedList(RejectedLinkedList *list);
void printRejectedDatabase(RejectedDatabase database);
void destroyRejectedList(RejectedLinkedList *link);
void destroyRejectedDatabase(RejectedDatabase database);
void addToRejectedDatabase(RejectedDatabase* database,
                           Node node,
                           int binNumber,
                           int charge);