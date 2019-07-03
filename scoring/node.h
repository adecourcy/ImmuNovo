typedef struct node_tag {
  int currentMiscleavage;
  int totalMiscleavage;
  float massTolerance;
  long mass;
  long long massScore;
  long aminoScore;
  long long totalScore;
  long long adjustedScore;
  char* peptideString;
  int pepLength;
  int charge;
  short singleCharge;
  short doubleCharge;
  short miscleavageSwap;
  int precursorError;
  int maxPeptideLength;
} Node;

Node createNode(int maxPeptideLength);
void destroyNode(Node node);
Node copyNode(Node node);
void addAmino(Node node, char* amino);
void printNode(Node node);
Node* branchPeptide(Node node,
                    AminoScores aminoScores,
                    AminoMasses aminoMasses,
                    AminoAcids aminoAcids);