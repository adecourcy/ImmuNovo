typedef struct spectrumNode_tag {
  long mass;
  long long score;
  int positionNumber;
  struct spectrumNode_tag* next;
} SpectrumNode;

typedef struct spectrumHash_tag {
  int size;
  struct spectrumNode_tag* array;
} SpectrumHash;

SpectrumHash newSpectrumTable(int size);
void addSpectrumNode(long mass,
                     long long score,
                     int positionNumber,
                     SpectrumHash *table,
                     int precision);
int getSpectrumHash(long mass, int size, int precision);
void printSpectrumTable(SpectrumHash* table);
void destroySpectrumNode(SpectrumNode* node);
void destroySpectrumHashTable(SpectrumHash* table);
SpectrumNode lookupMass(long mass,
                        int precision,
                        SpectrumHash table,
                        float massTolerance,
                        float *massToleranceResult);
float massToleranceCalculation(long calculatedMass,
                               long experimentalMass);