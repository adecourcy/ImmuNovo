typedef struct spectrumBin_tag {
  int size;
  struct head_tag *bin;
} SpectrumBin;

typedef struct viable_tag {
  int binSize;
  struct spectrumBin_tag spectrumBins;
  struct spectrumBin_tag spectrumDoubleBins;
  struct spectrumBin_tag miscleavageBins;
} ViablePeptides;


void printViablePeptides(ViablePeptides viablePeptides,
                         int debug);

SpectrumBin createMiscleavageBins(long largestSpectrum,
                                  int precision);

SpectrumBin createSpectrumBins(int spectrumLength);

ViablePeptides createViablePeptides(long largestSpectrum,
                                    int singleLength,
                                    int doubleLength,
                                    int binSize,
                                    int precision);

void destroyViablePeptides(ViablePeptides viablePeptides);

void destroySpectrumBin(SpectrumBin spectrumBin);

int getMiscleavageBinNumber(Node node, int precision);

void addToMiscleavageBin(ViablePeptides* viablePeptides,
                         Node node,
                         int precision);

long long calculateAdjustedScore(Node node,
                                 int averageMassTolerance,
                                 int maximumMassTolerance,
                                 float lowestPenalty);

void addToSpectrumBinHelper(Head *head,
                            RejectedDatabase *rejectedDatabase,
                            Node node,
                            int charge,
                            int binNumber,
                            int binSize);

int trySpectrumBin(SpectrumBin spectrumBin,
                   RejectedDatabase *rejectedDatabase,
                   SpectrumHash spectrumHash,
                   Node node,
                   long mass,
                   double BPenalty,
                   float massTolerance,
                   int averageMassTolerance,
                   int maximumMassTolerance,
                   float lowestPenalty,
                   int charge,
                   int binSize,
                   short bCompliment,
                   int precision,
                   short doubleCharge);

int addToSpectrumBin(ViablePeptides* viablePeptides,
                     RejectedDatabase *rejectedDatabase,
                     SpectrumHash spectrumHash,
                     SpectrumHash spectrumHashDouble,
                     Node node,
                     double BPenalty,
                     int precision,
                     float massTolerance,
                     int averageMassTolerance,
                     int maximumMassTolerance,
                     float lowestPenalty,
                     long peptideMass,
                     int charge,
                     long protonMass,
                     long lowestSpecMass);

long findBCompliment(long peptideMass, long nodeMass, long protonMass);

void addToViablePeptides(Node node,
                         ViablePeptides* viablePeptides,
                         RejectedDatabase* rejectedDatabase,
                         SpectrumHash spectrumHash,
                         SpectrumHash spectrumHashDouble,
                         int maxInstanceMiscleavage,
                         int maxTotalMiscleavage,
                         int maxCharge,
                         double BPenalty,
                         int precision,
                         float massTolerance,
                         int averageMassTolerance,
                         int maximumMassTolerance,
                         float lowestPenalty,
                         long peptideMass,
                         int charge,
                         long protonMass,
                         long lowestSpecMass);

ViablePeptides getNextGeneration(ViablePeptides* viablePeptides,
                                 RejectedDatabase* rejectedDatabase,
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
                                 long lowestSpecMass);