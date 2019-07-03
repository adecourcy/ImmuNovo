typedef struct aminoMasses_tag {
  int size;
  long *aminoAcid;
} AminoMasses;

typedef struct aminoAcids_tag {
  int size;
  char *array;
} AminoAcids;

typedef struct positionScores_tag {
  int size;
  long *position;
} PositionScores;

typedef struct aminoScores_tag {
  int size;
  struct positionScores_tag *aminoAcid;
} AminoScores;

typedef struct lengthScore_tag {
  int size;
  struct aminoScores_tag *aminoLength;
} LengthScore;

AminoMasses createAminoMasses(int size);
AminoAcids createAminoAcids(int size);
PositionScores createPositionScores(int size);
AminoScores createAminoScores(int numAcids, int pepLength);
LengthScore createLengthScores(int size);
void destroyAminoMasses(AminoMasses aminoMasses);
void destroyAminoAcids(AminoAcids AminoAcids);
void destroyPositionScores(PositionScores positionScores);
void destroyAminoScores(AminoScores aminoScores);
void destroyLengthScores(LengthScore lengthScores);
int getHash(char aminoAcid);
long getAminoMass(AminoMasses aminoMasses, char aminoAcid);
AminoMasses addAminoMass(AminoMasses aminoMasses, char aminoAcid, long mass);
long getAminoScore(AminoScores aminoScores, char aminoAcid, int position);
AminoScores addAminoScore(AminoScores aminoScores,
                          char aminoAcid,
                          int position,
                          long score);
