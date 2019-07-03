#include <stdio.h>
#include <stdlib.h>
#include "includes.h"
#include "inputFunctions.h"

// This function gets the next whitespace-separated token
void getNext(char* buffer, int* index, FILE* fp)
{
  int i;
  for (i = 0; i < 255; i++) {
    buffer[i] = fgetc(fp);
    if (buffer[i] == '\n' ||
        buffer[i] == ' ') {
      buffer[i] = '\0';
      break;
    } else if (buffer[i] == EOF) {
      break;
    }
  }
  
  if (i == 254) {
    printf("Output file name is too long. Must be under 255 characters\n");
    exit(0);
  }
  
  *index = i;
}

SpectrumHash fileToSpectrumHash(FILE* fileName, int precision)
{
  char buffer[255];
  int i;
  long mass;
  long long score;
  getNext(buffer, &i, fileName);
  SpectrumHash newSpectrumHash = newSpectrumTable((int) strtol(buffer, NULL, 10));
  for (int specPosition = 0; specPosition < newSpectrumHash.size; specPosition++) {
    getNext(buffer, &i, fileName);
    mass = strtol(buffer, NULL, 10);
    getNext(buffer, &i, fileName);
    score = strtoll(buffer, NULL, 10);
    addSpectrumNode(mass, score, specPosition, &newSpectrumHash, precision);
  }

  return newSpectrumHash;

}

AminoAcids fileToAminoAcids(FILE* file)
{
  char buffer[255];
  int i;
  int numAcids;
  getNext(buffer, &i, file);
  numAcids = (int) strtol(buffer, NULL, 10);
  AminoAcids aminoAcids = createAminoAcids(numAcids);
  for (int k = 0; k < numAcids; k++) {
    getNext(buffer, &i, file);
    aminoAcids.array[k] = buffer[0];
  }

  return aminoAcids;
}

AminoMasses fileToAminoMasses(FILE* file, AminoAcids aminoAcids)
{
  char buffer[255];
  int k;
  AminoMasses aminoMasses = createAminoMasses(aminoAcids.size);
  for (int i = 0; i < aminoAcids.size; i++) {
    getNext(buffer, &k, file);
    aminoMasses = addAminoMass(aminoMasses,
                               aminoAcids.array[i],
                               strtol(buffer, NULL, 10));
  }

  return aminoMasses;

}

AminoScores fileToAminoScores(FILE* file,
                              AminoAcids aminoAcids,
                              int pepLength)
{
  char buffer[255];
  int a;
  AminoScores aminoScores = createAminoScores(aminoAcids.size, pepLength);
  for (int i = 0; i < aminoAcids.size; i++) {
    for (int k = 0; k < pepLength; k++) {
      getNext(buffer, &a, file);
      aminoScores = addAminoScore(aminoScores,
                                  aminoAcids.array[i],
                                  k,
                                  strtol(buffer, NULL, 10));
    }
  }

  return aminoScores;
}