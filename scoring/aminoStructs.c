#include <stdlib.h>
#include "includes.h"

AminoMasses createAminoMasses(int size)
{
  AminoMasses aminos;
  aminos.size = size;
  long *array = (long*) GC_MALLOC(size * sizeof(long));
  //long *array = (long*) malloc(size * sizeof(long));
  aminos.aminoAcid = array;

  return aminos;
}

AminoAcids createAminoAcids(int size)
{
  AminoAcids aminos;
  aminos.size = size;
  char *array = (char*) GC_MALLOC(size * sizeof(char));
  //char *array = (char*) malloc(size * sizeof(char));
  aminos.array = array;

  return aminos;
}

PositionScores createPositionScores(int size)
{
  PositionScores positionScores;
  positionScores.size = size;
  long *array = (long*) GC_MALLOC(size * sizeof(long));
  //long *array = (long*) malloc(size * sizeof(long));
  positionScores.position = array;

  return positionScores;
}

AminoScores createAminoScores(int numAcids, int pepLength)
{
  AminoScores aminos;
  aminos.size = numAcids;
  PositionScores *array =
        (PositionScores*) GC_MALLOC(numAcids * sizeof(PositionScores));
        //(PositionScores*) malloc(numAcids * sizeof(PositionScores));
  aminos.aminoAcid = array;
  for (int i = 0; i < numAcids; i++) {
    aminos.aminoAcid[i] = createPositionScores(pepLength);
  }

  return aminos;
}

LengthScore createLengthScores(int size)
{
  LengthScore scores;
  scores.size = size;
  AminoScores *array = (AminoScores*) GC_MALLOC(size* sizeof(AminoScores));
  //AminoScores *array = (AminoScores*) malloc(size* sizeof(AminoScores));
  scores.aminoLength = array;

  return scores;
}

void destroyAminoMasses(AminoMasses aminoMasses)
{
  //free(aminoMasses.aminoAcid);
}

void destroyAminoAcids(AminoAcids aminoAcids)
{
  //free(aminoAcids.array);
}

void destroyPositionScores(PositionScores positionScores)
{
  //free(positionScores.position);
}

void destroyAminoScores(AminoScores aminoScores)
{ /*
  for (int i = 0; i < aminoScores.size; i++) {
    destroyPositionScores(aminoScores.aminoAcid[i]);
  }
  free(aminoScores.aminoAcid);
  */
}

void destroyLengthScores(LengthScore lengthScores)
{/*
  for (int i = 0; i < lengthScores.size; i++) {
    destroyAminoScores(lengthScores.aminoLength[i]);
  }
  free(lengthScores.aminoLength);
  */
}

int getHash(char aminoAcid)
{
  switch(aminoAcid) {
    case 'G':
      return 0;
    case 'A':
      return 1;
    case 'S':
      return 2;
    case 'P':
      return 3;
    case 'V':
      return 4;
    case 'T':
      return 5;
    case 'L':
      return 6;
    case 'I':
      return 7;
    case 'N':
      return 8;
    case 'D':
      return 9;
    case 'Q':
      return 10;
    case 'K':
      return 11;
    case 'E':
      return 12;
    case 'M':
      return 13;
    case 'H':
      return 14;
    case 'F':
      return 15;
    case 'R':
      return 16;
    case 'C':
      return 17;
    case 'Y':
      return 18;
    case 'W':
      return 19;
    case 'B':
      return 20;
    case 'J':
      return 21;
    case 'O':
      return 22;
  }
}

long getAminoMass(AminoMasses aminoMasses, char aminoAcid)
{
  return aminoMasses.aminoAcid[getHash(aminoAcid)];
}

AminoMasses addAminoMass(AminoMasses aminoMasses, char aminoAcid, long mass)
{
  aminoMasses.aminoAcid[getHash(aminoAcid)] = mass;
  return aminoMasses;
}

long getAminoScore(AminoScores aminoScores, char aminoAcid, int position)
{
  return aminoScores.aminoAcid[getHash(aminoAcid)].position[position];
}

AminoScores addAminoScore(AminoScores aminoScores,
                          char aminoAcid,
                          int position,
                          long score)
{
  aminoScores.aminoAcid[getHash(aminoAcid)].position[position] = score;
  return aminoScores;
}