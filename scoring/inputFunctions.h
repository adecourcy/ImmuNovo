void getNext(char* buffer, int* index, FILE* fp);
SpectrumHash fileToSpectrumHash(FILE* fileName, int precision);
AminoAcids fileToAminoAcids(FILE* file);
AminoMasses fileToAminoMasses(FILE* file, AminoAcids aminoAcids);
AminoScores fileToAminoScores(FILE* file,
                              AminoAcids aminoAcids,
                              int pepLength);