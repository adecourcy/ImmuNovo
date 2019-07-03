typedef struct results_tag {
  SpectrumBin singleResults;
  SpectrumBin doubleResults;
} Results;

Results createResults(int size);
Results findBestPeptides(ViablePeptides* viablePeptides,
                      AminoMasses aminoMasses,
                      AminoScores* aminoScoresArray,
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
                      long lowestSpecMass,
                      int minPeptideLength,
                      int maxPeptideLength,
                      long H20Mass,
                      long protonMass,
                      short yFirst);