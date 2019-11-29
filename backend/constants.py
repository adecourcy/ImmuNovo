from decimal import Decimal

PROTONMASS=Decimal("1.0072766")
H2OMASS=Decimal("18.01528")
NH3MASS=Decimal("17.03052")


# Constants for CSV Headers
SCORE_GLOBAL = 'Global Score'
SCORE_PSSM = 'PSSM Score'
SCORE_COMBINED = 'Combined Score'

PEPTIDE = 'Peptide'
PEPTIDE_GLOBAL_BEST = 'Peptide Global Best'
PEPTIDE_PSSM_BEST = 'Peptide PSSM Best'
PEPTIDE_COMBINED_BEST = 'Peptide Combined Best'

TITLE_PSSM = 'PSSM Title'
TITLE_SPECTRUM = 'Spectrum Title'

FILE_SPECTRUM = 'Spectrum File'


# Constants for CSV entries
NO_PEP = 'No Peptide Found'


# For use in analysis scripts
RESULT_IMMUNO = 'IMMUNOSCORES'
RESULTS_DECOY = 'DECOYSCORES'
RESULTS_DELTA = 'SCOREDELTA'
FDR = 'FDR'

MSGF = 'MSGF'