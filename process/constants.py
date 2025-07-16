RAW = "raw"
METRICS = "metrics"
LOG2_MEAN = "log2_mean"
ZERO_MAX = "0-max"
ANOVA = "ANOVA"
Q_VALUE = "q_value"
FISHER_G = "Fisher_G"
PROTEIN_ABUNDANCES = "protein_abundances"
ABUNDANCE_AVERAGE = "average"
NORMALISED = "normalised"
MEDIAN = "median"
MIN_MAX = "min-max"
LOG2_ARREST = "log2 arrest"
IMPUTED = "imputed"
P_VALUE = "p_value"
# TODO - should it be F_statistic, singular?
F_STATISTICS = "F_statistics"

# TODO - change all these to not have PROTEIN_
PROTEIN_OSCILLATION_ABUNDANCES = "protein_oscillation_abundances"
ABUNDANCE_AVERAGE = "abundance_average"
PROTEIN_INFO = "protein_info"
GENE_NAME = "gene_name"
PROTEIN_NAME = "protein_name"

POSITION_ABUNDANCES = "position_abundances"
PHOSPHORYLATION_ABUNDANCES = "phosphorylation_abundances"
PHOSPHO_REGRESSION = "phospho_regression"
PHOSPHORYLATION_SITE = "phosphorylation_site"
CURVE_FOLD_CHANGE = "curve_fold_change"
PROTEIN_PHOSPHO_CORRELATION = "protein-phospho-correlation"
# TODO - misspelled, but it's difficult to replace it in the fixture json files
PHOSPHO_PROTEIN_CFC_RATIO = "phosho-protein-cfc_ratio"
G_STATISTIC = "G_statistic"
FREQUENCY = "frequency"
KINASE_PREDICTION = "kinase_prediction"
PEPTIDE_SEQ = "peptide_seq"
# TODO - is this a good name? A bit vague.
PHOSPHO = "phospho"
# These are capitalised as they're fields found in the index_protein_names json file
PEPTIDE = 'Peptide'
PEPTOOLS_ANNOTATIONS = "PepTools_annotations"
CONSENSUS_MOTIF_MATCH = "consenus_motif_match"
KINASE_MOTIF_MATCH = "kinase_motif_match"

TOTAL_PROTEIN_INDEX_FILE = "total_protein_index.json"

DEFAULT_FISHER_STATS = {G_STATISTIC: 1, P_VALUE: 1, FREQUENCY: 1, Q_VALUE: 1}

PROTEIN_INFO_FIELDS = [
    "halflife_mean",
    "halflife_std",     
    "halflife_min",     
    "halflife_count",
    "relative_abundance_8h_count",
    "relative_abundance_8h_mean",   
    "relative_abundance_8h_std",    
    "mean_gene_effect",   
    "in_DRIVE_cancer_proteins",
    "in_CGC_cancer_proteins",
    "role_in_cancer",
    "tier",
    "cell_cycle",
    "mitotic_cell_cycle",
    "kinetochore",
    "spindle",
    "centriole",
    "replication fork",
    "G0_to_G1_transition",
    "G1/S_transition",
    "G2/M_transition",
    "S_phase",
    "transcription_factor",
    "kinase_domain_containing",
    "is_E3",
    "APC_complex",
    "dna_replication_machinery"
]

# StatisticTypes

PROTEIN_MEDIAN = "protein median"
PHOSPHO_MEDIAN = "phospho median"

# 0-max and log2_mean already defined above
# TODO - can probably shorten these names

ABUNDANCES_RAW = "abundances raw"
ABUNDANCES_IMPUTED = "abundances imputed"

ABUNDANCES_NORMALISED_LOG2_MEAN = "abundances normalised log2_mean"
ABUNDANCES_NORMALISED_MIN_MAX = "abundances normalised min-max"
ABUNDANCES_NORMALISED_MEDIAN = "abundances normalised median"
ABUNDANCES_NORMALISED_ZERO_MAX = "abundances normalised 0-max"
ABUNDANCES_NORMALISED_LOG2_ARREST ="abundances normalised log2 arrest"

PHOSPHO_REGRESSION_ZERO_MAX = "phospho_regression 0-max"
PHOSPHO_REGRESSION_LOG2_MEAN = "phospho_regression log2_mean"
PHOSPHO_REGRESSION_POSITION_ABUNDANCES_RAW = "phospho_regression position_abundances raw"
PHOSPHO_REGRESSION_POSITION_ABUNDANCES_IMPUTED = "phospho_regression position_abundances imputed"
PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_ZERO_MAX = "phospho_regression position_abundances normalised 0-max"
PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_MEDIAN = "phospho_regression position_abundances normalised median"
PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_MIN_MAX = "phospho_regression position_abundances normalised min-max"
PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_LOG2_MEAN = "phospho_regression position_abundances normalised log2_mean"
PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_LOG2_ARREST = "phospho_regression position_abundances normalised log2 arrest"
PROTEIN_OSCILLATION_ABUNDANCES_ZERO_MAX = "protein_oscillation_abundances 0-max"
PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN = "protein_oscillation_abundances log2_mean"

SL_SAMPLE_STAGE_NAME_1 = "RO sample 15 min"
SL_SAMPLE_STAGE_NAME_2 = "RO sample 30 min"
SL_SAMPLE_STAGE_NAME_3 = "Palb sample, 1h"
SL_SAMPLE_STAGE_NAME_4 = "RO sample 2h"
SL_SAMPLE_STAGE_NAME_5 = "Palb sample, 2.5h"
SL_SAMPLE_STAGE_NAME_6 = "Palb sample, 4h"
SL_SAMPLE_STAGE_NAME_7 = "RO sample 5h"
SL_SAMPLE_STAGE_NAME_8 = "Palb sample, 8h"
SL_SAMPLE_STAGE_NAME_9 = "Palb sample, 10h"
SL_SAMPLE_STAGE_NAME_10 = "Nocodozole"

ICR_ABUNDANCE_REP_1 = "abundance_rep_1"
ICR_ABUNDANCE_REP_2 = "abundance_rep_2"

PROTEIN_MAX_Q = 0.5
PHOSPHO_MAX_Q = 0.1