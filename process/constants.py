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
PEPTIDE = "Peptide"
PEPTOOLS_ANNOTATIONS = "PepTools_annotations"
CONSENSUS_MOTIF_MATCH = "consenus_motif_match"
KINASE_MOTIF_MATCH = "kinase_motif_match"

TOTAL_PROTEIN_INDEX_FILE = "total_protein_index.json"

DEFAULT_FISHER_STATS = {G_STATISTIC: 1, P_VALUE: 1, FREQUENCY: 1, Q_VALUE: 1}

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
ABUNDANCES_NORMALISED_LOG2_ARREST = "abundances normalised log2 arrest"

PHOSPHO_REGRESSION_ZERO_MAX = "phospho_regression 0-max"
PHOSPHO_REGRESSION_LOG2_MEAN = "phospho_regression log2_mean"
PHOSPHO_REGRESSION_POSITION_ABUNDANCES_RAW = (
    "phospho_regression position_abundances raw"
)
PHOSPHO_REGRESSION_POSITION_ABUNDANCES_IMPUTED = (
    "phospho_regression position_abundances imputed"
)
PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_ZERO_MAX = (
    "phospho_regression position_abundances normalised 0-max"
)
PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_MEDIAN = (
    "phospho_regression position_abundances normalised median"
)
PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_MIN_MAX = (
    "phospho_regression position_abundances normalised min-max"
)
PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_LOG2_MEAN = (
    "phospho_regression position_abundances normalised log2_mean"
)
PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_LOG2_ARREST = (
    "phospho_regression position_abundances normalised log2 arrest"
)
PROTEIN_OSCILLATION_ABUNDANCES_ZERO_MAX = "protein_oscillation_abundances 0-max"
PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN = "protein_oscillation_abundances log2_mean"

SL_SAMPLE_STAGE_NAME_1 = "G1"
SL_SAMPLE_STAGE_NAME_2 = "Early S"
SL_SAMPLE_STAGE_NAME_3 = "Mid S"
SL_SAMPLE_STAGE_NAME_4 = "Late S"
SL_SAMPLE_STAGE_NAME_5 = "Mid G2"
SL_SAMPLE_STAGE_NAME_6 = "Late G2"
SL_SAMPLE_STAGE_NAME_7 = "Prophase"
SL_SAMPLE_STAGE_NAME_8 = "Prometaphase arrest"
SL_SAMPLE_STAGE_NAME_9 = "Metaphase"
SL_SAMPLE_STAGE_NAME_10 = "Cytokinesis"

ICR_ABUNDANCE_REP_1 = "abundance_rep_1"
ICR_ABUNDANCE_REP_2 = "abundance_rep_2"

PROTEIN_MAX_Q = 0.05
PHOSPHO_MAX_Q = 0.01

PROJECT_SL = "SL"
PROJECT_ICR = "ICR"
PROJECT_ORIGINAL = "Original"

CCD_BATCH_SIZE = 400
