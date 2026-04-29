# Ported from Python CoNGA: tcr_distances.py (lines 11-17)
# and amino_acids.py

#' Gap penalty for V-region alignments
#'
#' Integer gap penalty applied when aligning V-region sequences.
#'
#' @format A length-1 integer.
#' @keywords internal
GAP_PENALTY_V_REGION <- 4L

#' Gap penalty for CDR3-region alignments
#'
#' Integer gap penalty applied when aligning CDR3-region sequences.
#'
#' @format A length-1 integer.
#' @keywords internal
GAP_PENALTY_CDR3_REGION <- 12L

#' Weight for V-region distance contribution
#'
#' Integer weight applied to V-region distances when computing combined TCR
#' distances.
#'
#' @format A length-1 integer.
#' @keywords internal
WEIGHT_V_REGION <- 1L

#' Weight for CDR3-region distance contribution
#'
#' Integer weight applied to CDR3-region distances when computing combined TCR
#' distances.
#'
#' @format A length-1 integer.
#' @keywords internal
WEIGHT_CDR3_REGION <- 3L

#' Gap character used in sequence alignments
#'
#' Single character representing a gap in aligned sequences.
#'
#' @format A length-1 character string.
#' @keywords internal
GAP_CHARACTER <- "."

#' Whether to trim CDR3 sequences before distance computation
#'
#' Logical flag controlling CDR3 trimming in TCR distance calculations.
#'
#' @format A length-1 logical.
#' @keywords internal
TRIM_CDR3S <- TRUE

#' Whether to align CDR3 sequences before distance computation
#'
#' Logical flag controlling CDR3 alignment in TCR distance calculations.
#'
#' @format A length-1 logical.
#' @keywords internal
ALIGN_CDR3S <- FALSE

#' Standard 20 amino acid single-letter codes
#'
#' A character vector of the 20 standard amino acids in alphabetical order,
#' ported from Python CoNGA \code{amino_acids.py}.
#'
#' @format A character vector of length 20.
#' @export
AMINO_ACIDS <- c(
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
)
