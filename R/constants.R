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
#'
#' @examples
#' AMINO_ACIDS
#' length(AMINO_ACIDS)  # 20
#'
#' @seealso \code{\link{bsd4_matrix}}, \code{\link{weighted_cdr3_distance}}
#' @export
AMINO_ACIDS <- c(
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
)


# ---------------------------------------------------------------------------
# BLOSUM62 substitution matrix (raw scores, NOT the BSD4 distance matrix)
# ---------------------------------------------------------------------------

#' Build the BLOSUM62 substitution matrix
#'
#' Constructs a 20x20 named integer matrix of BLOSUM62 scores from the
#' standard amino acid substitution table. Used for CDR3 alignment scoring
#' in sequence logo construction.
#'
#' @return A 20x20 named integer matrix.
#' @keywords internal
#' @noRd
.build_blosum62_matrix <- function() {
    aa <- AMINO_ACIDS
    # Upper-triangle values in AMINO_ACIDS order (A,C,D,E,F,G,H,I,K,L,
    # M,N,P,Q,R,S,T,V,W,Y). Diagonal + upper triangle = 210 entries.
    # We fill symmetrically.
    vals <- c(
    #   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
        4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2,  # A
            9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2,  # C
                6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3,  # D
                    5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2,  # E
                        6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3,  # F
                            6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3,  # G
                                8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2,  # H
                                    4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1,  # I
                                        5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2,  # K
                                            4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1,  # L
                                                5, -2, -2,  0, -1, -1, -1,  1, -1, -1,  # M
                                                    6, -2,  0,  0,  1,  0, -3, -4, -2,  # N
                                                        7, -1, -2, -1, -1, -2, -4, -3,  # P
                                                            5,  1,  0, -1, -2, -2, -1,  # Q
                                                                5, -1, -1, -3, -3, -2,  # R
                                                                    4,  1, -2, -3, -2,  # S
                                                                        5,  0, -2, -2,  # T
                                                                            4, -3, -1,  # V
                                                                               11,  2,  # W
                                                                                    7   # Y
    )
    mat <- matrix(0L, nrow = 20L, ncol = 20L, dimnames = list(aa, aa))
    k <- 1L
    for (i in seq_len(20L)) {
        for (j in i:20L) {
            mat[i, j] <- vals[k]
            mat[j, i] <- vals[k]
            k <- k + 1L
        }
    }
    mat
}

#' BLOSUM62 amino acid substitution matrix
#'
#' A 20x20 named integer matrix of BLOSUM62 substitution scores. Row and
#' column names are the 20 standard amino acids in \code{AMINO_ACIDS} order.
#' This is the raw substitution score matrix (positive = similar, negative =
#' dissimilar), NOT the BSD4 distance matrix used in TCRdist computation.
#'
#' @format A 20x20 named integer matrix.
#' @keywords internal
BLOSUM62 <- .build_blosum62_matrix()
