# Ported from Python CoNGA: conga/tcrdist/translation.py
# and conga/tcrdist/logo_tools.py (NUCLEOTIDE_CLASSES, BASE_PARTNER,
# reverse_complement)

#' IUPAC degenerate nucleotide classes (lowercase)
#'
#' Named list mapping each IUPAC nucleotide symbol (lowercase) to a character
#' string of its constituent standard bases (a, c, g, t). Matches the Python
#' \code{nucleotide_classes_lower_case} dictionary used by CoNGA.
#'
#' @format A named list of length 15.
#' @keywords internal
NUCLEOTIDE_CLASSES <- list(
    a = "a",
    c = "c",
    g = "g",
    t = "t",
    w = "at",
    s = "cg",
    k = "gt",
    m = "ac",
    y = "ct",
    r = "ag",
    b = "cgt",
    d = "agt",
    h = "act",
    v = "acg",
    n = "acgt"
)

#' Nucleotide complement mapping
#'
#' Named character vector mapping each nucleotide (and IUPAC ambiguity code and
#' gap character) to its Watson-Crick complement. Matches the Python
#' \code{base_partner} dictionary used by CoNGA.
#'
#' @format A named character vector of length 17.
#' @keywords internal
BASE_PARTNER <- c(
    a = "t", c = "g", g = "c", t = "a", n = "n",
    A = "T", C = "G", G = "C", T = "A", N = "N",
    R = "Y", Y = "R", S = "S", W = "W",
    K = "M", M = "K", "." = "."
)

#' Reverse complement of a nucleotide sequence
#'
#' Reverses a nucleotide sequence and replaces each base with its complement
#' using \code{BASE_PARTNER}. Supports standard bases (a/c/g/t, A/C/G/T),
#' IUPAC ambiguity codes, and the gap character (".").
#'
#' @param seq Character string. The nucleotide sequence to reverse-complement.
#' @return A character string of the same length as \code{seq}, reversed and
#'   complemented.
#' @examples
#' reverse_complement("ACGT")   # "ACGT"
#' reverse_complement("AACGT")  # "ACGTT"
#' @export
reverse_complement <- function(seq) {
    chars <- strsplit(seq, "", fixed = TRUE)[[1L]]
    rev_chars <- rev(chars)
    comp_chars <- BASE_PARTNER[rev_chars]
    paste0(comp_chars, collapse = "")
}

#' Extended genetic code including degenerate IUPAC codons
#'
#' A named character vector extending \code{GENETIC_CODE} with all degenerate
#' IUPAC codons formed from the 15 symbols in \code{NUCLEOTIDE_CLASSES}. For
#' each degenerate codon, all possible standard expansions are looked up in
#' \code{GENETIC_CODE}:
#' \itemize{
#'   \item If all expansions map to the same amino acid, that amino acid is
#'     assigned.
#'   \item If expansions map to different amino acids, \code{"X"} is assigned.
#' }
#' Standard 64 codons retain their original mappings. Total size is
#' \code{15^3 = 3375} entries (all possible combinations of the 15 IUPAC
#' symbols).
#'
#' Built once at package load time via \code{local()}.
#'
#' @format A named character vector of length 3375.
#' @keywords internal
EXTENDED_GENETIC_CODE <- local({
    bases_plus <- names(NUCLEOTIDE_CLASSES)
    extended   <- GENETIC_CODE  # start from the 64-codon standard code

    for (a in bases_plus) {
        a_bases <- strsplit(NUCLEOTIDE_CLASSES[[a]], "", fixed = TRUE)[[1L]]
        for (b in bases_plus) {
            b_bases <- strsplit(NUCLEOTIDE_CLASSES[[b]], "", fixed = TRUE)[[1L]]
            for (cc in bases_plus) {
                codon <- paste0(a, b, cc)
                if (!is.na(extended[codon])) next  # already present

                c_bases <- strsplit(NUCLEOTIDE_CLASSES[[cc]], "", fixed = TRUE)[[1L]]

                # Collect all AA translations for every standard expansion
                aas <- character(0L)
                for (a1 in a_bases) {
                    for (b1 in b_bases) {
                        for (c1 in c_bases) {
                            aas <- c(aas, GENETIC_CODE[[paste0(a1, b1, c1)]])
                        }
                    }
                }

                extended[codon] <- if (length(unique(aas)) == 1L) aas[[1L]] else "X"
            }
        }
    }

    extended
})

#' Translate a nucleotide sequence to a protein sequence
#'
#' Translates a nucleotide sequence into a single-letter amino acid string
#' using \code{EXTENDED_GENETIC_CODE}. Supports all IUPAC degenerate
#' nucleotide codes. Codons containing \code{"#"} are translated as
#' \code{"#"} (gap indicator used in some CoNGA input files). Codons not
#' found in the extended code are translated as \code{"X"}.
#'
#' @param seq Character string. The nucleotide sequence to translate. May
#'   contain standard or IUPAC-degenerate bases and may be any case (converted
#'   to lowercase internally after offset trimming).
#' @param frame Character string of length 1. Reading frame, e.g. \code{"+1"},
#'   \code{"+2"}, \code{"+3"} for forward strand, or \code{"-1"}, \code{"-2"},
#'   \code{"-3"} for reverse-complement strand. The sign determines strand and
#'   the absolute integer value (1, 2, or 3) determines the 0-based offset into
#'   the sequence before translation begins.
#' @return A character string of amino acids. Length is
#'   \code{floor((nchar(seq) - offset) / 3)}.
#' @examples
#' # Translate from frame +1 (no offset)
#' get_translation("ATGAAATTT", "+1")  # "MKF"
#'
#' # Frame +2 skips the first nucleotide
#' get_translation("AATGAAATTT", "+2")  # "MKF"
#' @export
get_translation <- function(seq, frame = "+1") {
    if (!is.character(frame) || length(frame) != 1L) {
        stop("'frame' must be a single character string, e.g. \"+1\" or \"-2\".")
    }
    strand <- substr(frame, 1L, 1L)
    if (!strand %in% c("+", "-")) {
        stop("'frame' must start with '+' or '-', e.g. \"+1\".")
    }

    if (strand == "-") {
        seq <- reverse_complement(seq)
    }

    offset <- abs(as.integer(frame)) - 1L
    seq <- tolower(substring(seq, offset + 1L))

    naa <- nchar(seq) %/% 3L
    if (naa == 0L) return("")

    # Extract codons as a character vector (vectorised)
    starts <- seq.int(1L, 3L * naa - 2L, by = 3L)
    ends   <- seq.int(3L, 3L * naa,      by = 3L)
    codons <- substring(seq, starts, ends)

    # Map codons: '#' -> '#'; known -> AA; unknown -> 'X'
    has_hash <- grepl("#", codons, fixed = TRUE)
    aas      <- EXTENDED_GENETIC_CODE[codons]
    aas[is.na(aas)]  <- "X"
    aas[has_hash]    <- "#"

    paste0(aas, collapse = "")
}
