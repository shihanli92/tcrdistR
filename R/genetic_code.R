# Ported from Python CoNGA: genetic_code.py

#' Standard genetic code (codon to amino acid)
#'
#' Named character vector mapping all 64 lowercase DNA codons to single-letter
#' amino acid codes. Stop codons map to \code{"*"}.
#'
#' @format A named character vector of length 64.
#' @keywords internal
GENETIC_CODE <- c(
    ttt = "F", ttc = "F", tta = "L", ttg = "L",
    tct = "S", tcc = "S", tca = "S", tcg = "S",
    tat = "Y", tac = "Y", taa = "*", tag = "*",
    tgt = "C", tgc = "C", tga = "*", tgg = "W",
    ctt = "L", ctc = "L", cta = "L", ctg = "L",
    cct = "P", ccc = "P", cca = "P", ccg = "P",
    cat = "H", cac = "H", caa = "Q", cag = "Q",
    cgt = "R", cgc = "R", cga = "R", cgg = "R",
    att = "I", atc = "I", ata = "I", atg = "M",
    act = "T", acc = "T", aca = "T", acg = "T",
    aat = "N", aac = "N", aaa = "K", aag = "K",
    agt = "S", agc = "S", aga = "R", agg = "R",
    gtt = "V", gtc = "V", gta = "V", gtg = "V",
    gct = "A", gcc = "A", gca = "A", gcg = "A",
    gat = "D", gac = "D", gaa = "E", gag = "E",
    ggt = "G", ggc = "G", gga = "G", ggg = "G"
)

#' Reverse genetic code (amino acid to codons)
#'
#' Named list mapping each single-letter amino acid code (and \code{"*"} for
#' stop codons) to a character vector of the corresponding lowercase DNA
#' codons.
#'
#' @format A named list of length 21 (20 amino acids + stop).
#' @keywords internal
REVERSE_GENETIC_CODE <- list(
    "*" = c("taa", "tag", "tga"),
    A   = c("gca", "gcc", "gcg", "gct"),
    C   = c("tgc", "tgt"),
    D   = c("gac", "gat"),
    E   = c("gaa", "gag"),
    F   = c("ttc", "ttt"),
    G   = c("gga", "ggc", "ggg", "ggt"),
    H   = c("cac", "cat"),
    I   = c("ata", "atc", "att"),
    K   = c("aaa", "aag"),
    L   = c("cta", "ctc", "ctg", "ctt", "tta", "ttg"),
    M   = c("atg"),
    N   = c("aac", "aat"),
    P   = c("cca", "ccc", "ccg", "cct"),
    Q   = c("caa", "cag"),
    R   = c("aga", "agg", "cga", "cgc", "cgg", "cgt"),
    S   = c("agc", "agt", "tca", "tcc", "tcg", "tct"),
    T   = c("aca", "acc", "acg", "act"),
    V   = c("gta", "gtc", "gtg", "gtt"),
    W   = c("tgg"),
    Y   = c("tac", "tat")
)

#' Degenerate codon representations per amino acid
#'
#' Named list mapping each single-letter amino acid code to a character vector
#' of IUPAC degenerate codon strings. Amino acids whose codons span multiple
#' two-letter prefixes have multiple entries (e.g., L = c("ctn", "ttr")).
#' Stop codons are excluded, matching the Python source.
#'
#' The degeneracy mapping used for the third nucleotide position:
#' \itemize{
#'   \item a -> a, c -> c, g -> g, t -> t
#'   \item ct -> y, ag -> r
#'   \item act -> h
#'   \item acgt -> n
#' }
#'
#' @format A named list of length 20.
#' @keywords internal
AA2DEGENERATE_CODONS <- local({
    # Partial degeneracy map (only what's needed, matching Python degmap)
    degmap <- list(
        "a"    = "a",
        "c"    = "c",
        "g"    = "g",
        "t"    = "t",
        "ct"   = "y",
        "ag"   = "r",
        "act"  = "h",
        "acgt" = "n"
    )

    result <- list()
    for (aa in names(REVERSE_GENETIC_CODE)) {
        if (aa == "*") next
        codons <- REVERSE_GENETIC_CODE[[aa]]
        # Group by first two nucleotides
        prefixes <- unique(substr(codons, 1, 2))
        degen_codons <- character(0)
        for (xy in sort(prefixes)) {
            third_nucs <- sort(unique(substr(
                codons[substr(codons, 1, 2) == xy], 3, 3
            )))
            third_key <- paste0(third_nucs, collapse = "")
            degen_codons <- c(degen_codons, paste0(xy, degmap[[third_key]]))
        }
        result[[aa]] <- degen_codons
    }
    result
})
