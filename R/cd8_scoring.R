# CD8 logistic regression scoring for paired TCRs.
#
# Ported from Python CoNGA: conga/cd8_scoring.py
#
# Scores each TCR using pre-trained logistic regression models on sorted
# CD4/CD8 bulk TCR data. Alpha and beta chain models are fit independently;
# their scores are averaged (each weighted 0.5).


# ---------------------------------------------------------------------------
# Internal: Load and cache CD8 logistic regression model parameters
# ---------------------------------------------------------------------------

#' Load and cache CD8 logistic regression model parameters
#'
#' Reads \code{inst/extdata/cd8_logreg_params_A.txt} and
#' \code{inst/extdata/cd8_logreg_params_B.txt}. Each file has 7 header lines
#' with integer parameters (key-value pairs), followed by tag-weight pairs for
#' the logistic regression model. The last tag is \code{BIAS}.
#'
#' @return A named list with elements \code{A} and \code{B}, each containing
#'   weights, gene indexers, and model parameters.
#' @keywords internal
.load_cd8_logreg_models <- function() {
    if (!is.null(.tcrdistR_env$cd8_logreg_models)) {
        return(.tcrdistR_env$cd8_logreg_models)
    }

    all_models <- list()

    for (ab in c("A", "B")) {
        fn <- sprintf("cd8_logreg_params_%s.txt", ab)
        filepath <- system.file("extdata", fn, package = "tcrdistR")

        if (!nzchar(filepath)) {
            # Development fallback (devtools::load_all)
            dev_path <- file.path(
                system.file(package = "tcrdistR"),
                "..", "..", "inst", "extdata", fn
            )
            if (file.exists(dev_path)) {
                filepath <- dev_path
            } else {
                stop(sprintf("CD8 model file '%s' not found.", fn))
            }
        }

        param_names <- c("window_size", "min_lenbin", "max_lenbin",
                         "NV", "NJ", "NL", "NC")
        model_params <- list()

        lines <- readLines(filepath, warn = FALSE)
        mtags <- character(0L)
        weights <- numeric(0L)

        for (line in lines) {
            fields <- strsplit(trimws(line), "\\s+")[[1L]]
            if (length(fields) != 2L) next
            tag <- fields[1L]
            value <- fields[2L]
            if (tag %in% param_names) {
                model_params[[tag]] <- as.integer(value)
            } else {
                mtags <- c(mtags, tag)
                weights <- c(weights, as.numeric(value))
            }
        }

        NV <- model_params[["NV"]]
        NJ <- model_params[["NJ"]]
        NL <- model_params[["NL"]]
        NC <- model_params[["NC"]]
        NTOT <- NV + NJ + NL + NC

        if (mtags[length(mtags)] != "BIAS") {
            stop("Expected last tag to be 'BIAS' in ", fn)
        }
        if (length(weights) != NTOT + 1L) {
            stop(sprintf(
                "Expected %d weights (NTOT+1) but got %d in %s",
                NTOT + 1L, length(weights), fn
            ))
        }

        # V-genes: first NV-1 tags (last V slot is UNK)
        vgenes <- mtags[seq_len(NV - 1L)]
        # J-genes: tags NV+1 to NV+NJ-1 (last J slot is UNK)
        jgenes <- mtags[seq(NV + 1L, NV + NJ - 1L)]

        # 0-based indexers (matching Python convention)
        vgene_indexer <- stats::setNames(seq(0L, length(vgenes) - 1L), vgenes)
        jgene_indexer <- stats::setNames(seq(0L, length(jgenes) - 1L), jgenes)

        model_params[["weights"]] <- weights
        model_params[["vgene_indexer"]] <- vgene_indexer
        model_params[["jgene_indexer"]] <- jgene_indexer

        all_models[[ab]] <- model_params
    }

    .tcrdistR_env$cd8_logreg_models <- all_models
    all_models
}


# ---------------------------------------------------------------------------
# Internal: TCR Feature Encoding
# ---------------------------------------------------------------------------

#' Encode a single TCR chain as a feature vector for CD8 logistic regression
#'
#' Creates a one-hot/k-hot feature vector encoding V-gene, J-gene, CDR3
#' length, and CDR3 amino acid composition at N-terminal, C-terminal, and
#' middle positions.
#'
#' @param vgene Character. V-gene name with allele.
#' @param jgene Character. J-gene name with allele.
#' @param cdr3 Character. CDR3 amino acid sequence.
#' @param model_params List. Model parameters from
#'   \code{.load_cd8_logreg_models}.
#' @return Numeric vector of length NTOT + 1 (features + bias).
#' @keywords internal
.encode_single_chain_tcr <- function(vgene, jgene, cdr3, model_params) {
    window_size <- model_params[["window_size"]]
    min_lenbin <- model_params[["min_lenbin"]]
    max_lenbin <- model_params[["max_lenbin"]]
    vgene_indexer <- model_params[["vgene_indexer"]]
    jgene_indexer <- model_params[["jgene_indexer"]]
    NV <- length(vgene_indexer) + 1L
    NJ <- length(jgene_indexer) + 1L
    NL <- max_lenbin - min_lenbin + 1L
    NC <- 20L * (2L * window_size + 1L)
    NTOT <- NV + NJ + NL + NC

    x <- numeric(NTOT + 1L)

    # V-gene: strip allele, look up 0-based index (UNK = NV-1)
    vgene_stripped <- sub("\\*.*$", "", vgene)
    iv <- if (vgene_stripped %in% names(vgene_indexer)) {
        vgene_indexer[[vgene_stripped]]
    } else {
        NV - 1L
    }
    x[iv + 1L] <- 1.0

    # J-gene: strip allele, look up 0-based index (UNK = NJ-1)
    jgene_stripped <- sub("\\*.*$", "", jgene)
    ij <- if (jgene_stripped %in% names(jgene_indexer)) {
        jgene_indexer[[jgene_stripped]]
    } else {
        NJ - 1L
    }
    x[NV + ij + 1L] <- 1.0

    # CDR3 length bin
    cdr3_len <- nchar(cdr3)
    lenbin <- if (cdr3_len <= min_lenbin) {
        min_lenbin
    } else if (cdr3_len >= max_lenbin) {
        max_lenbin
    } else {
        cdr3_len
    }
    il <- lenbin - min_lenbin
    x[NV + NJ + il + 1L] <- 1.0

    # CDR3 amino acid encoding: N-terminal + C-terminal windows + middle k-hot
    nterm <- min(window_size, cdr3_len %/% 2L)
    cterm <- min(window_size, cdr3_len - nterm)

    cdr3_chars <- strsplit(cdr3, "")[[1L]]

    for (i in seq(0L, window_size - 1L)) {
        # N-terminal window
        if (i < nterm) {
            aa <- cdr3_chars[i + 1L]
            aa_idx <- match(aa, AMINO_ACIDS) - 1L
            if (!is.na(aa_idx)) {
                pos <- NV + NJ + NL + 20L * i + aa_idx + 1L
                x[pos] <- x[pos] + 1.0
            }
        }
        # C-terminal window (reading from end)
        if (i < cterm) {
            aa <- cdr3_chars[cdr3_len - i]
            aa_idx <- match(aa, AMINO_ACIDS) - 1L
            if (!is.na(aa_idx)) {
                pos <- NV + NJ + NL + 20L * window_size + 20L * i + aa_idx + 1L
                x[pos] <- x[pos] + 1.0
            }
        }
    }

    # Middle portion: k-hot encoding
    if (nterm + 1L <= cdr3_len - cterm) {
        middle_chars <- cdr3_chars[seq(nterm + 1L, cdr3_len - cterm)]
        for (aa in middle_chars) {
            aa_idx <- match(aa, AMINO_ACIDS) - 1L
            if (!is.na(aa_idx)) {
                pos <- NV + NJ + NL + 40L * window_size + aa_idx + 1L
                x[pos] <- x[pos] + 1.0
            }
        }
    }

    # Bias term (last element)
    x[NTOT + 1L] <- 1.0

    x
}


# ---------------------------------------------------------------------------
# Public: CD8 Logistic Regression Scoring
# ---------------------------------------------------------------------------

#' CD8 logistic regression score for paired TCRs
#'
#' Scores each TCR using a logistic regression model trained on sorted CD4/CD8
#' bulk TCR data. The alpha-chain and beta-chain models are fit independently
#' and their scores are averaged (each weighted 0.5).
#'
#' The model encodes V-gene, J-gene, CDR3 length, and positional amino acid
#' features using one-hot and k-hot encodings, then computes a dot product
#' with pre-trained weights.
#'
#' @param tcr_df A data.frame with columns \code{va}, \code{ja}, \code{cdr3a},
#'   \code{vb}, \code{jb}, \code{cdr3b}. Gene names should include allele
#'   suffixes (e.g. \code{"TRAV1-2*01"}).
#' @param use_sigmoid Logical. If \code{TRUE}, apply the sigmoid function
#'   \eqn{1 / (1 + exp(-x))} to each chain's score before averaging.
#'   Default is \code{FALSE}.
#'
#' @return A numeric vector of CD8 logistic regression scores, one per row
#'   of \code{tcr_df}.
#'
#' @examples
#' \donttest{
#' tcr_df <- data.frame(
#'     va = "TRAV1-2*01", ja = "TRAJ33*01",
#'     cdr3a = "CAVMDSSYKLIF",
#'     vb = "TRBV6-4*01", jb = "TRBJ2-1*01",
#'     cdr3b = "CASSLAPGATNEKLFF",
#'     stringsAsFactors = FALSE
#' )
#' scores <- make_cd8_score_table_column(tcr_df)
#' }
#' @export
make_cd8_score_table_column <- function(tcr_df, use_sigmoid = FALSE) {
    if (!is.data.frame(tcr_df)) {
        stop("make_cd8_score_table_column: tcr_df must be a data.frame")
    }
    required_cols <- c("va", "ja", "cdr3a", "vb", "jb", "cdr3b")
    missing_cols <- setdiff(required_cols, names(tcr_df))
    if (length(missing_cols) > 0L) {
        stop("make_cd8_score_table_column: tcr_df is missing required ",
             "columns: ", paste(missing_cols, collapse = ", "))
    }

    models <- .load_cd8_logreg_models()
    n <- nrow(tcr_df)
    score_totals <- numeric(n)

    chain_cols <- list(
        A = c(v = "va", j = "ja", cdr3 = "cdr3a"),
        B = c(v = "vb", j = "jb", cdr3 = "cdr3b")
    )

    for (ab in c("A", "B")) {
        model <- models[[ab]]
        weights <- model[["weights"]]
        cols <- chain_cols[[ab]]

        tcr_matrix <- matrix(0.0, nrow = n, ncol = length(weights))
        for (i in seq_len(n)) {
            tcr_matrix[i, ] <- .encode_single_chain_tcr(
                tcr_df[[cols["v"]]][i],
                tcr_df[[cols["j"]]][i],
                tcr_df[[cols["cdr3"]]][i],
                model
            )
        }

        ab_scores <- as.numeric(tcr_matrix %*% weights)

        if (use_sigmoid) {
            ab_scores <- 1.0 / (1.0 + exp(-ab_scores))
        }

        score_totals <- score_totals + 0.5 * ab_scores
    }

    score_totals
}
