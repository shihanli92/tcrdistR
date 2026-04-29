#' DASH dataset: paired alpha-beta mouse TCRs across 7 epitopes
#'
#' A dataset of 1924 paired alpha-beta T-cell receptors from mice responding
#' to 7 viral epitopes, collected from 78 subjects. This is the benchmark
#' dataset from Dash et al. (2017) and is widely used for evaluating TCR
#' distance metrics.
#'
#' @format A data.frame with 1924 rows and 12 columns:
#' \describe{
#'   \item{subject}{Subject identifier (e.g., \code{"mouse_subject0050"}).}
#'   \item{epitope}{Epitope specificity: F2, M38, M45, m139, NP, PA, or PB1.}
#'   \item{count}{Clone count (number of cells observed for this clonotype).}
#'   \item{va}{V-alpha gene with allele (e.g., \code{"TRAV7-3*01"}).}
#'   \item{ja}{J-alpha gene with allele (e.g., \code{"TRAJ33*01"}).}
#'   \item{cdr3a}{CDR3-alpha amino acid sequence (e.g., \code{"CAVSLDSNYQLIW"}).}
#'   \item{cdr3a_nucseq}{CDR3-alpha nucleotide sequence.}
#'   \item{vb}{V-beta gene with allele (e.g., \code{"TRBV13-1*01"}).}
#'   \item{jb}{J-beta gene with allele (e.g., \code{"TRBJ2-3*01"}).}
#'   \item{cdr3b}{CDR3-beta amino acid sequence (e.g., \code{"CASSDFDWGGDAETLYF"}).}
#'   \item{cdr3b_nucseq}{CDR3-beta nucleotide sequence.}
#'   \item{clone_id}{Unique clone identifier.}
#' }
#'
#' @source Dash et al. (2017). Quantifiable predictive features define
#'   epitope-specific T cell receptor repertoires. \emph{Nature}, 547, 89--93.
#'   \doi{10.1038/nature22383}
#'
#' @examples
#' data(dash)
#' dim(dash)       # 1924 x 12
#' table(dash$epitope)
"dash"
