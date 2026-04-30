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


#' Human influenza TCR dataset from VDJdb
#'
#' A dataset of 2271 paired alpha-beta human T-cell receptors specific to
#' 4 influenza epitopes, compiled from the VDJdb database (September 2020
#' release). This dataset covers both MHC class I and class II restricted
#' responses from 63 subjects and is useful for benchmarking human TCR
#' distance calculations.
#'
#' @format A data.frame with 2271 rows and 13 columns:
#' \describe{
#'   \item{subject}{Subject identifier from the original study (may be
#'     \code{NA} if not annotated).}
#'   \item{epitope}{Epitope peptide sequence (e.g., \code{"GILGFVFTL"}).}
#'   \item{epitope_gene}{Source gene of the epitope (e.g., \code{"M1"}).}
#'   \item{mhc_a}{MHC alpha chain allele (e.g., \code{"HLA-A*02"}).}
#'   \item{mhc_b}{MHC beta chain allele (e.g., \code{"B2M"}).}
#'   \item{mhc_class}{MHC class: \code{"MHCI"} or \code{"MHCII"}.}
#'   \item{count}{Clone count (set to 1 for all entries).}
#'   \item{va}{V-alpha gene with allele (e.g., \code{"TRAV12-2*01"}).}
#'   \item{ja}{J-alpha gene with allele (e.g., \code{"TRAJ33*01"}).}
#'   \item{cdr3a}{CDR3-alpha amino acid sequence.}
#'   \item{vb}{V-beta gene with allele (e.g., \code{"TRBV19*01"}).}
#'   \item{jb}{J-beta gene with allele (e.g., \code{"TRBJ2-7*01"}).}
#'   \item{cdr3b}{CDR3-beta amino acid sequence.}
#' }
#'
#' @source VDJdb: \url{https://vdjdb.cdr3.net/}
#'
#'   Bagaev et al. (2020). VDJdb in 2019: database extension, new analysis
#'   infrastructure and a T-cell receptor motif compendium. \emph{Nucleic
#'   Acids Research}, 48(D1), D1057--D1062. \doi{10.1093/nar/gkz874}
#'
#' @examples
#' data(flu)
#' dim(flu)         # 2271 x 13
#' table(flu$epitope)
"flu"
