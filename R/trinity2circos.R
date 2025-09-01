#' Run Trinity-to-Circos preprocessing (bash helper)
#'
#' Calls the bash script included in the package to convert a Trinity-assembled
#' FASTA to a preBED-like table for Circos.
#'
#' @param input_fasta Character. Path to the Trinity FASTA file.
#' @param remove_tmp Logical. If `TRUE`, delete temporary files without prompting.
#' @param echo Logical. If `TRUE`, show bash output.
#'
#' @details
#' The helper script is installed under the package's `exec/` directory and
#' is located via `system.file("exec", "trinity2circos.sh", package = packageName())`.
#' Requires `bash`, `awk`, `sed`, `grep`, `sort` in `PATH`.
#'
#' @return (Invisibly) the exit status from `system2()` (0 = success).
#' @examples
#' \dontrun{
#' run_trinity2circos("Trinity.fasta", remove_tmp = FALSE)
#' }
#' @seealso [circosplot_degs]
#' @export
run_trinity2circos <- function(input_fasta, remove_tmp = FALSE, echo = TRUE) {
  stopifnot(is.character(input_fasta), length(input_fasta) == 1L, file.exists(input_fasta))

  # localizar bash
  bash <- Sys.which("bash")
  if (!nzchar(bash)) stop("Could not find 'bash' in PATH.")

  # localizar script dentro del paquete
  pkg <- utils::packageName()
  sh  <- system.file("exec", "trinity_nucl_to_preBED.sh", package = pkg)
  if (!nzchar(sh) || !file.exists(sh)) {
    stop("Bash helper not found in package exec/: ", sh)
  }

  # args para el script: -i <archivo>
  args <- c(shQuote(sh), "-i", shQuote(normalizePath(input_fasta)))

  # variable de entorno para limpieza no interactiva (tu bash debe respetarla)
  env <- character()
  if (isTRUE(remove_tmp)) env <- c(env, "RUN_REMOVE_TMP=1")

  status <- system2(command = bash,
                    args    = args,
                    stdout  = if (echo) "" else FALSE,
                    stderr  = if (echo) "" else FALSE,
                    env     = env)
  invisible(status)
}
