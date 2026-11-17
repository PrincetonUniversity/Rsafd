# Install the Imports/Depends/LinkingTo/Suggests of a DESCRIPTION file, plus any
# extra packages named on the command line. Fails if anything cannot be loaded.
#   Rscript scripts/install-deps.R path/to/DESCRIPTION [extra packages...]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || !file.exists(args[[1L]])) {
  stop("usage: Rscript install-deps.R <DESCRIPTION> [extra packages...]", call. = FALSE)
}
fields <- c("Depends", "Imports", "LinkingTo", "Suggests")
desc <- read.dcf(args[[1L]], fields = fields)
spec <- unlist(strsplit(paste(desc[!is.na(desc)], collapse = ","), ","))
pkgs <- trimws(sub("\\(.*$", "", spec))
base <- rownames(utils::installed.packages(priority = c("base", "recommended")))
pkgs <- setdiff(unique(c(pkgs[nzchar(pkgs)], args[-1L])), c("R", base))
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) utils::install.packages(missing)
load_error <- function(p) {
  tryCatch({ loadNamespace(p); NA_character_ }, error = function(e) conditionMessage(e))
}
errors <- vapply(pkgs, load_error, character(1))
failed <- errors[!is.na(errors)]
if (length(failed)) {
  for (p in names(failed)) message(p, ": ", failed[[p]])
  stop("could not install: ", paste(names(failed), collapse = ", "), call. = FALSE)
}
cat("dependencies ready:", paste(sort(pkgs), collapse = ", "), "\n")
