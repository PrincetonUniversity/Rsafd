# Send all plots to a throwaway device so tests run headless on every OS.
with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

# Names of objects in the user's workspace, used to detect global side effects.
global_names <- function() setdiff(ls(globalenv(), all.names = TRUE), ".Random.seed")
