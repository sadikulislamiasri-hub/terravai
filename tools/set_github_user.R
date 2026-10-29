# Replace the GITHUB_USER placeholder in all package files.
# Usage (from the package folder):  Rscript tools/set_github_user.R your-username
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !nzchar(args)) stop("Give your GitHub username as the only argument.")
files <- c("DESCRIPTION", "README.md", "inst/CITATION", "CITATION.cff")
for (f in files) {
  x <- readLines(f, warn = FALSE)
  writeLines(gsub("GITHUB_USER", args, x, fixed = TRUE), f)
}
message("Replaced GITHUB_USER with '", args, "' in: ", paste(files, collapse = ", "))
