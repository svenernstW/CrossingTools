# --- check_macos.R ---------------------------------------------------------
# One-shot setup for macOS CI checks + safe RcppParallel build flags
# Run from your package root:  source("check_macos.R")

message("\n=== CrossingTools: macOS CI + build setup ===\n")

pkg_root <- normalizePath(".", mustWork = TRUE)
desc_path <- file.path(pkg_root, "DESCRIPTION")
makevars_path <- file.path(pkg_root, "src", "Makevars")
makevars_win_path <- file.path(pkg_root, "src", "Makevars.win")
workflow_dir <- file.path(pkg_root, ".github", "workflows")
workflow_path <- file.path(workflow_dir, "R-CMD-check.yml")

stopifnot(file.exists(desc_path))

dir.create(file.path(pkg_root, "src"), showWarnings = FALSE, recursive = TRUE)
dir.create(workflow_dir, showWarnings = FALSE, recursive = TRUE)

backup <- function(path) {
  if (file.exists(path)) {
    file.copy(path, paste0(path, ".bak"), overwrite = TRUE)
  }
}

write_text <- function(path, txt) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(txt, con = path, useBytes = TRUE)
}

# -------------------- 1) Write GitHub Actions workflow ---------------------
workflow_yaml <- "
name: R-CMD-check

on:
  push:
  pull_request:

jobs:
  check:
    strategy:
      matrix:
        os: [macos-14, macos-13, ubuntu-latest, windows-latest]
        r: [release]
    runs-on: ${{ matrix.os }}

    steps:
      - uses: actions/checkout@v4

      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: ${{ matrix.r }}

      - uses: r-lib/actions/setup-pandoc@v2

      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::rcmdcheck
          needs: check

      - name: R CMD check
        run: |
          Rscript -e 'rcmdcheck::rcmdcheck(args=\"--as-cran\", error_on=\"warning\")'
"

backup(workflow_path)
write_text(workflow_path, workflow_yaml)
message("✓ Wrote workflow: .github/workflows/R-CMD-check.yml")

# -------------------- 2) Patch DESCRIPTION --------------------------------
desc <- readLines(desc_path, warn = FALSE)

has_field <- function(lines, field) any(grepl(paste0("^", field, ":"), lines, ignore.case = TRUE))

ensure_field_contains <- function(lines, field, items) {
  ix <- grep(paste0("^", field, ":"), lines, ignore.case = TRUE)
  if (length(ix) == 0) {
    # add new field at end
    value <- paste(items, collapse = ", ")
    return(c(lines, paste0(field, ": ", value)))
  } else {
    line <- lines[ix[1]]
    # preserve original capitalization of field label
    label <- sub(":.*$", "", line)
    rest <- sub("^[^:]*:\\s*", "", line)
    # split by comma, trim
    parts <- trimws(unlist(strsplit(rest, ",")))
    for (it in items) {
      if (!any(tolower(parts) == tolower(it))) {
        parts <- c(parts, it)
      }
    }
    new_line <- paste0(label, ": ", paste(unique(parts[nzchar(parts)]), collapse = ", "))
    lines[ix[1]] <- new_line
    return(lines)
  }
}

ensure_or_set <- function(lines, field, value) {
  ix <- grep(paste0("^", field, ":"), lines, ignore.case = TRUE)
  if (length(ix) == 0) {
    c(lines, paste0(field, ": ", value))
  } else {
    label <- sub(":.*$", "", lines[ix[1]])
    lines[ix[1]] <- paste0(label, ": ", value)
    lines
  }
}

backup(desc_path)

desc <- ensure_field_contains(desc, "LinkingTo", c("Rcpp", "RcppParallel"))
desc <- ensure_field_contains(desc, "Imports",   c("Rcpp", "RcppParallel"))
desc <- ensure_or_set(desc, "SystemRequirements", "C++17")

writeLines(desc, desc_path, useBytes = TRUE)
message("✓ Ensured DESCRIPTION has LinkingTo/Imports (Rcpp, RcppParallel) and SystemRequirements: C++17")

# -------------------- 3) Write Makevars / Makevars.win --------------------
makevars_txt <- '
CXX_STD = CXX17

# Ask RcppParallel for the right flags/libs at build time
PKG_CPPFLAGS = $(shell "${R_HOME}/bin/Rscript" -e "RcppParallel::RcppParallelCxxFlags()")
PKG_LIBS     = $(shell "${R_HOME}/bin/Rscript" -e "RcppParallel::RcppParallelLibs()")

# Keep these only if *your code* also uses OpenMP (harmless on macOS if unused)
PKG_CXXFLAGS += $(SHLIB_OPENMP_CXXFLAGS)
PKG_LIBS     += $(SHLIB_OPENMP_LDFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
'

backup(makevars_path);     write_text(makevars_path, makevars_txt)
backup(makevars_win_path); write_text(makevars_win_path, makevars_txt)
message("✓ Wrote src/Makevars and src/Makevars.win")

# -------------------- 4) Friendly next steps ------------------------------
cat("\nNext steps:\n",
    "1) In R, run:\n",
    "   Rcpp::compileAttributes(); devtools::document(); devtools::clean_dll(); devtools::load_all()\n",
    "   # Optional full check: devtools::check(args='--as-cran')\n\n",
    "2) Open GitHub Desktop:\n",
    "   - You should see these file changes staged:\n",
    "     * .github/workflows/R-CMD-check.yml\n",
    "     * DESCRIPTION (updated)\n",
    "     * src/Makevars\n",
    "     * src/Makevars.win\n",
    "   - Commit with a message like: 'Set up macOS CI + RcppParallel flags'\n",
    "   - Push to GitHub.\n\n",
    "3) On GitHub → Actions → R-CMD-check:\n",
    "   - Watch jobs for macos-14 (ARM) and macos-13 (Intel).\n",
    "   - If they pass, your package installs/loads on macOS.\n",
    "   - If a macOS job fails, click it and share the error log.\n\n",
    "Tip for users who still fail to install on Mac:\n",
    "   install.packages(c('Rcpp','RcppParallel'));\n",
    "   remove.packages('CrossingTools'); unlink(Sys.glob(file.path(.libPaths()[1],'00LOCK*')), recursive=TRUE, force=TRUE);\n",
    "   install.packages('CrossingTools', type='source')\n\n",
    "Done ✅\n")
# --------------------------------------------------------------------------
