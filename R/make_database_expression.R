#' Make a local SQL database of time-course data
#'
#' @param data A data frame with gene names in first column and time-series
#'    expression data in remaining columns.
#' @param path_to_csv Path to a CSV file with gene names in first column and
#'    time-series expression data in remaining columns. Default set to `NULL`
#' @param header Specifies if the file has a header or not.
#'    Either `TRUE` or `FALSE.` Default set to `TRUE`.
#' @param timepoints A numeric vector of time points that correspond to the
#'    time-series data.
#' @param expressed_timepoints_min Minimum number of timepoints that need to
#'    have at least 1 unit expression for a gene to be considered `expressed`.
#' @param cycle The light-dark conditions in which samples were collected.
#'    Either `LD` or `DD`. Default set to `LD`.
#' @param db_name The prefix to be used for the database tables created.
#'    By default, the db_name is set to `database`.
#'
#' @return Creates and connects to a local SQL database with four tables:
#' \enumerate{
#'    \item Input gene expression
#'    \item Expressed genes
#'    \item Log-transformed gene expression
#'    \item Z-scored gene expression
#' }
#'
#' @import dplyr RSQLite glue here utils
#'
#' @export
#'
make_database_expression <- function(data,
                                     path_to_csv = NULL,
                                     header = TRUE,
                                     timepoints,
                                     expressed_timepoints_min = 1,
                                     cycle = "LD",
                                     db_name = "database") {

  ## Load packages ----------
  # pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted)
  library(dplyr)

  ##-##-##-##-##-##-##-##-##-##-##-##-
  # SETUP
  ##-##-##-##-##-##-##-##-##-##-##-##-
  # Create local directory to store database
  if(!dir.exists(here::here(".databases/"))) {
    dir.create(here::here(".databases/"))
    cat("Created .databases/ in project root.")
  }

  # Create empty database
  db_path <- glue::glue(".databases/{db_name}.db")
  my.db <- RSQLite::dbConnect(
    RSQLite::SQLite(),
    here::here(
      db_path
    )
  )

  # Create table names
  tbl_names <- paste(
    db_name,
    c(
      "expression",
      "expressed_genes",
      "log2expression",
      "zscores"
    ),
    sep = "_"
  )

  if(length(RSQLite::dbListTables(my.db)) == 0) {
    # Load data
    if(!is.null(path_to_csv)) {
      # read data
      data <- read.csv(
        path_to_csv,
        header = header,
        stringsAsFactors = FALSE,
        na.strings = c(NA, "", " ")
      ) |>
        as_tibble()
    }

    # Save file to database
    RSQLite::dbWriteTable(
      my.db,
      tbl_names[1],
      data
    )
    cat("Table (", tbl_names[1], ") has been added to", db_path, "\n")

    # Expressed genes ---------------------------------------------------------

    ##-##-##-##-##-##-##-##-##-##-##-##-
    # list of all "Expressed" genes
    ##-##-##-##-##-##-##-##-##-##-##-##-
    expressed <- data |>
      mutate(
        n = rowSums(
          across(
            matches("^ZT|^CT"),
            ~ .x > 1
          ),
          na.rm = TRUE
        ),
        across(
          matches("^ZT|^CT"),
          ~ .x |> round(1)
        )
      ) |>
      mutate(
        expressed = if_else(
          n >= expressed_timepoints_min,
          "yes",
          "no"
        )
      ) |>
      select(
        gene_name,
        expressed
      )

    # Save file to databse
    RSQLite::dbWriteTable(
      my.db,
      tbl_names[2],
      expressed
    )
    cat("Table (", tbl_names[2], ") has been added to", db_path, "\n")

    # log2-expression --------------------------------------------------------

    ##-##-##-##-##-##-##-##-##-##-##-##-
    # log2-transformed expression data
    ##-##-##-##-##-##-##-##-##-##-##-##-
    log2.data <- data |>
      mutate(
        across(
          matches("^ZT|^CT"),
          ~ log2(.x + 1)
        )
      )

    # Save file to database
    RSQLite::dbWriteTable(
      my.db,
      tbl_names[3],
      log2.data
    )
    cat("Table (", tbl_names[3], ") has been added to", db_path, "\n")

    # zscore-log2-expression --------------------------------------------------

    ##-##-##-##-##-##-##-##-##-##-##-##-
    # z-score the log2-transformed exp. data
    ##-##-##-##-##-##-##-##-##-##-##-##-
    zscores.data <- log2.data |>
      tidyr::pivot_longer(
        cols = matches("^ZT|^CT"),
        names_to = "time",
        values_to = "exp"
      ) |>
      group_by(gene_name) |>
      mutate(
        exp = scale(exp)[,1]
      ) |>
      ungroup() |>
      tidyr::pivot_wider(
        names_from = "time",
        values_from = "exp"
      )

    # Save file to database
    RSQLite::dbWriteTable(
      my.db,
      tbl_names[4],
      zscores.data
    )
    cat("Table (", tbl_names[4], ") has been added to", db_path, "\n")
  } else {
    cat("Database exists. Returning link to the database.")
  }

  # export the connection
  my.db

}
