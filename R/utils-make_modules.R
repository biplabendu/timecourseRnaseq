
# A. Data prep ------------------------------------------------------------

#' Log-2 transform data
#'
#' @description
#' Log-2 transform all data columns at once.
#'
#' @note Applies `log2(.x + 1)`.
#'
#' @param data dataframe. One of the columns should contain unique gene
#'    (feature) IDs or names, and other columns should contain expression
#'    levels across timepoints.
#' @param id_column character of length one. Name of the column which
#'    contain unique gene (feature) IDs or names, and for which a
#'    log-transformation should not be attempted.
#' @param log2 logical. If set to FALSE, log-transformation will not be
#'    attempted. Default: `TRUE`.
#'
#' @return dataframe with values in all fields log2-transformed,
#'    except in the `id_column` field.
#'
#' @export
log2_transform_data <- function(data, id_column, log2) {
  chk::chk_data(data)
  chk::chk_logical(log2)

  if (log2) {
    cat("Applying log2-transformation...")
    chk::chk_character(id_column)
    chk::check_names(
      data,
      id_column
    )

    out <- data |>
      mutate(
        across(
          !all_of(id_column),
          ~ log2(.x + 1)
        )
      )
    cat("Done.")
    cat("\n")
  }

  out
}


#' Subset dataframe based on expression thresholds
#'
#' This function filters a data frame to include rows where the expression
#' level exceeds a specified threshold (`min_expression`) in at least a given
#' number of timepoints (`min_timepoints`). It allows specifying the column
#' representing unique IDs (default is `gene_name`).
#'
#' @inheritParams log2_transform_data
#' @param min_expression numeric of length one. Non-negative value specifying
#'    minimum expression required for a timepoint to be counted as _valid_.
#' @param min_timepoints integer of length one. Non-negative value specifying
#'    the minimum number of _valid_ timepoints, where the expression level must
#'    meet or exceed the value specified in `min_expression`.
#'
#' @return dataframe with rows filtered based on the expression threshold
#'   and timepoint count.
#'
#' @note
#' If duplicates exist in the `id_column`, will retain only the first
#' occurrence of each unique ID.
#'
#' @export
subset_data <- function(data,
                        id_column,
                        min_expression,
                        min_timepoints) {

  chk::chk_data(data)

  if (!is.null(min_expression) & !is.null(min_timepoints)) {
    cat("Subsetting data...")
    chk::chk_number(min_expression)
    chk::chk_integer(min_timepoints)
    chk::chk_range(
      min_timepoints, c(1, ncol(data) - 1)
    )
    chk::chk_character(id_column)
    chk::check_names(
      data,
      id_column
    )

    out <- data |>
      mutate(
        n_samples = rowSums(
          across(
            !all_of(id_column),
            ~ .x >= min_expression
          ),
          na.rm = TRUE
        )
      ) |>
      filter(
        n_samples >= min_timepoints
      ) |>
      select(
        - n_samples
      )

    cat("Done.")
    cat("\n")
    cat(
      glue::glue(
        "[ NOTE ]: After subsetting, { nrow(tmp_data) } of { nrow(data) } rows remain."
      )
    )
  }

  out
}


# B. WGCNA ----------------------------------------------------------------

transpose_data <- function(data) {
  datExpr <- as.data.frame(t(data[-1]))
  names(datExpr) <- data$gene_name
  rownames(datExpr) <- names(data)[-1]

  datExpr
}


calculate_gene_gene_sim <- function(data,
                                    method = "kendall",
                                    name = NULL,
                                    overwrite = FALSE,
                                    cache = FALSE,
                                    ...) {
  chk::chk_data(data)
  chk::chk_logical(overwrite)
  chk::chk_logical(cache)

  if(cache == FALSE) {
    cat("Running gene-gene similarity...")
    sim_matrix <- cor(
      datExpr,
      method = method,
      ...
    )
    cat("Done.")
    cat("\n")
  } else {
    chk::chk_scalar(name)

    # set up tmp workspace
    path = ".databases/.tmp/"
    sim_matrix_path = glue::glue(
      "{path}/01_sim_matrix_{name}.RDS"
    )

    if (!dir.exists(here::here(path))) {
      cat("Creating '.databases/.tmp/' for saving cache...")
      dir.create(here::here(path), recursive = TRUE)
      # this step takes time
      cat("Done.")
      cat("\n")
    }

    if (!file.exists(sim_matrix_path) | overwrite == TRUE) {
      cat("Running gene-gene similarity...")
      sim_matrix <- cor(
        datExpr,
        method = method,
        ...
      )
      cat("Done.")
      cat("\n")

      cat("Saving gene-gene similarity matrix...")
      saveRDS(
        sim_matrix,
        file = here::here(
          sim_matrix_path
        )
      )
      cat("Done.")
      cat("\n")
    } else {
      cat("Loading gene-gene similarity matrix from cache...")
      sim_matrix <- readRDS(
        file = here::here(
          sim_matrix_path
        )
      )
      cat("Done.")
    }
  }

  sim_matrix
}


analyze_network_topology <- function(data,
                                     max_power = 21,
                                     height = 0.9,
                                     plot = TRUE) {
  cat("Performing network topology analysis to pick
  soft-thresholding power...\n")

  # Choose a set of soft-thresholding powers
  powers = c(
    c(1:10),
    seq(
      from = 12,
      to = if_else(max_power < 12, 12, max_power),
      by = 3
    )
  )
  # # Call the network topology analysis function
  sft = WGCNA::pickSoftThreshold(
    data,
    powerVector = powers,
    verbose = 0
  )

  if (plot) {
    cat("\n")
    cat("Plotting resutls from the network topology analysis...")
    plot_network_topology(
      data = sft,
      powers = powers,
      height = height
    )
  }

  cat("Done.\n")
  cat(
    glue::glue(
      "[ NOTE, FIGURE ]: Red horizontal line indicates a signed R^2 of {height}"
    )
  )
  cat("\n\n")
  sft
}


perform_hclust <- function(data, plot = FALSE) {

  tree <- hclust(
    as.dist(data),
    method = "average"
  )

  if (plot) {
    # cat("Plotting the resulting clustering tree (dendrogram)")
    plot(
      tree,
      xlab="",
      sub="",
      main = "Hierarchical clustering (method = 'average')\non TOM-based dissimilarity (dissTOM)",
      labels = FALSE,
      hang = 0.04
    )
  }

  tree
}


create_modules <- function(tree,
                           dissTOM,
                           min_module_size,
                           verbose = 0,
                           data,
                           merge_cutoff_similarity = 0.9) {
  modules <- dynamicTreeCut::cutreeDynamic(
    dendro = tree,
    distM = dissTOM,
    method = "hybrid",
    verbose = verbose,
    deepSplit = 3,
    # see WGCNA for more info on tuning parameters
    pamRespectsDendro = FALSE,
    minClusterSize = min_module_size
  ) |>
    WGCNA::labels2colors()

  # Calculate eigengenes
  MEList = WGCNA::moduleEigengenes(
    data,
    colors = modules
  )
  MEs = MEList$eigengenes

  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs, method = "kendall");
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  # Plot the result
  plot(
    METree,
    main = "Initial classification into modules",
    xlab = glue::glue(
      "Red horizontal line shows modules that are ≥
        {merge_cutoff_similarity} similar."
    ),
    ylab = "Dissimilarity (height)",
    ylim = c(0, 0.9),
    sub = "MEDiss = 1-cor(MEs, method = 'kendall')"
  )
  abline(
    h = 1 - merge_cutoff_similarity,
    lwd = 4,
    col = "red",
    lty = 2
  )

  if (!is.null(merge_cutoff_similarity)) {
    cat(
      paste(
        "Merging modules that have a correlation ≥",
        merge_cutoff_similarity,
        "..."
      )
    )
    # Call an automatic merging function
    merged_modules <- WGCNA::mergeCloseModules(
      data,
      modules,
      cutHeight = 1 - merge_cutoff_similarity,
      verbose = 0
    )
    cat("Done.\n")

    # Visualize the merge
    # The merged module colors
    mergedColors = merged_modules$colors;
    # # Eigengenes of the new merged modules:
    cat(
      "[ NOTE, FIGURE ] Plotting identified clusters before and after merging."
    )
    cat("\n\n")
    WGCNA::plotDendroAndColors(
      geneTree,
      cbind(modules, mergedColors),
      main = "Cluster Dendogram (Dynamic Tree Cut)",
      c("Inital classification", "Merged dynamic"),
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = TRUE,
      guideHang = 0.05
    )

    writeLines("Module (cluster) size:")
    table(mergedColors) |> print()

    # return the merged modules
    mergedColors

  } else {
    writeLines("Initial classification into modules (clusters):")
    table(modules) |> print()

    modules
  }

}
