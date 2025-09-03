
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
log2_transform_data <- function(data,
                                id_column,
                                log2 = TRUE) {
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
    chk::chk_number(min_timepoints)
    chk::chk_lte(
      min_timepoints,
      ncol(data) - 1,
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
        "[ NOTE ]: After subsetting, { nrow(out) } of { nrow(data) } rows remain."
      )
    )
  }

  out
}


# B. WGCNA ----------------------------------------------------------------

#' Transpose data from a matrix to a data frame
#'
#' This function transposes the input data frame so that each column represent
#' genes, rows represent sequential time points.
#'
#' @inheritParams log2_transform_data
#'
#' @return dataframe, with rownames (samples or timepoints)
#'
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
      data,
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
        data,
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

#' Perform network topology analysis to pick soft-thresholding power.
#'
#' This function analyzes the network topology and uses it to determine optimal
#' soft-thresholding power. The results are also visualized as a plot.
#'
#' @inheritParams check_sample_quality
#' @param max_power Integer specifying the maximum possible value for the
#'    soft-thresholding power (default = 21)
#' @param height Real number specifying the signed R^2 of the
#'    correlation between sample means and the corresponding row means in the
#'    network topology matrix (default = 0.9)
#' @param plot Logical indicating whether to plot the results from the
#'    network topology analysis (default = TRUE).
#'    If true calls `plot_network_topology()`.
#'
#' @return listObject containing the optimized soft-thresholding power
#'
#' @export
#'
analyze_network_topology <- function(data,
                                     max_power = 21,
                                     height = 0.9,
                                     plot = TRUE) {
  chk::chk_range(height, c(0.5, 1))
  chk::chk_number(max_power)
  chk::chk_range(max_power, c(2, 30))

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

#' Perform hierarchical clustering.
#'
#' Performs average linkage clustering on the given data matrix.
#'
#' @param data Matrix of gene expression values
#' @param plot Logical indicating whether to display the resulting clustering
#' tree (default = FALSE)
#'
#' @return Hierarchical clustering tree object
#'
#' @keywords internal
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
                           merge_cutoff_similarity = 0.9,
                           plot = FALSE) {
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

  if (plot) {
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
      lty = 1
    )
  }

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
      tree,
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
    list(
      "colors" = mergedColors,
      "modules" = merged_modules
    )

  } else {
    writeLines("Initial classification into modules (clusters):")
    table(modules) |> print()

    modules
  }

}

create_modules_auto <- function(tree,
                                dissTOM,
                                min_module_size,
                                data,
                                merge_cutoff_similarity,
                                max_modules) {

  init_modules <- create_modules(
    tree = tree,
    dissTOM = dissTOM,
    data = data,
    merge_cutoff_similarity = merge_cutoff_similarity,
    min_module_size = min_module_size
  )

  n_modules <- ncol(init_modules$modules$newMEs)
  cutoff <- merge_cutoff_similarity

  while (n_modules > max_modules) {
    cutoff <- cutoff - 0.05
    merge <- create_modules(
      tree = tree,
      dissTOM = dissTOM,
      data = data,
      merge_cutoff_similarity = cutoff,
      min_module_size = min_module_size
    )
    n_modules <- ncol(merge$modules$newMEs)
  }

  cat(
    "\nCutoff used:", cutoff,
    "\nNumber of modules identified:", n_modules,
    "\n\n"
  )

  merge
}

tidy_modules <- function(merged_modules, mapping_tbl, data) {

  # Make a list that returns gene names for a given cluster
  module_genes <- list()

  # modules.to.exclude <- c(paste0("C",c(2,5,6,7,10:17,19)))
  modules.to.exclude <- c("")
  which.modules <- mapping_tbl %>%
    filter(!new_labels %in% modules.to.exclude) %>%
    pull(old_labels)
  which.labels <- mapping_tbl %>%
    filter(!new_labels %in% modules.to.exclude) %>%
    pull(new_labels)

  # Get the genes from each of the modules
  for (i in 1:length(which.modules)) {
    # which color
    mod.color = as.character(which.modules[[i]])
    # subset
    module_genes[[i]] <- names(data)[
      which(merged_modules == mod.color)
    ]
    names(module_genes)[[i]] <- as.character(
      which.labels[[i]]
    )
  }

  # # check the result | works
  # names(module_genes)
  # module_genes['C22']

  # [13 Dec 2021]
  # Save a csv with the module identity information for all genes used in building the GCN
  # make a dataframe with gene_name and module_identity
  for (i in 1:length(module_genes)){
    if (i == 1){
      dat_module_gene <- data.frame(
        gene_name = module_genes[[i]],
        module_identity = as.character(names(module_genes)[i])
      )
    }
    else{
      foo <- data.frame(
        gene_name = module_genes[[i]],
        module_identity = as.character(names(module_genes)[i])
      )
      dat_module_gene <- rbind(dat_module_gene, foo)
    }

  }

  dat_module_gene |>
    as_tibble() |>
    left_join(
      mapping_tbl,
      join_by(module_identity == new_labels)
    )
}

# C. INTERNAL ---------------------------------------------------------

#' Estimate minimum expression values across all columns except the ID column.
#'
#' This function calculates the minimum expression value in each column of
#' the input data, excluding the specified ID column.
#' The results are then converted to ceiling values.
#'
#' @inheritParams log2_transform_data
#'
#' @return numeric of length one
#'
#' @keywords internal
#'
estimate_min_expression <- function(data, id_column) {
  data |>
    select(
      !all_of(
        id_column
      )
    ) |>
    summarize(across(everything(), ~ mean(.x, na.rm = TRUE))) |>
    rowMeans() |>
    ceiling()
}

#' Check sample quality based on WGCNA gene expression data.
#'
#' This function uses the WGCNA package to check if all samples are of
#' good quality.
#'
#' @param data Data frame containing gene expression data, where
#' each column represent genes, rows represent sequential time points.
#' Intended for output of `transpose_data()`.
#'
#' @seealso [transpose_data()]
#'
#' @return stops if bad samples found, prints result in console
#'
#' @keywords internal
check_sample_quality <- function(data) {
  # check for bad samples
  gsg = WGCNA::goodSamplesGenes(data, verbose = 3);
  if(gsg$allOK == TRUE) {
    cat("All okay!")
  } else {
    paste("Bad sample(s) detected!")
    stop()
  }
}

#' Visualize log-transformed expression data across samples.
#'
#' This function plots log-transformed gene expression values for each sample,
#' with the x-axis representing log2 expression and the y-axis representing
#' density of a given expression across genes.
#'
#' @inheritParams check_sample_quality
#'
#' @return a ggplot
#' @keywords internal
plot_sample_expression <- function(data) {

  nSamples = nrow(data)

  writeLines("Visualizing the log-transformed data")

  data |>
    tibble::rownames_to_column(
      "sample"
    ) |>
    as_tibble() |>
    tidyr::pivot_longer(
      cols = !sample,
      names_to = "gene_id",
      values_to = "log2_fpkm"
    ) |>
    ggplot(
      aes(
        x = log2_fpkm,
        # color = sample,
        fill = sample
      )
    ) +
    geom_density(
      position = "stack"
    ) +
    theme_bw(20) +
    scale_fill_manual(
      values = viridis::viridis(nSamples)
    ) +
    labs(
      x = "Expression (log2)"
    ) +
    theme(
      legend.position = "bottom",
      legend.justification = "right"
    ) +
    guides(
      fill = guide_legend(
        nrow = 3,
        byrow=TRUE
      )
    )
}

#' Visualize network topology analysis results.
#'
#' This function plots the scale-free fit index and mean connectivity
#' as a function of the soft-thresholding power.
#'
#' @param data Object containing WGCNA network topology data
#' @param powers Vector specifying soft-thresholding powers used for analysis
#' @param height Real number specifying the signed R^2 threshold
#'
#' @return plot network topology analyses results
#' @keywords internal
#'
#' @seealso [analyze_network_topology()]
#'
plot_network_topology <- function(data,
                                  powers,
                                  height) {
  sft = data
  par(mfrow = c(1,2));
  cex1 = 1.1;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(
    sft$fitIndices[,1],
    -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)",
    ylab="Scale Free Topology\nModel Fit, (signed R^2)",
    type="n",
    main = paste("Scale independence")
  );
  # this line corresponds to using an R^2 cut-off of h
  abline(
    h = height,
    col="pink",
    lwd = 12,
    lty = 1
  );
  text(
    sft$fitIndices[,1],
    -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,
    cex=cex1,
    col="purple"
  );
  # Mean connectivity as a function of the soft-thresholding power
  plot(
    sft$fitIndices[,1],
    sft$fitIndices[,5],
    xlab="Soft Threshold (power)",
    ylab="Mean Connectivity",
    type="n",
    main = paste("Mean connectivity")
  )
  text(
    sft$fitIndices[,1],
    sft$fitIndices[,5],
    labels=powers,
    cex=cex1,
    col="purple"
  )
  par(mfrow = c(1,1))
}

#' Estimate optimal soft-thresholding power.
#'
#' Estimate the optimal soft-thresholding power based on SFT R.sq values.
#'
#' @param sft Object containing network topology data and fitted indices
#'
#' @return integer, between 6 and 16
#' @keywords internal
#'
#' @seealso [analyze_network_topology()]
#'
estimate_soft_power <- function(sft) {
  my_estimate <- sft$fitIndices |>
    filter(
      SFT.R.sq > 0.8 &
        SFT.R.sq <= 0.9,
      Power > 6 &
        Power < 16
    ) |>
    arrange(
      Power
    ) |>
    head(1) |>
    pull(Power)

  if_else(
    sft$powerEstimate < 16,
    max(sft$powerEstimate, my_estimate),
    my_estimate
  ) |>
    as.integer()
}

plot_adj_as_network <- function(matrix,
                                min_edge = 0.4,
                                node_size = 45,
                                node_fill_col = "yellow",
                                node_frame_col = NULL,
                                node_label_size = 1,
                                edge_size = 2,
                                layout = igraph::layout_in_circle,
                                ...) {
  set.seed(420)
  # get rid of low correlations
  matrix[matrix < min_edge] <- 0

  # build_network
  network <- igraph::graph_from_adjacency_matrix(
    matrix,
    mode = "upper",
    weighted = T,
    diag = F
  )

  igraph::V(network)$color <- node_fill_col

  igraph::V(network)$size <- node_size
  # igraph::V(network)$size <- igraph::degree(network, mode = "total")*node_size

  if (!is.null(node_label_size) & node_label_size > 0) {
    igraph::V(network)$label.color <- "black"
    igraph::V(network)$label.cex <- node_label_size
  }

  if (!is.null(node_frame_col)) {
    igraph::V(network)$frame.color <- node_frame_col
  }

  igraph::E(network)$width <- igraph::E(network)$weight*edge_size
  igraph::E(network)$color <- "black"

  writeLines("Visualizing a simplified representation of the circadian GCN")
  par(mfrow = c(1,1))
  plot(
    network,
    # size=100,
    # layout = igraph::layout.kamada.kawai,
    # layout = igraph::layout.fruchterman.reingold,
    # layout = igraph::layout.graphopt,
    layout = layout,
    # vertex.label=NA
    # vertex.shape="none"
    ...
  )

}

calculate_module_module_sim <- function(merged_modules, plot = TRUE) {

  cat("Calculating module-module similarity based
  on module-eigengene-expression...")
  # Calculate similarity of the eigen-genes
  sim_matrix_ME <- cor(
    merged_modules$newMEs,
    method = "kendall"
  )
  # calculate adj_matrix
  adj_matrix_ME <- WGCNA::adjacency.fromSimilarity(
    sim_matrix_ME,
    power=1,        # DO NOT power transform
    type='signed'
  ) |>
    as.matrix()
  cat("Done.\n")

  cat("Tidying module names...")
  ## CHANGE THE NAMES OF THE MODULES;
  module_ids <- data.frame(
    old_labels = rownames(adj_matrix_ME) %>%
      stringr::str_split("ME", 2) %>%
      sapply("[", 2) %>%
      as.character(),
    new_labels = paste0(
      "C",
      1:nrow(adj_matrix_ME)
    )
  )
  rownames(adj_matrix_ME) <- module_ids$new_labels
  colnames(adj_matrix_ME) <- module_ids$new_labels
  cat("Done.\n")

  if (plot) {
    cat("Plotting adjacency matrix for module-module similarity...\n")
    tryCatch({
      grDevices::dev.new()
      gplots::heatmap.2(
        t(adj_matrix_ME),
        col = viridis::inferno(100),
        trace = 'none',
        dendrogram = 'row',
        xlab = '',
        ylab = '',
        main = 'Adjacency matrix - MEs \n(module-module similarity)',
        density.info = 'none',
        revC = TRUE
      )
      trash <- grDevices::dev.off()
    }, error = function(e) {
      cat(
        "Plotting failed with error:", conditionMessage(e),
        "\nContinuing execution...\n"
      )
    })
  }

  list(
    "ME" = adj_matrix_ME,
    "mapping_tbl" = module_ids
  )

}
