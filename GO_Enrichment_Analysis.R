### GO Enrichment Analysis
### Performs comprehensive GO enrichment for UP and DOWN regulated proteins
### across BP, MF, and CC ontologies using weight01-fisher algorithm
### Generates hierarchy maps (igraph/ggraph) and exports results as CSV/PNG

library(topGO)
library(igraph)
library(ggraph)
library(GO.db)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)

# Named constants
MIN_PVALUE   <- 1e-100   # substitute for p-values reported as "<1e-XX" by topGO
MIN_EXPECTED <- 1e-6     # floor for Expected count to prevent division by zero

# ============================================================
# HELPER FUNCTIONS
# ============================================================

#' Return named list of parent GO IDs for a given ontology
#' @param ontology One of "BP", "MF", "CC"
get_parent_map <- function(ontology = c("BP", "MF", "CC")) {
  ontology <- match.arg(ontology)
  switch(ontology,
    BP = as.list(GOBPPARENTS),
    MF = as.list(GOMFPARENTS),
    CC = as.list(GOCCPARENTS)
  )
}

#' Return the root GO ID for a given ontology
#' @param ontology One of "BP", "MF", "CC"
get_root_id <- function(ontology = c("BP", "MF", "CC")) {
  ontology <- match.arg(ontology)
  switch(ontology,
    BP = "GO:0008150",
    MF = "GO:0003674",
    CC = "GO:0005575"
  )
}

#' Return all ancestor GO IDs for a set of GO IDs (BFS, bounded depth)
#' @param ids       Character vector of GO IDs (seeds)
#' @param parent_map Named list returned by get_parent_map()
#' @param max_depth  Maximum number of hops above the seeds to traverse.
#'   Default Inf (unlimited) keeps all existing call-sites compatible.
get_ancestors <- function(ids, parent_map, max_depth = Inf) {
  visited <- character(0)
  # Queue entries are two-element lists: list(id = <string>, depth = <int>)
  queue   <- lapply(ids, function(x) list(id = x, depth = 0L))

  while (length(queue) > 0) {
    entry   <- queue[[1]]
    queue   <- queue[-1]
    current <- entry$id
    depth   <- entry$depth

    if (current %in% visited) next
    visited <- c(visited, current)

    if (depth >= max_depth) next   # do not enqueue parents beyond the limit

    parents <- parent_map[[current]]
    parents <- parents[!is.na(parents) & parents != "all"]
    if (length(parents) > 0) {
      new_entries <- lapply(parents, function(p) list(id = p, depth = depth + 1L))
      queue <- c(queue, new_entries)
    }
  }

  setdiff(visited, ids)   # return ancestors only (not the seed ids)
}

#' Build an igraph subgraph for a set of GO IDs plus their ancestors
#' @param go_ids            Character vector of GO IDs to include
#' @param ontology          One of "BP", "MF", "CC"
#' @param selected_df       Optional data.frame with columns GO.ID, weight01Fisher,
#'   neg_log10_p, Annotated, Significant used to annotate nodes
#' @param max_ancestor_depth Maximum BFS hops above seed terms to traverse.
#'   Limits graph size by preventing deep traversal toward the ontology root.
build_go_subgraph <- function(go_ids, ontology = "BP", selected_df = NULL,
                               max_ancestor_depth = 4L) {
  parent_map <- get_parent_map(ontology)

  # Collect ancestors (bounded depth) so the hierarchy is connected
  ancestors  <- get_ancestors(go_ids, parent_map, max_depth = max_ancestor_depth)
  all_ids    <- unique(c(go_ids, ancestors))

  # Build edge list: parent -> child  (corrected direction)
  edges <- lapply(all_ids, function(id) {
    parents <- parent_map[[id]]
    parents <- parents[!is.na(parents) & parents != "all"]
    parents <- intersect(parents, all_ids)   # keep only nodes already in the set
    if (length(parents) == 0) return(NULL)
    data.frame(from = parents, to = id, stringsAsFactors = FALSE)
  })
  edges <- do.call(rbind, Filter(Negate(is.null), edges))

  if (is.null(edges) || nrow(edges) == 0) {
    g <- make_empty_graph(n = length(go_ids), directed = TRUE)
    V(g)$name        <- go_ids
    V(g)$selected    <- TRUE
    V(g)$neg_log10_p <- NA_real_
    V(g)$annotated   <- NA_integer_
    V(g)$significant <- NA_integer_
    V(g)$depth       <- 0L
    return(g)
  }

  # Pin vertex creation order via the `vertices` argument
  g <- graph_from_data_frame(
    edges,
    directed = TRUE,
    vertices = data.frame(name = all_ids, stringsAsFactors = FALSE)
  )

  # Node metadata — use match() to guard against row-reorder from select()
  go_terms_df <- tryCatch(
    AnnotationDbi::select(GO.db, keys = V(g)$name,
                          columns = "TERM", keytype = "GOID"),
    error = function(e) data.frame(GOID = V(g)$name, TERM = NA_character_,
                                   stringsAsFactors = FALSE)
  )
  term_idx      <- match(V(g)$name, go_terms_df$GOID)
  V(g)$term     <- go_terms_df$TERM[term_idx]
  V(g)$selected <- V(g)$name %in% go_ids

  # Enrichment data (only for seed/selected nodes)
  V(g)$neg_log10_p <- NA_real_
  V(g)$annotated   <- NA_integer_
  V(g)$significant <- NA_integer_

  if (!is.null(selected_df) && nrow(selected_df) > 0) {
    idx <- match(V(g)$name, selected_df$GO.ID)
    V(g)$neg_log10_p <- selected_df$neg_log10_p[idx]
    if ("Annotated"   %in% colnames(selected_df))
      V(g)$annotated   <- as.integer(selected_df$Annotated[idx])
    if ("Significant" %in% colnames(selected_df))
      V(g)$significant <- as.integer(selected_df$Significant[idx])
  }

  # Compute BFS depth from root(s) (in-degree == 0 after parent→child edges)
  roots <- which(degree(g, mode = "in") == 0)
  if (length(roots) == 0) roots <- 1L
  dist_mat  <- distances(g, v = roots, mode = "out")
  # For each node take the minimum distance from any root
  min_dists <- apply(dist_mat, 2, min)
  min_dists[!is.finite(min_dists)] <- 0L
  V(g)$depth <- as.integer(min_dists)

  g
}

#' ggraph hierarchy plot with depth-pinned y-axis and enrichment-aware encoding
#'
#' Produces a DAG with the GO root at the top (y = 0) and enriched leaf terms
#' at the bottom, matching the AEGIS / Reimand-2019 visual style:
#'   - Colour: continuous -log10(p) gradient (grey -> yellow -> orange -> red)
#'   - Size  : proportional to log1p(Annotated) for enriched nodes
#'   - Labels: top-N enriched terms (wrapped term name + p-value)
#'   - Arrows: downward, parent -> child
#'
#' @param g      igraph object returned by build_go_subgraph()
#' @param title  Plot title string
#' @param top_n  Number of top terms to label
plot_go_hierarchy <- function(g, title = "", top_n = 5) {
  if (vcount(g) == 0) {
    return(ggplot() +
             labs(title = title, subtitle = "No nodes to display") +
             theme_void())
  }

  # --- Node size: log1p(Annotated) for enriched; fixed small for ancestors ---
  node_sizes <- ifelse(
    V(g)$selected & !is.na(V(g)$annotated),
    log1p(V(g)$annotated) + 2,
    2
  )
  V(g)$node_size <- node_sizes

  # --- Identify top-N nodes to label (by neg_log10_p) ---
  top_ids <- character(0)
  if (any(!is.na(V(g)$neg_log10_p))) {
    top_ids <- V(g)$name[
      order(V(g)$neg_log10_p, decreasing = TRUE, na.last = TRUE)
    ][seq_len(min(top_n, sum(!is.na(V(g)$neg_log10_p))))]
  }

  # Safe loop-based label assignment (avoids ifelse NA-propagation)
  labels <- character(vcount(g))
  for (i in seq_len(vcount(g))) {
    if (V(g)$name[i] %in% top_ids) {
      term_txt  <- if (is.na(V(g)$term[i])) "" else str_wrap(V(g)$term[i], 22)
      p_val     <- V(g)$neg_log10_p[i]
      p_txt     <- if (is.na(p_val)) "" else sprintf("(p=%.1e)", 10^(-p_val))
      labels[i] <- paste0(term_txt, "\n", p_txt)
    }
  }
  V(g)$label <- labels

  # --- Sugiyama layout for x-ordering; override y with -depth ---
  lay <- tryCatch(
    create_layout(g, layout = "sugiyama"),
    error = function(e) create_layout(g, layout = "kk")
  )
  if ("depth" %in% vertex_attr_names(g)) {
    lay$y <- -V(g)$depth   # root at y = 0; deeper nodes at more negative y
  }

  ggraph(lay) +
    geom_edge_link(
      colour  = "grey55",
      arrow   = arrow(length = unit(2.5, "mm"), type = "closed"),
      end_cap = circle(4, "mm"),
      alpha   = 0.55
    ) +
    geom_node_point(
      aes(fill = neg_log10_p, size = node_size),
      shape  = 21,
      colour = "white",
      alpha  = 0.9
    ) +
    geom_node_text(
      aes(label = label),
      size         = 2.4,
      repel        = TRUE,
      max.overlaps = 25,
      colour       = "black",
      lineheight   = 0.9
    ) +
    scale_fill_gradientn(
      colours  = c("#CCCCCC", "#FFFFB2", "#FD8D3C", "#BD0026"),
      na.value = "#D9D9D9",
      name     = expression(-log[10](p)),
      guide    = guide_colourbar(barwidth = 8, barheight = 0.5,
                                 title.position = "top")
    ) +
    scale_size_continuous(range = c(2, 10), guide = "none") +
    labs(
      title    = title,
      subtitle = paste0(
        "Top ", top_n, " terms labelled  |  ",
        "Colour: -log\u2081\u2080(p)  |  Size \u221d log(Annotated + 1)"
      )
    ) +
    theme_void() +
    theme(
      plot.title       = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle    = element_text(size = 9,  hjust = 0.5, colour = "grey40"),
      legend.position  = "bottom",
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8)
    )
}

#' Build a gene -> GO mapping list for use as topGO's gene2GO argument
#' @param symbols Character vector of gene symbols
#' @param ontology One of "BP", "MF", "CC"
build_gene2go <- function(symbols, ontology = "BP") {
  # Map gene symbols -> Entrez IDs
  eg <- suppressMessages(
    AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys    = symbols,
      column  = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
  )
  eg <- eg[!is.na(eg)]

  # Retrieve GO terms for each gene via the gene-centric map
  go_map <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys    = unique(unname(eg)),
      columns = c("SYMBOL", "GO", "ONTOLOGY"),
      keytype = "ENTREZID"
    )
  ) %>%
    filter(ONTOLOGY == ontology) %>%
    group_by(SYMBOL) %>%
    summarise(go_ids = list(unique(GO)), .groups = "drop")

  gene2go <- setNames(go_map$go_ids, go_map$SYMBOL)

  # Ensure every input symbol is represented (even if empty)
  missing_syms <- setdiff(symbols, names(gene2go))
  if (length(missing_syms) > 0) {
    gene2go[missing_syms] <- lapply(missing_syms, function(x) character(0))
  }

  gene2go
}

# ============================================================
# UTILITY: parse GenTable into a clean tibble
# ============================================================

.make_table <- function(go_data_obj, result_obj, top_n = 25) {
  GenTable(go_data_obj, weight01Fisher = result_obj, topNodes = top_n) %>%
    as_tibble() %>%
    mutate(
      weight01Fisher = sapply(weight01Fisher, function(x) {
        if (is.character(x) && grepl("<", x)) return(MIN_PVALUE)
        suppressWarnings(as.numeric(x))
      }),
      Annotated      = as.numeric(Annotated),
      Significant    = as.numeric(Significant),
      Expected       = as.numeric(Expected),
      fold_enrichment = Significant / pmax(Expected, MIN_EXPECTED),
      neg_log10_p    = -log10(weight01Fisher)
    ) %>%
    filter(weight01Fisher < 0.05) %>%
    arrange(weight01Fisher)
}

# ============================================================
# MAIN ANALYSIS LOOP
# ============================================================

#' Run comprehensive GO enrichment for UP- and DOWN-regulated genes
#'
#' The function iterates over all requested ontologies and both expression
#' directions.  For each combination it:
#'   - creates a topGO result object
#'   - builds a GenTable tibble
#'   - builds a GO hierarchy igraph object
#'   - produces a ggraph hierarchy plot
#'   - exports a CSV and a PNG to output_dir
#'
#' Named global variables are created for interactive use (matching the
#' conventions required by downstream code):
#'   res_weight01_{ONT}_{DIR}, df_go_{ONT}_{DIR}_table,
#'   g_{ONT}_{DIR}, p_{ONT}_{DIR}
#'
#' The function invisibly returns a nested list with the same data so that
#' callers can use the results programmatically without relying on global state.
#'
#' @param df_go       tibble with columns SYMBOL, log2FC, adj.pvalue
#' @param output_dir  directory for CSV/PNG output (created if absent)
#' @param ontologies  character vector; subset of c("BP","MF","CC")
#' @param lfc_cutoff  |log2FC| threshold for gene selection (default 1)
#' @param top_nodes   max GO terms per direction for visualisation (default 25)
run_go_enrichment <- function(df_go,
                               output_dir = ".",
                               ontologies  = c("BP", "MF", "CC"),
                               lfc_cutoff  = 1,
                               top_nodes   = 25) {

  stopifnot(
    is.data.frame(df_go),
    all(c("SYMBOL", "log2FC", "adj.pvalue") %in% colnames(df_go))
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Gene universe: all measured symbols
  gene_universe <- unique(df_go$SYMBOL)

  # --- Significant gene sets by direction
  sig_up   <- df_go %>% filter(log2FC >=  lfc_cutoff) %>% pull(SYMBOL) %>% unique()
  sig_down <- df_go %>% filter(log2FC <= -lfc_cutoff) %>% pull(SYMBOL) %>% unique()

  cat(sprintf("Gene universe : %d genes\n", length(gene_universe)))
  cat(sprintf("Sig UP genes  : %d (log2FC >= %g)\n", length(sig_up),   lfc_cutoff))
  cat(sprintf("Sig DOWN genes: %d (log2FC <= -%g)\n\n", length(sig_down), lfc_cutoff))

  results <- list()

  for (ont in ontologies) {
    cat(sprintf("=== Ontology: %s ===\n", ont))

    # Build gene2GO mapping once per ontology
    gene2go <- build_gene2go(gene_universe, ontology = ont)

    for (direction in c("UP", "DOWN")) {
      sig_genes <- if (direction == "UP") sig_up else sig_down

      if (length(sig_genes) < 5) {
        cat(sprintf("  [%s %s] Skipping — fewer than 5 significant genes.\n",
                    ont, direction))
        next
      }

      cat(sprintf("  Running topGO for %s %s ...\n", ont, direction))

      # Build topGOdata and run test
      all_genes_fac <- factor(as.integer(gene_universe %in% sig_genes))
      names(all_genes_fac) <- gene_universe

      go_data <- tryCatch(
        new("topGOdata",
            ontology = ont,
            allGenes = all_genes_fac,
            geneSel  = function(x) x == 1,
            annot    = annFUN.gene2GO,
            gene2GO  = gene2go),
        error = function(e) {
          cat(sprintf("  [%s %s] topGOdata error: %s\n", ont, direction, e$message))
          NULL
        }
      )
      if (is.null(go_data)) next

      res_obj <- tryCatch(
        runTest(go_data, algorithm = "weight01", statistic = "fisher"),
        error = function(e) {
          cat(sprintf("  [%s %s] runTest error: %s\n", ont, direction, e$message))
          NULL
        }
      )
      if (is.null(res_obj)) next

      # Result table
      tbl <- tryCatch(
        .make_table(go_data, res_obj, top_n = top_nodes),
        error = function(e) {
          cat(sprintf("  [%s %s] GenTable error: %s\n", ont, direction, e$message))
          tibble()
        }
      )

      cat(sprintf("  [%s %s] %d significant GO terms\n",
                  ont, direction, nrow(tbl)))

      # --- Consistent variable names expected by the problem statement ---
      # res_weight01_{ONTOLOGY}_{DIRECTION}
      # df_go_{ONTOLOGY}_{DIRECTION}_table
      res_name <- sprintf("res_weight01_%s_%s",  ont, direction)
      tbl_name <- sprintf("df_go_%s_%s_table",   ont, direction)
      assign(res_name, res_obj,  envir = .GlobalEnv)
      assign(tbl_name, tbl,      envir = .GlobalEnv)

      # Export CSV
      csv_path <- file.path(output_dir,
                            sprintf("GO_%s_%s_enrichment.csv", ont, direction))
      write.csv(tbl, csv_path, row.names = FALSE)
      cat(sprintf("  [%s %s] CSV: %s\n", ont, direction, csv_path))

      # --- Build hierarchy graph ---
      if (nrow(tbl) == 0) next

      g <- tryCatch(
        build_go_subgraph(tbl$GO.ID, ontology = ont, selected_df = tbl),
        error = function(e) {
          cat(sprintf("  [%s %s] build_go_subgraph error: %s\n",
                      ont, direction, e$message))
          NULL
        }
      )
      if (is.null(g)) next

      # g_{ONTOLOGY}_{DIRECTION}
      g_name <- sprintf("g_%s_%s", ont, direction)
      assign(g_name, g, envir = .GlobalEnv)

      # Plot
      plot_title <- sprintf("GO Hierarchy - %s %s-regulated\n(%s)", ont, direction,
                            ifelse(direction == "UP",
                                   paste0("log2FC >= ", lfc_cutoff),
                                   paste0("log2FC <= -", lfc_cutoff)))

      p <- tryCatch(
        plot_go_hierarchy(g, title = plot_title, top_n = 5),
        error = function(e) {
          cat(sprintf("  [%s %s] plot error: %s\n", ont, direction, e$message))
          NULL
        }
      )
      if (is.null(p)) next

      # p_{ONTOLOGY}_{DIRECTION}
      p_name <- sprintf("p_%s_%s", ont, direction)
      assign(p_name, p, envir = .GlobalEnv)

      png_path <- file.path(output_dir,
                            sprintf("GO_%s_%s_hierarchy.png", ont, direction))
      tryCatch(
        ggsave(png_path, p, width = 14, height = 12, dpi = 300),
        error = function(e)
          cat(sprintf("  [%s %s] ggsave error: %s\n", ont, direction, e$message))
      )
      cat(sprintf("  [%s %s] PNG: %s\n", ont, direction, png_path))

      # Collect into return list
      results[[sprintf("%s_%s", ont, direction)]] <- list(
        result     = res_obj,
        table      = tbl,
        graph      = g,
        plot       = p
      )
    }
    cat("\n")
  }

  invisible(results)
}

# ============================================================
# EXAMPLE USAGE
# ============================================================
# The code below is commented out so that the script can be
# sourced/imported without running against actual data.
# Uncomment and adapt paths/variable names to run the analysis.

# # Load your data
# # df_go should be a tibble with columns: SYMBOL, log2FC, adj.pvalue
# # e.g.:
# # df_go <- read_csv("my_protein_data.csv")
#
# results <- run_go_enrichment(
#   df_go      = df_go,
#   output_dir = "GO_Results",
#   ontologies = c("BP", "MF", "CC"),
#   lfc_cutoff = 1,
#   top_nodes  = 25
# )
#
# # After running, the following objects are available globally:
# # res_weight01_BP_UP, res_weight01_BP_DOWN
# # res_weight01_MF_UP, res_weight01_MF_DOWN
# # res_weight01_CC_UP, res_weight01_CC_DOWN
# # df_go_BP_UP_table,  df_go_BP_DOWN_table
# # df_go_MF_UP_table,  df_go_MF_DOWN_table
# # df_go_CC_UP_table,  df_go_CC_DOWN_table
# # g_BP_UP,  g_BP_DOWN
# # g_MF_UP,  g_MF_DOWN
# # g_CC_UP,  g_CC_DOWN
# # p_BP_UP,  p_BP_DOWN
# # p_MF_UP,  p_MF_DOWN
# # p_CC_UP,  p_CC_DOWN
#
# # Print a summary of BP UP results:
# print(df_go_BP_UP_table)
