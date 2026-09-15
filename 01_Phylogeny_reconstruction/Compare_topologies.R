# Comparison_treefiles_topology_option.R
library(phytools)
library(ape)
library(viridis)

# 定义需要用黄色显示的物种名称（示例）
yellow_species <- c("1", "2")

get_tip_colors <- function(tree, species_list = yellow_species, highlight_color = "gold") {
  tip_colors <- rep("black", length(tree$tip.label))
  for (species in species_list) {
    matching_tips <- grep(species, tree$tip.label, value = FALSE)
    if (length(matching_tips) > 0) tip_colors[matching_tips] <- highlight_color
  }
  return(tip_colors)
}

# ---------------- Parameters ----------------
comparison_mode <- "group"   # "single" or "group"
work_dir <- "./"
ref_tree <- "./01_Phylogeny_reconstruction/ASTRAL4_HybSuite_HRS.bootstrap.rr.tre"
# right_tree <- "./01_Phylogeny_reconstruction/Final_species_tree_ASTRAL-IV.tre"
output_file <- "Final_ASTRAL-IV_comparison"
output_type <- "pdf"         # "png" or "pdf"
line_alpha <- 0.9           # numeric 0..1
single_title <- "Summary tree: MI+MO vs original 1to1"
outgroup_tips <- c("Sibiraea_angustata", "Petrophytum_caespitosum")  # vector or NULL
# Arabidopsis100: c("Sibiraea_angustata", "Petrophytum_caespitosum") 
# Angiosperms353: "Sibiraea_angustata"
cophylo_line_color <- "clade" # "single" or "clade"
show_branch_length <- FALSE  # FALSE => remove edge.length (topology only)
remove_outgroups <- FALSE    # logical: if true, will remove outgroup_tips in process
# Clades.tsv can be downloaded here:
# Angiosperms353:
# Arabidopsis100: 
clade_table <- "Clade.tsv" 

color_branches <- TRUE     # logical: if true, color tree branches according to clade.tsv colors
# --------------------------------------------------

if (comparison_mode == "group") {
  HRS_tree <- sub("_[^_]+\\.tre$", "_HRS.bootstrap.rr.tre", ref_tree)
  RLWP_tree <- sub("_[^_]+\\.tre$", "_RLWP.bootstrap.rr.tre", ref_tree)
  LS_tree <- sub("_[^_]+\\.tre$", "_LS.bootstrap.rr.tre", ref_tree)
  MI_tree <- sub("_[^_]+\\.tre$", "_MI.bootstrap.rr.tre", ref_tree)
  MO_tree <- sub("_[^_]+\\.tre$", "_MO.bootstrap.rr.tre", ref_tree)
  RT_tree <- sub("_[^_]+\\.tre$", "_RT.bootstrap.rr.tre", ref_tree)
  one_to_one_tree <- sub("_[^_]+\\.tre$", "_1to1.bootstrap.rr.tre", ref_tree)
  message("Derived group tree names:")
  print(list(HRS=HRS_tree, RLWP=RLWP_tree, LS=LS_tree, MI=MI_tree, MO=MO_tree, RT=RT_tree, one2one=one_to_one_tree))
}

setwd(work_dir)

# ----------------- process_tree (合并、安全) -----------------
process_tree <- function(tree_file,
                         outgroup_tips = NULL,
                         scale_factor = 0.5,
                         root_scale = 0.1,
                         ladderize_right = FALSE,
                         remove_outgroups = FALSE) {
  if (!file.exists(tree_file)) stop("tree_file not found: ", tree_file)
  tree <- read.tree(tree_file)
  if (remove_outgroups && !is.null(outgroup_tips) && length(outgroup_tips) > 0) {
    outgroup_in_tree <- intersect(tree$tip.label, outgroup_tips)
    if (length(outgroup_in_tree) > 0) tree <- drop.tip(tree, outgroup_in_tree)
  }
  tree <- ladderize(tree, right = ladderize_right)
  if (!is.null(tree$edge.length) && length(tree$edge.length) > 0) {
    bad_idx <- which(is.nan(tree$edge.length) | is.infinite(tree$edge.length))
    if (length(bad_idx) > 0) tree$edge.length[bad_idx] <- 1e-6
    tree$edge.length <- tree$edge.length * scale_factor
    root_edge <- which(tree$edge[,1] == Ntip(tree) + 1)
    if (length(root_edge) > 0) {
      tree$edge.length[root_edge] <- tree$edge.length[root_edge] * root_scale
    }
  }
  return(tree)
}

process_RT_tree <- function(tree_file, ...) {
  tree <- process_tree(tree_file, ..., remove_outgroups = FALSE)
  if (!is.null(tree$edge.length) && length(tree$edge.length) > 0) {
    tree$edge.length_original <- tree$edge.length
    tree$edge.length[is.nan(tree$edge.length)] <- 1e-6
    tree$edge.length[tree$edge.length == 0] <- 1e-6
  }
  return(tree)
}

# ----------------- make_edge_colors_from_clade -------------------
make_edge_colors_from_clade <- function(cophylo_obj,
                                        mapping_df,
                                        color_map,
                                        species_col = "Species",
                                        clade_col = "Clade",
                                        default_color = "#000000") {

  if (is.character(mapping_df)) {
    mapping_df <- read.table(mapping_df, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE)
  }

  species2clade <- setNames(
    as.character(mapping_df[[clade_col]]),
    mapping_df[[species_col]]
  )

  get_tree_edge_colors <- function(tree) {

    tip_clade <- species2clade[tree$tip.label]

    tip_colors <- color_map[tip_clade]

    tip_colors[is.na(tip_colors)] <- default_color

    edge_colors <- rep(default_color, nrow(tree$edge))

    for (i in seq_len(nrow(tree$edge))) {

      node <- tree$edge[i, 2]

      desc <- getDescendants(tree, node)

      tips <- desc[desc <= Ntip(tree)]

      if (length(tips) > 0) {

        clades <- tip_clade[tips]

        clades <- clades[!is.na(clades)]

        if (length(unique(clades)) == 1) {

          edge_colors[i] <- color_map[unique(clades)]

        } else {

          edge_colors[i] <- default_color

        }
      }
    }

    return(edge_colors)
  }

  left_tree  <- cophylo_obj$trees[[1]]
  right_tree <- cophylo_obj$trees[[2]]

  left_col  <- get_tree_edge_colors(left_tree)
  right_col <- get_tree_edge_colors(right_tree)

  return(list(left_col, right_col))
}

# ----------------- make_link_colors (校验 alpha) -----------------
make_link_colors <- function(cophylo_obj,
                             mapping_df,
                             color_map = NULL,
                             line_alpha = 0.7,
                             species_col = "Species",
                             clade_col = "Clade",
                             default_color = "#000000") {
  line_alpha_num <- as.numeric(line_alpha)
  if (is.na(line_alpha_num) || line_alpha_num < 0 || line_alpha_num > 1) {
    warning("line_alpha must be numeric between 0 and 1. Using 0.7 fallback.")
    line_alpha_num <- 0.7
  }
  if (is.character(mapping_df) && length(mapping_df) == 1 && file.exists(mapping_df)) {
    map_df <- read.table(mapping_df, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  } else if (is.data.frame(mapping_df)) {
    map_df <- mapping_df
  } else stop("mapping_df must be a data.frame or a valid filename")
  if (!all(c(species_col, clade_col) %in% colnames(map_df))) stop("mapping_df must contain columns: ", species_col, " and ", clade_col)

  okabe_ito_colors_dark <- c(
    "black"         = "#000000",
    "orange"        = "#C87500",
    "skyblue"       = "#2A7FBF",
    "bluishgreen"   = "#007B55",
    "yellow"        = "#C9B500",
    "blue"          = "#005EA2",
    "vermillion"    = "#B34100",
    "reddishpurple" = "#A23C80"
  )
  if (is.null(color_map)) {
    unique_clades <- sort(unique(map_df[[clade_col]]))
    palette_vals <- rep(unname(okabe_ito_colors_dark), length.out = length(unique_clades))
    color_map <- setNames(palette_vals, unique_clades)
  } else {
    if (!is.character(color_map) || is.null(names(color_map))) stop("color_map must be a named character vector where names are clade labels")
  }

  species2clade <- setNames(as.character(map_df[[clade_col]]), map_df[[species_col]])
  assoc <- cophylo_obj$assoc
  if (!is.null(assoc)) {
    if (is.numeric(assoc[,1])) {
      left_labels <- cophylo_obj$trees[[1]]$tip.label[assoc[,1]]
    } else left_labels <- as.character(assoc[,1])
  } else {
    left_labels <- cophylo_obj$trees[[1]]$tip.label
    warning("cophylo_obj$assoc not found; using left tree tip.label order for mapping.")
  }

  link_clades <- species2clade[left_labels]
  link_colors <- color_map[link_clades]
  bad <- which(is.na(link_colors) | link_colors == "")
  if (length(bad) > 0) {
    link_colors[bad] <- default_color
    warning(sprintf("Replaced %d missing/invalid colors with default_color '%s'. Example species: %s",
                    length(bad), default_color, paste(head(left_labels[bad], 5), collapse = ", ")))
  }

  safe_make_transparent <- function(cols, alpha) {
    cols2 <- as.character(cols)
    cols2[is.na(cols2) | cols2 == ""] <- default_color
    valid <- sapply(cols2, function(cl) {
      tryCatch({ grDevices::col2rgb(cl); TRUE }, error = function(e) FALSE)
    })
    cols2[!valid] <- default_color
    sapply(cols2, function(cl) make.transparent(cl, alpha))
  }

  final_colors <- safe_make_transparent(link_colors, line_alpha_num)
  names(final_colors) <- NULL
  return(final_colors)
}

# ------------- 新增函数：从 cophylo 对象中移除 edge.length -------------
remove_edge_lengths_from_cophylo <- function(cophylo_obj) {
  if (is.null(cophylo_obj) || !is.list(cophylo_obj) || is.null(cophylo_obj$trees)) return(cophylo_obj)
  for (i in seq_along(cophylo_obj$trees)) {
    tr <- cophylo_obj$trees[[i]]
    if (!is.null(tr$edge.length)) {
      tr$edge.length <- NULL
      # 另外也清空 cached x/y coordinates 以避免旧坐标基于枝长（若存在）
      tr$node.label <- tr$node.label  # 不变，但保留
    }
    cophylo_obj$trees[[i]] <- tr
  }
  return(cophylo_obj)
}
# -------------------------------------------------------------------

# 读取 clade 表
clade_df <- read.table(clade_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 预设 clade->color
okabe_ito_colors_dark <- c(
  "black"        = "#000000",
  "orange"       = "#1f78b4",
  "skyblue"      = "#8fce00",
  "bluishgreen"  = "#ffb000",
  "yellow"       = "#d18def",
  "reddishpurple"= "#70328b"
)
#  "#70328b",  # 顶生复伞房花序  明亮蓝
#  "#d18def",  # 侧生复伞房花序  青蓝（明显区分）
#  "#ffb000",  # 圆锥花序        金
#  "#8fce00",  # 伞房花序        绿
#  "#1f78b4",  # 伞形花序        蓝
#  "#c0bebe"
clade_color_map <- c(
  "1" = unname("#1f78b4"),
  "2" = unname("#8fce00"),
  "3" = unname("#70328b"),
  "4" = unname("#ffb000"),
  "5" = unname("#d18def"),
  "Outgroup" = unname("black")
)

# ---------------- main: single ----------------
if (comparison_mode == "single") {
  tree1 <- process_tree(ref_tree, outgroup_tips = outgroup_tips, remove_outgroups = remove_outgroups)
  if (grepl("RT", right_tree, ignore.case = FALSE)) {
    tree2 <- process_RT_tree(right_tree, outgroup_tips = outgroup_tips, remove_outgroups = remove_outgroups)
  } else {
    tree2 <- process_tree(right_tree, outgroup_tips = outgroup_tips, remove_outgroups = remove_outgroups)
  }

  cophylo_obj <- cophylo(tree1, tree2, rotate = TRUE)

  # 如果不显示枝长，就把 cophylo 对象里的边长移除（彻底删除）
  if (!show_branch_length) {
    cophylo_obj <- remove_edge_lengths_from_cophylo(cophylo_obj)
  }

  # open device once
  if (output_type == "png") {
    png(paste0(output_file, ".png"), width = 6600, height = 4500, res = 300)
  } else if (output_type == "pdf") {
    pdf(paste0(output_file, ".pdf"), width = 20, height = 18)
  }

  par(mar = c(0,0,0,0))
  use_edge_flag <- show_branch_length

  if (cophylo_line_color == "single") {
    plot.cophylo(cophylo_obj,
                 link.type = "curved", fsize = 1.2,
                 link.lwd = 3, tip.lty = 5, tip.lwd = 0.8, lwd = 2, link.lty = "solid",
                 tip.len = 0, pts = TRUE,
                 link.col = make.transparent("tomato", line_alpha),
                 use.edge.length = use_edge_flag)
  } else if (cophylo_line_color == "clade") {
    link_colors <- make_link_colors(cophylo_obj = cophylo_obj,
                                    mapping_df = clade_df,
                                    color_map = clade_color_map,
                                    line_alpha = line_alpha,
                                    species_col = "Species",
                                    clade_col = "Clade",
                                    default_color = "#000000")
    edge_colors <- make_edge_colors_from_clade(cophylo_obj = cophylo_obj,
                                               mapping_df = clade_df,
                                               color_map = clade_color_map,
                                               species_col = "Species",
                                               clade_col = "Clade",
                                               default_color = "#000000"
                                              )
    plot.cophylo(cophylo_obj,
                 link.type = "curved", fsize = 1.2,
                 link.lwd = 3, tip.lty = 5, tip.lwd = 0.8, lwd = 2, link.lty = "solid",
                 tip.len = 0, pts = TRUE,
                 link.col = link_colors,
                 use.edge.length = use_edge_flag)
  }

  # Format node labels to show only 2 decimal places
  left_node_labels <- cophylo_obj$trees[[1]]$node.label
  right_node_labels <- cophylo_obj$trees[[2]]$node.label
  
  # Convert to numeric and format to 2 decimal places, keeping non-numeric values as is
  format_node_labels <- function(labels) {
    if (is.null(labels)) return(NULL)
    sapply(labels, function(label) {
      if (is.na(label) || label == "") return(label)
      # Try to convert to numeric
      num_label <- as.numeric(label)
      if (!is.na(num_label)) {
        return(sprintf("%.2f", num_label))
      } else {
        return(label)
      }
    })
  }
  
  left_node_labels_formatted <- format_node_labels(left_node_labels)
  right_node_labels_formatted <- format_node_labels(right_node_labels)
  
  nodelabels.cophylo(text = left_node_labels_formatted, frame = "none",
                     adj = c(1.3, -0.4), cex = 1.0, which = "left", col = "black")
  nodelabels.cophylo(text = right_node_labels_formatted, frame = "none",
                     adj = c(-0.3, -0.4), cex = 1.0, which = "right", col = "black")
  title(main = single_title, cex.main = 2.0,
        col.main = "black", font.main = 2, line = -1.5)

  dev.off()
}

# ---------------- main: group ----------------
if (comparison_mode == "group") {
  ref_tree_processed <- process_tree(ref_tree, outgroup_tips = outgroup_tips, remove_outgroups = remove_outgroups)

  cophylo_list <- list()
  titles <- c()

  add_if_exists <- function(label, treefile) {
    if (!is.null(treefile) && treefile != ref_tree && file.exists(treefile)) {
      tproc <- if (grepl("RT", treefile, ignore.case = FALSE)) process_RT_tree(treefile) else process_tree(treefile)
      cophylo_list[[label]] <<- cophylo(ref_tree_processed, tproc, rotate = TRUE)
      titles <<- c(titles, label)
    }
  }

  add_if_exists("RLWP", RLWP_tree)
  add_if_exists("LS", LS_tree)
  add_if_exists("MI", MI_tree)
  add_if_exists("MO", MO_tree)
  add_if_exists("RT", RT_tree)
  add_if_exists("1to1", one_to_one_tree)

  n_plots <- length(cophylo_list)
  if (n_plots > 0) {
    if (n_plots <= 2) {
      nrow <- 1; ncol <- n_plots
    } else if (n_plots <= 4) {
      nrow <- 2; ncol <- 2
    } else {
      nrow <- 2; ncol <- 3
    }

    if (output_type == "png") {
      png(paste0(output_file, "_group.png"), width = 6600 * ncol / 2, height = 4800 * nrow / 2, res = 300)
    } else if (output_type == "pdf") {
      pdf(paste0(output_file, "_group.pdf"), width = 22 * ncol / 2, height = 20 * nrow / 2)
    }

    par(mfrow = c(nrow, ncol), mar = c(0,20,5,10))
    use_edge_flag <- show_branch_length

    for (i in seq_len(n_plots)) {
      cophylo_obj <- cophylo_list[[i]]
      # 如果不显示枝长，就移除 cophylo 对象内的边长
      if (!show_branch_length) cophylo_obj <- remove_edge_lengths_from_cophylo(cophylo_obj)

      plot_title <- paste("HRS vs", titles[i])
      if (cophylo_line_color == "single") {
        plot.cophylo(cophylo_obj,
                     link.type = "curved", fsize = 1.2,
                     link.lwd = 3, tip.lty = 5, tip.lwd = 0.8, lwd = 2, link.lty = "solid",
                     tip.len = 0, pts = TRUE,
                     link.col = make.transparent("tomato", line_alpha),
                     use.edge.length = use_edge_flag)
      } else if (cophylo_line_color == "clade") {
        link_colors <- make_link_colors(cophylo_obj = cophylo_obj,
                                        mapping_df = clade_df,
                                        color_map = clade_color_map,
                                        line_alpha = line_alpha,
                                        species_col = "Species",
                                        clade_col = "Clade",
                                        default_color = "#000000")
        if ( color_branches == TRUE ) {
          cat("Plotting color branches ...")
          edge_colors <- make_edge_colors_from_clade(cophylo_obj = cophylo_obj,
                                               mapping_df = clade_df,
                                               color_map = clade_color_map,
                                               species_col = "Species",
                                               clade_col = "Clade",
                                               default_color = "#000000"
                                              )
          plot.cophylo(cophylo_obj,
                 link.type = "curved", fsize = 1.2,
                 link.lwd = 2.5, tip.lty = 3, tip.lwd = 1.5, lwd = 4.5, link.lty = "solid",
                 tip.len = 0, pts = FALSE,
                 link.col = link_colors,
                 use.edge.length = use_edge_flag,
                 edge.col = list(left = edge_colors[[1]], right = edge_colors[[2]]))
        } else {
          plot.cophylo(cophylo_obj,
                 link.type = "curved", fsize = 1.5,
                 link.lwd = 3, tip.lty = 4, tip.lwd = 0.8, lwd = 2, link.lty = "solid",
                 tip.len = 0, pts = FALSE,
                 link.col = link_colors,
                 use.edge.length = use_edge_flag,
                 left.edge.col  = edge_colors[[1]],
                 right.edge.col = edge_colors[[2]])
        }
      }    

      # Format node labels to show only 2 decimal places
      left_node_labels <- cophylo_obj$trees[[1]]$node.label
      right_node_labels <- cophylo_obj$trees[[2]]$node.label
      
      # Convert to numeric and format to 2 decimal places, keeping non-numeric values as is
      format_node_labels <- function(labels) {
        if (is.null(labels)) return(NULL)
        sapply(labels, function(label) {
          if (is.na(label) || label == "") return(label)
          # Try to convert to numeric
          num_label <- as.numeric(label)
          if (!is.na(num_label)) {
            return(sprintf("%.2f", num_label))
          } else {
            return(label)
          }
        })
      }
      
      left_node_labels_formatted <- format_node_labels(left_node_labels)
      right_node_labels_formatted <- format_node_labels(right_node_labels)
      
      nodelabels.cophylo(text = left_node_labels_formatted, frame = "none",
                        adj = c(1.3, -0.4), cex = 1.0, which = "left", col = "black")
      nodelabels.cophylo(text = right_node_labels_formatted, frame = "none",
                        adj = c(-0.3, -0.4), cex = 1.0, which = "right", col = "black")
      title(main = plot_title, cex.main = 2.4, col.main = "black", font.main = 2, line = -1.5)
    }

    dev.off()
    cat("成功绘制了", n_plots, "个比较图到文件:", paste0(output_file, "_group.", output_type), "\n")
    cat("包含的比较:", paste(titles, collapse = ", "), "\n")
  } else {
    cat("没有找到有效的树文件进行比较\n")
  }
}
