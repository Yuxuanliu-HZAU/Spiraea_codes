# Compare topologies

## 01 Configure your parameters
Change the parameters on the beginning of the script [`Compare_topologies.R`]().
```
# ---------------- Parameters ----------------
comparison_mode <- "group"   # "single" or "group"
work_dir <- "./"
ref_tree <- "./01_Phylogeny_reconstruction/ASTRAL4_HybSuite_HRS.bootstrap.rr.tre"
output_file <- "Final_ASTRAL-IV_comparison"
output_type <- "pdf"  # output figures in "png" or "pdf" format
line_alpha <- 0.9     # numeric 0..1
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
```

## 02 Run `Compare_topologies.R`

```
Rscript Compare_topologies.R
```

After running this script, a figure showing the topology camparison among seven nuclear trees generated with different paralog-handling methods will be outputted.

In this study, we campare the topologies of seven nuclear trees genrated via HRS, RLWP, LS, MI, MO, RT, and 1to1, using HRS topology as reference.

The figure is like:
