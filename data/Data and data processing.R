# DATA Processing
# data collected from CMap -> LINCS L1000
# no data pre-processing if the data used is level 5


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("cmapR")


# Load the compound metadata file
compound_metadata <- read.delim("compoundinfo_beta.txt", header = TRUE, sep = "\t")

# Check the first few rows to understand the structure
head(compound_metadata)

# Filter compounds for any drug related to AML (using 'AML' as a keyword)
aml_compounds <- subset(compound_metadata, grepl("AML", cmap_name, ignore.case=TRUE) | 
                          grepl("AML", compound_aliases, ignore.case=TRUE))

# View the filtered compounds
head(aml_compounds)

# Filter for drugs with a specific mechanism of action (e.g., kinase inhibitors)
kinase_inhibitors <- subset(compound_metadata, grepl("kinase", moa, ignore.case=TRUE))

# View the filtered drugs
head(kinase_inhibitors)

# Extract the pert_id for the filtered AML-related compounds
subset_cid <- aml_compounds$pert_id

# View the selected perturbagen IDs
head(subset_cid)

# Extract the column IDs (cids) from the GCTX file
gctx_treatment <- "level5_beta_trt_cp_n720216x12328.gctx"
all_cid <- parse_gctx_ids(gctx_treatment, dim = "col")

# View the first few column IDs
head(all_cid)

library(cmapR)
ls("package:cmapR")

# Load only the column metadata (pert_id)
meta <- read_gctx_meta(gctx_treatment, dim = "col")

# Extract all column IDs (pert_id) from the metadata
all_cid <- meta$id
head(all_cid)

# Find matching cids between your subset and the GCTX file
matching_cid <- intersect(subset_cid, all_cid)
# ---- output ---- #
# [1] "ABY001_A375_XH:BRD-A61304759:0.625:24"
# [2] "ABY001_A375_XH:BRD-A61304759:0.625:3" 
# [3] "ABY001_A375_XH:BRD-A61304759:10:24"   
# [4] "ABY001_A375_XH:BRD-A61304759:10:3"    
# [5] "ABY001_A375_XH:BRD-A61304759:2.5:24"  
# [6] "ABY001_A375_XH:BRD-A61304759:2.5:3" 
# ---------------- #

# Check the matched cids
head(matching_cid)
# ---- output ---- #
# character(0)
# ---------------- #
# figure out what happened here





# Read the treatment file (drug-treated samples)
#gctx_treatment <- "level5_beta_trt_cp_n720216x12328.gctx"
#ds_treatment <- parse_gctx(gctx_treatment)

# Explore the data structure
#head(ds_treatment@mat)  # View the gene expression matrix

# Read the control file (untreated samples)
#gctx_control <- "level5_beta_ctl_n58022x12328.gctx"
#ds_control <- parse_gctx(gctx_control)

# Explore the control data
#head(ds_control@mat)

# Read the compound metadata
#compound_meta <- read.delim("compoundinfo_beta.txt", header = TRUE, sep = "\t")

# Explore the metadata
#head(compound_meta)

# Filter compound metadata for AML-related compounds
#aml_compounds <- subset(compound_meta, grepl("AML", Target_Disease, ignore.case=TRUE))

# Filter the expression data (assuming drug IDs in ds@rid match the compound metadata)
#aml_expression_data <- ds@mat[ds@rid %in% aml_compounds$pert_id, ]
