import pandas as pd
from cmapPy.pandasGEXpress.parse_gctx import parse

# Load the compound metadata file
compound_metadata = pd.read_csv("compoundinfo_beta.txt", sep="\t")
print("Loaded compound metadata:")
print(compound_metadata.head())

# Filter compounds for breast cancer relevance (modify as needed)
breast_cancer_compounds = compound_metadata[
    compound_metadata['cmap_name'].str.contains("breast", case=False, na=False) |
    compound_metadata['compound_aliases'].str.contains("breast", case=False, na=False)
]
print("Filtered compounds for breast cancer:")
print(breast_cancer_compounds.head())





# Include a list of known breast cancer drugs or related terms
breast_cancer_terms = ['tamoxifen', 'trastuzumab', 'aromatase', 'breast']

breast_cancer_compounds = compound_metadata[
    compound_metadata['cmap_name'].str.contains('|'.join(breast_cancer_terms), case=False, na=False) |
    compound_metadata['compound_aliases'].str.contains('|'.join(breast_cancer_terms), case=False, na=False)
]

print("Column names in compound metadata:", compound_metadata.columns)

print("Unique cmap_name values:", compound_metadata['cmap_name'].unique()[:10])
print("Unique compound_aliases values:", compound_metadata['compound_aliases'].unique()[:10])






# Extract unique pert_id for breast cancer compounds
subset_cid = breast_cancer_compounds['pert_id'].unique()
subset_cid = [id.strip().lower() for id in subset_cid]
print("Subset of Perturbagen IDs:", subset_cid[:5])

# Specify the path to the Level 5 GCTX file
gctx_treatment = "level5_beta_trt_cp_n720216x12328.gctx"

# Load the GCTX file
gctx_data = parse(gctx_treatment)

# Extract column metadata
gctx_metadata = gctx_data.col_metadata_df
print("GCTX column metadata:")
print(gctx_metadata.head())

# Extract and clean all column IDs
all_cid = gctx_metadata.index.values
cleaned_all_cid = [id.split(":")[1].strip().lower() if ":" in id else id.strip().lower() for id in all_cid]
print("Example cleaned column IDs:", cleaned_all_cid[:5])

# Find matching perturbagen IDs
matching_cid = set(subset_cid).intersection(set(cleaned_all_cid))
print("Matching Perturbagen IDs:", matching_cid)

if matching_cid:
    # Extract expression data for matching compounds
    expression_data = parse(gctx_treatment, cid=list(matching_cid))
    expression_df = pd.DataFrame(expression_data.data_df)
    print("Extracted expression data shape:", expression_df.shape)

    # Save to CSV
    expression_df.to_csv("filtered_breast_cancer_expression_data.csv", index=False)
    print("Filtered expression data saved.")
else:
    print("No matching perturbagen IDs found.")


actual_cid_in_gctx = set(gctx_metadata.index.values)
print(f"Does 'brd-k93754473' exist in GCTX metadata? {'brd-k93754473' in actual_cid_in_gctx}")

# Print all unique CIDs in the GCTX file
actual_cid_in_gctx = set(gctx_metadata.index.values)
print("All CIDs in GCTX metadata:", list(actual_cid_in_gctx)[:10])  # Print only the first 10 for brevity

