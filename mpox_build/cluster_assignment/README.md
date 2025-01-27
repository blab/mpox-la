## Extracting out local outbreak clusters from ML tree

### About the analysis
In order to conduct our multi-tree unstructured coalescent analysis, we need to first identify local outbreak clusters in LA via a maximum parsimony from our maximum likelihood phylogeny produced via the Nextstrain pipeline

### Running the analysis

1. First we need to format the nextstrain metadata in order to import it easier into matlab. To do so, run `./scripts/move_focus_area_column.sh metadata.tsv` which puts the column of interest into the desired column location and outputs `updated_metadata.tsv`. 

2. Open `scripts/la_mpox_parsimonyclusters.m` in matlab and update the file location for both the newick tree (found in the results folder of the nextstrain build) and the updated_metadata.tsv. This will output a file called `la_clusters.tsv`

3. Create the combined metadata and cluster infomation by running `scripts/combining_metadata_and_cluster_assignments.ipynb` and move the resulting `la_clusters_with_metadata.tsv` that is found in `clustering_results/` to `../../multitree_coalescent/data/` for further processing. 


### Figure S3: Sensitivity analysis of the relationship between background sequence proportion and number of idenified clusters

1. To test the sensitivity of our clustering to different proportions of background sequences as well as sequences within LA county, run `getParsimonyClusterComparison.m`
