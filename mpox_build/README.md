# Nextstrain repository for mpox virus

This repo is forked and customized from the original Nextstrain mpox build found here: [nextstrain.org/mpox](https://nextstrain.org/mpox).

## for this specific project:

This repo is customized for understanding mpox diversity and spread in LA County. 

To create the maximum likelihood phylogeny used in the manuscript, you can run: 

```
cd phylogenetic/
nextstrain build . --configfile my_configs/all_cladeII_config.yaml -p
nextstrain view .
```
This will create the nextstrain json, newick file, and alignment needed to continue the workflow. 