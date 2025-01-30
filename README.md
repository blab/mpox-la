# mpox_la

## Viral introductions and return to baseline sexual behaviors maintain low-level mpox incidence in Los Angeles County, USA, 2023-2024

### Abstract:

In 2022, mpox clade IIb disseminated around the world, causing outbreaks in more than 117 countries. Despite the decay of the 2022 epidemic and the expected accumulation of immunity within queer sexual networks, mpox continues to persist at low incidence in North America without extinction, raising concerns of future outbreaks. We combined phylodynamic inference and microsimulation modeling to understand the heterogeneous dynamics governing local mpox persistence in Los Angeles (LA) County from 2023–2024. Our Bayesian phylodynamic analysis revealed a time-varying pattern of viral importations into the county that seeded a heavy-tailed distribution of mpox outbreak clusters. Our phylodynamics-informed microsimulation model demonstrated that the persistent number of mpox cases in LA County can be explained by a combination of waves of viral introductions and a return to near-baseline sexual behaviors that were altered during the 2022 epidemic. Finally, our counterfactual scenario modeling showed that public health interventions that either promote increased isolation of symptomatic, infectious individuals or enact behavior-modifying campaigns during the periods with the highest viral importation intensity are both actionable and effective at curbing mpox cases. Our work highlights the heterogeneous factors that maintain present-day mpox dynamics in a large, urban US county and describes how to leverage these results to design timely and community-centered public health interventions.


-----------


This repository contains the analytic code needed to reproduce the phylodynamic results from the above paper. 

## Primary analysis outline

1. To start, begin with the folder [`mpox_build`](mpox_build/) to run the maximum likelihood analysis and create the temporally resolved phylogeny of clade II mpox. The ML phylogeny used in the manuscript can be found at https://nextstrain.org/groups/blab/mpox/allcladeIIseqs?c=focus_areas

2. Identify the local outbreak clusters for Los Angeles County using maximum parsimony found in [`cluster_assignment`](mpox_build/cluster_assignment/)

3. Once the clusters have been identified and the cluster and sequence metadata have been combined, switch to the [`multitree_coalescent`](multitree_coalescent/) folder which contains the scripts to create the BEAST2 XMLs that can be run using the provided [`jar file`](multitree_coalescent/mpox_la.jar). All the xmls used in the above analysis can be found in [`xmls`](multitree_coalescent/xmls/)

4. Analyze the results to create the manuscript figures using the scripts for each corresponding figure found in [`scripts`](scripts/)

5. To run simulations of mpox dynamics with superspreading which are part of our validation analysis for our multitree approach, you can find the necessary code in [`simulations`](simulations/)

6. of note, while the code for the microsimulation model used in our study is found on a seperate repo, the results of the microsim model used in our manuscript can be found under [`data`](data/) and used to recreate the main figures.

Of note, most of the xml and result files have been compressed. To decompress use the following format in the command line:
`xz --decompress --keep {file.name}`