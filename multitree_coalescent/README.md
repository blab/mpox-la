# Running multitree coalescent analysis

We can analyze local outbreak clusters of sequences in order to understand mpox transmission dynamics via a Bayesian coalescent approach that was previously developed in [Müller et al, 2021](https://www.science.org/doi/10.1126/scitranslmed.abf0202). 

## Analysis workflow

1. Once the local outbreak cluster information and metadata have been combined and moved to the [`data`](data/) subfolder, the xmls for this project can be created using the two matlab scripts found under the [`scripts`](scripts/) subfolder. Use `buildMulticoal_case_prior.m` to create the xmls for the skyline prior informed by case counts and `buildMulticoal_Skygrowth.m` for the skygrowth prior without case counts. 

2. To add the cases into your xmls with the case-based prior, you can use the script [`add_cases_to_empirical_xml.ipynb`](scripts/add_cases_to_emprical_xml.ipynb) to add those in. All the xmls used in the manuscript, with the various sensitivity analyses can be found in [`xmls`](xmls/)

3. The xmls can then be run using the custom jar file found in this folder using the following command: 
`java -jar mpox_la.jar {xml.name}`. you can also add in different memory specifications if you're analyzing a lot of sequences such as: `-Xms8g -Xmx16g`. The code for this custom jar file can be found here: https://github.com/miparedes/mab/tree/master 

4. After checking for convergence using [Tracer](https://www.beast2.org/tracer-2/), analyze the results to create the manuscript figures using the scripts for each corresponding figure found in [`scripts`](scripts/).

All results for main and supplementary analyses can be found in [`results`](results/)