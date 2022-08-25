# Evidence-of-Hybridisation-and-Introgression-in-Eucalyptus-Plantations-in-Ethiopia
Code for Samuel Rawlinson's Master's Project, "Evidence of Hybridisation and Introgression in Eucalyptus in Ethiopia", submitted to Imperial College London in August 2022 in partial fulfilment of the requirements of the degree of Master of Research Ecology, Evolution and Conservation

The "PCA.R" file takes an input of the vcf file of SNP data for the Eucalyptus population studies and the accompanying meta-data file. It uses these to conduct a PCA analysis.

The "Maps_for_Figure.R" file takes an inout of the meta-data file for the Eucalyptus population studies. It returns an output of a map on which the location of samples had been plotted.

The "STRUCTURE_input_and_MAF.R" file takes an input of the vcf file of SNP data for the Eucalyptus population studied. It uses the vcf file to produce a file in a format appropriate for input into STRUCTURE. It also contains the code for carrying out a minor allele frequency filter, which returns a list of loci which pass each filter.

The "Conservative Admixture Proportion Estimation.R" file takes an input of the q-value tables from the STRUCTURE output and summarises them, keeping the most conservative q-value estimates for each individual. It returns a file with the most conservative q-value estimates for each individual. It also returns a plot off these conservative q-values.

The "gghybrid_analysis.R" file takes a STRUCTURE input file as an input. It conductes all gghybrid analyses and filters for candidate adaptively introgressed loci based on the resultant genomic clines. It returns all gghybrid output files, a list of candidate SNPs in each transect of interest, as well as candidate SNPs which aresignificant in multiple transects.

The "gghybrid_input_random_sampling.R" file takes the gghybrid input file produced in "gghybrid_analysis.R" and outputs an excel file with a random list of 91 SNPs from said file.

The "plotting_pop_identity_against_altitude.R" file takes an input of the meta-data file, the output of "Conservative Admixture Proportion Estimation.R" and an esth file from gghybrid. It produces a model of hybrid index against altitude and plots it.

The "hybrid_value_model.R" file takes an input of the meta-data file, the output of "Conservative Admixture Proportion Estimation.R" and an esth file from gghybrid. It produces a model of admixture value against growth per year. It also produces a plot of average admixture value at each altitude.

The "Pair Comparison.R" file takes an input of the meta-data file, the output of "Conservative Admixture Proportion Estimation.R" and an esth file from gghybrid. It collates pairs of pure and admixed individuals from each site along transects of interest and conducts t-test to determine whether parental and admixed individuals differ in average growth per year.

The "Annotated Gene Freqeuency Analysis.R" file takes the input of a excel file which contains the details of gene annotation for candidate SNPs, and a second excel file which contains the same details for random subsets of SNPs. It produces plots for the number of SNPs mapping to genes and exons for each category, as well as plots showing the number of SNPs mapping to genes with different gene ontology annotations.
