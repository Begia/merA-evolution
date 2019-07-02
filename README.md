This repository contains code and files for the study:

"Swift Evolutionary Response of Microbes to Historical Rise in Anthropogenic Mercury in the Northern Hemisphere",
by Matti O. Ruuskanen, Stéphane Aris-Brosou, and Alexandre J. Poulain

## Included files

1. Dating files
	- All primary dating information collected in this study
	- Dating for Kevojarvi, Lake Hazen, and Pocket Lake cores was not completed in this study with 210-Pb dating
	- Folder: "dating"
	- Filenames: "Dating_profiles_*SITE*_*YEAR*.xlsx"
	- Data summarized in "Core_dating_raw_data.csv"

3. High-throughput sequence data processing
	- Scripts that were run on the CAC cluster to produce the primary data from Illumina read files
	- Folder: "sequence_processing"
	- Filenames 
		- Main pipelines: "merA_runs.sh" and "rpoB_runs.sh"
		- Files used in the scripts 
			- Summary of samples: "merA_samples.txt" and "rpoB_samples.txt"
			- Collection of sequences used in chimera detection: "merA_614_all_curated_filtered.fasta" and "rpoB_out_90.fa"
			- Hidden Markov Models used in QC of the sequences: "merA_614_all.hmm" and "RNA_pol_Rpb2_45.hmm"
			- Custom accessory programs used in the pipelines: "uc_to_otu_table.py", "transpose_table" and "swarm_construct_otu_table.R"

4. Analysis script and related input files
	- Main analysis script run with R and its input files
	- Folder: "analysis_and_input_files"
	- Filenames
		- Main pipeline: "analysis_script.R"
		- Files used in the script
			- Input files for analysis of both merA and rpoB:
				- Biom (variant table) files: "merA_swarm_otu_table_no_singletons.biom" and "rpoB_swarm_otu_table_no_singletons.biom"
				- Final variant sequences and phylogenetic trees constructed from them (from FastTree): "merA_trimal_v2.fasta/tre" and "rpoB_trimal_v2.fasta/tre"
			- Metadata mapping file: "mapping_corrected_concatenated.txt"
			- Combined tot-Hg measurements (also includes Canadian samples): "finland_sediment_totHg_overview.csv"
			- Combined dating results: "Core_dating_raw_data.csv"
			- Combined DNA concentration results: "FIN_CAN_DNA_conc_mixed.csv"
			- ddPCR data (4 files, two for each merA and glnA): "June_8th_qmerA_rows_*.csv" and "June_9th_glnA_rows_*.csv"
			- Accessory tables: "sample_coding_levels.csv" and "sampling_site_coordinates.csv"
			- Data from other studies for comparisons
				- Clackett et al. 2018; tree-ring tot-Hg data: "Clackett_Hg_data.csv"
				- Thienpoint et al. 2016; tot-Hg data from Pocket Lake for dating cross-comparison: "Pocket_lake_core_2014_THg_dating_data.csv"

5. Evolutionary analyses
	- Files related to the Bayesian analysis of molecular sequences using MCMC ran with BEAST 1.8.0
	- Folder: "evolutionary_analyses"
	- Filenames:
		- Input FASTA files with dated sequences for each core: "*GENE*_*SITE*_ref.fasta"
		- Xml parameter files for BEAST: "*GENE*_*SITE*_ref_25.xml"
		- Combined "target" trees from TreeAnnotator: "ref_25_*GENE*_*SITE*_combined.trees.tre"
		- Bayesian skyline data: "ref_25_*GENE*_*SITE*_skyline.txt"