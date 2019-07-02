#!/bin/bash
#SBATCH -J rpoB_runs
#SBATCH -c 10                             # Number of CPUS requested. If omitted, the default is 1 CPU.
#SBATCH --mem=100G                    # Memory requested in megabytes. If omitted, the default is 1024 MB.
#SBATCH -t 2-00:00:00      # How long will your job run for? If omitted, the default is 3 hours. days-hours:minutes:seconds

#enter the qiime1 environment with the correct versions of programs pre-installed
source activate qiime1

#run individual files with GNU Parallel
for sample in `awk '{print $1}' ~/matti/rpoB_2018/rpoB_samples.txt`
do 
	#1. join paired ends, overlap > 100 bp, highest p-value stringency
	( pear -f ~/matti/rpoB_2018/MI.M03555_0327.001.*.${sample}_rpoB_R1.fastq.gz -r ~/matti/rpoB_2018/MI.M03555_0327.001.*.${sample}_rpoB_R2.fastq.gz -o ~/matti/rpoB_2018/rpoB_${sample}_paired -p 0.0001 -v 100
	#2. check sequence quality with fastqc
	mkdir -p ~/matti/rpoB_2018/rpoB_fastqc_${sample}
	fastqc -o ~/matti/rpoB_2018/rpoB_fastqc_${sample}/ ~/matti/rpoB_2018/rpoB_${sample}_paired.assembled.fastq 
	#3. convert fastq to fasta and qual files
	convert_fastaqual_fastq.py -c fastq_to_fastaqual -f ~/matti/rpoB_2018/rpoB_${sample}_paired.assembled.fastq -o ~/matti/rpoB_2018/
	#4. check the mapping file
	validate_mapping_file.py -m ~/matti/rpoB_2018/rpoB_${sample}_mapping.txt -o ~/matti/rpoB_2018/ -b
	#5. split libraries according to the mapping file, remove both forward and reverse primers, barcode length 8, min length 150, quality filter to remove N and nucleotides with avg quality scores under 28 in a window of 2
	split_libraries.py -m ~/matti/rpoB_2018/rpoB_${sample}_mapping_corrected.txt -f ~/matti/rpoB_2018/rpoB_${sample}_paired.assembled.fna -q ~/matti/rpoB_2018/rpoB_${sample}_paired.assembled.qual -o ~/matti/rpoB_2018/rpoB_${sample}/ -z truncate_remove -b 0 -l 150 -x -s 28 -w 2
	#6. chimera detection and filtering with vsearch uchime
	vsearch --uchime_ref ~/matti/rpoB_2018/rpoB_${sample}/seqs.fna --db ~/matti/rpoB_2018/rpoB_out_90.fa --nonchimeras ~/matti/rpoB_2018/rpoB_${sample}/seqs_nonchimeras.fna --uchimeout ~/matti/rpoB_2018/rpoB_${sample}/chimeras_log.txt
	#7. translate the sequences in 3 frames for hmmsearch
	~/bin/EMBOSS-6.5.7/emboss/transeq -sequence ~/matti/rpoB_2018/rpoB_${sample}/seqs_nonchimeras.fna -outseq ~/matti/rpoB_2018/rpoB_${sample}/seqs_nonchimeras_aa.fna -frame=F -clean
	#8. remove translations with stop characters
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ~/matti/rpoB_2018/rpoB_${sample}/seqs_nonchimeras_aa.fna > ~/matti/rpoB_2018/rpoB_${sample}/seqs_nonchimeras_aa_linear.fna
	awk '!/^>/ { next } { getline seq } {if (seq !~ /X/) print $0 "\n" seq }' ~/matti/rpoB_2018/rpoB_${sample}/seqs_nonchimeras_aa_linear.fna > ~/matti/rpoB_2018/rpoB_${sample}/seqs_nonchimeras_aa.fna
	#9. query the RNA polymerase beta subunit external 1 domain HMM (from PF10385) against the sequences with hmmsearch
	hmmsearch -o ~/matti/rpoB_2018/rpoB_${sample}/seqs_nonchimeras_hits.txt --tblout ~/matti/rpoB_2018/rpoB_${sample}/seqs_nonchimeras_hits_table.txt --incE 1.0e-5 -E 1.0e-5 --noali ~/matti/rpoB_2018/RNA_pol_Rpb2_45.hmm ~/matti/rpoB_2018/rpoB_${sample}/seqs_nonchimeras_aa.fna
	#10. make lists of hits and subset the files
	awk -F "_" '{print ($1"_"$2)}' ~/matti/rpoB_2018/rpoB_${sample}/seqs_nonchimeras_hits_table.txt | grep "_" | sort | uniq > ~/matti/rpoB_2018/rpoB_${sample}/seqs_nonchimeras_hits.list
	faSomeRecords ~/matti/rpoB_2018/rpoB_${sample}/seqs_nonchimeras.fna ~/matti/rpoB_2018/rpoB_${sample}/seqs_nonchimeras_hits.list ~/matti/rpoB_2018/rpoB_${sample}/seqs_nonchimeras_hits.fna ) &
done
wait

#11. combine files and dereplicate with vsearch
mkdir -p ~/matti/rpoB_2018/rpoB_all
cat ~/matti/rpoB_2018/rpoB_*/seqs_nonchimeras_hits.fna > ~/matti/rpoB_2018/rpoB_all/seqs_nonchimeras_hits.fna
vsearch --derep_fulllength ~/matti/rpoB_2018/rpoB_all/seqs_nonchimeras_hits.fna --output ~/matti/rpoB_2018/rpoB_all/dereplicated.fna -uc ~/matti/rpoB_2018/rpoB_all/dereplicated.uc --sizeout
source activate anvio4
uc_to_OTU_table.py -c ~/matti/rpoB_2018/rpoB_all/dereplicated.uc -o ~/matti/rpoB_2018/rpoB_all/ -f ~/matti/rpoB_2018/rpoB_all/dereplicated.fna
source activate qiime1
awk -F ";|=" '/^>/{print $3}' ~/matti/rpoB_2018/rpoB_all/dereplicated.fna > ~/matti/rpoB_2018/rpoB_all/dereplicated_abundance.txt
awk 'FNR==NR{a[++i]=$1;next} /^>/{print $1 "_" a[++j]}!/^>/{print}' ~/matti/rpoB_2018/rpoB_all/dereplicated_abundance.txt ~/matti/rpoB_2018/rpoB_all/dereplicated_renamed.fasta > ~/matti/rpoB_2018/rpoB_all/dereplicated_abundance.fasta

#12. cluster sequences with Swarm v2 and make an OTU table with custom perl and R scripts
mkdir -p ~/matti/rpoB_2018/rpoB_all/swarm
swarm-2.1.9-linux-x86_64 -f -o ~/matti/rpoB_2018/rpoB_all/swarm/otu_map.txt -t 10 -s ~/matti/rpoB_2018/rpoB_all/swarm/statistics.txt ~/matti/rpoB_2018/rpoB_all/dereplicated_abundance.fasta
#transpose the otu table to a better format for the constructor script (perl script from https://stackoverflow.com/questions/1729824/an-efficient-way-to-transpose-a-file-in-bash) 
transpose_table ~/matti/rpoB_2018/rpoB_all/otu_table.txt > ~/matti/rpoB_2018/rpoB_all/otu_table_t.txt
sed -i -e '1s/Sample/OTU/' ~/matti/rpoB_2018/rpoB_all/otu_table_t.txt
swarm_construct_otu_table.R -i ~/matti/rpoB_2018/rpoB_all/swarm/otu_map.txt -t ~/matti/rpoB_2018/rpoB_all/otu_table_t.txt -o ~/matti/rpoB_2018/rpoB_all/swarm
awk -v OFS="\t" '$1=$1' ~/matti/rpoB_2018/rpoB_all/swarm/swarm_otu_table.txt > ~/matti/rpoB_2018/rpoB_all/swarm/swarm_otu_table_with_singletons.txt
biom convert -i ~/matti/rpoB_2018/rpoB_all/swarm/swarm_otu_table_with_singletons.txt -o ~/matti/rpoB_2018/rpoB_all/swarm_otu_table_with_singletons.biom --table-type="OTU table" --to-json

#13. remove singleton OTUs from the OTU table and subset the swarm fasta file to no singleton cluster seeds
awk '{for(i=2; i<=NF;i++) j+=$i; if(j==1){j=0; next}; print $0; j=0}' ~/matti/rpoB_2018/rpoB_all/swarm/swarm_otu_table.txt > ~/matti/rpoB_2018/rpoB_all/swarm/swarm_otu_table_no_singletons.tmp
awk -v OFS="\t" '$1=$1' ~/matti/rpoB_2018/rpoB_all/swarm/swarm_otu_table_no_singletons.tmp > ~/matti/rpoB_2018/rpoB_all/swarm/swarm_otu_table_no_singletons.txt
awk '{print $1}' ~/matti/rpoB_2018/rpoB_all/swarm/swarm_otu_table_no_singletons.txt > ~/matti/rpoB_2018/rpoB_all/swarm/swarm_cluster_seed_list.tmp
faSomeRecords ~/matti/rpoB_2018/rpoB_all/dereplicated_renamed.fasta ~/matti/rpoB_2018/rpoB_all/swarm/swarm_cluster_seed_list.tmp ~/matti/rpoB_2018/rpoB_all/swarm_cluster_seeds.fasta
rm ~/matti/rpoB_2018/rpoB_all/swarm/*.tmp

#14. convert OTU table to json
biom convert -i ~/matti/rpoB_2018/rpoB_all/swarm/swarm_otu_table_no_singletons.txt -o ~/matti/rpoB_2018/rpoB_all/swarm_otu_table_no_singletons.biom --table-type="OTU table" --to-json

#15. add outgroups for alignment
cat ~/matti/rpoB_2018/rpoB_all/swarm_cluster_seeds.fasta ~/matti/rpoB_2018/rpoB_outgroups_truncated.fa > ~/matti/rpoB_2018/rpoB_all/swarm_cluster_seeds_outgroups.fasta

#16. align the sequences with translatorX and refine the alignment with muscle
perl ~/.local/bin/translatorx_vLocal.pl -i ~/matti/rpoB_2018/rpoB_all/swarm_cluster_seeds_outgroups.fasta -o ~/matti/rpoB_2018/rpoB_all/rpoB_translatorx -c 11 -t T
muscle -in ~/matti/rpoB_2018/rpoB_all/rpoB_translatorx.nt_ali.fasta -out ~/matti/rpoB_2018/rpoB_all/rpoB_translatorx.nt_ali_refined.fasta -refine

#17. manually curate the alignment in JalView (MSA viewer)
#18. realign and refine a trimmed alignment
perl ~/.local/bin/translatorx_vLocal.pl -i ~/matti/rpoB_2018/rpoB_all/rpoB_trimmed.fasta -o ~/matti/rpoB_2018/rpoB_all/rpoB_trimmed_translatorx_mafft -c 11 -p F

#19. continue analysis in R from ~/matti/rpoB_2018/rpoB_all/rpoB_trimmed_translatorx_mafft.nt_ali.fasta