# cgMLST-Corynebacterium-striatum

# cgMLST-Corynebacterium-striatum

The objective of this repository is to describe how we created a cgMLST for C. striatum. This scheme was created with the ChewBBACA pipeline (link below). 

chewBBACA
To download ChewBBACA access the link:
https://github.com/B-UMMI/chewBBACA_tutorial

* Softwares and Downloads
* conda install ChewBBACA

Workflow used to create the scheme
* Step 1: Scheme Creation of cgMLST
* Step 2: Allele calling
* Step 3: Scheme Validation (Allele call) of cgMLST
* Step 4: Extracting the Matrix target genes
* Step 5: Minimum Spanning Tree (MST)
* Step 6: Graphical evaluation of the scheme

# Step 1: Schema Creation
Selection of complete genomes for schema creation
As of August 29, 2021, 271 publicly available genome sequences were available at the NCBI (National Center for Biotechnology Information) genome sequence repository (https://www.ncbi.nlm.nih.gov/assembly).
A list of 271 publicly available genomes used to create this protocol can be found in the file 271-publicly-available-genome.txt.
C. striatum WP1a (RefSeq assembly accession: GCF_004138065.1) was used by Prodigal algorithm as reference to recognize coding sequences (CDs). Prodigal generated the GCF_004138065.1.trn file at this step.

Command:
```
#create schema
chewBBACA.py CreateSchema -i publicly-available-genome --cpu 46 -o schema_seed --ptf GCF_004138065.1.trn
```
The above command uses 46 CPU and creates a preliminary scheme (wgMLST) in the schema_seed folder using the trained product folder GCF_004138065.1.trn that was generated using the reference genome WP1a (GCF_004138065.1) and 271 publicly available genome sequences. The wgMLST scheme generated contained 4934 loci based on the 271 complete genomes.

# Step 2: Allele calling
In this step the allele calling is performed using the resulting set of loci determined in step 1.

Command:
```
#run allelecall
chewBBACA.py AlleleCall -i publicly-available-genome -g schema_seed/ -o results_cg --cpu 46 --ptf GCF_004138065.1.trn
```
The allele calling used the default BLAST Score Ratio (BSR) threshold of 0.6.

# Step 2.1: Paralog detection
In this step genes considered paralogous from result of the allelecall (see above) are removed
Command:
```
#run remove genes
chewBBACA.py RemoveGenes -i results_cg/results_ 20220321T145123/results_alleles.tsv -g results_cg/results_ 20220321T145123/RepeatedLoci.txt -o alleleCallMatrix_cg
```
In this step, 66 genes were identified as possible paralogs and were removed from further analysis.

# Step 2.2: Genome Quality Control
In this step we define a Threshold for the scheme that limits the loss of loci targets defined in the previous steps per genome and excludes genomes considered to be of low quality due to significant loci absence.
With this analysis we define the percentage of loci that will constitute the scheme based on how many targets we want to keep in this phase. For example, 100%, 99.5%, 99% and 95% of the loci may present in the set of high-quality genomes. This is one of the main steps in defining the cgMLST scheme targets.
Command:
```
#run test genome quality
chewBBACA.py TestGenomeQuality -i alleleCallMatrix_cg.tsv -n 13 -t 300 -s 5
```
In this stage we chose the loci present in 95% of the publicly available genomes and the Threshold 15 that limited the loss of the loci in the genomes. In the Threshold 15 a set of 1917 loci were found to be present in 95% the analyzed genomes. In this Threshold (15) 35 genomes were removed due to loss of loci.
Command:
```
#run ExtractCgMLST
chewBBACA.py ExtractCgMLST -i alleleCallMatrix_cg.tsv -o cgMLST_15 --g GenomeRemoved15thr.txt  --t 0.95
```
This script selects all * loci * present in the selected * Threshold *. The value * t * is associated with the percentage of * loci * we want to be present in the set of genomes, for example: * t 1.0 * selects all * loci * present in the * Threshold * chosen in all genomes ie those present in 100% of genomes at that * Threshold *. Subsequently a cgMLST_15 folder receives the result of the allelic profile for each of the 1917 candidate * loci * (allelic profile matrix). The file in this folder (cgMLST.tsv) contains the allelic profile of 1917 selected * loci * and will be used to create the core gene list.

# Step 2.3: Creating the core gene list
This command selects all target genes from the "cgMLST.tsv" spreadsheet.
```
head -n 1 cgMLST.tsv > Genes_95%_Core_15.txt
```
This step generated the file Genes_95%_Core_15.txt. This list needs to be transposed so that each core gene name is reported in a single line:
Command:
```
#transpose table
datamash -W transpose < Genes_95%_Core_15.txt > Genes_Core_Al.txt 
```
This step generated the file > Genes_Core_Al.txt
You can see the list file with 1917 target genes at Genes_Core_Al.txt and for the subsequent steps we added the full path to each locus fasta file.
This list Genes_Core_Al.txt was then modified so that each name was preceeded by schema_seed:
Command:
```
tail -n+1 Genes_Core_Al.txt | cut -f2 | perl -pe 's: :\n:g' | sort -Vu | awk '{print("schema_seed/"$1)}' > list_genes_core.txt
```

# Step 3: Scheme Validation (Allele calling)
For the validation step we selected 30 C. striatum strains from a nosocomial outbreak and 2 epidemiologically unrelated outgroup C. striatum strains. The list of all the Validation genomes used can be found in the file Genomes_Validation.txt.
We this set of genomes (32 drafts genomes) we repeated the allele call step using only the 1917 candidate target genes.
Command:
```
chewBBACA.py AlleleCall -i ../Genomes_Validation/ --gl list_genes_core.txt -o ../results_all --cpu 46  -g ../schema_seed/schema_seed –ptf GCF_004138065.1.trn
```
The folder Genomes_Validation contains the 32 validation genomes used for validation of the scheme.

# Step 3.1: Concatenate the allelic profiles
The purpose of this step is to concatenate the matrix of the loci that defined the scheme and the matrix of the loci from the validation genomes. Thus, to concatenate the allelic profile matrix obtained from the creation of the scheme cgMLST_15/cgMLST.tsv with the matrix obtained for the validation genomes results_all/results_20220321T170302/results_alleles.tsv. The following command was used:
Command:
```
#create header
head -n 1 cgMLST_15/cgMLST.tsv > cgMLST_all.tsv
```
```
#concatenate
grep -v ^FILE cgMLST_15/cgMLST.tsv results_all/results_ 20220321T170302/results_alleles.tsv >> cgMLST_all.tsv
```
cgMLST_all.tsv file (cgMLST_all.tsv) contains the allelic profile of the 32 validation genomes.

# Step 3.2: Evaluation of genome quality
After concatenation, we used the TestGenomeQuality to assess the impact of each validation genome on candidate loci in order to exclude low quality validation genomes. In this step you may need to define a new Threshold, as well as a new value of the parameter t, because loci that remain after the filters are the ones that will constituted the final scheme.
```
chewBBACA.py TestGenomeQuality -i cgMLST_all.tsv -n 13 -t 300 -s 5
```

# Step 4: Extracting the Matrix loci
At this step we chose loci present in 99% (t0.99) of the validation genomes and the Threshold 110 to limit the loss of the loci in the genomes.
In Threshold 1776 a set of 268 loci were found to be present in 99% the validation genomes.
Using Threshold (110) 0 genome was removed due to absence of loci targets.
Command:
```
chewBBACA.py ExtractCgMLST -i cgMLST_all.tsv -o cgMLST_110 --t 0.99 
```
The folder with the output file can be found at: cgMLST_110. This folder contains four files "cgMLST.tsv; cgMLSTschema.txt; mdata_stats.tsv and Presence_Absence.tsv".
The cgMLST targets can be found at: cgMLST_110/cgMLSTschema.txt It contains the list of 1776 genes in the core genome defined as targets for this scheme. 

# Step 5: Minimum Spanning Tree
For the visualization of results, minimum spanning trees were buitl. Based on the allelic profiles obtained by the cgMLST scheme for each of the 268 genomes minimum spanning trees (MST) were constructed using the software GrapeTree (version 2.1.0) (https://github.com/achtman-lab/GrapeTree/releases) with parameters implemented in MSTree v2 ignoring missing values for the entire strain collection. The cgMLST_110/cgMLST.tsv file contains the allelic profile of the 268 genomes typed by cgMLST.

# Step 6: Graphical evaluation of the scheme
To assess the variability of the gene targets of cgMLST as well explore and evaluate the type and extent of allelic variation detected in each of the chosen loci. We run this script and graphically visualize the data in a series of html files.
Command:
```
chewBBACA.py SchemaEvaluator -i schema_seed/schema_seed   --cpu 6 -o rms
```
