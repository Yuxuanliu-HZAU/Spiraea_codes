# Chloroplast
This page documents the running code for reconstructing chloroplast phylogeny in our study.

## 01 Prepare public data and self-sequenced WGS raw data
The SRA accession list of public data has been recorded in our supporting information file (Table S4) and the self-sequenced WGS raw data has been released in [NGDC](https://ngdc.cncb.ac.cn/gsa/search?searchTerm=CRA049166) with the accession number **CRA049166**.

Download them and place them together in a directory before going forward to the following steps (we call this directory as `<input_data>` temporarily in this document).
The filenames shoud be `<sample_name>_1.fq.gz` and `<sample_name>_2.fq.gz` for paired-end data and `<sample_name>.fq.gz` for single-end data.

## 02 Prepare sample name list before running HybSuite

The sample name list used in Chloroplast dataset can be checked and downloaded [here](https://github.com/Yuxuanliu-HZAU/Spiraea_codes/blob/main/01_Phylogeny_reconstruction/Sample_list.Chloroplast.txt).

The details of how to prepare sample name list file can be found in HybSuite manual [here](https://yuxuanliu-hzau.github.io/HybSuite.docs/docs/tutorial/#1-the-sample-list-file).

## 03 Prepare the Chloroplast target loci file

We have released the Chloroplast target loci file (79 loci) used in our study [here](https://github.com/Yuxuanliu-HZAU/Spiraea_codes/blob/main/01_Phylogeny_reconstruction/Reference_Chloroplast_CDS.fasta). Download it before running HybSuite.

## 04 Run HybSuite

Assume that we have installed [HybSuite](https://github.com/Yuxuanliu-HZAU/HybSuite) successfully. Details of installing HybSuite can be found [here](https://yuxuanliu-hzau.github.io/HybSuite.docs/docs/installation/).

Then, run HybSuite with the command showed as follow:

```
hybsuite full_pipeline \
-input_list <sample_name_list> \
-input_data <input_data> \
-output_dir <out_dir> \
-nt 7 -process 8 -PH 1 -sp_tree 4 \
-seqs_min_length 50 -seqs_min_sample_coverage 0.1 \
-mafft_algorithm linsi -mafft_maxiterate 1000 \
-t <target_loci_file> -check TRUE
```
- `<sample_name_list>`: The sample list file.
- `<input_data>`: The directory containing all input public data and self-sequenced WGS data.
- `<out_dir>`: The output directory.
- `<target_loci_fasta>`: The Chloroplast target loci file (79 loci).
