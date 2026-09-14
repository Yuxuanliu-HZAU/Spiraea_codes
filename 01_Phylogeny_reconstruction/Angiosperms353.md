# Angiosperms353
This page documents the running code for reconstructing nuclear phylogeny with Angiosperms353 dataset in our study.

## 01 Prepare public data and self-sequenced WGS raw data
The SRA accession list of public data has been recorded in our supporting information file (Table S4) and the self-sequenced WGS raw data has been released in [NGDC](https://ngdc.cncb.ac.cn/gsa/search?searchTerm=CRA049166) with the accession number CRA049166.

Download them and place them together in a directory before going forward to the following steps (we call this directory as `<input_data>` temporarily in this document).
The filenames shoud be `<sample_name>_1.fq.gz` and `<sample_name>_2.fq.gz` for paired-end data and `<sample_name>.fq.gz` for single-end data.

## 02 Prepare sample name list before running HybSuite

The sample name list used in Angiosperms353 dataset can be checked and downloaded [here]().

The details of how to prepare sample name list file can be founded in HybSuite manual [here](https://yuxuanliu-hzau.github.io/HybSuite.docs/docs/tutorial/#1-the-sample-list-file).

## 03 Run HybSuite

Assume that we have installed [HybSuite](https://github.com/Yuxuanliu-HZAU/HybSuite) successfully. Details of installing HybSuite can be found [here](https://yuxuanliu-hzau.github.io/HybSuite.docs/docs/installation/).

Then, run HybSuite with the command showed as follow:
```


```
