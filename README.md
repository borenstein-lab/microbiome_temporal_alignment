# microbiome_tempoarl_alignment

Welocme! 

This repository contains the analysis code for: "Tempoarl Alignment of Microbiome Data" (Armoni & Borenstein, in preperation).

***Abstract***

A major challenge in working with longitudinal data when studying some temporal process is the fact that differences in pace and dynamics might overshadow similarities between processes. In the case of longitudinal microbiome data, this may hinder efforts to characterize common temporal trends across individuals or to harness temporal information to better understand the link between the microbiome and the host. One possible solution to this challenge lies in the field of ‘temporal alignment’ – an approach for optimally aligning longitudinal samples obtained from processes that may vary in pace. In this work we investigate the use of alignment-based analysis in the microbiome domain, focusing on microbiome data from infants in their first years of life. Our analyses center around two main use-cases: First, using the overall alignment score as a measure of the similarity between microbiome developmental trajectories, and showing that this measure can capture biological differences between individuals. Second, using the specific matching obtained between pairs of samples in the alignment to highlight changes in pace and temporal dynamics, showing that it can be utilized to predict the age of infants based on their microbiome and to uncover developmental delays. Combined, our findings serve as a proof-of-concept for the use of temporal alignment as an important and beneficial tool in future longitudinal microbiome studies. 

**Scripts and files overview**

*Notebooks:* 
each notebook is named according to the corresponding figure numbers it was used to generate.

*Code:*

* Config.py - the "config" object used to define and store configurations of the analysis.
* age_analysis.py - scripts for alignment-based age prediction.
* alignment.py - wraper functions and objects for the alignment procedure and result (alignment itself is calculated using the DTW library)
* data_tables.py - functions and objects used to store, process and filter the raw data. 
* distance.py - wraper functions for distance matrix calculations.
* experiments.py - utility functions for experiments conducted in the analysis. 
* interpoaltion.py - interpolation functions. 
* simulation.py - utility functions for simulating data (only "shuffle_data" used in paper)
