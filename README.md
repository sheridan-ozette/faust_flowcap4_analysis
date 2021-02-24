# FlowCAP IV FAUST Analysis

This repository contains scripts to reproduce the
FAUST analysis of the FlowCAP IV dataset.

The FAUST analysis modifies the FloReMi method reported in
[FloReMi: Flow density survival regression using minimal feature redundancy](https://pubmed.ncbi.nlm.nih.gov/26243673/).
The analysis here uses the pre-processing script `1_preprocessing.R` from
[https://github.com/SofieVG/FloReMi](https://github.com/SofieVG/FloReMi)
(modified to point to the local file system).
As in FloReMi, we also use a random survival forest for the prediction modeling here.

To reproduce the FAUST analysis of the FlowCAP IV data, run the scripts in numeric order.

Scripts expect that all FlowCAP IV data downloaded from flow repository to be stored in the analysis directory, in
a sub-directory called "FR-FCM-ZZ99".

The scripts will generate intermediate RDS files and sub-directories in the analysis directory needed to
run the analysis.