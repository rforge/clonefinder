---
author: Kevin R. Coombes
date: 12 October 2018
title: Auxiliary Material for CloneFinfer
---

# Overview
This directory contains a combination of material that either
(1) constructs data sets that are part of the package or (2) performs
some of the analyses described in our manuscript.

## Package Material

1. `chlens.txt` is a simple table holding the lengths In base pairs)
   of human chromosomes. It is based on build [UNKNOWN} and was
   initialy obtained from [UNKNOWN].
2. `makeChlens.R` is a script that reads the previous file, converts
   it to a useful form, and saves it in the `sysdata.rda` file so it
   can be used by the package. In particular, it is used by functions
   that simulate `Tumor` objects to create  realistic segmentation
   data.
3. `makeTestExamples.R` is a script that simulates copy number and
   mutation data used for regression testing of the algorithms during
   package development. 

## Analyses

1. `simgen.R` is the script that created the original simulated data
   sets used in our paper to compare the CloneSeeker algorithm to other
   algorithms intended to detect clones.
2. `running.R` analyzes all th simulated data.
3. `asswssing.R` compares the infered results feomthe the algorithms
   to the true properties of the simulated adta.
4. `clinical.R` applies th CloneSeeker alorithm to CLL data.
5. `cnv_survival.R` tsts the associations of subclones with clinical
   variables in CLL.
6. ``figures.R` generates the figures in the manuscript.
