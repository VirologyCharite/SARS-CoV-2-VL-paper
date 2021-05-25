[![DOI](https://zenodo.org/badge/348435588.svg)](https://zenodo.org/badge/latestdoi/348435588)

# Estimating infectiousness throughout SARS-CoV-2 infection course

This repository contains data and code used in the research for "Estimating
infectiousness throughout SARS-CoV-2 infection course" published in
_Science_ on May 25, 2021.

# Data

The following files are available in the `data` directory.

## viral-load-with-negatives.tsv.zip

A [zip](https://en.wikipedia.org/wiki/ZIP_(file_format))
[TAB-separated values](https://en.wikipedia.org/wiki/Tab-separated_values)
file containing 879,624 RT-PCR results.  The same data is also available in
compressed [bzip2](https://en.wikipedia.org/wiki/Bzip2) format, in
[viral-load-with-negatives.tsv.bz2](data/viral-load-with-negatives.tsv.bz2).

Each row has 13 fields, as follows

1. `log10Load`: log10 viral load. `-1` represents negative tests.
1. `Ct`: The RT-PCR [cycle threshold](https://duckduckgo.com/?q=pcr+cycle+threshold)
1. `PCR`: Either "LC480" or "T2" accordign to the RT-PCR system used.
1. `Date`: The "YYYY-MM-DD" date the sample arrived at the diagnostic
   facility.
1. `Age`: The subject's age on the day the RT-PCR was done, rounded to one
   decimal place (to maintain anonymity)
1. `TestCentre`: A secure one-way hash value for the test centre.
1. `TestCentreCategory`: The test centre category where the sample was
   obtained (see above).
1. `Gender`: Either "F" (female), "M" (male), or "U" (unknown)
1. `Onset`: Either null or a "YYYY-MM-DD" date of symptom onset.
1. `personHash`: A secure one-way hash value for the individual.
1. `PAMS1`: True or False, according to whether the first-positive RT-PCR of
   the person was done in a walk-in center.
1. `Hospitalized`: True or false, according to whether the subject was ever
   hospitalised when a positive RT-PCR was obtained.
1. `B117`: True or false, indicating whether an infection was from lineage
   B.1.1.7

Note that a person may have a series of leading negative tests. These are
included in the file for the purposes of making Table S1 in the paper,
which gives the detection rates for the various test centre categories.

## min-3-timeseries.json

Contains a [JSON](https://en.wikipedia.org/wiki/JSON) object with the data
from 4344 subjects who had RT-PCR tests on at least three different days
(with at least two tests being positive). The main data is "people", a list
of 4344 objects, each containing the following attributes:

* `personHash`: A secure one-way hash value for the individual.
* `gender`: Either "F" (female), "M" (male), or "U" (unknown)
* `PAMS1`: True or False, according to whether the first-positive RT-PCR of
  the person was done in a walk-in center.
* `hospitalized`: True or false, according to whether the subject was ever
  hospitalised when a positive RT-PCR was obtained.
* `onset`: Either null or a "YYYY-MM-DD" date of symptom onset.
* `B117`: True or false, indicating whether an infection was from lineage
  B.1.1.7

The following attributes are all lists, containing a value for each RT-PCR
result for the person:

* `viralLoad`: A floating point log10 viral load (copies / swab). A value of `0`
  indicates a negative test.
* `testName`: Either "LC480" or "T2" accordign to the RT-PCR system used.
* `date`: The "YYYY-MM-DD" date the sample arrived at the diagnostic
  facility.
* `testCentre`: A secure one-way hash value for the test centre.
* `testCentreCategory`: The test centre category where the sample was
  obtained (see below).
* `age`: The age of the subject on the day the RT-PCR was done. These are
  rounded to one decimal place to help ensure anonymity. In the paper we
  used the full-precision ages.

### Test centre categories

As in Table S1 in the paper, test centre category abbreviations are as follows:

    ?: Unknown
    AIR: Airport
    C19: COVID-19 testing centre
    CP: Company physician
    ED: Emergency department
    FM: Forensic medicine
    H: Hospital
    ICU: Intensive care unit
    IDW: Infectious diseases ward
    L: Labor
    LW: Labour ward
    OD: Outpatient department
    PHD: Public health department
    PRI: Prison
    RES: Aged residence
    SM: Sports medicine
    WD: Ward

## Culture_probability_data_B.1.1.7.xlsx

An Excel spreadsheet with data from the cell culturing trials. The three
columns should be self-explanatory.

## Culture_probability_data_wild_type.xlsx

An Excel spreadsheet with culture isolation data from the Ranawaka et
al. paper
[SARS-CoV-2 Virus Culture and Subgenomic RNA for Respiratory Specimens from Patients with Mild Coronavirus Disease](https://wwwnc.cdc.gov/eid/article/26/11/20-3219_article)
and the van Kampen et al. paper
[Duration and key determinants of infectious virus shedding in hospitalized patients with coronavirus disease-2019 (COVID-19)](https://www.nature.com/articles/s41467-020-20568-4). Columns should be self-explanatory.

## Culture_probability_data_wild_type_woelfel.xlsx

An Excel spreadsheet with culture isolation data from the Wölfel et
al. paper
[Virological assessment of hospitalized patients with COVID-2019](https://www.nature.com/articles/s41586-020-2196-x).
Columns should be self-explanatory.

# Code

The following files are available.

## ExtendedMethods.Rmd

Contains [R Markdown](https://rmarkdown.rstudio.com/) with a description
(and the [R](https://www.r-project.org/) code) for the statistical
analysis.

The R Markdown has been processed into HTML (in
[ExtendedMethods.html](ExtendedMethods.html)). Click
[here](http://htmlpreview.github.io/?https://github.com/VirologyCharite/SARS-CoV-2-VL-paper/blob/main/ExtendedMethods.html)
to view the file in your browser.

## utils.R

Contains various [R](https://www.r-project.org/) utility functions.
