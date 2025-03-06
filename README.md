# smartSim:  Simulation of splice aware single cell Smart-seq3 data

## Prerequisites

In order to use this package: you need folowing software:
* conda (https://www.anaconda.com/)

We provide a conda environment, which can be created with:
```bash
conda env create -f smartSim_env.yml
```
## Input

* reference gtf file
* reference fasta file
* aligned and barcoded Smart-seq3 data similar to the data you want to simulate (alternatively you can use [example_data/aligned](example_data/aligned))

## Usage
We provide a tutorial to simulate reads using the example data: [tutorial.R](tutorial.R)

