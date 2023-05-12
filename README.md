<!-- badges: start -->
[![Travis build status](https://travis-ci.org/bomeara/phydo.svg?branch=master)](https://travis-ci.org/bomeara/phydo)
<!-- badges: end -->

# phydo
Project to quickly provide an overview of a clade, as a starting point for human-driven analysis.

This is inspired by projects such as Rod Page's iSpecies or the Sanderson et al.'s phylota, as well as by the issues we frequently see in students first learning about phylogenetics.

If it works, someone can just ask for a report on a clade, and it will give info on diversity, location, phylogeny, conservation status, research on it, funding for its study, and conclusions from models (is it diversifying particularly fast, what traits correlate with each other, what is its climate envelope, etc.). Basically, automatically write a (bad, dull) article on a clade's evolution and ecology.

This is for our lab's latest hackathon on May 29, 2019.

Basic structure:

## Getting data

Data sources:

* gbif
* genbank
* paleobiodb
* EOL trait bank
* Open Tree of Life
* Datelife
* endangered species status
* NSF (funding)

## Processing data

* cleaning data
* converting into more usable forms (latitude and longitude -> biogeographic regions)

## Doing analyses

* diversification models: tree only
* diversification models: tree plus trait
* correlation: continuous traits
* univariate models for continuous and discrete (phylo signal, OU/BM, etc.)
* biogeographic models


## Creating report

Make a sensible, readable summary. Can we create a bot that creates basic papers on a clade? A web page with tables, figures, and downloadable data?
