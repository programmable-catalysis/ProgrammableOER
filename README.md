# Programmable OER

## Description
This repository contains the code used to generate the data in Dynamic Promotion of the Oxygen Evolution Reaction via Programmable Metal Oxides, [DOI: 10.26434/chemrxiv-2024-gs6zn](https://chemrxiv.org/engage/chemrxiv/article-details/65af381d66c13817290d5404)

If you use the code or data contained in this repository, please cite the corresponding preprint using the following citation:
```
@article{Gathmann2024,
author = {Sallye R. Gathmann and Christopher J. Bartel and Lars C. Grabow and Omar A. Abdelrahman and C. Daniel Frisbie and Paul J. Dauenhauer},
doi = {10.26434/chemrxiv-2024-gs6zn},
journal = {ChemRxiv},
title = {Dynamic Promotion of the Oxygen Evolution Reaction via Programmable Metal Oxides},
year = {2024}
}
```



## Getting Started
This code was written in Julia 1.8.4 and has not been tested for compatability with newer versions of Julia. If you run into trouble using the current Julia release, v1.8.4 can be donwloaded from the [old releases page](https://julialang.org/downloads/oldreleases/). 

Required packages: CSV, DataFrames, DifferentialEquations, NBInclude, PyPlot, Trapz (for DRC only: DiffEqSensitivity, ForwardDiff)
