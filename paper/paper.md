---
title: 'pyIMLi-Heat3D: An extension for pyGIMLi to calculate thermal fields in 3D'
tags:
  - Python
  - geophysics
  - steady-state-heat flow
  - thermal calculations
  - Heat Flow
authors:
  - name: Georg-Maximilian GHüttner
    orcid: 0009-0003-3644-4818
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Johannes Gehrig
    affiliation: "2"
  - name: Wolfgang Szwillus
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: "2"
  - name: Jörg Ebbing
    affiliation: "2"
affiliations:
 - name: Institute of Geosciences, Kiel University, Germany
   index: 1
 - name: AWI, Germany 
   index: 2
date: 00 December 2025
bibliography: paper.bib

---

# Summary
The pyIMLi-Heat3D extension allows to calculate the thermal field and geothermal heat flow solving the steady-state-heat equation in 3D.


# Statement of need

Subglacial heat flux is an important factor for understanding and predicting ice-sheet evolution.
Calculating geothermal heat flux from the solid Earth can be done using the steady-state-heat equation, which requires information of the layering and properties of the crust and mantle (REF). This is often done in 1D (e.g. Lösing et al. 2000) and empiricial corrections exists (Colgan et al. 2021) to account for the influence of subglacial topography on geothermal heat flux. Such empirical relations do not exist for other paranmeters affection heat flux, e.g. lateral variations in thermal conductivity or radiogenic heat production.
Instead of using an empirical correction, we present a code to calculate directly heat flux in 3D from a geophysical model, where we can consider both variations in the geometry of the layers and the thermal properties based on finite elements. Such a tool is for example needed to link solid Earth thermal models to glaciological models (e.g. Wolovick et al. 202?=) or ice sheet models(e.g. ISSM, Larour et al. 2012) 

- 1d vs general heat func

# Design

- build upon pyGIMLI for mesh generation and finite element solver
- varibable input number of diff lithosphgeric layers and hp distribution in correct coords

# Usage

- load data and set thermal parameters
- let funcs do their job -> see exmaple
- look at ouput

# Acknowledgements

# References

