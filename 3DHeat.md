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
    orcid: 0000-0000-0000-0000
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
date: 00 August 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary
The pyIMLi-Heat3D extension allows to calculate the thermal field and geothermal heat flow solving the steady-state-heat equation in 3D.


# Statement of need


Subglacial heat flux is an important factor for understanding and predicting ice-sheet evolution.
Calculating geothermal heat flux from the solid Earth can be done using the steady-state-heat equation, which requires information of the layering and properties of the crust and mantle (REF). This is often done in 1D (e.g. Lösing et al. 2000) and empiricial corrections exists (Colgan et al. 2021) to account for the influence of subglacial topography on geothermal heat flux. Such empirical relations do not exist for other paranmeters affection heat flux, e.g. lateral variations in thermal conductivity or radiogenic heat production.
Instead of using an empirical correction, we present a code to calculate directly heat flux in 3D from a geophysical model, where we can consider both variations in the geometry of the layers and the thermal properties based on finite elements. Such a tool is for example needed to link solid Earth thermal models to glaciological models (e.g. Wolovick et al. 202?=) or ice sheet models(e.g. ISSM, Larour et al. 2012) 
# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge ....

# References