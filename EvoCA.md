# EvoCA

initial prompt for Claude

## Introduction

This is code for an evolutionary cellular automata model.  It follows on two articles on the GeneLife model, found in Docs (which read).  In Genelife, a binary 2d cellular automaton is has cells containing 1 associated with a genome governing the local dynamics.  The result is a spatially inhomogeneous CA, with local rules that evolve.

The binary CA variable will be denoted $v_{\bm{x}}$

Evoca follows on framework of GeneLife idea, but with the following changes:

- Instead of genomes being associated with 1's on the lattice, every cell is 'alive', with a genome determining its local rule.  The genome is present and operative for both 0 and 1 local states.

- The local rule specification for each cell will be part of that cell's genome.  The rule space will be a generalization of the quarter totalistic rules discussed in the geneliffe papers.  Neighborhood size is 5x5. The rule table is a lookup table with entries for each possible value of the weighted sum $S_{\bm{x}}=\sum_{\bm{\delta}} v_{\bm{x}+\bm{\delta}}\lVert \bm{\delta} \rVert$.  So the local rule LUT at $\bm{x}$ will be a map from the sum to $\{0,1\}$, $R_{\bm{x}}(S_{\bm{x}})$
    - Question:  how many distinct values does the sum take for the 5x5 neighborhood?  This will determine the size of the LUT.  If the size is prohibitively large, we must reduce neighborhood size down from full 5x5.

- There is a resource field that we will call food, $F(\bm{x})$:  the food lattice the same size as the CA lattice.  $F(\bm{x})$ denotes food found at each spatial location, and is a float with values between zero and one.  Food is affected 
    - by a uniform regeneration, adding an increment food_inc to every cell, at every time step.  `float food_inc` is a global metaparameter.
    - by each living cells:  a living cell at location $(\bm{x})$ will eat food, transfering a fraction of the food present on the food lattice at $(\bm{x})$ to its private store of food, $f(\bm{x})$.  The fraction is determined by a function of the local configuration around the cells location, $C(\bm{x}) = \{(\bm{x}) + \bm{\delta}\}$, with $\bm{\delta}$ ranging over 25 nearest neighbors of the surrounding 5x5 neighborhood (including $(0,0)$).  The local configuration is compared with a fiducial local configuration private to each living cell $E(\bm{x})$, and part of the living cell's genome, evolvable upon reproduction of the organism.  We will begin with $c(\bm{x})$ patterns being symmetric with respect to reflections about horizontal and vertical midlines, and about the diagonals.  Eventually we may want to generalize to other symmetries or random fiducial configurations.
        - Question:  How many bits are needed to specify $c(\bm{x})$ with these symmetries; i.e. how many independent $c(\bm{x})$ are there?
    - Thus the amount of food transferred from $F(\bm{x})$ to $f(\bm{x})$, i.e. the cells mouthful, is $$
    M(\bm{x}) = \frac{m}{25} \sum_{\substack{a \in C(\bm{x})\\ b \in c(\bm{x})}} \delta(a,b) $$
    with $\delta(a,b) = 1$ if $a=b$, $0$ otherwise.  $F(\bm{x})$ is decremented by $M(\bm{x})$, $f(\bm{x})$ is incremented by $M(\bm{x})$. $m$ is a global metaparameter to scale the mouthful.

- Reproduction with will occur when $f(\bm{x})$ reaches a threshold,`float food_repro`, another global metaparameter.  The reproducing cell will look at all nearest neighbors (Moore neighborhood) and replace the cell with the lowest $f(\bm{x}')$, where $\bm{x}' \in \text{nbrhd}(\bm{x})$.  'Replacement' means copying all genetic information, including $r(\bm{x})$, $c(\bm{x})$.  Food of the parent cell is set to set to $f(\bm{x})/2$, and the food of the child cell is set to $f(\bm{x})/2$  the genetic information 

## Code layout

We want to have the CA lattice with genome information in C++ in the C directory, with CA updates implemented in C for speed.  All the C variables must be mapped to python, which will handle
- iteration of a time step, calling a C function to iterate across the entire CA lattice implementing each different local rule.
- fast display in a separate SDL2 window.
- display of widgets to display and set global metaparams.
