---
title: 'Musubi: Octree based Lattice-Boltzmann solver with support for multiple species'
tags:
  - Fortran
  - fluid dynamics
  - Lattice-Boltzmann method
  - Maxwell-Stefan equations
  - particle transport
  - distributed parallelism
authors:
  - name: Kannan Masilamani
    orcid: 0000-0002-3640-2154
    equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Kartik Jain
    orcid: 0000-0002-6540-9304
    equal-contrib: true
    affiliation: 2
  - name: Harald Klimach
    orcid: 0000-0002-6054-5681
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 3
affiliations:
 - name: German Aerospace Center (DLR), Germany
   index: 1
 - name: University of Twente, Netherlands
   index: 2
date: 26 November 2025
bibliography: paper.bib
---

# Summary

Musubi is a multi-level, parallel lattice Boltzmann solver and part of the APES suite.
It is working on an octree mesh that is linearized by a (Morton) space-filling curve and
uses efficient data structures allowing adaptive parallel simulations.

Musubi offers several collision kernels (BGK, MRT, HRR, Cumulant) and is designed to deal
with huge meshes (billions of lattices) efficiently.

It is written in Fortran, with language constructs from Fortran 2003.

# Statement of need

Highly resolved fluid simulations are an integral part in many scientifc application areas.
Due to the nonlinearity and the large amount of degrees of freedom to consider in these
problems, these simulations require significant computational resources, which typically
are only available in distributed parallel systems.
Musubi implements the lattice Boltzmann method (LBM) with a Message Passing Interface
(MPI) parallelization with a fully distributed handling of the mesh data, avoiding
bottlenecks on individual processors and enabling the scaling of the simulation to
hundreds of thousands of MPI processes.

# Acknowlegements

# References
