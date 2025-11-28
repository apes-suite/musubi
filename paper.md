---
title: 'Musubi: Octree based Lattice-Boltzmann solver for multi-pyhsics'
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
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Kartik Jain
    orcid: 0000-0002-6540-9304
    affiliation: 2
  - name: Manuel Hasert
    affiliation: 3
  - name: Harald Klimach
    orcid: 0000-0002-6054-5681
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
affiliations:
 - name: German Aerospace Center, DLR, Germany
   index: 1
   ror: 04bwf3e34
 - name: University of Twente, Netherlands
   index: 2
   ror: 006hf6230
 - name: Festo SE & Co. KG
   index: 3
   ror: 03ga8q162
date: 28 November 2025
bibliography: paper.bib
---

# Summary

Musubi is a multi-level, parallel lattice Boltzmann solver and part of the APES suite.
It is working on an octree mesh that is linearized by a (Morton) space-filling curve and
uses efficient data structures allowing adaptive parallel simulations.

Musubi is designed to deal with huge meshes (billions of lattices) and complex geometries
on large computing systems efficiently.
It can be used for a wide range of application areas from electrodialysis [@Masilamani:2020]
over biomedical problems [@Jain:2016] and aero-dynamic setups [@Spinelli:2024] to aero-acoustic
simulations [@Hasert:2013].

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

# The lattice Boltzmann method

The lattice Boltzmann method utilizes ideas of cellular automata and represents at its
core a basic two step algorithm.
The state of the fluid is represented by particle density functions (PDF) in a discrete
velocity field.
These PDFs reside on the lattices and are exchanged along the discrete velocity directions.
The two steps of the algorithm are now the streaming of the PDF information along
the velocity directions, followed by the so-called collision, computing the new
PDF on each lattice.
This modeling with discrete velocities also allows for a straight forward handling of
complicated wall boundaries, as a simple line intersection with the wall geometry
can be used to accurately model the surface.
Due to these properties the method has gained popularity in the field of computational
fluid dynamics over the last decades.
Other Open Source solvers that utilize this method are for example Palabos [@Palabos2020],
OpenLB [@olbPaper2021] and waLBerla [@BAUER2021478].

# The Musubi implementation

Musubi implements the lattice Boltzmann method in the form of kernels that can be
run on individual refinement levels of an octree mesh.
It is developed within the APES-Suite [@Klimach:2014] of simulation tools based on the central
Treelm library [@Klimach:2012vi] that provides the handling of this octree mesh on distributed
parallel systems.
The interpolation and transformation between the involved levels for the local
refinement are separated from the kernel, allowing for an implementation of the
respective methods without encumberment by the interpolation between the different
resolutions.
This method was described in detail in [@hasert:2013jc].
There are various collision schemes implemented (BGK, MRT, HRR, Cumulants) [@Spinelli:2023],
which can be used on a range for stencil configurations (discrete velocity directions).
It is also possible to consider the transport of particles and passive scalars in
the flow.

# Acknowlegements

This software has been written by many people over the years.
The individual authors can be found in each file with the respective copyright
statement.
Not appearing in the list of authors is Sabine Roller, who enabled the development
of this software in the first place and we are very grateful for this possibility.
We especially thank our fellow contributors to this code basis Jiaxing Qi [@Qi:2017],
Jens Zudrop [@Zudrop:2015], Simon Zimny [@Zimny:2015], Peter Vitt, Jana Gericke,
Tristan Vlogman [@Vlogman:2025], Mengyu Wang and many students.
The development of Musubi was partially funded by the German Federal Ministry of Education and Research
(Bundesministerium für Bildung und Forschung, BMBF) in the framework of the HPC software initiative in
the project HISEEM, by the European Commission in the Seventh Framework Programme in the area of
Virtual Physiological Human (THROMBUS project, ICT-2009.5.3, project reference 269966), by the
German Research School of Simulation Sciences, the University of Siegen, the University of Twente
and the German Aerospace Center, DLR.
We are grateful for the computing time provided by LRZ in Munich and by HLRS in Stuttgart who
also contributed performance evaluations within the POP project.

# References
