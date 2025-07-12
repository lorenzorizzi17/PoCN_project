# PoCN_project


This repository contains the two projects completed as the final exam for the *Physics of Complex Networks* course at the University of Padua:
- Project n.37 - Explosive percolation (Score: 0.5)
- Project n.45 - Constructing a graph of Europe's railway infrastructure (Score: 1.0)

## Project n.37 - Explosive percolation
This project focuses on the simulation and analysis of Achlioptas processes, which lead to *explosive percolation*, on various network topologies (ER-like, scale-free). By implementing these processes, we study how sudden transitions in connectivity emerge, contrasting them with classical percolation models.

#### Implementation details
The simulation is implemented in C++ and builds a percolation algorithm based on the Newman-Ziff method (which uses pointers, thus the choice of C++). Code is divided as follows:
- Folder `gen\`: Generative part of the algorithm
- Folder `proc\`: Data analysis performed using ROOT software (however, any data analysis oriented language will do the job)
Inside `gen\`:
- `include\`: header files (automatically included)
- `lib\`: a lib-like C++ file that needs to be linked with the main.cpp
- Four different "main" C++ files. Each one of those, when linked with `lib\graph.cpp` will generate data depending on its specific purpose

To boost further performances, OpenMP framework was used to achieve SPMD parallelism (so, when compiling, one should also link against OpenMP libraries)

Output data include:
- LCC size/ $\chi$ with respect to m
- Finite size cluster distribution around criticality
- Scaling of transition window

## Project n.45 - Railways in Europe
This project builds a railway network graph starting from EUROGEOMAP geospatial data, enabling a network-like representation of railway connectivity in Europe

#### Implementation details
The task was performed in R. In particular, a RStudio Markdown is available. Alternatively, two plain R scripts are also provided:
- `EURailNetwork.R`: Builds the railroad network for a specific country, given its country code. Outputs (in `data/`) two files for each country, one containing information on nodes-stations and one containing network connectivity information. It also creates a pdf file representing the graph
- `EURWholeEurope.R`: Similar as above, but considering data from all Europe