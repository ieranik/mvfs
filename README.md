# Overview
This repository contains the source code for the paper titled [The Maximum Visibility Facility Selection Query in Spatial Databases](https://drive.google.com/file/d/1DYdF3Uk-N_qtIRDsM4yWJiSLJB7zJQM9/view?usp=sharing). This paper was published in _ACM SIGSPATIAL_ back in 2019. This paper proposes exact and approximation algorithms for the _Maximum Visibility Facility Selection_ (MVFS) problem. The objective of the MVFS problem is to select $k$ out of $n$ given locations to place cameras in an occluded space, such that the visibility coverage of the free space provided by the $k$ cameras is maximized. 

In this paper, we point out the limitations of the existing solutions to the MVFS problem and devise clever techniques to overcome the limitations. First, most of the existing scalable solutions assume a discretized data space, i.e., the data-space is sampled to form a set of control points, and the visibility coverage is measured in terms of the number of visible control points instead of the actual area/volume of the visible region. To this end, we propose an exact in-memory algorithm that does not require the data-space to be discretized into control points. The key idea of our approach is to construct a novel visibility-based triangulation of the data-space, namely equivisibility triangulation, using which we accurately determine the area of the region visible from one or more facility locations. 

Second, existing exact solutions to the MVFS problem involve solving BIP formulations, and the number of binary variables in this formulation is high. As a result, the performance of these exact algorithms is not satisfactory even for small instances. In this paper, we model the MVFS problem as a graph problem, and employ a divide-and-conquer strategy to obtain the optimum solution. We determine a vertex separator (a subset of the vertices, removal of which along with the associated edges renders the graph disconnected) of the graph, recursively solve the smaller subgraphs split by the vertex separator, and finally merge the solutions obtained from the subgraphs to form the overall solution to the MVFS problem. Our proposed algorithm can solve larger instances of the MVFS problem in comparison with the existing exact solutions.

Finally, no existing approach offers a disk-based solution and thus cannot handle the scenarios where the obstacle set is too large to be stored in RAM. Thus existing solutions are not applicable in the context of big data. In this work, we propose the first external-memory algorithm for the MVFS problem. To attain scalability, we apply a greedy approximate technique to reduce the computational cost, use a disk-resident spatial data structure to accelerate obstacle retrieval, and employ a heuristic-driven best-first search technique to reduce the I/O overhead.

# Directory Tree 

* `datasets`: We use both synthetic and real-world data for the experiments. Synthetic data for the MVFS problem can be generated using the code in `dataset/gen_dataset.cpp`. We also conduct experiments to test our proposed algorithms in a real-world setting using the [Boston 3D dataset](https://www.bostonplans.org/3d-data-maps).

* `algorithms`
  * `exact.cpp`: This file contains two versions of our proposed exact algorithm for the MVFS problem: _basicExact_ and _efficientExact_. In the basic exact algorithm, we construct the equivisibility triangulation of the space, solve each connected component independently, and merge the results of all the components to obtain the solution to the MVFS problem. In the efficient exact algorithm, we accelerate the process of solving each component by decomposing it using vertex separators.   
  * `greedy.cpp`: This file contains my implementation of our proposed greedy algorithm, which, at each step, chooses the facility that maximizes the sum of the areas of the uncovered triangles. This greedy algorithm has been proved as the best-possible polynomial time approximation algorithm for the MVFS problem and has an approximation ratio of $1 - \frac{1}{e}$.
  * `scalable_greedy.cpp`: This file contains the implementation of our proposed greedy algorithm that achieves scalability by discretizing the data space and indexing the obstacles in a spatial data structure called [R-Tree](https://en.wikipedia.org/wiki/R-tree). 



# Links

* [The Maximum Visibility Facility Selection Query in Spatial Databases](https://drive.google.com/file/d/1DYdF3Uk-N_qtIRDsM4yWJiSLJB7zJQM9/view?usp=sharing)
* [R-Tree](https://en.wikipedia.org/wiki/R-tree)
* [Boston 3D Dataset](https://www.bostonplans.org/3d-data-maps)



