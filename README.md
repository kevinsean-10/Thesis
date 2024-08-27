This is the code repository for the Thesis (final project) by Kevin Sean Hans Lopulalan, completed during his studies at the Bandung Institute of Technology. The thesis, titled "Comparison Study of Methods for Finding All Roots of Nonlinear Equation Systems Using Differential Evolution," explores the performance differences among three optimization methods used to find solutions for multiple nonlinear equation systems. The findings of this research provide new insights into the characteristics of each method and how they should be utilized. If you are interested in reviewing the results of this thesis, please email at kevinsean.lopulalan@gmail.com.

In general, three methods were used:
1. **Clustering Differential Evolution (CDE):** This method is adapted from the root-finding approach discussed in Sidarto and Kania's (2015) paper, which uses Clustering Spiral Optimization. However, in this case, the optimization method employed is Differential Evolution.
2. **Hypercube Differential Evolution (HDE):** This method is adapted from the hypercube approach by Mastorakis (2005), which uses incremental search to find all roots in nonlinear equation systems.
3. **Repulsion-based Adaptive Differential Evolution (RADE):** This method, described by Gong et al. (2020), employs repulsion techniques, diversity preservation mechanisms, and adaptive techniques with DE to find solutions for multiple nonlinear equation systems.

There are three main folders in this repository. Each folder contains files named [xxx.m], [xxx_trial.m], and [xxx_simulation.m], where "xxx" represents one of CDE, HDE, or RADE. The [xxx.m] files are the main class files containing the algorithms. The [xxx_trial.m] files are used to test the algorithms in the [xxx.m] files. The [xxx_simulation.m] files are used to perform N simulations on the algorithm of the respective method.

If you are interested in viewing the results of this research, please email at kevinsean.lopulalan@gmail.com.

Thank you for your time.
