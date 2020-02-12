# Code-for-BGSRP
A Fast Algorithm for Recovery of Bandlimited Graph Signals Based on the Reproducing Kernel Hilbert Space

Copyright (C) 2020 Qian Zhang and Lihua Yang. All rights reserved.
Notice:  This repository is free software for scientific research, you can redistribute it and/or modify it. But not for commercial use.
This repository is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability of fitness for a particular purpose.


This repository is used to open the codes for paper "A Fast Algorithm for Recovery of Bandlimited Graph Signals Based on the Reproducing Kernel Hilbert Space" for BGSRP to make it easier for readers to understand the proposed algorithms. 
 
In this repository, the "gsp_BGSRP_recon.m" is the Matlab code for Algorithm 1 in paper "A Fast Algorithm for Recovery...". In order to call this code and subsequent codes, you may need to use the GSPBOX toolbox. 

For a given graph G (a structure defined in GSPBOX), a graph signal f. Only a partial values of f are known, let x0 the locations of known labeled points, let y0 the labels of known labeled points. By utilizing the function "y = gsp_BGSRP_recon(G,x0,y0,param)", you 
can obtain a reconstructed grpah signal y. Here param is a struct for optional parameters. For more details, please download and read the gsp_BGSRP_recon.m file. 

The "gsp_BGSRP_recon_exfig4.m" is the code of experiment conducted for "Fig.4" in paper. This experiment is conducted on sensor graphs by varying the number of vertices in the range N = 800, 1200,..., 4000. In this experiment, we construct 50 bandlimited original signals for each graph. The known labeled vertices and the bandwidth coefficient are fixed as l=400,n=400. The computational costs, TIMEs, are depicted as a box diagram in Fig.4. For more details, please download and read the gsp_BGSRP_recon_exfig4.m file. 

 
 

When you use these codes, please kindly cite following references:
1    N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
     ArXiv e-prints, Aug. 2014.

2     Qian Zhang, Chao Huang, Zhihua Yang, and Lihua Yang. A Fast 
      Algorithm for Recovery of Bandlimited Graph Signals Based on the 
      Reproducing Kernel Hilbert Space, IEEE Transactions on Signal 
      Processing, 2020. 
      
      and if needed
3     J. Paratte and L. Martin. Fast eigenspace approximation using random
      signals. arXiv preprint arXiv:1611.00938, 2016.



