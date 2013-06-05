miccai2013
==========

MICCAI 2013 code - Segmenting Multiple Overlapping Cervical Cells by Joint Level Set.


====================
   What's inside?
====================
1. Dataset: four real EDF Pap smear images and synthetic test/training images in the folder './ims'.
2. Matlab code for the paper.

====================
   Before running
====================
1. This code requires third part libraries for GMM density estimation, Quick Shift and LBP. You need to download them separately.
    1.1 Statistical Pattern Recognition Toolbox: http://cmp.felk.cvut.cz/cmp/software/stprtool/
    1.2 VLFeat: http://www.vlfeat.org/
    1.3 LBP: http://www.cse.oulu.fi/CMV/Downloads/LBPMatlab (download the file "lbp.m" and "getmapping.m")
    1.4 
2. This code is compatible with Matlab R2012b or after.


====================
    How to run?
====================
1. Invoke script "Run.m". The argument in line 7 to invoke the function Runner_inOne(...) refers to:
    1.1 'EDF' - the real Pap smear images
    1.2 'Test' - the synthetic test images
    1.3 'Train' - the synthetic training images
