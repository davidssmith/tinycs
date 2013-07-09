
# tinycs

*tinycs* is a minimal compressed sensing (CS) toolkit designed
to allow MR imaging scientists to design undersampled
acquisitions and reconstruct the resulting data with CS without
needing to be a CS expert.

Currently, TinyCS supports Cartesian geometries with a total
variation constrained reconstruction only.  If there is
sufficient interest, I can add additional acquisition designs
and sparsity constraints.

The Cartesian reconstruction is based on the split Bregman
code written by Tom Goldstein, originally available here:
<http://tag7.web.rice.edu/Split_Bregman.html>

