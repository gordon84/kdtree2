A refactored version of KDTREE2.

All necessary code is in the single Fortran file mod_kdtree2.f90 except for mod_types.f90 which specifies a real kind rt.

If you compile it with $FC mod_types.f90 mod_kdtree2.f90 test_kdtree2.f90 (or $FC \*.f90 a few times if you are lazy), the executable will run some basic examples.

License and original readme are at the top of mod_kdtree2.f90.

Original version is by Matthew B. Kennel, see https://arxiv.org/abs/physics/0408067
