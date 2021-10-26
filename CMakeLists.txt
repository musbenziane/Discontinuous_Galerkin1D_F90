cmake_minimum_required(VERSION 3.21)
project(DG1D Fortran)


enable_language(Fortran)
add_executable(DG1D main.f90 readmodelfiles1D.f90 lagrangeprime.f90 gll.f90 legendrep.f90 ricker.f90 lagrangep.f90 shapefunc.f90
        connectivity_matrix.f90 compute_flux.f90)



