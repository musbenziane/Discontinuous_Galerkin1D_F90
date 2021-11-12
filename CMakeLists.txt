cmake_minimum_required(VERSION 3.21)
project(DG1D Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS         "-fopenmp")
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native")
endif()


enable_language(Fortran)
add_executable(DG1D main.f90 readmodelfiles1D.f90 lagrangeprime.f90 legendrep.f90 ricker.f90 lagrangep.f90 shapefunc.f90
        connectivity_matrix.f90 compute_flux.f90 gll_library.f90 )



