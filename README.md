<h1> Disconitnuous Galerkin implementation in 1D for the Elastic Wave PDE in Cartesian structrured mesh </h1>
<br/>
Check Cmake version before compiling if using CMake <br/>

<pre>
    !####################################################################################################################
    ! DG1D: Noda Discontinuous Galerkin Method in 1D for the Elastic Wave Equation with Regular Mesh
    ! NOTE: This program is yet to be tested against analytical solutions.
    ! 
    ! Language: Fortran 90, with parralel impelementation using OpenMP API
    ! 
    ! Sources:  Igel 2017, Leveque 2002, Hesthaven & Warburton 2008
    ! 
    ! The code used for arbitrary GLL points and weights was created by M.I.T departement of engineering. Link is hereafter
    ! https://geodynamics.org/cig/doxygen/release/specfem3d/    | file name: gll_library.f90
    ! 
    ! This is part of the Numerical Modelling Workshop.
    ! 
    ! Supervisor: Pr. Emmanual Chaljub
    ! Author    : Mus Benziane
    !
    ! Input file: [Example]
    ! testing              ! model name prefix (model names: prefix_vp, prefix_rho)
    ! 4000                 ! xmax
    ! 6                    ! Polynomial order
    ! 800                  ! Number of elements
    ! 5                    ! Element size
    ! 3.                   ! Wavelet's peak frequency
    ! 0.000040             ! Time step
    ! 50000                ! Number of time steps
    ! 3                    ! Source location - element number
    ! 1                    ! Source location - collocation point
    ! 100                  ! Snapshot interval
    ! 3                    ! [1/2/3] 1: Free surface, 2: Rigid wall, 3: Periodic
    ! .6                   ! Att constant for sponge layer
    ! 
    ! -> Model files in C-Style binary floats [doubles]: Vs, Rho files are needed.
    !                                                  : For simple models, use create1Dmodel_files.f90
    ! 
    ! -> Outputs are created in OUTPUT/ if OUTPUT/ is not created by the user, the program will not handle it.
    !    Output files in OUTPUT directory:
    !
    !####################################################################################################################
</pre>
