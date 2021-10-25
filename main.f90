program DG1D
implicit none

    real(kind=8)                                 :: xmax, h, f0, dt, CFL, mindist, Ja, Jai
    real(kind=8), dimension(:), allocatable      :: xi, wi, v1D, rho1D, src, v1Dgll, rho1Dgll, mu,Z
    real(kind=8), dimension(:,:), allocatable    :: lprime, xgll, Minv, Ke
    real(kind=8), dimension(:,:,:), allocatable  :: Al, Ar, u, unew, k1, k2, flux
    integer                                      :: N, ne, nt, esrc, gsrc, isnap, ngll, i, j, it, el
    integer, dimension(:,:), allocatable         :: Cij
    character(len=40)                            :: filename, outname, filecheck


    filename="parameters.in"
    outname = "snapshots.bin"

    write(*,*) "##########################################"
    write(*,*) "######## Reading parameters file #########"
    write(*,*) "##########################################"

    print*,"Is the parameters input file (parameters.in) [Yes/no]"
    read(*,*) filecheck

    if (filecheck=="Yes" .or. filecheck=="yes" .or. filecheck=="y" .or. &
            filecheck=="Y") then
        write(*,*) "Reading simulation parameters..."

    elseif  (filecheck=="No" .or. filecheck=="no" .or. filecheck=="n" .or. &
            filecheck=="N") then
        write(*,*) "Enter simulation parameters text file name with extension"
        write(*,*) "40 characters max"
        read(*,*) filename

    else
        write(*,*) "Only: Yes/yes/Y/y & No/no/N/n are handled"
        write(*,*) "The program have been terminated, please star over"
        stop
    end if

    open (2, file=filename, status = 'old')
    read(2,*) xmax
    read(2,*) N
    read(2,*) ne
    read(2,*) h
    read(2,*) f0
    read(2,*) dt
    read(2,*) nt
    read(2,*) esrc
    read(2,*) gsrc
    read(2,*) isnap
    close(2)

    print*,"Maximum distance         -> ",xmax
    print*,"Polynomial order         -> ",N
    print*,"Number of elements       -> ",ne
    print*,"Element size             -> ",h
    print*,"Wavelet's peak frequency -> ",f0
    print*,"Time step                -> ",dt
    print*,"Number of time steps     -> ",nt
    print*,"Source location          -> ",esrc, gsrc
    print*,"Snapshot interval        -> ",isnap

    ngll = (N + 1) * ne
    Ja   = h / 2.
    Jai  = 1 / ja

    !##########################################
    !#####    Matrices allocation         #####
    !##########################################

    allocate(xi(N+1))                   ! GLL points
    allocate(wi(N+1))                   ! GLL Quadrature weights
    allocate(v1D(ne))                   ! 1D velocity model in elements
    allocate(rho1D(ne))                 ! Density velocity model in elements
    allocate(rho1Dgll(ngll))            ! 1D density model mapped
    allocate(v1Dgll(ngll))              ! 1D velocity mapped
    allocate(mu(ne))                  ! Shear modulus mapped
    allocate(Z(ne))                     ! Impepdences
    allocate(xgll(ne,N+1))              ! Array for global mapping
    allocate(lprime(N+1,N+1))           ! Dervatives of Lagrange polynomials
    allocate(src(nt))                   ! Source time function
    allocate(Cij(N+1,ne))               ! Connectivity matrix
    allocate(Minv(N+1,N+1))             ! Elemental mass matrix
    allocate(Ke(N+1,N+1))               ! Elemental stifness matrix
    allocate(Ar(ne,2,2),Al(ne,2,2))
    allocate(u(ne,N+1,2),unew(ne,N+1,2))! Solution fields
    allocate(k1(ne,N+1,2),k2(ne,N+1,2))
    allocate(flux(ne,N+1,2))


    !##########################################
    !#####      Calling subroutines       #####
    !##########################################

    call gll(N,xi,wi)                              ! Getting GLL points and weights
    call readmodelfiles1D(v1D, rho1D, ne)          ! Reading model files
    call connectivity_matrix(N,ne,Cij)             ! Getting connectivity matrix
    call shapefunc(N,h,ne, xgll)                   ! Global domain mapping
    !call mapmodel(N,ne, rho1D,v1D,rho1Dgll,v1Dgll)! Mapping models
    call lagrangeprime(N,lprime)                   ! Lagrange polynomials derivatives
    call ricker(nt,f0,dt,src)                      ! Source time function

    mu = (v1D**2) * rho1D

    write(*,*)"##########################################"
    write(*,*)"############### CFL Check ################"

    mindist = xgll(2,1) - xgll(1,1)
    CFL = (dt / mindist) * maxval(v1D(:))
    if (CFL > .8) then
        print"( a14,f6.3)"," Courant number is ",CFL
        print*,"Decrease time step, the program has been terminated"
        stop
    else
        print"(a14,f6.3)","Courant number is ",CFL
        print*,"Simulation is stable"
    end if
    write(*,*)"##########################################"

    !##########################################
    !####### Construct the mass matrix ########
    !##########################################

    do i=1,N+1
        Minv(i,i) = 1. / (wi(i) * Ja)
    end do

    !##########################################
    !##### Construct the stiffness matrix #####
    !##########################################

    do i=1,N+1
        do j=1,N+1
               Ke(i,j) = wi(j) * lprime(i,j)
        end do
    end do

    !##########################################
    !##### Construct the Flux matrices    #####
    !##########################################

    Z = rho1D * v1D

    do i=1,ne-2
        Ar(i,1,1) =  .5 * v1D(i);
        Ar(i,1,2) = -.5 * Z(i) * v1D(i);
        Ar(i,2,1) = -.5 * v1D(i) / Z(i);
        Ar(i,2,2) =  .5 * v1D(i);

        Al(i,1,1) = -.5 * v1D(i);
        Al(i,1,2) = -.5 * Z(i) * v1D(i);
        Al(i,2,1) = -.5 * v1D(i) / Z(i);
        Al(i,2,2) = -.5 * v1D(i);
    end do
    flux(1,1,:)     = MATMUL(RESHAPE(Ar(i,:,:),(/2,2/)),  RESHAPE((-u(i-1,N+1,:)),(/2/))) + &
                      MATMUL(RESHAPE(Al(i,:,:),(/2,2/)),  RESHAPE((-u(i  ,1,:)),(/2/)))

    write(*,*) "##########################################"
    write(*,*) "########### Begin time  loop  ############"
    write(*,*) "##########################################"



do it=1,nt

    u(esrc,gsrc,1) = src(it)   ! source injection - displacement component

    call compute_flux(ne,u,N,Al,Ar,flux)

    do el=2,ne-1
        k1(el,:,1)   = matmul(Minv,(-mu(el) * matmul(Ke, u(el,:,2))) - flux(el,:,1))
        k1(el,:,2)   = matmul(Minv,(-1. / rho1D(el)) * matmul(Ke,u(el,:,1)) - flux(el,:,2))
    end do

    do el=2,ne-2
        unew(el,:,1) = dt * matmul(Minv,(-mu(el) * matmul(Ke,u(el,:,2))) - &
                       flux(el,:,1)) + u(el,:,1)
        unew(el,:,2) = dt * matmul(Minv,(-1. / rho1D(el)) * matmul(Ke,u(el,:,1)) - &
                flux(el,:,2)) + u(el,:,2)
    end do

    call compute_flux(ne,unew,N,Al,Ar,flux)

    do el=2,ne-2
        k2(el,:,1)   = matmul(Minv,(-mu(el) * matmul(Ke, unew(el,:,2))) - flux(el,:,1))
        k2(el,:,2)   = matmul(Minv,(-1. / rho1D(el)) * matmul(Ke,unew(el,:,1)) - flux(el,:,2))
    end do

    unew = u + .5 * dt * (k1 + k2)

    u    = unew;

end do



end program
