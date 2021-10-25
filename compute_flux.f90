subroutine compute_flux(ne,u,N,Al,Ar,flux)
    integer, intent(in)          :: ne, N
    real(kind=8), intent(in)     :: u(ne,N+1,2), Al(ne,2,2), Ar(ne,2,2)
    real(kind=8), intent(out)    :: flux(ne,N+1,2)

    integer                      :: i

    !##########################################
    !#### Inside the computational domain  ####
    !##########################################



    !##########################################
    !#####  At the  domain's boundaries  ######
    !##########################################

end subroutine compute_flux

