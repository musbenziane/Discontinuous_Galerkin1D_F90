subroutine compute_analytical(nt,dt,esrc,gsrc,ircv,src,igreen,ne,N,xgll,c,asigma)
implicit none
integer                                             :: nt, i, N, esrc, gsrc, ircv, igreen, ne
real(kind=8)                                        :: xgll(ne,N+1), c, dt
real(kind=8), dimension(nt)                         :: greenfunc, src, asigma

do i=1,nt
    if (i*dt-ABS(xgll(esrc,gsrc)-xgll(igreen*ircv,CEILING(REAL((N+1)/2))))/c >=0) then
        greenfunc(i) = 1./(2*c) 
    else
        greenfunc(i)= 0
    end if
end do


call conv_fft(greenfunc,src,asigma)



end subroutine compute_analytical