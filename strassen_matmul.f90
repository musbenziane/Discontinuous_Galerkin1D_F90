recursive subroutine strassen_multiply(AA, BB, CC, nn)
    USE OMP
    integer, intent(in)   :: nn
    real(8), intent(in)  :: AA(nn,nn), BB(nn,nn)
    real(8), intent(out) :: CC(nn,nn)
    real(8), dimension(nn/2,nn/2) :: A11, A21, A12, A22, B11, B21, B12, B22
    real(8), dimension(nn/2,nn/2) :: Q1, Q2, Q3, Q4, Q5, Q6, Q7
    integer                       :: threshold=128

    if(iand(nn,1) /= 0 .OR. nn < threshold) then
        CC = matmul(AA,BB)
    else
        A11 = AA(1:nn/2,1:nn/2)
        A21 = AA(nn/2+1:nn,1:nn/2)
        A12 = AA(1:nn/2,nn/2+1:nn)
        A22 = AA(nn/2+1:nn,nn/2+1:nn)
        B11 = BB(1:nn/2,1:nn/2)
        B21 = BB(nn/2+1:nn,1:nn/2)
        B12 = BB(1:nn/2,nn/2+1:nn)
        B22 = BB(nn/2+1:nn,nn/2+1:nn)
        !$OMP TASK  SHARED(A11,A22,B11,B22,Q1,nn) ! FINAL(nn<128) MERGEABLE
        call strassen_multiply(A11+A22, B11+B22, Q1, nn/2)
        !$OMP END TASK
        !$OMP TASK  SHARED(A21,A22,B11,Q2,nn) ! FINAL(nn<128) MERGEABLE
        call strassen_multiply(A21+A22, B11, Q2, nn/2)
        !$OMP END TASK
        !$OMP TASK  SHARED(A11,B12,B22,Q3,nn) ! FINAL(nn<128) MERGEABLE
        call strassen_multiply(A11, B12-B22, Q3, nn/2)
        !$OMP END TASK
        !$OMP TASK  SHARED(A22,B11,B21,Q4,nn) ! FINAL(nn<128) MERGEABLE
        call strassen_multiply(A22, -B11+B21, Q4, nn/2)
        !$OMP END TASK
        !$OMP TASK  SHARED(A11,A12,B22,Q5,nn) ! FINAL(nn<128) MERGEABLE
        call strassen_multiply(A11+A12, B22, Q5, nn/2)
        !$OMP END TASK
        !$OMP TASK  SHARED(A11,A21,B11,B12,Q6,nn) ! FINAL(nn<128) MERGEABLE
        call strassen_multiply(-A11+A21, B11+B12, Q6, nn/2)
        !$OMP END TASK
        !$OMP TASK  SHARED(A12,A22,B21,B22,Q7,nn) ! FINAL(nn<128) MERGEABLE
        call strassen_multiply(A12-A22, B21+B22, Q7, nn/2)
        !$OMP END TASK
        !$OMP TASKWAIT
        CC(1:nn/2,1:nn/2) = Q1+Q4-Q5+Q7
        CC(nn/2+1:nn,1:nn/2) = Q2+Q4
        CC(1:nn/2,nn/2+1:nn) = Q3+Q5
        CC(nn/2+1:nn,nn/2+1:nn) = Q1+Q3-Q2+Q6
    end if
    return
end subroutine strassen_multiply