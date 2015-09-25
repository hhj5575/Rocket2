
!-------------------------------------------------------------------------------------
include 'mkl_dss.f90'





subroutine cholesky(A, b, x, ndim, num_col_b)


! Solve symmetric psd Ax = b  using LAPACK based on Cholesky method
!    
! ndim: dimension of A
! LAPACK: DPOTRF, DPOTRS

integer, intent(in) :: ndim, num_col_b
real, intent(in) :: A(ndim,ndim), b(ndim,num_col_b) 
real, intent(out) :: x(ndim,num_col_b)
integer :: info


x = b  ! solution x must be separated from b

!call DPOTRF('U', ndim, A, ndim, info)
!if (info == 0) then
!  call DPOTRS('U', ndim, num_col_b, A, ndim, x, ndim, info)
!  if (info < 0) stop 'cholesky: LAPACK DPOTRS error- illegal argument value'
!else
!  x = 0.0
!endif

call DPOSV('U', ndim, num_col_b, A, ndim, x, ndim, info)	
if (info /= 0) stop 'cholesky: LAPACK DPOSV error- illegal argument value'  
  
end subroutine cholesky






subroutine general_lu(A, b, x, ndim, num_col_b)


! Solve general Ax = b  using LAPACK
!    
! ndim: dimension of A
! LAPACK: DPOTRF, DPOTRS

integer, intent(in) :: ndim, num_col_b
real, intent(in) :: A(ndim,ndim), b(ndim,num_col_b) 
real, intent(out) :: x(ndim,num_col_b)
integer :: info, ipiv(ndim)


x = b  ! solution x must be separated from b

!call DPOTRF('U', ndim, A, ndim, info)
!if (info == 0) then
!  call DPOTRS('U', ndim, num_col_b, A, ndim, x, ndim, info)
!  if (info < 0) stop 'cholesky: LAPACK DPOTRS error- illegal argument value'
!else
!  x = 0.0
!endif


call DGESV(ndim, num_col_b, A, ndim, ipiv, x, ndim, info)	
if (info /= 0) stop 'general_lu: LAPACK DGESV error- illegal argument value'  
  
end subroutine general_lu







subroutine triple_mul_sparse(m, n, ndim, col_info, row_info, sparse, matrix, result)

! R = R + S'AS
! sparse(S): sparse matrix (m-by-n) in CSC format
! matrix(A): a matrix (m-by-m)
! result matrix(R): n-by-n
! ndim: size(length) of row_info & sparse arrays

implicit none

integer,intent(in) :: m, n, ndim, col_info(n+1), row_info(ndim)
real,intent(in) :: sparse(ndim), matrix(m,m)
real,intent(inout) :: result(n,n)

integer :: j, k, start_pt, end_pt
real :: vector(m)


do k = 1, n
  start_pt = col_info(k)
  end_pt = col_info(k+1) - 1
  if (end_pt >= start_pt) then
    vector = MATMUL(matrix(:,row_info(start_pt:end_pt)), sparse(start_pt:end_pt))  
    do j = 1, n
      start_pt = col_info(j)
      end_pt = col_info(j+1) - 1
      if (end_pt >= start_pt) then
        result(j,k) = result(j,k) + DOT_PRODUCT(sparse(start_pt:end_pt),vector(row_info(start_pt:end_pt)))
      endif
    end do
  endif
end do  ! k

end subroutine triple_mul_sparse





subroutine mul_tr_sparse_matrix(m, l, n, ndim, col_info, row_info, sparse, matrix, result)

! R = R + S'A
! sparse(S): sparse matrix (m-by-l) in CSC format
! S': transpose of S (l-by-m)
! matrix(A): a matrix (m-by-n)
! result matrix(R): l-by-n
! ndim: size(length) of row_info & sparse arrays

implicit none

integer,intent(in) :: m, l, n, ndim, col_info(l+1), row_info(ndim)
real,intent(in) :: sparse(ndim), matrix(m,n)
real,intent(inout) :: result(l,n)

integer :: j, k, start_pt, end_pt


do k = 1, n
  do j = 1, l
    start_pt = col_info(j)
    end_pt = col_info(j+1) - 1
    if (end_pt >= start_pt) then
      result(j,k) = result(j,k) + DOT_PRODUCT(sparse(start_pt:end_pt),matrix(row_info(start_pt:end_pt),k))
    endif
  end do
end do  ! k


end subroutine mul_tr_sparse_matrix






subroutine mul_sparse_matrix(m, l, n, ndim, col_info, row_info, sparse, matrix, result)

! R = SA
! sparse(S): sparse matrix (m-by-l) in CSC format
! matrix(A): a matrix (l-by-n)
! result matrix(R): m-by-n
! ndim: size(length) of row_info & sparse arrays

implicit none

integer,intent(in) :: m, l, n, ndim, col_info(l+1), row_info(ndim)
real,intent(in) :: sparse(ndim), matrix(l,n)
real,intent(out) :: result(m,n)

integer :: i, j, k, start_pt, end_pt
real :: factor


result  = 0.0

do k = 1, n
  do i = 1, l
    factor = matrix(i,k)
    start_pt = col_info(i)
    end_pt = col_info(i+1) - 1
    if (end_pt >= start_pt) then
      result(row_info(start_pt:end_pt),k) = result(row_info(start_pt:end_pt),k) + factor*sparse(start_pt:end_pt)   
    endif
  end do
end do  ! k

end subroutine mul_sparse_matrix







subroutine mul_matrix_sparse(m, l, n, ndim, matrix, col_info, row_info, sparse, result)

! R = AS
! matrix(A): a matrix (m-by-l)
! sparse(S): sparse matrix (l-by-n) in CSC format
! result matrix(R): m-by-n
! ndim: size(length) of row_info & sparse arrays

implicit none

integer,intent(in) :: m, l, n, ndim, col_info(n+1), row_info(ndim)
real,intent(in) :: sparse(ndim), matrix(m,l)
real,intent(out) :: result(m,n)

integer :: i, j, k, start_pt, end_pt
real :: factor

result  = 0.0

do k = 1, n
  start_pt = col_info(k)
  end_pt = col_info(k+1) - 1
  if (end_pt >= start_pt) then
    result(:,k) = MATMUL(matrix(:,row_info(start_pt:end_pt)), sparse(start_pt:end_pt))   
  endif
end do  ! k

end subroutine mul_matrix_sparse






real function calc_length(vec1, vec2)

! Calculate length for 2D
! BLAS Library: DNRM2

implicit none

real, intent(in) :: vec1(2), vec2(2)
real :: x(2)
real :: DNRM2

x = vec1 - vec2 
calc_length = DNRM2(2, x, 1) 

end function calc_length






real function calc_length_n(vec1, vec2, n)

! Calculate length for n-dimension
! BLAS Library: DNRM2

implicit none

integer, intent(in) :: n
real, intent(in) :: vec1(n), vec2(n)
real :: x(n)
real :: DNRM2

x = vec1 - vec2 
calc_length_n = DNRM2(n, x, 1) 

end function calc_length_n






subroutine spect_3d(str, prst)

implicit none

real, intent(in) :: str(6)
real, intent(out) :: prst(3)
integer :: i,nrot
real:: pvect(3,3), a(3,3)


do i = 1,3
  a(i,i) = str(i)
enddo
a(1,2) = str(4)
a(1,3) = str(6)
a(2,3) = str(5)
a(2,1) = a(1,2)
a(3,1) = a(1,3)
a(3,2) = a(2,3)

call spect_jacobi(a,3,3,prst,pvect,nrot)
call spect_eigsrt(prst,pvect,3,3)

end subroutine




subroutine spect_2d(str, prst)

implicit none
   
real, intent(in) :: str(4)
real, intent(out) :: prst(2)
integer :: i,nrot
real:: pvect(2,2), a(2,2)


do i = 1,2
  a(i,i) = str(i)
enddo
a(1,2) = str(4)
a(2,1) = a(1,2)

call spect_jacobi(a,2,2,prst,pvect,nrot)
call spect_eigsrt(prst,pvect,2,2)

end subroutine





subroutine spect_JACOBI(A,N,NP,D,V,NROT)

integer :: IP,IQ,NROT,N,NP,I,J
real :: SM,TRESH,G,H,T,THETA,C,S,TAU
real :: A(NP,NP),D(NP),V(NP,NP),B(10),Z(10)

integer, parameter :: MAX = 50
real, parameter :: TOL = 1.e0-12
      
	  
      DO 12 IP=1,N
        DO 11 IQ=1,N
          V(IP,IQ)=0.d0
11      CONTINUE
        V(IP,IP)=1.d0
12    CONTINUE
      DO 13 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.d0
13    CONTINUE
      NROT=0
      DO 24 I=1,MAX
        SM=0.d0
        DO 15 IP=1,N-1
          DO 14 IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
14        CONTINUE
15      CONTINUE
        IF(SM.EQ.0.)RETURN
        IF(I.LT.4)THEN
          TRESH=0.2d0*SM/N**2
        ELSE
          TRESH=0.d0
        ENDIF
        DO 22 IP=1,N-1
          DO 21 IQ=IP+1,N
            G=100.*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP))).AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ)=0.d0
            ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(ABS(H)+G.EQ.ABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5d0*H/A(IP,IQ)
                T=1.d0/(ABS(THETA)+SQRT(1.d0+THETA**2))
                IF(THETA.LT.0.d0)T=-T
              ENDIF
              C=1.d0/SQRT(1.d0+T**2)
              S=T*C
              TAU=S/(1.d0+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=0.d0
              DO 16 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
16            CONTINUE
              DO 17 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
17            CONTINUE
              DO 18 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
18            CONTINUE
              DO 19 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
19            CONTINUE
              NROT=NROT+1
            ENDIF
21        CONTINUE
22      CONTINUE
        DO 23 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.d0
23      CONTINUE
24    CONTINUE
      STOP 'JACOBI: Maximum Iteration !'

end subroutine





subroutine spect_EIGSRT(D,V,N,NP)

integer :: N,NP,I,J,K
real :: D(NP),V(NP,NP),P
      
	  
      DO 13 I=1,N-1
        K=I
        P=D(I)
        DO 11 J=I+1,N
          IF(D(J).GE.P)THEN
            K=J
            P=D(J)
          ENDIF
11      CONTINUE
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO 12 J=1,N
            P=V(J,I)
            V(J,I)=V(J,K)
            V(J,K)=P
12        CONTINUE
        ENDIF
13    CONTINUE

end subroutine






subroutine quicksort(flag, list, cutvalue, ndim, sorted, nsize, cutline)

! integer list sort algorithm

implicit none

integer, intent(in) :: flag, ndim, list(ndim), cutvalue
integer, intent(out) :: sorted(ndim), nsize, cutline

logical :: existence(ndim)
integer :: i, j, number
integer :: new_list(ndim)


! remove the repeated values
existence = .TRUE.
do i = 1, ndim-1
  if (existence(i)) then
    number = list(i)
    do j = (i+1), ndim
      if (number == list(j)) existence(j) = .FALSE.
    end do
  endif
end do

nsize = COUNT(existence)
new_list(1:nsize) = PACK(list,existence)

! sorting procedure
sorted = 0
cutline = 0
existence = .TRUE.
do i = 1, nsize
  j = MINLOC(new_list(1:nsize), 1, existence(1:nsize))
  sorted(i) = new_list(j)
  existence(j) = .FALSE.
  if (sorted(i) <= cutvalue) cutline = cutline + 1
end do

end subroutine quicksort
