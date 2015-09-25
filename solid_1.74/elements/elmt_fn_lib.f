subroutine bbar2m(phi,shp,shpr,dvol,xji,lint,nel,npm,hh,theta,shpbar)


!     Purpose: Compute mixed formulation for the volumetric response

!     Inputs:
!        phi(6,*)       - Pressure/volume functions
!        shp(3,16,*)    - Shape functions and derivatives at gauss points
!        shpr(16,*)     - Shape functions divided by radius or zero
!        vol(*)         - Volume elements at gauss points at t_n+1
!        xji(2,*)       - Jacobian determinant at gauss points
!        lint           - Number of quadrature points
!        nel            - Number of nodes on element (should be 8)
!        npm            - Number of mixed modes

!     Outputs:
!        hh(6,6)        - Reference configuration shape integrals (inverse)
!        theta(2,*)     - Mixed jacobian determinant for element
!        shpbar(2,16,*) - Mixed derivatives of shape functions.

implicit none

integer   lint,  nel,  npm,  i,  j,  k,  l
real*8    phi(6,*),shp(3,16,*),shpr(16,*),dvol(*),theta(2,*)
real*8    shpbar(2,16,*),gg(6,2,16),hh(6,6),ji(6,2),xji(2,*)
real*8    hj(6,2), hg(6,2,16), voln, h0, h1, h2


!     Constant pressure elements

if (npm.eq.1) then

  do j = 1,nel
    shpbar(1,j,1) = 0.0d0
    shpbar(2,j,1) = 0.0d0
  end do ! j
  hh(1,1) = 0.d0
  h1      = 0.d0
  h2      = 0.d0

  do l = 1,lint

!         H-array and D-array

    voln    = dvol(l) / xji(1,l)
    hh(1,1) = hh(1,1) + voln
    h1      = h1      + dvol(l)
    h2      = h2      + voln*xji(2,l)
    
!         G-array

    do j = 1,nel
      shpbar(1,j,1) = shpbar(1,j,1) + (shp(1,j,l) + shpr(j,l))* dvol(l)
      shpbar(2,j,1) = shpbar(2,j,1) + shp(2,j,l)* dvol(l)
    end do ! j
  end do ! l

!       Modify shpbar for B-bar type computations

  h0 = 1.d0/h1
  do j = 1,nel
    do i = 1,2
      shpbar(i,j,1) = shpbar(i,j,1)*h0
    end do ! i
  end do ! j

!       Average Jacobian

  hh(1,1)    = 1.d0 / hh(1,1)
  theta(1,1) = h1   * hh(1,1)
  theta(2,1) = h2   * hh(1,1)
  do l = 2,lint
    theta(1,l) = theta(1,1)
    theta(2,l) = theta(2,1)
    do j = 1,nel
      shpbar(1,j,l) = shpbar(1,j,1)
      shpbar(2,j,l) = shpbar(2,j,1)
    end do ! j
    
  end do ! l

!     Higher order elements

else

  do i = 1,npm
    do j = 1,nel
      gg(i,1,j) = 0.0d0
      gg(i,2,j) = 0.0d0
    end do ! j
    do j = 1,npm
      hh(j,i) = 0.0d0
    end do ! j
    ji(i,1) = 0.0d0
    ji(i,2) = 0.0d0
  end do ! i

!       Quadrature loop

  do l = 1,lint
    do j = 1,npm

      h1 = phi(j,l)*dvol(l)
      h0 = h1/xji(1,l)
      h2 = h1*xji(2,l)

!           Ji-array

      ji(j,1) = ji(j,1) + h1
      ji(j,2) = ji(j,2) + h1*xji(2,l)/xji(1,l)

!           H-array

      do i = 1,npm
        hh(i,j)    = hh(i,j)  + phi(i,l)*h0
      end do ! i

!           G-array

      do i = 1,nel
        gg(j,1,i) = gg(j,1,i) + (shp(1,i,l) + shpr(i,l))*h1
        gg(j,2,i) = gg(j,2,i) +  shp(2,i,l)*h1
      end do ! i
    end do ! j

  end do ! l

!       Invert H-array

  call invert(hh,npm,6)

  do j = 1,2
    do i = 1,npm
      hj(i,j) = 0.0d0
      do k = 1,npm
        hj(i,j) = hj(i,j) + hh(i,k)*ji(k,j)
      end do ! k
    end do ! i
  end do ! j

  do j = 1,nel
    do i = 1,npm
      hg(i,1,j) = 0.0d0
      hg(i,2,j) = 0.0d0
      do k = 1,npm
        hg(i,1,j) = hg(i,1,j) + hh(i,k)*gg(k,1,j)
        hg(i,2,j) = hg(i,2,j) + hh(i,k)*gg(k,2,j)
      end do ! k
    end do ! i
  end do ! j

  do l = 1,lint
    theta(1,l) = hj(1,1)
    theta(2,l) = hj(1,2)
    do k = 2,npm
      theta(1,l) = theta(1,l) + phi(k,l)*hj(k,1)
      theta(2,l) = theta(2,l) + phi(k,l)*hj(k,2)
    end do ! k
    h0         = 1.d0/theta(1,l)
    do j = 1,nel
      shpbar(1,j,l) = hg(1,1,j)
      shpbar(2,j,l) = hg(1,2,j)
      do k = 2,npm
        shpbar(1,j,l) = shpbar(1,j,l) + phi(k,l)*hg(k,1,j)
        shpbar(2,j,l) = shpbar(2,j,l) + phi(k,l)*hg(k,2,j)
      end do ! k
      shpbar(1,j,l) = h0*shpbar(1,j,l)
      shpbar(2,j,l) = h0*shpbar(2,j,l)
    end do ! j
  end do ! l

endif

end subroutine bbar2m






subroutine bbar3m(phi,shp,dvol,xji,lint,nel,npm,hh,theta,shpbar)


!     Purpose: Compute mixed formulation for the volumetric response

!     Inputs:
!        shp(4,27,*)    - Shape functions and derivatives at gauss points
!        vol(*)         - Volume elements at gauss points at t_n+1
!        xji(2,*)       - Jacobian determinant at gauss points
!        lint           - Number of quadrature points
!        nel            - Number of nodes on element
!        npm            - Number of pressure modes/element

!     Outputs:
!        hh(4,4)        - Reference configuration shape integrals (inverse)
!        theta(2,*)     - Mixed jacobian determinant for element
!        shpbar(3,27,*) - Mixed derivatives of shape functions.

implicit none

integer, intent(in) :: lint, nel, npm
real, intent(in) :: phi(4,*), shp(4,27,*), dvol(*), xji(2,*)
real, intent(out) :: hh(4,4), theta(2,*), shpbar(3,27,*)
integer :: i, j, k, l
real ::  gg(4,3,27), ji(4,3), hj(4,3), hg(4,3,27), voln, h0, h1, h2


!     Constant pressure elements

if (npm.eq.1) then

  do j = 1,nel
    do i = 1,3
      shpbar(i,j,1) = 0.0d0
    end do ! i
  end do ! j
  hh(1,1)  = 0.d0
  h1       = 0.0d0
  h2       = 0.0d0

  do l = 1,lint

!         H-array and D-array

    voln     = dvol(l)/xji(1,l)
    hh(1,1)  = hh(1,1)  + voln
    h1       = h1       + dvol(l)
    h2       = h2       + voln*xji(2,l)

!         G-array

    do j = 1,nel
      do i = 1,3
        shpbar(i,j,1) = shpbar(i,j,1) + shp(i,j,l) * dvol(l)
      end do ! i
    end do ! j
  end do ! l

!       Modify shpbar for B-bar type computations

  h0 = 1.d0/h1

  do j = 1,nel
    do i = 1,3
      shpbar(i,j,1) = shpbar(i,j,1)*h0
    end do ! i
  end do ! j

!       Average Jacobian at t_n+1, t_n

  hh(1,1)     = 1.d0 / hh(1,1)
  theta(1,1)  = h1  *  hh(1,1)
  theta(2,1)  = h2  *  hh(1,1)

  do l = 2,lint
    theta(1,l) = theta(1,1)
    theta(2,l) = theta(2,1)
    do j = 1,nel
      do i = 1,3
        shpbar(i,j,l) = shpbar(i,j,1)
      end do ! i
    end do ! j
  end do ! l

!     Higher order elements

else

  do i = 1,npm
    do j = 1,nel
      gg(i,1,j) = 0.0d0
      gg(i,2,j) = 0.0d0
      gg(i,3,j) = 0.0d0
    end do ! j
    do j = 1,npm
      hh(j,i) = 0.0d0
    end do ! j
    ji(i,1) = 0.0d0
    ji(i,2) = 0.0d0
  end do ! i

!       Quadrature loop

  do l = 1,lint
    do j = 1,npm
      h1 = phi(j,l)*dvol(l)
      h0 = h1/xji(1,l)
      h2 = h1*xji(2,l)

!           Ji-array

      ji(j,1) = ji(j,1) + h1
      ji(j,2) = ji(j,2) + h1*xji(2,l)/xji(1,l)

!           H-array

      do i = 1,npm
        hh(i,j)    = hh(i,j)  + phi(i,l)*h0
      end do ! i

!           G-array

      do i = 1,nel
        do k = 1,3
          gg(j,k,i) = gg(j,k,i) + shp(k,i,l)*h1
        end do ! k
      end do ! i
    end do ! j

  end do ! l

!       Invert H-array

  call invert(hh,npm,4)

  do j = 1,2
    do i = 1,npm
      hj(i,j) = 0.0d0
      do k = 1,npm
        hj(i,j) = hj(i,j) + hh(i,k)*ji(k,j)
      end do ! k
    end do ! i
  end do ! j

  do j = 1,nel
    do i = 1,npm
      hg(i,1,j) = 0.0d0
      hg(i,2,j) = 0.0d0
      hg(i,3,j) = 0.0d0
      do k = 1,npm
        hg(i,1,j) = hg(i,1,j) + hh(i,k)*gg(k,1,j)
        hg(i,2,j) = hg(i,2,j) + hh(i,k)*gg(k,2,j)
        hg(i,3,j) = hg(i,3,j) + hh(i,k)*gg(k,3,j)
      end do ! k
    end do ! i
  end do ! j

  do l = 1,lint
    theta(1,l) = hj(1,1)
    theta(2,l) = hj(1,2)
    do k = 2,npm
      theta(1,l) = theta(1,l) + phi(k,l)*hj(k,1)
      theta(2,l) = theta(2,l) + phi(k,l)*hj(k,2)
    end do ! k
    h0 = 1.d0/theta(1,l)
    do j = 1,nel
      do i = 1,3
        shpbar(i,j,l) = hg(1,i,j)
        do k = 2,npm
          shpbar(i,j,l) = shpbar(i,j,l) + phi(k,l)*hg(k,i,j)
        end do ! k
        shpbar(i,j,l) = h0*shpbar(i,j,l)
      end do ! k
    end do ! j
  end do ! l

endif

end subroutine bbar3m






subroutine dmatdx(dd,sigb,p_bar,p_mix)


!     Purpose: Compute finite deformation mixed D arrays

!     Inputs:
!         dd(7,7)       Modified constitutive array
!         sigb(6)       Constitutive stresses
!         p_bar         Constitutive pressure
!         p_mix         Mixed pressure

!     Outputs:
!         dd(7,7)       Material tangent terms

implicit none


real, intent(in) :: sigb(6), p_bar, p_mix
real, intent(inout) :: dd(7,7)
integer :: i, j
real :: sigb_d(6), fac1
real, parameter :: two3 = 2.0/3.0


!     Compute deviator stress

do i = 1,3
  sigb_d(i  ) = two3 *(sigb(i  ) - p_bar)
  sigb_d(i+3) = two3 * sigb(i+3)
end do ! i

!     D_11: B_matrix part

fac1 = p_mix - two3*p_bar
do j = 1,3
  do i = 1,3
    dd(i,j) = dd(i,j) + fac1
  end do ! i
end do ! j

do j = 1,6
  do i = 1,3
    dd(i,j) = dd(i,j) - sigb_d(j)
    dd(j,i) = dd(j,i) - sigb_d(j)
  end do ! i
end do ! j

fac1 = p_bar - p_mix
do j = 1,3
  dd(j  ,j  ) = dd(j  ,j  ) + fac1*2.d0
  dd(j+3,j+3) = dd(j+3,j+3) + fac1
end do ! j

!     D_12: Coupling matrix with

do j = 1,6
  dd(7,j) = dd(7,j) + sigb_d(j)
  dd(j,7) = dd(j,7) + sigb_d(j)
end do ! j

!     D_22: Volumetric part

dd(7,7) = dd(7,7) - p_bar/3.0

end subroutine dmatdx






subroutine fbar(f,xji,theta,lint)


!     Purpose: Form F-bar and left Cauchy-Green tensors

!     Inputs:
!        f(3,3)      - Deformation gradient
!        xji(2,*)    - Determinant of deformation gradient (J)
!        theta(2,*)  - Mixed determinant of deformation gradient
!        lint        - Number of quadrature points

!     Outputs:
!        f(3,3)      - Mixed deformation gradient

implicit none

integer, intent(in) :: lint
real, intent(in) :: xji(2,*), theta(2,*)
real, intent(inout) :: f(9,2,*)
integer :: i, l
real :: ration, ratio1
real, parameter :: one3 = 1.0/3.0


!     Compute mixed deformation gradient

do l = 1, lint
  ratio1 = abs(theta(1,l)/xji(1,l))**one3
  ration = abs(theta(2,l)/xji(2,l))**one3
  do i = 1,9
    f(i,1,l) = ratio1*f(i,1,l)
    f(i,2,l) = ration*f(i,2,l)
  end do ! i
end do ! l

end subroutine fbar







subroutine kine2d(pro_type,shp,xl,ul,f,finv,df,detf,ndm,ndf,nel,lint)

!      Purpose: Compute kinematic quantities for finite deformations
!                >>  Only for 2D

!      Inputs:
!         shp(3,16,*)   - Reference configuration shape functions
!         xl(ndm,*)     - Nodal reference coordinates
!         ul(ndf,nel,*) - Nodal displacements
!         ndm           - Number mesh dimensions
!         ndf           - Number dof/node
!         nel           - Number nodes/element
!         lint          - Number of quadrature points

!      Outputs:
!         f(3,3,2,*)    - Deformation gradient
!         finv(3,3,*)   - Inverse deformation gradient
!         df(3,3,*)     - Incremental deformation gradient
!         detf(2,*)     - Determinant of deformation gradient

implicit none

integer, intent(in) :: pro_type, ndm, ndf, nel, lint
real, intent(in) :: xl(ndm,*), ul(ndf,nel,*)
real, intent(out) :: f(3,3,2,*), finv(3,3,*), df(3,3,*), detf(2,*)
real, intent(inout) :: shp(3,16,*)

integer :: i, j ,k, l
real :: detfi, temp, xx1, uu1, du1


!     Deformation gradient at t_n+1 : F_n+1 = I + GRAD u_n+1

do l = 1,lint

  do i = 1,3
    do j = 1,3
      f(j,i,1,l) = 0.0
      df(j,i,l)  = 0.0
    end do ! j
  end do ! i
  do i = 1,2
    do j = 1,2
      do k = 1,nel
        f(i,j,1,l) = f(i,j,1,l) + ul(i,k,1)*shp(j,k,l)
        df(i,j,l ) = df(i,j,l ) + ul(i,k,2)*shp(j,k,l)
      end do ! k
    end do ! j
    f(i,i,1,l) = f(i,i,1,l) + 1.0
  end do ! i

!       Axisymmetry

  if(pro_type == 2) then  ! Axisymmetry case
    xx1 = 0.0
    uu1 = 0.0
    du1 = 0.0
    do k = 1,nel
      xx1 = xx1 + xl(1,k  )*shp(3,k,l)
      uu1 = uu1 + ul(1,k,1)*shp(3,k,l)
      du1 = du1 + ul(1,k,2)*shp(3,k,l)
    end do ! k
    df(3,3,l)  = du1/xx1
    f(3,3,1,l) = uu1/xx1 + 1.0

  else
    f(3,3,1,l) = 1.0
    f(3,3,2,l) = 1.0
    df(3,3,l)  = 0.0
  endif

!       Deformation gradient at t_n: F_n

  do i = 1,3
    do j = 1,3
      f(j,i,2,l) = f(j,i,1,l) - df(j,i,l)
    end do ! j
  end do ! i

!       Invert F

  detf(1,l) = f(1,1,1,l)*f(2,2,1,l) - f(1,2,1,l)*f(2,1,1,l)
  detf(2,l) = f(1,1,2,l)*f(2,2,2,l) - f(1,2,2,l)*f(2,1,2,l)
  
  detfi   =  1.0/detf(1,l)
  finv(1,1,l) =  f(2,2,1,l)*detfi
  finv(1,2,l) = -f(1,2,1,l)*detfi
  finv(1,3,l) =  0.0
  finv(2,1,l) = -f(2,1,1,l)*detfi
  finv(2,2,l) =  f(1,1,1,l)*detfi
  finv(2,3,l) =  0.0
  finv(3,3,l) =  1.0/f(3,3,1,l)

!       Determinants

  detf(1,l) = detf(1,l)*f(3,3,1,l)
  detf(2,l) = detf(2,l)*f(3,3,2,l)

  finv(3,1,l) = 0.0
  finv(3,2,l) = 0.0
  

!       Transform shape functions to current configuration

  do k = 1,nel
    temp       = finv(1,1,l)*shp(1,k,l) + finv(2,1,l)*shp(2,k,l)
    shp(2,k,l) = finv(1,2,l)*shp(1,k,l) + finv(2,2,l)*shp(2,k,l)
    shp(1,k,l) = temp
  end do ! k
  

end do ! l

end subroutine kine2d








subroutine kine3d(shp,ul,f,fi,df,detf,ndf,nel)

!     Purpose: Compute deformation gradient and its inverse at tn+1

!     Inputs:
!        shp(4,*)  - Shape functions
!        ul(ndf,*) - Nodal solution values
!        ndf       - Degrees of freedom / node
!        nel       - Number of element nodes

!     Outputs:
!        f(3,3,2)  - Deformation gradients
!        fi(3,3)   - Inverse deformation gradient
!        df(3,3)   - Incremental deformation gradient
!        detf(2)  - Determinant of deformation gradient

implicit none

integer, intent(in) :: ndf, nel
real, intent(in) :: ul(ndf,nel,*)
real, intent(out) :: f(3,3,2), fi(3,3), df(3,3), detf(2)
real, intent(inout) :: shp(4,*)

integer :: i, j, k
real :: deti, temp(3)


!     Compute compatible deformation gradient at t-n+1: F = I + GRAD u

do i = 1,3
  do j = 1,3
    f(i,j,1) = 0.0
    f(i,j,2) = 0.0
    df(i,j)  = 0.0
    do k = 1, nel
      f(i,j,1) = f(i,j,1) + ul(i,k,1)*shp(j,k)
      df(i,j)  = df(i,j)  + ul(i,k,2)*shp(j,k)
    end do ! k
    f(i,j,2) = f(i,j,1) - df(i,j)
  end do ! j
  f(i,i,1) = f(i,i,1) + 1.0
  f(i,i,2) = f(i,i,2) + 1.0
end do ! i

!     Invert F_n

detf(2) = f(1,1,2)*f(2,2,2)*f(3,3,2) + f(1,2,2)*f(2,3,2)*f(3,1,2) &
        + f(1,3,2)*f(2,1,2)*f(3,2,2) - f(3,1,2)*f(2,2,2)*f(1,3,2) &
        - f(3,2,2)*f(2,3,2)*f(1,1,2) - f(3,3,2)*f(2,1,2)*f(1,2,2)

!     Invert F_n+1

detf(1) = f(1,1,1)*f(2,2,1)*f(3,3,1) + f(1,2,1)*f(2,3,1)*f(3,1,1) &
        + f(1,3,1)*f(2,1,1)*f(3,2,1) - f(3,1,1)*f(2,2,1)*f(1,3,1) &
        - f(3,2,1)*f(2,3,1)*f(1,1,1) - f(3,3,1)*f(2,1,1)*f(1,2,1)

deti    = 1.d0/detf(1)
fi(1,1) = (f(2,2,1)*f(3,3,1) - f(3,2,1)*f(2,3,1))*deti
fi(1,2) =-(f(1,2,1)*f(3,3,1) - f(3,2,1)*f(1,3,1))*deti
fi(1,3) = (f(1,2,1)*f(2,3,1) - f(2,2,1)*f(1,3,1))*deti
fi(2,1) =-(f(2,1,1)*f(3,3,1) - f(3,1,1)*f(2,3,1))*deti
fi(2,2) = (f(1,1,1)*f(3,3,1) - f(3,1,1)*f(1,3,1))*deti
fi(2,3) =-(f(1,1,1)*f(2,3,1) - f(2,1,1)*f(1,3,1))*deti
fi(3,1) = (f(2,1,1)*f(3,2,1) - f(3,1,1)*f(2,2,1))*deti
fi(3,2) =-(f(1,1,1)*f(3,2,1) - f(3,1,1)*f(1,2,1))*deti
fi(3,3) = (f(1,1,1)*f(2,2,1) - f(2,1,1)*f(1,2,1))*deti

!     Push forward standard shape functions

do k = 1, nel
  do i = 1, 3
    temp(i) = fi(1,i)*shp(1,k) + fi(2,i)*shp(2,k) + fi(3,i)*shp(3,k)
  end do ! i
  do i = 1, 3
    shp(i,k) = temp(i)
  end do ! i
end do ! k

end subroutine kine3d







subroutine kine3m(shp,ul,f,fi,df,detf,ndf,nel,lint)



!      Purpose: Compute kinematic quantities for finite deformations

!      Inputs:
!         shp(4,27,*)   - Reference configuration shape functions
!         ul(ndf,nel,*) - Nodal displacements
!         ndf           - Number dof/node
!         nel           - Number nodes/element
!         lint          - Number of quadrature points

!      Outputs:
!         f(3,3,2,*)    - Deformation gradient
!         fi(3,3,*)     - Inverse deformation gradient
!         df(3,3,*)     - Incremental deformation gradient
!         detf(2,*)     - Determinant of deformation gradient

implicit none

integer, intent(in) :: ndf, nel, lint
real, intent(in) :: ul(ndf,nel,*)
real, intent(out) :: f(3,3,2,*), fi(3,3,*), df(3,3,*), detf(2,*)
real, intent(inout) :: shp(4,27,*)

integer :: i, j, k, l
real :: detfi, cc(3)


!     Deformation gradient at t_n+1 : F_n+1 = I + GRAD u_n+1

do l = 1,lint
  do i = 1,3
    do j = 1,3
      f(j,i,1,l) = 0.0d0
      df(j,i,l)  = 0.0d0
    end do ! j
  end do ! i
  do i = 1,3
    do j = 1,3
      do k = 1,nel
        f(i,j,1,l) = f(i,j,1,l) + ul(i,k,1)*shp(j,k,l)
        df(i,j,l ) = df(i,j,l ) + ul(i,k,2)*shp(j,k,l)
      end do ! k
    end do ! j
    f(i,i,1,l) = f(i,i,1,l) + 1.0d0
  end do ! i

!       Deformation gradient at t_n: F_n

  do i = 1,3
    do j = 1,3
      f(j,i,2,l) = f(j,i,1,l) - df(j,i,l)
    end do ! j
  end do ! i

!       Invert F

  do i = 1,2
    detf(i,l) = f(1,1,i,l)*f(2,2,i,l)*f(3,3,i,l)  &
              + f(1,2,i,l)*f(2,3,i,l)*f(3,1,i,l)  &
              + f(1,3,i,l)*f(2,1,i,l)*f(3,2,i,l)  &
              - f(3,1,i,l)*f(2,2,i,l)*f(1,3,i,l)  &
              - f(3,2,i,l)*f(2,3,i,l)*f(1,1,i,l)  &
              - f(3,3,i,l)*f(2,1,i,l)*f(1,2,i,l)
  end do ! i

  detfi     = 1.d0/detf(1,l)
  fi(1,1,l) = (f(2,2,1,l)*f(3,3,1,l)-f(3,2,1,l)*f(2,3,1,l))*detfi
  fi(1,2,l) = (f(3,2,1,l)*f(1,3,1,l)-f(1,2,1,l)*f(3,3,1,l))*detfi
  fi(1,3,l) = (f(1,2,1,l)*f(2,3,1,l)-f(2,2,1,l)*f(1,3,1,l))*detfi

  fi(2,1,l) = (f(2,3,1,l)*f(3,1,1,l)-f(3,3,1,l)*f(2,1,1,l))*detfi
  fi(2,2,l) = (f(3,3,1,l)*f(1,1,1,l)-f(1,3,1,l)*f(3,1,1,l))*detfi
  fi(2,3,l) = (f(1,3,1,l)*f(2,1,1,l)-f(2,3,1,l)*f(1,1,1,l))*detfi

  fi(3,1,l) = (f(2,1,1,l)*f(3,2,1,l)-f(3,1,1,l)*f(2,2,1,l))*detfi
  fi(3,2,l) = (f(3,1,1,l)*f(1,2,1,l)-f(1,1,1,l)*f(3,2,1,l))*detfi
  fi(3,3,l) = (f(1,1,1,l)*f(2,2,1,l)-f(2,1,1,l)*f(1,2,1,l))*detfi

!       Transform shape functions to current configuration

  do k = 1,nel
    do i = 1,3
      cc(i) = fi(1,i,l)*shp(1,k,l)  &
            + fi(2,i,l)*shp(2,k,l)  &
            + fi(3,i,l)*shp(3,k,l)
    end do ! i
    do i = 1,3
      shp(i,k,l) = cc(i)
    end do ! i
  end do ! k
end do ! l

end subroutine kine3m
