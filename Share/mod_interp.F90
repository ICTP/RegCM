!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_interp

  use mod_dynparam , only : ds
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_stdio
  use mod_message
  use mod_memutil
  use mod_constants
  use mod_earth

  implicit none

  private

  public :: bilinx , cressmcr , cressmdt , distwgtcr , distwgtdt
  public :: kernsmooth
  public :: interp1d

  real(rkx) :: alatmn , alatmx , alonmn , alonmx
  real(rkx) :: glatmn , glatmx , glonmn , glonmx

  integer(ik4) :: imxmn = 0
  integer(ik4) :: lcross = 0
  integer(ik4) :: ldot = 0

  real(rkx) , parameter :: deg720 = d_two*deg360
  real(rkx) , parameter :: missl = -9999.0_rkx
  real(rkx) , parameter :: missc = -9990.0_rkx
  integer(ik4) , parameter :: p_factor = 2

  real(rkx) , pointer , dimension(:,:) :: dc1xa , dc1xb , dc1xc , dc1xd
  real(rkx) , pointer , dimension(:,:) :: dd1xa , dd1xb , dd1xc , dd1xd
  integer(ik4) , pointer , dimension(:,:) :: ic1dl , ic1dr , ic1ul , ic1ur , &
                                       jc1dl , jc1dr , jc1ul , jc1ur
  integer(ik4) , pointer , dimension(:,:) :: id1dl , id1dr , id1ul , id1ur , &
                                       jd1dl , jd1dr , jd1ul , jd1ur
  logical :: lonwrap = .false.
  logical :: latpole = .false.

  interface interp1d
    module procedure interp1d_r4
    module procedure interp1d_r8
  end interface interp1d

  interface cressmcr
    module procedure cressmcr3d
    module procedure cressmcr2d
  end interface cressmcr

  interface bilinx
    module procedure bilinx_2d
    module procedure bilinx_3d
  end interface bilinx

  interface kernsmooth
    module procedure kernsmooth2
    module procedure kernsmooth3
  end interface kernsmooth

  contains

  ! Interpolates with splines with tension in one dimension.
  ! The spline is defined imposing that the second derivative is
  ! the average of second derivatives computed at the two adjacent points.
  ! At interval extremes the second derivative is assumed null.
  ! This subroutine also extrapolates out of the interval where the
  ! input funtion g is defined
  subroutine interp1d_r4(xi,g,xo,f,alfa,ex1,ex2)
    implicit none
    !  Input:  function g defined at irregular but strictly monotonic xi
    !  Output: f interpolated values at arbitrary coordinates xo
    real(rk4) , dimension(:) , intent(in) :: xi , xo , g
    real(rk4) , dimension(:) , intent(out) :: f
    ! alfa: spline tension parameter, comprised between 0 and 1:
    ! if alfa=1, pure linear interpolation; if alfa=0, pure spline
    real(rk4) , intent(in) :: alfa
    ! ex1: param. determining extrapolation for x < xi(1)
    ! ex2: param. determining extrapolation for x > xi(npi)
    ! if ex1=0 or ex2=0, constant value extrapolation is used at extreme
    ! if ex1=1 or ex2=1, linear extrapolation is used at corresponding extreme
    ! intermediate values of ex1 and ex2 give intermediate extrapolation values
    real(rk4) , intent(in) :: ex1
    real(rk4) , intent(in) :: ex2

    real(rk4) :: zeps , ximed , gmed , fmm , fpp , xmm , xpp
    real(rk4) :: fm , xm , fp , xp , delx , delxp , delxm
    real(rk4) :: delx1 , delx2 , delxs , delx1s , delx2s
    real(rk4) :: spl , clin

    real(rk4) , dimension(size(xi)) :: zi , zg
    integer(ik4) :: npi , npo
    integer(ik4) :: k , j , jj , ir

    npi = size(xi)
    npo = size(xo)

    if ( npi < 2 ) then
      write(stderr,*) 'Refusing to work: too few input points.'
      call die('interp1d')
    end if

    if ( size(g) /= npi .or. size(f) /= npo ) then
      write(stderr,*) 'Refusing to work: different size coordinate/values'
      call die('interp1d')
    end if

    if ( xi(1) >= xi(npi) ) then
      do k = 1 , npi
        zi(k) = xi(npi-k+1)
        zg(k) = g(npi-k+1)
      end do
    else
      zi(:) = xi(:)
      zg(:) = g(:)
    end if

    ! small deviation used to set apart interlaced coordinates
    zeps = (zi(npi) - zi(1)) * 1.e-6_rk4
    deinterlace: &
    do
      do k = 2 , npi
        if ( zi(k) <= zi(k-1) ) then
          ximed = 0.5_rk4 * (zi(k) + zi(k-1))
          zi(k-1) = ximed - zeps
          zi(k) = ximed + zeps
          gmed = 0.5_rk4 * (zg(k) + zg(k-1))
          zg(k-1) = gmed
          zg(k) = gmed
        end if
      end do

      do k = 2 , npi
        if ( zi(k) <= zi(k-1) ) then
          cycle deinterlace
        end if
      end do
      exit deinterlace
    end do deinterlace

    do j = 1 , npo
      ! 2 cases of extrapolation
      if ( xo(j) < zi(1) ) then
        f(j) = zg(1) + ex1*(zg(1)-zg(2))/(zi(1)-zi(2)) * (xo(j)-zi(1))
        cycle
      else if ( xo(j) >= zi(npi) ) then
        f(j) = zg(npi) + ex2*(zg(npi)-zg(npi-1))/(zi(npi)-zi(npi-1)) * &
               (xo(j)-zi(npi))
        cycle
      end if
      ir = 0
      ! ir is a reference index determining the interpolation interval
      ! The interpolation expression is applied also if xo = zi(j)
      do jj = 1 , npi
        if ( xo(j) >= zi(jj) ) ir = ir + 1
      end do
      if ( ir == 1 ) then
        fmm = 2.0_rk4 * zg(1) - zg(2)
        xmm = 2.0_rk4 * zi(1) - zi(2)
        fpp = zg(ir+2)
        xpp = zi(ir+2)
      else if ( ir == (npi-1) ) then
        fpp = 2.0_rk4 * zg(npi) - zg(npi-1)
        xpp = 2.0_rk4 * zi(npi) - zi(npi-1)
        fmm = zg(ir-1)
        xmm = zi(ir-1)
      else
        fmm = zg(ir-1)
        xmm = zi(ir-1)
        fpp = zg(ir+2)
        xpp = zi(ir+2)
      end if
      fm     = zg(ir)
      xm     = zi(ir)
      fp     = zg(ir+1)
      xp     = zi(ir+1)
      delx   = xp - xm
      delxp  = xpp - xp
      delxm  = xm - xmm
      delx1  = xo(j) - xm
      delx2  = xp - xo(j)
      delxs  = delx**2
      delx1s = delx1**2
      delx2s = delx2**2
      !  Spline contribution to interpolation
      spl = fm*(delx2/delx + delx1*delx2s/(delxs*delxm) - delx1s*     &
            delx2/((delx+delxp)*delxs)) + fp*(delx1/delx +            &
            delx1s*delx2/(delxs*delxp) - delx1*delx2s/((delx+delxm)*  &
            delxs)) - fmm * delx1*delx2s/((delx+delxm)*delx*delxm) -  &
            fpp * delx1s*delx2/((delx+delxp)*delx*delxp)
      !  Linear interpolation contribution
      clin = (fm*delx2 + fp*delx1)/delx
      !  Final interpolation combined using alfa
      f(j) = alfa*clin + (1.0_rk4-alfa)*spl
    end do
  end subroutine interp1d_r4

  subroutine interp1d_r8(xi,g,xo,f,alfa,ex1,ex2)
    implicit none
    real(rk8) , dimension(:) , intent(in) :: xi , xo , g
    real(rk8) , dimension(:) , intent(out) :: f
    real(rk8) , intent(in) :: alfa
    real(rk8) , intent(in) :: ex1
    real(rk8) , intent(in) :: ex2

    real(rk8) :: zeps , ximed , gmed , fmm , fpp , xmm , xpp
    real(rk8) :: fm , xm , fp , xp , delx , delxp , delxm
    real(rk8) :: delx1 , delx2 , delxs , delx1s , delx2s
    real(rk8) :: spl , clin

    real(rk8) , dimension(size(xi)) :: zi , zg
    integer(ik4) :: npi , npo
    integer(ik4) :: k , j , jj , ir

    npi = size(xi)
    npo = size(xo)

    if ( npi < 2 ) then
      write(stderr,*) 'Refusing to work: too few input points.'
      call die('interp1d')
    end if

    if ( size(g) /= npi .or. size(f) /= npo ) then
      write(stderr,*) 'Refusing to work: different size coordinate/values'
      call die('interp1d')
    end if

    if ( xi(1) >= xi(npi) ) then
      do k = 1 , npi
        zi(k) = xi(npi-k+1)
        zg(k) = g(npi-k+1)
      end do
    else
      zi(:) = xi(:)
      zg(:) = g(:)
    end if

    zeps = (zi(npi) - zi(1)) * 1.e-6_rk8
    deinterlace: &
    do
      do k = 2 , npi
        if ( zi(k) <= zi(k-1) ) then
          ximed = 0.5_rk8 * (zi(k) + zi(k-1))
          zi(k-1) = ximed - zeps
          zi(k) = ximed + zeps
          gmed = 0.5_rk8 * (zg(k) + zg(k-1))
          zg(k-1) = gmed
          zg(k) = gmed
        end if
      end do

      do k = 2 , npi
        if ( zi(k) <= zi(k-1) ) then
          cycle deinterlace
        end if
      end do
      exit deinterlace
    end do deinterlace

    do j = 1 , npo
      if ( xo(j) < zi(1) ) then
        f(j) = zg(1) + ex1*(zg(1)-zg(2))/(zi(1)-zi(2)) * (xo(j)-zi(1))
        cycle
      else if ( xo(j) >= zi(npi) ) then
        f(j) = zg(npi) + ex2*(zg(npi)-zg(npi-1))/(zi(npi)-zi(npi-1)) * &
               (xo(j)-zi(npi))
        cycle
      end if
      ir = 0
      do jj = 1 , npi
        if ( xo(j) >= zi(jj) ) ir = ir + 1
      end do
      if ( ir == 1 ) then
        fmm = 2.0_rk8 * zg(1) - zg(2)
        xmm = 2.0_rk8 * zi(1) - zi(2)
        fpp = zg(ir+2)
        xpp = zi(ir+2)
      else if ( ir == (npi-1) ) then
        fpp = 2.0_rk8 * zg(npi) - zg(npi-1)
        xpp = 2.0_rk8 * zi(npi) - zi(npi-1)
        fmm = zg(ir-1)
        xmm = zi(ir-1)
      else
        fmm = zg(ir-1)
        xmm = zi(ir-1)
        fpp = zg(ir+2)
        xpp = zi(ir+2)
      end if
      fm     = zg(ir)
      xm     = zi(ir)
      fp     = zg(ir+1)
      xp     = zi(ir+1)
      delx   = xp - xm
      delxp  = xpp - xp
      delxm  = xm - xmm
      delx1  = xo(j) - xm
      delx2  = xp - xo(j)
      delxs  = delx**2
      delx1s = delx1**2
      delx2s = delx2**2
      spl = fm*(delx2/delx + delx1*delx2s/(delxs*delxm) - delx1s*     &
            delx2/((delx+delxp)*delxs)) + fp*(delx1/delx +            &
            delx1s*delx2/(delxs*delxp) - delx1*delx2s/((delx+delxm)*  &
            delxs)) - fmm * delx1*delx2s/((delx+delxm)*delx*delxm) -  &
            fpp * delx1s*delx2/((delx+delxp)*delx*delxp)
      clin = (fm*delx2 + fp*delx1)/delx
      f(j) = alfa*clin + (1.0_rk8-alfa)*spl
    end do
  end subroutine interp1d_r8

  subroutine bilinx_2d(b3,b2,alon,alat,hlon,hlat,nlon,nlat,jx,iy)
    implicit none
    integer(ik4) , intent(in) :: iy , jx , nlat , nlon
    real(rkx) , dimension(jx,iy) , intent(in) :: alat , alon
    real(rkx) , dimension(nlon,nlat) , intent(in) :: b2
    real(rkx) , dimension(nlat) , intent(in) :: hlat
    real(rkx) , dimension(nlon) , intent(in) :: hlon
    real(rkx) , dimension(jx,iy) , intent(out) :: b3
    real(rkx) :: p1 , p2 , q1 , q2 , dlon , dlat
    integer(ik4) :: i , i1 , i2 , j , j1 , j2
    !
    ! PERFORMING BI-LINEAR INTERPOLATION USING 4 GRID POINTS FROM A
    ! BIGGER RECTANGULAR GRID TO A GRID DESCRIBED BY ALON AND ALAT OF
    ! GRID2. A POINT ON GRID2 IS TRAPPED WITHIN FOUR GRID POINTS ON
    ! GRID4.THE GRID POINTS ARE ALWAYS TO THE NORTH AND EAST OF THE
    ! TRAPPED POINT. THE ALGORITHM COMPUTES THE FRACTIONAL DISTANCES IN
    ! BOTH X AND Y DIRECTION OF THE TRAPPED GRID POINT AND USES THE
    ! INFORMATION AS WEIGHTING FACTORS IN THE INTERPOLATION.
    !
    ! B2(JX,IX,NLEV) IS THE INPUT FIELD ON REGULAR LAT/LON GRID.
    ! B3(JX,IX,NLEV) IS THE OUTPUT FIELD ON PROJECTED GRID
    ! HLON......LONGITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
    ! HLAT......LATITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
    !
    dlon = abs(min(hlon(2)-hlon(1),hlon(nlon)-hlon(nlon-1)))
    dlat = abs(min(hlat(2)-hlat(1),hlat(nlat)-hlat(nlat-1)))
    do i = 1 , iy
      do j = 1 , jx
        j1 = whereislon(nlon,alon(j,i),hlon)
        ! Assume global data here
        if ( j1 > nlon ) j1 = 1
        if ( j1 < 1 ) j1 = nlon
        j2 = j1 + 1
        ! Assume global data here
        if ( j2 > nlon ) j2 = 1
        i1 = whereislat(nlat,alat(j,i),hlat)
        i2 = i1 + 1
        if ( b2(j1,i1) < missc .or. b2(j2,i1) < missc .or. &
             b2(j1,i2) < missc .or. b2(j2,i2) < missc ) then
          b3(j,i) = missl
        else
          p1 = mod(alon(j,i)-hlon(j1),dlon)
          if ( p1 < 0.0_rkx ) p1 = dlon + p1
          q1 = abs(mod(alat(j,i)-hlat(i1),dlat))
          q2 = dlat - q1
          p2 = dlon - p1
          b3(j,i) = ((b2(j1,i1)*p2+b2(j2,i1)*p1)*q2 + &
                     (b2(j1,i2)*p2+b2(j2,i2)*p1)*q1)/(dlon*dlat)
        end if
      end do
    end do
  end subroutine bilinx_2d

  subroutine bilinx_3d(b3,b2,alon,alat,hlon,hlat,nlon,nlat,jx,iy,llev)
    implicit none
    integer(ik4) , intent(in) :: iy , jx , llev , nlat , nlon
    real(rkx) , dimension(jx,iy) , intent(in) :: alat , alon
    real(rkx) , dimension(nlon,nlat,llev) , intent(in) :: b2
    real(rkx) , dimension(jx,iy,llev) , intent(out) :: b3
    real(rkx) , dimension(nlat) , intent(in) :: hlat
    real(rkx) , dimension(nlon) , intent(in) :: hlon
    integer(ik4) :: l
    do l = 1 , llev
      call bilinx_2d(b3(:,:,l),b2(:,:,l),alon,alat,hlon,hlat,nlon,nlat,jx,iy)
    end do
  end subroutine bilinx_3d

  subroutine compwgt(alon,alat,glon,glat,d1xa,d1xb,d1xc,d1xd, &
                     i1dl,i1dr,i1ul,i1ur,j1dl,j1dr,j1ul,j1ur, &
                     jx,iy,nlon,nlat)
    integer(ik4) , intent(in) :: iy , jx , nlat , nlon
    real(rkx) , dimension(jx,iy) , intent(in) :: alat , alon
    real(rkx) , dimension(nlon,nlat) , intent(in) :: glat , glon
    real(rkx) , dimension(jx,iy) , intent(out) :: d1xa , d1xb , d1xc , d1xd
    integer(ik4) , dimension(jx,iy) , intent(out) :: i1dl , i1dr , &
      i1ul , i1ur , j1dl , j1dr , j1ul , j1ur
    real(rkx) :: dist , wa , wb , wc , wd , distx
    integer(ik4) :: i , j , m , mdl , mdr , mul , mur , n , ndl ,  &
               ndr , nul , nur , mx , nx
    logical :: lsouthnorth , l360

    lsouthnorth = ( glat(1,1) < glat(nlon,nlat) )
    l360 = any( abs(glon) > 180.0_rkx )
    do i = 1 , iy
      do j = 1 , jx
        ! Find nearest point
        distx = 1.0e+20_rkx
        mx = -1
        nx = -1
        do n = 1 , nlat
          do m = 1 , nlon
            dist = gcdist(ds,glat(m,n),glon(m,n),alat(j,i),alon(j,i))
            if ( dist < distx ) then
              distx = dist
              mx = m
              nx = n
            end if
          end do
        end do
        if ( topof(glat(mx,nx),alat(j,i)) ) then
          if ( rightof(glon(mx,nx),alon(j,i)) ) then
            ! mx,nx is top right
            mur = mx
            nur = nx
            mul = mx-1
            nul = nx
            mdr = mx
            ndr = nx-1
            mdl = mx-1
            ndl = nx-1
          else
            ! mx,nx is top left
            mur = mx+1
            nur = nx
            mul = mx
            nul = nx
            mdr = mx+1
            ndr = nx-1
            mdl = mx
            ndl = nx-1
          end if
        else
          if ( rightof(glon(mx,nx),alon(j,i)) ) then
            ! mx,nx is bottom right
            mur = mx
            nur = nx+1
            mul = mx-1
            nul = nx+1
            mdr = mx
            ndr = nx
            mdl = mx-1
            ndl = nx
          else
            ! mx,nx is bottom left
            mur = mx+1
            nur = nx+1
            mul = mx
            nul = nx+1
            mdr = mx+1
            ndr = nx
            mdl = mx
            ndl = nx
          end if
        end if
        if ( lonwrap ) then
          if ( mur > nlon ) mur = mur - nlon
          if ( mdr > nlon ) mdr = mdr - nlon
          if ( mul < 1 ) mul = mul + nlon
          if ( mdl < 1 ) mdl = mdl + nlon
        end if
        if ( latpole ) then
          if ( nul > nlat ) nul = nlat
          if ( nur > nlat ) nur = nlat
          if ( ndl < 1 ) ndl = 1
          if ( ndr < 1 ) ndr = 1
        end if
        if ( mur < 1 .or. mur > nlon .or. &
             mul < 1 .or. mul > nlon .or. &
             mdl < 1 .or. mdl > nlon .or. &
             mdr < 1 .or. mdr > nlon .or. &
             nur < 1 .or. nur > nlat .or. &
             nul < 1 .or. nul > nlat .or. &
             ndl < 1 .or. ndl > nlat .or. &
             ndr < 1 .or. ndr > nlat ) then
          write (stderr,*) 'LOGIC ERROR in locating point'
          write (stderr,*) i , j , mx , nx , nlon , nlat
          write (stderr,*) alon(j,i)
          write (stderr,*) alat(j,i)
          call die('compwgt')
        end if
        i1ur(j,i) = mur
        j1ur(j,i) = nur
        i1ul(j,i) = mul
        j1ul(j,i) = nul
        i1dr(j,i) = mdr
        j1dr(j,i) = ndr
        i1dl(j,i) = mdl
        j1dl(j,i) = ndl
        wa = gcdist(ds,glat(mur,nur),glon(mur,nur),alat(j,i),alon(j,i))
        wb = gcdist(ds,glat(mul,nul),glon(mul,nul),alat(j,i),alon(j,i))
        wc = gcdist(ds,glat(mdr,ndr),glon(mdr,ndr),alat(j,i),alon(j,i))
        wd = gcdist(ds,glat(mdl,ndl),glon(mdl,ndl),alat(j,i),alon(j,i))
        d1xa(j,i) = d_one/(wa**p_factor)
        d1xb(j,i) = d_one/(wb**p_factor)
        d1xc(j,i) = d_one/(wc**p_factor)
        d1xd(j,i) = d_one/(wd**p_factor)
      end do
    end do

    contains

      logical function rightof(a,b)
        implicit none
        real(rkx) , intent(in) :: a , b
        if ( l360 ) then
          if ( a > 180.0_rkx ) then
            rightof = ( ( a - 360.0_rkx ) > b )
          else
            rightof = ( a > b )
          end if
        else
          if ( b > 180.0_rkx ) then
            rightof = ( a > ( b - 360.0_rkx ) )
          else
            rightof = ( a > b )
          end if
        end if
      end function rightof

      logical function topof(a,b)
        implicit none
        real(rkx) , intent(in) :: a , b
        if ( lsouthnorth ) then
          topof = ( a > b )
        else
          topof = ( b > a )
        end if
      end function topof

  end subroutine compwgt

  subroutine dwgt(jx,iy,nlon,nlat,b2,b3,d1xa,d1xb,d1xc,d1xd, &
                  i1dl,i1dr,i1ul,i1ur,j1dl,j1dr,j1ul,j1ur)
    implicit none
    integer(ik4) , intent(in) :: nlon , nlat , jx , iy
    real(rkx) , dimension(nlon,nlat) , intent(in) :: b2
    real(rkx) , dimension(jx,iy) , intent(out) :: b3
    real(rkx) , dimension(jx,iy) , intent(in) :: d1xa , d1xb , d1xc , d1xd
    integer(ik4) , dimension(jx,iy) , intent(in) :: i1dl , i1dr , &
      i1ul , i1ur , j1dl , j1dr , j1ul , j1ur
    real(rkx) :: wa , wb , wc , wd , wg , vv
    real(rkx) , dimension(jx,iy) :: smth1 , smth2
    integer(ik4) :: ifound
    integer(ik4) :: i , j , mdl , mdr , mul , mur , ndl , ndr , nul , nur

    do i = 1 , iy
      do j = 1 , jx
        mur = i1ur(j,i)
        nur = j1ur(j,i)
        mul = i1ul(j,i)
        nul = j1ul(j,i)
        mdr = i1dr(j,i)
        ndr = j1dr(j,i)
        mdl = i1dl(j,i)
        ndl = j1dl(j,i)
        wa = d1xa(j,i)
        wb = d1xb(j,i)
        wc = d1xc(j,i)
        wd = d1xd(j,i)
        ifound = 0
        wg = d_zero
        vv = d_zero
        if ( b2(mur,nur) > missc ) then
          ifound = ifound + 1
          vv = vv + b2(mur,nur)*wa
          wg = wg + wa
        end if
        if ( b2(mul,nul) > missc ) then
          ifound = ifound + 1
          vv = vv + b2(mul,nul)*wb
          wg = wg + wb
        end if
        if ( b2(mdr,ndr) > missc ) then
          ifound = ifound + 1
          vv = vv + b2(mdr,ndr)*wc
          wg = wg + wc
        end if
        if ( b2(mdl,ndl) > missc ) then
          ifound = ifound + 1
          vv = vv + b2(mdl,ndl)*wd
          wg = wg + wd
        end if
        if ( ifound /= 0 ) then
          b3(j,i) = vv/wg
        else
         b3(j,i) = missl
        end if
      end do
    end do
    ! Smooth the field
    smth1(:,:) = b3(:,:)
    smth2(:,:) = b3(:,:)
    do i = 1 , iy
      do j = 2 , jx - 1
        if ( b3(j,i) > missl .and. b3(j+1,i) > missl .and. &
             b3(j-1,i) > missl ) then
          smth2(j,i) = d_rfour*(d_two*smth1(j,i)+smth1(j+1,i)+smth1(j-1,i))
        end if
      end do
    end do
    do i = 2 , iy - 1
      do j = 1 , jx
        if ( b3(j,i) > missl .and. b3(j,i+1) > missl .and. &
             b3(j,i-1) > missl ) then
          smth1(j,i) = d_rfour*(d_two*smth2(j,i)+smth2(j,i+1)+smth2(j,i-1))
        end if
      end do
    end do
    b3(:,:) = smth1(:,:)
  end subroutine dwgt

  subroutine distwgtcr(b3,b2,alon,alat,glon,glat,jx,iy,nlon,nlat)
    implicit none
    integer(ik4) :: iy , jx , nlat , nlon
    real(rkx) , dimension(jx,iy) :: alat , alon
    real(rkx) , dimension(jx,iy) :: b3
    real(rkx) , dimension(nlon,nlat) :: glat , glon
    real(rkx) , dimension(nlon,nlat) :: b2
    intent (in) alat , alon , b2 , glat , glon , iy , jx , nlat , nlon
    intent (out) b3
    !
    ! FIND THE FOUR CLOSEST POINTS TO THE GRID WE WANT TO HAVE VALUE,
    ! THEN DO THE AVERAGE OF THOSE FOUR POINTS WEIGHTED BY THE DISTANCE.
    ! THE CLOSEST ONE HAS BIGGEST WEIGHT.
    !
    ! B2(JX,IX,NLEV) IS THE INPUT FIELD ON PREVIOUS regular or
    ! irregular GRID. B3(JX,IX,NLEV) IS THE OUTPUT FIELD ON new
    ! (regular or irregular) GRID. GLON......LONGITUDE VALUES IN
    ! DEGREES OF THE INTERMEDIATE GRID4. GLAT......LATITUDE VALUES IN
    ! DEGREES OF THE INTERMEDIATE GRID4.
    !
    ! Find the Minimum and Maximum of GLON, GLAT, ALON and ALAT
    !
    if ( imxmn == 0 ) then
      glonmn = minval(glon(1,:))
      glonmx = maxval(glon(nlon,:))
      alonmx = maxval(alon)
      alonmn = minval(alon)
      glatmx = maxval(glat(:,1))
      glatmn = minval(glat(:,nlat))
      alatmx = maxval(alat)
      alatmn = minval(alat)
      if ( maxval(glon) - minval(glon) > 350.0_rkx ) lonwrap = .true.
      if ( maxval(glat) - minval(glat) > 170.0_rkx ) latpole = .true.
      write (stdout,*) 'GLONMN,ALONMN,ALONMX,GLONMX = '
      write (stdout,*) glonmn , alonmn , alonmx , glonmx
      write (stdout,*) 'GLATMN,ALATMN,ALATMX,GLATMX = '
      write (stdout,*) glatmn , alatmn , alatmx , glatmx
      imxmn = 1
    end if
    if ( lcross == 0 ) then
      write (stdout,*) 'FIRST TIME in CRESSMCR'
      write (stdout,*) 'Calculating weights.... (will take long time)'
      if (.not. associated(ic1dl)) call getmem2d(ic1dl,1,jx,1,iy,'ic1dl')
      if (.not. associated(ic1dr)) call getmem2d(ic1dr,1,jx,1,iy,'ic1dr')
      if (.not. associated(ic1ul)) call getmem2d(ic1ul,1,jx,1,iy,'ic1ul')
      if (.not. associated(ic1ur)) call getmem2d(ic1ur,1,jx,1,iy,'ic1ur')
      if (.not. associated(jc1dl)) call getmem2d(jc1dl,1,jx,1,iy,'jc1dl')
      if (.not. associated(jc1dr)) call getmem2d(jc1dr,1,jx,1,iy,'jc1dr')
      if (.not. associated(jc1ul)) call getmem2d(jc1ul,1,jx,1,iy,'jc1ul')
      if (.not. associated(jc1ur)) call getmem2d(jc1ur,1,jx,1,iy,'jc1ur')
      if (.not. associated(dc1xa)) call getmem2d(dc1xa,1,jx,1,iy,'dc1xa')
      if (.not. associated(dc1xb)) call getmem2d(dc1xb,1,jx,1,iy,'dc1xb')
      if (.not. associated(dc1xc)) call getmem2d(dc1xc,1,jx,1,iy,'dc1xc')
      if (.not. associated(dc1xd)) call getmem2d(dc1xd,1,jx,1,iy,'dc1xd')
      call compwgt(alon,alat,glon,glat,dc1xa,dc1xb,dc1xc,dc1xd, &
                   ic1dl,ic1dr,ic1ul,ic1ur,jc1dl,jc1dr,jc1ul,jc1ur, &
                   jx,iy,nlon,nlat)
      write (stdout,*) 'Done.'
      lcross = 1
    end if
    call dwgt(jx,iy,nlon,nlat,b2,b3,dc1xa,dc1xb,dc1xc,dc1xd, &
              ic1dl,ic1dr,ic1ul,ic1ur,jc1dl,jc1dr,jc1ul,jc1ur)
  end subroutine distwgtcr

  subroutine distwgtdt(b3,b2,alon,alat,glon,glat,jx,iy,nlon,nlat)
    implicit none
    integer(ik4) :: iy , jx , nlat , nlon
    real(rkx) , dimension(jx,iy) :: alat , alon
    real(rkx) , dimension(jx,iy) :: b3
    real(rkx) , dimension(nlon,nlat) :: glat , glon
    real(rkx) , dimension(nlon,nlat) :: b2
    intent (in) alat , alon , b2 , glat , glon , iy , jx , nlat , nlon
    intent (out) b3
    !
    ! FIND THE FOUR CLOSEST POINTS TO THE GRID WE WANT TO HAVE VALUE,
    ! THEN DO THE AVERAGE OF THOSE FOUR POINTS WEIGHTED BY THE DISTANCE.
    ! THE CLOSEST ONE HAS BIGGEST WEIGHT.
    !
    ! B2(JX,IX,NLEV) IS THE INPUT FIELD ON PREVIOUS regular or
    ! irregular GRID. B3(JX,IX,NLEV) IS THE OUTPUT FIELD ON new
    ! (regular or irregular) GRID. GLON......LONGITUDE VALUES IN
    ! DEGREES OF THE INTERMEDIATE GRID4. GLAT......LATITUDE VALUES IN
    ! DEGREES OF THE INTERMEDIATE GRID4.
    !
    ! Find the Minimum and Maximum of GLON, GLAT, ALON and ALAT
    !
    if ( imxmn == 0 ) then
      glonmn = minval(glon(1,:))
      glonmx = maxval(glon(nlon,:))
      alonmx = maxval(alon)
      alonmn = minval(alon)
      glatmx = maxval(glat(:,1))
      glatmn = minval(glat(:,nlat))
      alatmx = maxval(alat)
      alatmn = minval(alat)
      write (stdout,*) 'GLONMN,ALONMN,ALONMX,GLONMX= '
      write (stdout,*) glonmn , alonmn , alonmx , glonmx
      write (stdout,*) 'GLATMN,ALATMN,ALATMX,GLATMX= '
      write (stdout,*) glatmn , alatmn , alatmx , glatmx
      if ( glonmx - glonmn > 350.0_rkx ) lonwrap = .true.
      if ( glatmx - glatmn > 170.0_rkx ) latpole = .true.
      imxmn = 1
    end if
    if ( ldot == 0 ) then
      write (stdout,*) 'FIRST TIME in CRESSMDT'
      write (stdout,*) 'Calculating weights.... (will take long time)'
      if (.not. associated(id1dl)) call getmem2d(id1dl,1,jx,1,iy,'id1dl')
      if (.not. associated(id1dr)) call getmem2d(id1dr,1,jx,1,iy,'id1dr')
      if (.not. associated(id1ul)) call getmem2d(id1ul,1,jx,1,iy,'id1ul')
      if (.not. associated(id1ur)) call getmem2d(id1ur,1,jx,1,iy,'id1ur')
      if (.not. associated(jd1dl)) call getmem2d(jd1dl,1,jx,1,iy,'jd1dl')
      if (.not. associated(jd1dr)) call getmem2d(jd1dr,1,jx,1,iy,'jd1dr')
      if (.not. associated(jd1ul)) call getmem2d(jd1ul,1,jx,1,iy,'jd1ul')
      if (.not. associated(jd1ur)) call getmem2d(jd1ur,1,jx,1,iy,'jd1ur')
      if (.not. associated(dd1xa)) call getmem2d(dd1xa,1,jx,1,iy,'dd1xa')
      if (.not. associated(dd1xb)) call getmem2d(dd1xb,1,jx,1,iy,'dd1xb')
      if (.not. associated(dd1xc)) call getmem2d(dd1xc,1,jx,1,iy,'dd1xc')
      if (.not. associated(dd1xd)) call getmem2d(dd1xd,1,jx,1,iy,'dd1xd')
      call compwgt(alon,alat,glon,glat,dd1xa,dd1xb,dd1xc,dd1xd, &
                   id1dl,id1dr,id1ul,id1ur,jd1dl,jd1dr,jd1ul,jd1ur, &
                   jx,iy,nlon,nlat)
      write (stdout,*) 'Done.'
      ldot = 1
    end if
    call dwgt(jx,iy,nlon,nlat,b2,b3,dd1xa,dd1xb,dd1xc,dd1xd, &
              id1dl,id1dr,id1ul,id1ur,jd1dl,jd1dr,jd1ul,jd1ur)
  end subroutine distwgtdt

  subroutine cressmcr3d(b3,b2,alon,alat,glon,glat,jx,iy,nlon,nlat,nlev,nf)
    implicit none
    integer(ik4) :: iy , jx , nlat , nlev , nlon , nf
    real(rkx) , dimension(jx,iy) :: alat , alon
    real(rkx) , dimension(jx,iy,nlev*nf) :: b3
    real(rkx) , dimension(nlon,nlat) :: glat , glon
    real(rkx) , dimension(nlon,nlat,nlev*nf) :: b2
    intent (in) alat , alon , b2 , glat , glon , iy , jx , nlat ,     &
                nlev , nlon , nf
    intent (out) b3
    integer(ik4) :: k , l , kin
    do l = 1 , nf
      do k = 1 , nlev
         kin = (l-1)*nlev+k
         call distwgtcr(b3(:,:,kin),b2(:,:,kin),alon,alat, &
                           glon,glat,jx,iy,nlon,nlat)
      end do
    end do
  end subroutine cressmcr3d

  subroutine cressmcr2d(b3,b2,alon,alat,glon,glat,jx,iy,nlon,nlat)
    implicit none
    integer(ik4) :: iy , jx , nlat , nlon
    real(rkx) , dimension(jx,iy) :: alat , alon
    real(rkx) , dimension(jx,iy) :: b3
    real(rkx) , dimension(nlon,nlat) :: glat , glon
    real(rkx) , dimension(nlon,nlat) :: b2
    intent (in) alat , alon , b2 , glat , glon , iy , jx , nlat , nlon
    intent (out) b3

    call distwgtcr(b3(:,:),b2(:,:),alon,alat,glon,glat,jx,iy,nlon,nlat)
  end subroutine cressmcr2d

  subroutine cressmdt(b3,b2,alon,alat,glon,glat,jx,iy,nlon,nlat,nlev,nf)
    implicit none
    integer(ik4) :: iy , jx , nlat , nlev , nlon , nf
    real(rkx) , dimension(jx,iy) :: alat , alon
    real(rkx) , dimension(jx,iy,nlev*nf) :: b3
    real(rkx) , dimension(nlon,nlat) :: glat , glon
    real(rkx) , dimension(nlon,nlat,nlev*nf) :: b2
    intent (in) alat , alon , b2 , glat , glon , iy , jx , nlat ,     &
                nlev , nlon , nf
    intent (out) b3
    integer(ik4) :: k , l , kin
    do l = 1 , nf
      do k = 1 , nlev
         kin = (l-1)*nlev+k
         call distwgtdt(b3(:,:,kin),b2(:,:,kin),alon,alat, &
                           glon,glat,jx,iy,nlon,nlat)
      end do
    end do
  end subroutine cressmdt

  subroutine kernsmooth2(f,nx,ny,npass)
    implicit none
    integer(ik4) , intent(in) :: nx , ny , npass
    real(rkx) , dimension(nx,ny) , intent(inout) :: f
    integer(ik4) :: i , j , n
    real(rkx) , dimension(nx,ny) :: newf
    do n = 1 , npass
      do j = 2 , ny - 1
        do i = 2 , nx - 1
          newf(i,j) = (f(i+1,j-1) + f(i+1,j) + f(i+1,j+1) + &
            f(i,j-1) + f(i,j) * 4.0_rkx + f(i,j+1) + &
            f(i-1,j-1) + f(i-1,j) + f(i-1,j+1)) / 12.0_rkx
        end do
      end do
      do j = 2 , ny-1
        newf(1,j) = (f(1,j-1) + f(1,j) * 7.0_rkx + f(1,j+1) + &
          f(2,j-1) + f(2,j) + f(2,j+1)) / 12.0_rkx
        newf(nx,j) = (f(nx,j-1) + f(nx,j) * 7.0_rkx + f(nx,j+1) + &
          f(nx-1,j-1) + f(nx-1,j) + f(nx-1,j+1)) / 12.0_rkx
      end do
      do i = 2 , nx-1
        newf(i,1) = (f(i-1,1) + f(i,1) * 7.0_rkx + f(i+1,1) + &
          f(i-1,2) + f(i,2) + f(i+1,2)) / 12.0_rkx
        newf(i,ny) = (f(i-1,ny) + f(i,ny) * 7.0_rkx + f(i+1,ny) + &
          f(i-1,ny-1) + f(i,ny-1) + f(i+1,ny-1)) / 12.0_rkx
      end do
      newf(1,1) = f(1,1)
      newf(nx,1) = f(nx,1)
      newf(1,ny) = f(1,ny)
      newf(nx,ny) = f(nx,ny)
      f(:,:) = newf(:,:)
    end do
  end subroutine kernsmooth2

  subroutine kernsmooth3(f,nx,ny,nz,npass)
    implicit none
    integer(ik4) , intent(in) :: nx , ny , nz , npass
    real(rkx) , dimension(nx,ny,nz) , intent(inout) :: f
    integer(ik4) :: i , j , k , n
    real(rkx) , dimension(nx,ny) :: newf
    do n = 1 , npass
      do k = 1 , nz
        do j = 2 , ny - 1
          do i = 2 , nx - 1
            newf(i,j) = (f(i+1,j-1,k) + f(i+1,j,k) + f(i+1,j+1,k) + &
              f(i,j-1,k) + f(i,j,k) * 4.0_rkx + f(i,j+1,k) + &
              f(i-1,j-1,k) + f(i-1,j,k) + f(i-1,j+1,k)) / 12.0_rkx
          end do
        end do
        do j = 2 , ny-1
          newf(1,j) = (f(1,j-1,k) + f(1,j,k) * 7.0_rkx + f(1,j+1,k) + &
            f(2,j-1,k) + f(2,j,k) + f(2,j+1,k)) / 12.0_rkx
          newf(nx,j) = (f(nx,j-1,k) + f(nx,j,k) * 7.0_rkx + f(nx,j+1,k) + &
            f(nx-1,j-1,k) + f(nx-1,j,k) + f(nx-1,j+1,k)) / 12.0_rkx
        end do
        do i = 2 , nx-1
          newf(i,1) = (f(i-1,1,k) + f(i,1,k) * 7.0_rkx + f(i+1,1,k) + &
            f(i-1,2,k) + f(i,2,k) + f(i+1,2,k)) / 12.0_rkx
          newf(i,ny) = (f(i-1,ny,k) + f(i,ny,k) * 7.0_rkx + f(i+1,ny,k) + &
            f(i-1,ny-1,k) + f(i,ny-1,k) + f(i+1,ny-1,k)) / 12.0_rkx
        end do
        newf(1,1) = f(1,1,k)
        newf(nx,1) = f(nx,1,k)
        newf(1,ny) = f(1,ny,k)
        newf(nx,ny) = f(nx,ny,k)
        f(:,:,k) = newf(:,:)
      end do
    end do
  end subroutine kernsmooth3

  integer(ik4) function whereislon(nlon,lon,lonarr) result(jj)
    implicit none
    integer(ik4) , intent(in) :: nlon
    real(rkx) , intent(in) :: lon
    real(rkx) , dimension(nlon) , intent(in) :: lonarr
    real(rkx) , dimension(nlon) :: xlonarr
    real(rkx) :: xlon , dlon

    xlonarr(:) = lonarr(:)
    xlon = lon
    if ( xlon < 0.0_rkx ) xlon = xlon + 360.0_rkx
    if ( xlonarr(1) > 0.0_rkx .and. xlonarr(nlon) < 0.0_rkx ) then
      ! window crossing timeline
      xlonarr = mod(xlonarr+360.0_rkx,360.0_rkx)
      xlon = mod(xlon+360.0_rkx,360.0_rkx)
    else if ( xlonarr(1) > xlonarr(nlon) ) then
      ! window crossing greenwich
      where ( xlonarr > 180.0_rkx )
        xlonarr = xlonarr - 360.0_rkx
      end where
      if ( xlon > 180.0_rkx ) xlon = xlon - 360.0_rkx
    end if
    dlon = (xlonarr(nlon)-xlonarr(1))/real(nlon-1,rkx)
    jj = int((xlon-xlonarr(1))/dlon) + 1
  end function whereislon

  integer(ik4) function whereislat(nlat,lat,latarr) result(ii)
    implicit none
    real(rkx) , intent(in) :: lat
    integer(ik4) , intent(in) :: nlat
    real(rkx) , dimension(nlat) , intent(in) :: latarr
    real(rkx) , dimension(nlat) :: xlatarr
    real(rkx) :: xlat , dlat

    xlatarr(:) = latarr(:)
    xlat = lat
    dlat = (xlatarr(nlat)-xlatarr(1))/real(nlat-1,rkx)
    if ( dlat > 0 ) then
      ii = int((xlat-xlatarr(1))/dlat) + 1
    else
      ii = int((xlatarr(1)-xlat)/abs(dlat)) + 1
    end if
    if ( ii < 1 ) ii = 1
    if ( ii > nlat - 1 ) ii = nlat - 1
  end function whereislat

end module mod_interp

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
