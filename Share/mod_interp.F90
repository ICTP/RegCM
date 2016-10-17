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

  implicit none

  private

  public :: bilinx , cressmcr , cressmdt , distwgtcr , distwgtdt
  public :: gcdist , kernsmooth
  public :: global_domain , get_window

  real(rkx) :: alatmn , alatmx , alonmn , alonmx
  real(rkx) :: glatmn , glatmx , glonmn , glonmx

  integer(ik4) :: imxmn = 0
  integer(ik4) :: lcross = 0
  integer(ik4) :: ldot = 0

  real(rkx) , parameter :: deg720 = d_two*deg360
  real(rkx) , parameter :: missl = -9999.0_rkx
  real(rkx) , parameter :: missc = -9990.0_rkx
  real(rkx) , parameter :: mindist = 1.0e-6_rkx
  integer(ik4) , parameter :: p_factor = 2

  real(rkx) , pointer , dimension(:,:) :: dc1xa , dc1xb , dc1xc , dc1xd
  real(rkx) , pointer , dimension(:,:) :: dd1xa , dd1xb , dd1xc , dd1xd
  integer(ik4) , pointer , dimension(:,:) :: ic1dl , ic1dr , ic1ul , ic1ur , &
                                       jc1dl , jc1dr , jc1ul , jc1ur
  integer(ik4) , pointer , dimension(:,:) :: id1dl , id1dr , id1ul , id1ur , &
                                       jd1dl , jd1dr , jd1ul , jd1ur
  logical :: lonwrap = .false.
  logical :: latpole = .false.

  type global_domain
    integer(ik4) :: global_ni
    integer(ik4) :: global_nj
    integer(ik4) :: ntiles
    integer(ik4) , dimension(2) :: ni
    integer(ik4) , dimension(2) :: igstart
    integer(ik4) , dimension(2) :: igstop
    integer(ik4) :: nj
    integer(ik4) :: jgstart
    integer(ik4) :: jgstop
  end type global_domain

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
          p2 = dlon - p1
          q1 = mod(alat(j,i)-hlat(i1),dlat)
          if ( q1 < 0.0_rkx ) q1 = dlat + q1
          q2 = dlat - q1
          if ( i2 > i1 ) then
            if ( j2 > j1 ) then
              b3(j,i) = ((b2(j1,i1)*p2+b2(j2,i1)*p1)*q1 + &
                         (b2(j1,i2)*p2+b2(j2,i2)*p1)*q2)/(dlon*dlat)
            else
              b3(j,i) = ((b2(j1,i1)*p1+b2(j2,i1)*p2)*q1 + &
                         (b2(j1,i2)*p1+b2(j2,i2)*p2)*q2)/(dlon*dlat)
            end if
          else
            if ( j2 > j1 ) then
              b3(j,i) = ((b2(j1,i1)*p2+b2(j2,i1)*p1)*q2 + &
                         (b2(j1,i2)*p2+b2(j2,i2)*p1)*q1)/(dlon*dlat)
            else
              b3(j,i) = ((b2(j1,i1)*p1+b2(j2,i1)*p2)*q2 + &
                         (b2(j1,i2)*p1+b2(j2,i2)*p2)*q1)/(dlon*dlat)
            end if
          end if
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
            dist = gcdist(glat(m,n),glon(m,n),alat(j,i),alon(j,i))
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
        wa = gcdist(glat(mur,nur),glon(mur,nur),alat(j,i),alon(j,i))
        wb = gcdist(glat(mul,nul),glon(mul,nul),alat(j,i),alon(j,i))
        wc = gcdist(glat(mdr,ndr),glon(mdr,ndr),alat(j,i),alon(j,i))
        wd = gcdist(glat(mdl,ndl),glon(mdl,ndl),alat(j,i),alon(j,i))
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
      if ( glonmx - glonmn > 350.0_rkx ) lonwrap = .true.
      if ( glatmx - glatmn > 170.0_rkx ) latpole = .true.
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

  real(rkx) function gcdist_simple(lat1,lon1,lat2,lon2)
    implicit none
    real(rkx) , intent(in) :: lat1 , lon1 , lat2, lon2
    real(rkx) :: clat1 , slat1 , clat2 , slat2 , cdlon , crd
    clat1 = cos(lat1*degrad)
    slat1 = sin(lat1*degrad)
    clat2 = cos(lat2*degrad)
    slat2 = sin(lat2*degrad)
    cdlon = cos((lon1-lon2)*degrad)
    crd   = slat1*slat2+clat1*clat2*cdlon
    ! Have it in km to avoid numerical problems :)
    gcdist_simple = erkm*acos(crd)
  end function gcdist_simple

  real(rkx) function gcdist(lat1,lon1,lat2,lon2)
    implicit none
    real(rkx) , intent(in) :: lat1 , lon1 , lat2, lon2
    real(rkx) :: clat1 , slat1 , clat2 , slat2 , cdlon , sdlon
    real(rkx) :: y , x
    clat1 = cos(lat1*degrad)
    slat1 = sin(lat1*degrad)
    clat2 = cos(lat2*degrad)
    slat2 = sin(lat2*degrad)
    cdlon = cos((lon1-lon2)*degrad)
    sdlon = sin((lon1-lon2)*degrad)
    y = sqrt((clat2*sdlon)**2+(clat1*slat2-slat1*clat2*cdlon)**2)
    x = slat1*slat2+clat1*clat2*cdlon
    ! Have it in km to avoid numerical problems :)
    gcdist = max(erkm*atan2(y,x)/ds,mindist)
  end function gcdist

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

  ! Assumptions :
  !  GLAT is regular ordered array of latitudes in a GLOBAL grid
  !  GLON is regular ordered array of longitudes in a GLOBAL grid
  !  XLAT is local grid latitudes going South to North
  !  XLON is local grid longitudes  going West to East
  ! ALL LONGITUDES ARE in the range -180.0 <-> 180.0
  subroutine get_window(glat,glon,xlat,xlon,i_band,domain)
    implicit none
    real(rkx) , dimension(:) , intent(in) :: glat
    real(rkx) , dimension(:) , intent(in) :: glon
    real(rkx) , dimension(:,:) , intent(in) :: xlat
    real(rkx) , dimension(:,:) , intent(in) :: xlon
    integer(ik4) , intent(in) :: i_band
    type(global_domain) , intent(out) :: domain

    real(rkx) :: dlat , dlon
    real(rkx) , allocatable , dimension(:,:) :: xlon360
    real(rkx) :: maxlat
    real(rkx) :: minlat
    real(rkx) :: maxlon
    real(rkx) :: minlon
    integer :: gi , gj , xi , xj , l1 , l2 , i , j

    xi = size(xlon,1)
    xj = size(xlat,2)
    maxlat = maxval(xlat)
    minlat = minval(xlat)
    maxlon = maxval(xlon)
    minlon = minval(xlon)

    gi = size(glon)
    gj = size(glat)
    dlat = abs(glat(2) - glat(1))
    dlon = abs(glon(2) - glon(1))
    domain%global_ni = gi
    domain%global_nj = gj

    if ( i_band == 1 ) then
      domain%ntiles = 1
      domain%igstart(1) = 1
      domain%igstop(1) = gi
      domain%igstart(2) = 0
      domain%igstop(2) = 0
    else if ( glon(gi) > 350.0_rkx ) then
      ! Input data is     0 : 360 , xlon is -180 : 180
      allocate(xlon360(xi,xj))
      xlon360 = xlon
      where ( xlon < 0.0_rkx )
        xlon360 = 360.0_rkx + xlon
      end where
      if ( minval(xlon360(1,:)) < maxval(xlon360(xi,:)) ) then
        domain%ntiles = 1
        minlon = minval(xlon360)
        maxlon = maxval(xlon360)
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = int((maxlon-glon(1))/dlon) + 3
        domain%igstart(2) = 0
        domain%igstop(2) = 0
      else
        ! Cross Greenwich line
        minlon = minval(xlon360(1,:))
        maxlon = maxval(xlon(xi,:))
        domain%ntiles = 2
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = gi
        domain%igstart(2) = 1
        domain%igstop(2) = int((maxlon-glon(1))/dlon) + 3
      end if
      deallocate(xlon360)
    else
      ! Input Data is -180 : 180 , xlon is -180 : 180
      if ( xlon(1,xj)   <= xlon(xi,xj)   .and. &
           xlon(1,xj/2) <= xlon(xi,xj/2) .and. &
           xlon(1,1)    <= xlon(xi,1) ) then
        ! it is not crossing timeline
        domain%ntiles = 1
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = int((maxlon-glon(1))/dlon) + 3
        domain%igstart(2) = 0
        domain%igstop(2) = 0
      else
        domain%ntiles = 2
        minlon = 180.0_rkx
        do j = 1 , xj
          if ( xlon(1,j) > 0.0_rkx ) minlon = min(minlon,xlon(1,j))
        end do
        maxlon = -180.0_rkx
        do j = 1 , xj
          if ( xlon(xi,j) < 0.0_rkx ) maxlon = max(maxlon,xlon(xi,j))
        end do
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = gi
        domain%igstart(2) = 1
        domain%igstop(2) = int((maxlon-glon(1))/dlon) + 3
      end if
    end if
    if ( has_north_pole(xlat,xi/2) ) then
      ! North pole inside
      write(stdout,*) 'Correcting North Pole'
      l1 = int((minval(xlat(:,1))-glat(1))/dlat) - 2
      l2 = int((minval(xlat(:,xj))-glat(1))/dlat) - 2
      domain%jgstart = min(l1,l2)
      domain%jgstop = gj
      domain%ntiles = 1
      domain%igstart(1) = 1
      domain%igstop(1) = gi
      domain%igstart(2) = 0
      domain%igstop(2) = 0
    else if ( has_south_pole(xlat,xi/2) ) then
      ! South Pole inside
      write(stdout,*) 'Correcting South Pole'
      l1 = int((maxval(xlat(:,1))-glat(1))/dlat) + 3
      l2 = int((maxval(xlat(:,xj))-glat(1))/dlat) + 3
      domain%jgstart = 1
      domain%ntiles = 1
      domain%jgstop = max(l1,l2)
      domain%igstart(1) = 1
      domain%igstop(1) = gi
      domain%igstart(2) = 0
      domain%igstop(2) = 0
    else
      if ( glat(1) < glat(gj) ) then
        domain%jgstart = int((minlat-glat(1))/dlat) - 2
        domain%jgstop  = int((maxlat-glat(1))/dlat) + 3
      else
        domain%jgstart = int((minlat-glat(gj))/dlat) - 2
        domain%jgstop  = int((maxlat-glat(gj))/dlat) + 3
        domain%nj =  domain%jgstop - domain%jgstart + 1
        domain%jgstart = (gj+1)-(domain%jgstart+domain%nj-1)
        domain%jgstop = domain%jgstart + domain%nj - 1
      end if
    end if

    domain%ni = 0
    do i = 1 , domain%ntiles
      domain%igstart(i) = min(max(1,domain%igstart(i)),gi)
      domain%igstop(i) = max(min(domain%igstop(i),gi),1)
      domain%ni(i) =  domain%igstop(i) - domain%igstart(i) + 1
    end do
    domain%jgstart = min(max(1,domain%jgstart),gj)
    domain%jgstop = max(min(domain%jgstop,gj),1)
    domain%nj =  domain%jgstop - domain%jgstart + 1

    contains

      logical function has_north_pole(l,i)
        real(rkx) , intent(in) , dimension(:,:) :: l
        integer(ik4) , intent(in) :: i
        integer(ik4) :: j
        has_north_pole = .false.
        do j = 2 , size(l,2)
          if ( l(i,j) < l(i,j-1) ) then
            has_north_pole = .true.
            exit
          end if
        end do
      end function has_north_pole

      logical function has_south_pole(l,i)
        real(rkx) , intent(in) , dimension(:,:) :: l
        integer(ik4) , intent(in) :: i
        integer(ik4) :: j
        has_south_pole = .false.
        do j = 2 , size(l,2)
          if ( l(i,j) < l(i,j-1) ) then
            has_south_pole = .true.
            exit
          end if
        end do
      end function has_south_pole

  end subroutine get_window

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
    if ( xlonarr(1) > xlon ) then
       jj = 1
       return
    else if ( xlonarr(nlon) < xlon ) then
      jj = nlon
      return
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
      ii = int((xlat-xlatarr(1))/dlat)
    else
      ii = int((xlatarr(1)-xlat)/abs(dlat)) + 1
    end if
    if ( ii < 1 ) ii = 1
    if ( ii > nlat - 1 ) ii = nlat - 1
  end function whereislat

end module mod_interp

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
