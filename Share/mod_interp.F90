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

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_stdio
  use mod_message
  use mod_memutil
  use mod_constants

  implicit none

  private

  public :: bilinx , bilinx2 , cressmcr , cressmdt , distwgtcr , distwgtdt
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
  real(rkx) , parameter :: mindist = 0.001_rkx ! 1 meter
  real(rkx) , parameter :: p_factor = 2.00_rkx

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

  interface bilinx2
    module procedure bilinx2_2d
    module procedure bilinx2_3d
  end interface bilinx2

  interface kernsmooth
    module procedure kernsmooth2
    module procedure kernsmooth3
  end interface kernsmooth

  contains

  subroutine bilinx(fin,fout,lono,lato,loni,lati,nloni,nlati,jx,iy,nflds)
    implicit none
    integer(ik4) :: jx , iy , nflds , nlati , nloni
    real(rkx) , dimension(nloni,nlati,nflds) :: fin
    real(rkx) , dimension(nlati) :: lati
    real(rkx) , dimension(jx,iy) :: lato , lono
    real(rkx) , dimension(nloni) :: loni
    real(rkx) , dimension(jx,iy,nflds) :: fout
    intent (in) fin , jx , iy , lati , lato , loni , lono , nflds ,   &
                nlati , nloni
    intent (out) fout

    real(rkx) :: bas , lon360 , p , q , xsum , xind , yind
    integer(ik4) :: i , ip , ipp1 , j , jq , jqp1 , l
    logical :: lg
    !
    ! PERFORMING BI-LINEAR INTERPOLATION USING 4 GRID POINTS FROM A
    ! BIGGER RECTANGULAR GRID TO A GRID DESCRIBED BY XLONS AND XLATS OF
    ! GRID2. A POINT ON GRID2 IS TRAPPED WITHIN FOUR GRID POINTS ON
    ! GRID4.THE GRID POINTS ARE ALWAYS TO THE NORTH AND EAST OF THE
    ! TRAPPED POINT.. THE ALGORITHM COMPUTES THE FRACTIONAL DISTANCES
    ! IN BOTH X AND Y DIRECTION OF THE TRAPPED GRID POINT AND USES THE
    ! INFORMATION AS WEIGHTING FACTORS IN THE INTERPOLATION.
    ! THERE IS ONE LESS ROW AND COLUMN WHEN THE SCALAR FIELDS ARE
    ! INTERPOLATED BECAUSE XLATS AND XLONS ARE NOT DEFINED FOR
    ! THE CROSS POINTS IN THE MM4 MODEL.
    !
    ! IN(NLONI,NLATI,NFLDS)  IS THE INPUT FIELD ON REGULAR LAT/LON GRID.
    ! OUT(NLATO,NLONO,NFLDS) IS THE OUTPUT FIELD ON LAMBERT CONFORMAL
    ! GRID. LONI.....LONGITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
    ! LATI.....LATITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
    ! P.........EAST-WEST WEIGHTING FACTOR.
    ! Q.........NORTH-SOUTH WEIGHTING FACTOR.
    ! IP........GRID POINT LOCATION IN EAST-WEST OF TRAPPED GRID POINT.
    ! IQ........GRID POINT LOCATION IN NORTH-SOUTH OF TRAPPED GRID POINT.
    !
    ! Global dataset ?
    !
    lg = .true.
    do i = 1 , iy
      do j = 1 , jx
        yind = (((lato(j,i)-lati(1))/(lati(nlati)-lati(1))) * &
                   real(nlati-1,rkx))+d_one
        jq = int(yind)
        jq = max(jq,1)
        jqp1 = min(jq+1,nlati)
        q = yind - jq
        lon360 = lono(j,i)
        if ( lono(j,i) < deg00 ) lon360 = lono(j,i) + deg360
        xind = (((lon360-loni(1))/(loni(nloni)-loni(1))) * &
                    real(nloni-1,rkx))+d_one
        if ( xind < d_one .and. lg ) then
          ip = nloni
          ipp1 = 1
          p = xind
        else if ( (xind-nloni) > deg00 .and. lg ) then
          ip = nloni
          ipp1 = 1
          p = xind - nloni
        else
          ip = int(xind)
          ip = max(ip,1)
          ipp1 = min(ip+1,nloni)
          p = xind - ip
        end if
        do l = 1 , nflds
          xsum = d_zero
          bas = d_zero
          if ( fin(ip,jq,l) < missc .and. fin(ipp1,jq,l) < missc .and. &
               fin(ipp1,jqp1,l) < missc .and. fin(ip,jqp1,l) < missc ) then
            fout(j,i,l) = missl
          else
            if ( fin(ip,jq,l) > missc ) then
              xsum = xsum + (d_one-q)*(d_one-p)*fin(ip,jq,l)
              bas = bas + (d_one-q)*(d_one-p)
            end if
            if ( fin(ipp1,jq,l) > missc ) then
              xsum = xsum + (d_one-q)*p*fin(ipp1,jq,l)
              bas = bas + (d_one-q)*p
            end if
            if ( fin(ipp1,jqp1,l) > missc ) then
              xsum = xsum + q*p*fin(ipp1,jqp1,l)
              bas = bas + q*p
            end if
            if ( fin(ip,jqp1,l) > missc ) then
              xsum = xsum + q*(d_one-p)*fin(ip,jqp1,l)
              bas = bas + q*(d_one-p)
            end if
            if ( bas > dlowval ) then
              fout(j,i,l) = xsum/bas
            else
              fout(j,i,l) = missl
            end if
          end if
        end do
      end do
    end do
  end subroutine bilinx

  subroutine bilinx2_3d(b3,b2,alon,alat,hlon,hlat,nlon,nlat,jx,iy,llev)
    implicit none
    integer(ik4) :: iy , jx , llev , nlat , nlon
    real(rkx) , dimension(jx,iy) :: alat , alon
    real(rkx) , dimension(nlon,nlat,llev) :: b2
    real(rkx) , dimension(jx,iy,llev) :: b3
    real(rkx) , dimension(nlat) :: hlat
    real(rkx) , dimension(nlon) :: hlon
    intent (in) alat , alon , b2 , hlat , hlon , iy , jx , llev , nlat , nlon
    intent (out) b3
    real(rkx) :: ave , p1 , p2 , q1 , q2
    integer(ik4) :: i , i1 , i2 , ii , j , j1 , j2 , jj , l
    !
    ! PERFORMING BI-LINEAR INTERPOLATION USING 4 GRID POINTS FROM A
    ! BIGGER RECTANGULAR GRID TO A GRID DESCRIBED BY ALON AND ALAT OF
    ! GRID2. A POINT ON GRID2 IS TRAPPED WITHIN FOUR GRID POINTS ON
    ! GRID4.THE GRID POINTS ARE ALWAYS TO THE NORTH AND EAST OF THE
    ! TRAPPED POINT. THE ALGORITHM COMPUTES THE FRACTIONAL DISTANCES IN
    ! BOTH X AND Y DIRECTION OF THE TRAPPED GRID POINT AND USES THE
    ! INFORMATION AS WEIGHTING FACTORS IN THE INTERPOLATION.
    ! THERE IS ONE LESS ROW AND COLUMN WHEN THE SCALAR FIELDS ARE
    ! INTERPOLATED BECAUSE ALAT AND ALON ARE NOT DEFINED FOR
    ! THE CROSS POINTS IN THE RegCM MODEL.
    !
    ! B2(JX,IX,NLEV) IS THE INPUT FIELD ON REGULAR LAT/LON GRID.
    ! B3(JX,IX,NLEV) IS THE OUTPUT FIELD ON LAMBERT CONFORMAL GRID.
    ! HLON......LONGITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
    ! HLAT......LATITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
    ! P.........EAST-WEST WEIGHTING FACTOR.
    ! Q.........NORTH-SOUTH WEIGHTING FACTOR.
    ! IP........GRID POINT LOCATION IN EAST-WEST OF TRAPPED GRID POINT.
    ! IQ........GRID POINT LOCATION IN NORTH-SOUTH OF TRAPPED GRID POINT.
    !
    i1 = 0
    i2 = 0
    j1 = 0
    j2 = 0
    q1 = d_zero
    q2 = d_zero
    p1 = d_zero
    p2 = d_zero
    do i = 1 , iy
      do j = 1 , jx
        i1 = 1000
        do ii = 1 , nlon - 1
          if ( alon(j,i) >= hlon(ii) .and. alon(j,i) < hlon(ii+1) ) then
            p1 = alon(j,i) - hlon(ii)
            p2 = hlon(ii+1) - alon(j,i)
            i1 = ii
            i2 = ii + 1
            exit
          else if ( alon(j,i) >= hlon(ii)-deg360 .and. &
                    alon(j,i) < hlon(ii+1)-deg360 ) then
            p1 = alon(j,i) - (hlon(ii)-deg360)
            p2 = (hlon(ii+1)-deg360) - alon(j,i)
            i1 = ii
            i2 = ii + 1
            exit
          else if ( alon(j,i) >= hlon(ii)+deg360 .and. &
                    alon(j,i) < hlon(ii+1)+deg360 ) then
            p1 = alon(j,i) - (hlon(ii)+deg360)
            p2 = (hlon(ii+1)+deg360) - alon(j,i)
            i1 = ii
            i2 = ii + 1
            exit
          end if
        end do
        if ( alon(j,i) >= hlon(nlon) .and. alon(j,i) < hlon(1)+deg360 ) then
          p1 = alon(j,i) - hlon(nlon)
          p2 = (hlon(1)+deg360) - alon(j,i)
          i1 = nlon
          i2 = 1
        else if ( alon(j,i) >= hlon(nlon)+deg360 .and. &
                  alon(j,i) < hlon(1)+deg720 ) then
          p1 = alon(j,i) - (hlon(nlon)+deg360)
          p2 = (hlon(1)+deg720) - alon(j,i)
          i1 = nlon
          i2 = 1
        else if ( alon(j,i) >= hlon(nlon)-deg360 .and. &
                  alon(j,i) < hlon(1) ) then
          p1 = alon(j,i) - (hlon(nlon)-deg360)
          p2 = hlon(1) - alon(j,i)
          i1 = nlon
          i2 = 1
        else if ( alon(j,i) >= hlon(nlon)-deg720 .and. &
                  alon(j,i) < hlon(1)-deg360 ) then
          p1 = alon(j,i) - (hlon(nlon)-deg720)
          p2 = (hlon(1)-deg360) - alon(j,i)
          i1 = nlon
          i2 = 1
        end if
        if ( i1 == 1000 ) then
          call die('bilinx2','Could not find the right longitude',1)
        end if
        j1 = 1000
        do jj = 1 , nlat - 1
          if ( alat(j,i) >= hlat(jj) .and. alat(j,i) < hlat(jj+1) ) then
            q1 = alat(j,i) - hlat(jj)
            q2 = hlat(jj+1) - alat(j,i)
            j1 = jj
            j2 = jj + 1
            exit
          else if ( alat(j,i) <= hlat(1) ) then
            q1 = d_one
            q2 = d_one
            j1 = 1
            j2 = 1
            exit
          else if ( alat(j,i) >= hlat(nlat) ) then
            q1 = d_one
            q2 = d_one
            j1 = nlat
            j2 = nlat
          end if
        end do
        if ( j1 == 1000 ) then
          call die('bilinx2','Could not find the right latitude',1)
        end if
        if ( j1 > 0 .and. j1 < nlat ) then
          do l = 1 , llev
            b3(j,i,l) = ((b2(i1,j1,l)*p2+b2(i2,j1,l)*p1)*q2+(b2(i1,j2,l)* &
                          p2+b2(i2,j2,l)*p1)*q1)/(p1+p2)/(q1+q2)
          end do
        else if ( j1 == 0 ) then
          do l = 1 , llev
            ave = d_zero
            do ii = 1 , nlon
              ave = ave + b2(ii,1,l)
            end do
            ave = ave/real(nlon,rkx)
            b3(j,i,l) = ((ave*(p1+p2))*q2+(b2(i1,j2,l)*p2+b2(i2,j2,l)* &
                                       p1)*q1)/(p1+p2)/(q1+q2)
          end do
        else if ( j1 == nlat ) then
          do l = 1 , llev
            ave = d_zero
            do ii = 1 , nlon
              ave = ave + b2(ii,nlat,l)
            end do
            ave = ave/real(nlon,rkx)
            b3(j,i,l) = ((b2(i1,j1,l)*p2+b2(i2,j1,l)*p1)*q2+ &
                         (ave*(p1+p2))*q1)/(p1+p2)/(q1+q2)
          end do
        else
        end if
      end do
    end do
  end subroutine bilinx2_3d

  subroutine bilinx2_2d(b3,b2,alon,alat,hlon,hlat,nlon,nlat,jx,iy)
    implicit none
    integer(ik4) :: iy , jx , nlat , nlon
    real(rkx) , dimension(jx,iy) :: alat , alon
    real(rkx) , dimension(nlon,nlat) :: b2
    real(rkx) , dimension(jx,iy) :: b3
    real(rkx) , dimension(nlat) :: hlat
    real(rkx) , dimension(nlon) :: hlon
    intent (in) alat , alon , b2 , hlat , hlon , iy , jx , nlat , nlon
    intent (out) b3
    real(rkx) :: ave , p1 , p2 , q1 , q2
    integer(ik4) :: i , i1 , i2 , ii , j , j1 , j2 , jj
    !
    ! PERFORMING BI-LINEAR INTERPOLATION USING 4 GRID POINTS FROM A
    ! BIGGER RECTANGULAR GRID TO A GRID DESCRIBED BY ALON AND ALAT OF
    ! GRID2. A POINT ON GRID2 IS TRAPPED WITHIN FOUR GRID POINTS ON
    ! GRID4.THE GRID POINTS ARE ALWAYS TO THE NORTH AND EAST OF THE
    ! TRAPPED POINT. THE ALGORITHM COMPUTES THE FRACTIONAL DISTANCES IN
    ! BOTH X AND Y DIRECTION OF THE TRAPPED GRID POINT AND USES THE
    ! INFORMATION AS WEIGHTING FACTORS IN THE INTERPOLATION.
    ! THERE IS ONE LESS ROW AND COLUMN WHEN THE SCALAR FIELDS ARE
    ! INTERPOLATED BECAUSE ALAT AND ALON ARE NOT DEFINED FOR
    ! THE CROSS POINTS IN THE RegCM MODEL.
    !
    ! B2(JX,IX,NLEV) IS THE INPUT FIELD ON REGULAR LAT/LON GRID.
    ! B3(JX,IX,NLEV) IS THE OUTPUT FIELD ON LAMBERT CONFORMAL GRID.
    ! HLON......LONGITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
    ! HLAT......LATITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
    ! P.........EAST-WEST WEIGHTING FACTOR.
    ! Q.........NORTH-SOUTH WEIGHTING FACTOR.
    ! IP........GRID POINT LOCATION IN EAST-WEST OF TRAPPED GRID POINT.
    ! IQ........GRID POINT LOCATION IN NORTH-SOUTH OF TRAPPED GRID POINT.
    !
    i1 = 0
    i2 = 0
    j1 = 0
    j2 = 0
    q1 = d_zero
    q2 = d_zero
    p1 = d_zero
    p2 = d_zero
    do i = 1 , iy
      do j = 1 , jx
        i1 = 1000
        do ii = 1 , nlon - 1
          if ( alon(j,i) >= hlon(ii) .and. alon(j,i) < hlon(ii+1) ) then
            p1 = alon(j,i) - hlon(ii)
            p2 = hlon(ii+1) - alon(j,i)
            i1 = ii
            i2 = ii + 1
            exit
          else if ( alon(j,i) >= hlon(ii)-deg360 .and. &
                    alon(j,i) < hlon(ii+1)-deg360 ) then
            p1 = alon(j,i) - (hlon(ii)-deg360)
            p2 = (hlon(ii+1)-deg360) - alon(j,i)
            i1 = ii
            i2 = ii + 1
            exit
          else if ( alon(j,i) >= hlon(ii)+deg360 .and. &
                    alon(j,i) < hlon(ii+1)+deg360 ) then
            p1 = alon(j,i) - (hlon(ii)+deg360)
            p2 = (hlon(ii+1)+deg360) - alon(j,i)
            i1 = ii
            i2 = ii + 1
            exit
          end if
        end do
        if ( alon(j,i) >= hlon(nlon) .and. alon(j,i) < hlon(1)+deg360 ) then
          p1 = alon(j,i) - hlon(nlon)
          p2 = (hlon(1)+deg360) - alon(j,i)
          i1 = nlon
          i2 = 1
        else if ( alon(j,i) >= hlon(nlon)+deg360 .and. &
                  alon(j,i) < hlon(1)+deg720 ) then
          p1 = alon(j,i) - (hlon(nlon)+deg360)
          p2 = (hlon(1)+deg720) - alon(j,i)
          i1 = nlon
          i2 = 1
        else if ( alon(j,i) >= hlon(nlon)-deg360 .and. &
                  alon(j,i) < hlon(1) ) then
          p1 = alon(j,i) - (hlon(nlon)-deg360)
          p2 = hlon(1) - alon(j,i)
          i1 = nlon
          i2 = 1
        else if ( alon(j,i) >= hlon(nlon)-deg720 .and. &
                  alon(j,i) < hlon(1)-deg360 ) then
          p1 = alon(j,i) - (hlon(nlon)-deg720)
          p2 = (hlon(1)-deg360) - alon(j,i)
          i1 = nlon
          i2 = 1
        end if
        if ( i1 == 1000 ) then
          call die('bilinx2','Could not find the right longitude',1)
        end if
        j1 = 1000
        do jj = 1 , nlat - 1
          if ( alat(j,i) >= hlat(jj) .and. alat(j,i) < hlat(jj+1) ) then
            q1 = alat(j,i) - hlat(jj)
            q2 = hlat(jj+1) - alat(j,i)
            j1 = jj
            j2 = jj + 1
            exit
          else if ( alat(j,i) <= hlat(1) ) then
            q1 = d_one
            q2 = d_one
            j1 = 1
            j2 = 1
            exit
          else if ( alat(j,i) >= hlat(nlat) ) then
            q1 = d_one
            q2 = d_one
            j1 = nlat
            j2 = nlat
          end if
        end do
        if ( j1 == 1000 ) then
          call die('bilinx2','Could not find the right latitude',1)
        end if
        if ( j1 > 0 .and. j1 < nlat ) then
          b3(j,i) = ((b2(i1,j1)*p2+b2(i2,j1)*p1)*q2 + &
                     (b2(i1,j2)*p2+b2(i2,j2)*p1)*q1)/(p1+p2)/(q1+q2)
        else if ( j1 == 0 ) then
          ave = d_zero
          do ii = 1 , nlon
            ave = ave + b2(ii,1)
          end do
          ave = ave/real(nlon,rkx)
          b3(j,i) = ((ave*(p1+p2))*q2+(b2(i1,j2)*p2 + &
                      b2(i2,j2)*p1)*q1)/(p1+p2)/(q1+q2)
        else if ( j1 == nlat ) then
          ave = d_zero
          do ii = 1 , nlon
            ave = ave + b2(ii,nlat)
          end do
          ave = ave/real(nlon,rkx)
          b3(j,i) = ((b2(i1,j1)*p2+b2(i2,j1)*p1)*q2 + &
                     (ave*(p1+p2))*q1)/(p1+p2)/(q1+q2)
        end if
      end do
    end do
  end subroutine bilinx2_2d

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
    real(rkx) , dimension(2) :: dists
    integer(ik4) :: i , j , m , mdl , mdr , mul , mur , n , ndl ,  &
               ndr , nul , nur , mx , nx , mm , nn
    logical , dimension(4) :: q
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
        ! Find dists from 4 points around this
        !        2 | 1
        !      ---------
        !        3 | 4
        q(:) = .true.
        nn = nx
        mm = mx + 1
        if ( lonwrap ) then
          if ( mm > nlon ) mm = mm - nlon
        else
          if ( mm > nlon ) then
            write (stderr,*) 'EXCEEDING NLON'
            write (stderr,*) 'NEST DOMAIN TOO NEAR TO PARENT.'
            write (stderr,*) i , j , mm , nn
            write (stderr,*) alon(j,i)
            write (stderr,*) alat(j,i)
            call die('compwgt')
          end if
        end if
        dists(1) = gcdist(glat(mm,nn),glon(mm,nn),alat(j,i),alon(j,i))
        mm = mx - 1
        if ( lonwrap ) then
          if ( mm < 1 ) mm = nlon - mm
        else
          if ( mm < 1 ) then
            write (stderr,*) 'NLON < 1'
            write (stderr,*) 'NEST DOMAIN TOO NEAR TO PARENT.'
            write (stderr,*) i , j , mm , nn
            write (stderr,*) alon(j,i)
            write (stderr,*) alat(j,i)
            call die('compwgt')
          end if
        end if
        dists(2) = gcdist(glat(mm,nn),glon(mm,nn),alat(j,i),alon(j,i))
        if ( dists(1) > dists(2) ) then
          q(1) = .false.
          q(4) = .false.
        else
          q(2) = .false.
          q(3) = .false.
        end if
        mm = mx
        nn = nx + 1
        if ( latpole ) then
          if ( nn > nlat ) nn = nlat
        else
          if ( nn > nlat ) then
            write (stderr,*) 'EXCEEDING NLAT'
            write (stderr,*) 'NEST DOMAIN TOO NEAR TO PARENT.'
            write (stderr,*) i , j , mm , nn
            write (stderr,*) alon(j,i)
            write (stderr,*) alat(j,i)
            call die('compwgt')
          end if
        end if
        dists(1) = gcdist(glat(mm,nn),glon(mm,nn),alat(j,i),alon(j,i))
        nn = nx - 1
        if ( latpole ) then
          if ( nn < 1 ) nn = 1
        else
          if ( nn < 1 ) then
            write (stderr,*) 'NLAT < 1'
            write (stderr,*) 'NEST DOMAIN TOO NEAR TO PARENT.'
            write (stderr,*) i , j , mm , nn
            write (stderr,*) alon(j,i)
            write (stderr,*) alat(j,i)
            call die('compwgt')
          end if
        end if
        dists(2) = gcdist(glat(mm,nn),glon(mm,nn),alat(j,i),alon(j,i))
        if ( dists(1) > dists(2) ) then
          q(1) = .false.
          q(2) = .false.
        else
          q(3) = .false.
          q(4) = .false.
        end if
        if ( q(1) ) then
          mur = mx+1
          nur = nx+1
          mul = mx
          nul = nx+1
          mdr = mx+1
          ndr = nx
          mdl = mx
          ndl = nx
        else if ( q(2) ) then
          mur = mx
          nur = nx+1
          mul = mx-1
          nul = nx+1
          mdr = mx
          ndr = nx
          mdl = mx-1
          ndl = nx
        else if ( q(3) ) then
          mur = mx
          nur = nx
          mul = mx-1
          nul = nx
          mdr = mx
          ndr = nx-1
          mdl = mx-1
          ndl = nx-1
        else if ( q(4) ) then
          mur = mx+1
          nur = nx
          mul = mx
          nul = nx
          mdr = mx+1
          ndr = nx-1
          mdl = mx
          ndl = nx-1
        else
          write (stderr,*) 'LOGIC ERROR in locating point'
          write (stderr,*) i , j , mx , nx
          write (stderr,*) alon(j,i)
          write (stderr,*) alat(j,i)
          call die('compwgt')
        end if
        ! Check if in a global window again
        if ( lonwrap ) then
          if ( mur > nlon ) mur = nlon
          if ( mdr > nlon ) mdr = nlon
          if ( mul < 1 ) mul = 1
          if ( mdl < 1 ) mdl = 1
        end if
        if ( latpole ) then
          if ( nur > nlat ) nur = nlat
          if ( nul > nlat ) nul = nlat
          if ( ndr < 1 ) ndr = 1
          if ( ndr < 1 ) ndr = 1
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
        if ( wa > mindist ) then
          d1xa(j,i) = d_one/(wa**p_factor)
        else
          d1xa(j,i) = d_one
          d1xb(j,i) = d_zero
          d1xc(j,i) = d_zero
          d1xd(j,i) = d_zero
          cycle
        end if
        wb = gcdist(glat(mul,nul),glon(mul,nul),alat(j,i),alon(j,i))
        if ( wb > mindist ) then
          d1xb(j,i) = d_one/(wb**p_factor)
        else
          d1xa(j,i) = d_zero
          d1xb(j,i) = d_one
          d1xc(j,i) = d_zero
          d1xd(j,i) = d_zero
          cycle
        end if
        wc = gcdist(glat(mdr,ndr),glon(mdr,ndr),alat(j,i),alon(j,i))
        if ( wc > mindist ) then
          d1xc(j,i) = d_one/(wc**p_factor)
        else
          d1xa(j,i) = d_zero
          d1xb(j,i) = d_zero
          d1xc(j,i) = d_one
          d1xd(j,i) = d_zero
          cycle
        end if
        wd = gcdist(glat(mdl,ndl),glon(mdl,ndl),alat(j,i),alon(j,i))
        if ( wd > mindist ) then
          d1xd(j,i) = d_one/(wd**p_factor)
        else
          d1xa(j,i) = d_zero
          d1xb(j,i) = d_zero
          d1xc(j,i) = d_zero
          d1xd(j,i) = d_one
          cycle
        end if
      end do
    end do
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
    real(rkx) :: wa , wb , wc , wd , wg
    real(rkx) , dimension(4) :: vv
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
        vv(:) = d_zero
        if ( b2(mur,nur) > missc ) then
          ifound = ifound + 1
          vv(ifound) = b2(mur,nur)*wa
          wg = wg + wa
        end if
        if ( b2(mul,nul) > missc ) then
          ifound = ifound + 1
          vv(ifound) = b2(mul,nul)*wb
          wg = wg + wb
        end if
        if ( b2(mdr,ndr) > missc ) then
          ifound = ifound + 1
          vv(ifound) = b2(mdr,ndr)*wc
          wg = wg + wc
        end if
        if ( b2(mdl,ndl) > missc ) then
          ifound = ifound + 1
          vv(ifound) = b2(mdl,ndl)*wd
          wg = wg + wd
        end if
        if ( ifound /= 0 ) then
          b3(j,i) = sum(vv(1:ifound))/wg
        else
          b3(j,i) = missl
        end if
      end do
    end do
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
      glonmx = maxval(glon)
      glonmn = minval(glon)
      alonmx = maxval(alon)
      alonmn = minval(alon)
      glatmx = maxval(glat)
      glatmn = minval(glat)
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
      glonmx = maxval(glon)
      glonmn = minval(glon)
      alonmx = maxval(alon)
      alonmn = minval(alon)
      glatmx = maxval(glat)
      glatmn = minval(glat)
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
    gcdist_simple = earthrad*acos(crd)*d_r1000
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
    gcdist = earthrad*atan2(y,x)*d_r1000
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
    else if ( xlon(1,xj)   <= xlon(xi,xj)   .and. &
              xlon(1,xj/2) <= xlon(xi,xj/2) .and. &
              xlon(1,1)    <= xlon(xi,1) ) then
      ! it is not crossing timeline
      domain%ntiles = 1
      domain%igstart(1) = int((minlon-glon(1))/dlon) - 1
      domain%igstop(1) = int((maxlon-glon(1))/dlon) + 2
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
      domain%igstart(1) = int((minlon-glon(1))/dlon) - 1
      domain%igstop(1) = gi
      domain%igstart(2) = 1
      domain%igstop(2) = int((maxlon-glon(1))/dlon) + 2
    end if
    if ( has_north_pole(xlat,xi/2) ) then
      ! North pole inside
      l1 = int((minval(xlat(:,1))-glat(1))/dlat) - 1
      l2 = int((minval(xlat(:,xj))-glat(1))/dlat) - 1
      domain%jgstart = min(l1,l2)
      domain%jgstop = gj
    else if ( has_south_pole(xlat,xi/2) ) then
      ! South Pole inside
      l1 = int((maxval(xlat(:,1))-glat(1))/dlat) + 2
      l2 = int((maxval(xlat(:,xj))-glat(1))/dlat) + 2
      domain%jgstart = 1
      domain%jgstop = max(l1,l2)
    else
      domain%jgstart = int((minlat-glat(1))/dlat) - 1
      domain%jgstop  = int((maxlat-glat(1))/dlat) + 2
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

end module mod_interp

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
