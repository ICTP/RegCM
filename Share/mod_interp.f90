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

  private

  public :: bilinx , bilinx2 , cressmcr , cressmdt , distwgtcr , distwgtdt
  public :: gcdist

  real(rk8) :: alatmn , alatmx , alonmn , alonmx
  real(rk8) :: glatmn , glatmx , glonmn , glonmx

  integer(ik4) :: imxmn = 0
  integer(ik4) :: lcross = 0
  integer(ik4) :: ldot = 0

  real(rk8) , parameter :: deg720 = d_two*deg360
  real(rk8) , parameter :: missl = -9999.0D0
  real(rk8) , parameter :: missc = -9990.0D0

  real(rk8) , pointer , dimension(:,:) :: dc1xa , dc1xb , dc1xc , dc1xd
  real(rk8) , pointer , dimension(:,:) :: dd1xa , dd1xb , dd1xc , dd1xd
  integer(ik4) , pointer, dimension(:,:) :: ic1dl , ic1dr , ic1ul , ic1ur , &
                                       jc1dl , jc1dr , jc1ul , jc1ur
  integer(ik4) , pointer, dimension(:,:) :: id1dl , id1dr , id1ul , id1ur , &
                                       jd1dl , jd1dr , jd1ul , jd1ur
  logical :: lonwrap = .false. , latpole = .false.

  contains
!
!-----------------------------------------------------------------------
!
  subroutine bilinx(fin,fout,lono,lato,loni,lati,nloni,nlati,jx,iy,nflds)
  implicit none
!
  integer(ik4) :: jx , iy , nflds , nlati , nloni
  real(rk8) , dimension(nloni,nlati,nflds) :: fin
  real(rk8) , dimension(nlati) :: lati
  real(rk8) , dimension(jx,iy) :: lato , lono
  real(rk8) , dimension(nloni) :: loni
  real(rk8) , dimension(jx,iy,nflds) :: fout
  intent (in) fin , jx , iy , lati , lato , loni , lono , nflds ,   &
              nlati , nloni
  intent (out) fout
!
  real(rk8) :: bas , lon360 , p , q , xsum , xind , yind
  integer(ik4) :: i , ip , ipp1 , j , jq , jqp1 , l
  logical :: lg
!
!
!     PERFORMING BI-LINEAR INTERPOLATION USING 4 GRID POINTS FROM A
!     BIGGER RECTANGULAR GRID TO A GRID DESCRIBED BY XLONS AND XLATS OF
!     GRID2. A POINT ON GRID2 IS TRAPPED WITHIN FOUR GRID POINTS ON
!     GRID4.THE GRID POINTS ARE ALWAYS TO THE NORTH AND EAST OF THE
!     TRAPPED POINT.. THE ALGORITHM COMPUTES THE FRACTIONAL DISTANCES
!     IN BOTH X AND Y DIRECTION OF THE TRAPPED GRID POINT AND USES THE
!     INFORMATION AS WEIGHTING FACTORS IN THE INTERPOLATION.
!     THERE IS ONE LESS ROW AND COLUMN WHEN THE SCALAR FIELDS ARE
!     INTERPOLATED BECAUSE XLATS AND XLONS ARE NOT DEFINED FOR
!     THE CROSS POINTS IN THE MM4 MODEL.
!
!     IN(NLONI,NLATI,NFLDS)  IS THE INPUT FIELD ON REGULAR LAT/LON GRID.
!     OUT(NLATO,NLONO,NFLDS) IS THE OUTPUT FIELD ON LAMBERT CONFORMAL
!     GRID. LONI.....LONGITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
!     LATI.....LATITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
!     P.........EAST-WEST WEIGHTING FACTOR.
!     Q.........NORTH-SOUTH WEIGHTING FACTOR.
!     IP........GRID POINT LOCATION IN EAST-WEST OF TRAPPED GRID POINT.
!     IQ........GRID POINT LOCATION IN NORTH-SOUTH OF TRAPPED GRID
!               POINT.
!
!     Global dataset ?
! 
  lg = .true.
!
!
  do i = 1 , iy
    do j = 1 , jx
 
      yind = (((lato(j,i)-lati(1))/(lati(nlati)-lati(1)))*dble(nlati-1))+d_one
      jq = int(yind)
      jq = max0(jq,1)
      jqp1 = min0(jq+1,nlati)
      q = yind - jq
 
      lon360 = lono(j,i)
      if ( lono(j,i) < deg00 ) lon360 = lono(j,i) + deg360
      xind = (((lon360-loni(1))/(loni(nloni)-loni(1)))*dble(nloni-1))+d_one
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
        ip = max0(ip,1)
        ipp1 = min0(ip+1,nloni)
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
          fout(j,i,l) = xsum/bas
        end if
      end do
    end do
 
  end do
 
  end subroutine bilinx
!
!-----------------------------------------------------------------------
!
  subroutine bilinx2(b3,b2,alon,alat,hlon,hlat,nlon,nlat,jx,iy,llev)
  implicit none
!
  integer(ik4) :: iy , jx , llev , nlat , nlon
  real(rk8) , dimension(jx,iy) :: alat , alon
  real(rk8) , dimension(nlon,nlat,llev) :: b2
  real(rk8) , dimension(jx,iy,llev) :: b3
  real(rk8) , dimension(nlat) :: hlat
  real(rk8) , dimension(nlon) :: hlon
  intent (in) alat , alon , b2 , hlat , hlon , iy , jx , llev , nlat , nlon
  intent (out) b3
!
  real(rk8) :: ave , p1 , p2 , q1 , q2
  integer(ik4) :: i , i1 , i2 , ii , j , j1 , j2 , jj , l
!
!     PERFORMING BI-LINEAR INTERPOLATION USING 4 GRID POINTS FROM A
!     BIGGER RECTANGULAR GRID TO A GRID DESCRIBED BY ALON AND ALAT OF
!     GRID2. A POINT ON GRID2 IS TRAPPED WITHIN FOUR GRID POINTS ON
!     GRID4.THE GRID POINTS ARE ALWAYS TO THE NORTH AND EAST OF THE
!     TRAPPED POINT. THE ALGORITHM COMPUTES THE FRACTIONAL DISTANCES IN
!     BOTH X AND Y DIRECTION OF THE TRAPPED GRID POINT AND USES THE
!     INFORMATION AS WEIGHTING FACTORS IN THE INTERPOLATION.
!     THERE IS ONE LESS ROW AND COLUMN WHEN THE SCALAR FIELDS ARE
!     INTERPOLATED BECAUSE ALAT AND ALON ARE NOT DEFINED FOR
!     THE CROSS POINTS IN THE RegCM MODEL.
!
!     B2(JX,IX,NLEV) IS THE INPUT FIELD ON REGULAR LAT/LON GRID.
!     B3(JX,IX,NLEV) IS THE OUTPUT FIELD ON LAMBERT CONFORMAL GRID.
!     HLON......LONGITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
!     HLAT......LATITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
!     P.........EAST-WEST WEIGHTING FACTOR.
!     Q.........NORTH-SOUTH WEIGHTING FACTOR.
!     IP........GRID POINT LOCATION IN EAST-WEST OF TRAPPED GRID POINT.
!     IQ........GRID POINT LOCATION IN NORTH-SOUTH OF TRAPPED GRID
!     POINT.
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
      else if ( alon(j,i) >= hlon(nlon)-deg360 .and. alon(j,i) < hlon(1) ) then
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
          ave = ave/dble(nlon)
          b3(j,i,l) = ((ave*(p1+p2))*q2+(b2(i1,j2,l)*p2+b2(i2,j2,l)* &
                                     p1)*q1)/(p1+p2)/(q1+q2)
        end do
      else if ( j1 == nlat ) then
        do l = 1 , llev
          ave = d_zero
          do ii = 1 , nlon
            ave = ave + b2(ii,nlat,l)
          end do
          ave = ave/dble(nlon)
          b3(j,i,l) = ((b2(i1,j1,l)*p2+b2(i2,j1,l)*p1)*q2+ &
                       (ave*(p1+p2))*q1)/(p1+p2)/(q1+q2)
        end do
      else
      end if
    end do
  end do

  end subroutine bilinx2
!
  subroutine distwgtcr(b3,b2,alon,alat,glon,glat,jx,iy,nlon,nlat)
  implicit none
!
  integer(ik4) :: iy , jx , nlat , nlon
  real(rk8) , dimension(jx,iy) :: alat , alon
  real(rk8) , dimension(jx,iy) :: b3
  real(rk8) , dimension(nlon,nlat) :: glat , glon
  real(rk8) , dimension(nlon,nlat) :: b2
  intent (in) alat , alon , b2 , glat , glon , iy , jx , nlat , nlon
  intent (out) b3
!
  real(rk8) :: dist , dista , distb , distc , distd , distx
  real(rk8) , dimension(2) :: dists
  real(rk8) :: v1 , v2 , v3 , v4 , wg
  integer(ik4) :: i , j , m , mdl , mdr , mul , mur , n , ndl ,  &
             ndr , nul , nur , ifound , mx , nx , mm , nn
  logical , dimension(4) :: q
!
!     FIND THE FOUR CLOSEST POINTS TO THE GRID WE WANT TO HAVE VALUE,
!     THEN DO THE AVERAGE OF THOSE FOUR POINTS WEIGHTED BY THE DISTANCE.
!     THE CLOSEST ONE HAS BIGGEST WEIGHT.
!
!     B2(JX,IX,NLEV) IS THE INPUT FIELD ON PREVIOUS regular or
!     irregular GRID. B3(JX,IX,NLEV) IS THE OUTPUT FIELD ON new
!     (regular or irregular) GRID. GLON......LONGITUDE VALUES IN
!     DEGREES OF THE INTERMEDIATE GRID4. GLAT......LATITUDE VALUES IN
!     DEGREES OF THE INTERMEDIATE GRID4.
!
!
!     Find the Minimum and Maximum of GLON, GLAT, ALON and ALAT
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
    if ( glonmx - glonmn > 350.0D0 ) lonwrap   = .true.
    if ( glatmx - glatmn > 170.0D0 ) latpole = .true.
    imxmn = 1
  end if

  if ( lcross == 0 ) then
    if (.not. associated(ic1dl)) call getmem2d(ic1dl,1,jx,1,iy,'interp:ic1dl')
    if (.not. associated(ic1dr)) call getmem2d(ic1dr,1,jx,1,iy,'interp:ic1dr')
    if (.not. associated(ic1ul)) call getmem2d(ic1ul,1,jx,1,iy,'interp:ic1ul')
    if (.not. associated(ic1ur)) call getmem2d(ic1ur,1,jx,1,iy,'interp:ic1ur')
    if (.not. associated(jc1dl)) call getmem2d(jc1dl,1,jx,1,iy,'interp:jc1dl')
    if (.not. associated(jc1dr)) call getmem2d(jc1dr,1,jx,1,iy,'interp:jc1dr')
    if (.not. associated(jc1ul)) call getmem2d(jc1ul,1,jx,1,iy,'interp:jc1ul')
    if (.not. associated(jc1ur)) call getmem2d(jc1ur,1,jx,1,iy,'interp:jc1ur')
    if (.not. associated(dc1xa)) call getmem2d(dc1xa,1,jx,1,iy,'interp:dc1xa')
    if (.not. associated(dc1xb)) call getmem2d(dc1xb,1,jx,1,iy,'interp:dc1xb')
    if (.not. associated(dc1xc)) call getmem2d(dc1xc,1,jx,1,iy,'interp:dc1xc')
    if (.not. associated(dc1xd)) call getmem2d(dc1xd,1,jx,1,iy,'interp:dc1xd')
    write (stdout,*) 'FIRST TIME in CRESSMCR'
    write (stdout,*) 'Calculating weights.... (will take long time)'
    do i = 1 , iy
      do j = 1 , jx
        ! Find nearest point
        distx = 1.0D+08
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
            call die('cressmcr')
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
            call die('cressmcr')
          end if
        end if
        dists(2) = gcdist(glat(mm,nn),glon(mm,nn),alat(j,i),alon(j,i))
        if ( dists(1) > dists(2) ) then
          q(2) = .false.
          q(4) = .false.
        else
          q(1) = .false.
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
            call die('cressmcr')
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
            call die('cressmcr')
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
          mur = mx
          nur = nx+1
          mul = mx-1
          nul = nx-1
          mdr = mx
          ndr = nx
          mdl = mx-1
          ndl = nx
        else if ( q(2) ) then
          mur = mx+1
          nur = nx+1
          mul = mx
          nul = nx+1
          mdr = mx+1
          ndr = nx
          mdl = mx
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
          call die('cressmcr')
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
        dista = gcdist(glat(mur,nur),glon(mur,nur),alat(j,i),alon(j,i))
        distb = gcdist(glat(mul,nul),glon(mul,nul),alat(j,i),alon(j,i))
        distc = gcdist(glat(mdr,ndr),glon(mdr,ndr),alat(j,i),alon(j,i))
        distd = gcdist(glat(mdl,ndl),glon(mdl,ndl),alat(j,i),alon(j,i))
        ic1ur(j,i) = mur
        jc1ur(j,i) = nur
        ic1ul(j,i) = mul
        jc1ul(j,i) = nul
        ic1dr(j,i) = mdr
        jc1dr(j,i) = ndr
        ic1dl(j,i) = mdl
        jc1dl(j,i) = ndl
        dc1xa(j,i) = (d_one/dista)**2
        dc1xb(j,i) = (d_one/distb)**2
        dc1xc(j,i) = (d_one/distc)**2
        dc1xd(j,i) = (d_one/distd)**2
      end do
    end do
    write (stdout,*) 'Done.'
    lcross = 1
  end if

  b3(:,:) = missl

  do i = 1 , iy
    do j = 1 , jx
      mur = ic1ur(j,i)
      nur = jc1ur(j,i)
      mul = ic1ul(j,i)
      nul = jc1ul(j,i)
      mdr = ic1dr(j,i)
      ndr = jc1dr(j,i)
      mdl = ic1dl(j,i)
      ndl = jc1dl(j,i)
      dista = dc1xa(j,i)
      distb = dc1xb(j,i)
      distc = dc1xc(j,i)
      distd = dc1xd(j,i)
 
      ifound = 0
      wg = d_zero
      v1 = d_zero
      v2 = d_zero
      v3 = d_zero
      v4 = d_zero
      if ( b2(mur,nur) > missc ) then
        v1 = b2(mur,nur)*dista
        ifound = 1
        wg = wg + dista
      end if
      if ( b2(mul,nul) > missc ) then
        v2 = b2(mul,nul)*distb
        ifound = 1
        wg = wg + distb
      end if
      if ( b2(mdr,ndr) > missc ) then
        v3 = b2(mdr,ndr)*distc
        ifound = 1
        wg = wg + distc
      end if
      if ( b2(mdl,ndl) > missc ) then
        v4= b2(mdl,ndl)*distd
        ifound = 1
        wg = wg + distd
      end if
      if ( ifound /= 0 ) b3(j,i) = v1/wg+v2/wg+v3/wg+v4/wg
    end do
  end do
  end subroutine distwgtcr
!
  subroutine distwgtdt(b3,b2,alon,alat,glon,glat,jx,iy,nlon,nlat)
  implicit none
!
  integer(ik4) :: iy , jx , nlat , nlon
  real(rk8) , dimension(jx,iy) :: alat , alon
  real(rk8) , dimension(jx,iy) :: b3
  real(rk8) , dimension(nlon,nlat) :: glat , glon
  real(rk8) , dimension(nlon,nlat) :: b2
  intent (in) alat , alon , b2 , glat , glon , iy , jx , nlat , nlon
  intent (out) b3
!
  real(rk8) :: dist , dista , distb , distc , distd , distx
  real(rk8) , dimension(2) :: dists
  real(rk8) :: v1 , v2 , v3 , v4 , wg
  integer(ik4) :: i , j , m , mdl , mdr , mul , mur , n , ndl ,  &
             ndr , nul , nur , ifound , mx , nx , mm , nn
  logical , dimension(4) :: q
!
!     FIND THE FOUR CLOSEST POINTS TO THE GRID WE WANT TO HAVE VALUE,
!     THEN DO THE AVERAGE OF THOSE FOUR POINTS WEIGHTED BY THE DISTANCE.
!     THE CLOSEST ONE HAS BIGGEST WEIGHT.
!
!     B2(JX,IX,NLEV) IS THE INPUT FIELD ON PREVIOUS regular or
!     irregular GRID. B3(JX,IX,NLEV) IS THE OUTPUT FIELD ON new
!     (regular or irregular) GRID. GLON......LONGITUDE VALUES IN
!     DEGREES OF THE INTERMEDIATE GRID4. GLAT......LATITUDE VALUES IN
!     DEGREES OF THE INTERMEDIATE GRID4.
!
!
!     Find the Minimum and Maximum of GLON, GLAT, ALON and ALAT
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
    if ( glonmx - glonmn > 350.0D0 ) lonwrap   = .true.
    if ( glatmx - glatmn > 170.0D0 ) latpole = .true.
    imxmn = 1
  end if

  if ( ldot == 0 ) then
    if (.not. associated(id1dl)) call getmem2d(id1dl,1,jx,1,iy,'interp:id1dl')
    if (.not. associated(id1dr)) call getmem2d(id1dr,1,jx,1,iy,'interp:id1dr')
    if (.not. associated(id1ul)) call getmem2d(id1ul,1,jx,1,iy,'interp:id1ul')
    if (.not. associated(id1ur)) call getmem2d(id1ur,1,jx,1,iy,'interp:id1ur')
    if (.not. associated(jd1dl)) call getmem2d(jd1dl,1,jx,1,iy,'interp:jd1dl')
    if (.not. associated(jd1dr)) call getmem2d(jd1dr,1,jx,1,iy,'interp:jd1dr')
    if (.not. associated(jd1ul)) call getmem2d(jd1ul,1,jx,1,iy,'interp:jd1ul')
    if (.not. associated(jd1ur)) call getmem2d(jd1ur,1,jx,1,iy,'interp:jd1ur')
    if (.not. associated(dd1xa)) call getmem2d(dd1xa,1,jx,1,iy,'interp:dd1xa')
    if (.not. associated(dd1xb)) call getmem2d(dd1xb,1,jx,1,iy,'interp:dd1xb')
    if (.not. associated(dd1xc)) call getmem2d(dd1xc,1,jx,1,iy,'interp:dd1xc')
    if (.not. associated(dd1xd)) call getmem2d(dd1xd,1,jx,1,iy,'interp:dd1xd')
    write (stdout,*) 'FIRST TIME in CRESSMDT'
    write (stdout,*) 'Calculating weights.... (will take long time)'
    do i = 1 , iy
      do j = 1 , jx
        ! Find nearest point
        distx = 1.0D+08
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
            call die('cressmdt')
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
            call die('cressmdt')
          end if
        end if
        dists(2) = gcdist(glat(mm,nn),glon(mm,nn),alat(j,i),alon(j,i))
        if ( dists(1) > dists(2) ) then
          q(2) = .false.
          q(4) = .false.
        else
          q(1) = .false.
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
            call die('cressmdt')
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
            call die('cressmdt')
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
          mur = mx
          nur = nx+1
          mul = mx-1
          nul = nx-1
          mdr = mx
          ndr = nx
          mdl = mx-1
          ndl = nx
        else if ( q(2) ) then
          mur = mx+1
          nur = nx+1
          mul = mx
          nul = nx+1
          mdr = mx+1
          ndr = nx
          mdl = mx
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
          call die('cressmdt')
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
        dista = gcdist(glat(mur,nur),glon(mur,nur),alat(j,i),alon(j,i))
        distb = gcdist(glat(mul,nul),glon(mul,nul),alat(j,i),alon(j,i))
        distc = gcdist(glat(mdr,ndr),glon(mdr,ndr),alat(j,i),alon(j,i))
        distd = gcdist(glat(mdl,ndl),glon(mdl,ndl),alat(j,i),alon(j,i))
        id1ur(j,i) = mur
        jd1ur(j,i) = nur
        id1ul(j,i) = mul
        jd1ul(j,i) = nul
        id1dr(j,i) = mdr
        jd1dr(j,i) = ndr
        id1dl(j,i) = mdl
        jd1dl(j,i) = ndl
        dd1xa(j,i) = (d_one/dista)**2
        dd1xb(j,i) = (d_one/distb)**2
        dd1xc(j,i) = (d_one/distc)**2
        dd1xd(j,i) = (d_one/distd)**2
      end do
    end do
    write (stdout,*) 'Done.'
    ldot = 1
  end if

  b3(:,:) = missl

  do i = 1 , iy
    do j = 1 , jx
      mur = id1ur(j,i)
      nur = jd1ur(j,i)
      mul = id1ul(j,i)
      nul = jd1ul(j,i)
      mdr = id1dr(j,i)
      ndr = jd1dr(j,i)
      mdl = id1dl(j,i)
      ndl = jd1dl(j,i)
      dista = dd1xa(j,i)
      distb = dd1xb(j,i)
      distc = dd1xc(j,i)
      distd = dd1xd(j,i)
 
      ifound = 0
      wg = d_zero
      v1 = d_zero
      v2 = d_zero
      v3 = d_zero
      v4 = d_zero
      if ( b2(mur,nur) > missc ) then
        v1 = b2(mur,nur)*dista
        ifound = 1
        wg = wg + dista
      end if
      if ( b2(mul,nul) > missc ) then
        v2 = b2(mul,nul)*distb
        ifound = 1
        wg = wg + distb
      end if
      if ( b2(mdr,ndr) > missc ) then
        v3 = b2(mdr,ndr)*distc
        ifound = 1
        wg = wg + distc
      end if
      if ( b2(mdl,ndl) > missc ) then
        v4= b2(mdl,ndl)*distd
        ifound = 1
        wg = wg + distd
      end if
      if ( ifound /= 0 ) b3(j,i) = v1/wg+v2/wg+v3/wg+v4/wg
    end do
  end do
  end subroutine distwgtdt
!
!-----------------------------------------------------------------------
!
  subroutine cressmcr(b3,b2,alon,alat,glon,glat,jx,iy,nlon,nlat,nlev,nf)
  implicit none
!
  integer(ik4) :: iy , jx , nlat , nlev , nlon , nf
  real(rk8) , dimension(jx,iy) :: alat , alon
  real(rk8) , dimension(jx,iy,nlev*nf) :: b3
  real(rk8) , dimension(nlon,nlat) :: glat , glon
  real(rk8) , dimension(nlon,nlat,nlev*nf) :: b2
  intent (in) alat , alon , b2 , glat , glon , iy , jx , nlat ,     &
              nlev , nlon , nf
  intent (out) b3
!
  integer(ik4) :: k , l , kin
!
  do l = 1 , nf
    do k = 1 , nlev
       kin = (l-1)*nlev+k
       call distwgtcr(b3(:,:,kin),b2(:,:,kin),alon,alat, &
                         glon,glat,jx,iy,nlon,nlat)
    end do
  end do
 
  end subroutine cressmcr
!
!-----------------------------------------------------------------------
!
  subroutine cressmdt(b3,b2,alon,alat,glon,glat,jx,iy,nlon,nlat,nlev,nf)
  implicit none
!
  integer(ik4) :: iy , jx , nlat , nlev , nlon , nf
  real(rk8) , dimension(jx,iy) :: alat , alon
  real(rk8) , dimension(jx,iy,nlev*nf) :: b3
  real(rk8) , dimension(nlon,nlat) :: glat , glon
  real(rk8) , dimension(nlon,nlat,nlev*nf) :: b2
  intent (in) alat , alon , b2 , glat , glon , iy , jx , nlat ,     &
              nlev , nlon , nf
  intent (out) b3
!
  integer(ik4) :: k , l , kin
!
  do l = 1 , nf
    do k = 1 , nlev
       kin = (l-1)*nlev+k
       call distwgtdt(b3(:,:,kin),b2(:,:,kin),alon,alat, &
                         glon,glat,jx,iy,nlon,nlat)
    end do
  end do
 
  end subroutine cressmdt
!
  real(rk8) function gcdist(lat1,lon1,lat2,lon2)
    implicit none
    real(rk8) , intent(in) :: lat1 , lon1 , lat2, lon2
    real(rk8) :: clat1 , slat1 , clat2 , slat2 , cdlon , crd

    clat1 = cos(lat1*degrad)
    slat1 = sin(lat1*degrad)
    clat2 = cos(lat2*degrad)
    slat2 = sin(lat2*degrad)
    cdlon = cos((lon1-lon2)*degrad)
    crd   = slat1*slat2+clat1*clat2*cdlon
    gcdist = earthrad*acos(crd)
  end function gcdist
!
end module mod_interp
