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
 
  use m_realkinds
  use m_die
  use m_stdio

  use mod_constants , only : degrad

  real(sp) :: alatmn , alatmx , alonmn , alonmx , glatmn , glatmx , &
              glonmn , glonmx
  real(sp) :: dlatmax , dlonmax
  integer :: imxmn , lcross , ldot

  real(sp) , allocatable , dimension(:,:) :: dc1xa , dc1xb , dc1xc ,&
                                            dc1xd , dc1xt
  integer , allocatable, dimension(:,:) :: ic1dl , ic1dr , ic1ul ,  &
                              ic1ur , jc1dl , jc1dr , jc1ul , jc1ur
  real(sp) , allocatable , dimension(:,:) :: dd1xa , dd1xb , dd1xc ,&
                                             dd1xd , dd1xt
  integer , allocatable, dimension(:,:) :: id1dl , id1dr , id1ul ,  &
                              id1ur , jd1dl , jd1dr , jd1ul , jd1ur

  data imxmn /0/
  data lcross /0/
  data ldot /0/
  contains
!
!-----------------------------------------------------------------------
!
  subroutine bilinx(fin,fout,lono,lato,loni,lati,nloni,nlati,iy,jx,nflds)
  implicit none
!
  integer :: iy , jx , nflds , nlati , nloni
  real(sp) , dimension(nloni,nlati,nflds) :: fin
  real(sp) , dimension(nlati) :: lati
  real(sp) , dimension(iy,jx) :: lato , lono
  real(sp) , dimension(nloni) :: loni
  real(sp) , dimension(iy,jx,nflds) :: fout
  intent (in) fin , iy , jx , lati , lato , loni , lono , nflds ,   &
              nlati , nloni
  intent (out) fout
!
  real(sp) :: bas , lon360 , p , q , xsum , xind , yind
  integer :: i , ip , ipp1 , j , jq , jqp1 , l
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
  do j = 1 , jx
    do i = 1 , iy
 
      yind = (((lato(i,j)-lati(1))/(lati(nlati)-lati(1)))*float(nlati-1))+1.
      jq = int(yind)
      jq = max0(jq,1)
      jqp1 = min0(jq+1,nlati)
      q = yind - jq
 
      lon360 = lono(i,j)
      if ( lono(i,j) < 0. ) lon360 = lono(i,j) + 360.
      xind = (((lon360-loni(1))/(loni(nloni)-loni(1)))*float(nloni-1))+1.
      if ( xind < 1.0 .and. lg ) then
        ip = nloni
        ipp1 = 1
        p = xind
      else if ( (xind-nloni) > 0.0 .and. lg ) then
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
        xsum = 0.0
        bas = 0.0
        if ( fin(ip,jq,l) < -9990.0 .and. fin(ipp1,jq,l) < -9990.0 .and. &
             fin(ipp1,jqp1,l) < -9990.0 .and. fin(ip,jqp1,l) < -9990.0 ) then
          fout(i,j,l) = -9999.
        else
          if ( fin(ip,jq,l) > -9990.0 ) then
            xsum = xsum + (1.-q)*(1.-p)*fin(ip,jq,l)
            bas = bas + (1.-q)*(1.-p)
          end if
          if ( fin(ipp1,jq,l) > -9990.0 ) then
            xsum = xsum + (1.-q)*p*fin(ipp1,jq,l)
            bas = bas + (1.-q)*p
          end if
          if ( fin(ipp1,jqp1,l) > -9990.0 ) then
            xsum = xsum + q*p*fin(ipp1,jqp1,l)
            bas = bas + q*p
          end if
          if ( fin(ip,jqp1,l) > -9990.0 ) then
            xsum = xsum + q*(1.-p)*fin(ip,jqp1,l)
            bas = bas + q*(1.-p)
          end if
          fout(i,j,l) = xsum/bas
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
  integer :: iy , jx , llev , nlat , nlon
  real(sp) , dimension(jx,iy) :: alat , alon
  real(sp) , dimension(nlon,nlat,llev) :: b2
  real(sp) , dimension(jx,iy,llev) :: b3
  real(sp) , dimension(nlat) :: hlat
  real(sp) , dimension(nlon) :: hlon
  intent (in) alat , alon , b2 , hlat , hlon , iy , jx , llev , nlat , nlon
  intent (out) b3
!
  real(sp) :: ave , p1 , p2 , q1 , q2
  integer :: i , i1 , i2 , ii , j , j1 , j2 , jj , l
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
  q1 = 0.0
  q2 = 0.0
  p1 = 0.0
  p2 = 0.0

  do j = 1 , iy
    do i = 1 , jx
 
      i1 = 1000
      do ii = 1 , nlon - 1
        if ( alon(i,j) >= hlon(ii) .and. alon(i,j) < hlon(ii+1) ) then
          p1 = alon(i,j) - hlon(ii)
          p2 = hlon(ii+1) - alon(i,j)
          i1 = ii
          i2 = ii + 1
          exit
        else if ( alon(i,j) >= hlon(ii)-360 .and. &
                  alon(i,j) < hlon(ii+1)-360. ) then
          p1 = alon(i,j) - (hlon(ii)-360.)
          p2 = (hlon(ii+1)-360.) - alon(i,j)
          i1 = ii
          i2 = ii + 1
          exit
        else if ( alon(i,j) >= hlon(ii)+360 .and. &
                  alon(i,j) < hlon(ii+1)+360. ) then
          p1 = alon(i,j) - (hlon(ii)+360.)
          p2 = (hlon(ii+1)+360.) - alon(i,j)
          i1 = ii
          i2 = ii + 1
          exit
        end if
      end do
      if ( alon(i,j) >= hlon(nlon) .and. alon(i,j) < hlon(1)+360. ) then
        p1 = alon(i,j) - hlon(nlon)
        p2 = (hlon(1)+360.) - alon(i,j)
        i1 = nlon
        i2 = 1
      else if ( alon(i,j) >= hlon(nlon)+360 .and. &
                alon(i,j) < hlon(1)+720. ) then
        p1 = alon(i,j) - (hlon(nlon)+360.)
        p2 = (hlon(1)+720.) - alon(i,j)
        i1 = nlon
        i2 = 1
      else if ( alon(i,j) >= hlon(nlon)-360 .and. alon(i,j) < hlon(1) ) then
        p1 = alon(i,j) - (hlon(nlon)-360.)
        p2 = hlon(1) - alon(i,j)
        i1 = nlon
        i2 = 1
      else if ( alon(i,j) >= hlon(nlon)-720 .and. &
                alon(i,j) < hlon(1)-360. ) then
        p1 = alon(i,j) - (hlon(nlon)-720.)
        p2 = (hlon(1)-360.) - alon(i,j)
        i1 = nlon
        i2 = 1
      end if
      if ( i1 == 1000 ) then
        call die('cressmdt','Could not find the right longitude',1)
      end if
      j1 = 1000
      do jj = 1 , nlat - 1
        if ( alat(i,j) >= hlat(jj) .and. alat(i,j) < hlat(jj+1) ) then
          q1 = alat(i,j) - hlat(jj)
          q2 = hlat(jj+1) - alat(i,j)
          j1 = jj
          j2 = jj + 1
          exit
        else if ( alat(i,j) <= hlat(1) ) then
          q1 = 1.0
          q2 = 1.0
          j1 = 1
          j2 = 1
          exit
        else if ( alat(i,j) >= hlat(nlat) ) then
          q1 = 1.0
          q2 = 1.0
          j1 = nlat
          j2 = nlat
        end if
      end do
      if ( j1 == 1000 ) then
        call die('cressmdt','Could not find the right latitude',1)
      end if
      if ( j1 > 0 .and. j1 < nlat ) then
        do l = 1 , llev
          b3(i,j,l) = ((b2(i1,j1,l)*p2+b2(i2,j1,l)*p1)*q2+(b2(i1,j2,l)* &
                        p2+b2(i2,j2,l)*p1)*q1)/(p1+p2)/(q1+q2)
        end do
      else if ( j1 == 0 ) then
        do l = 1 , llev
          ave = 0.0
          do ii = 1 , nlon
            ave = ave + b2(ii,1,l)
          end do
          ave = ave/float(nlon)
          b3(i,j,l) = ((ave*(p1+p2))*q2+(b2(i1,j2,l)*p2+b2(i2,j2,l)* &
                                     p1)*q1)/(p1+p2)/(q1+q2)
        end do
      else if ( j1 == nlat ) then
        do l = 1 , llev
          ave = 0.0
          do ii = 1 , nlon
            ave = ave + b2(ii,nlat,l)
          end do
          ave = ave/float(nlon)
          b3(i,j,l) = ((b2(i1,j1,l)*p2+b2(i2,j1,l)*p1)*q2+ &
                       (ave*(p1+p2))*q1)/(p1+p2)/(q1+q2)
        end do
      else
      end if
    end do
  end do

  end subroutine bilinx2
!
!-----------------------------------------------------------------------
!
  subroutine cressmcr(b3,b2,alon,alat,glon,glat,jx,iy,nlon,nlat,nlev,nf)
  implicit none
!
  integer :: iy , jx , nlat , nlev , nlon , nf
  real(sp) , dimension(jx,iy) :: alat , alon
  real(sp) , dimension(jx,iy,nlev*nf) :: b3
  real(sp) , dimension(nlon,nlat) :: glat , glon
  real(sp) , dimension(nlon,nlat,nlev*nf) :: b2
  intent (in) alat , alon , b2 , glat , glon , iy , jx , nlat ,     &
              nlev , nlon , nf
  intent (out) b3
!
  real(sp) :: aaa , dist , dista , distb , distc , distd
  integer :: i , j , k , l , m , mdl , mdr , mul , mur , n , ndl ,  &
             ndr , nul , nur , kin
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
    dlatmax = (glatmx-glatmn)/nlat * 2.
    dlonmax = (glonmx-glonmn)/nlon * 2.
    write (stdout,*) 'GLONMN,ALONMN,ALONMX,GLONMX = '
    write (stdout,*) glonmn , alonmn , alonmx , glonmx
    write (stdout,*) 'GLATMN,ALATMN,ALATMX,GLATMX = '
    write (stdout,*) glatmn , alatmn , alatmx , glatmx
    imxmn = 1
  end if

  if ( lcross == 0 ) then
    if (.not. allocated(ic1dl)) allocate (ic1dl(jx,iy))
    if (.not. allocated(ic1dr)) allocate (ic1dr(jx,iy))
    if (.not. allocated(ic1ul)) allocate (ic1ul(jx,iy))
    if (.not. allocated(ic1ur)) allocate (ic1ur(jx,iy))
    if (.not. allocated(jc1dl)) allocate (jc1dl(jx,iy))
    if (.not. allocated(jc1dr)) allocate (jc1dr(jx,iy))
    if (.not. allocated(jc1ul)) allocate (jc1ul(jx,iy))
    if (.not. allocated(jc1ur)) allocate (jc1ur(jx,iy))
    if (.not. allocated(dc1xa)) allocate (dc1xa(jx,iy))
    if (.not. allocated(dc1xb)) allocate (dc1xb(jx,iy))
    if (.not. allocated(dc1xc)) allocate (dc1xc(jx,iy))
    if (.not. allocated(dc1xd)) allocate (dc1xd(jx,iy))
    if (.not. allocated(dc1xt)) allocate (dc1xt(jx,iy))
    write (stdout,*) 'FIRST TIME in CRESSMCR'
    write (stdout,*) 'Calculating weights.... (will take long time)'
    do j = 1 , iy
      do i = 1 , jx
 
        mur = -1000
        nur = -1000
        mul = -1000
        nul = -1000
        mdr = -1000
        ndr = -1000
        mdl = -1000
        ndl = -1000
 
        dista = 1.E8
        distb = 1.E8
        distc = 1.E8
        distd = 1.E8
        do n = 2 , nlat
          do m = 2 , nlon
            if ( (glon(m,n) >= alon(i,j) .and.                      &
                  glon(m,n)-alon(i,j) < dlonmax ) .and.             &
                 (glat(m,n) >= alat(i,j) .and.                      &
                  glat(m,n)-alat(i,j) < dlatmax ) ) then
              aaa = real((dble(glon(m,n)-alon(i,j)) *               &
                    dcos(dble(glat(m,n)+alat(i,j))*degrad))**2.0D0 +&
                    dble(glat(m,n)-alat(i,j))**2.0D0)
              if ( dista > aaa ) then
                dista = aaa
                mur = m
                nur = n
              end if
            end if
            if ( (glon(m,n) < alon(i,j) .and.                       &
                  alon(i,j)-glon(m,n) < dlonmax ) .and.             &
                 (glat(m,n) >= alat(i,j) .and.                      &
                  glat(m,n)-alat(i,j) < dlatmax ) ) then
              aaa = real((dble(glon(m,n)-alon(i,j)) *               &
                    dcos(dble(glat(m,n)+alat(i,j))*degrad))**2.0D0 +&
                    dble(glat(m,n)-alat(i,j))**2.0D0)
              if ( distb > aaa ) then
                distb = aaa
                mul = m
                nul = n
              end if
            end if
            if ( (glon(m,n) >= alon(i,j) .and.                      &
                  glon(m,n)-alon(i,j) < dlonmax ) .and.             &
                 (glat(m,n) < alat(i,j) .and.                       &
                  alat(i,j)-glat(m,n) < dlatmax ) ) then
              aaa = real((dble(glon(m,n)-alon(i,j)) *               &
                    dcos(dble(glat(m,n)+alat(i,j))*degrad))**2.0D0 +&
                    dble(glat(m,n)-alat(i,j))**2.0D0)
              if ( distc > aaa ) then
                distc = aaa
                mdr = m
                ndr = n
              end if
            end if
            if ( (glon(m,n) < alon(i,j) .and.                       &
                  alon(i,j)-glon(m,n) < dlonmax ) .and.             &
                 (glat(m,n) < alat(i,j) .and.                       &
                  alat(i,j)-glat(m,n) < dlatmax ) ) then
              aaa = real((dble(glon(m,n)-alon(i,j)) *               &
                    dcos(dble(glat(m,n)+alat(i,j))*degrad))**2.0D0 +&
                    dble(glat(m,n)-alat(i,j))**2.0D0)
              if ( distd > aaa ) then
                distd = aaa
                mdl = m
                ndl = n
              end if
            end if
          end do
        end do

        if ( mur < 0. .or. nur < 0. .or. mul < 0. .or. nul < 0. .or.&
             mdr < 0. .or. ndr < 0. .or. mdl < 0. .or. ndl < 0 )    &
        then
          write (stderr,*) 'NEST DOMAIN TOO NEAR TO PARENT.'
          write (stderr,*) mur , nur , mdr , ndr
          write (stderr,*) mul , nul , mdl , ndl
          write (stderr,*) i , j
          write (stderr,*) alon(i,j)
          write (stderr,*) alat(i,j)
          call die('cressmcr')
        end if

        dist = amin1(dista,distb,distc,distd)
 
        ic1ur(i,j) = mur
        jc1ur(i,j) = nur
        ic1ul(i,j) = mul
        jc1ul(i,j) = nul
        ic1dr(i,j) = mdr
        jc1dr(i,j) = ndr
        ic1dl(i,j) = mdl
        jc1dl(i,j) = ndl
        dc1xt(i,j) = dist
        dc1xa(i,j) = dista
        dc1xb(i,j) = distb
        dc1xc(i,j) = distc
        dc1xd(i,j) = distd

        do l = 1 , nf
          do k = 1 , nlev
            kin = (l-1)*nlev+k
            if ( dist > 0.000001 ) then
              b3(i,j,kin) = (b2(mur,nur,kin)/dista+b2(mul,nul,kin)  &
                            /distb+b2(mdr,ndr,kin)                  &
                            /distc+b2(mdl,ndl,kin)/distd)           &
                            /(1./dista+1./distb+1./distc+1./distd)
            else if ( dist == dista ) then
              b3(i,j,kin) = b2(mur,nur,kin)
            else if ( dist == distb ) then
              b3(i,j,kin) = b2(mul,nul,kin)
            else if ( dist == distc ) then
              b3(i,j,kin) = b2(mdr,ndr,kin)
            else if ( dist == distd ) then
              b3(i,j,kin) = b2(mdl,ndl,kin)
            else
            end if
          end do
        end do
      end do
    end do
    write (stdout,*) 'Done.'
    lcross = 1
  else
    do j = 1 , iy
      do i = 1 , jx
 
        mur = ic1ur(i,j)
        nur = jc1ur(i,j)
        mul = ic1ul(i,j)
        nul = jc1ul(i,j)
        mdr = ic1dr(i,j)
        ndr = jc1dr(i,j)
        mdl = ic1dl(i,j)
        ndl = jc1dl(i,j)
        dist = dc1xt(i,j)
        dista = dc1xa(i,j)
        distb = dc1xb(i,j)
        distc = dc1xc(i,j)
        distd = dc1xd(i,j)
 
        do l = 1 , nf
          do k = 1 , nlev
            kin = (l-1)*nlev+k
            if ( dist > 0.000001 ) then
              b3(i,j,kin) = (b2(mur,nur,kin)/dista+b2(mul,nul,kin)  &
                            /distb+b2(mdr,ndr,kin)                  &
                            /distc+b2(mdl,ndl,kin)/distd)           &
                            /(1./dista+1./distb+1./distc+1./distd)
            else if ( dist == dista ) then
              b3(i,j,kin) = b2(mur,nur,kin)
            else if ( dist == distb ) then
              b3(i,j,kin) = b2(mul,nul,kin)
            else if ( dist == distc ) then
              b3(i,j,kin) = b2(mdr,ndr,kin)
            else if ( dist == distd ) then
              b3(i,j,kin) = b2(mdl,ndl,kin)
            else
            end if
          end do
        end do
      end do
    end do
  end if
 
  end subroutine cressmcr
!
!-----------------------------------------------------------------------
!
  subroutine cressmdt(b3,b2,alon,alat,glon,glat,jx,iy,nlon,nlat,nlev,nf)
  implicit none
!
  integer :: iy , jx , nlat , nlev , nlon , nf
  real(sp) , dimension(jx,iy) :: alat , alon
  real(sp) , dimension(jx,iy,nlev*nf) :: b3
  real(sp) , dimension(nlon,nlat) :: glat , glon
  real(sp) , dimension(nlon,nlat,nlev*nf) :: b2
  intent (in) alat , alon , b2 , glat , glon , iy , jx , nlat ,     &
              nlev , nlon , nf
  intent (out) b3
!
  real(sp) :: aaa , dist , dista , distb , distc , distd
  integer :: i , j , k , l , m , mdl , mdr , mul , mur , n , ndl ,  &
             ndr , nul , nur , kin
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
    dlatmax = (glatmx-glatmn)/nlat * 2.
    dlonmax = (glonmx-glonmn)/nlon * 2.
    write (stdout,*) 'GLONMN,ALONMN,ALONMX,GLONMX = '
    write (stdout,*) glonmn , alonmn , alonmx , glonmx
    write (stdout,*) 'GLATMN,ALATMN,ALATMX,GLATMX = '
    write (stdout,*) glatmn , alatmn , alatmx , glatmx
    imxmn = 1
  end if

  if ( ldot == 0 ) then
    if (.not. allocated(id1dl)) allocate (id1dl(jx,iy))
    if (.not. allocated(id1dr)) allocate (id1dr(jx,iy))
    if (.not. allocated(id1ul)) allocate (id1ul(jx,iy))
    if (.not. allocated(id1ur)) allocate (id1ur(jx,iy))
    if (.not. allocated(jd1dl)) allocate (jd1dl(jx,iy))
    if (.not. allocated(jd1dr)) allocate (jd1dr(jx,iy))
    if (.not. allocated(jd1ul)) allocate (jd1ul(jx,iy))
    if (.not. allocated(jd1ur)) allocate (jd1ur(jx,iy))
    if (.not. allocated(dd1xa)) allocate (dd1xa(jx,iy))
    if (.not. allocated(dd1xb)) allocate (dd1xb(jx,iy))
    if (.not. allocated(dd1xc)) allocate (dd1xc(jx,iy))
    if (.not. allocated(dd1xd)) allocate (dd1xd(jx,iy))
    if (.not. allocated(dd1xt)) allocate (dd1xt(jx,iy))
    write (stdout,*) 'FIRST TIME in CRESSMDT'
    write (stdout,*) 'Calculating weights.... (will take long time)'
    do j = 1 , iy
      do i = 1 , jx
 
        mur = 1000
        nur = 1000
        mul = 1000
        nul = 1000
        mdr = 1000
        ndr = 1000
        mdl = 1000
        ndl = 1000
 
        dista = 1.E8
        distb = 1.E8
        distc = 1.E8
        distd = 1.E8
        do n = 2 , nlat
          do m = 2 , nlon
            if ( (glon(m,n) >= alon(i,j) .and.                      &
                  glon(m,n)-alon(i,j) < dlonmax ) .and.             &
                 (glat(m,n) >= alat(i,j) .and.                      &
                  glat(m,n)-alat(i,j) < dlatmax ) ) then
              aaa = real((dble(glon(m,n)-alon(i,j)) *               &
                    dcos(dble(glat(m,n)+alat(i,j))*degrad))**2.0D0 +&
                    dble(glat(m,n)-alat(i,j))**2.0D0)
              if ( dista > aaa ) then
                dista = aaa
                mur = m
                nur = n
              end if
            end if
            if ( (glon(m,n) < alon(i,j) .and.                       &
                  alon(i,j)-glon(m,n) < dlonmax ) .and.             &
                 (glat(m,n) >= alat(i,j) .and.                      &
                  glat(m,n)-alat(i,j) < dlatmax ) ) then
              aaa = real((dble(glon(m,n)-alon(i,j)) *               &
                    dcos(dble(glat(m,n)+alat(i,j))*degrad))**2.0D0 +&
                    dble(glat(m,n)-alat(i,j))**2.0D0)
              if ( distb > aaa ) then
                distb = aaa
                mul = m
                nul = n
              end if
            end if
            if ( (glon(m,n) >= alon(i,j) .and.                      &
                  glon(m,n)-alon(i,j) < dlonmax ) .and.             &
                 (glat(m,n) < alat(i,j) .and.                       &
                  alat(i,j)-glat(m,n) < dlatmax ) ) then
              aaa = real((dble(glon(m,n)-alon(i,j)) *               &
                    dcos(dble(glat(m,n)+alat(i,j))*degrad))**2.0D0 +&
                    dble(glat(m,n)-alat(i,j))**2.0D0)
              if ( distc > aaa ) then
                distc = aaa
                mdr = m
                ndr = n
              end if
            end if
            if ( (glon(m,n) < alon(i,j) .and.                       &
                  alon(i,j)-glon(m,n) < dlonmax ) .and.             &
                 (glat(m,n) < alat(i,j) .and.                       &
                  alat(i,j)-glat(m,n) < dlatmax ) ) then
              aaa = real((dble(glon(m,n)-alon(i,j)) *               &
                    dcos(dble(glat(m,n)+alat(i,j))*degrad))**2.0D0 +&
                    dble(glat(m,n)-alat(i,j))**2.0D0)
              if ( distd > aaa ) then
                distd = aaa
                mdl = m
                ndl = n
              end if
            end if
          end do
        end do

        if ( mur < 0. .or. nur < 0. .or. mul < 0. .or. nul < 0. .or.&
             mdr < 0. .or. ndr < 0. .or. mdl < 0. .or. ndl < 0 )    &
        then
          write (stderr,*) 'NEST DOMAIN TOO NEAR TO PARENT.'
          write (stderr,*) mur , nur , mdr , ndr
          write (stderr,*) mul , nul , mdl , ndl
          write (stderr,*) i , j
          write (stderr,*) alon(i,j)
          write (stderr,*) alat(i,j)
          call die('cressmdt')
        end if

        dist = amin1(dista,distb,distc,distd)
 
        id1ur(i,j) = mur
        jd1ur(i,j) = nur
        id1ul(i,j) = mul
        jd1ul(i,j) = nul
        id1dr(i,j) = mdr
        jd1dr(i,j) = ndr
        id1dl(i,j) = mdl
        jd1dl(i,j) = ndl
        dd1xt(i,j) = dist
        dd1xa(i,j) = dista
        dd1xb(i,j) = distb
        dd1xc(i,j) = distc
        dd1xd(i,j) = distd
        do l = 1 , nf
          do k = 1 , nlev
            kin = (l-1)*nlev+k
            if ( dist > 0.000001 ) then
              b3(i,j,kin) = (b2(mur,nur,kin)/dista+b2(mul,nul,kin)  &
                            /distb+b2(mdr,ndr,kin)                  &
                            /distc+b2(mdl,ndl,kin)/distd)           &
                            /(1./dista+1./distb+1./distc+1./distd)
            else if ( dist == dista ) then
              b3(i,j,kin) = b2(mur,nur,kin)
            else if ( dist == distb ) then
              b3(i,j,kin) = b2(mul,nul,kin)
            else if ( dist == distc ) then
              b3(i,j,kin) = b2(mdr,ndr,kin)
            else if ( dist == distd ) then
              b3(i,j,kin) = b2(mdl,ndl,kin)
            else
            end if
          end do
        end do
      end do
    end do
    ldot = 1
    write (stdout,*) 'Done.'
  else
    do j = 1 , iy
      do i = 1 , jx
 
        mur = id1ur(i,j)
        nur = jd1ur(i,j)
        mul = id1ul(i,j)
        nul = jd1ul(i,j)
        mdr = id1dr(i,j)
        ndr = jd1dr(i,j)
        mdl = id1dl(i,j)
        ndl = jd1dl(i,j)
        dist = dd1xt(i,j)
        dista = dd1xa(i,j)
        distb = dd1xb(i,j)
        distc = dd1xc(i,j)
        distd = dd1xd(i,j)
 
        do l = 1 , nf
          do k = 1 , nlev
            kin = (l-1)*nlev+k
            if ( dist > 0.000001 ) then
              b3(i,j,kin) = (b2(mur,nur,kin)/dista+b2(mul,nul,kin)  &
                            /distb+b2(mdr,ndr,kin)                  &
                            /distc+b2(mdl,ndl,kin)/distd)           &
                            /(1./dista+1./distb+1./distc+1./distd)
            else if ( dist == dista ) then
              b3(i,j,kin) = b2(mur,nur,kin)
            else if ( dist == distb ) then
              b3(i,j,kin) = b2(mul,nul,kin)
            else if ( dist == distc ) then
              b3(i,j,kin) = b2(mdr,ndr,kin)
            else if ( dist == distd ) then
              b3(i,j,kin) = b2(mdl,ndl,kin)
            else
            end if
          end do
        end do
      end do
    end do
  end if
 
  end subroutine cressmdt
!
end module mod_interp
