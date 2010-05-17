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

      contains
!
!-----------------------------------------------------------------------
!
      subroutine bilinx(fin,fout,lono,lato,loni,lati,nloni,nlati,iy,jx, &
                      & nflds)
      implicit none
!
! Dummy arguments
!
      integer :: iy , jx , nflds , nlati , nloni
      real(4) , dimension(nloni,nlati,nflds) :: fin
      real(4) , dimension(nlati) :: lati
      real(4) , dimension(iy,jx) :: lato , lono
      real(4) , dimension(nloni) :: loni
      real(4) , dimension(iy,jx,nflds) :: fout
      intent (in) fin , iy , jx , lati , lato , loni , lono , nflds ,   &
                & nlati , nloni
      intent (out) fout
!
! Local variables
!
      real(4) :: bas , lon360 , p , q , xsum , xind , yind
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
 
          yind = (((lato(i,j)-lati(1))/(lati(nlati)-lati(1)))           &
               & *float(nlati-1)) + 1.
          jq = int(yind)
          jq = max0(jq,1)
          jqp1 = min0(jq+1,nlati)
          q = yind - jq
 
          lon360 = lono(i,j)
          if ( lono(i,j)<0. ) lon360 = lono(i,j) + 360.
          xind = (((lon360-loni(1))/(loni(nloni)-loni(1)))              &
               & *float(nloni-1)) + 1.
          if ( xind<1.0 .and. lg ) then
            ip = nloni
            ipp1 = 1
            p = xind
          else if ( (xind-nloni)>0.0 .and. lg ) then
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
            if ( fin(ip,jq,l)<-9990.0 .and. fin(ipp1,jq,l)<-9990.0 .and.&
               & fin(ipp1,jqp1,l)<-9990.0 .and. fin(ip,jqp1,l)<-9990.0 )&
               & then
              fout(i,j,l) = -9999.
            else
              if ( fin(ip,jq,l)>-9990.0 ) then
                xsum = xsum + (1.-q)*(1.-p)*fin(ip,jq,l)
                bas = bas + (1.-q)*(1.-p)
              end if
              if ( fin(ipp1,jq,l)>-9990.0 ) then
                xsum = xsum + (1.-q)*p*fin(ipp1,jq,l)
                bas = bas + (1.-q)*p
              end if
              if ( fin(ipp1,jqp1,l)>-9990.0 ) then
                xsum = xsum + q*p*fin(ipp1,jqp1,l)
                bas = bas + q*p
              end if
              if ( fin(ip,jqp1,l)>-9990.0 ) then
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
! Dummy arguments
!
      integer :: iy , jx , llev , nlat , nlon
      real(4) , dimension(jx,iy) :: alat , alon
      real(4) , dimension(nlon,nlat,llev) :: b2
      real(4) , dimension(jx,iy,llev) :: b3
      real(4) , dimension(nlat) :: hlat
      real(4) , dimension(nlon) :: hlon
      intent (in) alat , alon , b2 , hlat , hlon , iy , jx , llev ,     &
                & nlat , nlon
      intent (out) b3
!
! Local variables
!
      real(4) :: ave , p1 , p2 , q1 , q2
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

      do j = 1 , iy
        do i = 1 , jx
 
          i1 = 1000
          do ii = 1 , nlon - 1
            if ( alon(i,j)>=hlon(ii) .and. alon(i,j)<hlon(ii+1) ) then
              p1 = alon(i,j) - hlon(ii)
              p2 = hlon(ii+1) - alon(i,j)
              i1 = ii
              i2 = ii + 1
              go to 20
            else if ( alon(i,j)>=hlon(ii)-360 .and. alon(i,j)<hlon(ii+1)&
                    & -360. ) then
              p1 = alon(i,j) - (hlon(ii)-360.)
              p2 = (hlon(ii+1)-360.) - alon(i,j)
              i1 = ii
              i2 = ii + 1
              go to 20
            else if ( alon(i,j)>=hlon(ii)+360 .and. alon(i,j)<hlon(ii+1)&
                    & +360. ) then
              p1 = alon(i,j) - (hlon(ii)+360.)
              p2 = (hlon(ii+1)+360.) - alon(i,j)
              i1 = ii
              i2 = ii + 1
              go to 20
            else
            end if
          end do
          if ( alon(i,j)>=hlon(nlon) .and. alon(i,j)<hlon(1)+360. ) then
            p1 = alon(i,j) - hlon(nlon)
            p2 = (hlon(1)+360.) - alon(i,j)
            i1 = nlon
            i2 = 1
          else if ( alon(i,j)>=hlon(nlon)+360 .and. alon(i,j)<hlon(1)   &
                  & +720. ) then
            p1 = alon(i,j) - (hlon(nlon)+360.)
            p2 = (hlon(1)+720.) - alon(i,j)
            i1 = nlon
            i2 = 1
          else if ( alon(i,j)>=hlon(nlon)-360 .and. alon(i,j)<hlon(1) ) &
                  & then
            p1 = alon(i,j) - (hlon(nlon)-360.)
            p2 = hlon(1) - alon(i,j)
            i1 = nlon
            i2 = 1
          else if ( alon(i,j)>=hlon(nlon)-720 .and. alon(i,j)<hlon(1)   &
                  & -360. ) then
            p1 = alon(i,j) - (hlon(nlon)-720.)
            p2 = (hlon(1)-360.) - alon(i,j)
            i1 = nlon
            i2 = 1
          else
          end if
 20       continue
          if ( i1==1000 ) stop 'Could not find the right longitute'
          j1 = 1000
          do jj = 1 , nlat - 1
            if ( alat(i,j)>=hlat(jj) .and. alat(i,j)<hlat(jj+1) ) then
              q1 = alat(i,j) - hlat(jj)
              q2 = hlat(jj+1) - alat(i,j)
              j1 = jj
              j2 = jj + 1
            else if ( alat(i,j)<hlat(1) ) then
              if ( alat(i,j)>=-90. ) then
                q1 = alat(i,j) + 90.
                q2 = hlat(1) - alat(i,j)
                j1 = 0
                j2 = 1
              end if
            else if ( alat(i,j)>hlat(nlat) ) then
              if ( alat(i,j)<=90. ) then
                q1 = alat(i,j) - hlat(nlat)
                q2 = 90. - alat(i,j)
                j1 = nlat
                j2 = nlat + 1
              end if
            else
            end if
          end do
          if ( j1==1000 ) stop 'Could not find the right latitude'
          if ( j1>0 .and. j1<nlat ) then
            do l = 1 , llev
              b3(i,j,l) = ((b2(i1,j1,l)*p2+b2(i2,j1,l)*p1)*q2+(b2(i1,j2,&
                        & l)*p2+b2(i2,j2,l)*p1)*q1)/(p1+p2)/(q1+q2)
            end do
          else if ( j1==0 ) then
            do l = 1 , llev
              ave = 0.0
              do ii = 1 , nlon
                ave = ave + b2(ii,1,l)
              end do
              ave = ave/float(nlon)
              b3(i,j,l) = ((ave*(p1+p2))*q2+(b2(i1,j2,l)*p2+b2(i2,j2,l)*&
                        & p1)*q1)/(p1+p2)/(q1+q2)
            end do
          else if ( j1==nlat ) then
            do l = 1 , llev
              ave = 0.0
              do ii = 1 , nlon
                ave = ave + b2(ii,nlat,l)
              end do
              ave = ave/float(nlon)
              b3(i,j,l) = ((b2(i1,j1,l)*p2+b2(i2,j1,l)*p1)*q2+(ave*(p1+ &
                        & p2))*q1)/(p1+p2)/(q1+q2)
            end do
          else
          end if
        end do
      end do

      end subroutine bilinx2
!
!-----------------------------------------------------------------------
!
      subroutine cressmcr(b3,b2,alon,alat,glon,glat,jx,iy,i1ur,i1ul,    &
                        & i1dr,i1dl,j1ur,j1ul,j1dr,j1dl,d1xt,d1xa,d1xb, &
                        & d1xc,d1xd,nlon,nlat,nlev)
      use mod_constants , only : mathpi
      use mod_mxncom
      implicit none
!
! Dummy arguments
!
      integer :: iy , jx , nlat , nlev , nlon
      real(4) , dimension(jx,iy) :: alat , alon , d1xa , d1xb , d1xc ,  &
                               & d1xd , d1xt
      real(4) , dimension(nlon,nlat,nlev,5) :: b2
      real(4) , dimension(jx,iy,nlev,5) :: b3
      real(4) , dimension(nlon,nlat) :: glat , glon
      integer , dimension(jx,iy) :: i1dl , i1dr , i1ul , i1ur , j1dl ,  &
                                  & j1dr , j1ul , j1ur
      intent (in) alat , alon , b2 , glat , glon , iy , jx , nlat ,     &
                & nlev , nlon
      intent (out) b3
      intent (inout) d1xa , d1xb , d1xc , d1xd , d1xt , i1dl , i1dr ,   &
                   & i1ul , i1ur , j1dl , j1dr , j1ul , j1ur
!
! Local variables
!
      real(4) :: aaa , dist , dista , distb , distc , distd
      integer :: i , j , k , l , m , mdl , mdr , mul , mur , n , ndl ,  &
               & ndr , nul , nur
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
      if ( imxmn==0 ) then
        glonmx = -361.
        glonmn = 361.
        do n = 1 , nlat
          do m = 1 , nlon
            if ( glonmx<glon(m,n) ) glonmx = glon(m,n)
            if ( glonmn>glon(m,n) ) glonmn = glon(m,n)
          end do
        end do
        alonmx = -361.
        alonmn = 361.
        do j = 1 , iy
          do i = 1 , jx
            if ( alonmx<alon(i,j) ) alonmx = alon(i,j)
            if ( alonmn>alon(i,j) ) alonmn = alon(i,j)
          end do
        end do
        glatmx = -91.
        glatmn = 91.
        do n = 1 , nlat
          do m = 1 , nlon
            if ( glatmx<glat(m,n) ) glatmx = glat(m,n)
            if ( glatmn>glat(m,n) ) glatmn = glat(m,n)
          end do
        end do
        alatmx = -91.
        alatmn = 91.
        do j = 1 , iy
          do i = 1 , jx
            if ( alatmx<alat(i,j) ) alatmx = alat(i,j)
            if ( alatmn>alat(i,j) ) alatmn = alat(i,j)
          end do
        end do
        write (*,*) 'GLONMN,ALONMN,ALONMX,GLONMX = '
        write (*,*) glonmn , alonmn , alonmx , glonmx
        write (*,*) 'GLATMN,ALATMN,ALATMX,GLATMX = '
        write (*,*) glatmn , alatmn , alatmx , glatmx
        imxmn = 1
      end if
      if ( lcross==0 ) then
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
                if ( (glon(m,n)>=alon(i,j) .and. glon(m,n)-alon(i,j)    &
                   & <10.) .and.                                        &
                   & (glat(m,n)>=alat(i,j) .and. glat(m,n)-alat(i,j)    &
                   & <10.) ) then
                  aaa = ((glon(m,n)-alon(i,j))                          &
                      & *cos((glat(m,n)+alat(i,j))/360.*mathpi))        &
                      & **2 + (glat(m,n)-alat(i,j))**2
                  if ( dista>aaa ) then
                    dista = aaa
                    mur = m
                    nur = n
                  end if
                end if
                if ( (glon(m,n)<alon(i,j) .and. alon(i,j)-glon(m,n)<10.)&
                   & .and.                                              &
                   & (glat(m,n)>=alat(i,j) .and. glat(m,n)-alat(i,j)    &
                   & <10.) ) then
                  aaa = ((glon(m,n)-alon(i,j))                          &
                      & *cos((glat(m,n)+alat(i,j))/360.*mathpi))        &
                      & **2 + (glat(m,n)-alat(i,j))**2
                  if ( distb>aaa ) then
                    distb = aaa
                    mul = m
                    nul = n
                  end if
                end if
                if ( (glon(m,n)>=alon(i,j) .and. glon(m,n)-alon(i,j)    &
                   & <10.) .and.                                        &
                   & (glat(m,n)<alat(i,j) .and. alat(i,j)-glat(m,n)<10.)&
                   & ) then
                  aaa = ((glon(m,n)-alon(i,j))                          &
                      & *cos((glat(m,n)+alat(i,j))/360.*mathpi))        &
                      & **2 + (glat(m,n)-alat(i,j))**2
                  if ( distc>aaa ) then
                    distc = aaa
                    mdr = m
                    ndr = n
                  end if
                end if
                if ( (glon(m,n)<alon(i,j) .and. alon(i,j)-glon(m,n)<10.)&
                   & .and.                                              &
                   & (glat(m,n)<alat(i,j) .and. alat(i,j)-glat(m,n)<10.)&
                   & ) then
                  aaa = ((glon(m,n)-alon(i,j))                          &
                      & *cos((glat(m,n)+alat(i,j))/360.*mathpi))        &
                      & **2 + (glat(m,n)-alat(i,j))**2
                  if ( distd>aaa ) then
                    distd = aaa
                    mdl = m
                    ndl = n
                  end if
                end if
              end do
            end do
            dist = amin1(dista,distb,distc,distd)
 
            i1ur(i,j) = mur
            j1ur(i,j) = nur
            i1ul(i,j) = mul
            j1ul(i,j) = nul
            i1dr(i,j) = mdr
            j1dr(i,j) = ndr
            i1dl(i,j) = mdl
            j1dl(i,j) = ndl
            d1xt(i,j) = dist
            d1xa(i,j) = dista
            d1xb(i,j) = distb
            d1xc(i,j) = distc
            d1xd(i,j) = distd
            do l = 1 , 4
              do k = 1 , nlev
                if ( dist>0.000001 ) then
                  b3(i,j,k,l) = (b2(mur,nur,k,l)/dista+b2(mul,nul,k,l)  &
                              & /distb+b2(mdr,ndr,k,l)                  &
                              & /distc+b2(mdl,ndl,k,l)/distd)           &
                              & /(1./dista+1./distb+1./distc+1./distd)
                else if ( dist==dista ) then
                  b3(i,j,k,l) = b2(mur,nur,k,l)
                else if ( dist==distb ) then
                  b3(i,j,k,l) = b2(mul,nul,k,l)
                else if ( dist==distc ) then
                  b3(i,j,k,l) = b2(mdr,ndr,k,l)
                else if ( dist==distd ) then
                  b3(i,j,k,l) = b2(mdl,ndl,k,l)
                else
                end if
              end do
            end do
          end do
        end do
        lcross = 1
      else
        do j = 1 , iy
          do i = 1 , jx
 
            mur = i1ur(i,j)
            nur = j1ur(i,j)
            mul = i1ul(i,j)
            nul = j1ul(i,j)
            mdr = i1dr(i,j)
            ndr = j1dr(i,j)
            mdl = i1dl(i,j)
            ndl = j1dl(i,j)
            dist = d1xt(i,j)
            dista = d1xa(i,j)
            distb = d1xb(i,j)
            distc = d1xc(i,j)
            distd = d1xd(i,j)
 
            do l = 1 , 4
              do k = 1 , nlev
                if ( dist>0.000001 ) then
                  b3(i,j,k,l) = (b2(mur,nur,k,l)/dista+b2(mul,nul,k,l)  &
                              & /distb+b2(mdr,ndr,k,l)                  &
                              & /distc+b2(mdl,ndl,k,l)/distd)           &
                              & /(1./dista+1./distb+1./distc+1./distd)
                else if ( dist==dista ) then
                  b3(i,j,k,l) = b2(mur,nur,k,l)
                else if ( dist==distb ) then
                  b3(i,j,k,l) = b2(mul,nul,k,l)
                else if ( dist==distc ) then
                  b3(i,j,k,l) = b2(mdr,ndr,k,l)
                else if ( dist==distd ) then
                  b3(i,j,k,l) = b2(mdl,ndl,k,l)
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
      subroutine cressmdt(b3,b2,alon,alat,glon,glat,jx,iy,i1ur,i1ul,    &
                        & i1dr,i1dl,j1ur,j1ul,j1dr,j1dl,d1xt,d1xa,d1xb, &
                        & d1xc,d1xd,nlon,nlat,nlev)
      use mod_constants , only : mathpi
      use mod_mxncom
      implicit none
!
! Dummy arguments
!
      integer :: iy , jx , nlat , nlev , nlon
      real(4) , dimension(jx,iy) :: alat , alon , d1xa , d1xb , d1xc ,  &
                               & d1xd , d1xt
      real(4) , dimension(nlon,nlat,nlev,5) :: b2
      real(4) , dimension(jx,iy,nlev,5) :: b3
      real(4) , dimension(nlon,nlat) :: glat , glon
      integer , dimension(jx,iy) :: i1dl , i1dr , i1ul , i1ur , j1dl ,  &
                                  & j1dr , j1ul , j1ur
      intent (in) alat , alon , b2 , glat , glon , iy , jx , nlat ,     &
                & nlev , nlon
      intent (out) b3
      intent (inout) d1xa , d1xb , d1xc , d1xd , d1xt , i1dl , i1dr ,   &
                   & i1ul , i1ur , j1dl , j1dr , j1ul , j1ur
!
! Local variables
!
      real(4) :: aaa , dist , dista , distb , distc , distd
      integer :: i , j , k , l , m , mdl , mdr , mul , mur , n , ndl ,  &
               & ndr , nul , nur
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
      if ( imxmn==0 ) then
        glonmx = -361.
        glonmn = 361.
        do n = 1 , nlat
          do m = 1 , nlon
            if ( glonmx<glon(m,n) ) glonmx = glon(m,n)
            if ( glonmn>glon(m,n) ) glonmn = glon(m,n)
          end do
        end do
        alonmx = -361.
        alonmn = 361.
        do j = 1 , iy
          do i = 1 , jx
            if ( alonmx<alon(i,j) ) alonmx = alon(i,j)
            if ( alonmn>alon(i,j) ) alonmn = alon(i,j)
          end do
        end do
        glatmx = -91.
        glatmn = 91.
        do n = 1 , nlat
          do m = 1 , nlon
            if ( glatmx<glat(m,n) ) glatmx = glat(m,n)
            if ( glatmn>glat(m,n) ) glatmn = glat(m,n)
          end do
        end do
        alatmx = -91.
        alatmn = 91.
        do j = 1 , iy
          do i = 1 , jx
            if ( alatmx<alat(i,j) ) alatmx = alat(i,j)
            if ( alatmn>alat(i,j) ) alatmn = alat(i,j)
          end do
        end do
        write (*,*) 'GLONMN,ALONMN,ALONMX,GLONMX = '
        write (*,*) glonmn , alonmn , alonmx , glonmx
        write (*,*) 'GLATMN,ALATMN,ALATMX,GLATMX = '
        write (*,*) glatmn , alatmn , alatmx , glatmx
        imxmn = 1
      end if
      if ( ldot==0 ) then
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
                if ( (glon(m,n)>=alon(i,j) .and. glon(m,n)-alon(i,j)    &
                   & <10.) .and.                                        &
                   & (glat(m,n)>=alat(i,j) .and. glat(m,n)-alat(i,j)    &
                   & <10.) ) then
                  aaa = ((glon(m,n)-alon(i,j))                          &
                      & *cos((glat(m,n)+alat(i,j))/360.*mathpi))        &
                      & **2 + (glat(m,n)-alat(i,j))**2
                  if ( dista>aaa ) then
                    dista = aaa
                    mur = m
                    nur = n
                  end if
                end if
                if ( (glon(m,n)<alon(i,j) .and. alon(i,j)-glon(m,n)<10.)&
                   & .and.                                              &
                   & (glat(m,n)>=alat(i,j) .and. glat(m,n)-alat(i,j)    &
                   & <10.) ) then
                  aaa = ((glon(m,n)-alon(i,j))                          &
                      & *cos((glat(m,n)+alat(i,j))/360.*mathpi))        &
                      & **2 + (glat(m,n)-alat(i,j))**2
                  if ( distb>aaa ) then
                    distb = aaa
                    mul = m
                    nul = n
                  end if
                end if
                if ( (glon(m,n)>=alon(i,j) .and. glon(m,n)-alon(i,j)    &
                   & <10.) .and.                                        &
                   & (glat(m,n)<alat(i,j) .and. alat(i,j)-glat(m,n)<10.)&
                   & ) then
                  aaa = ((glon(m,n)-alon(i,j))                          &
                      & *cos((glat(m,n)+alat(i,j))/360.*mathpi))        &
                      & **2 + (glat(m,n)-alat(i,j))**2
                  if ( distc>aaa ) then
                    distc = aaa
                    mdr = m
                    ndr = n
                  end if
                end if
                if ( (glon(m,n)<alon(i,j) .and. alon(i,j)-glon(m,n)<10.)&
                   & .and.                                              &
                   & (glat(m,n)<alat(i,j) .and. alat(i,j)-glat(m,n)<10.)&
                   & ) then
                  aaa = ((glon(m,n)-alon(i,j))                          &
                      & *cos((glat(m,n)+alat(i,j))/360.*mathpi))        &
                      & **2 + (glat(m,n)-alat(i,j))**2
                  if ( distd>aaa ) then
                    distd = aaa
                    mdl = m
                    ndl = n
                  end if
                end if
              end do
            end do
            dist = amin1(dista,distb,distc,distd)
 
            i1ur(i,j) = mur
            j1ur(i,j) = nur
            i1ul(i,j) = mul
            j1ul(i,j) = nul
            i1dr(i,j) = mdr
            j1dr(i,j) = ndr
            i1dl(i,j) = mdl
            j1dl(i,j) = ndl
            d1xt(i,j) = dist
            d1xa(i,j) = dista
            d1xb(i,j) = distb
            d1xc(i,j) = distc
            d1xd(i,j) = distd
            do l = 1 , 3
              do k = 1 , nlev
                if ( dist>0.000001 ) then
                  b3(i,j,k,l) = (b2(mur,nur,k,l)/dista+b2(mul,nul,k,l)  &
                              & /distb+b2(mdr,ndr,k,l)                  &
                              & /distc+b2(mdl,ndl,k,l)/distd)           &
                              & /(1./dista+1./distb+1./distc+1./distd)
                else if ( dist==dista ) then
                  b3(i,j,k,l) = b2(mur,nur,k,l)
                else if ( dist==distb ) then
                  b3(i,j,k,l) = b2(mul,nul,k,l)
                else if ( dist==distc ) then
                  b3(i,j,k,l) = b2(mdr,ndr,k,l)
                else if ( dist==distd ) then
                  b3(i,j,k,l) = b2(mdl,ndl,k,l)
                else
                end if
              end do
            end do
          end do
        end do
        ldot = 1
      else
        do j = 1 , iy
          do i = 1 , jx
 
            mur = i1ur(i,j)
            nur = j1ur(i,j)
            mul = i1ul(i,j)
            nul = j1ul(i,j)
            mdr = i1dr(i,j)
            ndr = j1dr(i,j)
            mdl = i1dl(i,j)
            ndl = j1dl(i,j)
            dist = d1xt(i,j)
            dista = d1xa(i,j)
            distb = d1xb(i,j)
            distc = d1xc(i,j)
            distd = d1xd(i,j)
 
            do l = 1 , 3
              do k = 1 , nlev
                if ( dist>0.000001 ) then
                  b3(i,j,k,l) = (b2(mur,nur,k,l)/dista+b2(mul,nul,k,l)  &
                              & /distb+b2(mdr,ndr,k,l)                  &
                              & /distc+b2(mdl,ndl,k,l)/distd)           &
                              & /(1./dista+1./distb+1./distc+1./distd)
                else if ( dist==dista ) then
                  b3(i,j,k,l) = b2(mur,nur,k,l)
                else if ( dist==distb ) then
                  b3(i,j,k,l) = b2(mul,nul,k,l)
                else if ( dist==distc ) then
                  b3(i,j,k,l) = b2(mdr,ndr,k,l)
                else if ( dist==distd ) then
                  b3(i,j,k,l) = b2(mdl,ndl,k,l)
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
