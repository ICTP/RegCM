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

      subroutine bilinxo(fin,loni,lati,nloni,nlati,fout,lono,lato,iy,jx,&
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
 
!     POINT.
!     write(*,*)'BILINX',IN(201,140,1),IN(200,141,1)
!     write(*,*)'BILINX',LONI(20),LATI(35)
 
      do j = 1 , jx
        do i = 1 , iy
 
          yind = (((lato(i,j)-lati(1))/(lati(nlati)-lati(1)))           &
               & *float(nlati-1)) + 1.
          jq = int(yind)
          jqp1 = min0(jq+1,nlati)
          q = yind - jq
 
          lon360 = lono(i,j)
 
!         IF(LONO(I,J).LT.0. .and. INVNTRY .eq. 'edgar') THEN
!         IF(LONO(I,J).LT. 0. )THEN
!     write(*,*)'NOHA',LONO(I,J)
          if ( loni(1)>0 .and. lono(i,j)<0. ) lon360 = lono(i,j) + 360.
          xind = (((lon360-loni(1))/(loni(nloni)-loni(1)))              &
               & *float(nloni-1)) + 1.
          ip = int(xind)
          ipp1 = min0(ip+1,nloni)
          p = xind - ip
 
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
!         write(*,*)'test',IN(IP,JQ,1),OUT(I,J,1)
        end do
      end do
 
      end subroutine bilinxo
