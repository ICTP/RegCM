      subroutine bilinx(fin,loni,lati,nloni,nlati,fout,lono,lato,iy,jx, &
                      & nflds)
      implicit none
!
! Dummy arguments
!
      integer :: iy , jx , nflds , nlati , nloni
      real , dimension(nloni,nlati,nflds) :: fin
      real , dimension(nlati) :: lati
      real , dimension(iy,jx) :: lato , lono
      real , dimension(nloni) :: loni
      real , dimension(iy,jx,nflds) :: fout
      intent (in) fin , iy , jx , lati , lato , loni , lono , nflds ,   &
                & nlati , nloni
      intent (out) fout
!
! Local variables
!
      real :: bas , lon360 , p , q , xsum , xind , yind
      integer :: i , ip , ipp1 , j , jq , jqp1 , l
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
 
!     POINT.
 
      do j = 1 , jx
        do i = 1 , iy
 
          yind = (((lato(i,j)-lati(1))/(lati(nlati)-lati(1)))           &
               & *float(nlati-1)) + 1.
          jq = int(yind)
          jqp1 = min0(jq+1,nlati)
          q = yind - jq
 
          lon360 = lono(i,j)
          if ( lono(i,j)<0. ) lon360 = lono(i,j) + 360.
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
        end do
 
      end do
 
      end subroutine bilinx
