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

module mod_ecwcp

  use mod_dynparam
  use mod_memutil
  use mod_stdio
  use mod_realkinds
  use mod_message

  private

  integer , parameter :: nlev = 15 , jlat = 64 , ilon = 128

  real(sp) , dimension(jlat) :: hlat
  real(sp) , dimension(ilon) :: hlon
  real(sp) , dimension(nlev) :: sigma1 , sigmar

  real(sp) , dimension(ilon,jlat,nlev) :: w1

  real(sp) , target , dimension(ilon,jlat,nlev*3) :: b2
  real(sp) , target , dimension(ilon,jlat,nlev*2) :: d2
  real(sp) , pointer , dimension(:,:,:) :: b3
  real(sp) , pointer , dimension(:,:,:) :: d3

  real(sp) , pointer , dimension(:,:,:) :: t1 , q1 , h1
  real(sp) , pointer , dimension(:,:,:) :: u1 , v1
  real(sp) , pointer , dimension(:,:,:) :: t3 , q3 , h3
  real(sp) , pointer , dimension(:,:,:) :: u3 , v3

  public :: getecwcp , headerec

  contains

  subroutine getecwcp(idate)
  use mod_grid
  use mod_write
  use mod_interp , only : bilinx2
  use mod_vertint
  use mod_hgt
  use mod_humid
  use mod_mksst
  use mod_uvrot
  use mod_vectutil
  implicit none
!
  type(rcm_time_and_date) , intent(in) :: idate
!
  character(12) , dimension(12,5) :: finm
  integer :: i , j , k , nrec
  integer :: year , month , day , hour
  logical :: there
!
  data finm/'ECT421993JAN' , 'ECT421993FEB' , 'ECT421993MAR' ,  &
            'ECT421993APR' , 'ECT421993MAY' , 'ECT421993JUN' ,  &
            'ECT421993JUL' , 'ECT421993AUG' , 'ECT421993SEP' ,  &
            'ECT421993OCT' , 'ECT421993NOV' , 'ECT421993DEC' ,  &
            'ECT421994JAN' , 'ECT421994FEB' , 'ECT421994MAR' ,  &
            'ECT421994APR' , 'ECT421994MAY' , 'ECT421994JUN' ,  &
            'ECT421994JUL' , 'ECT421994AUG' , 'ECT421994SEP' ,  &
            'ECT421994OCT' , 'ECT421994NOV' , 'ECT421994DEC' ,  &
            'ECT421995JAN' , 'ECT421995FEB' , 'ECT421995MAR' ,  &
            'ECT421995APR' , 'ECT421995MAY' , 'ECT421995JUN' ,  &
            'ECT421995JUL' , 'ECT421995AUG' , 'ECT421995SEP' ,  &
            'ECT421995OCT' , 'ECT421995NOV' , 'ECT421995DEC' ,  &
            'ECT421996JAN' , 'ECT421996FEB' , 'ECT421996MAR' ,  &
            'ECT421996APR' , 'ECT421996MAY' , 'ECT421996JUN' ,  &
            'ECT421996JUL' , 'ECT421996AUG' , 'ECT421996SEP' ,  &
            'ECT421996OCT' , 'ECT421996NOV' , 'ECT421996DEC' ,  &
            'ECT421997JAN' , 'ECT421997FEB' , 'ECT421997MAR' ,  &
            'ECT421997APR' , 'ECT421997MAY' , 'ECT421997JUN' ,  &
            'ECT421997JUL' , 'ECT421997AUG' , 'ECT421997SEP' ,  &
            'ECT421997OCT' , 'ECT421997NOV' , 'ECT421997DEC'/
!
 
  call split_idate(idate,year,month,day,hour)
!
  inquire (file=trim(inpglob)//'/ECWCRP/'// &
           finm(month,year-1992), exist=there)
  if ( .not.there ) then
    call die('getecwcp',trim(inpglob)//'/ECWCRP/'// &
             finm(month,year-1992)//' is not available',1)
  end if
  open (63,file=trim(inpglob)//'/ECWCRP/'// &
        finm(month,year-1992),  &
        form='unformatted',recl=ilon*jlat*ibyte,access='direct')
  nrec = (((day-1)*4+hour)/6)*(nlev*6+1)
  nrec = nrec + 1
  do k = 1 , nlev
    nrec = nrec + 1
    read (63,rec=nrec) ((h1(i,j,k),i=1,ilon),j=1,jlat)
  end do
  do k = 1 , nlev
    nrec = nrec + 1
    read (63,rec=nrec) ((t1(i,j,k),i=1,ilon),j=1,jlat)
  end do
  do k = 1 , nlev
    nrec = nrec + 1
    read (63,rec=nrec) ((u1(i,j,k),i=1,ilon),j=1,jlat)
  end do
  do k = 1 , nlev
    nrec = nrec + 1
    read (63,rec=nrec) ((v1(i,j,k),i=1,ilon),j=1,jlat)
  end do
  do k = 1 , nlev
    nrec = nrec + 1
    read (63,rec=nrec) ((w1(i,j,k),i=1,ilon),j=1,jlat)
  end do
  do k = 1 , nlev
    nrec = nrec + 1
    read (63,rec=nrec) ((q1(i,j,k),i=1,ilon),j=1,jlat)
  end do
  write (stdout,*) 'READ IN fields at DATE:' , tochar(idate) , ' from ' , &
              finm(month,year-1992)
  close (21)
!
!     HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
  call bilinx2(b3,b2,xlon,xlat,hlon,hlat,ilon,jlat,jx,iy,nlev*3)
  call bilinx2(d3,d2,dlon,dlat,hlon,hlat,ilon,jlat,jx,iy,nlev*2)
!
!     ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
  call uvrot4(u3,v3,dlon,dlat,clon,clat,xcone,jx,iy,nlev,plon,plat,iproj)
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!     V E R T I C A L   I N T E R P O L A T I O N
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!     ******           NEW CALCULATION OF P* ON RegCM TOPOGRAPHY.
  call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,nlev)
 
  call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
  if(i_band == 1) then
     call p1p2_band(b3pd,ps4,jx,iy)
  else
     call p1p2(b3pd,ps4,jx,iy)
  endif
 
!
!     F0    DETERMINE SURFACE TEMPS ON RegCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
  call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,nlev)
 
  call readsst(ts4,idate)

!     F3     INTERPOLATE U, V, T, AND Q.
  call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,nlev)
  call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,nlev)
!
  call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,iy,kz,nlev)
 
  call humid1(t3,q3,100.0,0.0D0,sigma1,jx,iy,nlev)
  call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,iy,kz,nlev)
  call humid2(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
!
!     F4     DETERMINE H
  call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,dsigma,jx,iy,kz)
!
  end subroutine getecwcp
!
  subroutine headerec
  implicit none
!
  integer :: i , k , kr
!
  hlat(1) = -87.8638
  hlat(2) = -85.0965
  hlat(3) = -82.3129
  hlat(4) = -79.5256
  hlat(5) = -76.7369
  hlat(6) = -73.9475
  hlat(7) = -71.1578
  hlat(8) = -68.3678
  hlat(9) = -65.5776
  hlat(10) = -62.7874
  hlat(11) = -59.9970
  hlat(12) = -57.2066
  hlat(13) = -54.4162
  hlat(14) = -51.6257
  hlat(15) = -48.8352
  hlat(16) = -46.0447
  hlat(17) = -43.2542
  hlat(18) = -40.4636
  hlat(19) = -37.6731
  hlat(20) = -34.8825
  hlat(21) = -32.0919
  hlat(22) = -29.3014
  hlat(23) = -26.5108
  hlat(24) = -23.7202
  hlat(25) = -20.9296
  hlat(26) = -18.1390
  hlat(27) = -15.3484
  hlat(28) = -12.5578
  hlat(29) = -9.76715
  hlat(30) = -6.97653
  hlat(31) = -4.18592
  hlat(32) = -1.39531
  hlat(33) = 1.39531
  hlat(34) = 4.18592
  hlat(35) = 6.97653
  hlat(36) = 9.76715
  hlat(37) = 12.5578
  hlat(38) = 15.3484
  hlat(39) = 18.1390
  hlat(40) = 20.9296
  hlat(41) = 23.7202
  hlat(42) = 26.5108
  hlat(43) = 29.3014
  hlat(44) = 32.0919
  hlat(45) = 34.8825
  hlat(46) = 37.6731
  hlat(47) = 40.4636
  hlat(48) = 43.2542
  hlat(49) = 46.0447
  hlat(50) = 48.8352
  hlat(51) = 51.6257
  hlat(52) = 54.4162
  hlat(53) = 57.2066
  hlat(54) = 59.9970
  hlat(55) = 62.7874
  hlat(56) = 65.5776
  hlat(57) = 68.3678
  hlat(58) = 71.1578
  hlat(59) = 73.9475
  hlat(60) = 76.7369
  hlat(61) = 79.5256
  hlat(62) = 82.3129
  hlat(63) = 85.0965
  hlat(64) = 87.8638
  do i = 1 , ilon
    hlon(i) = float(i-1)*2.8125
  end do
!
  sigmar(1) = .01
  sigmar(2) = .03
  sigmar(3) = .05
  sigmar(4) = .07
  sigmar(5) = .1
  sigmar(6) = .15
  sigmar(7) = .2
  sigmar(8) = .25
  sigmar(9) = .3
  sigmar(10) = .4
  sigmar(11) = .5
  sigmar(12) = .7
  sigmar(13) = .85
  sigmar(14) = .925
  sigmar(15) = 1.0
!
!     CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
  do k = 1 , nlev
    kr = nlev - k + 1
    sigma1(k) = sigmar(kr)
  end do

  call getmem3d(b3,1,jx,1,iy,1,nlev*3,'mod_ecwcp:b3')
  call getmem3d(d3,1,jx,1,iy,1,nlev*2,'mod_ecwcp:d3')
  
!     Set up pointers

  u3 => d3(:,:,1:nlev)
  v3 => d3(:,:,nlev+1:2*nlev)
  t3 => b3(:,:,1:nlev)
  h3 => b3(:,:,nlev+1:2*nlev)
  q3 => b3(:,:,2*nlev+1:3*nlev)
  u1 => d2(:,:,1:nlev)
  v1 => d2(:,:,nlev+1:2*nlev)
  t1 => b2(:,:,1:nlev)
  h1 => b2(:,:,nlev+1:2*nlev)
  q1 => b2(:,:,2*nlev+1:3*nlev)

  end subroutine headerec

end module mod_ecwcp
