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

module mod_erahi

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_message
  use mod_stdio
  use mod_grid
  use mod_write
  use mod_interp
  use mod_vertint
  use mod_hgt
  use mod_humid
  use mod_mksst
  use mod_uvrot
  use mod_vectutil

  private

  integer(ik4) , parameter :: nlev1 = 60
  integer(ik4) , parameter :: nlats = 160
  integer(ik4) , parameter :: nlons = 320
  integer(ik4) , parameter :: nlev2 = 18

  real(rk8) , dimension(nlev1+1) :: ak , bk
  real(rk8) , dimension(nlev2) :: pplev , sigma1 , sigmar
  real(rk8) , dimension(nlats) :: slat
  real(rk8) , dimension(nlons) :: slon

  real(rk4) , dimension(nlons,nlats) :: iobuf
  real(rk8) , dimension(nlons,nlats) :: lsm , ps2 , zs2
  real(rk8) , dimension(nlons,nlats,nlev1) :: q2 , t2 , u2 , v2
  real(rk8) , dimension(nlons,nlats,nlev1) :: pp3d , z1

  real(rk8) , target , dimension(nlons,nlats,nlev2*3) :: b2
  real(rk8) , target , dimension(nlons,nlats,nlev2*2) :: d2
  real(rk8) , pointer , dimension(:,:,:) :: b3
  real(rk8) , pointer , dimension(:,:,:) :: d3

  real(rk8) , pointer , dimension(:,:,:) :: tp , qp , hp
  real(rk8) , pointer , dimension(:,:,:) :: up , vp
  real(rk8) , pointer , dimension(:,:,:) :: t3 , q3 , h3
  real(rk8) , pointer , dimension(:,:,:) :: u3 , v3

  public :: geterahi , headerehi

  contains

  subroutine geterahi(idate)
  implicit none
!
  type(rcm_time_and_date) , intent(in) :: idate
!
  character(len=256) :: finame
  integer(ik4) :: i , j , k , nrec
  logical :: there
  real(rk8) :: slonmax , slonmin , xlonmax , xlonmin
  integer :: hireclen
!
  if ( idate == globidate1 ) then
    xlonmin = 400.
    xlonmax = -400.
    do j = 1 , iy
      do i = 1 , jx
        if ( xlon(i,j) < xlonmin ) xlonmin = xlon(i,j)
        if ( xlon(i,j) > xlonmax ) xlonmax = xlon(i,j)
      end do
    end do
    write (stdout,*) 'XLONMIN,XLONMAX = ' , xlonmin , xlonmax
    slonmin = 400.
    slonmax = -400.
    do i = 1 , nlons
      if ( slon(i) < slonmin ) slonmin = slon(i)
      if ( slon(i) > slonmax ) slonmax = slon(i)
    end do
    write (stdout,*) 'SLONMIN,SLONMAX = ' , slonmin , slonmax
  end if
  write (finame,99001) trim(inpglob),pthsep,'ERAHI',pthsep,toint10(idate)
  inquire (file=finame,exist=there)
  if ( .not. there ) then
    call die('ERAHI', trim(finame)//' is not available',1)
  end if
  inquire(iolength=hireclen) iobuf
  open (61,file=finame,form='unformatted',recl=hireclen,action='read', &
        access='direct',status='old')
  nrec = 0
  nrec = nrec + 1
  read (61,rec=nrec) iobuf
  zs2 = iobuf
  nrec = nrec + 1
  read (61,rec=nrec) iobuf
  lsm = iobuf
  nrec = nrec + 1
  read (61,rec=nrec) iobuf
  ps2 = iobuf
  do j = 1 , nlats
    do i = 1 , nlons
      if ( lsm(i,j) < 0.5 ) zs2(i,j) = 0.000
    end do
  end do
  do k = 1 , nlev1
    nrec = nrec + 1
    read (61,rec=nrec) iobuf
    t2(:,:,k) = iobuf
  end do
  do k = 1 , nlev1
    nrec = nrec + 1
    read (61,rec=nrec) iobuf
    q2(:,:,k) = iobuf
  end do
  do k = 1 , nlev1
    nrec = nrec + 1
    read (61,rec=nrec) iobuf
    u2(:,:,k) = iobuf
  end do
  do k = 1 , nlev1
    nrec = nrec + 1
    read (61,rec=nrec) iobuf
    v2(:,:,k) = iobuf
  end do

  write (stdout,*) 'READ IN fields at DATE:' , tochar(idate)
  do k = 1 , nlev1
    do j = 1 , nlats
      do i = 1 , nlons
        if ( ps2(i,j) > -9995. ) then
          pp3d(i,j,k) = ps2(i,j)*0.5*(bk(k)+bk(k+1))   &
                        + 0.5*(ak(k)+ak(k+1))
        else
          pp3d(i,j,k) = -9999.0
        end if
      end do
    end do
  end do
!
!     to calculate Heights on sigma surfaces.
  call htsig(t2,z1,pp3d,ps2,zs2,nlons,nlats,nlev1)
!
!     to interpolate H,U,V,T,Q and QC
!     1. For Heights
  call height(hp,z1,t2,ps2,pp3d,zs2,nlons,nlats,nlev1,pplev,nlev2)
!     2. For Zonal and Meridional Winds
  call intlin(up,u2,ps2,pp3d,nlons,nlats,nlev1,pplev,nlev2)
  call intlin(vp,v2,ps2,pp3d,nlons,nlats,nlev1,pplev,nlev2)
!     3. For Temperatures
  call intlog(tp,t2,ps2,pp3d,nlons,nlats,nlev1,pplev,nlev2)
!     4. For Moisture qva & qca
  call humid1fv(t2,q2,pp3d,nlons,nlats,nlev1)
  call intlin(qp,q2,ps2,pp3d,nlons,nlats,nlev1,pplev,nlev2)
!
!     HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
  call bilinx2(b3,b2,xlon,xlat,slon,slat,nlons,nlats,jx,iy,nlev2*3)
  call bilinx2(d3,d2,dlon,dlat,slon,slat,nlons,nlats,jx,iy,nlev2*2)
!
!     ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
  call uvrot4(u3,v3,dlon,dlat,clon,clat,xcone,jx,iy,nlev2,plon,plat,iproj)
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!     V E R T I C A L   I N T E R P O L A T I O N
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!HH:  CHANGE THE VERTICAL ORDER.
  call top2btm(t3,jx,iy,nlev2)
  call top2btm(q3,jx,iy,nlev2)
  call top2btm(h3,jx,iy,nlev2)
  call top2btm(u3,jx,iy,nlev2)
  call top2btm(v3,jx,iy,nlev2)
!HH:OVER
!
!     ******           NEW CALCULATION OF P* ON RegCM TOPOGRAPHY.
  call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,nlev2)

  call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
  if(i_band == 1) then
     call p1p2_band(b3pd,ps4,jx,iy)
  else
     call p1p2(b3pd,ps4,jx,iy)
  endif
!
!     F0    DETERMINE SURFACE TEMPS ON RegCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
  call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,nlev2)

  call readsst(ts4,idate)

!     F3     INTERPOLATE U, V, T, AND Q.
  call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,nlev2)
  call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,nlev2)
!
  call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,iy,kz,nlev2)

  call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,iy,kz,nlev2)
  call humid2(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
!
!     F4     DETERMINE H
  call hydrost(h4,t4,topogm,ps4,ptop,sigma2,jx,iy,kz)
!
!
99001 format (a,a,a,a,'EHI_',i10)
!
  end subroutine geterahi

!
  subroutine headerehi
  implicit none
!
  integer(ik4) :: i , k , kr
!
  slat(1) = -89.142
  slat(2) = -88.029
  slat(3) = -86.911
  slat(4) = -85.791
  slat(5) = -84.670
  slat(6) = -83.549
  slat(7) = -82.428
  slat(8) = -81.307
  slat(9) = -80.185
  slat(10) = -79.064
  slat(11) = -77.943
  slat(12) = -76.821
  slat(13) = -75.700
  slat(14) = -74.578
  slat(15) = -73.457
  slat(16) = -72.336
  slat(17) = -71.214
  slat(18) = -70.093
  slat(19) = -68.971
  slat(20) = -67.850
  slat(21) = -66.728
  slat(22) = -65.607
  slat(23) = -64.485
  slat(24) = -63.364
  slat(25) = -62.242
  slat(26) = -61.121
  slat(27) = -60.000
  slat(28) = -58.878
  slat(29) = -57.757
  slat(30) = -56.635
  slat(31) = -55.514
  slat(32) = -54.392
  slat(33) = -53.271
  slat(34) = -52.149
  slat(35) = -51.028
  slat(36) = -49.906
  slat(37) = -48.785
  slat(38) = -47.663
  slat(39) = -46.542
  slat(40) = -45.420
  slat(41) = -44.299
  slat(42) = -43.177
  slat(43) = -42.056
  slat(44) = -40.934
  slat(45) = -39.813
  slat(46) = -38.691
  slat(47) = -37.570
  slat(48) = -36.448
  slat(49) = -35.327
  slat(50) = -34.205
  slat(51) = -33.084
  slat(52) = -31.962
  slat(53) = -30.841
  slat(54) = -29.719
  slat(55) = -28.598
  slat(56) = -27.476
  slat(57) = -26.355
  slat(58) = -25.234
  slat(59) = -24.112
  slat(60) = -22.991
  slat(61) = -21.869
  slat(62) = -20.748
  slat(63) = -19.626
  slat(64) = -18.505
  slat(65) = -17.383
  slat(66) = -16.262
  slat(67) = -15.140
  slat(68) = -14.019
  slat(69) = -12.897
  slat(70) = -11.776
  slat(71) = -10.654
  slat(72) = -9.533
  slat(73) = -8.411
  slat(74) = -7.290
  slat(75) = -6.168
  slat(76) = -5.047
  slat(77) = -3.925
  slat(78) = -2.804
  slat(79) = -1.682
  slat(80) = -0.561
  slat(81) = 0.561
  slat(82) = 1.682
  slat(83) = 2.804
  slat(84) = 3.925
  slat(85) = 5.047
  slat(86) = 6.168
  slat(87) = 7.290
  slat(88) = 8.411
  slat(89) = 9.533
  slat(90) = 10.654
  slat(91) = 11.776
  slat(92) = 12.897
  slat(93) = 14.019
  slat(94) = 15.140
  slat(95) = 16.262
  slat(96) = 17.383
  slat(97) = 18.505
  slat(98) = 19.626
  slat(99) = 20.748
  slat(100) = 21.869
  slat(101) = 22.991
  slat(102) = 24.112
  slat(103) = 25.234
  slat(104) = 26.355
  slat(105) = 27.476
  slat(106) = 28.598
  slat(107) = 29.719
  slat(108) = 30.841
  slat(109) = 31.962
  slat(110) = 33.084
  slat(111) = 34.205
  slat(112) = 35.327
  slat(113) = 36.448
  slat(114) = 37.570
  slat(115) = 38.691
  slat(116) = 39.813
  slat(117) = 40.934
  slat(118) = 42.056
  slat(119) = 43.177
  slat(120) = 44.299
  slat(121) = 45.420
  slat(122) = 46.542
  slat(123) = 47.663
  slat(124) = 48.785
  slat(125) = 49.906
  slat(126) = 51.028
  slat(127) = 52.149
  slat(128) = 53.271
  slat(129) = 54.392
  slat(130) = 55.514
  slat(131) = 56.635
  slat(132) = 57.757
  slat(133) = 58.878
  slat(134) = 60.000
  slat(135) = 61.121
  slat(136) = 62.242
  slat(137) = 63.364
  slat(138) = 64.485
  slat(139) = 65.607
  slat(140) = 66.728
  slat(141) = 67.850
  slat(142) = 68.971
  slat(143) = 70.093
  slat(144) = 71.214
  slat(145) = 72.336
  slat(146) = 73.457
  slat(147) = 74.578
  slat(148) = 75.700
  slat(149) = 76.821
  slat(150) = 77.943
  slat(151) = 79.064
  slat(152) = 80.185
  slat(153) = 81.307
  slat(154) = 82.428
  slat(155) = 83.549
  slat(156) = 84.670
  slat(157) = 85.791
  slat(158) = 86.911
  slat(159) = 88.029
  slat(160) = 89.142

  do i = 1 , nlons
    slon(i) = float(i-1)*1.125
  end do

  pplev(1) = 30.
  pplev(2) = 50.
  pplev(3) = 70.
  pplev(4) = 100.
  pplev(5) = 150.
  pplev(6) = 200.
  pplev(7) = 250.
  pplev(8) = 300.
  pplev(9) = 350.
  pplev(10) = 420.
  pplev(11) = 500.
  pplev(12) = 600.
  pplev(13) = 700.
  pplev(14) = 780.
  pplev(15) = 850.
  pplev(16) = 920.
  pplev(17) = 960.
  pplev(18) = 1000.

  do k = 1 , nlev2
    sigmar(k) = pplev(k)*0.001
  end do
!HH:OVER
!     CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
  do k = 1 , nlev2
    kr = nlev2 - k + 1
    sigma1(k) = sigmar(kr)
  end do

  ak(1) = 0.00000000
  ak(2) = 0.20000000
  ak(3) = 0.38425343
  ak(4) = 0.63647804
  ak(5) = 0.95636963
  ak(6) = 1.34483307
  ak(7) = 1.80584351
  ak(8) = 2.34779053
  ak(9) = 2.98495789
  ak(10) = 3.73971924
  ak(11) = 4.64618134
  ak(12) = 5.75651001
  ak(13) = 7.13218079
  ak(14) = 8.83660522
  ak(15) = 10.94834717
  ak(16) = 13.56474609
  ak(17) = 16.80640259
  ak(18) = 20.82273926
  ak(19) = 25.79888672
  ak(20) = 31.96421631
  ak(21) = 39.60291504
  ak(22) = 49.06708496
  ak(23) = 60.18019531
  ak(24) = 73.06631348
  ak(25) = 87.65053711
  ak(26) = 103.76126953
  ak(27) = 120.77446289
  ak(28) = 137.75325195
  ak(29) = 153.79805664
  ak(30) = 168.19474609
  ak(31) = 180.45183594
  ak(32) = 190.27695313
  ak(33) = 197.55109375
  ak(34) = 202.22205078
  ak(35) = 204.29863281
  ak(36) = 203.84480469
  ak(37) = 200.97402344
  ak(38) = 195.84330078
  ak(39) = 188.64750000
  ak(40) = 179.61357422
  ak(41) = 168.99468750
  ak(42) = 157.06447266
  ak(43) = 144.11124023
  ak(44) = 130.43218750
  ak(45) = 116.32758789
  ak(46) = 102.09500977
  ak(47) = 88.02356445
  ak(48) = 74.38803223
  ak(49) = 61.44314941
  ak(50) = 49.41778320
  ak(51) = 38.50913330
  ak(52) = 28.87696533
  ak(53) = 20.63779785
  ak(54) = 13.85912598
  ak(55) = 8.55361755
  ak(56) = 4.67333588
  ak(57) = 2.10393890
  ak(58) = 0.65889244
  ak(59) = 0.07367743
  ak(60) = 0.00000000
  ak(61) = 0.00000000

  bk(1) = 0.00000000
  bk(2) = 0.00000000
  bk(3) = 0.00000000
  bk(4) = 0.00000000
  bk(5) = 0.00000000
  bk(6) = 0.00000000
  bk(7) = 0.00000000
  bk(8) = 0.00000000
  bk(9) = 0.00000000
  bk(10) = 0.00000000
  bk(11) = 0.00000000
  bk(12) = 0.00000000
  bk(13) = 0.00000000
  bk(14) = 0.00000000
  bk(15) = 0.00000000
  bk(16) = 0.00000000
  bk(17) = 0.00000000
  bk(18) = 0.00000000
  bk(19) = 0.00000000
  bk(20) = 0.00000000
  bk(21) = 0.00000000
  bk(22) = 0.00000000
  bk(23) = 0.00000000
  bk(24) = 0.00000000
  bk(25) = 0.00007582
  bk(26) = 0.00046139
  bk(27) = 0.00181516
  bk(28) = 0.00508112
  bk(29) = 0.01114291
  bk(30) = 0.02067788
  bk(31) = 0.03412116
  bk(32) = 0.05169041
  bk(33) = 0.07353383
  bk(34) = 0.09967469
  bk(35) = 0.13002251
  bk(36) = 0.16438432
  bk(37) = 0.20247594
  bk(38) = 0.24393314
  bk(39) = 0.28832296
  bk(40) = 0.33515489
  bk(41) = 0.38389215
  bk(42) = 0.43396294
  bk(43) = 0.48477158
  bk(44) = 0.53570992
  bk(45) = 0.58616841
  bk(46) = 0.63554746
  bk(47) = 0.68326861
  bk(48) = 0.72878581
  bk(49) = 0.77159661
  bk(50) = 0.81125343
  bk(51) = 0.84737492
  bk(52) = 0.87965691
  bk(53) = 0.90788388
  bk(54) = 0.93194032
  bk(55) = 0.95182151
  bk(56) = 0.96764523
  bk(57) = 0.97966272
  bk(58) = 0.98827010
  bk(59) = 0.99401945
  bk(60) = 0.99763012
  bk(61) = 1.00000000

  call getmem3d(b3,1,jx,1,iy,1,nlev2*3,'mod_erahi:b3')
  call getmem3d(d3,1,jx,1,iy,1,nlev2*2,'mod_erahi:b3')

!     Set up pointers

  tp => b2(:,:,1:nlev2)
  qp => b2(:,:,nlev2+1:2*nlev2)
  hp => b2(:,:,2*nlev2+1:3*nlev2)
  up => d2(:,:,1:nlev2)
  vp => d2(:,:,nlev2+1:2*nlev2)
  t3 => b3(:,:,1:nlev2)
  q3 => b3(:,:,nlev2+1:2*nlev2)
  h3 => b3(:,:,2*nlev2+1:3*nlev2)
  u3 => d3(:,:,1:nlev2)
  v3 => d3(:,:,nlev2+1:2*nlev2)

  end subroutine headerehi

end module mod_erahi
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
