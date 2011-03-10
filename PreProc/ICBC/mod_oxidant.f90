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

module mod_oxidant

  use mod_dynparam
  use mod_grid
  use mod_wrtoxd
  use mod_interp
  use mod_date
  use m_realkinds
  use m_die
  use m_mall
  use m_zeit

  private
!
  integer, parameter :: ilon = 128, jlat = 64, ilev = 26, itime = 12

  integer :: nyear , month , nday , nhour
  integer :: k , l
  integer :: k0
   
  real(sp), dimension(ilon) ::  t42lon
  real(sp), dimension(jlat) ::  t42lat
  real(sp), dimension(ilev) ::  t42hyam , t42hybm
  real(sp), dimension(ilon,jlat,itime) ::  t42ps
  real(sp), dimension(ilon,ilev,jlat,itime) ::  t42oh , t42ho2 , &
                                        t42o3 , t42no3 , t42h2o2
  real(sp), dimension(ilon,jlat) :: poxid_2
  real(sp), dimension(ilon,jlat,ilev) :: oh2 , ho22 , o32 , no32 , h2o22
  real(sp), allocatable, dimension(:,:)   :: poxid_3
  real(sp), allocatable, dimension(:,:,:) :: oh3 , ho23 , o33 , no33 , h2o23
  real(sp) :: prcm , pmpi , pmpj
  integer :: ncid , istatus

  public  :: headermozart , getmozart , freemozart

  contains

  subroutine headermozart
    use netcdf
    implicit none
    integer :: i , ierr

    allocate(poxid_3(jx,iy), stat=ierr)
    if (ierr /= 0) call die('headermozart','allocate poxid_3',ierr)
    call mall_mci(poxid_3,'mod_oxidant')
    allocate(oh3(jx,iy,ilev), stat=ierr)
    if (ierr /= 0) call die('headermozart','allocate oh3',ierr)
    call mall_mci(oh3,'mod_oxidant')
    allocate(ho23(jx,iy,ilev), stat=ierr)
    if (ierr /= 0) call die('headermozart','allocate ho23',ierr)
    call mall_mci(ho23,'mod_oxidant')
    allocate(o33(jx,iy,ilev), stat=ierr)
    if (ierr /= 0) call die('headermozart','allocate o33',ierr)
    call mall_mci(o33,'mod_oxidant')
    allocate(no33(jx,iy,ilev), stat=ierr)
    if (ierr /= 0) call die('headermozart','allocate no33',ierr)
    call mall_mci(no33,'mod_oxidant')
    allocate(h2o23(jx,iy,ilev), stat=ierr)
    if (ierr /= 0) call die('headermozart','allocate h2o23',ierr)
    call mall_mci(h2o23,'mod_oxidant')

    istatus = nf90_open(trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
                      'oxid_3d_64x128_L26_c030722.nc', nf90_nowrite,ncid)
    if ( istatus /= nf90_noerr ) then
      call die('headermozart','Cannot open input file '// &
               'oxid_3d_64x128_L26_c030722.nc',1,         &
               nf90_strerror(istatus),istatus)
    end if

    do i = 1 , ilon
      t42lon(i) = float(i-1)*2.8125
    end do

    t42lat(1)  = -87.86380
    t42lat(2)  = -85.09653
    t42lat(3)  = -82.31291
    t42lat(4)  = -79.52560
    t42lat(5)  = -76.73690
    t42lat(6)  = -73.94752
    t42lat(7)  = -71.15775
    t42lat(8)  = -68.36776
    t42lat(9)  = -65.57761
    t42lat(10) = -62.78735
    t42lat(11) = -59.99702
    t42lat(12) = -57.20663
    t42lat(13) = -54.41620
    t42lat(14) = -51.62573
    t42lat(15) = -48.83524
    t42lat(16) = -46.04473
    t42lat(17) = -43.25420
    t42lat(18) = -40.46365
    t42lat(19) = -37.67309
    t42lat(20) = -34.88252
    t42lat(21) = -32.09195
    t42lat(22) = -29.30136
    t42lat(23) = -26.51077
    t42lat(24) = -23.72017
    t42lat(25) = -20.92957
    t42lat(26) = -18.13897
    t42lat(27) = -15.34836
    t42lat(28) = -12.55776
    t42lat(29) = -9.767145
    t42lat(30) = -6.976533
    t42lat(31) = -4.185921
    t42lat(32) = -1.395307
    t42lat(33) = 1.395307
    t42lat(34) = 4.185921
    t42lat(35) = 6.976533
    t42lat(36) = 9.767145
    t42lat(37) = 12.55776
    t42lat(38) = 15.34836
    t42lat(39) = 18.13897
    t42lat(40) = 20.92957
    t42lat(41) = 23.72017
    t42lat(42) = 26.51077
    t42lat(43) = 29.30136
    t42lat(44) = 32.09195
    t42lat(45) = 34.88252
    t42lat(46) = 37.67309
    t42lat(47) = 40.46365
    t42lat(48) = 43.25420
    t42lat(49) = 46.04473
    t42lat(50) = 48.83524
    t42lat(51) = 51.62573
    t42lat(52) = 54.41620
    t42lat(53) = 57.20663
    t42lat(54) = 59.99702
    t42lat(55) = 62.78735
    t42lat(56) = 65.57761
    t42lat(57) = 68.36776
    t42lat(58) = 71.15775
    t42lat(59) = 73.94752
    t42lat(60) = 76.73690
    t42lat(61) = 79.52560
    t42lat(62) = 82.31291
    t42lat(63) = 85.09653
    t42lat(64) = 87.86380

    t42hyam(1)  =  3.5446379
    t42hyam(2)  =  7.3888130
    t42hyam(3)  = 13.967210
    t42hyam(4)  = 23.944629
    t42hyam(5)  = 37.230290
    t42hyam(6)  = 53.114600
    t42hyam(7)  = 70.059143
    t42hyam(8)  = 77.912569
    t42hyam(9)  = 76.607011
    t42hyam(10) = 75.071082
    t42hyam(11) = 73.264152
    t42hyam(12) = 71.138389
    t42hyam(13) = 68.637542
    t42hyam(14) = 65.695412
    t42hyam(15) = 62.234160
    t42hyam(16) = 58.162171
    t42hyam(17) = 53.371690
    t42hyam(18) = 47.735929
    t42hyam(19) = 41.105751
    t42hyam(20) = 33.305701
    t42hyam(21) = 24.968440
    t42hyam(22) = 17.095910
    t42hyam(23) = 10.214710
    t42hyam(24) =  4.8031751
    t42hyam(25) =  1.2606850
    t42hyam(26) =  0.0000000

    t42hybm(1)  = 0.0000000
    t42hybm(2)  = 0.0000000
    t42hybm(3)  = 0.0000000
    t42hybm(4)  = 0.0000000
    t42hybm(5)  = 0.0000000
    t42hybm(6)  = 0.0000000
    t42hybm(7)  = 0.0000000
    t42hybm(8)  = 7.5265500E-03 
    t42hybm(9)  = 2.3907680E-02
    t42hybm(10) = 4.3179251E-02
    t42hybm(11) = 6.5851226E-02
    t42hybm(12) = 9.2523657E-02 
    t42hybm(13) = 0.1239024
    t42hybm(14) = 0.1608178
    t42hybm(15) = 0.2042469
    t42hybm(16) = 0.2553391
    t42hybm(17) = 0.3154463
    t42hybm(18) = 0.3861593
    t42hybm(19) = 0.4693495
    t42hybm(20) = 0.5672184
    t42hybm(21) = 0.6718278
    t42hybm(22) = 0.7706061
    t42hybm(23) = 0.8569460
    t42hybm(24) = 0.9248458
    t42hybm(25) = 0.9692941
    t42hybm(26) = 0.9925561
  end subroutine headermozart

  subroutine getmozart(idate)
    use netcdf
    implicit none
!
    integer :: i , j , k , k0
    integer :: idate
    integer, dimension(10) ::  istart , icount

    call zeit_ci('getmozart')

    do i = 4 , 10
      istart(i) = 0
    end do

    do i = 1 , 3
      istart(i) = 1
    end do

    icount(1) = ilon
    icount(2) = jlat
    icount(3) = itime
!C  Retrieve data for Variable 'PS'
    istatus = nf90_get_var(ncid,13,t42PS,istart,icount)

!C  Retrieve data for Variable 'OH', 'HO2', 'O3', 'NO3', 'H2O2'
    istart(4) = 1
    icount(1) = ilon
    icount(2) = ilev
    icount(3) = jlat
    icount(4) = itime

    istatus = nf90_get_var(ncid,14,t42OH,istart,icount)
    istatus = nf90_get_var(ncid,15,t42HO2,istart,icount)
    istatus = nf90_get_var(ncid,16,t42O3,istart,icount)
    istatus = nf90_get_var(ncid,17,t42NO3,istart,icount)
    istatus = nf90_get_var(ncid,18,t42H2O2,istart,icount)

    call split_idate(idate,nyear,month,nday,nhour)

    do j = 1 , jlat
      do i = 1 , ilon
         poxid_2(i,j) = T42PS(i,j,month)*0.01
      end do
    end do


    if ((month == 1 .or. month == 3. .or. month == 5 .or. month == 7 .or. &
         month == 8 .or. month == 10) .and. nday >= 16) then
      if (month == 1) then
        do k = 1 , ilev
          do j = 1 , jlat
            do i = 1 , ilon
              oh2(i,j,k) = T42OH(i,k,j,month)*(1.-float(nday-16)/30.)    &
                          +T42OH(i,k,j,month+1)*(float(nday-16)/30.)  
              ho22(i,j,k) = T42HO2(i,k,j,month)*(1.-float(nday-16)/30.)   &
                          +T42HO2(i,k,j,month+1)*(float(nday-16)/30.)  
              o32(i,j,k) = T42O3(i,k,j,month)*(1.-float(nday-16)/30.)    &
                          +T42O3(i,k,j,month+1)*(float(nday-16)/30.)   
              no32(i,j,k) = T42NO3(i,k,j,month)*(1.-float(nday-16)/30.)   &
                          +T42NO3(i,k,j,month+1)*(float(nday-16)/30.)  
              h2o22(i,j,k) = T42H2O2(i,k,j,month)*(1.-float(nday-16)/30.) &
                           +T42H2O2(i,k,j,month+1)*(float(nday-16)/30.)  
            end do
          end do
        end do
      else if(month == 3 .or. month == 5 .or.   &
              month == 8 .or. month == 10) then                 
        do k = 1 , ilev
          do j = 1 , jlat
            do i = 1 , ilon
              oh2(i,j,k) = T42OH(i,k,j,month)*(1.-float(nday-16)/31.)    &
                          +T42OH(i,k,j,month+1)*(float(nday-16)/31.)   
              ho22(i,j,k) = T42HO2(i,k,j,month)*(1.-float(nday-16)/31.)   &
                          +T42HO2(i,k,j,month+1)*(float(nday-16)/31.)  
              o32(i,j,k) = T42O3(i,k,j,month)*(1.-float(nday-16)/31.)    &
                          +T42O3(i,k,j,month+1)*(float(nday-16)/31.)   
              no32(i,j,k) = T42NO3(i,k,j,month)*(1.-float(nday-16)/31.)   &
                          +T42NO3(i,k,j,month+1)*(float(nday-16)/31.)  
              h2o22(i,j,k) = T42H2O2(i,k,j,month)*(1.-float(nday-16)/31.) &
                           +T42H2O2(i,k,j,month+1)*(float(nday-16)/31.)  
            end do
          end do
        end do
      else if (month == 7) then
        do k = 1 , ilev
          do j = 1 , jlat
            do i = 1 , ilon
              oh2(i,j,k) = T42OH(i,k,j,month)*(1.-float(nday-16)/32.)    &
                          +T42OH(i,k,j,month+1)*(float(nday-16)/32.)   
              ho22(i,j,k) = T42HO2(i,k,j,month)*(1.-float(nday-16)/32.)   &
                          +T42HO2(i,k,j,month+1)*(float(nday-16)/32.)  
              o32(i,j,k) = T42O3(i,k,j,month)*(1.-float(nday-16)/32.)    &
                          +T42O3(i,k,j,month+1)*(float(nday-16)/32.)   
              no32(i,j,k) = T42NO3(i,k,j,month)*(1.-float(nday-16)/32.)   &
                          +T42NO3(i,k,j,month+1)*(float(nday-16)/32.)  
              h2o22(i,j,k) = T42H2O2(i,k,j,month)*(1.-float(nday-16)/32.) &
                           +T42H2O2(i,k,j,month+1)*(float(nday-16)/32.)  
            end do
          end do
        end do
      end if
    else if (month == 2 .and. nday >= 15) then
      do k = 1 , ilev
        do j = 1 , jlat
          do i = 1 , ilon
            oh2(i,j,k) = T42OH(i,k,j,month)*(1.-float(nday-16)/30.)    &
                        +T42OH(i,k,j,month+1)*(float(nday-16)/30.)     
            ho22(i,j,k) = T42HO2(i,k,j,month)*(1.-float(nday-16)/30.)   &
                        +T42HO2(i,k,j,month+1)*(float(nday-16)/30.)    
            o32(i,j,k) = T42O3(i,k,j,month)*(1.-float(nday-16)/30.)    &
                        +T42O3(i,k,j,month+1)*(float(nday-16)/30.)     
            no32(i,j,k) = T42NO3(i,k,j,month)*(1.-float(nday-16)/30.)   &
                        +T42NO3(i,k,j,month+1)*(float(nday-16)/30.)    
            h2o22(i,j,k) = T42H2O2(i,k,j,month)*(1.-float(nday-16)/30.) &
                         +T42H2O2(i,k,j,month+1)*(float(nday-16)/30.)   
          end do
        end do
      end do
    else if((month == 4. .or. month == 6 .or. month == 9 .or.  &
             month == 11) .and. nday >= 16) then
      do k = 1 , ilev
        do j = 1 , jlat
          do i = 1 ,ilon
            oh2(i,j,k) = T42OH(i,k,j,month)*(1.-float(nday-16)/31.)    &
                        +T42OH(i,k,j,month+1)*(float(nday-16)/31.)     
            ho22(i,j,k) = T42HO2(i,k,j,month)*(1.-float(nday-16)/31.)   &
                        +T42HO2(i,k,j,month+1)*(float(nday-16)/31.)
            o32(i,j,k) = T42O3(i,k,j,month)*(1.-float(nday-16)/31.)    &
                        +T42O3(i,k,j,month+1)*(float(nday-16)/31.)
            no32(i,j,k) =  T42NO3(i,k,j,month)*(1.-float(nday-16)/31.)  &
                         +T42NO3(i,k,j,month+1)*(float(nday-16)/31.)
            h2o22(i,j,k) = T42H2O2(i,k,j,month)*(1.-float(nday-16)/31.) &
                         +T42H2O2(i,k,j,month+1)*(float(nday-16)/31.)
          end do
        end do
      end do
    else if (month == 12. .and. nday >= 16) then
      do k = 1 , ilev
        do j = 1 , jlat
          do i = 1 , ilon
            oh2(i,j,k) = T42OH(i,k,j,12)*(1.-float(nday-16)/32.)       &
                        +T42OH(i,k,j,1)*(float(nday-16)/32.)
            ho22(i,j,k) = T42HO2(i,k,j,12)*(1.-float(nday-16)/32.)      &
                        +T42HO2(i,k,j,1)*(float(nday-16)/32.)
            o32(i,j,k) = T42O3(i,k,j,12)*(1.-float(nday-16)/32.)       &
                        +T42O3(i,k,j,1)*(float(nday-16)/32.)
            no32(i,j,k) = T42NO3(i,k,j,12)*(1.-float(nday-16)/32.)      &
                        +T42NO3(i,k,j,1)*(float(nday-16)/32.)
            h2o22(i,j,k) = T42H2O2(i,k,j,12)*(1.-float(nday-16)/32.)    &
                         +T42H2O2(i,k,j,1)*(float(nday-16)/32.)
          end do
        end do
      end do
    else if ((month == 3 .or. month == 5 .or. month == 7 .or.    &
              month == 8 .or. month == 10 .or. month == 12)      &
              .and. nday < 16) then
      if (month == 3) then
        do k = 1 , ilev
          do j = 1 , jlat
            do i = 1 , ilon
              oh2(i,j,k) = T42OH(i,k,j,month-1)*(float(16-nday)/30.)    &
                          +T42OH(i,k,j,month)*(1.-float(16-nday)/30.)
              ho22(i,j,k) = T42HO2(i,k,j,month-1)*(float(16-nday)/30.)   &
                          +T42HO2(i,k,j,month)*(1.-float(16-nday)/30.)  
              o32(i,j,k) = T42O3(i,k,j,month-1)*(float(16-nday)/30.)    &
                          +T42O3(i,k,j,month)*(1.-float(16-nday)/30.)
              no32(i,j,k) = T42NO3(i,k,j,month-1)*(float(16-nday)/30.)   &
                          +T42NO3(i,k,j,month)*(1.-float(16-nday)/30.)
              h2o22(i,j,k) = T42H2O2(i,k,j,month-1)*(float(16-nday)/30.) &
                           +T42H2O2(i,k,j,month)*(1.-float(16-nday)/30.)
            end do
          end do
        end do
      else if (month == 5 .or. month == 7 .or. month == 10 .or. &
               month == 12) then
        do k = 1 , ilev
          do j = 1 , jlat
            do i = 1 , ilon
              oh2(i,j,k) = T42OH(i,k,j,month-1)*(float(16-nday)/31.)    &
                          +T42OH(i,k,j,month)*(1.-float(16-nday)/31.)
              ho22(i,j,k) = T42HO2(i,k,j,month-1)*(float(16-nday)/31.)   &
                          +T42HO2(i,k,j,month)*(1.-float(16-nday)/31.)
              o32(i,j,k) = T42O3(i,k,j,month-1)*(float(16-nday)/31.)    &
                          +T42O3(i,k,j,month)*(1.-float(16-nday)/31.)
              no32(i,j,k) = T42NO3(i,k,j,month-1)*(float(16-nday)/31.)   &
                          +T42NO3(i,k,j,month)*(1.-float(16-nday)/31.)
              h2o22(i,j,k) = T42H2O2(i,k,j,month-1)*(float(16-nday)/31.) &
                           +T42H2O2(i,k,j,month)*(1.-float(16-nday)/31.)
            end do
          end do
        end do
      else if (month == 8) then
        do k = 1 , ilev
          do j = 1 , jlat
            do i = 1 , ilon
              oh2(i,j,k) = T42OH(i,k,j,month-1)*(float(16-nday)/32.)     &
                          +T42OH(i,k,j,month)*(1.-float(16-nday)/32.)
              ho22(i,j,k) = T42HO2(i,k,j,month-1)*(float(16-nday)/32.)    &
                          +T42HO2(i,k,j,month)*(1.-float(16-nday)/32.)
              o32(i,j,k) = T42O3(i,k,j,month-1)*(float(16-nday)/32.)     &
                          +T42O3(i,k,j,month)*(1.-float(16-nday)/32.)
              no32(i,j,k) = T42NO3(i,k,j,month-1)*(float(16-nday)/32.)    &
                          +T42NO3(i,k,j,month)*(1.-float(16-nday)/32.)
              h2o22(i,j,k) = T42H2O2(i,k,j,month-1)*(float(16-nday)/32.)  &
                           +T42H2O2(i,k,j,month)*(1.-float(16-nday)/32.)
            end do
          end do
        end do
      end if
    else if (month == 2 .and. nday < 15) then
      do k = 1 , ilev
        do j = 1 , jlat
          do i = 1 , ilon
            oh2(i,j,k) = T42OH(i,k,j,month-1)*(float(15-nday)/30.)      &
                        +T42OH(i,k,j,month)*(1.-float(15-nday)/30.)    
            ho22(i,j,k) = T42HO2(i,k,j,month-1)*(float(15-nday)/30.)     &
                        +T42HO2(i,k,j,month)*(1.-float(15-nday)/30.)
            o32(i,j,k) = T42O3(i,k,j,month-1)*(float(15-nday)/30.)      &
                        +T42O3(i,k,j,month)*(1.-float(15-nday)/30.)
            no32(i,j,k) = T42NO3(i,k,j,month-1)*(float(15-nday)/30.)     &
                        +T42NO3(i,k,j,month)*(1.-float(15-nday)/30.)
            h2o22(i,j,k) = T42H2O2(i,k,j,month-1)*(float(15-nday)/30.)   &
                         +T42H2O2(i,k,j,month)*(1.-float(15-nday)/30.)
          end do
        end do
      end do
    else if ((month == 4 .or. month == 6 .or. month == 9 .or.      &
              month == 11) .and. nday < 16) then
      do k = 1 , ilev
        do j = 1 , jlat
          do i = 1 , ilon
            oh2(i,j,k) = T42OH(i,k,j,month-1)*(float(16-nday)/31.)      &
                        +T42OH(i,k,j,month)*(1.-float(16-nday)/31.)
            ho22(i,j,k) = T42HO2(i,k,j,month-1)*(float(16-nday)/31.)     &
                        +T42HO2(i,k,j,month)*(1.-float(16-nday)/31.)
            o32(i,j,k) = T42O3(i,k,j,month-1)*(float(16-nday)/31.)      &
                        +T42O3(i,k,j,month)*(1.-float(16-nday)/31.)
            no32(i,j,k) = T42NO3(i,k,j,month-1)*(float(16-nday)/31.)     &
                        +T42NO3(i,k,j,month)*(1.-float(16-nday)/31.)
            h2o22(i,j,k) = T42H2O2(i,k,j,month-1)*(float(16-nday)/31.)   &
                         +T42H2O2(i,k,j,month)*(1.-float(16-nday)/31.)
          end do
        end do
      end do
    else if (month == 1 .and. nday < 16) then
      do k = 1 , ilev
        do j = 1 , jlat
          do i = 1 , ilon
            oh2(i,j,k) = T42OH(i,k,j,12)*(float(16-nday)/32.)   &
                       + T42OH(i,k,j,1)*(1.-float(16-nday)/32.)
            ho22(i,j,k) = T42HO2(i,k,j,12)*(float(16-nday)/32.)   &
                       +T42HO2(i,k,j,1)*(1.-float(16-nday)/32.)
            o32(i,j,k) = T42O3(i,k,j,12)*(float(16-nday)/32.)   &
                       + T42O3(i,k,j,1)*(1.-float(16-nday)/32.)
            no32(i,j,k) = T42NO3(i,k,j,12)*(float(16-nday)/32.)   &
                       +T42NO3(i,k,j,1)*(1.-float(16-nday)/32.)
            h2o22(i,j,k) = T42H2O2(i,k,j,12)*(float(16-nday)/32.) &
                      +T42H2O2(i,k,j,1)*(1.-float(16-nday)/32.)
          end do
        end do
      end do
    end if

    call bilinx2(oh3,oh2,xlon,xlat,t42lon,t42lat,ilon,jlat,iy,jx,ilev)
    call bilinx2(ho23,ho22,xlon,xlat,t42lon,t42lat,ilon,jlat,iy,jx,ilev)
    call bilinx2(o33,o32,xlon,xlat,t42lon,t42lat,ilon,jlat,iy,jx,ilev)
    call bilinx2(no33,no32,xlon,xlat,t42lon,t42lat,ilon,jlat,iy,jx,ilev)
    call bilinx2(h2o23,h2o22,xlon,xlat,t42lon,t42lat,ilon,jlat,iy,jx,ilev)
    call bilinx2(poxid_3,poxid_2,xlon,xlat,t42lon,t42lat,ilon,jlat,iy,jx,1)
    do i = 1 , iy
      do j = 1 , jx
        do l = 1 , kz
          prcm = ((poxid_3(j,i)*0.1-real(ptop))*sigma2(l)+real(ptop))*10.

          k0 = -1
          do k = ilev , 1 , -1
            pmpi = poxid_3(j,i)*t42hybm(k)+t42hyam(k)
            k0 = k
            if (prcm > pmpi) exit
          end do
          if (k0 == ilev) then
            pmpi = poxid_3(j,i)*t42hybm(ilev)+t42hyam(ilev)
            oh4(j,i,l) = oh3(j,i,ilev)+(oh3(j,i,ilev)-oh3(j,i,ilev-1))         &
                        *(prcm-pmpi)/(pmpi-pmpj)
            ho24(j,i,l) = ho23(j,i,ilev)+(ho23(j,i,ilev)-ho23(j,i,ilev-1))     &
                         *(prcm-pmpi)/(pmpi-pmpj)
            o34(j,i,l) = o33(j,i,ilev)+(o33(j,i,ilev)-o33(j,i,ilev-1))         &
                        *(prcm-pmpi)/(pmpi-pmpj)
            no34(j,i,l) = no33(j,i,ilev)+(no33(j,i,ilev)-no33(j,i,ilev-1))     &
                         *(prcm-pmpi)/(pmpi-pmpj)
            h2o24(j,i,l) = h2o23(j,i,ilev)+(h2o23(j,i,ilev)-h2o23(j,i,ilev-1)) &
                          *(prcm-pmpi)/(pmpi-pmpj)
          else if (k0 >= 1) then
            pmpi = poxid_3(j,i)*t42hybm(k0+1)+t42hyam(k0+1)
            oh4(j,i,l) = (oh3(j,i,k0+1)*(prcm-pmpj)            &
                          +oh3(j,i,k0)*(pmpi-prcm))/(pmpi-pmpj)
            ho24(j,i,l) = (ho23(j,i,k0+1)*(prcm-pmpj)          &
                           +ho23(j,i,k0)*(pmpi-prcm))/(pmpi-pmpj)
            o34(j,i,l) = (o33(j,i,k0+1)*(prcm-pmpj)            &
                          +o33(j,i,k0)*(pmpi-prcm))/(pmpi-pmpj)
            no34(j,i,l) = (no33(j,i,k0+1)*(prcm-pmpj)          &
                           +no33(j,i,k0)*(pmpi-prcm))/(pmpi-pmpj)
            h2o24(j,i,l) = (h2o23(j,i,k0+1)*(prcm-pmpj)        &
                            +h2o23(j,i,k0)*(pmpi-prcm))/(pmpi-pmpj)
          endif
        end do
      end do
    end do
    call writeox(idate)
    call zeit_co('getmozart')
  end subroutine getmozart

  subroutine freemozart
    use netcdf
    implicit none
    istatus = nf90_close(ncid)
    if ( istatus /= nf90_noerr ) then
      call die('freemozart','Cannot close input file',1, &
               nf90_strerror(istatus),istatus)
    end if
    call mall_mco(poxid_3,'mod_oxidant')
    deallocate(poxid_3)
    call mall_mco(oh3,'mod_oxidant')
    deallocate(oh3)
    call mall_mco(ho23,'mod_oxidant')
    deallocate(ho23)
    call mall_mco(o33,'mod_oxidant')
    deallocate(o33)
    call mall_mco(no33,'mod_oxidant')
    deallocate(no33)
    call mall_mco(h2o23,'mod_oxidant')
    deallocate(h2o23)
  end subroutine freemozart

end module mod_oxidant
