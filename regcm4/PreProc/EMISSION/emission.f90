      program emission
!#####################################################################
!#    This program read the emission species from different emission
!# inventories, and interpolat the data to model grid.
!# all data written in NetCDF format, so the program use the
!# netcdf library.
!# the following is the definition of the variables inside.
!#####################################################################
! the program depends on three subroutine
! 1- JULIAN2 to transform the DATE format (e.g 2004020100) to its original
!    components of year, month, day and hour
! 2- GRIDDEFF to extract some data related to the model GRID from
!    DOMAIN.INFO file and to creat AERO.ctl file
! 3- reademission to read emission in netcdf format from the data base
!    and write it to AERO.dat file
!######################################################################
! AER1(NLONS1,NLATS1):the array of requested element
!      NLONS1=720, NLATS1=360
! AER2(NLONS2,NLATS2):the arrray of requested element
!      NLONS2=360, NLATS2=180
! AERMM(IY,JX) : the output array after interpolation to model GRID
 
! LONI1(NLONS1) : the array of longitude inside netcdf file
! LONI2(NLONS2) : the array of longitude inside netcdf file
! LATI1(NLATS1) : the array of latitude inside netcdf file
! LATI2(NLATS2) : the array of latitude inside netcdf file
! we retrieve many variable from Input/DOMAIN.INFO file
! XLON(IY,JX)   : the array of model longitude written in DOMAIN.INFO
! XLAT(IY,JX)   : the array of model latitude  written in DOMAIN.INFO
! TRUELATL, TRUELATH : written in DOMAIN.INFO
! record : to trace the record number in the AERO.dat output file
! julday,julncep,jul1900 : to transform the data format(e.g 2005030100) to
! its original components of year, month, day, hour
! iyr1, imo1, idy1, ihr1, ayr1 for IDATE1
! iyr2, imo2, idy2, ihr2, ayr2 for IDATE2
! TREC : the total number of records per element
! (e.g if we run for 5 month, then there are 5 record for each elements)
 
! EMA=1  for  anthropogenic emission
! EMB=2  for  biogenic emission
 
! a   for RETRO inventory
! b   for POET inventory
! c   for GFED inventory
! d   for EDGAR inventory
 
! INVNTRY1='RETRO'
! INVNTRY2='POET'
! INVNTRY3='GFED'
! INVNTRY4='EDGAR'
 
! EMSRC1=1 : read anthrobogenic emission
! EMSRC2=2 : read biogenic emission
 
! NSPC1A   : number of Anthropogenic species from RETRO
! NSPC1B   : number of Biogenic species from RETRO
! NSPC2    : number of Anthropogenic or biogenic species from POET
! NSPC3    : number of Biogenic species from GFEDV2
! NSPC4A   : number of Anthropogenic species from EDGAR
! NSPC4B   : number of Biogenic species from EDGAR
 
! DR1      : the directory of RETRO data
! DR2      : the directory of POET data
! DR3      : the directory of GFED data
! DR4      : the directory of EDGAR data
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mod_param
      use mod_emission
      implicit none
!
! Local variables
!
      real(4) , dimension(nlons1,nlats1) :: aer1
      real(4) , dimension(nlons2,nlats2) :: aer2
      real(4) , dimension(nlats1) :: lati1
      real(4) , dimension(nlats2) :: lati2
      real(4) , dimension(nlons1) :: loni1
      real(4) , dimension(nlons2) :: loni2

      integer :: ayr1 , ayr2 , dt , idy1 , idy2 , ihr1 , ihr2 , imo1 ,  &
               & imo2 , iyr1 , iyr2 , jul1900 , julday , julncep ,      &
               & mm , mm1 , mm2 , mm22 , mm3 , mm4 , mm5 , mone ,       &
               & mons , record , trec , yr , yr1 , yr2 , yr22 ,  yr3 ,  &
               & yr4 , yr5
      logical :: there
 
      open (30,file='namelist.input',status='old')
!     READ(30,NML=SPECIES)
 
      inquire (file='../../Input/AERO_new.dat',exist=there)
      if ( there ) call unlink('../../Input/AERO.dat')
 
!     ******    ON WHAT RegCM GRID ARE AEROSOL DESIRED?
      open (10,file='../ICBC/fort.10',form='unformatted',               &
          & recl=ix*jx*ibyte,access='direct',status='old',err=100)
!     &   ,convert='big_endian' ,status='unknown',ERR=4830)
      inquire (file='../../Input/AERO.dat',exist=there)
      if ( there ) call unlink('../../Input/oxid.dat')
!     OPEN(55,file='../../Input/oxid.dat',form='unformatted'
!     &       ,recl=IY*JX*ibyte,access='direct')
 
      dt = 1
      call julian2(idate1,dt,julday,julncep,jul1900,iyr1,imo1,idy1,ihr1,&
                 & ayr1)
      call julian2(idate2,dt,julday,julncep,jul1900,iyr2,imo2,idy2,ihr2,&
                 & ayr2)
 
      write (*,*) 'START YEAR=' , iyr1 , 'END YEAR=' , iyr2
      write (*,*) 'START MONTH=' , imo1 , 'END MONTHE=' , imo2
      trec = (iyr2-iyr1+1)*12 - (imo1-1) - (12-imo2)
      write (*,*) 'TOTAL RECORD NUMBER per element=' , trec
 
      call griddef(iyr2,trec)
 
      record = 0
 
      do yr = iyr1 , iyr2  !loop over years
        mons = 1
        mone = 12
        if ( yr==iyr1 ) mons = imo1
        if ( yr==iyr2 ) mone = imo2
        do mm = mons , mone !loop over months
          yr1 = yr
          yr2 = yr
          yr22 = yr
          yr3 = yr
          yr4 = yr
          yr5 = yr
          mm1 = mm
          mm2 = mm
          mm22 = mm
          mm3 = mm
          mm4 = mm
          mm5 = mm
          select case (aertyp)
          case ('AER01D0','AER01D1')
!           in this case we read biogenic from EDGAR, RETRO, GFED
            write (*,*) 'HI BIOGENIC'
            if ( iretro ) call reademission(nlats1,nlons1,nlvls1,'a',   &
                &'DR1',emb,nspc1b,aer1,loni1,lati1,yr1,mm1,record)
            if ( ipoet ) call reademission(nlats2,nlons2,nlvls2,'b',    &
                &'DR2',emb,nspc2b,aer2,loni2,lati2,yr2,mm2,record)
 
            if ( yr>2000 ) then
              if ( igfed ) call reademission(nlats2,nlons2,nlvls2,'c',  &
                  &'DR3',emb,nspc3,aer2,loni2,lati2,yr3,mm3,record)
            end if
            if ( iedgar ) call reademission(nlats2,nlons2,nlvls2,'d',   &
                &'DR4',emb,nspc4b,aer2,loni2,lati2,yr5,mm5,record)
            if ( imozart ) call reademission(nlats2,nlons2,nlvls2,'m',  &
                &'DR5',emb,nspc5b,aer2,loni2,lati2,yr5,mm5,record)
          case ('AER10D0','AER10D1')
!           in this case we read anthropogenic from RETRO, EDGAR
            write (*,*) 'HI ANTHROPOGENIC'
            if ( iretro ) call reademission(nlats1,nlons1,nlvls1,'a',   &
                &'DR1',ema,nspc1a,aer1,loni1,lati1,yr1,mm1,record)
            if ( ipoet ) call reademission(nlats2,nlons2,nlvls2,'b',    &
                &'DR2',ema,nspc2a,aer2,loni2,lati2,yr2,mm2,record)
            if ( iedgar ) call reademission(nlats2,nlons2,nlvls2,'d',   &
                &'DR4',ema,nspc4a,aer2,loni2,lati2,yr5,mm5,record)
            if ( imozart ) call reademission(nlats2,nlons2,nlvls2,'m',  &
                &'DR5',emb,nspc5a,aer2,loni2,lati2,yr5,mm5,record)
          case ('AER11D0','AER11D1')
            write (*,*) 'HI BOTH'
!--------------------Anthropogenic data ----------------------
            if ( iretro ) call reademission(nlats1,nlons1,nlvls1,'a',   &
                &'DR1',ema,nspc1a,aer1,loni1,lati1,yr1,mm1,record)
 
            if ( ipoet ) then
              write (*,*) 'hi poet' , yr2 , mons
              call reademission(nlats2,nlons2,nlvls2,'b','DR2',ema,     &
                              & nspc2a,aer2,loni2,lati2,yr2,mm2,record)
              write (*,*) 'by poet' , yr2 , mons
            end if
            if ( iedgar ) then
              write (*,*) 'HI EDGAR'
              call reademission(nlats2,nlons2,nlvls2,'d','DR4',ema,     &
                              & nspc4a,aer2,loni2,lati2,yr5,mm5,record)
            end if
            if ( imozart ) call reademission(nlats2,nlons2,nlvls2,'m',  &
                &'DR5',ema,nspc5a,aer2,loni2,lati2,yr5,mm5,record)
!-------------------Biogenic data--------------------------------
            write (*,*) '---------------yy' , yr1 , yr
            if ( iretro ) call reademission(nlats1,nlons1,nlvls1,'a',   &
                &'DR1',emb,nspc1b,aer1,loni1,lati1,yr1,mm1,record)
            if ( ipoet ) call reademission(nlats2,nlons2,nlvls2,'b',    &
                &'DR2',emb,nspc2b,aer2,loni2,lati2,yr22,mm22,record)
 
            if ( yr>2000 ) then
              if ( igfed ) call reademission(nlats2,nlons2,nlvls2,'c',  &
                  &'DR3',emb,nspc3,aer2,loni2,lati2,yr3,mm3,record)
            end if
            if ( iedgar ) call reademission(nlats2,nlons2,nlvls2,'d',   &
                &'DR4',emb,nspc4b,aer2,loni2,lati2,yr5,mm5,record)
            if ( imozart ) call reademission(nlats2,nlons2,nlvls2,'m',  &
                &'DR5',emb,nspc5b,aer2,loni2,lati2,yr5,mm5,record)
          case default
          end select
        end do     !end of month's
      end do   !end of years
 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
      stop 99999

      print * , 'ERROR OPENING AEROSOL FILE'
      stop '4810 IN PROGRAM AEROSOL'
 100  continue
      print * , 'ERROR OPENING DOMAIN HEADER FILE'
      stop '4830 IN PROGRAM RDSST'
      end program emission
