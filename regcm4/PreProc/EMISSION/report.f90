      program report
      use mod_regcm_param , only : ntr ,  nbin
      use mod_emission
      implicit none
!
! Local variables
!
      character(8) , dimension(ntr) :: chtrname
      real , dimension(nbin) :: dustbsiz
!
!     This program produce a report have all information about
!     emission species and transported species
!
      namelist /species/ ele_retroa , ele_retrob , ele_poeta ,          &
      & ele_poetb , ele_gfed , ele_edgara , ele_edgarb , ele_mozrta ,   &
      & ele_mozrtb
      namelist /transport/ chtrname
      namelist /dust_bins/ dustbsiz

      open (10,file='namelist.report',status='old')
      read (10,nml=species)
      read (10,nml=transport)
      read (10,nml=dust_bins)
 
      write (*,*) ele_retroa
      write (*,*) ele_edgara
      write (*,*) chtrname
      write (*,*) dustbsiz

      end program report
