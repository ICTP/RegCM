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

module mod_runparams

  use mod_realkinds
  use mod_date
  use mod_dynparam
  use mod_memutil

  implicit none

  integer , public , parameter :: nqx = 2
  integer , public , parameter :: iqv = 1
  integer , public , parameter :: iqc = 2

  type(rcm_time_and_date) , save :: idate0 , idate1 , idate2

  type(rcm_time_and_date) , save :: idatex
  integer :: xyear , xmonth , xday , xhour

  type(rcm_time_and_date) , save :: bdydate1 , bdydate2

  type(rcm_time_interval) , save :: intmdl
  type(rcm_time_interval) , save :: intbdy

  real(dp) :: declin , deltmx
  real(dp) :: xbctime
  real(dp) :: calday , twodt

  ! Step counter. Is zero at idate0, always increasing, never reset.
  integer(8) :: ktau
  ! Final number of step for THIS run
  integer(8) :: mtau
  ! How many steps for an hour (updates date fields Y m d H)
  integer(8) :: khour
  ! Output k values for I/O operations.
  integer(8) :: katm , krad , kche , ksav , kdbg , kbdy , ksrf
  ! Seconds counter in between boundary conditions read
  integer(8) :: nbdytime
  ! Step counters to activate surface and radiation schemes
  integer(8) :: ntsrf , ntrad
  ! Model timestep in seconds (real and integer)
  integer(8) :: ntsec
  real(dp) :: dtsec
  ! Internal count for how many SRF outputs every LAK output
  integer :: klak
  ! Internal count for how many SRF outputs per day
  integer(8) :: ksts , kstsoff
  !
  ! Cumulus scheme index
  integer :: icup
  ! Closure index for Grell
  integer :: igcc
  ! Boundary layer index
  integer :: ibltyp
  ! Lake model activation index
  integer :: lakemod
  ! Diurnal cycle SST index
  integer :: idcsst
  ! Sea Ice scheme index
  integer :: iseaice
  ! Seasonal albedo for desert index
  integer :: idesseas
  ! Ocean model switch indexes
  integer :: iocnrough , iocnflx , iocncpl
  ! Radiation switch controls
  integer :: idirect , iemiss
!
  real(dp) :: dt , dt2 , dtbdys
  real(dp) :: dx , dx2 , dx4 , dx8 , dx16 , dxsq
  real(dp) :: c200 , rdxsq , dtsrf , dtabem , dtrad , cpldt
  real(dp) :: xkhmax , xkhz

  integer :: iboudy , ichem , ipgf , ipptls , cplexvars , cplinterp

  logical :: ifrest , rfstrt , doing_restart , cplbdysmooth
  logical :: lband

  integer :: kchi , kclo , kcmd , cpldbglevel
!
  real(dp) :: akht1 , akht2

  real(dp) , pointer , dimension(:) :: dtau
  real(dp) , pointer , dimension(:) :: hsigma , dsigma , qcon
  real(dp) , pointer , dimension(:) :: sigma
  real(dp) , pointer , dimension(:,:) :: twt

  character(len=8) :: scenario

  integer, private  :: ierr 
  real(dp) , private :: total_allocation_size

  data total_allocation_size /d_zero/
  data doing_restart /.false./

  contains

  subroutine allocate_mod_runparams
    implicit none
    call getmem1d(hsigma,1,kz,'mod_runparams:hsigma')
    call getmem1d(dsigma,1,kz,'mod_runparams:dsigma')
    call getmem1d(qcon,1,kz,'mod_runparams:qcon')
    call getmem1d(sigma,1,kzp1,'mod_runparams:sigma')
    call getmem2d(twt,1,kz,1,2,'mod_runparams:twt')
    call getmem1d(dtau,1,nsplit,'mod_runparams:nsplit')
  end subroutine allocate_mod_runparams
!
  logical function iswater(a)
    real(dp) , intent(in) :: a
    iswater = .false.
    if (a > 13.5D0 .and. a < 15.5D0) iswater = .true.
  end function
!
  logical function isocean(a)
    real(dp) , intent(in) :: a
    isocean = .false.
    if (a > 14.5D0 .and. a < 15.5D0) isocean = .true.
  end function
!
  logical function islake(a)
    real(dp) , intent(in) :: a
    islake = .false.
    if (a > 13.5D0 .and. a < 14.5D0) islake = .true.
  end function
!
  subroutine cdumlbcs(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (iboudy)
      case(0)
       write (cdum,'(a)') 'Fixed'
      case(1)
       write (cdum,'(a)') 'Relaxation, linear technique'
      case(2)
       write (cdum,'(a)') 'Time-dependent'
      case(3)
       write (cdum,'(a)') 'Time and inflow/outflow dependent'
      case(4)
       write (cdum,'(a)') 'Sponge (Perkey & Kreitzberg, MWR 1976)'
      case(5)
       write (cdum,'(a)') 'Relaxation, exponential technique'
      case default 
       write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumlbcs

  subroutine cdumcums(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (icup)
      case(1)
       write (cdum,'(a)') 'Kuo'
      case(2)
       write (cdum,'(a)') 'Grell'
      case(3)
       write (cdum,'(a)') 'Betts-Miller (1986)'
      case(4)
       write (cdum,'(a)') 'Emanuel (1991)'
      case(5)
        write (cdum,'(a)') 'Tiedtke (1986)'
      case(98)
        write (cdum,'(a)') 'Grell over ocean, Emanuel (1991) over land'
      case(99)
        write (cdum,'(a)') 'Emanuel (1991) over ocean, Grell over land'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumcums

  subroutine cdumcumcl(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (igcc)
      case(1)
        write (cdum,'(a)') 'Arakawa & Schubert (1974)'
      case(2)
        write (cdum,'(a)') 'Fritsch & Chappell (1980)'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumcumcl

  subroutine cdumpbl(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (ibltyp)
      case(0)
        write (cdum,'(a)') 'Frictionless'
      case(1)
        write (cdum,'(a)') 'Holtslag PBL (Holtslag, 1990)'
      case(2)
        write (cdum,'(a)') 'UW PBL (Bretherton and McCaa, 2004)'
      case(99)
        write (cdum,'(a)') 'Holtslag PBL, with UW in diag. mode'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumpbl

  subroutine cdummoist(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (ipptls)
      case(1)
        write (cdum,'(a)') 'Explicit moisture (SUBEX; Pal et al 2000)'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdummoist

  subroutine cdumocnflx(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (iocnflx)
      case(1)
        write (cdum,'(a)') 'Use BATS1e Monin-Obukhov'
      case(2)
        write (cdum,'(a)') 'Zeng et al (1998)'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumocnflx

  subroutine cdumpgfs(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (ipgf)
      case(0)
        write (cdum,'(a)') 'Use full fields'
      case(1)
        write (cdum,'(a)') 'Hydrostatic deduction with perturbation temperature'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumpgfs

  subroutine cdumemiss(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (iemiss)
      case(0)
        write (cdum,'(a)') 'No'
      case(1)
        write (cdum,'(a)') 'Yes'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumemiss

  subroutine cdumlakes(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (lakemod)
      case(0)
        write (cdum,'(a)') 'No'
      case(1)
        write (cdum,'(a)') 'Yes'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumlakes

  subroutine cdumchems(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (ichem)
      case(0)
        write (cdum,'(a)') 'Not active'
      case(1)
        write (cdum,'(a)') 'Active'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumchems

  subroutine cdumdcsst(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (idcsst)
      case(0)
        write (cdum,'(a)') 'Not active'
      case(1)
        write (cdum,'(a)') 'Active'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumdcsst

  subroutine cdumseaice(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (iseaice)
      case(0)
        write (cdum,'(a)') 'Not active'
      case(1)
        write (cdum,'(a)') 'Active'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumseaice

  subroutine cdumdesseas(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (idesseas)
      case(0)
        write (cdum,'(a)') 'Not active'
      case(1)
        write (cdum,'(a)') 'Active'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumdesseas

end module mod_runparams
