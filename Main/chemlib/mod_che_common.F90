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

module mod_che_common

  use mod_realkinds
  use mod_dynparam
  use mod_mppparam
  use mod_memutil
  use mod_mpmessage
  use mod_che_param
  use mod_che_species
!
  public

  logical :: lcband

  character(len=8)   :: chemsimtype 
  integer , parameter :: nbin = 4
  integer , parameter :: sbin = 2
  integer , parameter :: maxntr = 40 

  ! chemistry nameliste option
  integer , public :: ichcumtra , ichdrdepo , ichremcvc , &
                      ichremlsc , ichsursrc, ichsolver
  ! usefull flags
  integer  :: iaerosol , igaschem

  ! tracer variables

  real(dp) , pointer , dimension(:,:,:,:) :: chi
  real(dp) , pointer , dimension(:,:,:,:) :: chic , chiten, chemten
!
  real(dp) , pointer , dimension(:,:,:,:) :: chemsrc
  real(dp) , pointer , dimension(:,:,:,:) :: chia , chib
  real(dp) , pointer , dimension(:,:,:) :: srclp2
  real(dp) , pointer , dimension(:,:,:) :: ddsfc , dtrace , wdcvc , &
                                           wdlsc , wxaq , wxsg , drydepv
  real(dp) , pointer , dimension(:,:,:) :: taucld
  real(dp) , pointer , dimension(:,:,:,:) :: chemall
  integer , pointer , dimension(:,:) :: kcumtop , kcumbot , cveg2d
!
  character(len=5) , pointer , dimension(:) :: chtrname
!
  real(dp) , pointer , dimension(:,:) :: chtrdpv

  ! evap of l-s precip (see mod_precip.f90; [kg_h2o/kg_air/s)
  ! cum h2o vapor tendency for cum precip (kg_h2o/kg_air/s)
  real(dp) , pointer , dimension(:,:,:) :: chevap , checum

  real(dp) , pointer , dimension(:,:) :: chtrsize
  real(dp) , pointer , dimension(:) :: chtrsol

  real(dp), pointer , dimension(:,:,:) :: cchifxuw
!
  integer , pointer , dimension(:) :: isslt , icarb , idust
!
  real(dp) , pointer , dimension(:,:,:) :: cemtr , cemtrac , remdrd
  real(dp) , pointer , dimension(:,:,:,:) :: remcvc , remlsc , &
                                             rxsaq1 , rxsaq2 , rxsg


!*****************************************************************************
! INTERFACE VARIABLES  for chemistry / regcm 
!   the pointer targets are defined in mod_che_interface
!*****************************************************************************

  real(dp) :: ccalday , crdxsq
  real(dp) , pointer , dimension(:,:,:,:) ::chib3d
  real(dp) , pointer , dimension(:,:,:) :: ctb3d , cqvb3d , cubx3d , cvbx3d , &
         crhob3d , cqcb3d , cfcc , cza , cdzq , ccldfra , crembc , cremrat
  real(dp) , pointer , dimension(:,:) :: cpsb , ctg , clndcat , cht , cssw2da, &
         cvegfrac , csol2d , csdeltk2d , csdelqk2d , ctwt , cuvdrag
  real(dp) , pointer , dimension(:) :: hlev , cdsigma , canudg
  real(dp) , pointer , dimension(:,:) :: czen
  real(dp) :: chptop

  type cbound_area
    logical :: dotflag
    logical :: havebound
    logical , pointer , dimension(:,:) :: bsouth
    logical , pointer , dimension(:,:) :: bnorth
    logical , pointer , dimension(:,:) :: beast
    logical , pointer , dimension(:,:) :: bwest
    integer :: ns , nn , ne , nw
    integer :: nsp
    integer , pointer , dimension(:,:) :: ibnd
  end type cbound_area

  type(cbound_area) :: cba_cr

  contains

    subroutine allocate_mod_che_common(ichem)
      implicit none

      integer , intent(in) :: ichem

      if ( ichem == 1 ) lch = .true.

      if ( lch ) then

        call getmem4d(chia,jce1-ma%jbl2,jce2+ma%jbr2, &
                      ice1-ma%ibb2,ice2+ma%ibt2,1,kz,1,ntr,'che_common:chia')
        call getmem4d(chib,jce1-ma%jbl2,jce2+ma%jbr2, &
                      ice1-ma%ibb2,ice2+ma%ibt2,1,kz,1,ntr,'che_common:chib')
        call getmem4d(chi,jce1-ma%jbl1,jce2+ma%jbr1, &
                      ice1-ma%ibb1,ice2+ma%ibt1,1,kz,1,ntr,'che_common:chi')
        call getmem4d(chi,jce1,jce2,ice1,ice2,1,kz,1,ntr,'che_common:chic')
        call getmem4d(chiten,jce1,jce2,ice1,ice2,1,kz,1,ntr,'che_common:chiten')
        call getmem4d(chemten,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:chemten')

        call getmem4d(chemsrc,jce1,jce2,ice1,ice2, &
                      1,mpy,1,ntr,'mod_che_common:chemsrc')

        call getmem3d(cemtr,1,iy,1,jxp,1,ntr,'mod_che_common:cemtr')
        call getmem3d(cemtrac,1,iy,1,jxp,1,ntr,'mod_che_common:cemtrac')
        call getmem3d(remdrd,1,iy,1,jxp,1,ntr,'mod_che_common:remdrd')
        call getmem4d(remcvc,1,iy,1,kz,1,jxp,1,ntr,'mod_che_common:remcvc')
        call getmem4d(remlsc,1,iy,1,kz,1,jxp,1,ntr,'mod_che_common:remlsc')
        call getmem4d(rxsaq1,1,iy,1,kz,1,jxp,1,ntr,'mod_che_common:rxsaq1')
        call getmem4d(rxsaq2,1,iy,1,kz,1,jxp,1,ntr,'mod_che_common:rxsaq2')
        call getmem4d(rxsg,1,iy,1,kz,1,jxp,1,ntr,'mod_che_common:rxsg')

        call getmem1d(chtrsol,1,ntr,'mod_che_common:chtrsol')
        call getmem2d(chtrdpv,1,ntr,1,2,'mod_che_common:chtrdpv')
        call getmem1d(idust,1,nbin,'mod_che_common:idust')
        call getmem1d(isslt,1,sbin,'mod_che_common:isslt')
        call getmem1d(icarb,1,5,'mod_che_common:icarb')
        call getmem2d(chtrsize,1,nbin,1,2,'mod_che_common:chtrsize')

        call getmem4d(chemall,1,jxp,1,iy,1,kz,1,totsp,'mod_che_common:chemall')
     
        call getmem3d(taucld,1,iy,1,kz,1,jxp,'mod_che_common:taucld')
        call getmem3d(srclp2,1,iy,1,jxp,1,ntr,'mod_che_common:srclp2')
        call getmem3d(ddsfc,1,iy,1,jxp,1,ntr,'mod_che_common:ddsfc')
        call getmem3d(dtrace,1,iy,1,jxp,1,ntr,'mod_che_common:dtrace')
        call getmem3d(drydepv,1,iy,1,jxp,1,ntr,'mod_che_common:drydepv')
        call getmem3d(wdcvc,1,iy,1,jxp,1,ntr,'mod_che_common:wdcvc')
        call getmem3d(wdlsc,1,iy,1,jxp,1,ntr,'mod_che_common:wdlsc')
        call getmem3d(wxaq,1,iy,1,jxp,1,ntr,'mod_che_common:wxaq')
        call getmem3d(wxsg,1,iy,1,jxp,1,ntr,'mod_che_common:wxsg')
        call getmem3d(chevap,1,jxp,1,iy,1,kz,'mod_che_common:chevap')
        call getmem3d(checum,1,jxp,1,iy,1,kz,'mod_che_common:checum')

        call getmem3d(cchifxuw,1,jxp,1,iy,1,ntr,'mod_che_common:cchifxuw')

      end if
    end subroutine allocate_mod_che_common
!

    subroutine chem_config

      implicit none
      ! Define here the possible types of simulation and fix the dimension
      ! of relevant tracer dimension and parameters 
      igaschem = 0
      iaerosol = 0
      if ( chemsimtype(1:4) == 'DUST' ) then
        ntr = nbin
        allocate(chtrname(nbin))
        chtrname(1:ntr)(1:5) = (/'DUST1','DUST2','DUST3','DUST4'/)
        iaerosol = 1
        write (aline,*) 'DUST simulation , used tracers: ', chtrname(:)
        call say
      else if ( chemsimtype(1:4) == 'SSLT' ) then 
        ntr = sbin
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:5) = (/'SSLT1','SSLT2'/)
        iaerosol = 1
        write (aline,*) 'SSLT simulation , used tracers: ', chtrname(:)
        call say
      else if ( chemsimtype(1:4) == 'CARB' ) then 
        ntr = 4
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:5) = (/'BC_HL','BC_HB','OC_HB','OC_HL'/)
        iaerosol = 1
        write (aline,*) 'CARB simulation , used tracers: ', chtrname(:)
        call say
      else if ( chemsimtype(1:4) == 'SULF' ) then 
        ntr = 2
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:5) = (/'SO2  ','SO4  '/)
        iaerosol = 1
        write (aline,*) 'SULF simulation , used tracers: ', chtrname(:)
        call say
      else if ( chemsimtype(1:4) == 'SUCA' ) then 
        ntr = 6
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:5) = (/'SO2  ','SO4  ','BC_HL', &
                                 'BC_HB','OC_HB','OC_HL' /)
        iaerosol = 1
        write (aline,*) 'SUCA simulation , used tracers: ', chtrname(:)
        call say
      else if ( chemsimtype(1:4) == 'AERO' ) then 
        ntr = 12 
        allocate(chtrname(ntr))
        iaerosol = 1
        chtrname(1:ntr)(1:5) = (/'SO2  ','SO4  ','BC_HL','BC_HB', &
                                 'OC_HB','OC_HL','DUST1','DUST2', &
                                 'DUST3','DUST4','SSLT1','SSLT2' /)
        write (aline,*) 'AERO simulation , used tracers: ', chtrname(:)
        call say
      else if ( chemsimtype(1:4) == 'CBMZ' ) then 
        ntr = 25
        allocate(chtrname(ntr))      
        chtrname(1:ntr)(1:5) = (/'SO2  ','SO4  ','DMS  ','O3   ', &
                                 'NO2  ','NO   ','CO   ','H2O2 ', &
                                 'HNO3 ','N2O5 ','HCHO ','ALD2 ', &
                                 'ISOP ','C2H6 ','PAR  ','ACET ', &
                                 'MOH  ','OLT  ','OLI  ','TOLUE', &
                                 'XYL  ','ETHE ','PAN  ','CH4  ', &
                                 'NH3  '  /)
        igaschem = 1
        write (aline,*) 'CBMZ gas-phase + sulfate simulation , ', &
                        'used tracers: ', chtrname(:)
        call say
      else 
        write (aline,*) 'Not a valid chemtype simulation : STOP !'
        call say
        call fatal(__FILE__,__LINE__,'INVALID CHEM CONFIGURATION')
      end if
    end subroutine chem_config

end module mod_che_common
