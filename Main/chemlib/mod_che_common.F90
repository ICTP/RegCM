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

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_mppparam
  use mod_runparams
  use mod_memutil
  use mod_stdio
  use mod_mpmessage
  use mod_che_param
  use mod_che_species
  use mod_che_indices
  use mod_cbmz_global , only : xr , xrin , xrout , c
  use mod_cbmz_parameters , only : nfix

  implicit none

  public

  real(rkx) :: cfdout

  integer(ik4) , parameter :: sbin = 2
  integer(ik4) , parameter :: cbin = 12
  !
  ! Only one cover type per grid cell for now
  !
  integer(ik4) , parameter :: luc = 1

  ! tracer indices :
  type tracer
    integer(ik4) , pointer , dimension(:) :: index
    integer(ik4) , pointer , dimension(:) :: indcbmz
    integer(ik4) , pointer , dimension(:) :: indchbdy
    real(rkx)    , pointer , dimension(:) :: mw
  end type tracer

  type(tracer) trac

  ! tracer variables

  real(rkx) , pointer , dimension(:,:,:,:) :: chi
  real(rkx) , pointer , dimension(:,:,:,:) :: chiten , chemten

  real(rkx) , pointer , dimension(:,:,:) :: chemsrc, tmpsrc
  real(rkx) , pointer , dimension(:,:,:) :: chemsrcbb , chemsrcan
  real(rkx) , pointer , dimension(:,:,:,:) :: chia , chib , chemt
  real(rkx) , pointer , dimension(:,:,:) :: dtrace , wdwout , wdrout , ddv_out
  real(rkx) , pointer , dimension(:,:,:) :: drydepv , cfmz

  real(rkx) , pointer , dimension(:,:,:,:) :: chemall , jphoto
  integer(ik4) , pointer , dimension(:,:) :: kcumtop , kcumbot , cveg2d

  real(rkx) , pointer , dimension(:,:)   :: chtrsize
  real(rkx) , pointer , dimension(:)     :: chtrsol

  real(rkx), pointer , dimension(:,:,:)  :: chifxuw , cccn

  integer(ik4) , pointer , dimension(:) :: isslt , icarb , idust
  integer(ik4) , pointer , dimension(:,:) :: imine
  integer(ik4) , parameter :: nphoto = 56

  real(rkx) , pointer , dimension(:,:,:) :: convcldfra , cemtrac , remdrd

  ! diagnostic
  real(rkx) , pointer , dimension(:,:,:,:) :: washout , rainout
  real(rkx) , pointer , dimension(:,:,:,:) :: chemdiag , cadvhdiag , &
          cadvvdiag , cdifhdiag , cconvdiag , cbdydiag , ctbldiag ,  &
          cseddpdiag
  real(rkx) , pointer , dimension(:,:,:,:) :: cemisdiag

  !*****************************************************************************
  ! INTERFACE VARIABLES  for chemistry / regcm
  !   the pointer targets are defined in mod_che_interface
  !*****************************************************************************

  real(rkx) , pointer , dimension(:,:,:,:) ::chib3d
  real(rkx) , pointer , dimension(:,:,:,:) :: cqxb3d
  real(rkx) , pointer , dimension(:,:) :: bndp0 , bndp1
  real(rkx) , pointer , dimension(:,:,:) :: tvirt0 , tvirt1
  real(rkx) , pointer , dimension(:,:,:) :: ctb3d , cubx3d , cvbx3d ,  &
         crhob3d , cpb3d , cpf3d , cfcc , cza , czq , cdzq , ccldfra , &
         crembc , cremrat ,  cconvpr , crhb3d , cdrydepflx , cwetdepflx
  real(rkx) , pointer , dimension(:,:) :: cpsb , ctg , ctga , clndcat , cht , &
         cssw2da , cvegfrac , cxlai2d , csol2d , csdeltk2d , csdelqk2d , &
         custar , csfracv2d , csfracb2d , csfracs2d , cxlat , crainc ,  &
         cps2d , cps0 , cptrop, cw10m, cdlat,cdlon
  real(rkx) , pointer , dimension(:,:) :: psbb0 , psbb1 , crho2d
  real(rkx) , pointer , dimension(:,:) :: czen
  real(rkx) , pointer , dimension(:,:,:,:) :: ctaucld
#if defined CLM45
  ! Tracer mask that uses MEGAN indices
  integer(ik4) , pointer , dimension(:) :: bvoc_trmask
  real(rkx) , pointer , dimension(:,:,:) :: cvoc_em_clm
  real(rkx) , pointer , dimension(:,:,:) :: cdustflx_clm
  real(rkx) , pointer , dimension(:,:,:) :: csw_vol
  real(rkx) , pointer , dimension(:,:,:) :: ctsoi
#endif

  contains

  subroutine allocate_mod_che_common
    implicit none

    if ( ichem /= 1 ) return

    call getmem1d(trac%index,1,ntr,'mod_che_common:trac%index')
    call getmem1d(trac%indcbmz,1,ntr,'mod_che_common:trac%indcbmz')
    call getmem1d(trac%mw,1,ntr,'mod_che_common:trac%mw')
    call getmem1d(trac%indchbdy,1,ntr,'mod_che_common:trac%indchbdy')

    if ( igaschem == 1 .and. ichsolver > 0 ) then
      call getmem4d(chemten,jci1,jci2, &
                    ici1,ici2,1,kz,1,ntr,'che_common:chemten')
    end if

    call getmem3d(chemsrc,jce1,jce2,ice1,ice2, &
                  1,ntr,'mod_che_common:chemsrc')
    call getmem3d(chemsrcbb,jce1,jce2,ice1,ice2, &
                  1,ntr,'mod_che_common:chemsrcbb')
    call getmem3d(chemsrcan,jce1,jce2,ice1,ice2, &
                  1,ntr,'mod_che_common:chemsrcan')
    call getmem3d(tmpsrc,jce1,jce2,ice1,ice2, &
                  1,ntr,'mod_che_common:tmpsrc')
    call getmem3d(chifxuw,jci1,jci2,ici1,ici2, &
                  1,ntr,'mod_che_common:chifxuw')
    call getmem3d(convcldfra,jci1,jci2,ici1,ici2, &
                  1,kz,'mod_che_common:convcldfra')
    call getmem4d(rainout,jci1,jci2,ici1,ici2,1,kz,1,ntr, &
                  'che_common:rainout')
    call getmem4d(washout,jci1,jci2,ici1,ici2,1,kz,1,ntr, &
                  'che_common:washout')
    call getmem3d(remdrd,jci1,jci2,ici1,ici2,1,ntr,'che_common:remdrd')
    call getmem1d(chtrsol,1,ntr,'mod_che_common:chtrsol')
    call getmem1d(idust,1,nbin,'mod_che_common:idust')
    call getmem1d(isslt,1,sbin,'mod_che_common:isslt')
    call getmem1d(icarb,1,cbin,'mod_che_common:icarb')
    call getmem2d(chtrsize,1,nbin,1,2,'mod_che_common:chtrsize')

    if ( nmine > 0 ) then
      call getmem2d(imine,1,nbin,1,nmine,'mod_che_common:imine')
    end if
    if ( igaschem == 1 .and. ichsolver > 0 ) then
      call getmem4d(chemall,jci1,jci2,ici1,ici2, &
                    1,kz,1,totsp,'mod_che_common:chemall')
      call getmem4d(jphoto,jci1,jci2,ici1,ici2,1,kz, &
                    1,nphoto,'che_common:jphoto')
      call getmem1d(xr,1,totsp,'che_common:xr')
      call getmem1d(xrin,1,totsp,'che_common:xrin')
      call getmem1d(xrout,1,totsp,'che_common:xrout')
      call getmem1d(c,1,totsp+nfix,'che_common:c')
    end if

    call getmem3d(dtrace,jci1,jci2,ici1,ici2,1,ntr,'che_common:dtrace')
    call getmem3d(wdrout,jci1,jci2,ici1,ici2,1,ntr,'che_common:wdrout')
    call getmem3d(wdwout,jci1,jci2,ici1,ici2,1,ntr,'che_common:wdwout')

    call getmem3d(cemtrac,jci1,jci2,ici1,ici2,1,ntr,'che_common:cemtrac')
    call getmem3d(drydepv,jci1,jci2,ici1,ici2,1,ntr,'che_common:drydepv')
    call getmem3d(ddv_out,jci1,jci2,ici1,ici2,1,ntr,'che_common:ddv_out')

    if ( ichdiag > 0 ) then
      call getmem4d(chemdiag,jci1,jci2, &
                    ici1,ici2,1,kz,1,ntr,'che_common:chemdiag')
      call getmem4d(cadvhdiag,jci1,jci2, &
                    ici1,ici2,1,kz,1,ntr,'che_common:cadvhdiag')
      call getmem4d(cadvvdiag,jci1,jci2, &
                    ici1,ici2,1,kz,1,ntr,'che_common:cadvvdiag')
      call getmem4d(cdifhdiag,jci1,jci2, &
                    ici1,ici2,1,kz,1,ntr,'che_common:cdifhdiag')
      call getmem4d(cconvdiag,jci1,jci2, &
                    ici1,ici2,1,kz,1,ntr,'che_common:cconvdiag')
      call getmem4d(ctbldiag,jci1,jci2, &
                    ici1,ici2,1,kz,1,ntr,'che_common:ctbldiag')
      call getmem4d(cbdydiag,jci1,jci2, &
                    ici1,ici2,1,kz,1,ntr,'che_common:cbdydiag')
      call getmem4d(cseddpdiag,jci1,jci2, &
                    ici1,ici2,1,kz,1,ntr,'che_common:cseddpdiag')
      call getmem4d(cemisdiag,jci1,jci2, &
                    ici1,ici2,1,kz,1,ntr,'che_common:cemisdiag')
    end if
#if defined CLM45 || (defined CLM && defined VOC)
    call getmem1d(bvoc_trmask,1,ntr,'mod_che_common:bvoc_trmask')
#endif
  end subroutine allocate_mod_che_common

  subroutine chem_config
    implicit none
    ! Define here the possible types of simulation and fix the dimension
    ! of relevant tracer dimension and parameters
    ntr = 0
    nbin = 0
    nmine = 0
    igaschem = 0
    iaerosol = 0
    iisoropia = 0
    ioxclim = 0
    if ( chemsimtype(1:4) == 'DUST' ) then
      nbin = 4
      ntr = nbin
      allocate(chtrname(nbin))
      chtrname(1:ntr)(1:6) = ['DUST01','DUST02','DUST03','DUST04']
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'DUST simulation'
    else if ( chemsimtype(1:4) == 'DU12' ) then
      nbin = 12
      ntr = nbin
      allocate(chtrname(nbin))
      chtrname(1:ntr)(1:6) = ['DUST01','DUST02','DUST03','DUST04',&
                              'DUST05','DUST06','DUST07','DUST08',&
                              'DUST09','DUST10','DUST11','DUST12' ]
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'DUST 12 bins simulation'
    else if (chemsimtype(1:4) == 'MINE') then
      nmine = 13
      nbin = 4
      ntr =  nbin * (nmine + 1)
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = ['DUST01', 'DUST02', 'DUST03', 'DUST04' ,&
                              'IRON01', 'IRON02', 'IRON03', 'IRON04' ,&
                              'HEMT01', 'HEMT02', 'HEMT03', 'HEMT04' ,&
                              'CALC01', 'CALC02', 'CALC03', 'CALC04' ,&
                              'GOTH01', 'GOTH02', 'GOTH03', 'GOTH04' ,&
                              'CHLR01', 'CHLR02', 'CHLR03', 'CHLR04' ,&
                              'FLDS01', 'FLDS02', 'FLDS03', 'FLDS04' ,&
                              'QRTZ01', 'QRTZ02', 'QRTZ03', 'QRTZ04' ,&
                              'SMEC01', 'SMEC02', 'SMEC03', 'SMEC04' ,&
                              'VRMC01', 'VRMC02', 'VRMC03', 'VRMC04' ,&
                              'GYPS01', 'GYPS02', 'GYPS03', 'GYPS04' ,&
                              'MICA01', 'MICA02', 'MICA03', 'MICA04' ,&
                              'KALO01', 'KALO02', 'KALO03', 'KALO04' ,&
                              'ILIT01', 'ILIT02', 'ILIT03', 'ILIT04' ]
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'MINERALS 70 bins simulation'
    else if ( chemsimtype(1:4) == 'SSLT' ) then
      ntr = sbin
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = ['SSLT01','SSLT02']
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'SSLT simulation'
    else if ( chemsimtype(1:4) == 'DUSS' ) then
      nbin = 4
      ntr = sbin + nbin
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = ['DUST01','DUST02','DUST03','DUST04', &
                              'SSLT01','SSLT02']
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'DUSS simulation'
    else if ( chemsimtype(1:4) == 'FIRE' ) then
      ntr = 2
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = ['SM1  ','SM2  ']
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'FIRE simulation'
    else if ( chemsimtype(1:4) == 'CARB' ) then
      if ( ismoke == 1 ) then
        ntr = 6
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:6) = ['BC_HL ','BC_HB ','OC_HL ', &
                                 'OC_HB ','SM1   ', 'SM2   ']
      else
        ntr = 4
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:6) = ['BC_HL ','BC_HB ','OC_HL ','OC_HB ']
      end if
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'CARB simulation'
    else if ( chemsimtype(1:4) == 'SULF' ) then
      ntr = 2
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = ['SO2   ','SO4   ']
      iaerosol = 1
      ioxclim  = 1
      if ( myid == italk ) write(stdout,*) 'SULF simulation'
    else if ( chemsimtype(1:4) == 'SUCA' ) then
      ntr = 6
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = ['BC_HL ','BC_HB ','OC_HL ','OC_HB ', &
                              'SO2   ','SO4   ']
      iaerosol = 1
      ioxclim  = 1
      if ( myid == italk ) write(stdout,*) 'SUCA simulation'
    else if ( chemsimtype(1:4) == 'SUCE' ) then
      ntr = 10
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = ['BC_HB ','BC_HL ','BC_HL1','BC_HL2', &
                              'OC_HB ','OC_HL ','OC_HL1','OC_HL2', &
                              'SO2   ','SO4   ']
      iaerosol = 1
      ioxclim  = 1
      if ( myid == italk ) write(stdout,*) 'SUCE simulation'
    else if ( chemsimtype(1:4) == 'AERO' ) then
      nbin = 4
      iaerosol = 1
      ioxclim = 1
      if ( ismoke == 1 ) then
        ntr = 14
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:6) = ['BC_HL ','BC_HB ','OC_HL ','OC_HB ', &
                                'SO2   ','SO4   ','DUST01','DUST02', &
                                'DUST03','DUST04','SSLT01','SSLT02', &
                                'SM1   ','SM2   ' ]
      else
        ntr = 12
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:6) = ['BC_HL ','BC_HB ','OC_HL ','OC_HB ', &
                                'SO2   ','SO4   ','DUST01','DUST02', &
                                'DUST03','DUST04','SSLT01','SSLT02' ]
      end if
      if ( myid == italk ) write(stdout,*) 'AERO simulation'
    else if ( chemsimtype(1:4) == 'AERE' ) then
      nbin = 4
      iaerosol = 1
      ioxclim = 1
      ntr = 16
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = ['BC_HB ','BC_HL ','BC_HL1','BC_HL2', &
                              'OC_HB ','OC_HL ','OC_HL1','OC_HL2', &
                              'SO2   ','SO4   ','DUST01','DUST02', &
                              'DUST03','DUST04','SSLT01','SSLT02' ]
      if ( myid == italk ) write(stdout,*) 'AERE simulation'
    else if ( chemsimtype(1:4) == 'DCCB' ) then
      nbin = 4
      if ( ismoke == 1 ) then
        ntr = 52
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:6) = ['NO    ','NO2   ','N2O5  ','HNO2  ',&
                                'HNO3  ','HNO4  ','O3    ','H2O2  ',&
                                'CO    ','SO2   ','DMS   ','H2SO4 ',&
                                'CH4   ','C2H6  ','PAR   ','CH3OH ',&
                                'HCHO  ','ALD2  ','AONE  ','ETH   ',&
                                'OLET  ','OLEI  ','TOL   ','XYL   ',&
                                'ISOP  ','ONIT  ','PAN   ','HCOOH ',&
                                'RCOOH ','CH3OOH','ETHOOH','ROOH  ',&
                                'MGLY  ','ISOPRD','ISOPN ','OPEN  ',&
                                'CRES  ','NH3   ','DUST01','DUST02',&
                                'DUST03','DUST04','BC_HL ','BC_HB ',&
                                'OC_HL ','OC_HB ','SSLT01','SSLT02',&
                                'ANO3  ','ANH4  ','SM1   ','SM2   ' ]
      else
        ntr = 50
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:6) = ['NO    ','NO2   ','N2O5  ','HNO2  ',&
                                'HNO3  ','HNO4  ','O3    ','H2O2  ',&
                                'CO    ','SO2   ','DMS   ','H2SO4 ',&
                                'CH4   ','C2H6  ','PAR   ','CH3OH ',&
                                'HCHO  ','ALD2  ','AONE  ','ETH   ',&
                                'OLET  ','OLEI  ','TOL   ','XYL   ',&
                                'ISOP  ','ONIT  ','PAN   ','HCOOH ',&
                                'RCOOH ','CH3OOH','ETHOOH','ROOH  ',&
                                'MGLY  ','ISOPRD','ISOPN ','OPEN  ',&
                                'CRES  ','NH3   ','DUST01','DUST02',&
                                'DUST03','DUST04','BC_HL ','BC_HB ',&
                                'OC_HL ','OC_HB ','SSLT01','SSLT02',&
                                'ANO3  ','ANH4  ']
      end if
      iaerosol = 1
      igaschem = 1
      iisoropia = 1
      if ( myid == italk ) write(stdout,*) 'DCCB simulation'
    else if ( chemsimtype(1:4) == 'CBMZ' ) then
      ! This does not include any aerosol(NH3) or monoterpens(APIN, LIMO)
      ntr = 37
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = ['NO    ','NO2   ','N2O5  ','HNO2  ',&
                              'HNO3  ','HNO4  ','O3    ','H2O2  ',&
                              'CO    ','SO2   ','DMS   ','H2SO4 ',&
                              'CH4   ','C2H6  ','PAR   ','CH3OH ',&
                              'HCHO  ','ALD2  ','AONE  ','ETH   ',&
                              'OLET  ','OLEI  ','TOL   ','XYL   ',&
                              'ISOP  ','ONIT  ','PAN   ','HCOOH ',&
                              'RCOOH ','CH3OOH','ETHOOH','ROOH  ',&
                              'MGLY  ','ISOPRD','ISOPN ','OPEN  ',&
                              'CRES  ']
      igaschem = 1
      if ( myid == italk ) write(stdout,*) 'CBMZ simulation'
    else if ( chemsimtype(1:6) == 'POLLEN' ) then
      ntr = 1
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = ['POLLEN' ]
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'POLLEN simulation'
    else
      if ( myid == italk ) then
        write (stderr,*) 'Not a valid chemtype simulation : STOP !'
        write (stderr,*) 'Valid simulations are : ' , &
           'DUST DU12 SSLT DUSS CARB SULF SUCA SUCE AERO CBMZ DCCB POLLEN'
      end if
      call fatal(__FILE__,__LINE__,'INVALID CHEM CONFIGURATION')
    end if

  end subroutine chem_config

end module mod_che_common
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
