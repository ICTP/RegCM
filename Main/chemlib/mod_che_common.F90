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
  use mod_mpmessage
  use mod_che_param
  use mod_che_species
  use mod_che_indices
  use mod_cbmz_global , only : xr , xrin , xrout , c
  use mod_cbmz_parameters , only : nfix

  implicit none

  public

  integer(ik4) , parameter :: sbin = 2
  integer(ik4) , parameter :: maxntr = 50

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
  real(rkx) , pointer , dimension(:,:,:,:) :: chic , chiten , chiten0 , chemten

  real(rkx) , pointer , dimension(:,:,:) :: chemsrc, tmpsrc
  real(rkx) , pointer , dimension(:,:,:,:) :: chia , chib
  real(rkx) , pointer , dimension(:,:,:) :: dtrace , wdwout , &
                                           wdrout , wxaq , wxsg , ddv_out
  real(rkx) , pointer , dimension(:,:,:) :: drydepv

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
  real(rkx) , pointer , dimension(:,:,:,:) :: washout , rainout , &
                                             rxsaq1 , rxsaq2 , rxsg
  real(rkx) , pointer , dimension(:,:,:,:) :: chemdiag , cadvhdiag , &
          cadvvdiag , cdifhdiag , cconvdiag , cbdydiag , ctbldiag ,  &
          cseddpdiag
  real(rkx) , pointer , dimension(:,:,:) :: cemisdiag


!*****************************************************************************
! INTERFACE VARIABLES  for chemistry / regcm
!   the pointer targets are defined in mod_che_interface
!*****************************************************************************

  real(rkx) , pointer , dimension(:,:,:,:) ::chib3d
  real(rkx) , pointer , dimension(:,:,:,:) :: cqxb3d
  real(rkx) , pointer , dimension(:,:,:) :: bndt0 , bndt1 , bndq0 , bndq1
  real(rkx) , pointer , dimension(:,:) :: bndp0 , bndp1
  real(rkx) , pointer , dimension(:,:) :: sp0 , sp1
  real(rkx) , pointer , dimension(:,:,:) :: ctb3d , cubx3d , cvbx3d ,  &
         crhob3d , cpb3d , cpf3d , cfcc , cza , czq , cdzq , ccldfra , &
         crembc , cremrat ,  cconvpr , crhb3d , cdrydepflx , cwetdepflx , tvirt
  real(rkx) , pointer , dimension(:,:) :: cpsb , ctg , ctga , clndcat , cht , &
         cssw2da , cvegfrac , cxlai2d , csol2d , csdeltk2d , csdelqk2d , &
         cuvdrag , csfracv2d , csfracb2d , csfracs2d , cxlat , crainc ,  &
         cps2d , cps0 , cptrop
  real(rkx) , pointer , dimension(:,:) :: psbb0 , psbb1
  real(rkx) , pointer , dimension(:,:) :: czen
  real(rkx) , pointer , dimension(:,:,:,:) :: ctaucld
#if defined CLM45
  ! Tracer mask that uses MEGAN indices
  integer(ik4) , pointer , dimension(:) :: bvoc_trmask
  real(rkx) , pointer , dimension(:,:,:) :: cvoc_em_clm
  real(rkx) , pointer , dimension(:,:,:) :: cdustflx_clm
#endif
#if (defined CLM && defined VOC)
  ! Tracer mask that uses MEGAN indices
  integer(ik4) , pointer , dimension(:) :: bvoc_trmask
  real(rkx) , pointer , dimension(:,:) :: cvoc_em0
  real(rkx) , pointer , dimension(:,:) :: cvoc_em1
  real(rkx) , pointer , dimension(:,:) :: cvoc_em2
#endif
#if defined CLM
  real(rkx) , pointer , dimension(:,:,:) :: cdep_vels
#endif

  contains

  subroutine allocate_mod_che_common(isladvec)
    implicit none
    integer(ik4) , intent(in) :: isladvec

    if ( ichem == 1 ) then

      call getmem1d(trac%index,1,ntr,'mod_che_common:trac%index')
      call getmem1d(trac%indcbmz,1,ntr,'mod_che_common:trac%indcbmz')
      call getmem1d(trac%mw,1,ntr,'mod_che_common:trac%mw')
      call getmem1d(trac%indchbdy,1,ntr,'mod_che_common:trac%indchbdy')

      call getmem4d(chia,jce1ga,jce2ga,ice1ga,ice2ga, &
                         1,kz,1,ntr,'che_common:chia')
      if ( isladvec == 1 ) then
        call getmem4d(chib,jce1sl,jce2sl,ice1sl,ice2sl, &
                           1,kz,1,ntr,'che_common:chib')
      else
        call getmem4d(chib,jce1gb,jce2gb,ice1gb,ice2gb, &
                           1,kz,1,ntr,'che_common:chib')
      end if
      call getmem4d(chi,jce1gb,jce2gb,ice1gb,ice2gb, &
                        1,kz,1,ntr,'che_common:chi')
      if ( idynamic == 2 ) then
        call getmem3d(tvirt,jce1,jce2,ice1,ice2,1,kz,'che_common:tvirt')
        call getmem2d(sp0,jce1,jce2,ice1,ice2,'che_common:sp0')
        call getmem2d(sp1,jce1,jce2,ice1,ice2,'che_common:sp1')
      end if
      call getmem4d(chic,jce1,jce2,ice1,ice2,1,kz,1,ntr,'che_common:chic')
      call getmem4d(chiten,jce1,jce2,ice1,ice2,1,kz,1,ntr,'che_common:chiten')
      call getmem4d(chemten,jce1,jce2, &
                    ice1,ice2,1,kz,1,ntr,'che_common:chemten')
      call getmem3d(chemsrc,jce1,jce2,ice1,ice2, &
                    1,ntr,'mod_che_common:chemsrc')
      call getmem3d(tmpsrc,jce1,jce2,ice1,ice2, &
                    1,ntr,'mod_che_common:tmpsrc')
      call getmem3d(chifxuw,jci1,jci2,ici1,ici2, &
                    1,ntr,'mod_che_common:chifxuw')
      call getmem3d(convcldfra,jci1,jci2,ici1,ici2, &
                    1,kz,'mod_che_common:convcldfra')
      call getmem4d(rxsg,jce1,jce2,ice1,ice2,1,kz,1,ntr, &
                    'che_common:rxsg')
      call getmem4d(rxsaq1,jce1,jce2,ice1,ice2,1,kz,1,ntr, &
                    'che_common:rxsaq1')
      call getmem4d(rxsaq2,jce1,jce2,ice1,ice2,1,kz,1,ntr, &
                    'che_common:rxsaq2')
      call getmem4d(rainout,jce1,jce2,ice1,ice2,1,kz,1,ntr, &
                    'che_common:rainout')
      call getmem4d(washout,jce1,jce2,ice1,ice2,1,kz,1,ntr, &
                    'che_common:washout')
      call getmem3d(remdrd,jce1,jce2,ice1,ice2,1,ntr,'che_common:remdrd')

      call getmem1d(chtrsol,1,ntr,'mod_che_common:chtrsol')
      call getmem1d(idust,1,nbin,'mod_che_common:idust')
      call getmem1d(isslt,1,sbin,'mod_che_common:isslt')
      call getmem1d(icarb,1,7,'mod_che_common:icarb')
      call getmem2d(chtrsize,1,nbin,1,2,'mod_che_common:chtrsize')
      call getmem2d(imine,1,nbin,1,nmine,'mod_che_common:imine')

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

      call getmem3d(dtrace,jce1,jce2,ice1,ice2,1,ntr,'che_common:dtrace')
      call getmem3d(wdrout,jce1,jce2,ice1,ice2,1,ntr,'che_common:wdrout')
      call getmem3d(wdwout,jce1,jce2,ice1,ice2,1,ntr,'che_common:wdwout')

      call getmem3d(wxsg,jce1,jce2,ice1,ice2,1,ntr,'che_common:wxsg')
      call getmem3d(wxaq,jce1,jce2,ice1,ice2,1,ntr,'che_common:wxaq')
      call getmem3d(cemtrac,jce1,jce2,ice1,ice2,1,ntr,'che_common:cemtrac')
      call getmem3d(drydepv,jce1,jce2,ice1,ice2,1,ntr,'che_common:drydepv')
      call getmem3d(ddv_out,jce1,jce2,ice1,ice2,1,ntr,'che_common:ddv_out')

      if ( ichdiag > 0 ) then
        call getmem4d(chiten0,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:chiten0')
        call getmem4d(chemdiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:chemdiag')
        call getmem4d(cadvhdiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:cadvhdiag')
        call getmem4d(cadvvdiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:cadvvdiag')
        call getmem4d(cdifhdiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:cdifhdiag')
        call getmem4d(cconvdiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:cconvdiag')
        call getmem4d(ctbldiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:ctbldiag')
        call getmem4d(cbdydiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:cbdydiag')
        call getmem4d(cseddpdiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:cseddpdiag')
        call getmem3d(cemisdiag,jce1,jce2, &
                      ice1,ice2,1,ntr,'che_common:cemisdiag')
      end if
#if defined CLM45 || (defined CLM && defined VOC)
      call getmem1d(bvoc_trmask,1,ntr,'mod_che_common:bvoc_trmask')
#endif
    end if
  end subroutine allocate_mod_che_common
!
  subroutine chem_config
    implicit none
    ! Define here the possible types of simulation and fix the dimension
    ! of relevant tracer dimension and parameters
    igaschem = 0
    iaerosol = 0
    iisoropia = 0
    if ( chemsimtype(1:4) == 'DUST' ) then
      nbin = 4
      ntr = nbin
      allocate(chtrname(nbin))
      chtrname(1:ntr)(1:6) = (/'DUST01','DUST02','DUST03','DUST04'/)
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'DUST simulation'
    else if ( chemsimtype(1:4) == 'DU12' ) then
      nbin = 12
      ntr = nbin
      allocate(chtrname(nbin))
      chtrname(1:ntr)(1:6) = (/'DUST01','DUST02','DUST03','DUST04',&
                               'DUST05','DUST06','DUST07','DUST08',&
                               'DUST09','DUST10','DUST11','DUST12' /)
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'DUST 12 bins simulation'
    else if (chemsimtype(1:4) == 'MINE') then
      nmine = 3 
      nbin = 4
      ntr =  nbin * (nmine + 1)
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'DUST01','DUST02','DUST03','DUST04',&
                               'IRON01', 'IRON02', 'IRON03', 'IRON04',&
                               'HEMT01', 'HEMT02', 'HEMT03', 'HEMT04',&
                               'CALC01', 'CALC02', 'CALC03', 'CALC04' /)
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'MINERALS 4  bins simulation'
    else if ( chemsimtype(1:4) == 'SSLT' ) then
      ntr = sbin
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'SSLT01','SSLT02'/)
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'SSLT simulation'
    else if ( chemsimtype(1:4) == 'DUSS' ) then
      nbin = 4
      ntr = sbin + nbin
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'DUST01','DUST02','DUST03','DUST04', &
                               'SSLT01','SSLT02'/)
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'DUSS simulation'
    else if ( chemsimtype(1:4) == 'FIRE' ) then
      ntr = 2
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'SM1  ','SM2  '/)
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'FIRE simulation'
    else if ( chemsimtype(1:4) == 'CARB' ) then
      if ( ismoke == 1 ) then
        ntr = 6
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:6) = (/'BC_HL ','BC_HB ','OC_HL ', &
                                 'OC_HB ','SM1   ', 'SM2   '/)
      else
        ntr = 4
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:6) = (/'BC_HL ','BC_HB ','OC_HL ','OC_HB '/)
      end if
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'CARB simulation'
    else if ( chemsimtype(1:4) == 'SULF' ) then
      ntr = 2
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'SO2   ','SO4   '/)
      iaerosol = 1
      ioxclim  = 1
      if ( myid == italk ) write(stdout,*) 'SULF simulation'
    else if ( chemsimtype(1:4) == 'SUCA' ) then
      ntr = 6
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'BC_HL ','BC_HB ','OC_HL ','OC_HB ', &
                               'SO2   ','SO4   '/)
      iaerosol = 1
      ioxclim  = 1
      if ( myid == italk ) write(stdout,*) 'SUCA simulation'
    else if ( chemsimtype(1:4) == 'AERO' ) then
      nbin = 4
      ntr = 12
      allocate(chtrname(ntr))
      iaerosol = 1
      ioxclim  = 1
      chtrname(1:ntr)(1:6) = (/'BC_HL ','BC_HB ','OC_HL ','OC_HB ', &
                               'SO2   ','SO4   ','DUST01','DUST02', &
                               'DUST03','DUST04','SSLT01','SSLT02' /)
      if ( myid == italk ) write(stdout,*) 'AERO simulation'
    else if ( chemsimtype(1:4) == 'DCCB' ) then
      nbin = 4
      ntr = 50
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'NO    ','NO2   ','N2O5  ','HNO2  ',&
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
                               'ANO3  ','ANH4  ' /)
      iaerosol = 1
      igaschem = 1
      iisoropia = 1
      if ( myid == italk ) write(stdout,*) 'DCCB simulation'
    else if ( chemsimtype(1:4) == 'CBMZ' ) then
      ! This does not include any aerosol(NH3) or monoterpens(APIN, LIMO)
      ntr = 37
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'NO    ','NO2   ','N2O5  ','HNO2  ',&
                               'HNO3  ','HNO4  ','O3    ','H2O2  ',&
                               'CO    ','SO2   ','DMS   ','H2SO4 ',&
                               'CH4   ','C2H6  ','PAR   ','CH3OH ',&
                               'HCHO  ','ALD2  ','AONE  ','ETH   ',&
                               'OLET  ','OLEI  ','TOL   ','XYL   ',&
                               'ISOP  ','ONIT  ','PAN   ','HCOOH ',&
                               'RCOOH ','CH3OOH','ETHOOH','ROOH  ',&
                               'MGLY  ','ISOPRD','ISOPN ','OPEN  ',&
                               'CRES  '/)
      igaschem = 1
      if ( myid == italk ) write(stdout,*) 'CBMZ simulation'
    else if ( chemsimtype(1:6) == 'POLLEN' ) then
      ntr = 1
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'POLLEN' /)
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'POLLEN simulation'
    else
      if ( myid == italk ) then
        write (stderr,*) 'Not a valid chemtype simulation : STOP !'
        write (stderr,*) 'Valid simulations are : ' , &
           'DUST DU12 SSLT DUSS CARB SULF SUCA AERO CBMZ DCCB POLLEN'
      end if
      call fatal(__FILE__,__LINE__,'INVALID CHEM CONFIGURATION')
    end if

  end subroutine chem_config

end module mod_che_common
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
