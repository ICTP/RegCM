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
!
  public

  integer(ik4) , parameter :: nbin = 4
  integer(ik4) , parameter :: sbin = 2
  integer(ik4) , parameter :: maxntr = 50 

  ! tracer variables

  real(rk8) , pointer , dimension(:,:,:,:) :: chi
  real(rk8) , pointer , dimension(:,:,:,:) :: chic , chiten , chiten0 , chemten 
!
  real(rk8) , pointer , dimension(:,:,:) :: chemsrc, tmpsrc
  real(rk8) , pointer , dimension(:,:,:,:) :: chia , chib
  real(rk8) , pointer , dimension(:,:,:) :: srclp2
  real(rk8) , pointer , dimension(:,:,:) :: dtrace , wdcvc , &
                                           wdlsc , wxaq , wxsg , ddv_out
  real(rk8) , pointer , dimension(:,:,:) :: drydepv

  real(rk8) , pointer , dimension(:,:,:,:) :: chemall , jphoto
  integer(ik4) , pointer , dimension(:,:) :: kcumtop , kcumbot , cveg2d
!
  real(rk8) , pointer , dimension(:,:)   :: chtrsize
  real(rk8) , pointer , dimension(:)     :: chtrsol

  real(rk8), pointer , dimension(:,:,:)  :: cchifxuw
!
  integer(ik4) , pointer , dimension(:) :: isslt , icarb , idust
  integer(ik4) , parameter :: nphoto = 56
!
  real(rk8) , pointer , dimension(:,:,:) :: convcldfra , cemtrac , remdrd

!diagnostic
  real(rk8) , pointer , dimension(:,:,:,:) :: remcvc , remlsc , &
                                             rxsaq1 , rxsaq2 , rxsg
  real(rk8) , pointer , dimension(:,:,:,:) :: chemdiag , cadvhdiag , &
          cadvvdiag , cdifhdiag , cconvdiag , cbdydiag , ctbldiag ,  &
          cseddpdiag , cemisdiag 


!*****************************************************************************
! INTERFACE VARIABLES  for chemistry / regcm 
!   the pointer targets are defined in mod_che_interface
!*****************************************************************************

  real(rk8) , pointer , dimension(:,:,:,:) ::chib3d
  real(rk8) , pointer , dimension(:,:,:,:) :: cqxb3d
  real(rk8) , pointer , dimension(:,:,:) :: ctb3d , cubx3d , cvbx3d , &
         crhob3d , cfcc , cza , cdzq , ccldfra , crembc , cremrat ,  &
         cconvpr , crhb3d
  real(rk8) , pointer , dimension(:,:) :: cpsb , ctg , clndcat , cht , &
         cssw2da , cvegfrac , csol2d , csdeltk2d , csdelqk2d ,         &
         cuvdrag , csfracv2d , csfracb2d , csfracs2d , cxlat , crainc
  real(rk8) , pointer , dimension(:,:) :: psbb0 , psbb1
  real(rk8) , pointer , dimension(:,:) :: czen
  real(rk8) , pointer , dimension(:,:,:,:) :: ctaucld
#if (defined CLM)
#if (defined VOC)
  ! Tracer mask that uses MEGAN indices
  integer(ik4) , pointer , dimension(:) :: bvoc_trmask
  real(rk8) , pointer , dimension(:,:) :: cvoc_em
#endif
  real(rk8) , pointer , dimension(:,:,:) :: cdep_vels
#endif

  contains

    subroutine allocate_mod_che_common(isladvec)
      implicit none
      integer(ik4) , intent(in) :: isladvec

      if ( ichem == 1 ) then
        call getmem4d(chia,jce1-ma%jbl2,jce2+ma%jbr2, &
                      ice1-ma%ibb2,ice2+ma%ibt2,1,kz,1,ntr,'che_common:chia')
        if ( isladvec == 1 ) then
          call getmem4d(chib,jce1-ma%jbl4,jce2+ma%jbr4, &
                        ice1-ma%ibb4,ice2+ma%ibt4,1,kz,1,ntr,'che_common:chib')
          call getmem4d(chi,jce1-ma%jbl4,jce2+ma%jbr4, &
                        ice1-ma%ibb4,ice2+ma%ibt4,1,kz,1,ntr,'che_common:chi')
        else
          call getmem4d(chib,jce1-ma%jbl2,jce2+ma%jbr2, &
                        ice1-ma%ibb2,ice2+ma%ibt2,1,kz,1,ntr,'che_common:chib')
          call getmem4d(chi,jce1-ma%jbl1,jce2+ma%jbr1, &
                        ice1-ma%ibb1,ice2+ma%ibt1,1,kz,1,ntr,'che_common:chi')
        end if
        call getmem4d(chic,jce1,jce2,ice1,ice2,1,kz,1,ntr,'che_common:chic')
        call getmem4d(chiten,jce1,jce2,ice1,ice2,1,kz,1,ntr,'che_common:chiten')
        call getmem4d(chemten,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:chemten')

        call getmem3d(chemsrc,jce1,jce2,ice1,ice2, &
                      1,ntr,'mod_che_common:chemsrc')

        call getmem3d(tmpsrc,jce1,jce2,ice1,ice2, &
                      1,ntr,'mod_che_common:tmpsrc')

        call getmem3d(cchifxuw,jci1,jci2,ici1,ici2, &
                      1,ntr,'mod_che_common:cchifxuw')
       
        call getmem3d(convcldfra,jci1,jci2,ici1,ici2, &
                      1,kz,'mod_che_common:convcldfra')
        
        call getmem4d(rxsg,jce1,jce2,ice1,ice2,1,kz,1,ntr, &
                      'che_common:rxsg')
        call getmem4d(rxsaq1,jce1,jce2,ice1,ice2,1,kz,1,ntr, &
                      'che_common:rxsaq1')
        call getmem4d(rxsaq2,jce1,jce2,ice1,ice2,1,kz,1,ntr, &
                      'che_common:rxsaq2')
        call getmem4d(remlsc,jce1,jce2,ice1,ice2,1,kz,1,ntr, &
                      'che_common:remlsc')
        call getmem4d(remcvc,jce1,jce2,ice1,ice2,1,kz,1,ntr, &
                      'che_common:remcvc')
        call getmem3d(remdrd,jce1,jce2,ice1,ice2,1,ntr,'che_common:remdrd')

        call getmem1d(chtrsol,1,ntr,'mod_che_common:chtrsol')
        call getmem1d(idust,1,nbin,'mod_che_common:idust')
        call getmem1d(isslt,1,sbin,'mod_che_common:isslt')
        call getmem1d(icarb,1,5,'mod_che_common:icarb')
        call getmem2d(chtrsize,1,nbin,1,2,'mod_che_common:chtrsize')

        call getmem4d(chemall,jci1,jci2,ici1,ici2, &
                      1,kz,1,totsp,'mod_che_common:chemall')
        call getmem3d(srclp2,jci1,jci2,ici1,ici2,1,ntr,'mod_che_common:srclp2')
        call getmem4d(jphoto,jci1,jci2,ici1,ici2,1,kz, &
          1,nphoto,'che_common:jphoto')

        call getmem3d(dtrace,jce1,jce2,ice1,ice2,1,ntr,'che_common:dtrace')
        call getmem3d(wdlsc,jce1,jce2,ice1,ice2,1,ntr,'che_common:wdlsc')
        call getmem3d(wdcvc,jce1,jce2,ice1,ice2,1,ntr,'che_common:wdcvc')
 
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
          call getmem4d(cemisdiag,jce1,jce2, &
                        ice1,ice2,1,kz,1,ntr,'che_common:cemisdiag')
        end if
#if (defined VOC && defined CLM)
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
        ntr = nbin
        allocate(chtrname(nbin))
        chtrname(1:ntr)(1:6) = (/'DUST01','DUST02','DUST03','DUST04'/)
        iaerosol = 1
        if ( myid == italk ) write(stdout,*) 'DUST simulation'
      else if ( chemsimtype(1:4) == 'SSLT' ) then 
        ntr = sbin
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:6) = (/'SSLT01','SSLT02'/)
        iaerosol = 1
        if ( myid == italk ) write(stdout,*) 'SSLT simulation'
      else if ( chemsimtype(1:4) == 'DUSS' ) then
        ntr = sbin + nbin
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:6) = (/'DUST01','DUST02','DUST03','DUST04', &
                                 'SSLT01','SSLT02'/)
        iaerosol = 1
        if ( myid == italk ) write(stdout,*) 'DUSS simulation'
       else if ( chemsimtype(1:4) == 'CARB' ) then 
        ntr = 4
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:6) = (/'BC_HL ','BC_HB ','OC_HL ','OC_HB '/)
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
        ntr = 12 
        allocate(chtrname(ntr))
        iaerosol = 1
        ioxclim  = 1
        chtrname(1:ntr)(1:6) = (/'BC_HL ','BC_HB ','OC_HL ','OC_HB ', &
                                 'SO2   ','SO4   ','DUST01','DUST02', &
                                 'DUST03','DUST04','SSLT01','SSLT02' /)
        if ( myid == italk ) write(stdout,*) 'AERO simulation'
      else if ( chemsimtype(1:4) == 'DCCB' ) then 
        ntr = 49
        allocate(chtrname(ntr))      
        chtrname(1:ntr)(1:6) = (/'SO2   ','SO4   ','NH3   ','O3    ', &
                                 'NO2   ','NO    ','CO    ','H2O2  ', &
                                 'HNO3  ','N2O5  ','HCHO  ','ALD2  ', &
                                 'ISOP  ','C2H6  ','PAR   ','ACET  ', &
                                 'MOH   ','OLT   ','OLI   ','TOLUE ', &
                                 'XYL   ','ETHE  ','PAN   ','CH4   ', &
                                 'MGLY  ','CRES  ','OPEN  ','ISOPRD', &
                                 'ONIT  ','HCOOH ','RCOOH ','CH3OOH', &
                                 'ETHOOH','ROOH  ','HONO  ','HNO4  ', &
                                 'XO2   ','DUST01','DUST02','DUST03', &
                                 'DUST04','BC_HL ','BC_HB ','OC_HL ', &
                                 'OC_HB ','SSLT01','SSLT02','ANO3  ', &
                                 'ANH4  ' /)
        iaerosol = 1
        igaschem = 1
        iisoropia = 1
        if ( myid == italk ) write(stdout,*) 'DCCB simulation'
      else if ( chemsimtype(1:4) == 'CBMZ' ) then 
        ntr = 37
        allocate(chtrname(ntr))      
        chtrname(1:ntr)(1:6) = (/'SO2   ','SO4   ','NH3   ','O3    ', &
                                 'NO2   ','NO    ','CO    ','H2O2  ', &
                                 'HNO3  ','N2O5  ','HCHO  ','ALD2  ', &
                                 'ISOP  ','C2H6  ','PAR   ','ACET  ', &
                                 'MOH   ','OLT   ','OLI   ','TOLUE ', &
                                 'XYL   ','ETHE  ','PAN   ','CH4   ', &
                                 'MGLY  ','CRES  ','OPEN  ','ISOPRD', &
                                 'ONIT  ','HCOOH ','RCOOH ','CH3OOH', &
                                 'ETHOOH','ROOH  ','HONO  ','HNO4  ', &
                                 'XO2   ' /)
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
             'DUST SSLT DUSS CARB SULF SUCA AERO CBMZ DCCB POLLEN'
        end if
        call fatal(__FILE__,__LINE__,'INVALID CHEM CONFIGURATION')
      end if

    end subroutine chem_config

end module mod_che_common
