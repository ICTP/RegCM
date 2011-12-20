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
!
module mod_che_ncio
!
  use mod_realkinds
  use mod_dynparam
  use mod_mpmessage
  use mod_che_indices
  use mod_che_common
  use netcdf
!
  private
!
  public :: read_texture , read_aerosol , read_emission, recc
  public :: prepare_chem_out, init_mod_che_ncio, writerec_che2
!
  integer :: istatus
  integer :: recc

   integer , parameter :: n_chevar = 18
   integer, parameter :: n_optvar = 9
   integer, dimension(:), allocatable :: ncche     
   integer , dimension(n_chevar) :: ichevar
   integer, dimension(n_optvar) ::ioptvar 
   integer, dimension(9) :: idims 
 integer :: icherec, ioptrec

 type(rcm_time_and_date) , save :: icherefdate

        real(8) :: rpt,cfd

        integer :: o_is
        integer :: o_ie
        integer :: o_js
        integer :: o_je
        integer :: o_ni
        integer :: o_nj
         integer :: o_nz
        logical :: lwrap 

        real(4) , dimension(:) , allocatable :: hsigma
        real(4) , dimension(:,:) , allocatable :: ioxlat
        real(4) , dimension(:,:) , allocatable :: ioxlon
        real(4) , dimension(:,:) , allocatable :: iotopo
 
        real(4) , dimension(:,:,:) , allocatable :: dumio
      

        real(4) , dimension(2) :: latrange
        real(4) , dimension(2) :: lonrange

  contains


 subroutine init_mod_che_ncio(lband)
          implicit none
          logical , intent(in) :: lband
          character(3) :: sbstring
 

          if (lband) then
            o_is = 2
            o_ie = iy-1
            o_js = 1
            o_je = jx
            o_ni = iy-2
            o_nj = jx
            o_nz = kz
            lwrap = .true.
          else
            o_is = 2
            o_ie = iy-1
            o_js = 2
            o_je = jx-1
            o_ni = iy-2
            o_nj = jx-2
            o_nz = kz
            lwrap = .false.
          end if

          allocate(hsigma(o_nz))
          allocate(ioxlat(o_nj,o_ni))
          allocate(ioxlon(o_nj,o_ni))
          allocate(iotopo(o_nj,o_ni))
          allocate(dumio(o_nj,o_ni,o_nz))
          
        end subroutine init_mod_che_ncio



  subroutine read_texture(nats,texture)
    implicit none

    integer , intent(in) :: nats
    real(dp) , dimension(iy,jx,nats) , intent(out) :: texture

    integer :: ivarid
    integer :: i , j , n
    integer :: idmin
    integer , dimension(3) :: istart , icount
    character(256) :: dname
    real(sp), dimension(jx,iy) ::  toto

    dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
    istatus = nf90_open(dname, nf90_nowrite, idmin)
    call check_ok(__FILE__,__LINE__, &
                  'Error Opening Domain file '//trim(dname),'DOMAIN FILE OPEN')
    istatus = nf90_inq_varid(idmin, 'texture_fraction', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable texture_fraction miss', &
                  'DOMAIN FILE')
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = iy
    icount(1) = jx
    do n = 1 , nats
      istart(3) = n
      istatus = nf90_get_var(idmin, ivarid, toto, istart, icount)
      call check_ok(__FILE__,__LINE__,'Variable texture_frac read error', &
                    'DOMAIN FILE')
      do j = 1 , jx
        do i = 1 , iy
          texture(i,j,n) = dble(toto(j,i))*0.01D0
          if (texture(i,j,n)<d_zero) texture(i,j,n)=d_zero
        end do
      end do
    end do
    istatus = nf90_close(idmin)
    call check_ok(__FILE__,__LINE__,'Domain file close error','DOMAIN FILE')
  end subroutine read_texture

  subroutine read_aerosol(chtrname,chemsrc)
    implicit none

    character(256) :: aername
    character(5) , dimension(ntr) , intent(in) :: chtrname
    real(dp) , dimension(iy,jx,12,ntr) , intent(out) :: chemsrc

    integer :: ncid , ivarid
    real(sp) , dimension(jx,iy) :: toto
    character(5) :: aerctl
    integer , dimension(3) :: istart , icount
    integer :: itr , i , j , m

    aername = trim(dirglob)//pthsep//trim(domname)//'_AERO.nc'
    istatus = nf90_open(aername, nf90_nowrite, ncid)
    call check_ok(__FILE__,__LINE__, &
         'Error Opening Aerosol file '//trim(aername),'AEROSOL FILE OPEN')

    do itr = 1 , ntr
      aerctl = chtrname(itr)
      write (aline, *) itr , aerctl
      call say
      if ( aerctl(1:4) /= 'DUST') then
        if ( aerctl(1:3) == 'SO2' ) then
          if ( aertyp(4:4) == '1' ) then
            istatus = nf90_inq_varid(ncid, 'so2', ivarid)
            call check_ok(__FILE__,__LINE__, &
                          'Variable so2 miss','AEROSOL FILE')
            istatus = nf90_get_var(ncid, ivarid, toto)
            call check_ok(__FILE__,__LINE__, &
                          'Variable so2 read error','AEROSOL FILE')
            do m = 1 , 12
              do j = 1 , jx
                do i = 1 , iy
                  chemsrc(i,j,m,itr) = dble(toto(j,i))
                end do
              end do
            end do
          end if
          if ( aertyp(5:5) == '1' ) then
            istatus = nf90_inq_varid(ncid, 'so2_monthly', ivarid)
            call check_ok(__FILE__,__LINE__, &
                          'Variable so2_mon miss','AEROSOL FILE')
            istart(1) = 1
            istart(2) = 1
            icount(1) = jx
            icount(2) = iy
            icount(3) = 1
            do m = 1 , 12
              istart(3) = m
              istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
              call check_ok(__FILE__,__LINE__, &
                            'Variable so2_mon read err','AEROSOL FILE')
              do j = 1 , jx
                do i = 1 , iy
                  chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + dble(toto(j,i))
                end do
              end do
            end do
          end if
        else if ( aerctl(1:2) == 'BC' ) then
          if ( aertyp(4:4) == '1' ) then
            istatus = nf90_inq_varid(ncid, 'bc', ivarid)
            call check_ok(__FILE__,__LINE__, &
                          'Variable bc miss','AEROSOL FILE')
            istatus = nf90_get_var(ncid, ivarid, toto)
            call check_ok(__FILE__,__LINE__, &
                          'Variable bc read error','AEROSOL FILE')
            do m = 1 , 12
              do j = 1 , jx
                do i = 1 , iy
                  chemsrc(i,j,m,itr) = dble(toto(j,i))
                end do
              end do
            end do
          end if
          if ( aertyp(5:5) == '1' ) then
            istatus = nf90_inq_varid(ncid, 'bc_monthly', ivarid)
            call check_ok(__FILE__,__LINE__, &
                          'Variable bc_mon miss','AEROSOL FILE')
            istart(1) = 1
            istart(2) = 1
            icount(1) = jx
            icount(2) = iy
            icount(3) = 1
            do m = 1 , 12
              istart(3) = m
              istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
              call check_ok(__FILE__,__LINE__, &
                            'Variable bc_mon read err','AEROSOL FILE')
              do j = 1 , jx
                do i = 1 , iy
                  chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + dble(toto(j,i))
                end do
              end do
            end do
          end if
        else if ( aerctl(1:2) == 'OC' ) then
          if ( aertyp(4:4) == '1' ) then
            istatus = nf90_inq_varid(ncid, 'oc', ivarid)
            call check_ok(__FILE__,__LINE__, &
                          'Variable oc miss','AEROSOL FILE')
            istatus = nf90_get_var(ncid, ivarid, toto)
            call check_ok(__FILE__,__LINE__, &
                          'Variable oc read error','AEROSOL FILE')
            do m = 1 , 12
              do j = 1 , jx
                do i = 1 , iy
                  chemsrc(i,j,m,itr) = dble(toto(j,i))
                end do
              end do
            end do
          end if
          if ( aertyp(5:5) == '1' ) then
            istatus = nf90_inq_varid(ncid, 'oc_monthly', ivarid)
            call check_ok(__FILE__,__LINE__, &
                          'Variable oc_mon miss','AEROSOL FILE')
            istart(1) = 1
            istart(2) = 1
            icount(1) = jx
            icount(2) = iy
            icount(3) = 1
            do m = 1 , 12
              istart(3) = m
              istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
              call check_ok(__FILE__,__LINE__, &
                            'Variable oc_mon read err','AEROSOL FILE')
              do j = 1 , jx
                do i = 1 , iy
                  chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + dble(toto(j,i))
                end do
              end do
            end do
          end if
        end if
      end if
    end do

    istatus = nf90_close(ncid)
    call check_ok(__FILE__,__LINE__, &
                  'Error Close Aerosol file '//trim(aername), &
                  'AEROSOL FILE CLOSE')

  end subroutine read_aerosol

  subroutine read_emission(lmonth,echemsrc)
    implicit none

    integer , intent(in) :: lmonth

    real(dp) , dimension(iy,jx,12,ntr) , intent(inout) :: echemsrc

    character(256) :: aername
 
    integer :: ncid 
    integer , dimension(3) :: istart , icount
    integer :: itr , ivarid
 
! FAB: remember for now, we have 1 emission file containing all monthly emission for the whole simulation period
! change that in the future. Also lmonth is not really necessary here, but KEEP THIS DIMENSION FOR HIGHER TEMPORAL RESOLUTION INVENTORIES 
!all aggregations / lumping should in the future be done in emission preproc 


    aername = trim(dirglob)//pthsep//trim(domname)//'_CHEMISS.nc'

    print*, 'Opening ch. emission file : RETURN !!!!', aername
    !CARE FAB shut down the emission reading temporarily 

    return

    istatus = nf90_open(aername, nf90_nowrite, ncid)
    call check_ok(__FILE__,__LINE__, &
                  'Error Opening chem emissiom file '//trim(aername), &
                  'CHE EMISS FILE OPEN ERROR')

    !*** intialized in start_chem
    !*** Advice record counter
    recc = recc + 1
 


    istart(1) = 1
    istart(2) = 1
    istart(3) = recc
    icount(1) = jx
    icount(2) = iy
    icount(3) = 1

    ! NO emission                  
    if ( ino /= 0 ) then
      call rvar(ncid,istart,icount,ino,lmonth,echemsrc, &
                'a_NO',.false.,'bio_nox')
    end if
    ! CO emission
    if ( ico /= 0 ) then
      call rvar(ncid,istart,icount,ico,lmonth,echemsrc, &
                'a_CO',.false.,'bio_co','o_co')
    print*, 'FAB emis testco','ico', maxval(echemsrc)
    end if
    ! HCHO emission                  
    if ( ihcho /= 0 ) then
      call rvar(ncid,istart,icount,ihcho,lmonth,echemsrc,'a_HCHO',.false.)
    end if
    ! ACET emission                  
    if ( iacet /= 0 ) then
      call rvar(ncid,istart,icount,iacet,lmonth,echemsrc, &
                'a_ACET',.false.,'bio_acet')
    end if
    ! SO2 emission
    if ( iso2 /= 0 ) then
      call rvar(ncid,istart,icount,iso2,lmonth,echemsrc, &
                'a_SO2',.false.,'b_so2')
    end if
    ! CH4
    if ( ich4 /= 0 ) then
      call rvar(ncid,istart,icount,ich4,lmonth,echemsrc, &
                'a_ch4',.false.,'bio_ch4')
    end if
    ! Ethane
    if ( ic2h6 /= 0 ) then
      call rvar(ncid,istart,icount,ic2h6,lmonth,echemsrc, &
                'a_ETHANE',.false.,'bio_c2h6','o_c2h6')
    end if
    ! PAR
    if ( ipar /= 0 ) then
      call rvar(ncid,istart,icount,ipar,lmonth,echemsrc, &
                'a_c3h8',.false.,'a_butane','bio_c3h8','o_c3h8')
    end if
    ! Ethene
    if ( iethe /= 0 ) then
      call rvar(ncid,istart,icount,iethe,lmonth,echemsrc, &
                'a_ETHENE',.false.,'bio_c2h4','o_c2h4')
    end if
    ! Termenal Alkene
    if ( iolt /= 0 ) then
      call rvar(ncid,istart,icount,iolt,lmonth,echemsrc, &
                'a_PROPENE',.false.,'bio_c3h6')
    end if
    ! Internal Alkene
    if ( ioli /= 0 ) then
      call rvar(ncid,istart,icount,ioli,lmonth,echemsrc,'a_BIGENE',.true.)
    end if
    ! Isoprene
    if ( iisop /= 0 ) then
      call rvar(ncid,istart,icount,iisop,lmonth,echemsrc,'bio_isop',.false.)
    end if
    ! Toluene
    if ( itolue /= 0 ) then
      istatus = nf90_inq_varid(ncid, 'a_XYLENE', ivarid)
      if ( istatus == nf90_noerr ) then
        call rvar(ncid,istart,icount,itolue,lmonth,echemsrc, &
                  'a_TOLUENE',.true.,'b_TOLUENE')
      end if
    end if
    ! Xylene
    if ( ixyl /= 0 ) then
      istatus = nf90_inq_varid(ncid, 'a_TOLUENE', ivarid)
      if ( istatus == nf90_noerr ) then
        call rvar(ncid,istart,icount,ixyl,lmonth,echemsrc, &
                  'a_XYLENE',.true.)
      else
        call rvar(ncid,istart,icount,ixyl,lmonth,echemsrc, &
                  'a_TOLUENE',.true.,'b_XYLENE')
      end if
    end if
    ! Acetaldehyde
    if ( iald2 /= 0 ) then
      call rvar(ncid,istart,icount,iald2,lmonth,echemsrc,'a_ALD2',.false.)
    end if
    ! Methanol + Ethanol
    if ( imoh /= 0 ) then
      call rvar(ncid,istart,icount,imoh,lmonth,echemsrc, &
                'a_MOH',.false.,'b_ch3oh','bio_ch3oh')
    end if           
    ! DMS
    if ( idms /= 0 ) then
      call rvar(ncid,istart,icount,idms,lmonth,echemsrc,'o_DMS',.false.)
    end if
    ! OC and BC anthropogenic + biomass burning
    if ( ibchb /= 0 ) then
      call rvar(ncid,istart,icount,ibchb,lmonth,echemsrc, &
                'b_BC',.false.,'a_BC')
    end if
    if ( iochb /= 0 ) then
      call rvar(ncid,istart,icount,iochb,lmonth,echemsrc, &
                'b_OC',.false.,'a_OC')
    end if

 
    istatus = nf90_close(ncid)
    call check_ok(__FILE__,__LINE__, &
                  'Error Closing Chem emission file '//trim(aername), &
                  'CH EMISS FILE CLOSE ERROR')

  end subroutine read_emission

  subroutine rvar(ncid,istart,icount,ind,lmonth,echemsrc,cna,lh,cnb,cnc,cnd)
    integer , intent(in) :: ncid
    integer , dimension(3) , intent(in) :: istart , icount
    integer , intent(in) :: lmonth
    real(dp) , dimension(iy,jx,12,ntr) , intent(inout) :: echemsrc
    logical , intent(in) :: lh
    character(len=*) , intent(in) :: cna
    character(len=*) , intent(in) , optional :: cnb
    character(len=*) , intent(in) , optional :: cnc
    character(len=*) , intent(in) , optional :: cnd
    integer :: ivarid 
    real(sp) , dimension(jx,iy) :: toto
    integer :: i , j

    istatus = nf90_inq_varid(ncid, cna, ivarid)     
    call check_ok(__FILE__,__LINE__, &
                  'Variable '//cna//' miss','CHEM_EMISS FILE')
    istatus = nf90_get_var(ncid,ivarid,toto)
    call check_ok(__FILE__,__LINE__, &
                  'Variable '//cna//' read err','CHEM_EMISS FILE')
    if ( lh ) then  ! half of lumped Aromatics
      do j = 1 , jx
        do i = 1 , iy
          echemsrc(i,j,lmonth,ind) = d_half*toto(j,i)
        end do
      end do
    else
      do j = 1 , jx
        do i = 1 , iy
          echemsrc(i,j,lmonth,ind) = toto(j,i)
        end do
      end do
    end if
    if ( present(cnb) ) then
      istatus = nf90_inq_varid(ncid, cnb, ivarid)
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cnb//' miss','CHEM_EMISS FILE')
      istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cnb//' read err','CHEM_EMISS FILE')
      do j = 1 , jx
        do i = 1 , iy
          echemsrc(i,j,lmonth,ind) = toto(j,i) + echemsrc(i,j,lmonth,ind)
        end do
      end do
    end if
    if ( present(cnc) ) then
      istatus = nf90_inq_varid(ncid, cnc, ivarid)
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cnc//' miss','CHEM_EMISS FILE')
      istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cnc//' read err','CHEM_EMISS FILE')
      do j = 1 , jx
        do i = 1 , iy
          echemsrc(i,j,lmonth,ind) = toto(j,i) + echemsrc(i,j,lmonth,ind)
        end do
      end do
    end if
    if ( present(cnd) ) then
      istatus = nf90_inq_varid(ncid, cnd, ivarid)
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cnd//' miss','CHEM_EMISS FILE')
      istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cnd//' read err','CHEM_EMISS FILE')
      do j = 1 , jx
        do i = 1 , iy
          echemsrc(i,j,lmonth,ind) = toto(j,i) + echemsrc(i,j,lmonth,ind)
        end do
      end do
    end if
  end subroutine rvar


!------------------------------------------------------------------------------      
!       IROUTINE: prepare_chem_out
!       SUBROUTINE INTERFACE:

        subroutine prepare_chem_out(idate, ifrest)
    
          use netcdf
!     !DESCRIPTION:
!     prepare the dimensions variables, write global attributes
!     define the chemistry variables


          implicit none
           type(rcm_time_and_date) , intent(in) :: idate
          logical, intent(in) :: ifrest

          integer  :: itr
          character(128)::cdum
          character(8) ::chevarnam
          character(64) :: title
          character(32) :: fbname , csdate
          character(64) :: cmethodmax , cmethodmin
          character(16) :: fterr
          character(256) :: ofname , history
          integer , dimension(8) :: tvals
          real(4) :: hptop , rdum1
          real(4) , dimension(2) :: trlat , rdum2
          real(4) , dimension(iysg) :: yiy
          real(4) , dimension(jxsg) :: xjx
          integer :: ncid,ivarid
          integer , dimension(3) :: izvar
          integer , dimension(2) :: ivvar
          integer , dimension(4) :: isrvvar
          integer , dimension(5) :: illtpvar
          integer :: itvar , iyy , im , id , ih , i , j,k
          integer :: l1,l2,ibin,jbin

          integer , dimension(9) :: tyx
          integer , dimension(9) :: tzyx
          

          real(4) , dimension(jx,iy) :: sp2d
          real(4) , dimension(iy,jx) :: xlat,xlon,topo
          real(4) , dimension(kz)    :: ppp


          if (.not.allocated(ncche)) then
                  allocate( ncche(ntr))
!                  ncche(:) = -1
          end if

ibin = 0
jbin = 0
! tracer loop , since we are generating one output per tracer 
          do itr = 1, ntr
!           
            ncid = ncche(itr)
     
            chevarnam =  chtrname(itr)

            if( chtrname (itr)== 'DUST')  then
             ibin = ibin+1
             write( chevarnam(5:6),'(I1)') ibin
            end if
            if( chtrname(itr) == 'SSLT')  then
             jbin = jbin+1
             write( chevarnam(5:6),'(I1)') jbin
            end if
 

            title = 'ICTP Regional Climatic model V4  '  &
                     //chevarnam//' output'
            icherefdate = idate
            icherec = 1

!          call close_chem(ncid,itr)

    write (fterr, '(a3,a)') chevarnam, ' FILE'
    write (fbname,'(a,a,i10)') trim(chevarnam), '.', toint10(idate)
    ofname = trim(dirout)//pthsep//trim(domname)// &
             '_'//trim(fbname)//'.nc'
 !   ctime = tochar(cordex_refdate)
    write (aline, *) 'Opening new output file ', trim(ofname)
    call say

#ifdef NETCDF4_HDF5
    istatus = nf90_create(ofname, &
              ior(ior(nf90_clobber,nf90_hdf5),nf90_classic_model),ncid)
#else
    istatus = nf90_create(ofname, nf90_clobber, ncid)
#endif


            ncche(itr) = ncid


!Start Global Attributes
    istatus = nf90_put_att(ncid, nf90_global, 'title', title)
    call check_ok(__FILE__,__LINE__,'Error add title', fterr)
    istatus = nf90_put_att(ncid, nf90_global, 'institution', 'ICTP')
    call check_ok(__FILE__,__LINE__,'Error add institution', fterr)
    istatus = nf90_put_att(ncid, nf90_global, 'source', &
               'RegCM Model '//'SVN_REV'//' simulation output')
    call check_ok(__FILE__,__LINE__,'Error add source', fterr)
    istatus = nf90_put_att(ncid, nf90_global, 'Conventions', 'CF-1.4')
    call check_ok(__FILE__,__LINE__,'Error add Conventions', fterr)
    call date_and_time(values=tvals)
    write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)') &
         tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,       &
         tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,             &
         ' : Created by RegCM model'
    istatus = nf90_put_att(ncid, nf90_global, 'history', history)
    call check_ok(__FILE__,__LINE__,'Error add history', fterr)
    istatus = nf90_put_att(ncid, nf90_global, 'references', &
               'http://eforge.escience-lab.org/gf/project/regcm')
    call check_ok(__FILE__,__LINE__,'Error add references', fterr)
    istatus = nf90_put_att(ncid, nf90_global, 'experiment', domname)
    call check_ok(__FILE__,__LINE__,'Error add experiment', fterr)
    istatus = nf90_put_att(ncid, nf90_global, 'projection', iproj)
    call check_ok(__FILE__,__LINE__,'Error add projection', fterr)
    if (iproj == 'LAMCON') then
      istatus = nf90_put_att(ncid, nf90_global, &
                   'grid_mapping_name', 'lambert_conformal_conic')
      call check_ok(__FILE__,__LINE__,'Error add grid_mapping_name',fterr)
    else if (iproj == 'POLSTR') then
      istatus = nf90_put_att(ncid, nf90_global, &
                   'grid_mapping_name', 'stereographic')
      call check_ok(__FILE__,__LINE__,'Error add grid_mapping_name',fterr)
    else if (iproj == 'NORMER') then
      istatus = nf90_put_att(ncid, nf90_global, &
                   'grid_mapping_name', 'mercator')
      call check_ok(__FILE__,__LINE__,'Error add grid_mapping_name',fterr)
    else if (iproj == 'ROTMER') then
      istatus = nf90_put_att(ncid, nf90_global, &
            'grid_mapping_name', 'rotated_latitude_longitude')
      call check_ok(__FILE__,__LINE__,'Error add grid_mapping_name',fterr)
    end if
    istatus = nf90_put_att(ncid, nf90_global, 'grid_size_in_meters', ds*d_1000)
    call check_ok(__FILE__,__LINE__,'Error add gridsize', fterr)
    istatus = nf90_put_att(ncid, nf90_global, &
                           'latitude_of_projection_origin', clat)
    call check_ok(__FILE__,__LINE__,'Error add clat', fterr)
    istatus = nf90_put_att(ncid, nf90_global,   &
                 'longitude_of_projection_origin', clon)
    call check_ok(__FILE__,__LINE__,'Error add clon', fterr)
    istatus = nf90_put_att(ncid, nf90_global,   &
                 'longitude_of_central_meridian', clon)
    call check_ok(__FILE__,__LINE__,'Error add gmtllon', fterr)
    if (iproj == 'ROTMER') then
      istatus = nf90_put_att(ncid, nf90_global, &
                   'grid_north_pole_latitude', plat)
      call check_ok(__FILE__,__LINE__,'Error add plat', fterr)
      istatus = nf90_put_att(ncid, nf90_global, &
                   'grid_north_pole_longitude', plon)
      call check_ok(__FILE__,__LINE__,'Error add plon', fterr)
    else if (iproj == 'LAMCON') then
      trlat(1) = real(truelatl)
      trlat(2) = real(truelath)
      istatus = nf90_put_att(ncid, nf90_global, 'standard_parallel', trlat)
      call check_ok(__FILE__,__LINE__,'Error add truelat', fterr)
    else if (iproj == 'NORMER') then
      istatus = nf90_put_att(ncid, nf90_global, 'standard_parallel', clat)
      call check_ok(__FILE__,__LINE__,'Error add truelat', fterr)
    else if (iproj == 'POLSTR') then
      trlat(1) = 1.0
      istatus = nf90_put_att(ncid, nf90_global, &
                 'scale_factor_at_projection_origin', trlat(1:1))
      call check_ok(__FILE__,__LINE__,'Error add scfac', fterr)
    end if
!
!         ADD RUN PARAMETERS
!
!FAB    istatus = nf90_put_att(ncid, nf90_global, 'model_IPCC_scenario', scenario)
!    call check_ok(__FILE__,__LINE__,'Error add scenario', fterr)

!!$                        !!!Je skip ce bloc pour l'instant !!!
!!$    call cdumlbcs
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_boundary_conditions' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add lbcs', fterr)
!!$    call cdumcums
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_cumulous_convection_scheme' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add icup', fterr)
!!$    if (icup == 2 .or. icup == 99 .or. icup == 98) then
!!$      call cdumcumcl
!!$      istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_convective_closure_assumption' , trim(cdum))
!!$      call check_ok(__FILE__,__LINE__,'Error add igcc', fterr)
!!$    end if
!!$    call cdumpbl
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_boundary_layer_scheme' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add ibltyp', fterr)
!!$    call cdummoist
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_moist_physics_scheme' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add ipptls', fterr)
!!$    call cdumocnflx
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_ocean_flux_scheme' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add iocnflx', fterr)
!!$    call cdumpgfs
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_pressure_gradient_force_scheme' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add ipgf', fterr)
!!$    call cdumemiss
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_use_emission_factor' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add iemiss', fterr)
!!$    call cdumlakes
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_use_lake_model' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add lakemod', fterr)
!!$    call cdumchems
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_chemistry' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add ichem', fterr)
!!$    call cdumdcsst
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_diurnal_cycle_sst' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add dcsst', fterr)
!!$    call cdumseaice
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_seaice_effect' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add seaice', fterr)
!!$    call cdumdesseas
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_seasonal_desert_albedo_effect' , trim(cdum))
!!$
!!$    call check_ok(__FILE__,__LINE__,'Error add desseas', fterr)



    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_simulation_initial_start' , tochar(globidate1))
    call check_ok(__FILE__,__LINE__,'Error add globidate1', fterr)

!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_simulation_start' , tochar(idate1))
!!$    call check_ok(__FILE__,__LINE__,'Error add idate1', fterr)
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_simulation_expected_end' , tochar(idate2))
!!$    call check_ok(__FILE__,__LINE__,'Error add idate2', fterr)


    if (ifrest) then
      cdum = 'Yes'
    else
      cdum = 'No'
    end if
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_simulation_is_a_restart' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add ifrest', fterr)


!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_timestep_in_seconds' , dt)
!!$    call check_ok(__FILE__,__LINE__,'Error add dt', fterr)
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_timestep_in_minutes_solar_rad_calc' , dtrad)
!!$    call check_ok(__FILE__,__LINE__,'Error add dtrad', fterr)
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_timestep_in_seconds_bats_calc' , dtsrf)
!!$    call check_ok(__FILE__,__LINE__,'Error add dtsrf', fterr)
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_timestep_in_hours_radiation_calc' , dtabem)
!!$    call check_ok(__FILE__,__LINE__,'Error add dtabem', fterr)
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_timestep_in_hours_boundary_input' , ibdyfrq)
!!$    call check_ok(__FILE__,__LINE__,'Error add dtabem', fterr)


!End of Global Attributes
!
!         ADD DIMENSIONS
!
      istatus = nf90_def_dim(ncid, 'iy', o_ni, idims(2))
      call check_ok(__FILE__,__LINE__,'Error create dim iy', fterr)
      istatus = nf90_def_dim(ncid, 'jx', o_nj, idims(1))
      call check_ok(__FILE__,__LINE__,'Error create dim jx', fterr)
    istatus = nf90_def_dim(ncid, 'time', nf90_unlimited, idims(3))
    call check_ok(__FILE__,__LINE__,'Error create dim time', fterr)
    istatus = nf90_def_dim(ncid, 'kz', kz, idims(4))
    call check_ok(__FILE__,__LINE__,'Error create dim kz', fterr)

!
!         OUT TYPE DEPENDENT DIMENSIONS
!
    istatus = nf90_def_var(ncid, 'sigma', nf90_float, idims(4), izvar(1))
    call check_ok(__FILE__,__LINE__,'Error add var sigma', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'standard_name', &
                         'atmosphere_sigma_coordinate')
    call check_ok(__FILE__,__LINE__,'Error add sigma standard_name', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'long_name', &
                         'Sigma at model layers')
    call check_ok(__FILE__,__LINE__,'Error add sigma long_name', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'units', '1')
    call check_ok(__FILE__,__LINE__,'Error add sigma units', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'axis', 'Z')
    call check_ok(__FILE__,__LINE__,'Error add sigma axis', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'positive', 'down')
    call check_ok(__FILE__,__LINE__,'Error add sigma positive', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'formula_terms',  &
                         'sigma: sigma ps: ps ptop: ptop')
    call check_ok(__FILE__,__LINE__,'Error add sigma formula_terms', fterr)
   
    istatus = nf90_def_var(ncid, 'ptop', nf90_float, varid=izvar(2))
    call check_ok(__FILE__,__LINE__,'Error add var ptop', fterr)
    istatus = nf90_put_att(ncid, izvar(2), 'standard_name', 'air_pressure')
    call check_ok(__FILE__,__LINE__,'Error add ptop standard_name', fterr)
    istatus = nf90_put_att(ncid, izvar(2), 'long_name', 'Pressure at model top')
    call check_ok(__FILE__,__LINE__,'Error add ptop long_name', fterr)
    istatus = nf90_put_att(ncid, izvar(2), 'units', 'hPa')
    call check_ok(__FILE__,__LINE__,'Error add ptop units', fterr)
    istatus = nf90_def_var(ncid, 'iy', nf90_float, idims(2), ivvar(1))
    call check_ok(__FILE__,__LINE__,'Error add var iy', fterr)
    istatus = nf90_put_att(ncid, ivvar(1), 'standard_name', &
                           'projection_y_coordinate')
    call check_ok(__FILE__,__LINE__,'Error add iy standard_name', fterr)
    istatus = nf90_put_att(ncid, ivvar(1), 'long_name', &
                           'y-coordinate in Cartesian system')
    call check_ok(__FILE__,__LINE__,'Error add iy long_name', fterr)
    istatus = nf90_put_att(ncid, ivvar(1), 'units', 'km')
    call check_ok(__FILE__,__LINE__,'Error add iy units', fterr)
    istatus = nf90_def_var(ncid, 'jx', nf90_float, idims(1), ivvar(2))
    call check_ok(__FILE__,__LINE__,'Error add var jx', fterr)
    istatus = nf90_put_att(ncid, ivvar(2), 'standard_name', &
                         'projection_x_coordinate')
    call check_ok(__FILE__,__LINE__,'Error add jx standard_name', fterr)
    istatus = nf90_put_att(ncid, ivvar(2), 'long_name', &
                         'x-coordinate in Cartesian system')
    call check_ok(__FILE__,__LINE__,'Error add jx long_name', fterr)
    istatus = nf90_put_att(ncid, ivvar(2), 'units', 'km')
    call check_ok(__FILE__,__LINE__,'Error add jx units', fterr)
    istatus = nf90_def_var(ncid, 'xlat', nf90_float, idims(1:2), illtpvar(1))
    call check_ok(__FILE__,__LINE__,'Error add var xlat', fterr)
    istatus = nf90_put_att(ncid, illtpvar(1), 'standard_name', 'latitude')
    call check_ok(__FILE__,__LINE__,'Error add xlat standard_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(1), 'long_name', &
                         'Latitude at cross points')
    call check_ok(__FILE__,__LINE__,'Error add xlat long_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(1), 'units', 'degrees_north')
    call check_ok(__FILE__,__LINE__,'Error add xlat units', fterr)
    istatus = nf90_def_var(ncid, 'xlon', nf90_float, idims(1:2), illtpvar(2))
    call check_ok(__FILE__,__LINE__,'Error add var xlon', fterr)
    istatus = nf90_put_att(ncid, illtpvar(2), 'standard_name', 'longitude')
    call check_ok(__FILE__,__LINE__,'Error add xlon standard_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(2), 'long_name', &
                         'Longitude at cross points')
    call check_ok(__FILE__,__LINE__,'Error add xlon long_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(2), 'units', 'degrees_east')
    call check_ok(__FILE__,__LINE__,'Error add xlon units', fterr)
    istatus = nf90_def_var(ncid, 'topo', nf90_float, idims(1:2), illtpvar(3))
    call check_ok(__FILE__,__LINE__,'Error add var topo', fterr)
    istatus = nf90_put_att(ncid, illtpvar(3), 'standard_name', &
                         'surface_altitude')
    call check_ok(__FILE__,__LINE__,'Error add topo standard_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(3), 'long_name',     &
                         'Domain surface elevation')
    call check_ok(__FILE__,__LINE__,'Error add topo long_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(3), 'units', 'm')
    call check_ok(__FILE__,__LINE__,'Error add topo units', fterr)
    istatus = nf90_put_att(ncid, illtpvar(3), 'coordinates', 'xlat xlon')
    call check_ok(__FILE__,__LINE__,'Error add topo coord', fterr)
    istatus = nf90_put_att(ncid, illtpvar(3), 'grid_mapping', 'rcm_map')
    call check_ok(__FILE__,__LINE__,'Error add topo grid_mapping', fterr)
    istatus = nf90_def_var(ncid, 'mask', nf90_float, idims(1:2), illtpvar(4))
    call check_ok(__FILE__,__LINE__,'Error add var mask', fterr)
    istatus = nf90_put_att(ncid, illtpvar(4), 'standard_name', 'landmask')
    call check_ok(__FILE__,__LINE__,'Error add mask standard_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(4), 'long_name',     &
                         'Domain land/ocean mask')
    call check_ok(__FILE__,__LINE__,'Error add mask long_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(4), 'units', '1')
    call check_ok(__FILE__,__LINE__,'Error add mask units', fterr)
    istatus = nf90_put_att(ncid, illtpvar(4), 'coordinates', 'xlat xlon')
    call check_ok(__FILE__,__LINE__,'Error add mask coord', fterr)
    istatus = nf90_put_att(ncid, illtpvar(4), 'grid_mapping', 'rcm_map')
    call check_ok(__FILE__,__LINE__,'Error add mask grid_mapping', fterr)
    istatus = nf90_def_var(ncid, 'time', nf90_double, idims(3:3), itvar)
    call check_ok(__FILE__,__LINE__,'Error add var time', fterr)
    istatus = nf90_put_att(ncid, itvar, 'standard_name', 'time')
    call check_ok(__FILE__,__LINE__,'Error add time standard_name', fterr)
    istatus = nf90_put_att(ncid, itvar, 'long_name', 'time')
    call check_ok(__FILE__,__LINE__,'Error add time long_name', fterr)
    istatus = nf90_put_att(ncid, itvar, 'calendar', calstr(idate%calendar))
    call check_ok(__FILE__,__LINE__,'Error add time calendar', fterr)
!    istatus = nf90_put_att(ncid, itvar, 'units', 'hours since '//ctime)
!    call check_ok(__FILE__,__LINE__,'Error add time units', fterr)
   
    istatus = nf90_def_var(ncid, 'ps', nf90_float, idims(1:3), illtpvar(5))
    call check_ok(__FILE__,__LINE__,'Error add var ps', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'standard_name', &
                         'surface_air_pressure')
    call check_ok(__FILE__,__LINE__,'Error add ps standard_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'long_name', 'Surface pressure')
    call check_ok(__FILE__,__LINE__,'Error add ps long_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'units', 'hPa')
    call check_ok(__FILE__,__LINE__,'Error add ps units', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'coordinates', 'xlat xlon')
    call check_ok(__FILE__,__LINE__,'Error add ps coord', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'grid_mapping', 'rcm_map')
    call check_ok(__FILE__,__LINE__,'Error add ps grid_mapping', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'cell_methods', 'time: point')
    call check_ok(__FILE__,__LINE__,'Error add ps cell_methods', fterr)

    tyx = (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/)
    tzyx = (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/)
         

            ichevar = -1
            ichevar(1) = itvar
            ichevar(2) = illtpvar(4)

            call ch_addvara(ncid,chevarnam,chevarnam, &
                'atmosphere_mixing_ratio_of_tracer', &
                'Tracers mixing ratios','kg kg-1', &
                tzyx,.false.,ichevar(3))

            call ch_addvara(ncid,chevarnam,'wetdep_ls_flx', &
                'wet_deposition_from_large_scale_precip', &
                'Wet deposition LS','mg/m2/d', &
                tyx,.false.,ichevar(4))

            call ch_addvara(ncid,chevarnam,'wetdep_conv_flx', &
                'wet_deposition_from_convective_precip', &
                'Wet deposition CONV','mg/m2/d', &
                tyx,.false.,ichevar(5))

            call ch_addvara(ncid,chevarnam,'drydep_flx', &
                'dry_deposition', &
                'Dry deposition rate','mg/m2/d', &
                tyx,.false.,ichevar(6))

            call ch_addvara(ncid,chevarnam,'emiss_flx', &
                'surface_emission_rate', &
                'Emission rate','mg/m2/d', &
                tyx,.false.,ichevar(7))

             call ch_addvara(ncid,chevarnam,'drydep_vel', &
                'dry deposition velocity', &
                'dr. dep. vel','m.s-1', &
                tyx,.false.,ichevar(8))

    istatus = nf90_enddef(ncid)
    call check_ok(__FILE__,__LINE__,'Error End Definitions NetCDF output',fterr)

    istatus = nf90_put_var(ncid, izvar(1), hsigma)
    call check_ok(__FILE__,__LINE__,'Error var sigma write', fterr)
    hptop = real(ptop*d_10)
    istatus = nf90_put_var(ncid, izvar(2), hptop)
    call check_ok(__FILE__,__LINE__,'Error var ptop write', fterr)
  
      yiy(1) = -real((dble(o_ni-1)/2.0D0)*ds)
      xjx(1) = -real((dble(o_nj-1)/2.0D0)*ds)
      do i = 2 , o_ni
        yiy(i) = yiy(i-1)+real(ds)
      end do
      do j = 2 , o_nj
        xjx(j) = xjx(j-1)+real(ds)
      end do
      istatus = nf90_put_var(ncid, ivvar(1), yiy(1:o_ni))
      call check_ok(__FILE__,__LINE__,'Error var iy write', fterr)
      istatus = nf90_put_var(ncid, ivvar(2), xjx(1:o_nj))
      call check_ok(__FILE__,__LINE__,'Error var jx write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(1), ioxlat)
      call check_ok(__FILE__,__LINE__,'Error var xlat write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(2), ioxlon)
      call check_ok(__FILE__,__LINE__,'Error var xlon write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(3), iotopo)
      call check_ok(__FILE__,__LINE__,'Error var topo write', fterr)
!!$      istatus = nf90_put_var(ncid, illtpvar(4), iomask)
!!$      call check_ok(__FILE__,__LINE__,'Error var mask write', fterr)
!!$  
      ncche(itr) = ncid

! end of tracer loop       
       end do      

        end subroutine prepare_chem_out
!=====================================================================
        subroutine ch_addvara(ncid,ctype,vname,vst,vln,vuni,idims,lmiss, &
                           ivar)
          use netcdf
          implicit none
          integer , intent(in) :: ncid
          character(3) , intent(in) :: ctype
          character(len=*) , intent(in) :: vname
          character(len=*) , intent(in) :: vst , vln , vuni
          integer , dimension(5) , intent(in) :: idims
          logical , intent(in) :: lmiss
          integer , intent(out) :: ivar
          character(64) :: cdum
          real(4) , parameter :: fillv = -1E+34
          integer :: i , ndims
       
          ndims = 0
          do i = 1 , 5
            if (idims(i) > 0) ndims = ndims+1
          end do

          cdum = vname
          istatus = nf90_def_var(ncid, cdum, nf90_float, &
                             &   idims(1:ndims), ivar)
          call check_ok(__FILE__,__LINE__,'Error adding variable '//vname, &
                        ctype//' FILE ERROR')
#ifdef NETCDF4_HDF5
          istatus = nf90_def_var_deflate(ncid, ivar, 1, 1, 9)
        call  check_ok(__FILE__,__LINE__,'Error setting deflate on variable '//vname, &
                        ctype//' FILE ERROR')
#endif
          cdum = vst
          istatus = nf90_put_att(ncid, ivar, 'standard_name', cdum)
        call  check_ok(__FILE__,__LINE__,'Error adding '//vname//' standard_name', &
                        ctype//' FILE ERROR')
          cdum = vln
          istatus = nf90_put_att(ncid, ivar, 'long_name', cdum)
       call  check_ok(__FILE__,__LINE__,'Error adding '//vname//' long_name', &
                        ctype//' FILE ERROR')
          cdum = vuni
          istatus = nf90_put_att(ncid, ivar, 'units', cdum)
       call  check_ok(__FILE__,__LINE__,'Error adding '//vname//' units', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, ivar, 'coordinates', &
                            &  'xlat xlon')
       call  check_ok(__FILE__,__LINE__,'Error adding '//vname//' coordinates', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, ivar, 'grid_mapping', &
                            &  'rcm_map')
      call   check_ok(__FILE__,__LINE__,'Error adding '//vname//' grid_mapping', &
                        ctype//' FILE ERROR')
          if (lmiss) then
            istatus = nf90_put_att(ncid, ivar, '_FillValue', &
                              &  fillv)
         call  check_ok(__FILE__,__LINE__,'Error adding '//vname//' coordinates', &
                          ctype//' FILE ERROR')
          end if
        end subroutine ch_addvara

!==================================================================================
        subroutine writerec_che2(nx, ny, nnx, nny, nz, nt, chia, wdlsc, &
                                 wdcvc, ddsfc, cemtrac, drydepv, ps, idate)


!     DESCRIPTION
!     write the chemistry field to the corresbonding netcdf file                                
          use netcdf

          implicit none
          
type(rcm_time_and_date) , intent(in) :: idate
          integer , intent(in) :: nx , ny , nnx , nny , nz , nt 
          real(8) , dimension(iy,kz,jx,nt) , intent(in) :: chia
          real(8) , dimension(iy,jx) , intent(in) :: ps
          real(8) , dimension(iy,jx,nt), intent(in) :: wdlsc, wdcvc, ddsfc,  &
                                             &         cemtrac, drydepv

          integer :: n , k
          integer , dimension(5) :: istart , icount
                    
          real(8) , dimension(1) :: nctime
          type(rcm_time_interval) :: tdif
          character(len=36) :: ctime
          real(8) :: cfd2

          if (nx < o_nj .or. ny < o_ni .or. nz > o_nz) then
            write (6,*) 'Error writing record on CHE file'
            write (6,*) 'Expecting layers ', o_nz, 'x', o_nj, 'x', o_ni
            write (6,*) 'Got layers       ', nz, 'x', nx, 'x', ny
            call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
          end if


         do n=1,ntr        

       !     write (ctime, '(i10)') idate
       
          istart(1) = icherec
          icount(1) = 1

          ctime = tochar(idate)
          tdif = idate-icherefdate
          nctime(1) = tohours(tdif)

          istatus = nf90_put_var(ncche(n), ichevar(1), nctime, &
                                 istart(1:1), icount(1:1))
         call  check_ok(__FILE__,__LINE__,'Error writing itime '//ctime, 'CHE FILE ERROR')

          dumio(:,:,1) = transpose(ps(o_is:o_ie,o_js:o_je)+rpt)*10.0
          istart(3) = icherec
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj
          istatus = nf90_put_var(ncche(n), ichevar(2), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
         call  check_ok(__FILE__,__LINE__,'Error writing ps at '//ctime, 'CHE FILE ERROR')


            istart(4) = icherec
            istart(3) = 1
            istart(2) = 1
            istart(1) = 1
            icount(5) = 1
            icount(4) = 1
            icount(3) = o_nz
            icount(2) = o_ni
            icount(1) = o_nj

           !*** tracer concentration
            do k = 1 , nz
              dumio(:,:,k) = transpose(chia(o_is:o_ie,nz-k+1,o_js:o_je,n) / &
                                       ps(o_is:o_ie,o_js:o_je))
            end do
            istatus = nf90_put_var(ncche(n), ichevar(3), &
                                 dumio, istart, icount)

           call  check_ok(__FILE__,__LINE__,'Error writing '//chtrname(n)//' at '//ctime,&
                          'CHE FILE ERROR')


          
           !*** output 2-D gas chem fields
            istart(4) = 1
            istart(3) = icherec
            istart(2) = 1
            istart(1) = 1
            icount(4) = 1
            icount(3) = 1
            icount(2) = o_ni
            icount(1) = o_nj

! accumulated quantities between two output steps are converted to deposition/emission mean rate (mg /m2/per day)  
           cfd = 24./chemfrq
           cfd2 = dtche / (chemfrq *3600.)

          !*** wet deposition from large-scale precip
            dumio(:,:,1) = transpose(wdlsc(o_is:o_ie,o_js:o_je,n))*cfd
            istatus = nf90_put_var(ncche(n), ichevar(4), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
           call  check_ok(__FILE__,__LINE__,'Error writing wet dep LS at '//ctime,&
                          'CHE FILE ERROR')

          !*** wet deposition from convective precip
            dumio(:,:,1) = transpose(wdcvc(o_is:o_ie,o_js:o_je,n))*cfd
            istatus = nf90_put_var(ncche(n), ichevar(5), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
           call  check_ok(__FILE__,__LINE__,'Error writing wet dep CONV at '//ctime,&
                          'CHE FILE ERROR')

          !*** dry deposition
            dumio(:,:,1) = transpose(ddsfc(o_is:o_ie,o_js:o_je,n))*cfd
            istatus = nf90_put_var(ncche(n), ichevar(6), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
           call  check_ok(__FILE__,__LINE__,'Error writing dry dep '//ctime, 'CHE FILE ERROR')

          !*** emission rates
            dumio(:,:,1) = transpose(cemtrac(o_is:o_ie,o_js:o_je,n))*cfd
            istatus = nf90_put_var(ncche(n), ichevar(7), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
           call  check_ok(__FILE__,__LINE__,'Error writing emission rate '//ctime, 'CHE FILE ERROR')

          !*** dry dep vel 
            dumio(:,:,1) = transpose(drydepv(o_is:o_ie,o_js:o_je,n))*cfd2
            istatus = nf90_put_var(ncche(n), ichevar(8), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
           call  check_ok(__FILE__,__LINE__,'Error writing dr.dep.vel '//ctime, 'CHE FILE ERROR')

            istatus = nf90_sync(ncche(n))
            call check_ok(__FILE__,__LINE__,'Error sync at '//ctime, 'CHE FILE ERROR')

         end do     !main species looop

         icherec = icherec + 1


   end subroutine writerec_che2


    

  subroutine check_ok(f,l,m1,mf)
    implicit none
    character(*) , intent(in) :: f, m1 , mf
    integer , intent(in) :: l
    if (istatus /= nf90_noerr) then
      write (6,*) trim(m1)
      write (6,*) nf90_strerror(istatus)
      call fatal(f,l,trim(mf))
    end if
  end subroutine check_ok

end module mod_che_ncio
