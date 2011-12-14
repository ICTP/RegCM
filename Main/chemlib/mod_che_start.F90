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

  module mod_che_start 

    use mod_dynparam
    use mod_che_common
!    use mod_mpmessage
    use mod_che_indices
    use mod_che_bdyco
    use mod_che_wetdep
    use mod_che_carbonaer
    use mod_che_ncio

    
    implicit none

    public  :: start_chem

    contains

!----------------------------------------------------------------------------------------------------

subroutine start_chem (ifrest)
!  use mpi
  implicit none
 
  logical , intent(in) :: ifrest
  character(256)       :: namelistfile,prgname
  integer              :: ierr,err
  integer              :: itr,i,k,j,ibin,jbin,kbin

! je sais pas a quoi ca sert
!!$  if ( myid.eq.0 ) then
!!$
!!$  call getarg(0, prgname)
!!$  call getarg(1, namelistfile)
!!$  call init_mod_ncio2(.false.)
!!$
!!$  end if

!A : Intialise chemistry tracer indices         

        iso2 = 0
        iso4 = 0
        idms = 0
        imsa = 0
        ibchl = 0
        ibchb = 0
        iochl = 0
        iochb = 0


        io3       =  0
        ino       =  0
        ino2       =  0
        ino3      =  0
        ioh       =  0
        iho2      =  0
        ih2o2     =  0
        ihno2     =  0
        ihno3     =  0
        ihno4     =  0
        isulf     =  0
        ih2so4    =  0
        ihono     =  0
        in2o5     =  0
        ihc       =  0
        ihcr      =  0
        ic2h4     =  0
        ico       =  0
        ihcho     =  0
        iald2     =  0
        iethe     =  0
        ic2h6     =  0
        ic3h8     =  0
        iisop     =  0
        itolue    =  0
        ixyl      =  0
        inh3      =  0
        ipan      =  0
        irooh     =  0
        iacet     =  0
        ibenz     =  0
        inox      =  0
        ihox      =  0
        isox      =  0
        ich4      =  0
        ieoh      =  0
        imoh      =  0
        iaco2     =  0
        ico2      =  0
        in2o      =  0
        ipar      =  0
        iolt      =  0
        ioli      =  0

        !abt *** For initializing megan tracer biogenic voc mask  
        !    *** Mask not equal to zero when any MEGAN species is
        !    *** defined as a tracer within regcm.in  
        !    *** If not equal to zero then that compound will be
        !    *** used as a surface emission from MEGAN and not
        !    *** from inventory (see surftnd.F for use)
#if (defined VOC)
        if(ibvoc.eq.1) then
           if(.not.allocated(bvoc_trmask))  &
     &        allocate(bvoc_trmask(ntr))
           bvoc_trmask(:) = 0
        endif
#endif


ibin=0
jbin=0
kbin=0



        print*,'startchem',chtrname

        do itr = 1 , ntr

           if ( chtrname(itr).eq.'SO2' ) iso2 = itr
         
           if (chtrname(itr).eq. 'DMS')  idms = itr

           if ( chtrname(itr).eq.'SO4' ) then
             ! sulfate index is added to carb vector for treatment in drydep and wetdep
             ! sulfate effective diameter and bin is taken equal to ochl
                kbin = kbin + 1
                iso4 = itr
                icarb(kbin) = itr
                carbsiz(kbin,1) = reffochl - 0.1
                carbsiz(kbin,2) = reffochl + 0.1
           end if
           if ( chtrname(itr).eq.'BC_HL' ) then 
                kbin = kbin + 1
                ibchl = itr
                icarb(kbin) = itr
                carbsiz(kbin,1) = reffbchl - 0.06
                carbsiz(kbin,2) = reffbchl + 0.06
           end if
           if ( chtrname(itr).eq.'BC_HB' ) then 
                  kbin = kbin + 1
                  ibchb = itr
                  icarb(kbin) = itr
                  carbsiz(kbin,1) = reffbc - 0.01
                  carbsiz(kbin,2) = reffbc + 0.01
           end if
           if ( chtrname(itr).eq.'OC_HL' ) then
                   kbin = kbin + 1
                   iochl = itr
                   icarb(kbin) = itr
                   carbsiz(kbin,1) = reffochl - 0.1
                   carbsiz(kbin,2) = reffochl + 0.1
           end if
           if ( chtrname(itr).eq.'OC_HB' ) then
                   kbin = kbin + 1
                   iochb = itr
                   icarb(kbin) = itr 
                   carbsiz(kbin,1) = reffoc - 0.07
                   carbsiz(kbin,2) = reffoc + 0.07
           end if


           if ( chtrname(itr)(1:4).eq. 'DUST') then
              ibin = ibin + 1
              idust(ibin) = itr
           end if

           if ( chtrname(itr)(1:4).eq. 'SSLT') then
            jbin = jbin + 1
              isslt(jbin) = itr
           end if


           if  (chtrname(itr).eq. 'O3'    ) io3       = itr
           if  (chtrname(itr).eq. 'NO'    ) ino       = itr
           if  (chtrname(itr).eq. 'NO2'   ) ino2      = itr
           if  (chtrname(itr).eq. 'NO3'   ) ino3      = itr
           if  (chtrname(itr).eq. 'OH'    ) ioh       = itr
           if  (chtrname(itr).eq. 'HO2'   ) iho2      = itr
           if  (chtrname(itr).eq. 'H2O2'  ) ih2o2     = itr
           if  (chtrname(itr).eq. 'HNO2'  ) ihno2     = itr
           if  (chtrname(itr).eq. 'HNO3'  ) ihno3     = itr
           if  (chtrname(itr).eq. 'HNO4'  ) ihno4     = itr
           if  (chtrname(itr).eq. 'SULF'  ) isulf     = itr
           if  (chtrname(itr).eq. 'SO4'   ) iso4      = itr
           if  (chtrname(itr).eq. 'H2SO4' ) ih2so4    = itr
           if  (chtrname(itr).eq. 'HONO'  ) ihono     = itr
           if  (chtrname(itr).eq. 'N2O5'  ) in2o5     = itr
           if  (chtrname(itr).eq. 'HC'    ) ihc       = itr
           if  (chtrname(itr).eq. 'HCR'   ) ihcr      = itr
           if  (chtrname(itr).eq. 'C2H4'  ) ic2h4     = itr
           if  (chtrname(itr).eq. 'OLT'   ) iolt      = itr
           if  (chtrname(itr).eq. 'OLI'   ) ioli      = itr
           if  (chtrname(itr).eq. 'ALK4'  ) ialk4     = itr
           if  (chtrname(itr).eq. 'ALK7'  ) ialk7     = itr
           if  (chtrname(itr).eq. 'CO'    ) ico       = itr
           if  (chtrname(itr).eq. 'HCHO'  ) ihcho     = itr
           if  (chtrname(itr).eq. 'ALD2'  ) iald2     = itr
           if  (chtrname(itr).eq. 'ETHE'  ) iethe     = itr
           if  (chtrname(itr).eq. 'C2H6'  ) ic2h6     = itr
           if  (chtrname(itr).eq. 'C3H8'  ) ic3h8     = itr
           if  (chtrname(itr).eq. 'ISOP'  ) iisop     = itr
           if  (chtrname(itr).eq. 'TOLUE' ) itolue    = itr
           if  (chtrname(itr).eq. 'XYL'   ) ixyl      = itr
           if  (chtrname(itr).eq. 'NH3'   ) inh3      = itr
           if  (chtrname(itr).eq. 'PAN'   ) ipan      = itr
           if  (chtrname(itr).eq. 'ROOH'  ) irooh     = itr
           if  (chtrname(itr).eq. 'ACET'  ) iacet     = itr
           if  (chtrname(itr).eq. 'BENZ'  ) ibenz     = itr
           if  (chtrname(itr).eq. 'CH4'   ) ich4      = itr
           if  (chtrname(itr).eq. 'MOH'   ) imoh      = itr
           if  (chtrname(itr).eq. 'EOH'   ) ieoh      = itr
           if  (chtrname(itr).eq. 'ACO2'  ) iaco2     = itr
           if  (chtrname(itr).eq. 'CO2'   ) ico2      = itr
           if  (chtrname(itr).eq. 'DMS'   ) idms      = itr
           if  (chtrname(itr).eq. 'NOX'   ) inox      = itr
           if  (chtrname(itr).eq. 'HOX'   ) ihox      = itr
           if  (chtrname(itr).eq. 'SOX'   ) isox      = itr
           if  (chtrname(itr).eq. 'PAR'   ) ipar      = itr


        !abt *** Check to make sure SO4 is not defined twice as SULF or SO4 in          
        !    *** regcm.in namelist.  If both are defined then STOP                      
             if(iso4.ne.0 .and. isulf.ne.0) then
              write(*,*)"******* ERROR: Defined both SO4 and SULF"
              write(*,*)"*******        in regcm.in              "
              write(*,*)"*******        must choose one b/c they "
              write(*,*)"*******        both represent Sulfate   "
              stop
           end if


        !abt *** Added below to determine which MEGAN biogenic emission species         
        !    *** will be passed to the gas phase mechanism                              
        !    *** commented out lines correspond to species not advected but potentially
        !    *** used in chemistry mechanism.  Uncomment to give potential to advect    

!!$
!!$#if (defined VOC)
!!$           if(ibvoc.eq.1) then
!!$             if  (chtrname(itr).eq. 'ISOP'  ) bvoc_trmask(itr) = 1
!!$             if  (chtrname(itr).eq. 'APIN'  ) bvoc_trmask(itr) = 7
!!$             if  (chtrname(itr).eq. 'LIMO'  ) bvoc_trmask(itr) = 4
!!$          end if
!!$#endif
!!$        !abt above added 

        end do


print*, 'After startchem', icarb, isslt,idust


!!$        !*** abt added for wet deposition scheme
!        if(.not.allocated(chevap)) allocate(chevap(iy,kz))
!        if(.not.allocated(checum)) allocate(checum(iy,kz))

  !*** Initialize record read counter for CH EMISSI (see mod_che_ncio.F90)
  recc = 0
  
  open( 26,file='TUVGRID2', status='old')
  open( 25,file='REACTION.DAT_CBMZ', status='old')  

  call chem_initial_state(ifrest)

 print*, 'aprese chem_initial'



! FAB Traiter le prbleme du restart apres
!  call regchem


 end subroutine start_chem
!----------------------------------------------------------------


    subroutine chem_initial_state(ifrest)

!!$#ifndef IBM
!!$    use mpi
!!$#else
!!$      include 'mpif.h'
!!$#endif


       use mod_che_indices

    implicit none

    logical , intent(in) :: ifrest
    integer    :: i,j,k,ierr
    integer    :: icbc_date

    real(8) , dimension(iy,jxp) :: psdot


! FAB : TEMPORARY CUT the netcdf reading part to test the core

!!$  call allocate_species_bc
!!$  call allocate_mod_chem_bdycon
!!$  call allocate_mod_chem_mppio(.false.)
!!$
!!$
!!$   if(myid .eq. 0) then
!!$ 
!!$    write(*,*)'START CHEM INITIALIZATION',idate1,globidate1,myid
!!$ 
!!$    ndate0 = idate1
!!$    ndate1 = ndate0
!!$
!!$
!!$    if (ndate0.eq.globidate1 .or.               &
!!$            (((ndate0/10000)*100+1)*100 .eq.    &
!!$            ((globidate1/10000)*100+1)*100 ) ) then
!!$            icbc_date = globidate1
!!$       else
!!$            icbc_date = ((ndate0/10000)*100+1)*100
!!$    end if        
!!$
!!$
!!$
!!$   end if
!!$
!!$
!!$
!!$!FAB compatibility with standard chem mode
!!$!( igasphase==0 does not open and read oxydant fields and initial /boundary conditions are set to zero) 
!!$
!!$     if (igaschem==0) then !*** abt added: only for simple sulfate chemistry
!!$
!!$      if (myid ==0 ) then
!!$
!!$     call open_oxcl(icbc_date)
!!$      
!!$     call read_oxcl(ndate0,                            &
!!$           ohc0_io,ho2c0_io,o3c0_io,no3c0_io, h2o2c0_io) 
!!$
!!$
!!$      do j = 1 , jx
!!$         do k = 1 , kz
!!$            do i = 1 , iy
!!$               savch_0(i,k       ,j)  = ohc0_io(i,k,j)
!!$               savch_0(i,kz    +k,j)  = ho2c0_io(i,k,j)
!!$               savch_0(i,kz*2  +k,j)  = o3c0_io(i,k,j)
!!$               savch_0(i,kz*3  +k,j)  = no3c0_io(i,k,j)
!!$               savch_0(i,kz*4  +k,j)  = h2o2c0_io(i,k,j)
!!$            end do
!!$         end do
!!$      end do
!!$
!!$      end if
!!$
!!$!       Start transmission of data to other processors
!!$! nb here use already existing savch* arrays (dimensionned to 25), even if 5 variables are relevant( maybe more in the future)
!!$      call mpi_scatter(savch_0,iy*kz*25*jxp,mpi_real8,         &
!!$                       savch0, iy*kz*25*jxp,mpi_real8,         &
!!$                       0,mpi_comm_world,ierr)
!!$
!!$
!!$      do j = 1 , jendl
!!$         do k = 1 , kz
!!$            do i = 1 , iy
!!$               ohc0(i,k,j)        = savch0(i,      k,j)
!!$               ho2c0(i,k,j)        = savch0(i,kz  +k,j)
!!$               o3c0(i,k,j)       = savch0(i,kz*2  +k,j)
!!$               no3c0(i,k,j)      = savch0(i,kz*3  +k,j)
!!$               h2o2c0(i,k,j)      = savch0(i,kz*4 +k,j)
!!$            end do
!!$         end do
!!$      end do
!!$           
!!$     end if !end of igaschem condition
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$      if (igaschem==1) then 
!!$       
!!$      if ( myid.eq.0 ) then
!!$
!!$             call open_chbc(icbc_date)
!!$
!!$             call read_chbc(ndate0,                           &
!!$              o3b0_io,nob0_io,no2b0_io,hno3b0_io,n2o5b0_io     &
!!$             ,h2o2b0_io,ch4b0_io,cob0_io,hchob0_io,ch3ohb0_io  &
!!$             ,c2h5ohb0_io,c2h4b0_io,c2h6b0_io,ch3chob0_io      &
!!$             ,c3h6b0_io,c3h8b0_io,ch3coch3b0_io,bigeneb0_io    &
!!$             ,bigalkb0_io,isopb0_io,tolueb0_io,panb0_io        &
!!$             ,so2b0_io,dmsb0_io,so4b0_io)
!!$
!!$      do j = 1 , jx
!!$         do k = 1 , kz
!!$            do i = 1 , iy
!!$               savch_0(i,k       ,j)  = o3b0_io(i,k,j)
!!$               savch_0(i,kz    +k,j)  = nob0_io(i,k,j)
!!$               savch_0(i,kz*2  +k,j)  = no2b0_io(i,k,j)
!!$               savch_0(i,kz*3  +k,j)  = hno3b0_io(i,k,j)
!!$               savch_0(i,kz*4  +k,j)  = n2o5b0_io(i,k,j)
!!$               savch_0(i,kz*5  +k,j)  = h2o2b0_io(i,k,j)
!!$               savch_0(i,kz*6  +k,j)  = ch4b0_io(i,k,j)
!!$               savch_0(i,kz*7  +k,j)  = cob0_io(i,k,j)
!!$               savch_0(i,kz*8  +k,j)  = hchob0_io(i,k,j)
!!$               savch_0(i,kz*9  +k,j)  = ch3ohb0_io(i,k,j)
!!$               savch_0(i,kz*10  +k,j) = c2h5ohb0_io(i,k,j)
!!$               savch_0(i,kz*11  +k,j) = c2h4b0_io(i,k,j)
!!$               savch_0(i,kz*12  +k,j) = c2h6b0_io(i,k,j)
!!$               savch_0(i,kz*13  +k,j) = ch3chob0_io(i,k,j)
!!$               savch_0(i,kz*14  +k,j) = c3h6b0_io(i,k,j)
!!$               savch_0(i,kz*15  +k,j) = c3h8b0_io(i,k,j)
!!$               savch_0(i,kz*16  +k,j) = ch3coch3b0_io(i,k,j)
!!$               savch_0(i,kz*17  +k,j) = bigeneb0_io(i,k,j)
!!$               savch_0(i,kz*18  +k,j) = bigalkb0_io(i,k,j)
!!$               savch_0(i,kz*19  +k,j) = isopb0_io(i,k,j)
!!$               savch_0(i,kz*20  +k,j) = tolueb0_io(i,k,j)
!!$               savch_0(i,kz*21  +k,j) = panb0_io(i,k,j)
!!$               savch_0(i,kz*22  +k,j) = so2b0_io(i,k,j)
!!$               savch_0(i,kz*23  +k,j) = dmsb0_io(i,k,j)
!!$               savch_0(i,kz*24  +k,j) = so4b0_io(i,k,j)
!!$
!!$            end do
!!$         end do
!!$      end do
!!$
!!$     end if
!!$
!!$
!!$!       Start transmission of data to other processors
!!$      call mpi_scatter(savch_0,iy*kz*25*jxp,mpi_real8,         &
!!$                       savch0, iy*kz*25*jxp,mpi_real8,         &
!!$                       0,mpi_comm_world,ierr)
!!$
!!$
!!$      do j = 1 , jendl
!!$         do k = 1 , kz
!!$            do i = 1 , iy
!!$               o3b0(i,k,j)        = savch0(i,       k,j)
!!$               nob0(i,k,j)        = savch0(i,kz    +k,j)
!!$               no2b0(i,k,j)       = savch0(i,kz*2  +k,j)
!!$               hno3b0(i,k,j)      = savch0(i,kz*3  +k,j)
!!$               n2o5b0(i,k,j)      = savch0(i,kz*4  +k,j)
!!$               h2o2b0(i,k,j)      = savch0(i,kz*5  +k,j)
!!$               ch4b0(i,k,j)       = savch0(i,kz*6  +k,j)
!!$               cob0(i,k,j)        = savch0(i,kz*7  +k,j)
!!$               hchob0(i,k,j)      = savch0(i,kz*8  +k,j)
!!$               ch3ohb0(i,k,j)     = savch0(i,kz*9  +k,j)
!!$               c2h5ohb0(i,k,j)    = savch0(i,kz*10  +k,j)
!!$               c2h4b0(i,k,j)      = savch0(i,kz*11  +k,j)
!!$               c2h6b0(i,k,j)      = savch0(i,kz*12  +k,j)
!!$               ch3chob0(i,k,j)    = savch0(i,kz*13  +k,j)
!!$               c3h6b0(i,k,j)      = savch0(i,kz*14  +k,j)
!!$               c3h8b0(i,k,j)      = savch0(i,kz*15  +k,j)
!!$               ch3coch3b0(i,k,j)  = savch0(i,kz*16  +k,j)
!!$               bigeneb0(i,k,j)    = savch0(i,kz*17  +k,j)
!!$               bigalkb0(i,k,j)    = savch0(i,kz*18  +k,j)
!!$               isopb0(i,k,j)      = savch0(i,kz*19  +k,j)
!!$               tolueb0(i,k,j)     = savch0(i,kz*20  +k,j)
!!$               panb0(i,k,j)       = savch0(i,kz*21  +k,j)
!!$               so2b0(i,k,j)       = savch0(i,kz*22  +k,j)
!!$               dmsb0(i,k,j)       = savch0(i,kz*23  +k,j)
!!$               so4b0(i,k,j)       = savch0(i,kz*24  +k,j)
!!$            end do
!!$         end do
!!$      end do
!!$
!!$!
!!$!=======================================================================
!!$!       Couple pressure u,v,t,q
!!$!=======================================================================
!!$!
!!$!       this routine determines p(.) from p(x) by a 4-point
!!$!       interpolation. on the x-grid, a p(x) point outside the grid
!!$!       domain is assumed to satisfy p(0,j)=p(1,j); p(iy,j)=p(iym1,j);
!!$!       and similarly for the i's.
!!$
!!$      do j = jbegin , jendx
!!$         do i = 2 , iym1
!!$            psdot(i,j) = 0.25*(ps0(i,j)   + ps0(i-1,j) +      &
!!$                         ps0(i,j-1) + ps0(i-1,j-1))
!!$         end do
!!$      end do
!!$      do j = jbegin , jendx
!!$         psdot(1,j) = 0.5*(ps0(1,j)+ps0(1,j-1))
!!$         psdot(iy,j) = 0.5*(ps0(iym1,j)+ps0(iym1,j-1))
!!$      end do
!!$
!!$
!!$! intialise the chib0 tracer : correspondance betwwen bdy tables and chib
!!$
!!$      do k = 1 , kz
!!$         do j = 1 , jendl
!!$            do i = 1 , iy
!!$          
!!$                if(io3 > 0)    chib0(i,k,j,io3)     =  o3b0(i,k,j)*psdot(i,j)
!!$                if(ino2 > 0)   chib0(i,k,j,ino2)    =  no2b0(i,k,j)*psdot(i,j)
!!$                if(ino > 0)    chib0(i,k,j,ino)     =  nob0(i,k,j)*psdot(i,j)
!!$                if(ico > 0)    chib0(i,k,j,ico)     =  cob0(i,k,j)*psdot(i,j)
!!$                if(ih2o2 > 0)  chib0(i,k,j,ih2o2)   =  h2o2b0(i,k,j)*psdot(i,j)
!!$                if(ihno3 > 0)  chib0(i,k,j,ihno3)   =  hno3b0(i,k,j)*psdot(i,j)
!!$                if(in2o5 > 0)  chib0(i,k,j,in2o5)   =  n2o5b0(i,k,j)*psdot(i,j)
!!$                if(iso2 > 0)   chib0(i,k,j,iso2)    =  so2b0(i,k,j)*psdot(i,j)
!!$                if(iso4 > 0)   chib0(i,k,j,iso4)    =  so4b0(i,k,j)*psdot(i,j)
!!$                if(idms > 0)   chib0(i,k,j,idms)    =  dmsb0(i,k,j)*psdot(i,j)
!!$                if(ihcho > 0)  chib0(i,k,j,ihcho)   =  hchob0(i,k,j)*psdot(i,j)
!!$                if(iald2 > 0)  chib0(i,k,j,iald2)   =  ch3chob0(i,k,j)*psdot(i,j)
!!$                if(iisop > 0)  chib0(i,k,j,iisop)   =  isopb0(i,k,j)*psdot(i,j)
!!$                if(ic2h6 > 0)  chib0(i,k,j,ic2h6)   =  c2h6b0(i,k,j)*psdot(i,j)
!!$                if(ic3h8 > 0)  chib0(i,k,j,ic3h8)   =  c3h8b0(i,k,j)*psdot(i,j)
!!$                if(ialk4 > 0)  chib0(i,k,j,ialk4)   =  bigalkb0(i,k,j)*psdot(i,j)
!!$                if(ialk7 > 0)  chib0(i,k,j,ialk7)   =  bigalkb0(i,k,j)*psdot(i,j)
!!$                if(iacet > 0)  chib0(i,k,j,iacet)   =  ch3coch3b0(i,k,j)*psdot(i,j)
!!$                if(imoh > 0)   chib0(i,k,j,imoh)    =  ch3ohb0(i,k,j)*psdot(i,j)
!!$                if(ieoh > 0)   chib0(i,k,j,ieoh)    =  c2h5ohb0(i,k,j)*psdot(i,j)
!!$                if(ioli > 0)   chib0(i,k,j,ioli)    =  bigeneb0(i,k,j)*psdot(i,j)
!!$                if(iolt > 0)   chib0(i,k,j,iolt)    =  c3h6b0(i,k,j)*psdot(i,j)
!!$                if(itolue > 0) chib0(i,k,j,itolue)  =  tolueb0(i,k,j)*psdot(i,j)
!!$                if(ich4 > 0)   chib0(i,k,j,ich4)    =  ch4b0(i,k,j)*psdot(i,j)
!!$                if(ipan > 0)   chib0(i,k,j,ipan)    =  panb0(i,k,j)*psdot(i,j)
!!$                if(iethe > 0)  chib0(i,k,j,iethe)   =  c2h4b0(i,k,j)*psdot(i,j)
!!$            end do
!!$         end do
!!$      end do
!!$
!!$
!!$ end if !end test igaschem
!!$



! also initialise chia and chib at this stage
!if restart is not active 
! note if igaschem ==0 , aerosol BC are set to chib0 = 0  for now 


!!!!! FAB: temporay for the core testing 
chib0 = d_zero
chib1 = d_zero
!!!!!!!!

if (.not.ifrest ) then
         do j=1,jendl
         do i=1,iy
         do k=1,kz
              chia(i,k,j,:) =   chib0(i,k,j,:) 
              chib(i,k,j,:) =   chib0(i,k,j,:)
         end do
         end do
         end do

end if






   end subroutine chem_initial_state




  end module mod_che_start
