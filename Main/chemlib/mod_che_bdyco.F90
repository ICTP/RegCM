!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or
!    modify
!    it under the terms of the GNU General Public License as
!    published by
!    the Free Software Foundation, either version 3 of the
!    License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty
!    of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public
!    License
!    along with ICTP RegCM.  If not, see
!    <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_che_bdyco

use mod_realkinds
  use mod_dynparam
use mod_memutil

  use mod_service
  use mod_mpmessage
  use mod_che_common
  use mod_che_mppio
  use mod_che_ncio
  use mod_che_species

  private

  public :: allocate_mod_che_bdyco , chem_bdyin , chem_bdyval

#ifndef BAND
  public :: chieb, chiebt , chiwb , chiwbt
#endif
  public :: chisb , chisbt , chinb , chinbt
  public :: chib0 , chib1,ichbdy2trac

  real(dp) ,  pointer, dimension(:,:,:,:) :: chib0 , chib1
!
#ifndef BAND
  real(dp) ,  pointer, dimension(:,:,:,:) :: chieb , chiebt , &
                                                 chiwb , chiwbt
#endif
  real(dp) , pointer, dimension(:,:,:,:) :: chinb , chinbt , &
                                                 chisb , chisbt
 
  integer, pointer, dimension(:) ::ichbdy2trac
   

  contains

!
  subroutine allocate_mod_che_bdyco
    implicit none

    character (len=64) :: subroutine_name='allocate_mod_che_bdyco'
    integer :: idindx=0
!
!    call time_begin(subroutine_name,idindx)

    call getmem4d(chib0,1,iy,1,kz,1,jxp,1,ntr,'mod_che_bdyco:chib0')
    call getmem4d(chib1,1,iy,1,kz,1,jxp,1,ntr,'mod_che_bdyco:chib1')
#ifndef BAND
    call getmem4d(chieb,1,iy,1,kz,0,jxp+1,1,ntr,'mod_che_bdyco:chieb')
    call getmem4d(chiebt,1,iy,1,kz,0,jxp+1,1,ntr,'mod_che_bdyco:chiebt')
    call getmem4d(chiwb,1,iy,1,kz,0,jxp+1,1,ntr,'mod_che_bdyco:chiwb')
    call getmem4d(chiwbt,1,iy,1,kz,0,jxp+1,1,ntr,'mod_che_bdyco:chiwbt')
#endif
    call getmem4d(chinb,1,nspgx,1,kz,0,jxp+1,1,ntr,'mod_che_bdyco:chinb')
    call getmem4d(chinbt,1,nspgx,1,kz,0,jxp+1,1,ntr,'mod_che_bdyco:chinbt')
    call getmem4d(chisb,1,nspgx,1,kz,0,jxp+1,1,ntr,'mod_che_bdyco:chisb')
    call getmem4d(chisbt,1,nspgx,1,kz,0,jxp+1,1,ntr,'mod_che_bdyco:chinst')

    call getmem1d(ichbdy2trac,1,25,'mod_che_bdyco:ichbdytrac')

!    call time_end(subroutine_name,idindx)

  end subroutine allocate_mod_che_bdyco
!!$!
!!$!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$!                                                                     c
!!$!     this subroutine reads in the boundary conditions.               c
!!$!                                                                     c
!!$!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine chem_bdyin (dtbdys, bdydate1 , bdydate2)

    use mod_che_indices

#ifndef IBM
    use mpi
#endif
#ifdef IBM
    include 'mpif.h'
#endif
implicit none

    type(rcm_time_and_date) ,intent(in) :: bdydate1 , bdydate2

    real(dp), intent(in) :: dtbdys
    integer :: i , j , k ,n, nn , nnb , mmrec, itr
    integer ::  nxeb , nxwb, ierr
#ifndef BAND
    integer :: nkk
#endif
    real(8) , dimension(iy,jxp) :: psdot , tdum
    character (len=64) :: subroutine_name='chem_bdyin'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
!!$    if ( dabs(xtime).gt.0.0001 ) return
!!$!
 

! FAB : le bloc de lecture est commente pour l'instant
!!$
!!$    if ( myid .eq. 0 .and. igaschem==0)then
!!$
!!$    mmrec = oxcl_search(ndate1)
!!$
!!$    print*, 'mmrec',mmrec
!!$      if (mmrec < 0) then
!!$     call open_oxcl(imonfirst(ndate1))
!!$      print*,'after open_oxcl'
!!$      end if
!!$
!!$
!!$       call read_oxcl(ndate1,ohc0_io,ho2c0_io,o3c0_io,no3c0_io, h2o2c0_io)
!!$
!!$    do j = 1 , jx
!!$       do k = 1 , kz
!!$          do i = 1 , iy
!!$             savch_0(i,k       ,j)  = ohc0_io(i,k,j)
!!$             savch_0(i,kz    +k,j)  = ho2c0_io(i,k,j)
!!$             savch_0(i,kz*2  +k,j)  = o3c0_io(i,k,j)
!!$             savch_0(i,kz*3  +k,j)  = no3c0_io(i,k,j)
!!$             savch_0(i,kz*4  +k,j)  = h2o2c0_io(i,k,j)
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    print*, 'OXCLIM READING', maxval(h2o2c0_io), minval(h2o2c0_io)
!!$
!!$    end if
!!$
!!$!       Start transmission of data to other processors
!!$! nb here use already existing savch* arrays (dimensionned to 25), even if 5 variables are relevant( maybe more in the future)
!!$    call mpi_scatter(savch_0,iy*kz*25*jxp,mpi_real8,         &
!!$                     savch0, iy*kz*25*jxp,mpi_real8,         &
!!$                     0,mpi_comm_world,ierr)
!!$
!!$
!!$    do j = 1 , jendl
!!$       do k = 1 , kz
!!$          do i = 1 , iy
!!$             ohc0(i,k,j)       = savch0(i,      k,j)
!!$             ho2c0(i,k,j)      = savch0(i,kz  +k,j)
!!$             o3c0(i,k,j)       = savch0(i,kz*2  +k,j)
!!$             no3c0(i,k,j)      = savch0(i,kz*3  +k,j)
!!$             h2o2c0(i,k,j)     = savch0(i,kz*4 +k,j)
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    print*, 'after', myid, maxval(h2o2c0), minval(h2o2c0)
!!$
!!$
!!$savch0  = 0.
!!$savch_0 =0.
!!$
!!$
!!$!  call allocate_species_bc
!!$!  call allocate_mod_chem_mppio(.false.)
!!$
      if(myid .eq. 0 ) then

    mmrec = chbc_search(bdydate2)
    write(*,*)'CH_BDYIN -----',mmrec,bdydate2
    if (mmrec < 0) then
            call open_chbc(bdydate2)
    end if
    
     call read_chbc(chebdy_io)  

   do n = 1, 25
      do j = 1 , jx
        do k = 1 , kz
           do i = 1 , iy
             savch_0(i,kz*(n-1) + k       ,j)  = chebdy_io(i,k,j,n)
!!$             savch_0(i,k       ,j)  = o3b1_io(i,k,j)
!!$             savch_0(i,kz    +k,j)  = nob1_io(i,k,j)
!!$             savch_0(i,kz*2  +k,j)  = no2b1_io(i,k,j)
!!$             savch_0(i,kz*3  +k,j)  = hno3b1_io(i,k,j)
!!$             savch_0(i,kz*4  +k,j)  = n2o5b1_io(i,k,j)
!!$             savch_0(i,kz*5  +k,j)  = h2o2b1_io(i,k,j)
!!$             savch_0(i,kz*6  +k,j)  = ch4b1_io(i,k,j)
!!$             savch_0(i,kz*7  +k,j)  = cob1_io(i,k,j)
!!$             savch_0(i,kz*8  +k,j)  = hchob1_io(i,k,j)
!!$             savch_0(i,kz*9  +k,j)  = ch3ohb1_io(i,k,j)
!!$             savch_0(i,kz*10  +k,j) = c2h5ohb1_io(i,k,j)
!!$             savch_0(i,kz*11  +k,j) = c2h4b1_io(i,k,j)
!!$             savch_0(i,kz*12  +k,j) = c2h6b1_io(i,k,j)
!!$             savch_0(i,kz*13  +k,j) = ch3chob1_io(i,k,j)
!!$             savch_0(i,kz*14  +k,j) = c3h6b1_io(i,k,j)
!!$             savch_0(i,kz*15  +k,j) = c3h8b1_io(i,k,j)
!!$             savch_0(i,kz*16  +k,j) = ch3coch3b1_io(i,k,j)
!!$             savch_0(i,kz*17  +k,j) = bigeneb1_io(i,k,j)
!!$             savch_0(i,kz*18  +k,j) = bigalkb1_io(i,k,j)
!!$             savch_0(i,kz*19  +k,j) = isopb1_io(i,k,j)
!!$             savch_0(i,kz*20  +k,j) = tolueb1_io(i,k,j)
!!$             savch_0(i,kz*21  +k,j) = panb1_io(i,k,j)
!!$             savch_0(i,kz*22  +k,j) = so2b1_io(i,k,j)
!!$             savch_0(i,kz*23  +k,j) = dmsb1_io(i,k,j)
!!$             savch_0(i,kz*24  +k,j) = so4b1_io(i,k,j)
          end do
       end do
    end do
    end do
!!$
     end if
    call mpi_scatter(savch_0,iy*kz*25*jxp,mpi_real8,      &
                     savch0, iy*kz*25*jxp,mpi_real8,      &
                     0,mpi_comm_world,ierr)
!!$
!!$       print*,' CIAO ,',myid, size(no2b1,3), size(savch0,3),jendl,jbegin,jendx,jxp
!!$

    do n=1,25 
    do j = 1 , jendl
    do k = 1 , kz
      do i = 1 , iy

        
             chebdy(i,k,j,n)        = savch0(i, kz*(n-1) + k,j)
!!$             nob1(i,k,j)        = savch0(i,kz    +k,j)
!!$             no2b1(i,k,j)       = savch0(i,kz*2  +k,j)
!!$             hno3b1(i,k,j)      = savch0(i,kz*3  +k,j)
!!$             n2o5b1(i,k,j)      = savch0(i,kz*4  +k,j)
!!$             h2o2b1(i,k,j)      = savch0(i,kz*5  +k,j)
!!$             ch4b1(i,k,j)       = savch0(i,kz*6  +k,j)
!!$             cob1(i,k,j)        = savch0(i,kz*7  +k,j)
!!$             hchob1(i,k,j)      = savch0(i,kz*8  +k,j)
!!$             ch3ohb1(i,k,j)     = savch0(i,kz*9  +k,j)
!!$             c2h5ohb1(i,k,j)    = savch0(i,kz*10  +k,j)
!!$             c2h4b1(i,k,j)      = savch0(i,kz*11  +k,j)
!!$             c2h6b1(i,k,j)      = savch0(i,kz*12  +k,j)
!!$             ch3chob1(i,k,j)    = savch0(i,kz*13  +k,j)
!!$             c3h6b1(i,k,j)      = savch0(i,kz*14  +k,j)
!!$             c3h8b1(i,k,j)      = savch0(i,kz*15  +k,j)
!!$             ch3coch3b1(i,k,j)  = savch0(i,kz*16  +k,j)
!!$             bigeneb1(i,k,j)    = savch0(i,kz*17  +k,j)
!!$             bigalkb1(i,k,j)    = savch0(i,kz*18  +k,j)
!!$             isopb1(i,k,j)      = savch0(i,kz*19  +k,j)
!!$             tolueb1(i,k,j)     = savch0(i,kz*20  +k,j)
!!$             panb1(i,k,j)       = savch0(i,kz*21  +k,j)
!!$             so2b1(i,k,j)       = savch0(i,kz*22  +k,j)
!!$             dmsb1(i,k,j)       = savch0(i,kz*23  +k,j)
!!$             so4b1(i,k,j)       = savch0(i,kz*24  +k,j)
          end do
       end do
    end do

     print*, 'apres transfert', myid, n,maxval(chebdy(:,:,:,n)) 

    end do


!!$    do j = jbegin , jendx
!!$       do i = 2 , iym1
!!$           psdot(i,j) = 0.25*(ps1(i,j)   + ps1(i-1,j) +  &
!!$                        ps1(i,j-1) + ps1(i-1,j-1))
!!$       end do
!!$    end do
!!$    do i = 2 , iym1
!!$       if ( myid.eq.0 ) psdot(i,1) = 0.5*(ps1(i,1)+ps1(i-1,1))
!!$       if ( myid.eq.nproc-1 )                             &
!!$           psdot(i,jendl) = 0.5*(ps1(i,jendx)+ps1(i-1,jendx))
!!$    end do
!!$    do j = jbegin , jendx
!!$       psdot(1,j) = 0.5*(ps1(1,j)+ps1(1,j-1))
!!$       psdot(iy,j) = 0.5*(ps1(iym1,j)+ps1(iym1,j-1))
!!$    end do
!!$
!!$    if ( myid.eq.0 ) then
!!$         psdot(1,1) = ps1(1,1)
!!$         psdot(iy,1) = ps1(iym1,1)
!!$    end if
!!$
!!$    if ( myid.eq.nproc-1 ) then
!!$       psdot(1,jendl) = ps1(1,jendx)
!!$       psdot(iy,jendl) = ps1(iym1,jendx)
!!$    end if
!!$    write(*,*)' BDYIN---ZZZ--',maxval(o3b1(:,1,:)),maxval(psdot(:,:)),myid
!!$!=======================================================================
!!$!       Couple pressure u,v,t,q
!!$
!!$! fab modif
!!$! intialise the chib1 tracer : correspondance betwwen bdy tables and chib
!!$
   do n=1,25
    do k = 1 , kz
       do j = 1 , jendl
          do i = 1 , iy
              if(ichbdy2trac(n) > 0) chib1(i,k,j,ichbdy2trac(n)) = chebdy(i,k,j,n)*cpsb(j,i)
          end do
       end do
    end do
   end do




!-----compute boundary conditions for p*chi
! first define nxeb and nxwb
!

!!$
#ifdef BAND
    nxwb=0
    nxeb=0
#else
    if ( nspgx.le.jxp ) then
      nxwb = nspgx
    else
      nkk = nspgx/jxp
      if ( nspgx.eq.nkk*jxp ) then
        nxwb = jxp
      else
        nxwb = nspgx - nkk*jxp
      end if
    end if
    if ( nxwb+myid*jxp.gt.nspgx ) then
      nxwb = 0
    else if ( nxwb+myid*jxp.lt.nspgx ) then
      nxwb = jxp
    end if
!
    if ( nspgx.le.jxp-1 ) then
      nxeb = nspgx
    else
      nkk = (nspgx-jxp+1)/jxp
      if ( (nspgx-jxp+1).eq.nkk*jxp ) then
        nxeb = jxp
      else
        nxeb = (nspgx-jxp+1) - nkk*jxp
      end if
    end if
    if ( jxm1-(myid*jxp+jxp-nxeb).gt.nspgx ) then
      nxeb = 0
    else if ( jxm1-(myid*jxp+jxp-nxeb).lt.nspgx ) then
      nxeb = min(jendx,jxp)
    end if

#endif

!
#ifndef BAND
    do nn = 1 , nxwb
      do k = 1 , kz
        do i = 1 , iym1

          chiwb(i,k,nn,:) = chib0(i,k,nn,:)
          chiwbt(i,k,nn,:) = (chib1(i,k,nn,:)-chib0(i,k,nn,:))/dtbdys

        end do
      end do
    end do
    do nn = 1 , nxeb
      nnb = min(jendx,jxp) - nn + 1
      do k = 1 , kz
        do i = 1 , iym1
          chieb(i,k,nn,:) = chib0(i,k,nnb,:)
          chiebt(i,k,nn,:) = (chib1(i,k,nnb,:)-chib0(i,k,nnb,:))/dtbdys

        end do
      end do
    end do
#endif
    do nn = 1 , nspgx
      nnb = iym1 - nn + 1
      do k = 1 , kz
        do j = 1 , jendx
          chinb(nn,k,j,:) = chib0(nnb,k,j,:)
          chisb(nn,k,j,:) = chib0(nn,k,j,:)
          chinbt(nn,k,j,:) = (chib1(nnb,k,j,:)-chib0(nnb,k,j,:))/dtbdys
          chisbt(nn,k,j,:) = (chib1(nn,k,j,:)-chib0(nn,k,j,:))/dtbdys
        end do
      end do
    end do
! prepare for next bdy step
    do k = 1 , kz
      do j = 1 , jendl
        do i = 1 , iy
          chib0(i,k,j,:) = chib1(i,k,j,:)
        end do
      end do
    end do
  end subroutine chem_bdyin



  subroutine chem_bdyval(xt,iexec,nbdytime,dtbdys,ktau,ifrest)
    use mod_che_indices
!
    implicit none
!
    integer :: iexec
    integer(8)::nbdytime, ktau
    real(8) :: xt, dtbdys
    logical :: ifrest
    intent (in) iexec , xt, nbdytime,dtbdys,ktau,ifrest
!
    real(8) :: chix , chix1 , chix2 , dtb , vavg
    integer :: itr , j , k
#ifndef BAND
    integer :: i
    real(8) :: uavg
#endif
    character (len=50) :: subroutine_name='chem_bdyval'
    integer :: idindx=0
!
!    call time_begin(subroutine_name,idindx)
!
!*********************************************************************
!*****fill up the boundary value for xxb variables from xxa variables:
!     if this subroutine is called for the first time, this part
!     shall be skipped.
!
    if ( iexec.ne.1 ) then
!
!      if ( ichem.eq.1 ) then
!-----for p*chi (tracers)
        do itr = 1 , ntr
          do k = 1 , kz
#ifndef BAND
            do i = 1 , iym1
              if ( myid.eq.0 ) chib(i,k,1,itr) = chia(i,k,1,itr)
              if ( myid.eq.nproc-1 ) chib(i,k,jendx,itr)              &
                 & = chia(i,k,jendx,itr)
            end do
#endif
            do j = jbegin , jendm
              chib(1,k,j,itr) = chia(1,k,j,itr)
              chib(iym1,k,j,itr) = chia(iym1,k,j,itr)
            end do
          end do
        end do
!      end if
!chem2_
!
!
    end if      !end if(iexec.ne.1) test
!**********************************************************************
!*****compute the boundary values for xxa variables:
!
!-----compute the time interval for boundary tendency:
!
    dtb = xt
    if ( nbdytime == 0 .and. ktau > 0 ) then
      dtb = dtbdys
    end if

!
!
!-----set boundary values for p*t:
!-----set boundary values for p*qv:
!

    if (iexec.eq.1 .and. ifrest)  return


!    if ( iboudy.eq.0 ) then
!.....fixed boundary conditions:
!    end if
!


!.....time-dependent boundary conditions:
! for chemistry relaxation towrds
! time dependant boundary conditions is considered

    do k = 1 , kz

#ifndef BAND
      do i = 1 , iym1
        if ( myid.eq.0 ) chia(i,k,1,:) = chiwb(i,k,1,:) + dtb*chiwbt(i,k,1,:)
        if ( myid.eq.nproc-1 ) chia(i,k,jendx,:) = chieb(i,k,1,:)        &
           & + dtb*chiebt(i,k,1,:)
      end do
#endif
      do j = jbegin , jendm
        chia(1,k,j,:) = chisb(1,k,j,:) + dtb*chisbt(1,k,j,:)
        chia(iym1,k,j,:) = chinb(1,k,j,:) + dtb*chinbt(1,k,j,:)
      end do
    end do

!!$    if ( iboudy.eq.3 .or. iboudy.eq.4 ) then
!!$!
!!$!-----determine boundary values depends on inflow/outflow:0
!!$!overwrite the dfault case : keep it maybe for certainspecific  tracers in the future
!!$
!!$!    if ( ichem.eq.1 ) then
!!$!chem2
!!$
!!$!----add tracer bc's
!!$!
!!$      do itr = 1 , ntr
!!$        do k = 1 , kz
!!$#ifndef BAND
!!$!
!!$!.....west  boundary:
!!$!
!!$          if ( myid.eq.0 ) then
!!$            do i = 1 , iym1
!!$              chix1 = chia(i,k,1,itr)/sps1%ps(i,1)
!!$              chix2 = chia(i,k,2,itr)/sps1%ps(i,2)
!!$              uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
!!$              if ( uavg.ge.0. ) then
!!$                chix = chix1
!!$              else
!!$                chix = chix2
!!$              end if
!!$              chia(i,k,1,itr) = chix*sps1%ps(i,1)
!!$            end do
!!$          end if
!!$!
!!$!.....east  boundary:
!!$!
!!$          if ( myid.eq.nproc-1 ) then
!!$            do i = 1 , iym1
!!$              chix1 = chia(i,k,jendx,itr)/sps1%ps(i,jendx)
!!$              chix2 = chia(i,k,jendm,itr)/sps1%ps(i,jendm)
!!$              uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
!!$              if ( uavg.lt.0. ) then
!!$                chix = chix1
!!$
!!$              else
!!$                chix = chix2
!!$              end if
!!$              chia(i,k,jendx,itr) = chix*sps1%ps(i,jendx)
!!$            end do
!!$          end if
!!$#endif
!!$!
!!$!.....south boundary:
!!$!
!!$          do j = jbegin , jendm
!!$            chix1 = chia(1,k,j,itr)/sps1%ps(1,j)
!!$            chix2 = chia(2,k,j,itr)/sps1%ps(2,j)
!!$            vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
!!$            if ( vavg.ge.0. ) then
!!$              chix = chix1
!!$! fab test
!!$            else
!!$              chix = chix2
!!$            end if
!!$            chia(1,k,j,itr) = chix*sps1%ps(1,j)
!!$          end do
!!$!
!!$!.....north boundary:
!!$!
!!$          do j = jbegin , jendm
!!$            chix1 = chia(iym1,k,j,itr)/sps1%ps(iym1,j)
!!$            chix2 = chia(iym2,k,j,itr)/sps1%ps(iym2,j)
!!$            vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
!!$            if ( vavg.lt.0. ) then
!!$              chix = chix1
!!$            else
!!$              chix = chix2
!!$            end if
!!$            chia(iym1,k,j,itr) = chix*sps1%ps(iym1,j)
!!$          end do
!!$        end do
!!$      end do
!!$!chem2_
!!$!   end if

!!$    end if      !end if(iboudy.eq.3.or.4) test
!!$!
!    call time_end(subroutine_name,idindx)


  end subroutine chem_bdyval
!!$!

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                     c
  !     these subroutines apply relaxation boundary conditions to the   c
  !     tendency term - xpten.                                          c
  !                                                                     c
  !     ip    : is the number of slices affected by nudging.            c
  !                                                                     c
  !     xt    : is the time in minutes for variable "psb".              c
  !                                                                     c
  !     fcoef : are the coefficients for the newtonian term.            c
  !                                                                     c
  !     gcoef : are the coefficients for the diffusion term.            c
  !                                                                     c
  !     xpten : is the tendency calculated from the model.              c
  !                                                                     c
  !     peb, pwb, pss, pnb : are the observed boundary values           c
  !                   on east, west, south, and north boundaries.       c
  !                                                                     c
  !     pebt, pwbt, psbt, pnbt : are the large-scale or observed        c
  !             tendencies at east, west, south, and north boundaries.  c
  !                                                                     c
  !     psb    : is the variable at tau-1.                               c
  !                                                                     c
  !     ie = iy, je = jx for dot-point variables.                       c
  !     ie = iym1, je = jxm1 for cross-point variables.                 c
  !                                                                     c
  !     j    : is the j'th slice of the tendency to be adjusted.        c
  !     ibdy : type of boundary condition relaxation, 1=linear        c
  !              5 = exponential                                        c
  !                                                                     c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine nudge_chi(ip,fcoef,gcoef,xt,ften,j,ibdy)
    !
    implicit none
    !
    real(8) :: fcoef , gcoef , xt
    integer :: ibdy , ip , j
    real(8) , dimension(iy,kz,ntr) :: ften
    intent (in) fcoef , gcoef , ibdy , ip , j , xt
    intent (inout) ften
    !
    real(8) :: dtb , fcx ,gcx
    real(8), dimension(ntr)::    fls0 , fls1 , fls2 , fls3 , fls4
    integer :: i , ii , k

#ifndef BAND
    integer :: ibeg , iend , jj , jsls, jwb , jeb
#endif


    !
    cHARACTER (len=50) :: subroutine_name='nudge_chi'
    integer :: idindx=0
    !
!    call time_begin(subroutine_name,idindx)


#ifdef BAND
    !----------------------------------------------------------------------
    !
    !FAB: care here 
    dtb = xt
    !
    !-----determine which relaxation method to use:linear/expon.
    !
    if ( ibdy.eq.1 ) then
       !
    print*,'LINEAR NUDGING NOT IMPLEMENTED for CHEM fields !!'

    else if ( ibdy.eq.5 ) then

       !----------use exponential method
       !------interior j slices:
       do i = 2 , ip
          ii = iym1 - i + 1
          do k = 1 , kz
             fcx = fcoef*xfune(i,k)
             gcx = gcoef*xfune(i,k)
             !........south boundary:
             fls0(:) = (chisb(i,k,j,:)+dtb*chisbt(i,k,j,:)) -chib(i,k,j,:)
             fls1(:) = (chisb(i,k,j-1,:)+dtb*chisbt(i,k,j-1,:)) -chib(i,k,j-1,:)
             fls2(:) = (chisb(i,k,j+1,:)+dtb*chisbt(i,k,j+1,:)) -chib(i,k,j+1,:)
             fls3(:) = (chisb(i-1,k,j,:)+dtb*chisbt(i-1,k,j,:)) -chib(i-1,k,j,:)
             fls4(:) = (chisb(i+1,k,j,:)+dtb*chisbt(i+1,k,j,:)) -chib(i+1,k,j,:)
             ften(i,k,:) = ften(i,k,:) + fcx*fls0(:) -                        &
                  & gcx*crdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-4.*fls0(:))
             !........north boundary:
             fls0(:) = (chinb(i,k,j,:)+dtb*chinbt(i,k,j,:)) -chib(ii,k,j,:)
             fls1(:) = (chinb(i,k,j-1,:)+dtb*chinbt(i,k,j-1,:))-chib(ii,k,j-1,:)
             fls2(:) = (chinb(i,k,j+1,:)+dtb*chinbt(i,k,j+1,:))-chib(ii,k,j+1,:)
             fls3(:) = (chinb(i-1,k,j,:)+dtb*chinbt(i-1,k,j,:))-chib(ii-1,k,j,:)
             fls4(:) = (chinb(i+1,k,j,:)+dtb*chinbt(i+1,k,j,:))-chib(ii+1,k,j,:)
             ften(ii,k,:) = ften(ii,k,:) + fcx*fls0(:) -                      &
                  & gcx*crdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-4.*fls0(:))
          end do
       end do

    else
    end if

#else ! fin ifdef band
    !----------------------------------------------------------------------
    !
    dtb = xt

    jsls = j + myid*jxp
    jj = jx - jsls
    if ( jj.le.ip ) jsls = jj
    jwb = jsls
    if ( jwb.gt.jxp ) jwb = mod(jwb,jxp)
    if ( jwb.eq.0 ) jwb = jxp
    if ( myid.eq.nproc-1 ) then
       jeb = jsls
    else
       jeb = jsls + 1
    end if
    if ( jeb.gt.jxp ) jeb = mod(jeb,jxp)
    if ( jeb.eq.0 ) jeb = jxp

    !
    !-----determine which relaxation method to use:linear/expon.
    !
    if ( ibdy.eq.1 ) then
       !
     print*,'CHEM LINEAR NUDGING NOT IMPLEMENTED !!!'
       !
    else if ( ibdy.eq.5 ) then

       !----------use exponential method

       if ( jsls.gt.ip ) then
          !------interior j slices:
          do i = 2 , ip
             ii = iym1 - i + 1
             do k = 1 , kz
                fcx = fcoef*xfune(i,k)
                gcx = gcoef*xfune(i,k)
!!$                !........south boundary:
                fls0(:) = (chisb(i,k,j,:)+dtb*chisbt(i,k,j,:)) -chib(i,k,j,:)
                fls1(:) = (chisb(i,k,j-1,:)+dtb*chisbt(i,k,j-1,:)) -chib(i,k,j-1,:)
                fls2(:) = (chisb(i,k,j+1,:)+dtb*chisbt(i,k,j+1,:)) -chib(i,k,j+1,:)
                fls3(:) = (chisb(i-1,k,j,:)+dtb*chisbt(i-1,k,j,:)) -chib(i-1,k,j,:)
                fls4(:) = (chisb(i+1,k,j,:)+dtb*chisbt(i+1,k,j,:)) -chib(i+1,k,j,:)
                ften(i,k,:) = ften(i,k,:) + fcx*fls0(:) -                        &
                     & gcx*crdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-4.*fls0(:))
                !........north boundary:
                fls0(:) = (chinb(i,k,j,:)+dtb*chinbt(i,k,j,:)) -chib(ii,k,j,:)
                fls1(:) = (chinb(i,k,j-1,:)+dtb*chinbt(i,k,j-1,:))-chib(ii,k,j-1,:)
                fls2(:) = (chinb(i,k,j+1,:)+dtb*chinbt(i,k,j+1,:))-chib(ii,k,j+1,:)
                fls3(:) = (chinb(i-1,k,j,:)+dtb*chinbt(i-1,k,j,:))-chib(ii-1,k,j,:)
                fls4(:) = (chinb(i+1,k,j,:)+dtb*chinbt(i+1,k,j,:))-chib(ii+1,k,j,:)
                ften(ii,k,:) = ften(ii,k,:) + fcx*fls0(:) -                      &
                     & gcx*crdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-4.*fls0(:))

           end do
          end do
          !
       else if ( jsls.le.ip ) then
          !------east or west boundary slices:
          ibeg = 2
          iend = iym1 - 1
          if ( jsls.gt.2 ) then
               do i = 2 , jsls - 1
                ii = iym1 - i + 1
                do k = 1 , kz
                   fcx = fcoef*xfune(i,k)
                   gcx = gcoef*xfune(i,k)
                   !.........south boundary:
                   fls0(:) = (chisb(i,k,j,:)+dtb*chisbt(i,k,j,:)) -chib(i,k,j,:)
                   fls1(:) = (chisb(i,k,j-1,:)+dtb*chisbt(i,k,j-1,:))-chib(i,k,j-1,:)
                   fls2(:) = (chisb(i,k,j+1,:)+dtb*chisbt(i,k,j+1,:))-chib(i,k,j+1,:)
                   fls3(:) = (chisb(i-1,k,j,:)+dtb*chisbt(i-1,k,j,:))-chib(i-1,k,j,:)
                   fls4(:) = (chisb(i+1,k,j,:)+dtb*chisbt(i+1,k,j,:))-chib(i+1,k,j,:)
                   ften(i,k,:) = ften(i,k,:) + fcx*fls0(:) -                      &
                        & gcx*crdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-4.*fls0(:))
                   !.........north boundary:
                   fls0(:) = (chinb(i,k,j,:)+dtb*chinbt(i,k,j,:)) -chib(ii,k,j,:)
                   fls1(:) = (chinb(i,k,j-1,:)+dtb*chinbt(i,k,j-1,:)) - &
                        chib(ii,k,j-1,:)
                   fls2(:) = (chinb(i,k,j+1,:)+dtb*chinbt(i,k,j+1,:)) - &
                        chib(ii,k,j+1,:)
                   fls3(:) = (chinb(i-1,k,j,:)+dtb*chinbt(i-1,k,j,:)) - &
                        chib(ii-1,k,j,:)
                   fls4(:) = (chinb(i+1,k,j,:)+dtb*chinbt(i+1,k,j,:)) - &
                        chib(ii+1,k,j,:)
                   ften(ii,k,:) = ften(ii,k,:) + fcx*fls0(:) -                    &
                        & gcx*crdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-4.*fls0(:))
                end do
             end do
             ibeg = jsls
             iend = iym1 - jsls + 1
          end if
          !
          if ( jj.gt.ip ) then
             !-------west-boundary slice:
             do k = 1 , kz
                fcx = fcoef*xfune(jsls,k)
                gcx = gcoef*xfune(jsls,k)
                do i = ibeg , iend

                   fls0(:) = (chiwb(i,k,jwb,:)+dtb*chiwbt(i,k,jwb,:)) -chib(i,k,j,:)
                   fls1(:) = (chiwb(i-1,k,jwb,:)+dtb*chiwbt(i-1,k,jwb,:))             &
                        & -chib(i-1,k,j,:)
                   fls2(:) = (chiwb(i+1,k,jwb,:)+dtb*chiwbt(i+1,k,jwb,:))             &
                        & -chib(i+1,k,j,:)
                   fls3(:) = (chiwb(i,k,jwb-1,:)+dtb*chiwbt(i,k,jwb-1,:))             &
                        & -chib(i,k,j-1,:)
                   fls4(:) = (chiwb(i,k,jwb+1,:)+dtb*chiwbt(i,k,jwb+1,:))             &
                        & -chib(i,k,j+1,:)

                   ften(i,k,:) = ften(i,k,:) + fcx*fls0(:) -                      &
                        & gcx*crdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-4.*fls0(:))
                end do
             end do
          else if ( jj.le.ip ) then
             !-------east-boundary slice:
             do k = 1 , kz
                fcx = fcoef*xfune(jsls,k)
                gcx = gcoef*xfune(jsls,k)
                do i = ibeg , iend

                   fls0(:) = (chieb(i,k,jeb,:)+dtb*chiebt(i,k,jeb,:)) -chib(i,k,j,:)
                   fls1(:) = (chieb(i-1,k,jeb,:)+dtb*chiebt(i-1,k,jeb,:))             &
                        & -chib(i-1,k,j,:)
                   fls2(:) = (chieb(i+1,k,jeb,:)+dtb*chiebt(i+1,k,jeb,:))             &
                        & -chib(i+1,k,j,:)
                   fls3(:) = (chieb(i,k,jeb-1,:)+dtb*chiebt(i,k,jeb-1,:))             &
                        & -chib(i,k,j-1,:)
                   fls4(:) = (chieb(i,k,jeb+1,:)+dtb*chiebt(i,k,jeb+1,:))             &
                        & -chib(i,k,j+1,:)

                   ften(i,k,:) = ften(i,k,:) + fcx*fls0(:) -                      &
                        & gcx*crdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-4.*fls0(:))
                end do
             end do
          else
          end if
       else
       end if
       !
    else
    end if
#endif
!    call time_end(subroutine_name,idindx)
  end subroutine nudge_chi
  !

      function xfune(mm,kk)
        implicit none
        real(8) :: xfune
        integer , intent(in) :: mm , kk
        xfune = exp(-dble(mm-2)/canudg(kk))
        end function xfune



end module mod_che_bdyco
