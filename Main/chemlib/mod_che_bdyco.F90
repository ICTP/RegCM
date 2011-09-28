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

  use mod_dynparam
  use mod_service
  use mod_message
  use mod_memutil
  use mod_che_common
  use mod_che_mppio
  use mod_che_species

  private

  public :: allocate_mod_che_bdyco , chem_bdyin , chem_bdyval

#ifndef BAND
  public :: chieb, chiebt , chiwb , chiwbt
#endif
  public :: chisb , chisbt , chinb , chinbt
  public :: chib0 , chib1

  real(dp) , allocatable , dimension(:,:,:,:) :: chib0 , chib1
!
#ifndef BAND
  real(dp) , allocatable , dimension(:,:,:,:) :: chieb , chiebt , &
                                                 chiwb , chiwbt
#endif
  real(dp) , allocatable , dimension(:,:,:,:) :: chinb , chinbt , &
                                                 chisb , chisbt

  contains

!
  subroutine allocate_mod_chem_bdycon
    implicit none

    character (len=64) :: subroutine_name='allocate_mod_che_bdyco'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)

    call getmem4d(chib0,1,iy,1,kz,1,jxp,1,totsp,'che_bdyco:chib0')
    call getmem4d(chib1,1,iy,1,kz,1,jxp,1,totsp,'che_bdyco:chib1')
#ifndef BAND
    call getmem4d(chieb,1,iy,1,kz,0,jxp+1,1,totsp,'che_bdyco:chieb')
    call getmem4d(chiebt,1,iy,1,kz,0,jxp+1,1,totsp,'che_bdyco:chiebt')
    call getmem4d(chiwb,1,iy,1,kz,0,jxp+1,1,totsp,'che_bdyco:chiwb')
    call getmem4d(chiwbt,1,iy,1,kz,0,jxp+1,1,totsp,'che_bdyco:chiwbt')
#endif
    call getmem4d(chinb,1,nspgx,1,kz,0,jxp+1,1,totsp,'che_bdyco:chinb')
    call getmem4d(chinbt,1,nspgx,1,kz,0,jxp+1,1,totsp,'che_bdyco:chinbt')
    call getmem4d(chisb,1,nspgx,1,kz,0,jxp+1,1,totsp,'che_bdyco:chisb')
    call getmem4d(chisbt,1,nspgx,1,kz,0,jxp+1,1,totsp,'che_bdyco:chinst')

    call time_end(subroutine_name,idindx)

  end subroutine allocate_mod_chem_bdycon
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine reads in the boundary conditions.               c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine chem_bdyin

    use mod_che_indices
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif

    real(dp) :: dtbdys
    integer :: i , j , k , nn , nnb , mmrec, itr
    integer ::  nxeb , nxwb, ierr
#ifndef BAND
    integer :: nkk
#endif
    real(8) , dimension(iy,jxp) :: psdot , tdum
    integer :: n
    character (len=64) :: subroutine_name='chem_bdyin'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
    if ( dabs(xtime).gt.0.0001 ) return
!
    dtbdys = ibdyfrq*60.*60.

    if ( myid .eq. 0 .and. igaschem==0)then

    mmrec = oxcl_search(ndate1)

    print*, 'mmrec',mmrec
      if (mmrec < 0) then
     call open_oxcl(imonfirst(ndate1))
      print*,'after open_oxcl'
      end if


       call read_oxcl(ndate1,ohc0_io,ho2c0_io,o3c0_io,no3c0_io, h2o2c0_io)

    do j = 1 , jx
       do k = 1 , kz
          do i = 1 , iy
             savch_0(i,k       ,j)  = ohc0_io(i,k,j)
             savch_0(i,kz    +k,j)  = ho2c0_io(i,k,j)
             savch_0(i,kz*2  +k,j)  = o3c0_io(i,k,j)
             savch_0(i,kz*3  +k,j)  = no3c0_io(i,k,j)
             savch_0(i,kz*4  +k,j)  = h2o2c0_io(i,k,j)
          end do
       end do
    end do

    print*, 'OXCLIM READING', maxval(h2o2c0_io), minval(h2o2c0_io)

    end if

!       Start transmission of data to other processors
! nb here use already existing savch* arrays (dimensionned to 25), even if 5 variables are relevant( maybe more in the future)
    call mpi_scatter(savch_0,iy*kz*25*jxp,mpi_real8,         &
                     savch0, iy*kz*25*jxp,mpi_real8,         &
                     0,mpi_comm_world,ierr)


    do j = 1 , jendl
       do k = 1 , kz
          do i = 1 , iy
             ohc0(i,k,j)       = savch0(i,      k,j)
             ho2c0(i,k,j)      = savch0(i,kz  +k,j)
             o3c0(i,k,j)       = savch0(i,kz*2  +k,j)
             no3c0(i,k,j)      = savch0(i,kz*3  +k,j)
             h2o2c0(i,k,j)     = savch0(i,kz*4 +k,j)
          end do
       end do
    end do

    print*, 'after', myid, maxval(h2o2c0), minval(h2o2c0)


savch0  = 0.
savch_0 =0.


!  call allocate_species_bc
!  call allocate_mod_chem_mppio(.false.)

    if(myid .eq. 0 .and. igaschem==1) then

!      call addhours(ndate1, ibdyfrq)
    mmrec = chbc_search(ndate1)
    write(*,*)'BDYIN -----',mmrec,ndate1
    if (mmrec < 0) then
            call open_chbc(imonfirst(ndate1))
    end if


       call read_chbc(ndate1,                                 &
            o3b1_io,nob1_io,no2b1_io,hno3b1_io,n2o5b1_io      &
           ,h2o2b1_io,ch4b1_io,cob1_io,hchob1_io,ch3ohb1_io   &
           ,c2h5ohb1_io,c2h4b1_io,c2h6b1_io,ch3chob1_io       &
           ,c3h6b1_io,c3h8b1_io,ch3coch3b1_io,bigeneb1_io     &
           ,bigalkb1_io,isopb1_io,tolueb1_io,panb1_io         &
           ,so2b1_io,dmsb1_io,so4b1_io)


    do j = 1 , jx
       do k = 1 , kz
          do i = 1 , iy
             savch_0(i,k       ,j)  = o3b1_io(i,k,j)
             savch_0(i,kz    +k,j)  = nob1_io(i,k,j)
             savch_0(i,kz*2  +k,j)  = no2b1_io(i,k,j)
             savch_0(i,kz*3  +k,j)  = hno3b1_io(i,k,j)
             savch_0(i,kz*4  +k,j)  = n2o5b1_io(i,k,j)
             savch_0(i,kz*5  +k,j)  = h2o2b1_io(i,k,j)
             savch_0(i,kz*6  +k,j)  = ch4b1_io(i,k,j)
             savch_0(i,kz*7  +k,j)  = cob1_io(i,k,j)
             savch_0(i,kz*8  +k,j)  = hchob1_io(i,k,j)
             savch_0(i,kz*9  +k,j)  = ch3ohb1_io(i,k,j)
             savch_0(i,kz*10  +k,j) = c2h5ohb1_io(i,k,j)
             savch_0(i,kz*11  +k,j) = c2h4b1_io(i,k,j)
             savch_0(i,kz*12  +k,j) = c2h6b1_io(i,k,j)
             savch_0(i,kz*13  +k,j) = ch3chob1_io(i,k,j)
             savch_0(i,kz*14  +k,j) = c3h6b1_io(i,k,j)
             savch_0(i,kz*15  +k,j) = c3h8b1_io(i,k,j)
             savch_0(i,kz*16  +k,j) = ch3coch3b1_io(i,k,j)
             savch_0(i,kz*17  +k,j) = bigeneb1_io(i,k,j)
             savch_0(i,kz*18  +k,j) = bigalkb1_io(i,k,j)
             savch_0(i,kz*19  +k,j) = isopb1_io(i,k,j)
             savch_0(i,kz*20  +k,j) = tolueb1_io(i,k,j)
             savch_0(i,kz*21  +k,j) = panb1_io(i,k,j)
             savch_0(i,kz*22  +k,j) = so2b1_io(i,k,j)
             savch_0(i,kz*23  +k,j) = dmsb1_io(i,k,j)
             savch_0(i,kz*24  +k,j) = so4b1_io(i,k,j)
          end do
       end do
    end do

    end if
    call mpi_scatter(savch_0,iy*kz*25*jxp,mpi_real8,      &
                     savch0, iy*kz*25*jxp,mpi_real8,      &
                     0,mpi_comm_world,ierr)

       print*,' CIAO ,',myid, size(no2b1,3), size(savch0,3),jendl,jbegin,jendx,jxp

    do j = 1 , jendl
    do k = 1 , kz
          do i = 1 , iy
             o3b1(i,k,j)        = savch0(i,       k,j)
             nob1(i,k,j)        = savch0(i,kz    +k,j)
             no2b1(i,k,j)       = savch0(i,kz*2  +k,j)
             hno3b1(i,k,j)      = savch0(i,kz*3  +k,j)
             n2o5b1(i,k,j)      = savch0(i,kz*4  +k,j)
             h2o2b1(i,k,j)      = savch0(i,kz*5  +k,j)
             ch4b1(i,k,j)       = savch0(i,kz*6  +k,j)
             cob1(i,k,j)        = savch0(i,kz*7  +k,j)
             hchob1(i,k,j)      = savch0(i,kz*8  +k,j)
             ch3ohb1(i,k,j)     = savch0(i,kz*9  +k,j)
             c2h5ohb1(i,k,j)    = savch0(i,kz*10  +k,j)
             c2h4b1(i,k,j)      = savch0(i,kz*11  +k,j)
             c2h6b1(i,k,j)      = savch0(i,kz*12  +k,j)
             ch3chob1(i,k,j)    = savch0(i,kz*13  +k,j)
             c3h6b1(i,k,j)      = savch0(i,kz*14  +k,j)
             c3h8b1(i,k,j)      = savch0(i,kz*15  +k,j)
             ch3coch3b1(i,k,j)  = savch0(i,kz*16  +k,j)
             bigeneb1(i,k,j)    = savch0(i,kz*17  +k,j)
             bigalkb1(i,k,j)    = savch0(i,kz*18  +k,j)
             isopb1(i,k,j)      = savch0(i,kz*19  +k,j)
             tolueb1(i,k,j)     = savch0(i,kz*20  +k,j)
             panb1(i,k,j)       = savch0(i,kz*21  +k,j)
             so2b1(i,k,j)       = savch0(i,kz*22  +k,j)
             dmsb1(i,k,j)       = savch0(i,kz*23  +k,j)
             so4b1(i,k,j)       = savch0(i,kz*24  +k,j)
          end do
       end do
    end do
    do j = jbegin , jendx
       do i = 2 , iym1
           psdot(i,j) = 0.25*(ps1(i,j)   + ps1(i-1,j) +  &
                        ps1(i,j-1) + ps1(i-1,j-1))
       end do
    end do
    do i = 2 , iym1
       if ( myid.eq.0 ) psdot(i,1) = 0.5*(ps1(i,1)+ps1(i-1,1))
       if ( myid.eq.nproc-1 )                             &
           psdot(i,jendl) = 0.5*(ps1(i,jendx)+ps1(i-1,jendx))
    end do
    do j = jbegin , jendx
       psdot(1,j) = 0.5*(ps1(1,j)+ps1(1,j-1))
       psdot(iy,j) = 0.5*(ps1(iym1,j)+ps1(iym1,j-1))
    end do

    if ( myid.eq.0 ) then
         psdot(1,1) = ps1(1,1)
         psdot(iy,1) = ps1(iym1,1)
    end if

    if ( myid.eq.nproc-1 ) then
       psdot(1,jendl) = ps1(1,jendx)
       psdot(iy,jendl) = ps1(iym1,jendx)
    end if
    write(*,*)' BDYIN---ZZZ--',maxval(o3b1(:,1,:)),maxval(psdot(:,:)),myid
!=======================================================================
!       Couple pressure u,v,t,q

! fab modif
! intialise the chib1 tracer : correspondance betwwen bdy tables and chib

    do k = 1 , kz
       do j = 1 , jendl
          do i = 1 , iy
              if(io3 > 0) chib1(i,k,j,io3)     =  o3b1(i,k,j)*psdot(i,j)
              if(ino > 0) chib1(i,k,j,ino)     =  nob1(i,k,j)*psdot(i,j)
              if(ino2 > 0) chib1(i,k,j,ino2)    =  no2b1(i,k,j)*psdot(i,j)
              if(ihno3 > 0) chib1(i,k,j,ihno3)   =  hno3b1(i,k,j)*psdot(i,j)
              if(in2o5 > 0) chib1(i,k,j,in2o5)   =  n2o5b1(i,k,j)*psdot(i,j)
              if(ih2o2 > 0) chib1(i,k,j,ih2o2)   =  h2o2b1(i,k,j)*psdot(i,j)
              if(ich4 > 0) chib1(i,k,j,ich4)    =  ch4b1(i,k,j)*psdot(i,j)
              if(ico > 0) chib1(i,k,j,ico)     =  cob1(i,k,j)*psdot(i,j)
              if(ihcho > 0) chib1(i,k,j,ihcho)   =  hchob1(i,k,j)*psdot(i,j)
              if(imoh > 0) chib1(i,k,j,imoh)    =  ch3ohb1(i,k,j)*psdot(i,j)
              if(ieoh > 0) chib1(i,k,j,ieoh)    =  c2h5ohb1(i,k,j)*psdot(i,j)
              if(iethe > 0) chib1(i,k,j,iethe)   =  c2h4b1(i,k,j)*psdot(i,j)
              if(ic2h6 > 0) chib1(i,k,j,ic2h6)   =  c2h6b1(i,k,j)*psdot(i,j)
              if(iald2 > 0) chib1(i,k,j,iald2)   =  ch3chob1(i,k,j)*psdot(i,j)
              if(iolt > 0) chib1(i,k,j,iolt)    =  c3h6b1(i,k,j)*psdot(i,j)
              if(ic3h8 > 0) chib1(i,k,j,ic3h8)   =  c3h8b1(i,k,j)*psdot(i,j)
              if(iacet > 0) chib1(i,k,j,iacet)   =  ch3coch3b1(i,k,j)*psdot(i,j)
              if(ioli > 0) chib1(i,k,j,ioli)    =  bigeneb1(i,k,j)*psdot(i,j)
              if(iisop > 0) chib1(i,k,j,iisop)   =  isopb1(i,k,j)*psdot(i,j)
              if(itolue > 0) chib1(i,k,j,itolue)  =  tolueb1(i,k,j)*psdot(i,j)
              if(ipan > 0) chib1(i,k,j,ipan)    =  panb1(i,k,j)*psdot(i,j)
              if(iso2 > 0) chib1(i,k,j,iso2)    =  so2b1(i,k,j)*psdot(i,j)
              if(idms > 0) chib1(i,k,j,idms)    =  dmsb1(i,k,j)*psdot(i,j)
              if(iso4 > 0) chib1(i,k,j,iso4)    =  so4b1(i,k,j)*psdot(i,j)

          end do
       end do
    end do


! fab modif ici
!-----compute boundary conditions for p*chi
! first define nxeb and nxwb
!


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
    do k = 1 , kz
      do j = 1 , jendl
        do i = 1 , iy
          chib0(i,k,j,:) = chib1(i,k,j,:)
        end do
      end do
    end do
  end subroutine chem_bdyin

  subroutine chem_bdyval(xt,iexec)
    use mod_che_indices
!
    implicit none
!
    integer :: iexec
    real(8) :: xt
    intent (in) iexec , xt
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
    call time_begin(subroutine_name,idindx)
!
!*********************************************************************
!*****fill up the boundary value for xxb variables from xxa variables:
!     if this subroutine is called for the first time, this part
!     shall be skipped.

!gaffe je considere defaut mpp

!
    if ( iexec.ne.1 ) then
!
      if ( ichem.eq.1 ) then
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
      end if
!chem2_
!
!
    end if      !end if(iexec.ne.1) test
!**********************************************************************
!*****compute the boundary values for xxa variables:
!
!-----compute the time interval for boundary tendency:
!
    dtb = xt*60.
    if ( dabs(xt).lt.0.00001 .and. ldatez.gt.idate0 )                 &
       & dtb = ibdyfrq*60.*60.
!
!
!-----set boundary values for p*t:
!-----set boundary values for p*qv:
!

    if (iexec.eq.1 .and. ifrest)  return


    if ( iboudy.eq.0 ) then
!.....fixed boundary conditions:
    end if
!
!.....time-dependent boundary conditions:
!
!
! FAB modif chemistry time dependant boundary conditions
! default case here

    do k = 1 , kz
!fabtest remplace chia par chib
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

! for the
    if ( iboudy.eq.3 .or. iboudy.eq.4 ) then
!
!-----determine boundary values depends on inflow/outflow:
!overwrite the dfault case : keep it maybe for certainspecific  tracers in the future

    if ( ichem.eq.1 ) then
!chem2

!----add tracer bc's
!
      do itr = 1 , ntr
        do k = 1 , kz
#ifndef BAND
!
!.....west  boundary:
!
          if ( myid.eq.0 ) then
            do i = 1 , iym1
              chix1 = chia(i,k,1,itr)/sps1%ps(i,1)
              chix2 = chia(i,k,2,itr)/sps1%ps(i,2)
              uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
              if ( uavg.ge.0. ) then
                chix = chix1
              else
                chix = chix2
              end if
              chia(i,k,1,itr) = chix*sps1%ps(i,1)
            end do
          end if
!
!.....east  boundary:
!
          if ( myid.eq.nproc-1 ) then
            do i = 1 , iym1
              chix1 = chia(i,k,jendx,itr)/sps1%ps(i,jendx)
              chix2 = chia(i,k,jendm,itr)/sps1%ps(i,jendm)
              uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
              if ( uavg.lt.0. ) then
                chix = chix1

              else
                chix = chix2
              end if
              chia(i,k,jendx,itr) = chix*sps1%ps(i,jendx)
            end do
          end if
#endif
!
!.....south boundary:
!
          do j = jbegin , jendm
            chix1 = chia(1,k,j,itr)/sps1%ps(1,j)
            chix2 = chia(2,k,j,itr)/sps1%ps(2,j)
            vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
            if ( vavg.ge.0. ) then
              chix = chix1
! fab test
            else
              chix = chix2
            end if
            chia(1,k,j,itr) = chix*sps1%ps(1,j)
          end do
!
!.....north boundary:
!
          do j = jbegin , jendm
            chix1 = chia(iym1,k,j,itr)/sps1%ps(iym1,j)
            chix2 = chia(iym2,k,j,itr)/sps1%ps(iym2,j)
            vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
            if ( vavg.lt.0. ) then
              chix = chix1
            else
              chix = chix2
            end if
            chia(iym1,k,j,itr) = chix*sps1%ps(iym1,j)
          end do
        end do
      end do
!chem2_
    end if
!!$
    end if      !end if(iboudy.eq.3.or.4) test
!
!
!!$!
    call time_end(subroutine_name,idindx)
  end subroutine chem_bdyval
!
end module mod_che_bdyco
