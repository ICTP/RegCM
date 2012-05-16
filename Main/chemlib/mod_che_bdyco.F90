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
  use mod_mppparam  
  use mod_service
  use mod_mpmessage
  use mod_che_common
  use mod_che_mppio
  use mod_che_ncio 
  use mod_che_species
  use mod_che_indices
  use mod_che_emission
  use mod_mppparam  

  private

  public :: allocate_mod_che_bdyco , chem_bdyin , chem_bdyval
  public :: nudge_chi , setup_che_bdycon
  public :: che_init_bdy , chib0 , chib1 , chibt , ichbdy2trac , chebdy

  type(rcm_time_and_date) , save :: chbdydate1 , chbdydate2

  real(dp) , pointer , dimension(:,:,:,:) :: chib0 , chib1 , chibt , chebdy
  real(dp) , pointer , dimension(:,:) :: cefc , cegc
  integer , pointer , dimension(:) :: ichbdy2trac
   
  integer :: cnbdm

  contains

!
  subroutine allocate_mod_che_bdyco
    implicit none
    cnbdm = max(nspgx,nspgd)
    call getmem4d(chib0,jde1-ma%jbl1,jde2+ma%jbr1, &
                        ide1-ma%ibb1,ide2+ma%ibt1, &
                        1,kz,1,ntr,'mod_che_bdyco:chib0')
    call getmem4d(chib1,jde1-ma%jbl1,jde2+ma%jbr1, &
                        ide1-ma%ibb1,ide2+ma%ibt1, &
                        1,kz,1,ntr,'mod_che_bdyco:chib1')
    call getmem4d(chibt,jde1-ma%jbl1,jde2+ma%jbr1, &
                        ide1-ma%ibb1,ide2+ma%ibt1, &
                        1,kz,1,ntr,'mod_che_bdyco:chibt')
    call getmem4d(chebdy,jde1-ma%jbl1,jde2+ma%jbr1, &
                        ide1-ma%ibb1,ide2+ma%ibt1, &
                        1,kz,1,ntr,'mod_che_bdyco:chebdy')
    call getmem1d(ichbdy2trac,1,25,'mod_che_bdyco:ichbdytrac')
    call getmem2d(cefc,1,cnbdm,1,kz,'bdycon:fcx')
    call getmem2d(cegc,1,cnbdm,1,kz,'bdycon:fcx')
  end subroutine allocate_mod_che_bdyco

  subroutine che_init_bdy(idate1,intbdy,dtbdys,ifrest)
    implicit none
    logical :: ifrest
    integer :: datefound , i , j , k , ierr,n
    real (dp) :: dtbdys
    character(len=32) :: appdat
    type (rcm_time_and_date) :: idate1, chbc_date
    type(rcm_time_interval)::  intbdy
    
    character (len=64) :: subroutine_name='che_init_bdy'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)

    chbdydate1 = idate1
    chbdydate2 = idate1

    if ( myid == 0 ) then
      if ( chbdydate1 == globidate1 ) then
        chbc_date = chbdydate1
      else
        chbc_date = monfirst(chbdydate1)
      end if

      call open_chbc(chbc_date)

      datefound = chbc_search(chbdydate1)
      if (datefound < 0) then
        !
        ! Cannot run without initial conditions
        !
        appdat = tochar(chbdydate2)
        call fatal(__FILE__,__LINE__,'CHBC for '//appdat//' not found')
      end if

      call read_chbc(chebdy_in)
      chebdy_io0 =d_zero
      do n=1,size(ichbdy2trac)
             if(ichbdy2trac(n) > 0) chebdy_io0(:,:,:,ichbdy2trac(n)) = chebdy_in(:,:,:,n)
      end do

      appdat = tochar(chbdydate1)
      if ( .not. ifrest ) then
        write (6,*) 'READY ICCH  DATA for ', appdat
      else
        write (6,*) 'READY BCCH DATA for ', appdat
      end if

      chbdydate2 = chbdydate2 + intbdy
 
      write (6,*) 'SEARCH CHBC data for ', toint10(chbdydate2)
      datefound = chbc_search(chbdydate2)
      if (datefound < 0) then
        call open_chbc(monfirst(chbdydate2))
        datefound = chbc_search(chbdydate2)
        if (datefound < 0) then
          appdat = tochar(chbdydate2)
          call fatal(__FILE__,__LINE__,'CHBC for '//appdat//' not found')
        end if
      end if
      call read_chbc(chebdy_in )
      chebdy_io1 = d_zero
      do n=1,size(ichbdy2trac)
             if(ichbdy2trac(n) > 0) chebdy_io1(:,:,:,ichbdy2trac(n)) = chebdy_in(:,:,:,n)
      end do


      write (6,*) 'READY  CHBC from     ' , &
            toint10(chbdydate1) , ' to ' , toint10(chbdydate2)

    end if

    call date_bcast(chbdydate2,0,mycomm,ierr)
    chbdydate1 = chbdydate2
    !
    ! Send each processor its computing slice
    !
    call deco1_scatter(chebdy_io0,chebdy, &
                       jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
  

    do n = 1 , ntr
  
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              chib0(j,i,k,n) = chebdy(j,i,k,n)*cpsb(j,i)
            end do
          end do
        end do
  
    end do
    call deco1_exchange_left(chib0,1,ice1,ice2,1,kz,1,ntr)
    call deco1_exchange_right(chib0,1,ice1,ice2,1,kz,1,ntr)
  

    !
    ! Repeat fot T2
    !
    call deco1_scatter(chebdy_io1,chebdy, &
                       jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
    do n = 1 , ntr

        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              chib1(j,i,k,n) = chebdy(j,i,k,n)*cpsb(j,i)
            end do
          end do
        end do

    end do
    call deco1_exchange_left(chib1,1,ice1,ice2,1,kz,1,ntr)
    call deco1_exchange_right(chib1,1,ice1,ice2,1,kz,1,ntr)
    !
    ! Calculate time varying component
    !
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          chibt(j,i,k,:) = (chib1(j,i,k,:)-chib0(j,i,k,:))/dtbdys
        end do
      end do
    end do
    call deco1_exchange_left(chibt,1,ice1,ice2,1,kz,1,ntr)
    call deco1_exchange_right(chibt,1,ice1,ice2,1,kz,1,ntr)

    call time_end(subroutine_name,idindx)
  end subroutine che_init_bdy

  subroutine chem_bdyin(dtbdys,intbdy)
    implicit none
    type(rcm_time_interval) :: intbdy
    real(dp) , intent(in) :: dtbdys
    integer :: i , j , k , n , mmrec
    character(len=32) :: appdat
    integer :: lyear , lmonth , lday , lhour

    character (len=64) :: subroutine_name='chem_bdyin'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
  
    chbdydate2 = chbdydate2 + intbdy
    call split_idate(chbdydate2,lyear,lmonth,lday,lhour)

    chib0(:,:,:,:) = chib1(:,:,:,:)

    if ( myid == 0 ) then
      write (6,'(a,i10)') 'SEARCH BC data for ', toint10(chbdydate2)
      mmrec = chbc_search(chbdydate2)
      if (mmrec < 0) then
        call open_chbc(monfirst(chbdydate2))
        mmrec = chbc_search(chbdydate2)
        if (mmrec < 0) then
          appdat = tochar(chbdydate2)
          call fatal(__FILE__,__LINE__,'chBC for '//appdat//' not found')
        end if
      end if
      call read_chbc(chebdy_in)
      chebdy_io1 = d_zero
      do n = 1 , size(ichbdy2trac)
        if ( ichbdy2trac(n) > 0 ) then
          chebdy_io1(:,:,:,ichbdy2trac(n)) = chebdy_in(:,:,:,n)
        end if
      end do
    end if
 
    call deco1_scatter(chebdy_io1,chebdy, &
                       jcross1,jcross2,icross1,icross2,1,kz,1,ntr)

    do n = 1 , ntr

        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              chib1(j,i,k,n) = chebdy(j,i,k,n)*cpsb(j,i)
            end do
          end do
        end do

    end do
    call deco1_exchange_left(chib1,1,ice1,ice2,1,kz,1,ntr)
    call deco1_exchange_right(chib1,1,ice1,ice2,1,kz,1,ntr)

    do k=1,kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          chibt(j,i,k,:) = (chib1(j,i,k,:)-chib0(j,i,k,:))/dtbdys
        end do
      end do
    end do
    call deco1_exchange_left(chibt,1,ice1,ice2,1,kz,1,ntr)
    call deco1_exchange_right(chibt,1,ice1,ice2,1,kz,1,ntr)

    call chem_emission(lyear,lmonth,lday,lhour)

    call time_end(subroutine_name,idindx)
  end subroutine chem_bdyin

  subroutine chem_bdyval(xt,ktau)
    implicit none
    integer(8) , intent(in) :: ktau
    real(dp) , intent(in) :: xt
!
    integer :: itr , j , k , i
    character (len=64) :: subroutine_name='chem_bdyval'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
    if ( ktau > 1 ) then
      !
      ! West boundary
      !
      if ( ma%hasleft ) then
        do itr = 1 , ntr
          do k = 1 , kz
            do i = ici1 , ici2
              chib(jce1,i,k,itr) = chia(jce1,i,k,itr)
            end do
          end do
        end do
      end if
      !
      ! East boundary
      !
      if ( ma%hasright ) then
        do itr = 1 , ntr
          do k = 1 , kz
            do i = ici1 , ici2
              chib(jce2,i,k,itr) = chia(jce2,i,k,itr)
            end do
          end do
        end do
      end if
      !
      ! North and South boundaries
      !
      if ( ma%hasbottom ) then
        do itr = 1 , ntr
          do k = 1 , kz
            do j = jce1 , jce2
              chib(j,ice1,k,itr) = chia(j,ice1,k,itr)
            end do
          end do
        end do
      end if
      if ( ma%hastop ) then
        do itr = 1 , ntr
          do k = 1 , kz
            do j = jce1 , jce2
              chib(j,ice2,k,itr) = chia(j,ice2,k,itr)
            end do
          end do
        end do
      end if
    end if  !end if (ktau > 1) test

!.....time-dependent boundary conditions:
! for chemistry relaxation towrds
! time dependant boundary conditions is considered
!    if ( iboudy.eq.0 ) then
!.....fixed boundary conditions:
!    end if
!
!.....time-dependent boundary conditions:
! for chemistry relaxation towrds
! time dependant boundary conditions is considered

    if ( ma%hasleft ) then
      do k = 1 , kz
        do i = ici1 , ici2
          chia(jce1,i,k,:) = chib0(jce1,i,k,:) + xt*chibt(jce1,i,k,:)
        end do
      end do
    end if
    if ( ma%hasright ) then
      do k = 1 , kz
        do i = ici1 , ici2
          chia(jce2,i,k,:) = chib0(jce2,i,k,:) + xt*chibt(jce2,i,k,:)
        end do
      end do
    end if
    if ( ma%hasbottom ) then
      do k = 1 , kz
        do j = jce1 , jce2
          chia(j,ice1,k,:) = chib0(j,ice1,k,:) + xt*chibt(j,ice1,k,:)
        end do
      end do
    end if
    if ( ma%hastop ) then
      do k = 1 , kz
        do j = jce1 , jce2
          chia(j,ice2,k,:) = chib0(j,ice2,k,:) + xt*chibt(j,ice2,k,:)
        end do
      end do
    end if
    call time_end(subroutine_name,idindx)

  end subroutine chem_bdyval

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
  subroutine nudge_chi(nk,cba,xt,f,ften)
    implicit none
    integer , intent(in) :: nk
    real(dp) , intent(in) :: xt
    real(dp) , pointer , intent(in) , dimension(:,:,:,:) :: f
    type(cbound_area), intent(in) :: cba   
    real(dp) , pointer , intent(inout) , dimension(:,:,:,:) :: ften
!
    real(dp) :: xf , xg
    real(dp), dimension(ntr) :: fls0 , fls1 , fls2 , fls3 , fls4 

    integer :: i , j , k , ib , i1 , i2 , j1 , j2
    character (len=64) :: subroutine_name='nudge3d'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
    if ( cba%dotflag ) then
      i1 = idi1
      i2 = idi2
      j1 = jdi1
      j2 = jdi2
    else
      i1 = ici1
      i2 = ici2
      j1 = jci1
      j2 = jci2
    end if
!
    if ( cba%ns /= 0 ) then
      do k = 1 , nk
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. cba%bsouth(j,i) ) cycle
            ib = cba%ibnd(j,i)
            xf = cefc(ib,k)
            xg = cegc(ib,k)
            fls0(:) = (chib0(j,i,k,:)  +xt*chibt(j,i,k,:))   - f(j,i,k,:)
            fls1(:) = (chib0(j-1,i,k,:)+xt*chibt(j-1,i,k,:)) - f(j-1,i,k,:)
            fls2(:) = (chib0(j+1,i,k,:)+xt*chibt(j+1,i,k,:)) - f(j+1,i,k,:)
            fls3(:) = (chib0(j,i-1,k,:)+xt*chibt(j,i-1,k,:)) - f(j,i-1,k,:)
            fls4(:) = (chib0(j,i+1,k,:)+xt*chibt(j,i+1,k,:)) - f(j,i+1,k,:)
            ften(j,i,k,:) = ften(j,i,k,:) + xf*fls0(:) - &
                    xg*crdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-d_four*fls0(:))
          end do
        end do
      end do
    end if
    if ( cba%nn /= 0 ) then
      do k = 1 , nk
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. cba%bnorth(j,i) ) cycle
            ib = cba%ibnd(j,i)
            xf = cefc(ib,k)
            xg = cegc(ib,k)
            fls0(:) = (chib0(j,i,k,:)  +xt*chibt(j,i,k,:))   - f(j,i,k,:)
            fls1(:) = (chib0(j-1,i,k,:)+xt*chibt(j-1,i,k,:)) - f(j-1,i,k,:)
            fls2(:) = (chib0(j+1,i,k,:)+xt*chibt(j+1,i,k,:)) - f(j+1,i,k,:)
            fls3(:) = (chib0(j,i-1,k,:)+xt*chibt(j,i-1,k,:)) - f(j,i-1,k,:)
            fls4(:) = (chib0(j,i+1,k,:)+xt*chibt(j,i+1,k,:)) - f(j,i+1,k,:)
            ften(j,i,k,:) = ften(j,i,k,:) + xf*fls0(:) - &
                     xg*crdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-d_four*fls0(:))
          end do
        end do
      end do
    end if
    if ( cba%nw /= 0 ) then
      do k = 1 , nk
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. cba%bwest(j,i) ) cycle
            ib = cba%ibnd(j,i)
            xf = cefc(ib,k)
            xg = cegc(ib,k)
            fls0(:) = (chib0(j,i,k,:)  +xt*chibt(j,i,k,:))   - f(j,i,k,:)
            fls1(:) = (chib0(j,i-1,k,:)+xt*chibt(j,i-1,k,:)) - f(j,i-1,k,:)
            fls2(:) = (chib0(j,i+1,k,:)+xt*chibt(j,i+1,k,:)) - f(j,i+1,k,:)
            fls3(:) = (chib0(j-1,i,k,:)+xt*chibt(j-1,i,k,:)) - f(j-1,i,k,:)
            fls4(:) = (chib0(j+1,i,k,:)+xt*chibt(j+1,i,k,:)) - f(j+1,i,k,:)
            ften(j,i,k,:) = ften(j,i,k,:) + xf*fls0(:) - &
                     xg*crdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-d_four*fls0(:))
          end do
        end do
      end do
    end if
    if ( cba%ne /= 0 ) then
      do k = 1 , nk
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. cba%beast(j,i) ) cycle
            ib = cba%ibnd(j,i)
            xf = cefc(ib,k)
            xg = cegc(ib,k)
            fls0(:) = (chib0(j,i,k,:)  +xt*chibt(j,i,k,:))   - f(j,i,k,:)
            fls1(:) = (chib0(j,i-1,k,:)+xt*chibt(j,i-1,k,:)) - f(j,i-1,k,:)
            fls2(:) = (chib0(j,i+1,k,:)+xt*chibt(j,i+1,k,:)) - f(j,i+1,k,:)
            fls3(:) = (chib0(j-1,i,k,:)+xt*chibt(j-1,i,k,:)) - f(j-1,i,k,:)
            fls4(:) = (chib0(j+1,i,k,:)+xt*chibt(j+1,i,k,:)) - f(j+1,i,k,:)
            ften(j,i,k,:) = ften(j,i,k,:) + xf*fls0(:) -  &
                     xg*crdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-d_four*fls0(:))
          end do
        end do
      end do
    end if

    call time_end(subroutine_name,idindx)
  end subroutine nudge_chi
!
  subroutine setup_che_bdycon
    implicit none
    integer :: n , k
    real(dp) :: fnudge , gnudge
    !
    ! Specify the coefficients for nudging boundary conditions:
    !
    fnudge = 0.1D0/ ( 2 * dtche)
    gnudge = (1.0D0/crdxsq/dtche)/50.0D0
    do k = 1 , kz
      do n = 2 , cnbdm-1
        cefc(n,k) = fnudge*xfune(n,k)
        cegc(n,k) = gnudge*xfune(n,k)
      end do
    end do
  end subroutine setup_che_bdycon
!
  function xfune(mm,kk)
    implicit none
    real(dp) :: xfune
    integer , intent(in) :: mm , kk
    xfune = dexp(-dble(mm-2)/canudg(kk))
  end function xfune
!
end module mod_che_bdyco
