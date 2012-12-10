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

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_mppparam  
  use mod_runparams
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
  public :: che_init_bdy , chib0 , chib1 , chibt , ichbdy2trac , oxcl

  type(rcm_time_and_date) , save :: chbdydate1 , chbdydate2

  real(rk8) , pointer , dimension(:,:,:,:) :: chib0 , chib1 , chibt , oxcl
  real(rk8) , pointer , dimension(:,:) :: cefc , cegc
  integer(ik4) , pointer , dimension(:) :: ichbdy2trac
  !
  ! Boundary conditions arrays
  !
  real(rk8) , pointer , dimension(:,:,:,:) :: chebdy

  integer(ik4) , parameter :: max_input_tracers = 50
  integer(ik4) , parameter :: noxcl = 5

  type cbound_area
    logical :: havebound
    logical , pointer , dimension(:,:) :: bsouth
    logical , pointer , dimension(:,:) :: bnorth
    logical , pointer , dimension(:,:) :: beast
    logical , pointer , dimension(:,:) :: bwest
    integer(ik4) :: ns , nn , ne , nw
    integer(ik4) :: nsp
    integer(ik4) , pointer , dimension(:,:) :: ibnd
  end type cbound_area
  public :: cbound_area
  type(cbound_area) , public :: cba

  contains
!
  subroutine allocate_mod_che_bdyco
    implicit none
    call getmem1d(ichbdy2trac,1,max_input_tracers,'che_bdyco:ichbdytrac')
    call getmem4d(chebdy,jce1,jce2,ice1,ice2,1,kz, &
                  1,max_input_tracers,'che_bdyco:chebdy')
    call getmem4d(chib0,jde1-ma%jbl1,jde2+ma%jbr1, &
                        ide1-ma%ibb1,ide2+ma%ibt1, &
                        1,kz,1,ntr,'mod_che_bdyco:chib0')
    call getmem4d(chib1,jde1-ma%jbl1,jde2+ma%jbr1, &
                        ide1-ma%ibb1,ide2+ma%ibt1, &
                        1,kz,1,ntr,'mod_che_bdyco:chib1')
    call getmem4d(chibt,jde1-ma%jbl1,jde2+ma%jbr1, &
                        ide1-ma%ibb1,ide2+ma%ibt1, &
                        1,kz,1,ntr,'mod_che_bdyco:chibt')
    call getmem2d(cefc,1,nspgx,1,kz,'che_bdyco:fcx')
    call getmem2d(cegc,1,nspgx,1,kz,'che_bdyco:fcx')
    if ( ioxclim == 1 ) then
      call getmem4d(oxcl,jde1-ma%jbl1,jde2+ma%jbr1, &
                         ide1-ma%ibb1,ide2+ma%ibt1, &
                         1,kz,1,noxcl,'che_bdyco:oxcl')
    end if 
  end subroutine allocate_mod_che_bdyco

  subroutine che_init_bdy
    implicit none
    integer(ik4) :: datefound , i , j , k , n, after
    character(len=32) :: appdat
    type (rcm_time_and_date) :: chbc_date
    integer(ik4) :: lyear , lmonth , lday , lhour
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'che_init_bdy'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    chbdydate1 = idate1
    chbdydate2 = idate1

    if ( ichsursrc == 1 ) then
      call split_idate(chbdydate1,lyear,lmonth,lday,lhour)
      call chem_emission(lyear,lmonth,lday,lhour)
    end if

    if ( ichebdy == 1 ) then 

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
        appdat = tochar(chbdydate1)
        call fatal(__FILE__,__LINE__,'CHBC for '//appdat//' not found')
      end if

      call read_chbc(chebdy)

      chib0 = d_zero
      after = 0
      do n = 1 , size(ichbdy2trac)
        if ( ichbdy2trac(n) > 0 ) then
          do k = 1 , kz
            do i = ice1 , ice2
              do j = jce1 , jce2
                chib0(j,i,k,ichbdy2trac(n)) = chebdy(j,i,k,n)
              end do
            end do
          end do
          after = after + 1
        end if
      end do

      ! handle oxidant climatology 
      if ( ioxclim == 1 ) then 
        do n = 1 , noxcl
          do k = 1 , kz
            do i = ice1 , ice2
              do j = jce1 , jce2
                oxcl(j,i,k,n) = chebdy(j,i,k,after+n)
              end do
            end do
          end do
        end do
      end if 
      
      if ( myid == italk ) then
        appdat = tochar(chbdydate1)
        if ( .not. ifrest ) then
          write(stdout,*) 'READY ICCH DATA for ', appdat
        else
          write(stdout,*) 'READY BCCH DATA for ', appdat
        end if
      end if

      chbdydate2 = chbdydate2 + intbdy
 
      if ( myid == italk ) then
        write (stdout,*) 'SEARCH CHBC data for ', toint10(chbdydate2)
      end if

      datefound = chbc_search(chbdydate2)
      if (datefound < 0) then
        call open_chbc(monfirst(chbdydate2))
        datefound = chbc_search(chbdydate2)
        if (datefound < 0) then
          appdat = tochar(chbdydate2)
          call fatal(__FILE__,__LINE__,'CHBC for '//appdat//' not found')
        end if
      end if

      call read_chbc(chebdy)

      chib1 = d_zero
      do n = 1 , size(ichbdy2trac)
        if ( ichbdy2trac(n) > 0 ) then
          do k = 1 , kz
            do i = ice1 , ice2
              do j = jce1 , jce2
                chib1(j,i,k,ichbdy2trac(n)) = chebdy(j,i,k,n)
              end do
            end do
          end do
        end if
      end do

      if ( myid == italk ) then
        write (stdout,*) 'READY  CHBC from     ' , &
            toint10(chbdydate1) , ' to ' , toint10(chbdydate2)
      end if

      chbdydate1 = chbdydate2

      !
      ! Couple with pstar
      !
      do n = 1 , ntr
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              chib0(j,i,k,n) = chib0(j,i,k,n)*psbb0(j,i)
            end do
          end do
        end do
      end do
      call exchange(chib0,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      !
      ! Repeat fot T2
      !
      do n = 1 , ntr
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              chib1(j,i,k,n) = chib1(j,i,k,n)*psbb1(j,i)
            end do
          end do
        end do
      end do
      call exchange(chib1,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
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
      call exchange(chibt,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)

      ! handle oxc lima

      if ( ioxclim == 1 ) then 
        call exchange(oxcl,1,jce1,jce2,ice1,ice2,1,kz,1,noxcl)
      end if
    else
      chbdydate2 = chbdydate2 + intbdy
      chbdydate1 = chbdydate2
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine che_init_bdy

  subroutine chem_bdyin
    implicit none
    integer(ik4) :: i , j , k , n , datefound, after
    character(len=32) :: appdat
    integer(ik4) :: lyear , lmonth , lday , lhour
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'chem_bdyin'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
  
    chib0(:,:,:,:) = chib1(:,:,:,:)

    if ( ichsursrc == 1 ) then
      call split_idate(chbdydate1,lyear,lmonth,lday,lhour)
      call chem_emission(lyear,lmonth,lday,lhour)
    end if

    chbdydate2 = chbdydate2 + intbdy

    if ( ichebdy == 1 ) then

      if ( myid == italk ) then
        write (stdout,*) 'SEARCH CHBC data for ', toint10(chbdydate2)
      end if

      datefound = chbc_search(chbdydate2)
      if (datefound < 0) then
        call open_chbc(monfirst(chbdydate2))
        datefound = chbc_search(chbdydate2)
        if (datefound < 0) then
          appdat = tochar(chbdydate2)
          call fatal(__FILE__,__LINE__,'CHBC for '//appdat//' not found')
        end if
      end if

      call read_chbc(chebdy)

      chib1 = d_zero
      after = 0  
      do n = 1 , size(ichbdy2trac)
        if ( ichbdy2trac(n) > 0 ) then
          do k = 1 , kz
            do i = ice1 , ice2
              do j = jce1 , jce2
                chib1(j,i,k,ichbdy2trac(n)) = chebdy(j,i,k,n)
              end do
            end do
          end do
          after = after + 1
        end if
      end do
      if ( ioxclim == 1 ) then
        do n = 1 , noxcl
          do k = 1 , kz
            do i = ice1 , ice2
              do j = jce1 , jce2
                oxcl(j,i,k,n) = chebdy(j,i,k,after+n)
              end do
            end do
          end do
        end do
      end if

      do n = 1 , ntr
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              chib1(j,i,k,n) = chib1(j,i,k,n)*psbb1(j,i)
            end do
          end do
        end do
      end do
      call exchange(chib1,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            chibt(j,i,k,:) = (chib1(j,i,k,:)-chib0(j,i,k,:))/dtbdys
          end do
        end do
      end do
      call exchange(chibt,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)

      ! handle oxidant climatology 
      if ( ioxclim == 1 ) then 
        call exchange(oxcl,1,jce1,jce2,ice1,ice2,1,kz,1,noxcl)
      end if

      if ( myid == italk ) then
        write (stdout,*) 'READY  CHBC from     ' , &
            toint10(chbdydate1) , ' to ' , toint10(chbdydate2)
      end if
    else
      chbdydate1 = chbdydate2
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine chem_bdyin

  subroutine chem_bdyval(xt)
    implicit none
    real(rk8) , intent(in) :: xt
!
    integer(ik4) :: itr , j , k , i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'chem_bdyval'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( ktau > 1 ) then
      !
      ! West boundary
      !
      if ( ma%has_bdyleft ) then
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
      if ( ma%has_bdyright ) then
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
      if ( ma%has_bdybottom ) then
        do itr = 1 , ntr
          do k = 1 , kz
            do j = jce1 , jce2
              chib(j,ice1,k,itr) = chia(j,ice1,k,itr)
            end do
          end do
        end do
      end if
      if ( ma%has_bdytop ) then
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

    if ( ma%has_bdyleft ) then
      do k = 1 , kz
        do i = ici1 , ici2
          chia(jce1,i,k,:) = chib0(jce1,i,k,:) + xt*chibt(jce1,i,k,:)
        end do
      end do
    end if
    if ( ma%has_bdyright ) then
      do k = 1 , kz
        do i = ici1 , ici2
          chia(jce2,i,k,:) = chib0(jce2,i,k,:) + xt*chibt(jce2,i,k,:)
        end do
      end do
    end if
    if ( ma%has_bdybottom ) then
      do k = 1 , kz
        do j = jce1 , jce2
          chia(j,ice1,k,:) = chib0(j,ice1,k,:) + xt*chibt(j,ice1,k,:)
        end do
      end do
    end if
    if ( ma%has_bdytop ) then
      do k = 1 , kz
        do j = jce1 , jce2
          chia(j,ice2,k,:) = chib0(j,ice2,k,:) + xt*chibt(j,ice2,k,:)
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
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
  subroutine nudge_chi(nk,xt,f,ften)
    implicit none
    integer(ik4) , intent(in) :: nk
    real(rk8) , intent(in) :: xt
    real(rk8) , pointer , intent(in) , dimension(:,:,:,:) :: f
    real(rk8) , pointer , intent(inout) , dimension(:,:,:,:) :: ften
!
    real(rk8) :: xf , xg
    real(rk8), dimension(ntr) :: fls0 , fls1 , fls2 , fls3 , fls4 

    integer(ik4) :: i , j , k , ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'nudge_chi'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
!
    if ( cba%ns /= 0 ) then
      do k = 1 , nk
        do i = ici1 , ici2
          do j = jci1 , jci2
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
                    xg*rdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-d_four*fls0(:))
          end do
        end do
      end do
    end if
    if ( cba%nn /= 0 ) then
      do k = 1 , nk
        do i = ici1 , ici2
          do j = jci1 , jci2
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
                     xg*rdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-d_four*fls0(:))
          end do
        end do
      end do
    end if
    if ( cba%nw /= 0 ) then
      do k = 1 , nk
        do i = ici1 , ici2
          do j = jci1 , jci2
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
                     xg*rdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-d_four*fls0(:))
          end do
        end do
      end do
    end if
    if ( cba%ne /= 0 ) then
      do k = 1 , nk
        do i = ici1 , ici2
          do j = jci1 , jci2
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
                     xg*rdxsq*(fls1(:)+fls2(:)+fls3(:)+fls4(:)-d_four*fls0(:))
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine nudge_chi
!
  subroutine setup_che_bdycon
    implicit none
    integer(ik4) :: n , k
    real(rk8) :: fnudge , gnudge
    !
    ! Specify the coefficients for nudging boundary conditions:
    !
    fnudge = 0.1D0/dt2
    gnudge = (dxsq/dt)/50.0D0
    do k = 1 , kz
      do n = 2 , nspgx-1
        cefc(n,k) = fnudge*xfune(n,k)
        cegc(n,k) = gnudge*xfune(n,k)
      end do
    end do
  end subroutine setup_che_bdycon
!
  function xfune(mm,kk)
    implicit none
    real(rk8) :: xfune
    integer(ik4) , intent(in) :: mm , kk
    xfune = dexp(-dble(mm-2)/anudg(kk))
  end function xfune
!
end module mod_che_bdyco
