!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_che_bdyco

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_date
  use mod_memutil
  use mod_nhinterp
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
  use mod_zita

  implicit none

  private

  public :: allocate_mod_che_bdyco, chem_bdyin, chem_bdyval
  public :: nudge_chi, setup_che_bdycon
  public :: che_init_bdy, chib0, chib1, chibt, ichbdy2trac, oxcl

  type(rcm_time_and_date), save :: chbdydate1, chbdydate2

  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chib0, chib1, chibt, oxcl
  real(rkx), pointer, contiguous, dimension(:,:) :: cefc, cegc
  integer(ik4), pointer, contiguous, dimension(:) :: ichbdy2trac
  !
  ! Boundary conditions arrays
  !
  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chebdy
  real(rkx), pointer, contiguous, dimension(:,:,:) :: fg

  integer(ik4), parameter :: max_input_tracers = 50
  integer(ik4), parameter :: noxcl = 5

  type cbound_area
    logical :: havebound
    logical, pointer, contiguous, dimension(:,:) :: bsouth
    logical, pointer, contiguous, dimension(:,:) :: bnorth
    logical, pointer, contiguous, dimension(:,:) :: beast
    logical, pointer, contiguous, dimension(:,:) :: bwest
    integer(ik4) :: ns, nn, ne, nw
    integer(ik4) :: nsp
    integer(ik4), pointer, contiguous, dimension(:,:) :: ibnd
  end type cbound_area
  public :: cbound_area
  type(cbound_area), public :: cba

  interface nudge_chi
    module procedure nudge_chiten
    module procedure monudgechi
  end interface nudge_chi

  interface chem_bdyval
    module procedure chem_bdyval_coupled
    module procedure chem_bdyval_uncoupled
  end interface chem_bdyval

  contains

  subroutine allocate_mod_che_bdyco
    implicit none
    call getmem1d(ichbdy2trac,1,max_input_tracers,'che_bdyco:ichbdytrac')
    call getmem4d(chib0,jde1ga,jde2ga,ide1ga,ide2ga, &
                        1,kz,1,ntr,'mod_che_bdyco:chib0')
    if ( ioxclim == 1 ) then
      call getmem4d(oxcl,jde1ga,jde2ga,ide1ga,ide2ga, &
                         1,kz,1,noxcl,'che_bdyco:oxcl')
    end if
    if ( ichebdy == 1 ) then
      call getmem4d(chebdy,jce1,jce2,ice1,ice2,1,kz, &
                    1,max_input_tracers,'che_bdyco:chebdy')
      call getmem4d(chib1,jde1ga,jde2ga,ide1ga,ide2ga, &
                          1,kz,1,ntr,'mod_che_bdyco:chib1')
      call getmem4d(chibt,jde1ga,jde2ga,ide1ga,ide2ga, &
                          1,kz,1,ntr,'mod_che_bdyco:chibt')
    end if
    call getmem2d(cefc,1,nspgx,1,kz,'che_bdyco:fcx')
    call getmem2d(cegc,1,nspgx,1,kz,'che_bdyco:fcx')
    call getmem3d(fg,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'che_bdyco:fg')
  end subroutine allocate_mod_che_bdyco

  subroutine che_init_bdy
    implicit none
    integer(ik4) :: datefound, i, j, k, n, after
    character(len=32) :: appdat
    type (rcm_time_and_date) :: chbc_date
    integer(ik4) :: lyear, lmonth, lday, lhour
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'che_init_bdy'
    integer(ik4), save :: idindx = 0
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
        call fatal(__FILE__,__LINE__, &
                   'CHBC for '//appdat//' not found')
      end if

      call read_chbc(chebdy)

      ! Interpolate to non-hydrostatic levels
      if ( idynamic == 2 ) then
        call nhinterp(ice1,ice2,jce1,jce2,kz,max_input_tracers, &
                      hsigma,sigma,chebdy,tvirt0,bndp0,cps0)
      end if
      if ( idynamic == 3 ) then
        call zita_interp(jce1,jce2,ice1,ice2,kz,max_input_tracers, &
                         chebdy,cza,tvirt0,hsigma,bndp0,0)
      end if

      chib0 = d_zero
      after = 0
      do n = 1, size(ichbdy2trac)
        if ( ichbdy2trac(n) > 0 ) then
          do k = 1, kz
            do i = ice1, ice2
              do j = jce1, jce2
                chib0(j,i,k,ichbdy2trac(n)) = chebdy(j,i,k,n)
              end do
            end do
          end do
          after = after + 1
        end if
      end do

      ! handle oxidant climatology
      if ( ioxclim == 1 ) then
        do n = 1, noxcl
          do k = 1, kz
            do i = ice1, ice2
              do j = jce1, jce2
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
          write(stdout,*) 'READY CHBC DATA for ', appdat
        end if
      end if

      chbdydate2 = chbdydate2 + intbdy

      if ( myid == italk ) then
        write (stdout,*) 'SEARCH CHBC data for ', tochar10(chbdydate2)
      end if

      datefound = chbc_search(chbdydate2)
      if (datefound < 0) then
        call open_chbc(monfirst(chbdydate2))
        datefound = chbc_search(chbdydate2)
        if (datefound < 0) then
          appdat = tochar(chbdydate2)
          call fatal(__FILE__,__LINE__, &
                     'CHBC for '//appdat//' not found')
        end if
      end if

      call read_chbc(chebdy)

      if ( idynamic == 2 ) then
        call nhinterp(ice1,ice2,jce1,jce2,kz,max_input_tracers, &
                      hsigma,sigma,chebdy,tvirt1,bndp1,cps0)
      end if
      if ( idynamic == 3 ) then
        call zita_interp(jce1,jce2,ice1,ice2,kz,max_input_tracers, &
                         chebdy,cza,tvirt1,hsigma,bndp1,0)
      end if

      chib1 = d_zero
      do n = 1, size(ichbdy2trac)
        if ( ichbdy2trac(n) > 0 ) then
          do k = 1, kz
            do i = ice1, ice2
              do j = jce1, jce2
                chib1(j,i,k,ichbdy2trac(n)) = chebdy(j,i,k,n)
              end do
            end do
          end do
        end if
      end do

      if ( myid == italk ) then
        write (stdout,*) 'READY  CHBC from     ', &
            tochar10(chbdydate1), ' to ', tochar10(chbdydate2)
      end if

      chbdydate1 = chbdydate2

      if ( idynamic /= 3 ) then
        !
        ! Couple with pstar
        !
        do n = 1, ntr
          do k = 1, kz
            do i = ice1, ice2
              do j = jce1, jce2
                chib0(j,i,k,n) = chib0(j,i,k,n)*psbb0(j,i)
              end do
            end do
          end do
        end do
        !
        ! Repeat fot T2
        !
        do n = 1, ntr
          do k = 1, kz
            do i = ice1, ice2
              do j = jce1, jce2
                chib1(j,i,k,n) = chib1(j,i,k,n)*psbb1(j,i)
              end do
            end do
          end do
        end do
      end if
      call exchange(chib0,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      call exchange(chib1,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      !
      ! Calculate time varying component
      !
      do n = 1, ntr
        do k = 1, kz
          do i = ice1ga, ice2ga
            do j = jce1ga, jce2ga
              chibt(j,i,k,n) = (chib1(j,i,k,n)-chib0(j,i,k,n))/dtbdys
            end do
          end do
        end do
      end do

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
    integer(ik4) :: i, j, k, n, datefound, after
    character(len=32) :: appdat
    integer(ik4) :: lyear, lmonth, lday, lhour
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'chem_bdyin'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( ichsursrc == 1 ) then
      call split_idate(chbdydate1,lyear,lmonth,lday,lhour)
      call chem_emission(lyear,lmonth,lday,lhour)
    end if

    chbdydate2 = chbdydate2 + intbdy

    if ( ichebdy == 1 ) then

      chib0(:,:,:,:) = chib1(:,:,:,:)

      if ( myid == italk ) then
        write (stdout,*) 'SEARCH CHBC data for ', tochar10(chbdydate2)
      end if

      datefound = chbc_search(chbdydate2)
      if (datefound < 0) then
        call open_chbc(monfirst(chbdydate2))
        datefound = chbc_search(chbdydate2)
        if (datefound < 0) then
          appdat = tochar(chbdydate2)
          call fatal(__FILE__,__LINE__, &
                     'CHBC for '//appdat//' not found')
        end if
      end if

      call read_chbc(chebdy)

      if ( idynamic == 2 ) then
        call nhinterp(ice1,ice2,jce1,jce2,kz,max_input_tracers, &
                      hsigma,sigma,chebdy,tvirt1,bndp1,cps0)
      end if
      if ( idynamic == 3 ) then
        call zita_interp(jce1,jce2,ice1,ice2,kz,max_input_tracers, &
                         chebdy,cza,tvirt1,hsigma,bndp1,0)
      end if

      chib1 = d_zero
      after = 0
      do n = 1, size(ichbdy2trac)
        if ( ichbdy2trac(n) > 0 ) then
          do k = 1, kz
            do i = ice1, ice2
              do j = jce1, jce2
                chib1(j,i,k,ichbdy2trac(n)) = chebdy(j,i,k,n)
              end do
            end do
          end do
          after = after + 1
        end if
      end do
      if ( ioxclim == 1 ) then
        do n = 1, noxcl
          do k = 1, kz
            do i = ice1, ice2
              do j = jce1, jce2
                oxcl(j,i,k,n) = chebdy(j,i,k,after+n)
              end do
            end do
          end do
        end do
      end if

      if ( idynamic /= 3 ) then
        do n = 1, ntr
          do k = 1, kz
            do i = ice1, ice2
              do j = jce1, jce2
                chib1(j,i,k,n) = chib1(j,i,k,n)*psbb1(j,i)
              end do
            end do
          end do
        end do
      end if
      call exchange(chib1,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      do n = 1, ntr
        do k = 1, kz
          do i = ice1ga, ice2ga
            do j = jce1ga, jce2ga
              chibt(j,i,k,n) = (chib1(j,i,k,n)-chib0(j,i,k,n))/dtbdys
            end do
          end do
        end do
      end do
      ! handle oxidant climatology
      if ( ioxclim == 1 ) then
        call exchange(oxcl,1,jce1,jce2,ice1,ice2,1,kz,1,noxcl)
      end if

      if ( myid == italk ) then
        write (stdout,*) 'READY  CHBC from     ', &
            tochar10(chbdydate1), ' to ', tochar10(chbdydate2)
      end if
      chbdydate1 = chbdydate2
    else
      chbdydate1 = chbdydate2
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine chem_bdyin

  subroutine chem_bdyval_uncoupled(u,v)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: u, v
    real(rkx) :: xt, windavg, trint
    integer(ik4) :: itr, i, j, k, n
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'chem_bdyval_uncoupled'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    xt = xbctime + dt
    if ( ichebdy == 0 ) then
      !
      ! flux dependent bdy
      ! west boundary:
      !
      if ( ma%has_bdyleft ) then
        do itr = 1, ntr
          do k = 1, kz
            do i = ice1, ice2
              trint = chemt(jci1,i,k,itr)
              windavg = u(jde1,i,k) - u(jdi1,i,k)
              if ( windavg < d_zero ) then
                chemt(jce1,i,k,itr) = trint
              else
                chemt(jce1,i,k,itr) = d_zero
              end if
            end do
          end do
        end do
      end if
      !
      ! east boundary:
      !
      if ( ma%has_bdyright ) then
        do itr = 1, ntr
          do k = 1, kz
            do i = ice1, ice2
              trint = chemt(jci2,i,k,itr)
              windavg = u(jde2,i,k) - u(jdi2,i,k)
              if ( windavg > d_zero ) then
                chemt(jce2,i,k,itr) = trint
              else
                chemt(jce2,i,k,itr) = d_zero
              end if
            end do
          end do
        end do
      end if
      !
      ! south boundary:
      !
      if ( ma%has_bdybottom ) then
        do itr = 1, ntr
          do k = 1, kz
            do j = jci1, jci2
              trint = chemt(j,ici1,k,itr)
              windavg = v(j,ide1,k) - v(j,idi1,k)
              if ( windavg < d_zero ) then
                chemt(j,ice1,k,itr) = trint
              else
                chemt(j,ice1,k,itr) = d_zero
              end if
            end do
          end do
        end do
      end if
      !
      ! north boundary:
      !
      if ( ma%has_bdytop ) then
        do itr = 1, ntr
          do k = 1, kz
            do j = jci1, jci2
              trint = chemt(j,ici2,k,itr)
              windavg = v(j,ide2,k) - v(j,idi2,k)
              if ( windavg > d_zero ) then
                chemt(j,ice2,k,itr) = trint
              else
                chemt(j,ice2,k,itr) = d_zero
              end if
            end do
          end do
        end do
      end if
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

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
      do n = 1, ntr
        do k = 1, kz
          do i = ici1, ici2
            chemt(jce1,i,k,n) = chib0(jce1,i,k,n) + xt*chibt(jce1,i,k,n)
          end do
        end do
      end do
    end if
    if ( ma%has_bdyright ) then
      do n = 1, ntr
        do k = 1, kz
          do i = ici1, ici2
            chemt(jce2,i,k,n) = chib0(jce2,i,k,n) + xt*chibt(jce2,i,k,n)
          end do
        end do
      end do
    end if
    if ( ma%has_bdybottom ) then
      do n = 1, ntr
        do k = 1, kz
          do j = jce1, jce2
            chemt(j,ice1,k,n) = chib0(j,ice1,k,n) + xt*chibt(j,ice1,k,n)
          end do
        end do
      end do
    end if
    if ( ma%has_bdytop ) then
      do n = 1, ntr
        do k = 1, kz
          do j = jce1, jce2
            chemt(j,ice2,k,n) = chib0(j,ice2,k,n) + xt*chibt(j,ice2,k,n)
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine chem_bdyval_uncoupled

  subroutine chem_bdyval_coupled(psa,wue,wui,eue,eui,nve,nvi,sve,svi)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: psa
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: wue, wui
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: eue, eui
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: nve, nvi
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: sve, svi
    real(rkx) :: xt, windavg, trint
    integer(ik4) :: itr, i, j, k, n
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'chem_bdyval_coupled'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    xt = xbctime + dt
    if ( rcmtimer%integrating( ) ) then
      !
      ! West boundary
      !
      if ( ma%has_bdyleft ) then
        do itr = 1, ntr
          do k = 1, kz
            do i = ice1, ice2
              chib(jce1,i,k,itr) = chia(jce1,i,k,itr)
            end do
          end do
        end do
      end if
      !
      ! East boundary
      !
      if ( ma%has_bdyright ) then
        do itr = 1, ntr
          do k = 1, kz
            do i = ice1, ice2
              chib(jce2,i,k,itr) = chia(jce2,i,k,itr)
            end do
          end do
        end do
      end if
      !
      ! North and South boundaries
      !
      if ( ma%has_bdybottom ) then
        do itr = 1, ntr
          do k = 1, kz
            do j = jci1, jci2
              chib(j,ice1,k,itr) = chia(j,ice1,k,itr)
            end do
          end do
        end do
      end if
      if ( ma%has_bdytop ) then
        do itr = 1, ntr
          do k = 1, kz
            do j = jci1, jci2
              chib(j,ice2,k,itr) = chia(j,ice2,k,itr)
            end do
          end do
        end do
      end if
    end if ! end if not start

    if ( ichebdy == 0 ) then
      !
      ! flux dependent bdy
      ! west boundary:
      !
      if ( ma%has_bdyleft ) then
        do itr = 1, ntr
          do k = 1, kz
            do i = ice1, ice2
              trint = chia(jci1,i,k,itr)/psa(jci1,i)
              windavg = wue(i,k) + wue(i+1,k) + wui(i,k) + wui(i+1,k)
              if ( windavg < d_zero ) then
                chia(jce1,i,k,itr) = trint*psa(jce1,i)
              else
                chia(jce1,i,k,itr) = d_zero
              end if
            end do
          end do
        end do
      end if
      !
      ! east boundary:
      !
      if ( ma%has_bdyright ) then
        do itr = 1, ntr
          do k = 1, kz
            do i = ice1, ice2
              trint = chia(jci2,i,k,itr)/psa(jci2,i)
              windavg = eue(i,k) + eue(i+1,k) + eui(i,k) + eui(i+1,k)
              if ( windavg > d_zero ) then
                chia(jce2,i,k,itr) = trint*psa(jce2,i)
              else
                chia(jce2,i,k,itr) = d_zero
              end if
            end do
          end do
        end do
      end if
      !
      ! south boundary:
      !
      if ( ma%has_bdybottom ) then
        do itr = 1, ntr
          do k = 1, kz
            do j = jci1, jci2
              trint = chia(j,ici1,k,itr)/psa(j,ici1)
              windavg = sve(j,k) + sve(j+1,k) + svi(j,k) + svi(j+1,k)
              if ( windavg < d_zero ) then
                chia(j,ice1,k,itr) = trint*psa(j,ice1)
              else
                chia(j,ice1,k,itr) = d_zero
              end if
            end do
          end do
        end do
      end if
      !
      ! north boundary:
      !
      if ( ma%has_bdytop ) then
        do itr = 1, ntr
          do k = 1, kz
            do j = jci1, jci2
              trint = chia(j,ici2,k,itr)/psa(j,ici2)
              windavg = nve(j,k) + nve(j+1,k) + nvi(j,k) + nvi(j+1,k)
              if ( windavg > d_zero ) then
                chia(j,ice2,k,itr) = trint*psa(j,ice2)
              else
                chia(j,ice2,k,itr) = d_zero
              end if
            end do
          end do
        end do
      end if
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

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
      do n = 1, ntr
        do k = 1, kz
          do i = ici1, ici2
            chia(jce1,i,k,n) = chib0(jce1,i,k,n) + xt*chibt(jce1,i,k,n)
          end do
        end do
      end do
    end if
    if ( ma%has_bdyright ) then
      do n = 1, ntr
        do k = 1, kz
          do i = ici1, ici2
            chia(jce2,i,k,n) = chib0(jce2,i,k,n) + xt*chibt(jce2,i,k,n)
          end do
        end do
      end do
    end if
    if ( ma%has_bdybottom ) then
      do n = 1, ntr
        do k = 1, kz
          do j = jce1, jce2
            chia(j,ice1,k,n) = chib0(j,ice1,k,n) + xt*chibt(j,ice1,k,n)
          end do
        end do
      end do
    end if
    if ( ma%has_bdytop ) then
      do n = 1, ntr
        do k = 1, kz
          do j = jce1, jce2
            chia(j,ice2,k,n) = chib0(j,ice2,k,n) + xt*chibt(j,ice2,k,n)
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine chem_bdyval_coupled

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
  subroutine nudge_chiten(f,ften)
    implicit none
    real(rkx), pointer, contiguous, intent(in), dimension(:,:,:,:) :: f
    real(rkx), pointer, contiguous, intent(inout), dimension(:,:,:,:) :: ften

    real(rkx) :: xt, xf, xg
    real(rkx) :: fls0, fls1, fls2, fls3, fls4

    integer(ik4) :: i, j, k, n, ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'nudge_chiten'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( ichebdy == 0 ) then
      if ( cba%ns /= 0 ) then
        do n = 1, ntr
          do k = 1, kz
            do i = ici1, ici2
              do j = jci1, jci2
                if ( .not. cba%bsouth(j,i) ) cycle
                ib = cba%ibnd(j,i)
                xf = cefc(ib,k)
                xg = cegc(ib,k)
                fls0 = f(j-1,ici1,k,n) - f(j,i,k,n)
                fls1 = f(j-1,ici1,k,n) - f(j-1,i,k,n)
                fls2 = f(j+1,ici1,k,n) - f(j+1,i,k,n)
                fls3 = f(j,ice1,k,n)   - f(j,i-1,k,n)
                fls4 = f(j,ici1+1,k,n) - f(j,i+1,k,n)
                ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                        xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end do
    end if
    if ( cba%nn /= 0 ) then
      do n = 1, ntr
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              if ( .not. cba%bnorth(j,i) ) cycle
              ib = cba%ibnd(j,i)
              xf = cefc(ib,k)
              xg = cegc(ib,k)
              fls0 = f(j,ici2,k,n) - f(j,i,k,n)
              fls1 = f(j-1,ici2,k,n) - f(j-1,i,k,n)
              fls2 = f(j+1,ici2,k,n) - f(j+1,i,k,n)
              fls3 = f(j,ici2-1,k,n) - f(j,i-1,k,n)
              fls4 = f(j,ice2,k,n) - f(j,i+1,k,n)
              ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                       xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end do
    end if
    if ( cba%nw /= 0 ) then
      do n = 1, ntr
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              if ( .not. cba%bwest(j,i) ) cycle
              ib = cba%ibnd(j,i)
              xf = cefc(ib,k)
              xg = cegc(ib,k)
              fls0 = f(jci2,i,k,n) - f(j,i,k,n)
              fls1 = f(jci2-1,i-1,k,n) - f(j,i-1,k,n)
              fls2 = f(jci2,i+1,k,n) - f(j,i+1,k,n)
              fls3 = f(jci2-1,i,k,n) - f(j-1,i,k,n)
              fls4 = f(jce2,i,k,n) - f(j+1,i,k,n)
              ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                       xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end do
    end if
    if ( cba%ne /= 0 ) then
      do n = 1, ntr
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              if ( .not. cba%beast(j,i) ) cycle
              ib = cba%ibnd(j,i)
              xf = cefc(ib,k)
              xg = cegc(ib,k)
              fls0 = f(jci1,i,k,n) - f(j,i,k,n)
              fls1 = f(jci1,i-1,k,n) - f(j,i-1,k,n)
              fls2 = f(jci1,i+1,k,n) - f(j,i+1,k,n)
              fls3 = f(jce1,i,k,n) - f(j-1,i,k,n)
              fls4 = f(jci1+1,i,k,n) - f(j+1,i,k,n)
              ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 -  &
                       xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end do
    end if
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    xt = xbctime + dt
    if ( cba%ns /= 0 ) then
      do n = 1, ntr
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              if ( .not. cba%bsouth(j,i) ) cycle
              ib = cba%ibnd(j,i)
              xf = cefc(ib,k)
              xg = cegc(ib,k)
              fls0 = (chib0(j,i,k,n)  +xt*chibt(j,i,k,n))   - f(j,i,k,n)
              fls1 = (chib0(j-1,i,k,n)+xt*chibt(j-1,i,k,n)) - f(j-1,i,k,n)
              fls2 = (chib0(j+1,i,k,n)+xt*chibt(j+1,i,k,n)) - f(j+1,i,k,n)
              fls3 = (chib0(j,i-1,k,n)+xt*chibt(j,i-1,k,n)) - f(j,i-1,k,n)
              fls4 = (chib0(j,i+1,k,n)+xt*chibt(j,i+1,k,n)) - f(j,i+1,k,n)
              ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                      xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end do
    end if
    if ( cba%nn /= 0 ) then
      do n = 1, ntr
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              if ( .not. cba%bnorth(j,i) ) cycle
              ib = cba%ibnd(j,i)
              xf = cefc(ib,k)
              xg = cegc(ib,k)
              fls0 = (chib0(j,i,k,n)  +xt*chibt(j,i,k,n))   - f(j,i,k,n)
              fls1 = (chib0(j-1,i,k,n)+xt*chibt(j-1,i,k,n)) - f(j-1,i,k,n)
              fls2 = (chib0(j+1,i,k,n)+xt*chibt(j+1,i,k,n)) - f(j+1,i,k,n)
              fls3 = (chib0(j,i-1,k,n)+xt*chibt(j,i-1,k,n)) - f(j,i-1,k,n)
              fls4 = (chib0(j,i+1,k,n)+xt*chibt(j,i+1,k,n)) - f(j,i+1,k,n)
              ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                       xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end do
    end if
    if ( cba%nw /= 0 ) then
      do n = 1, ntr
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              if ( .not. cba%bwest(j,i) ) cycle
              ib = cba%ibnd(j,i)
              xf = cefc(ib,k)
              xg = cegc(ib,k)
              fls0 = (chib0(j,i,k,n)  +xt*chibt(j,i,k,n))   - f(j,i,k,n)
              fls1 = (chib0(j,i-1,k,n)+xt*chibt(j,i-1,k,n)) - f(j,i-1,k,n)
              fls2 = (chib0(j,i+1,k,n)+xt*chibt(j,i+1,k,n)) - f(j,i+1,k,n)
              fls3 = (chib0(j-1,i,k,n)+xt*chibt(j-1,i,k,n)) - f(j-1,i,k,n)
              fls4 = (chib0(j+1,i,k,n)+xt*chibt(j+1,i,k,n)) - f(j+1,i,k,n)
              ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                       xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end do
    end if
    if ( cba%ne /= 0 ) then
      do n = 1, ntr
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              if ( .not. cba%beast(j,i) ) cycle
              ib = cba%ibnd(j,i)
              xf = cefc(ib,k)
              xg = cegc(ib,k)
              fls0 = (chib0(j,i,k,n)  +xt*chibt(j,i,k,n))   - f(j,i,k,n)
              fls1 = (chib0(j,i-1,k,n)+xt*chibt(j,i-1,k,n)) - f(j,i-1,k,n)
              fls2 = (chib0(j,i+1,k,n)+xt*chibt(j,i+1,k,n)) - f(j,i+1,k,n)
              fls3 = (chib0(j-1,i,k,n)+xt*chibt(j-1,i,k,n)) - f(j-1,i,k,n)
              fls4 = (chib0(j+1,i,k,n)+xt*chibt(j+1,i,k,n)) - f(j+1,i,k,n)
              ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 -  &
                       xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine nudge_chiten

  subroutine monudgechi(f)
    implicit none
    real(rkx), pointer, contiguous, intent(inout), dimension(:,:,:,:) :: f

    real(rkx) :: xt, xf, xg
    real(rkx) :: fls0, fls1, fls2, fls3, fls4

    integer(ik4) :: i, j, k, n, ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'monudgechi'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( ichebdy == 0 ) then
      do n = 1, ntr
        do concurrent ( j = jce1ga:jce2ga, i = ice1ga:ice2ga, k = 1:kz )
          fg(j,i,k) = f(j,i,k,n)
        end do
        if ( cba%ns /= 0 ) then
          do k = 1, kz
            do i = ici1, ici2
              do j = jci1, jci2
                if ( .not. cba%bsouth(j,i) ) cycle
                ib = cba%ibnd(j,i)
                xf = cefc(ib,k)
                xg = cegc(ib,k)
                fls0 = fg(j-1,ici1,k) - fg(j,i,k)
                fls1 = fg(j-1,ici1,k) - fg(j-1,i,k)
                fls2 = fg(j+1,ici1,k) - fg(j+1,i,k)
                fls3 = fg(j,ice1,k)   - fg(j,i-1,k)
                fls4 = fg(j,ici1+1,k) - fg(j,i+1,k)
                f(j,i,k,n) = f(j,i,k,n) + xf*fls0 - &
                        xg*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          end do
        end if
        if ( cba%nn /= 0 ) then
          do k = 1, kz
            do i = ici1, ici2
              do j = jci1, jci2
                if ( .not. cba%bnorth(j,i) ) cycle
                ib = cba%ibnd(j,i)
                xf = cefc(ib,k)
                xg = cegc(ib,k)
                fls0 = fg(j,ici2,k) - fg(j,i,k)
                fls1 = fg(j-1,ici2,k) - fg(j-1,i,k)
                fls2 = fg(j+1,ici2,k) - fg(j+1,i,k)
                fls3 = fg(j,ici2-1,k) - fg(j,i-1,k)
                fls4 = fg(j,ice2,k) - fg(j,i+1,k)
                f(j,i,k,n) = f(j,i,k,n) + xf*fls0 - &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          end do
        end if
        if ( cba%nw /= 0 ) then
          do k = 1, kz
            do i = ici1, ici2
              do j = jci1, jci2
                if ( .not. cba%bwest(j,i) ) cycle
                ib = cba%ibnd(j,i)
                xf = cefc(ib,k)
                xg = cegc(ib,k)
                fls0 = fg(jci2,i,k) - fg(j,i,k)
                fls1 = fg(jci2-1,i-1,k) - fg(j,i-1,k)
                fls2 = fg(jci2,i+1,k) - fg(j,i+1,k)
                fls3 = fg(jci2-1,i,k) - fg(j-1,i,k)
                fls4 = fg(jce2,i,k) - fg(j+1,i,k)
                f(j,i,k,n) = f(j,i,k,n) + xf*fls0 - &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          end do
        end if
        if ( cba%ne /= 0 ) then
          do k = 1, kz
            do i = ici1, ici2
              do j = jci1, jci2
                if ( .not. cba%beast(j,i) ) cycle
                ib = cba%ibnd(j,i)
                xf = cefc(ib,k)
                xg = cegc(ib,k)
                fls0 = fg(jci1,i,k) - fg(j,i,k)
                fls1 = fg(jci1,i-1,k) - fg(j,i-1,k)
                fls2 = fg(jci1,i+1,k) - fg(j,i+1,k)
                fls3 = fg(jce1,i,k) - fg(j-1,i,k)
                fls4 = fg(jci1+1,i,k) - fg(j+1,i,k)
                f(j,i,k,n) = f(j,i,k,n) + xf*fls0 -  &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          end do
        end if
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    xt = xbctime + dt
    do n = 1, ntr
      do concurrent ( j = jce1ga:jce2ga, i = ice1ga:ice2ga, k = 1:kz )
        fg(j,i,k) = (chib0(j,i,k,n) + xt*chibt(j,i,k,n)) - f(j,i,k,n)
      end do
      if ( cba%ns /= 0 ) then
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              if ( .not. cba%bsouth(j,i) ) cycle
              ib = cba%ibnd(j,i)
              xf = cefc(ib,k)
              xg = cegc(ib,k)
              fls0 = fg(j,i,k)
              fls1 = fg(j-1,i,k)
              fls2 = fg(j+1,i,k)
              fls3 = fg(j,i-1,k)
              fls4 = fg(j,i+1,k)
              f(j,i,k,n) = f(j,i,k,n) + xf*fls0 - &
                      xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( cba%nn /= 0 ) then
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              if ( .not. cba%bnorth(j,i) ) cycle
              ib = cba%ibnd(j,i)
              xf = cefc(ib,k)
              xg = cegc(ib,k)
              fls0 = fg(j,i,k)
              fls1 = fg(j-1,i,k)
              fls2 = fg(j+1,i,k)
              fls3 = fg(j,i-1,k)
              fls4 = fg(j,i+1,k)
              f(j,i,k,n) = f(j,i,k,n) + xf*fls0 - &
                       xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( cba%nw /= 0 ) then
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              if ( .not. cba%bwest(j,i) ) cycle
              ib = cba%ibnd(j,i)
              xf = cefc(ib,k)
              xg = cegc(ib,k)
              fls0 = fg(j,i,k)
              fls1 = fg(j-1,i,k)
              fls2 = fg(j+1,i,k)
              fls3 = fg(j,i-1,k)
              fls4 = fg(j,i+1,k)
              f(j,i,k,n) = f(j,i,k,n) + xf*fls0 - &
                       xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( cba%ne /= 0 ) then
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              if ( .not. cba%beast(j,i) ) cycle
              ib = cba%ibnd(j,i)
              xf = cefc(ib,k)
              xg = cegc(ib,k)
              fls0 = fg(j,i,k)
              fls1 = fg(j-1,i,k)
              fls2 = fg(j+1,i,k)
              fls3 = fg(j,i-1,k)
              fls4 = fg(j,i+1,k)
              f(j,i,k,n) = f(j,i,k,n) + xf*fls0 -  &
                       xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine monudgechi

  subroutine setup_che_bdycon
    implicit none
    integer(ik4) :: n, k
    real(rkx) :: fnudge, gnudge
    real(rkx), dimension(kz) :: anudgh
    !
    ! Specify the coefficients for nudging boundary conditions:
    !
    if ( bdy_nm > d_zero ) then
      fnudge = bdy_nm
    else
      fnudge = 0.1_rkx/dt2
    end if
    if ( bdy_dm > d_zero ) then
      gnudge = bdy_dm
    else
      gnudge = 0.02_rkx/dt2
    end if
    call exponential_nudging(anudgh)
    do k = 1, kz
      do n = 2, nspgx-1
        cefc(n,k) = fnudge*xfune(n,anudgh(k))
        cegc(n,k) = gnudge*xfune(n,anudgh(k))
      end do
    end do
    contains
      pure real(rkx) function xfune(mm,an)
        implicit none
        integer(ik4), intent(in) :: mm
        real(rkx), intent(in) :: an
        xfune = exp(-real(mm-2,rkx)/an)
      end function xfune
  end subroutine setup_che_bdycon

end module mod_che_bdyco

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
