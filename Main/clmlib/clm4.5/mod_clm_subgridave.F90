module mod_clm_subgridave
  !
  ! Utilities to perfrom subgrid averaging
  !
  use mod_realkinds
  use mod_intkinds
  use mod_stdio
  use mod_mpmessage
  use mod_clm_type
  use mod_clm_varcon , only : spval , isturb
  use mod_clm_varcon , only : icol_roof , icol_sunwall , icol_shadewall
  use mod_clm_varcon , only : icol_road_perv , icol_road_imperv
  use mod_clm_varpar , only : max_pft_per_col
  use mod_clm_varpar , only : max_pft_per_lu
  use mod_clm_varpar,  only : max_pft_per_gcell
  use mod_clm_varcon , only : max_lunit
  use mod_clm_varcon , only : istsoil, istice, istdlak, istslak, istwet, &
                              isturb, istcrop, max_lunit, spval

  implicit none

  private

  save

  public :: p2c   ! Perfrom an average from pfts to columns
  public :: p2l   ! Perfrom an average from pfts to landunits
  public :: p2g   ! Perfrom an average from pfts to gridcells
  public :: c2l   ! Perfrom an average from columns to landunits
  public :: c2g   ! Perfrom an average from columns to gridcells
  public :: l2g   ! Perfrom an average from landunits to gridcells

  interface p2c
    module procedure p2c_1d
    module procedure p2c_2d
    module procedure p2c_1d_filter
    module procedure p2c_2d_filter
  end interface
  interface p2l
    module procedure p2l_1d
    module procedure p2l_2d
  end interface
  interface p2g
    module procedure p2g_1d
    module procedure p2g_2d
  end interface
  interface c2l
    module procedure c2l_1d
    module procedure c2l_2d
  end interface
  interface c2g
    module procedure c2g_1d
    module procedure c2g_2d
  end interface
  interface l2g
    module procedure l2g_1d
    module procedure l2g_2d
  end interface

  ! - I believe that scale_p2c, scale_c2l and scale_l2g should be included in
  !   the sumwt accumulations (e.g., sumwt = sumwt + wtgcell * scale_p2c *
  !   scale_c2l * scale_l2g), but that requires some more thought to (1) make
  !   sure that is correct, and (2) make sure it doesn't break the urban
  !   scaling. (See also my notes in create_scale_l2g_lookup.)
  ! - Once that is done, you could use a scale of 0, avoiding the need for
  !   the use of spval and the special checks that requires.
  ! - Currently, there is a lot of repeated code to calculate scale_c2l.
  !   This should be cleaned up.
  ! - At a minimum, should collect the repeated code into a subroutine to
  !   eliminate this repitition
  ! - The best thing might be to use a lookup array, as is done for scale_l2g

  contains

  !
  ! Perfrom subgrid-average from pfts to columns.
  ! Averaging is only done for points that are not equal to "spval".
  !
  subroutine p2c_1d(lbp,ubp,lbc,ubc,parr,carr,p2c_scale_type)
    implicit none
    integer(ik4) , intent(in) :: lbp , ubp  ! beginning and ending pft
    integer(ik4) , intent(in) :: lbc , ubc  ! beginning and ending column
    real(rk8) , intent(in) , dimension(lbp:ubp) :: parr  ! pft array
    real(rk8) , intent(out) , dimension(lbc:ubc) :: carr ! column array
    character(len=*) , intent(in) :: p2c_scale_type ! scale type
    integer(ik4) :: p , c , iind   ! indices
    ! scale factor for column->landunit mapping
    real(rk8) , dimension(lbp:ubp) :: scale_p2c
    logical  :: found  ! temporary for error check
    real(rk8) , dimension(lbc:ubc) :: sumwt  ! sum of weights
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer , dimension(:) :: pactive
    ! weight of pft relative to column
    real(rk8) , pointer , dimension(:) :: wtcol
    ! column index of corresponding pft
    integer(ik4) , pointer , dimension(:) :: pcolumn
    ! number of pfts in column
    integer(ik4) , pointer , dimension(:) :: npfts
    ! initial pft index in column
    integer(ik4) , pointer , dimension(:) :: pfti

    pactive  => clm3%g%l%c%p%active
    wtcol    => clm3%g%l%c%p%wtcol
    pcolumn  => clm3%g%l%c%p%column
    npfts    => clm3%g%l%c%npfts
    pfti     => clm3%g%l%c%pfti

    if (p2c_scale_type == 'unity') then
      do p = lbp,ubp
        scale_p2c(p) = 1.0D0
      end do
    else
      write(stderr,*)'p2c_1d error: scale type ',p2c_scale_type,' not supported'
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if

    carr(lbc:ubc) = spval
    sumwt(lbc:ubc) = 0.D0
    do p = lbp , ubp
      if ( pactive(p) .and. wtcol(p) /= 0.D0 ) then
        if ( parr(p) /= spval ) then
          c = pcolumn(p)
          if (sumwt(c) == 0.D0) carr(c) = 0.D0
          carr(c) = carr(c) + parr(p) * scale_p2c(p) * wtcol(p)
          sumwt(c) = sumwt(c) + wtcol(p)
        end if
      end if
    end do
    found = .false.
    do c = lbc , ubc
      if ( sumwt(c) > 1.0D0 + 1.D-4 ) then
        found = .true.
        iind = c
      else if ( sumwt(c) /= 0.D0 ) then
        carr(c) = carr(c)/sumwt(c)
      end if
    end do
    if (found) then
      write(stderr,*)'p2c error: sumwt is greater than 1.0 at c= ',iind
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if
  end subroutine p2c_1d
  !
  ! Perfrom subgrid-average from landunits to gridcells.
  ! Averaging is only done for points that are not equal to "spval".
  !
  subroutine p2c_2d(lbp,ubp,lbc,ubc,num2d,parr,carr,p2c_scale_type)
    implicit none
    integer(ik4) , intent(in) :: lbp , ubp  ! beginning and ending pft
    integer(ik4) , intent(in) :: lbc , ubc  ! beginning and ending column
    integer(ik4) , intent(in) :: num2d      ! size of second dimension
    real(rk8) , intent(in) , dimension(lbp:ubp,num2d) :: parr  ! pft array
    real(rk8) , intent(out) , dimension(lbc:ubc,num2d) :: carr ! column array
    character(len=*) , intent(in) :: p2c_scale_type ! scale type
    integer(ik4) :: j , p , c , iind         ! indices
    ! scale factor for column->landunit mapping
    real(rk8) , dimension(lbp:ubp) :: scale_p2c
    logical :: found                  ! temporary for error check
    ! sum of weights
    real(rk8) , dimension(lbc:ubc) :: sumwt(lbc:ubc)
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer , dimension(:) :: pactive
    ! weight of pft relative to column
    real(rk8) , pointer , dimension(:) :: wtcol
    ! column index of corresponding pft
    integer(ik4) , pointer , dimension(:) :: pcolumn
    ! number of pfts in column
    integer(ik4) , pointer , dimension(:) :: npfts
    ! initial pft index in column
    integer(ik4) , pointer , dimension(:) :: pfti

    pactive  => clm3%g%l%c%p%active
    wtcol    => clm3%g%l%c%p%wtcol
    pcolumn  => clm3%g%l%c%p%column
    npfts    => clm3%g%l%c%npfts
    pfti     => clm3%g%l%c%pfti

    if (p2c_scale_type == 'unity') then
      do p = lbp,ubp
        scale_p2c(p) = 1.0D0
      end do
    else
      write(stderr,*)'p2c_2d error: scale type ',p2c_scale_type,' not supported'
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if

    carr(:,:) = spval
    do j = 1 , num2d
      sumwt(:) = 0.D0
      do p = lbp , ubp
        if (pactive(p) .and. wtcol(p) /= 0.D0) then
          if (parr(p,j) /= spval) then
            c = pcolumn(p)
            if (sumwt(c) == 0.D0) carr(c,j) = 0.D0
            carr(c,j) = carr(c,j) + parr(p,j) * scale_p2c(p) * wtcol(p)
            sumwt(c) = sumwt(c) + wtcol(p)
          end if
        end if
      end do
      found = .false.
      do c = lbc , ubc
        if (sumwt(c) > 1.0D0 + 1.D-4) then
          found = .true.
          iind = c
        else if (sumwt(c) /= 0.D0) then
          carr(c,j) = carr(c,j)/sumwt(c)
        end if
      end do
      if (found) then
        write(stderr,*) &
             'p2c_2d error: sumwt is greater than 1.0 at c= ',iind,' lev= ',j
        call fatal(__FILE__,__LINE__,'clm is now stopping')
      end if
    end do
  end subroutine p2c_2d
  !
  ! perform pft to column averaging for single level pft arrays
  !
  subroutine p2c_1d_filter(numfc,filterc,pftarr,colarr)
    implicit none
    integer(ik4) , intent(in) :: numfc
    integer(ik4) , intent(in) , dimension(numfc) :: filterc
    real(rk8) , pointer , dimension(:) :: pftarr
    real(rk8) , pointer , dimension(:) :: colarr
    integer(ik4) :: fc , c , p  ! indices
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer , dimension(:) :: pactive
    integer(ik4) , pointer , dimension(:) :: npfts
    integer(ik4) , pointer , dimension(:) :: pfti
    integer(ik4) , pointer , dimension(:) :: pftf
    real(rk8) , pointer , dimension(:) :: wtcol

    pactive => clm3%g%l%c%p%active
    npfts   => clm3%g%l%c%npfts
    pfti    => clm3%g%l%c%pfti
    pftf    => clm3%g%l%c%pftf
    wtcol   => clm3%g%l%c%p%wtcol

    do fc = 1,numfc
      c = filterc(fc)
      colarr(c) = 0.D0
      do p = pfti(c), pftf(c)
        if (pactive(p)) colarr(c) = colarr(c) + pftarr(p) * wtcol(p)
      end do
    end do
  end subroutine p2c_1d_filter
  !
  ! perform pft to column averaging for multi level pft arrays
  !
  subroutine p2c_2d_filter(lev,numfc,filterc,pftarr,colarr)
    implicit none
    integer(ik4) , intent(in) :: lev
    integer(ik4) , intent(in) :: numfc
    integer(ik4) , intent(in) , dimension(numfc) :: filterc
    real(rk8) , pointer , dimension(:,:) :: pftarr
    real(rk8) , pointer , dimension(:,:) :: colarr
    integer(ik4) :: fc , c , p , j    ! indices
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer , dimension(:) :: pactive
    integer(ik4) , pointer , dimension(:) :: npfts
    integer(ik4) , pointer , dimension(:) :: pfti
    integer(ik4) , pointer , dimension(:) :: pftf
    real(rk8) , pointer , dimension(:) :: wtcol

    pactive => clm3%g%l%c%p%active
    npfts   => clm3%g%l%c%npfts
    pfti    => clm3%g%l%c%pfti
    pftf    => clm3%g%l%c%pftf
    wtcol   => clm3%g%l%c%p%wtcol

    do j = 1,lev
      do fc = 1,numfc
        c = filterc(fc)
        colarr(c,j) = 0.D0
        do p = pfti(c), pftf(c)
          if (pactive(p)) colarr(c,j) = colarr(c,j) + pftarr(p,j) * wtcol(p)
        end do
      end do
    end do
  end subroutine p2c_2d_filter
  !
  ! Perfrom subgrid-average from pfts to landunits
  ! Averaging is only done for points that are not equal to "spval".
  !
  subroutine p2l_1d(lbp,ubp,lbc,ubc,lbl,ubl,parr,larr, &
                    p2c_scale_type,c2l_scale_type)
    implicit none
    integer(ik4) , intent(in) :: lbp , ubp ! beginning and ending pft indices
    integer(ik4) , intent(in) :: lbc , ubc ! beginning and ending column indices
    integer(ik4) , intent(in) :: lbl , ubl ! beginning and ending ldunit indices
    real(rk8) , intent(in) , dimension(lbp:ubp) :: parr  ! input column array
    real(rk8) , intent(out) , dimension(lbl:ubl) :: larr ! output landunit array
    ! scale factor type for averaging
    character(len=*) , intent(in) :: p2c_scale_type
    ! scale factor type for averaging
    character(len=*) , intent(in) :: c2l_scale_type
    integer(ik4) :: p , c , l , iind       ! indices
    logical :: found                            ! temporary for error check
    real(rk8) , dimension(lbl:ubl) :: sumwt     ! sum of weights
    ! scale factor for pft->column mapping
    real(rk8) , dimension(lbc:ubc) :: scale_p2c
    ! scale factor for column->landunit mapping
    real(rk8) , dimension(lbc:ubc) :: scale_c2l
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer , dimension(:) :: pactive
    ! weight of pft relative to landunit
    real(rk8) , pointer , dimension(:) :: wtlunit
    ! column of corresponding pft
    integer(ik4) , pointer , dimension(:) :: pcolumn
    ! landunit of corresponding pft
    integer(ik4) , pointer , dimension(:) :: plandunit
    ! number of pfts in landunit
    integer(ik4) , pointer , dimension(:) :: npfts
    ! initial pft index in landunit
    integer(ik4) , pointer , dimension(:) :: pfti
    ! landunit of corresponding column
    integer(ik4) , pointer , dimension(:) :: clandunit
    integer(ik4) , pointer , dimension(:) :: ctype      ! column type
    integer(ik4) , pointer , dimension(:) :: ltype      ! landunit type
    ! urban canyon height to width ratio
    real(rk8) , pointer , dimension(:) :: canyon_hwr

    pactive    => clm3%g%l%c%p%active
    canyon_hwr => clm3%g%l%canyon_hwr
    ltype      => clm3%g%l%itype
    ctype      => clm3%g%l%c%itype
    clandunit  => clm3%g%l%c%landunit
    wtlunit    => clm3%g%l%c%p%wtlunit
    pcolumn    => clm3%g%l%c%p%column
    plandunit  => clm3%g%l%c%p%landunit
    npfts      => clm3%g%l%npfts
    pfti       => clm3%g%l%pfti

    if ( c2l_scale_type == 'unity' ) then
      do c = lbc , ubc
        scale_c2l(c) = 1.0D0
      end do
    else if ( c2l_scale_type == 'urbanf' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else if ( c2l_scale_type == 'urbans' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l))/(2.0D0*canyon_hwr(l)+1.0D0)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l))/(2.0D0*canyon_hwr(l)+1.0D0)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0 / (2.0D0*canyon_hwr(l) + 1.0D0)
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else
      write(stderr,*)'p2l_1d error: scale type ',c2l_scale_type,' not supported'
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if

    if ( p2c_scale_type == 'unity' ) then
      do p = lbp , ubp
        scale_p2c(p) = 1.0D0
      end do
    else
      write(stderr,*)'p2l_1d error: scale type ',p2c_scale_type,' not supported'
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if

    larr(:) = spval
    sumwt(:) = 0.D0
    do p = lbp,ubp
      if ( pactive(p) .and. wtlunit(p) /= 0.D0 ) then
        c = pcolumn(p)
        if ( parr(p) /= spval .and. scale_c2l(c) /= spval ) then
          l = plandunit(p)
          if (sumwt(l) == 0.D0) larr(l) = 0.D0
          larr(l) = larr(l) + parr(p) * scale_p2c(p) * scale_c2l(c) * wtlunit(p)
          sumwt(l) = sumwt(l) + wtlunit(p)
        end if
      end if
    end do
    found = .false.
    do l = lbl , ubl
      if ( sumwt(l) > 1.0D0 + 1.D-4 ) then
        found = .true.
        iind = l
      else if ( sumwt(l) /= 0.D0 ) then
        larr(l) = larr(l)/sumwt(l)
      end if
    end do
    if (found) then
      write(stderr,*)'p2l_1d error: sumwt is greater than 1.0 at l= ',iind
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if
  end subroutine p2l_1d
  !
  ! Perfrom subgrid-average from pfts to landunits
  ! Averaging is only done for points that are not equal to "spval".
  !
  subroutine p2l_2d(lbp,ubp,lbc,ubc,lbl,ubl,num2d,parr,larr, &
                    p2c_scale_type,c2l_scale_type)
    implicit none
    integer(ik4) , intent(in) :: lbp , ubp ! beginning and ending pft indices
    integer(ik4) , intent(in) :: lbc , ubc ! beginning and ending column indices
    integer(ik4) , intent(in) :: lbl , ubl ! beginning and ending ldunit indices
    integer(ik4) , intent(in)  :: num2d ! size of second dimension
    real(rk8) , intent(in) , dimension(lbp:ubp,num2d) :: parr ! input pft array
    ! output gridcell array
    real(rk8) , intent(out) , dimension(lbl:ubl,num2d) :: larr
    ! scale factor type for averaging
    character(len=*) , intent(in) :: p2c_scale_type
    ! scale factor type for averaging
    character(len=*) , intent(in) :: c2l_scale_type
    integer(ik4) :: j , p , c , l , iind ! indices
    logical :: found ! temporary for error check
    real(rk8) , dimension(lbl:ubl) :: sumwt ! sum of weights
    ! scale factor for pft->column mapping
    real(rk8) , dimension(lbc:ubc) :: scale_p2c
    ! scale factor for column->landunit mapping
    real(rk8) , dimension(lbc:ubc) :: scale_c2l
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer , dimension(:) :: pactive
    ! weight of pft relative to landunit
    real(rk8) , pointer , dimension(:) :: wtlunit
    ! column of corresponding pft
    integer(ik4) , pointer , dimension(:) :: pcolumn
    ! landunit of corresponding pft
    integer(ik4) , pointer , dimension(:) :: plandunit
    ! number of pfts in landunit
    integer(ik4) , pointer , dimension(:) :: npfts
    ! initial pft index in landunit
    integer(ik4) , pointer , dimension(:) :: pfti
    ! landunit of corresponding column
    integer(ik4) , pointer , dimension(:) :: clandunit
    integer(ik4) , pointer , dimension(:) :: ctype ! column type
    integer(ik4) , pointer , dimension(:) :: ltype ! landunit type
    ! urban canyon height to width ratio
    real(rk8) , pointer , dimension(:) :: canyon_hwr

    pactive    => clm3%g%l%c%p%active
    canyon_hwr => clm3%g%l%canyon_hwr
    ltype      => clm3%g%l%itype
    clandunit  => clm3%g%l%c%landunit
    ctype      => clm3%g%l%c%itype
    wtlunit   => clm3%g%l%c%p%wtlunit
    pcolumn   => clm3%g%l%c%p%column
    plandunit => clm3%g%l%c%p%landunit
    npfts     => clm3%g%l%npfts
    pfti      => clm3%g%l%pfti

    if ( c2l_scale_type == 'unity' ) then
      do c = lbc , ubc
        scale_c2l(c) = 1.0D0
      end do
    else if ( c2l_scale_type == 'urbanf' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else if ( c2l_scale_type == 'urbans' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l)) / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l)) / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0 / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else
      write(stderr,*)'p2l_2d error: scale type ',c2l_scale_type,' not supported'
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if

    if ( p2c_scale_type == 'unity' ) then
      do p = lbp , ubp
        scale_p2c(p) = 1.0D0
      end do
    else
      write(stderr,*)'p2l_2d error: scale type ',p2c_scale_type,' not supported'
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if

    larr(:,:) = spval
    do j = 1,num2d
      sumwt(:) = 0.D0
      do p = lbp , ubp
        if ( pactive(p) .and. wtlunit(p) /= 0.D0 ) then
          c = pcolumn(p)
          if ( parr(p,j) /= spval .and. scale_c2l(c) /= spval ) then
            l = plandunit(p)
            if ( sumwt(l) == 0.D0 ) larr(l,j) = 0.D0
            larr(l,j) = larr(l,j) + &
                        parr(p,j) * scale_p2c(p) * scale_c2l(c) * wtlunit(p)
            sumwt(l) = sumwt(l) + wtlunit(p)
          end if
        end if
      end do
      found = .false.
      do l = lbl , ubl
        if ( sumwt(l) > 1.0D0 + 1.D-4 ) then
          found = .true.
          iind = l
        else if ( sumwt(l) /= 0.D0 ) then
          larr(l,j) = larr(l,j)/sumwt(l)
        end if
      end do
      if (found) then
        write(stderr,*) &
            'p2l_2d error: sumwt is greater than 1.0 at l= ',iind,' j= ',j
        call fatal(__FILE__,__LINE__,'clm is now stopping')
      end if
    end do
  end subroutine p2l_2d
  !
  ! Perfrom subgrid-average from pfts to gridcells.
  ! Averaging is only done for points that are not equal to "spval".
  !
  subroutine p2g_1d(lbp,ubp,lbc,ubc,lbl,ubl,lbg,ubg,parr,garr, &
                    p2c_scale_type,c2l_scale_type,l2g_scale_type)
    implicit none
    integer(ik4) , intent(in) :: lbp , ubp ! beginning and ending pft indices
    integer(ik4) , intent(in) :: lbc , ubc ! beginning and ending column indices
    integer(ik4) , intent(in) :: lbl , ubl ! beginning and ending ldunit indices
    integer(ik4) , intent(in) :: lbg , ubg ! beginning and ending gdcell indices
    real(rk8) , intent(in) , dimension(lbp:ubp) :: parr  ! input pft array
    real(rk8) , intent(out) , dimension(lbg:ubg) :: garr ! output gridcell array
    ! scale factor type for averaging
    character(len=*) , intent(in) :: p2c_scale_type
    ! scale factor type for averaging
    character(len=*) , intent(in) :: c2l_scale_type
    ! scale factor type for averaging
    character(len=*) , intent(in) :: l2g_scale_type
    integer(ik4) :: p , c , l , g , iind       ! indices
    logical :: found                  ! temporary for error check
    real(rk8) , dimension(lbp:ubp) :: scale_p2c    ! scale factor
    real(rk8) , dimension(lbc:ubc) :: scale_c2l    ! scale factor
    real(rk8) , dimension(lbl:ubl) :: scale_l2g    ! scale factor
    real(rk8) , dimension(lbg:ubg) :: sumwt        ! sum of weights
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer , dimension(:) :: pactive
    ! weight of pfts relative to gridcells
    real(rk8) , pointer , dimension(:) :: wtgcell
    ! column of corresponding pft
    integer(ik4) , pointer , dimension(:) :: pcolumn
    ! landunit of corresponding pft
    integer(ik4) , pointer , dimension(:) :: plandunit
    ! gridcell of corresponding pft
    integer(ik4) , pointer , dimension(:) :: pgridcell
    ! number of pfts in gridcell
    integer(ik4) , pointer , dimension(:) :: npfts
    ! initial pft index in gridcell
    integer(ik4) , pointer , dimension(:) :: pfti
    ! column type
    integer(ik4) , pointer , dimension(:) :: ctype
    ! landunit of corresponding column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! landunit type
    integer(ik4) , pointer , dimension(:) :: ltype
    ! urban canyon height to width ratio
    real(rk8) , pointer , dimension(:) :: canyon_hwr

    pactive    => clm3%g%l%c%p%active
    canyon_hwr => clm3%g%l%canyon_hwr
    ltype      => clm3%g%l%itype
    clandunit  => clm3%g%l%c%landunit
    ctype      => clm3%g%l%c%itype
    wtgcell   => clm3%g%l%c%p%wtgcell
    pcolumn   => clm3%g%l%c%p%column
    pgridcell => clm3%g%l%c%p%gridcell
    plandunit => clm3%g%l%c%p%landunit
    npfts     => clm3%g%npfts
    pfti      => clm3%g%pfti

    call build_scale_l2g(l2g_scale_type,lbl,ubl,scale_l2g)

    if ( c2l_scale_type == 'unity' ) then
      do c = lbc , ubc
        scale_c2l(c) = 1.0D0
      end do
    else if ( c2l_scale_type == 'urbanf' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else if ( c2l_scale_type == 'urbans' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l)) / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l)) / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0 / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else
      write(stderr,*)'p2g_1d error: scale type ',c2l_scale_type,' not supported'
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if

    if ( p2c_scale_type == 'unity' ) then
      do p = lbp , ubp
        scale_p2c(p) = 1.0D0
      end do
    else
      write(stderr,*)'p2g_1d error: scale type ',c2l_scale_type,' not supported'
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if

    garr(:) = spval
    sumwt(:) = 0.D0
    do p = lbp , ubp
      if ( pactive(p) .and. wtgcell(p) /= 0.D0 ) then
        c = pcolumn(p)
        l = plandunit(p)
        if ( parr(p) /= spval .and. &
             scale_c2l(c) /= spval .and. scale_l2g(l) /= spval ) then
          g = pgridcell(p)
          if ( sumwt(g) == 0.D0 ) garr(g) = 0.D0
          garr(g) = garr(g) + &
              parr(p) * scale_p2c(p) * scale_c2l(c) * scale_l2g(l) * wtgcell(p)
          sumwt(g) = sumwt(g) + wtgcell(p)
        end if
      end if
    end do
    found = .false.
    do g = lbg , ubg
      if ( sumwt(g) > 1.0D0 + 1.D-4 ) then
        found = .true.
        iind = g
      else if ( sumwt(g) /= 0.D0 ) then
        garr(g) = garr(g)/sumwt(g)
      end if
    end do
    if ( found ) then
      write(stderr,*)'p2g_1d error: sumwt is greater than 1.0 at g= ',iind
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if
  end subroutine p2g_1d
  !
  ! Perfrom subgrid-average from pfts to gridcells.
  ! Averaging is only done for points that are not equal to "spval".
  !
  subroutine p2g_2d(lbp,ubp,lbc,ubc,lbl,ubl,lbg,ubg,num2d,parr,garr, &
                    p2c_scale_type,c2l_scale_type,l2g_scale_type)
    implicit none
    integer(ik4) , intent(in) :: lbp , ubp ! beginning and ending pft indices
    integer(ik4) , intent(in) :: lbc , ubc ! beginning and ending column indices
    integer(ik4) , intent(in) :: lbl , ubl ! beginning and ending ldunit indices
    integer(ik4) , intent(in) :: lbg , ubg ! beginning and ending gdcell indices
    integer(ik4) , intent(in) :: num2d ! size of second dimension
    real(rk8) , intent(in) , dimension(lbp:ubp,num2d) :: parr  ! input pft array
    ! output gridcell array
    real(rk8) , intent(out) , dimension(lbg:ubg,num2d) :: garr
    ! scale factor type for averaging
    character(len=*) , intent(in) :: p2c_scale_type
    ! scale factor type for averaging
    character(len=*) , intent(in) :: c2l_scale_type
    ! scale factor type for averaging
    character(len=*) , intent(in) :: l2g_scale_type
    integer(ik4) :: j , p , c , l , g , iind     ! indices
    logical :: found                  ! temporary for error check
    real(rk8) , dimension(lbp:ubp) :: scale_p2c   ! scale factor
    real(rk8) , dimension(lbc:ubc) :: scale_c2l   ! scale factor
    real(rk8) , dimension(lbl:ubl) :: scale_l2g   ! scale factor
    real(rk8) , dimension(lbg:ubg) :: sumwt       ! sum of weights
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer , dimension(:) :: pactive
    ! weight of pfts relative to gridcells
    real(rk8) , pointer , dimension(:) :: wtgcell
    ! column of corresponding pft
    integer(ik4) , pointer , dimension(:) :: pcolumn
    ! landunit of corresponding pft
    integer(ik4) , pointer , dimension(:) :: plandunit
    ! gridcell of corresponding pft
    integer(ik4) , pointer , dimension(:) :: pgridcell
    ! number of pfts in gridcell
    integer(ik4) , pointer , dimension(:) :: npfts
    ! initial pft index in gridcell
    integer(ik4) , pointer , dimension(:) :: pfti
    ! landunit of corresponding column
    integer(ik4) , pointer , dimension(:) :: clandunit
    integer(ik4) , pointer , dimension(:) :: ctype ! column type
    integer(ik4) , pointer , dimension(:) :: ltype ! landunit type
    ! urban canyon height to width ratio
    real(rk8) , pointer , dimension(:) :: canyon_hwr

    pactive      => clm3%g%l%c%p%active
    canyon_hwr   => clm3%g%l%canyon_hwr
    ltype        => clm3%g%l%itype
    clandunit    => clm3%g%l%c%landunit
    ctype        => clm3%g%l%c%itype
    wtgcell      => clm3%g%l%c%p%wtgcell
    pcolumn      => clm3%g%l%c%p%column
    pgridcell    => clm3%g%l%c%p%gridcell
    plandunit    => clm3%g%l%c%p%landunit
    npfts        => clm3%g%npfts
    pfti         => clm3%g%pfti

    call build_scale_l2g(l2g_scale_type,lbl,ubl,scale_l2g)

    if ( c2l_scale_type == 'unity' ) then
      do c = lbc , ubc
        scale_c2l(c) = 1.0D0
      end do
    else if ( c2l_scale_type == 'urbanf' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else if ( c2l_scale_type == 'urbans' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l)) / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l)) / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0 / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else
      write(stderr,*)'p2g_2d error: scale type ',c2l_scale_type,' not supported'
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if

    if ( p2c_scale_type == 'unity' ) then
      do p = lbp , ubp
        scale_p2c(p) = 1.0D0
      end do
    else
      write(stderr,*)'p2g_2d error: scale type ',c2l_scale_type,' not supported'
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if

    garr(:,:) = spval
    do j = 1 , num2d
      sumwt(:) = 0.D0
      do p = lbp , ubp
        if ( pactive(p) .and. wtgcell(p) /= 0.D0 ) then
          c = pcolumn(p)
          l = plandunit(p)
          if ( parr(p,j) /= spval .and. &
               scale_c2l(c) /= spval .and. scale_l2g(l) /= spval ) then
            g = pgridcell(p)
            if ( sumwt(g) == 0.D0 ) garr(g,j) = 0.D0
            garr(g,j) = garr(g,j) + &
                 parr(p,j)*scale_p2c(p)*scale_c2l(c)*scale_l2g(l)*wtgcell(p)
            sumwt(g) = sumwt(g) + wtgcell(p)
          end if
        end if
      end do
      found = .false.
      do g = lbg , ubg
        if ( sumwt(g) > 1.0D0 + 1.D-4 ) then
          found = .true.
          iind = g
        else if ( sumwt(g) /= 0.D0 ) then
          garr(g,j) = garr(g,j)/sumwt(g)
        end if
      end do
      if (found) then
        write(stderr,*) &
             'p2g_2d error: sumwt gt 1.0 at g/sumwt = ',iind,sumwt(iind)
        call fatal(__FILE__,__LINE__,'clm is now stopping')
      end if
    end do
  end subroutine p2g_2d
  !
  ! Perfrom subgrid-average from columns to landunits
  ! Averaging is only done for points that are not equal to "spval".
  !
  subroutine c2l_1d(lbc,ubc,lbl,ubl,carr,larr,c2l_scale_type)
    implicit none
    integer(ik4) , intent(in) :: lbc , ubc ! beginning and ending column indices
    integer(ik4) , intent(in) :: lbl , ubl ! beginning and ending ldunit indices
    real(rk8) , intent(in) , dimension(lbc:ubc) :: carr  ! input column array
    real(rk8) , intent(out) , dimension(lbl:ubl) :: larr ! output landunit array
    ! scale factor type for averaging
    character(len=*) , intent(in) :: c2l_scale_type
    integer(ik4) :: c , l , iind ! indices
    logical  :: found  ! temporary for error check
    ! scale factor for column->landunit mapping
    real(rk8) , dimension(lbc:ubc) :: scale_c2l
    real(rk8) , dimension(lbl:ubl) :: sumwt ! sum of weights
    ! true=>do computations on this column (see reweightMod for details)
    logical , pointer , dimension(:) :: cactive
    ! weight of landunits relative to gridcells
    real(rk8) , pointer , dimension(:) :: wtlunit
    ! gridcell of corresponding column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! number of columns in landunit
    integer(ik4) , pointer , dimension(:) :: ncolumns
    ! initial column index in landunit
    integer(ik4) , pointer , dimension(:) :: coli
    integer(ik4) , pointer , dimension(:) :: ctype ! column type
    integer(ik4) , pointer , dimension(:) :: ltype ! landunit type
    ! urban canyon height to width ratio
    real(rk8) , pointer , dimension(:) :: canyon_hwr

    cactive    => clm3%g%l%c%active
    ctype      => clm3%g%l%c%itype
    ltype      => clm3%g%l%itype
    canyon_hwr => clm3%g%l%canyon_hwr
    wtlunit    => clm3%g%l%c%wtlunit
    clandunit  => clm3%g%l%c%landunit
    ncolumns   => clm3%g%l%ncolumns
    coli       => clm3%g%l%coli

    if ( c2l_scale_type == 'unity' ) then
      do c = lbc , ubc
        scale_c2l(c) = 1.0D0
      end do
    else if ( c2l_scale_type == 'urbanf' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else if ( c2l_scale_type == 'urbans' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l)) / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l)) / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0 / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else
      write(stderr,*)'c2l_1d error: scale type ',c2l_scale_type,' not supported'
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if

    larr(:) = spval
    sumwt(:) = 0.D0
    do c = lbc , ubc
      if ( cactive(c) .and. wtlunit(c) /= 0.D0 ) then
        if ( carr(c) /= spval .and. scale_c2l(c) /= spval ) then
          l = clandunit(c)
          if (sumwt(l) == 0.D0) larr(l) = 0.D0
          larr(l) = larr(l) + carr(c) * scale_c2l(c) * wtlunit(c)
          sumwt(l) = sumwt(l) + wtlunit(c)
        end if
      end if
    end do
    found = .false.
    do l = lbl , ubl
      if ( sumwt(l) > 1.0D0 + 1.D-4 ) then
        found = .true.
        iind = l
      else if ( sumwt(l) /= 0.D0 ) then
        larr(l) = larr(l)/sumwt(l)
      end if
    end do
    if (found) then
      write(stderr,*)'c2l_1d error: sumwt is greater than 1.0 at l= ',iind
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if
  end subroutine c2l_1d
  !
  ! Perfrom subgrid-average from columns to landunits
  ! Averaging is only done for points that are not equal to "spval".
  !
  subroutine c2l_2d(lbc,ubc,lbl,ubl,num2d,carr,larr,c2l_scale_type)
    implicit none
    integer(ik4) , intent(in) :: lbc , ubc ! beginning and ending column indices
    integer(ik4) , intent(in) :: lbl , ubl ! beginning and ending ldunit indices
    integer(ik4) , intent(in) :: num2d     ! size of second dimension
    ! input column array
    real(rk8) , intent(in) , dimension(lbc:ubc,num2d) :: carr
    ! output landunit array
    real(rk8) , intent(out) , dimension(lbl:ubl,num2d) :: larr
    ! scale factor type for averaging
    character(len=*) , intent(in) :: c2l_scale_type
    integer(ik4) :: j , l , c , iind         ! indices
    logical :: found ! temporary for error check
    ! scale factor for column->landunit mapping
    real(rk8) , dimension(lbc:ubc) :: scale_c2l
    real(rk8) , dimension(lbl:ubl) :: sumwt ! sum of weights
    ! true=>do computations on this column (see reweightMod for details)
    logical , pointer , dimension(:) :: cactive
    ! weight of column relative to landunit
    real(rk8) , pointer , dimension(:) :: wtlunit
    ! landunit of corresponding column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! number of columns in landunit
    integer(ik4) , pointer , dimension(:) :: ncolumns
    ! initial column index in landunit
    integer(ik4) , pointer , dimension(:) :: coli
    integer(ik4) , pointer , dimension(:) :: ctype ! column type
    integer(ik4) , pointer , dimension(:) :: ltype ! landunit type
    ! urban canyon height to width ratio
    real(rk8) , pointer , dimension(:) :: canyon_hwr

    cactive    => clm3%g%l%c%active
    ctype      => clm3%g%l%c%itype
    ltype      => clm3%g%l%itype
    canyon_hwr => clm3%g%l%canyon_hwr
    wtlunit    => clm3%g%l%c%wtlunit
    clandunit  => clm3%g%l%c%landunit
    ncolumns   => clm3%g%l%ncolumns
    coli       => clm3%g%l%coli

    if ( c2l_scale_type == 'unity' ) then
      do c = lbc , ubc
        scale_c2l(c) = 1.0D0
      end do
    else if ( c2l_scale_type == 'urbanf' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else if ( c2l_scale_type == 'urbans' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l)) / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l)) / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0 / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else
      write(stderr,*)'c2l_2d error: scale type ',c2l_scale_type,' not supported'
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if

    larr(:,:) = spval
    do j = 1 , num2d
      sumwt(:) = 0.D0
      do c = lbc , ubc
        if ( cactive(c) .and. wtlunit(c) /= 0.D0 ) then
          if ( carr(c,j) /= spval .and. scale_c2l(c) /= spval ) then
            l = clandunit(c)
            if ( sumwt(l) == 0.D0 ) larr(l,j) = 0.D0
            larr(l,j) = larr(l,j) + carr(c,j) * scale_c2l(c) * wtlunit(c)
            sumwt(l) = sumwt(l) + wtlunit(c)
          end if
        end if
      end do
      found = .false.
      do l = lbl , ubl
        if ( sumwt(l) > 1.0D0 + 1.D-4 ) then
          found = .true.
          iind = l
        else if ( sumwt(l) /= 0.D0 ) then
          larr(l,j) = larr(l,j)/sumwt(l)
        end if
      end do
      if (found) then
        write(stderr,*) &
            'c2l_2d error: sumwt is greater than 1.0 at l= ',iind,' lev= ',j
        call fatal(__FILE__,__LINE__,'clm is now stopping')
      end if
    end do
  end subroutine c2l_2d
  !
  ! Perfrom subgrid-average from columns to gridcells.
  ! Averaging is only done for points that are not equal to "spval".
  !
  subroutine c2g_1d(lbc,ubc,lbl,ubl,lbg,ubg,carr,garr, &
                    c2l_scale_type,l2g_scale_type)
    implicit none
    integer(ik4) , intent(in) :: lbc , ubc ! beginning and ending column indices
    integer(ik4) , intent(in) :: lbl , ubl ! beginning and ending ldunit indices
    integer(ik4) , intent(in) :: lbg , ubg ! beginning and ending ldunit indices
    real(rk8) , intent(in) , dimension(lbc:ubc) :: carr  ! input column array
    real(rk8) , intent(out) , dimension(lbg:ubg) :: garr ! output gridcell array
    ! scale factor type for averaging
    character(len=*) , intent(in) :: c2l_scale_type
    ! scale factor type for averaging
    character(len=*) , intent(in) :: l2g_scale_type
    integer(ik4) :: c , l , g , iind  ! indices
    logical :: found ! temporary for error check
    real(rk8) , dimension(lbc:ubc) :: scale_c2l    ! scale factor
    real(rk8) , dimension(lbl:ubl) :: scale_l2g    ! scale factor
    real(rk8) , dimension(lbg:ubg) :: sumwt        ! sum of weights
    ! true=>do computations on this column (see reweightMod for details)
    logical , pointer , dimension(:) :: cactive
    ! weight of columns relative to gridcells
    real(rk8) , pointer , dimension(:) :: wtgcell
    ! landunit of corresponding column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! gridcell of corresponding column
    integer(ik4) , pointer , dimension(:) :: cgridcell
    ! number of columns in gridcell
    integer(ik4) , pointer , dimension(:) :: ncolumns
    ! initial column index in gridcell
    integer(ik4) , pointer , dimension(:) :: coli
    integer(ik4) , pointer , dimension(:) :: ctype ! column type
    integer(ik4) , pointer , dimension(:) :: ltype ! landunit type
    ! urban canyon height to width ratio
    real(rk8) , pointer , dimension(:) :: canyon_hwr

    cactive    => clm3%g%l%c%active
    ctype      => clm3%g%l%c%itype
    ltype      => clm3%g%l%itype
    canyon_hwr => clm3%g%l%canyon_hwr
    wtgcell    => clm3%g%l%c%wtgcell
    clandunit  => clm3%g%l%c%landunit
    cgridcell  => clm3%g%l%c%gridcell
    ncolumns   => clm3%g%ncolumns
    coli       => clm3%g%coli

    call build_scale_l2g(l2g_scale_type,lbl,ubl,scale_l2g)

    if ( c2l_scale_type == 'unity' ) then
      do c = lbc , ubc
        scale_c2l(c) = 1.0D0
      end do
    else if ( c2l_scale_type == 'urbanf' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else if ( c2l_scale_type == 'urbans' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l)) / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l)) / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0 / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else
     write(stderr,*)'c2l_1d error: scale type ',c2l_scale_type,' not supported'
     call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if

    garr(:) = spval
    sumwt(:) = 0.D0
    do c = lbc , ubc
      if ( cactive(c) .and. wtgcell(c) /= 0.D0 ) then
        l = clandunit(c)
        if ( carr(c) /= spval .and. &
             scale_c2l(c) /= spval .and. scale_l2g(l) /= spval ) then
          g = cgridcell(c)
          if ( sumwt(g) == 0.D0 ) garr(g) = 0.D0
          garr(g) = garr(g) + carr(c) * scale_c2l(c) * scale_l2g(l) * wtgcell(c)
          sumwt(g) = sumwt(g) + wtgcell(c)
        end if
      end if
    end do
    found = .false.
    do g = lbg , ubg
      if ( sumwt(g) > 1.0D0 + 1.D-4 ) then
        found = .true.
        iind = g
      else if ( sumwt(g) /= 0.D0 ) then
        garr(g) = garr(g)/sumwt(g)
      end if
    end do
    if (found) then
      write(stderr,*)'c2g_1d error: sumwt is greater than 1.0 at g= ',iind
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if
  end subroutine c2g_1d
  !
  ! Perfrom subgrid-average from columns to gridcells.
  ! Averaging is only done for points that are not equal to "spval".
  !
  subroutine c2g_2d(lbc,ubc,lbl,ubl,lbg,ubg,num2d,carr,garr, &
                    c2l_scale_type,l2g_scale_type)
    implicit none
    integer(ik4) , intent(in) :: lbc , ubc ! beginning and ending column indices
    integer(ik4) , intent(in) :: lbl , ubl ! beginning and ending launit indices
    integer(ik4) , intent(in) :: lbg , ubg ! beginning and ending gdcell indices
    integer(ik4) , intent(in) :: num2d ! size of second dimension
    ! input column array
    real(rk8) , intent(in) , dimension(lbc:ubc,num2d) :: carr
    ! output gridcell array
    real(rk8) , intent(out) , dimension(lbg:ubg,num2d) :: garr
    ! scale factor type for averaging
    character(len=*) , intent(in) :: c2l_scale_type
    ! scale factor type for averaging
    character(len=*) , intent(in) :: l2g_scale_type
    integer(ik4) :: j , c , g , l , iind  ! indices
    logical :: found                  ! temporary for error check
    real(rk8) , dimension(lbc:ubc) :: scale_c2l  ! scale factor
    real(rk8) , dimension(lbl:ubl) :: scale_l2g  ! scale factor
    real(rk8) , dimension(lbg:ubg) :: sumwt      ! sum of weights
    ! true=>do computations on this column (see reweightMod for details)
    logical , pointer , dimension(:) :: cactive
    ! weight of columns relative to gridcells
    real(rk8) , pointer , dimension(:) :: wtgcell
    ! landunit of corresponding column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! gridcell of corresponding column
    integer(ik4) , pointer , dimension(:) :: cgridcell
    ! number of columns in gridcell
    integer(ik4) , pointer , dimension(:) :: ncolumns
    ! initial column index in gridcell
    integer(ik4) , pointer , dimension(:) :: coli
    integer(ik4) , pointer , dimension(:) :: ctype ! column type
    integer(ik4) , pointer , dimension(:) :: ltype ! landunit type
    ! urban canyon height to width ratio
    real(rk8) , pointer , dimension(:) :: canyon_hwr

    cactive    => clm3%g%l%c%active
    ctype      => clm3%g%l%c%itype
    ltype      => clm3%g%l%itype
    canyon_hwr => clm3%g%l%canyon_hwr
    wtgcell    => clm3%g%l%c%wtgcell
    clandunit  => clm3%g%l%c%landunit
    cgridcell  => clm3%g%l%c%gridcell
    ncolumns   => clm3%g%ncolumns
    coli       => clm3%g%coli

    call build_scale_l2g(l2g_scale_type,lbl,ubl,scale_l2g)

    if ( c2l_scale_type == 'unity' ) then
      do c = lbc , ubc
        scale_c2l(c) = 1.0D0
      end do
    else if ( c2l_scale_type == 'urbanf' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = 3.0D0 * canyon_hwr(l)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else if ( c2l_scale_type == 'urbans' ) then
      do c = lbc , ubc
        l = clandunit(c)
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_sunwall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l)) / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_shadewall ) then
            scale_c2l(c) = (3.0D0 * canyon_hwr(l)) / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_road_perv .or. &
                    ctype(c) == icol_road_imperv ) then
            scale_c2l(c) = 3.0D0 / (2.D0*canyon_hwr(l) + 1.D0)
          else if ( ctype(c) == icol_roof ) then
            scale_c2l(c) = 1.0D0
          end if
        else
          scale_c2l(c) = 1.0D0
        end if
      end do
    else
      write(stderr,*)'c2g_2d error: scale type ',c2l_scale_type,' not supported'
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if

    garr(:,:) = spval
    do j = 1 , num2d
      sumwt(:) = 0.D0
      do c = lbc , ubc
        if ( cactive(c) .and. wtgcell(c) /= 0.D0 ) then
          l = clandunit(c)
          if ( carr(c,j) /= spval .and. &
               scale_c2l(c) /= spval .and. scale_l2g(l) /= spval) then
            g = cgridcell(c)
            if ( sumwt(g) == 0.D0) garr(g,j ) = 0.D0
            garr(g,j) = garr(g,j) + &
               carr(c,j) * scale_c2l(c) * scale_l2g(l) * wtgcell(c)
            sumwt(g) = sumwt(g) + wtgcell(c)
          end if
        end if
      end do
      found = .false.
      do g = lbg, ubg
        if ( sumwt(g) > 1.0D0 + 1.D-4 ) then
          found = .true.
          iind = g
        else if ( sumwt(g) /= 0.D0 ) then
          garr(g,j) = garr(g,j)/sumwt(g)
        end if
      end do
      if ( found ) then
        write(stderr,*)'c2g_2d error: sumwt is greater than 1.0 at g= ',iind
        call fatal(__FILE__,__LINE__,'clm is now stopping')
      end if
    end do
  end subroutine c2g_2d
  !
  ! Perfrom subgrid-average from landunits to gridcells.
  ! Averaging is only done for points that are not equal to "spval".
  !
  subroutine l2g_1d(lbl,ubl,lbg,ubg,larr,garr,l2g_scale_type)
    implicit none
    integer(ik4) , intent(in) :: lbl , ubl ! sub landunit indices
    integer(ik4) , intent(in) :: lbg , ubg ! gridcell indices
    real(rk8) , intent(in) , dimension(lbl:ubl) :: larr  ! input landunit array
    real(rk8) , intent(out) , dimension(lbg:ubg) :: garr ! output gridcell array
    ! scale factor type for averaging
    character(len=*) , intent(in) :: l2g_scale_type
    integer(ik4) :: l , g , iind ! indices
    logical :: found ! temporary for error check
    real(rk8) , dimension(lbl:ubl) :: scale_l2g     ! scale factor
    real(rk8) , dimension(lbg:ubg) :: sumwt         ! sum of weights
    ! true=>do computations on this landunit (see reweightMod for details)
    logical , pointer , dimension(:) :: lactive
    ! weight of landunits relative to gridcells
    real(rk8) , pointer , dimension(:) :: wtgcell
    ! gridcell of corresponding landunit
    integer(ik4) , pointer , dimension(:) :: lgridcell
    ! number of landunits in gridcell
    integer(ik4) , pointer , dimension(:) :: nlandunits
    ! initial landunit index in gridcell
    integer(ik4) , pointer , dimension(:) :: luni

    lactive    => clm3%g%l%active
    wtgcell    => clm3%g%l%wtgcell
    lgridcell  => clm3%g%l%gridcell
    nlandunits => clm3%g%nlandunits
    luni       => clm3%g%luni

    call build_scale_l2g(l2g_scale_type,lbl,ubl,scale_l2g)

    garr(:) = spval
    sumwt(:) = 0.D0
    do l = lbl , ubl
      if ( lactive(l) .and. wtgcell(l) /= 0.D0 ) then
        if ( larr(l) /= spval .and. scale_l2g(l) /= spval ) then
          g = lgridcell(l)
          if ( sumwt(g) == 0.D0 ) garr(g) = 0.D0
          garr(g) = garr(g) + larr(l) * scale_l2g(l) * wtgcell(l)
          sumwt(g) = sumwt(g) + wtgcell(l)
        end if
      end if
    end do
    found = .false.
    do g = lbg , ubg
      if ( sumwt(g) > 1.0D0 + 1.D-4 ) then
        found = .true.
        iind = g
      else if ( sumwt(g) /= 0.D0 ) then
        garr(g) = garr(g)/sumwt(g)
      end if
    end do
    if ( found ) then
      write(stderr,*)'l2g_1d error: sumwt is greater than 1.0 at g= ',iind
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if
  end subroutine l2g_1d
  !
  ! Perfrom subgrid-average from landunits to gridcells.
  ! Averaging is only done for points that are not equal to "spval".
  !
  subroutine l2g_2d(lbl,ubl,lbg,ubg,num2d,larr,garr,l2g_scale_type)
    implicit none
    integer(ik4) , intent(in) :: lbl , ubl ! beginning and ending column indices
    integer(ik4) , intent(in) :: lbg , ubg ! beginning and ending gdcell indices
    integer(ik4) , intent(in) :: num2d     ! size of second dimension
    ! input landunit array
    real(rk8) , intent(in) , dimension(lbl:ubl,num2d) :: larr
    ! output gridcell array
    real(rk8) , intent(out) , dimension(lbg:ubg,num2d) :: garr
    ! scale factor type for averaging
    character(len=*) , intent(in) :: l2g_scale_type
    integer(ik4) :: j , g , l , iind ! indices
    logical :: found ! temporary for error check
    real(rk8) , dimension(lbl:ubl) :: scale_l2g     ! scale factor
    real(rk8) , dimension(lbg:ubg) :: sumwt         ! sum of weights
    ! true=>do computations on this landunit (see reweightMod for details)
    logical , pointer , dimension(:) :: lactive
    ! weight of landunits relative to gridcells
    real(rk8) , pointer , dimension(:) :: wtgcell
    ! gridcell of corresponding landunit
    integer(ik4) , pointer , dimension(:) :: lgridcell
    ! number of landunits in gridcell
    integer(ik4) , pointer , dimension(:) :: nlandunits
    ! initial landunit index in gridcell
    integer(ik4) , pointer , dimension(:) :: luni

    lactive   => clm3%g%l%active
    wtgcell   => clm3%g%l%wtgcell
    lgridcell => clm3%g%l%gridcell
    nlandunits => clm3%g%nlandunits
    luni       => clm3%g%luni

    call build_scale_l2g(l2g_scale_type,lbl,ubl,scale_l2g)

    garr(:,:) = spval
    do j = 1 , num2d
      sumwt(:) = 0.D0
      do l = lbl , ubl
        if ( lactive(l) .and. wtgcell(l) /= 0.D0 ) then
          if ( larr(l,j) /= spval .and. scale_l2g(l) /= spval ) then
            g = lgridcell(l)
            if ( sumwt(g) == 0.D0 ) garr(g,j) = 0.D0
            garr(g,j) = garr(g,j) + larr(l,j) * scale_l2g(l) * wtgcell(l)
            sumwt(g) = sumwt(g) + wtgcell(l)
          end if
        end if
      end do
      found = .false.
      do g = lbg , ubg
        if ( sumwt(g) > 1.0D0 + 1.D-4 ) then
          found = .true.
          iind= g
        else if ( sumwt(g) /= 0.D0 ) then
          garr(g,j) = garr(g,j)/sumwt(g)
        end if
      end do
      if ( found ) then
        write(stderr,*) &
           'l2g_2d error: sumwt is greater than 1.0 at g= ',iind,' lev= ',j
        call fatal(__FILE__,__LINE__,'clm is now stopping')
      end if
    end do
  end subroutine l2g_2d
  !
  ! Fill the scale_l2g(lbl:ubl) array with appropriate values for the
  ! given l2g_scale_type.
  ! This array can later be used to scale each landunit in forming
  ! grid cell averages.
  !
  subroutine build_scale_l2g(l2g_scale_type,lbl,ubl,scale_l2g)
     implicit none
     ! scale factor type for averaging
     character(len=*) , intent(in) :: l2g_scale_type
     integer(ik4) , intent(in) :: lbl , ubl ! column indices
     real(rk8) , intent(out) , dimension(lbl:ubl) :: scale_l2g ! scale factor
     ! scale factor for each landunit type
     real(rk8) , dimension(max_lunit) :: scale_lookup
     integer(ik4) :: l                       ! index
     integer(ik4) , pointer , dimension(:) :: ltype ! landunit type

     ltype      => clm3%g%l%itype

     call create_scale_l2g_lookup(l2g_scale_type,scale_lookup)

     do l = lbl , ubl
       scale_l2g(l) = scale_lookup(ltype(l))
     end do
  end subroutine build_scale_l2g
  !
  ! Create a lookup array, scale_lookup(1..max_lunit), which gives the
  ! scale factor for each landunit type depending on l2g_scale_type
  !
  subroutine create_scale_l2g_lookup(l2g_scale_type,scale_lookup)
    implicit none
    ! scale factor type for averaging
    character(len=*) , intent(in) :: l2g_scale_type
    ! scale factor for each landunit type
    real(rk8) , dimension(max_lunit) , intent(out) :: scale_lookup

    ! ------------ WJS (10-14-11): IMPORTANT GENERAL NOTES ------------
    !
    ! Since scale_l2g is not currently included in the sumwt accumulations,
    ! you need to be careful about the scale values you use. Values of 1 and
    ! spval are safe (including having multiple landunits with value 1), but
    ! only use other values if you know what you are doing!
    ! For example, using a value of 0 is NOT the correct way to exclude a
    ! landunit from the average, because the normalization will be done
    ! incorrectly in this case: instead, use spval to exclude a landunit from
    ! the average. Similarly, using a value of 2 is NOT the correct way to
    ! give a landunit double relative weight in general, because the
    ! normalization won't be done correctly in this case, either.
    !
    ! In the longer-term, I believe that the correct solution to this problem
    ! is to include scale_l2g (and the other scale factors) in the sumwt
    ! accumulations (e.g., sumwt = sumwt + wtgcell * scale_p2c * scale_c2l
    ! * scale_l2g), but that requires some more thought to (1) make sure
    ! that is correct, and (2) make sure it doesn't break the urban scaling.
    !
    ! -----------------------------------------------------------------

    ! Initialize scale_lookup to spval for all landunits.
    ! Thus, any landunit that keeps the default value will be excluded
    ! from grid cell averages.
    scale_lookup(:) = spval

    if ( l2g_scale_type == 'unity' ) then
      scale_lookup(:) = 1.0D0
    else if ( l2g_scale_type == 'veg' ) then
      scale_lookup(istsoil) = 1.0D0
      scale_lookup(istcrop) = 1.0D0
    else if ( l2g_scale_type == 'ice' ) then
      scale_lookup(istice) = 1.0D0
    else if ( l2g_scale_type == 'nonurb' ) then
      scale_lookup(:) = 1.0D0
      scale_lookup(isturb) = spval
    else if ( l2g_scale_type == 'lake' ) then
      scale_lookup(istdlak) = 1.0D0
    else
      write(stderr,*) &
          'scale_l2g_lookup_array error: scale type ',l2g_scale_type, &
          ' not supported'
      call fatal(__FILE__,__LINE__,'clm is now stopping')
    end if
  end subroutine create_scale_l2g_lookup

end module mod_clm_subgridave
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
