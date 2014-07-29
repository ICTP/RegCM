module mod_clm_initgridcells
  !
  ! Initializes sub-grid mapping for each land grid cell
  !
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_mpmessage
  use mod_dynparam
  use mod_mppparam
  use mod_clm_varsur , only : wtxy, vegxy
  use mod_clm_varcon , only : rpi

  implicit none

  private

  save

  public :: initGridcells ! initialize sub-grid gridcell mapping

  contains
  !
  ! Initialize sub-grid mapping and allocates space for derived type hierarchy.
  ! For each land gridcell determine landunit, column and pft properties.
  !
  subroutine initGridcells ()
    use mod_clm_type , only : clm3 , gridcell_type , landunit_type , &
                             column_type , pft_type
    use mod_clm_domain , only : ldomain
    use mod_clm_decomp , only : get_proc_global , get_proc_bounds
    use mod_clm_varcon , only : istsoil , istice , istwet , istdlak , &
            isturb , udens_tbd , udens_hd , udens_md
    use mod_clm_varcon , only : istcrop
    use mod_clm_subgrid , only : subgrid_get_gcellinfo
    use mod_clm_surfrd , only : crop_prog

    implicit none
    integer(ik4) :: li , ci , pi , gdc ! indices
    integer(ik4) :: ltype ! landunit type
    integer(ik4) :: numg    ! total number of gridcells across all processors
    integer(ik4) :: numl    ! total number of landunits across all processors
    integer(ik4) :: numc    ! total number of columns across all processors
    integer(ik4) :: nump    ! total number of pfts across all processors
    integer(ik4) :: begg , endg  ! local beg/end gridcells gdc
    integer(ik4) :: begl , endl  ! local beg/end landunits
    integer(ik4) :: begc , endc  ! local beg/end columns
    integer(ik4) :: begp , endp  ! local beg/end pfts
    logical :: my_gcell   ! is gdc gridcell on my pe
    integer(ik4) :: nwtxy ! wtxy cell index

    type(gridcell_type) , pointer :: gptr ! pointer to gridcell derived subtype
    type(landunit_type) , pointer :: lptr ! pointer to landunit derived subtype
    type(column_type) , pointer :: cptr   ! pointer to column derived subtype
    type(pft_type) , pointer :: pptr      ! pointer to pft derived subtype

    ! Set pointers into derived types for this module

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Get total global number of grid cells, landunits, columns and pfts

    call get_proc_global(numg,numl,numc,nump)
    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! For each land gridcell on global grid determine landunit,
    ! column and pft properties

    li    = begl-1
    ci    = begc-1
    pi    = begp-1

    if ( crop_prog )then
      ltype = istcrop
    else
      ltype = istsoil
    end if

    !----- Set clm3 variables -----
    do gdc = begg , endg

      nwtxy = gdc

      my_gcell = .false.
      if (gdc >= begg .and. gdc <= endg) then
        my_gcell = .true.
      end if

      ! Determine naturally vegetated landunit

      call set_landunit_veg_compete(               &
            ltype=istsoil, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

      ! Determine crop landunit

      call set_landunit_crop_noncompete(           &
            ltype=ltype, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

      ! Determine urban tall building district landunit

      call set_landunit_urban( &
!           ltype=isturb, wtxy=wtxy, vegxy=vegxy,   &
            ltype=isturb, udenstype=udens_tbd, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

      ! Determine urban high density landunit

      call set_landunit_urban( &
!           ltype=isturb, wtxy=wtxy, vegxy=vegxy,   &
            ltype=isturb, udenstype=udens_hd, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

      ! Determine urban medium density landunit

      call set_landunit_urban( &
!           ltype=isturb, wtxy=wtxy, vegxy=vegxy,   &
            ltype=isturb, udenstype=udens_md, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

      ! Determine lake, wetland and glacier landunits

      call set_landunit_wet_ice_lake(              &
            ltype=istdlak, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

      call set_landunit_wet_ice_lake(              &
            ltype=istwet, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

      call set_landunit_wet_ice_lake(              &
            ltype=istice, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

      ! Make ice sheet masks

      gptr%gris_mask(gdc) = 0.D0
      gptr%gris_area(gdc) = 0.D0
      gptr%aais_mask(gdc) = 0.D0
      gptr%aais_area(gdc) = 0.D0

      ! Greenland mask
      if ( (ldomain%latc(gdc) >  58. .and. ldomain%latc(gdc) <= 67.  .and.   &
            ldomain%lonc(gdc) > 302. .and. ldomain%lonc(gdc) < 330.) .or.    &
           (ldomain%latc(gdc) >  67. .and. ldomain%latc(gdc) <= 70. .and.    &
            ldomain%lonc(gdc) > 300. .and. ldomain%lonc(gdc) < 345.)  .or.   &
           (ldomain%latc(gdc) >  70. .and. ldomain%latc(gdc) <= 75. .and.    &
            ldomain%lonc(gdc) > 295. .and. ldomain%lonc(gdc) < 350.) .or.    &
           (ldomain%latc(gdc) >  75. .and. ldomain%latc(gdc) <= 79. .and.    &
            ldomain%lonc(gdc) > 285. .and. ldomain%lonc(gdc) < 350.) .or.    &
           (ldomain%latc(gdc) >  79. .and. ldomain%latc(gdc) <  85. .and.    &
            ldomain%lonc(gdc) > 290. .and. ldomain%lonc(gdc) < 355.) ) then
        gptr%gris_mask(gdc) = 1.0D0
      else if (ldomain%latc(gdc) < -60.) then
        gptr%aais_mask(gdc) = 1.0D0
      end if  ! Greenland or Antarctic grid cell

      ! Set clm3 lats/lons

      if (my_gcell) then
        gptr%latdeg(gdc) = ldomain%latc(gdc)
        gptr%londeg(gdc) = ldomain%lonc(gdc)
        gptr%lat(gdc)    = gptr%latdeg(gdc) * rpi/180.D0
        gptr%lon(gdc)    = gptr%londeg(gdc) * rpi/180.D0
        gptr%area(gdc)   = ldomain%area(gdc)
      end if

    end do

    ! Fill in subgrid datatypes

    call clm_ptrs_compdown()
    call clm_ptrs_check()

  end subroutine initGridcells
  !
  ! Assumes the part of the subgrid pointing up has been set.  Fills
  ! in the data pointing down.  Up is p_c, p_l, p_g, c_l, c_g, and l_g.
  !
  ! This algorithm assumes all indices are monotonically increasing.
  !
  ! Algorithm works as follows.  The p, c, and l loops march through
  ! the full arrays (nump, numc, and numl) checking the "up" indexes.
  ! As soon as the "up" index of the current (p,c,l) cell changes relative
  ! to the previous (p,c,l) cell, the *i array will be set to point down
  ! to that cell.  The *f array follows the same logic, so it's always the
  ! last "up" index from the previous cell when an "up" index changes.
  !
  ! For example, a case where p_c(1:4) = 1 and p_c(5:12) = 2.  This
  ! subroutine will set c_pi(1) = 1, c_pf(1) = 4, c_pi(2) = 5, c_pf(2) = 12.
  !
  subroutine clm_ptrs_compdown()
    use mod_clm_type, only : clm3 , gridcell_type , landunit_type , &
                        column_type , pft_type
    use mod_clm_decomp , only : get_proc_bounds
    implicit none
    ! beg/end glcp
    integer(ik4) :: begg , endg , begl , endl , begc , endc , begp , endp
    integer(ik4) :: l , c , p  ! loop counters
    ! tracks g,l,c,p indexes in arrays
    integer(ik4) :: curg , curl , curc
    ! pointer to gridcell derived subtype
    type(gridcell_type) , pointer :: gptr
    ! pointer to landunit derived subtype
    type(landunit_type) , pointer :: lptr
    ! pointer to column derived subtype
    type(column_type) , pointer :: cptr
    ! pointer to pft derived subtype
    type(pft_type) , pointer :: pptr

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    !--- Set the current c,l,g (curc, curl, curg) to zero for initialization,
    !---   these indices track the current "up" index.
    !--- Take advantage of locality of g/l/c/p cells
    !--- Loop p through full local begp:endp length
    !--- Separately check the p_c, p_l, and p_g indexes for a change in
    !---   the "up" index.
    !--- If there is a change, verify that the current c,l,g is within the
    !---   valid range, and set c_pi, l_pi, or g_pi to that current c,l,g
    !--- Constantly update the c_pf, l_pf, and g_pf array.  When the
    !---   g, l, c index changes, the *_pf array will be set correctly
    !--- Do the same for cols setting c_li, c_gi, c_lf, c_gf and
    !---   lunits setting l_gi, l_gf.

    curc = 0
    curl = 0
    curg = 0
    do p = begp , endp
      if (pptr%column(p) /= curc) then
        curc = pptr%column(p)
        if (curc < begc .or. curc > endc) then
          write(stderr,*) 'clm_ptrs_compdown ERROR: pcolumn ',p,curc,begc,endc
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
        cptr%pfti(curc) = p
      end if
      cptr%pftf(curc) = p
      cptr%npfts(curc) = cptr%pftf(curc) - cptr%pfti(curc) + 1
      if (pptr%landunit(p) /= curl) then
        curl = pptr%landunit(p)
        if (curl < begl .or. curl > endl) then
          write(stderr,*) 'clm_ptrs_compdown ERROR: plandunit ',p,curl,begl,endl
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
        lptr%pfti(curl) = p
      end if
      lptr%pftf(curl) = p
      lptr%npfts(curl) = lptr%pftf(curl) - lptr%pfti(curl) + 1
      if (pptr%gridcell(p) /= curg) then
        curg = pptr%gridcell(p)
        if (curg < begg .or. curg > endg) then
          write(stderr,*) 'clm_ptrs_compdown ERROR: pgridcell ',p,curg,begg,endg
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
        gptr%pfti(curg) = p
      end if
      gptr%pftf(curg) = p
      gptr%npfts(curg) = gptr%pftf(curg) - gptr%pfti(curg) + 1
    end do

    curg = 0
    curl = 0
    do c = begc , endc
      if (cptr%landunit(c) /= curl) then
        curl = cptr%landunit(c)
        if (curl < begl .or. curl > endl) then
          write(stderr,*) 'clm_ptrs_compdown ERROR: clandunit ',c,curl,begl,endl
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
        lptr%coli(curl) = c
      end if
      lptr%colf(curl) = c
      lptr%ncolumns(curl) = lptr%colf(curl) - lptr%coli(curl) + 1
      if (cptr%gridcell(c) /= curg) then
        curg = cptr%gridcell(c)
        if (curg < begg .or. curg > endg) then
          write(stderr,*) 'clm_ptrs_compdown ERROR: cgridcell ',c,curg,begg,endg
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
        gptr%coli(curg) = c
      end if
      gptr%colf(curg) = c
      gptr%ncolumns(curg) = gptr%colf(curg) - gptr%coli(curg) + 1
    end do

    curg = 0
    do l = begl , endl
      if (lptr%gridcell(l) /= curg) then
        curg = lptr%gridcell(l)
        if (curg < begg .or. curg > endg) then
          write(stderr,*) 'clm_ptrs_compdown ERROR: lgridcell ',l,curg,begg,endg
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
        gptr%luni(curg) = l
      end if
      gptr%lunf(curg) = l
      gptr%nlandunits(curg) = gptr%lunf(curg) - gptr%luni(curg) + 1
    end do
  end subroutine clm_ptrs_compdown
  !
  ! Checks and writes out a summary of subgrid data
  !
  subroutine clm_ptrs_check()
    use mod_clm_type , only : clm3 , gridcell_type , landunit_type , &
                        column_type , pft_type
    use mod_clm_decomp , only : get_proc_bounds
    implicit none
    type(gridcell_type) , pointer :: gptr ! pointer to gridcell derived subtype
    type(landunit_type) , pointer :: lptr ! pointer to landunit derived subtype
    type(column_type) , pointer :: cptr   ! pointer to column derived subtype
    type(pft_type) , pointer :: pptr      ! pointer to pft derived subtype
    ! beg/end indices
    integer(ik4) :: begg , endg , begl , endl , begc , endc , begp , endp
    integer(ik4) :: g , l , c , p ! loop counters
    logical :: error ! error flag

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    if (myid==italk) write(stdout,*) ' '
    if (myid==italk) write(stdout,*) '---clm_ptrs_check:'
    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    !--- check index ranges ---
    error = .false.
    if (minval(gptr%luni) < begl .or. maxval(gptr%luni) > endl) error=.true.
    if (minval(gptr%lunf) < begl .or. maxval(gptr%lunf) > endl) error=.true.
    if (minval(gptr%coli) < begc .or. maxval(gptr%coli) > endc) error=.true.
    if (minval(gptr%colf) < begc .or. maxval(gptr%colf) > endc) error=.true.
    if (minval(gptr%pfti) < begp .or. maxval(gptr%pfti) > endp) error=.true.
    if (minval(gptr%pftf) < begp .or. maxval(gptr%pftf) > endp) error=.true.
    if (error) then
      write(stderr,*) '   clm_ptrs_check: g index ranges - ERROR'
      write(stderr,*)'minval,beg,maxval,end'
      write(stderr,*) minval(gptr%luni),begl,maxval(gptr%luni),endl
      write(stderr,*) minval(gptr%lunf),begl,maxval(gptr%lunf),endl
      write(stderr,*) minval(gptr%coli),begc,maxval(gptr%coli),endc
      write(stderr,*) minval(gptr%colf),begc,maxval(gptr%colf),endc
      write(stderr,*) minval(gptr%pfti),begp,maxval(gptr%pfti),endp
      write(stderr,*) minval(gptr%pftf),begp,maxval(gptr%pftf),endp
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if
    if (myid==italk) write(stdout,*) '   clm_ptrs_check: g index ranges - OK'

    error = .false.
    if (minval(lptr%gridcell) < begg .or. &
        maxval(lptr%gridcell) > endg) error=.true.
    if (minval(lptr%coli) < begc .or. maxval(lptr%coli) > endc) error=.true.
    if (minval(lptr%colf) < begc .or. maxval(lptr%colf) > endc) error=.true.
    if (minval(lptr%pfti) < begp .or. maxval(lptr%pfti) > endp) error=.true.
    if (minval(lptr%pftf) < begp .or. maxval(lptr%pftf) > endp) error=.true.
    if (error) then
      write(stderr,*) '   clm_ptrs_check: l index ranges - ERROR'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if
    if (myid==italk) write(stdout,*) '   clm_ptrs_check: l index ranges - OK'

    error = .false.
    if (minval(cptr%gridcell) < begg .or. &
        maxval(cptr%gridcell) > endg) error=.true.
    if (minval(cptr%landunit) < begl .or. &
        maxval(cptr%landunit) > endl) error=.true.
    if (minval(cptr%pfti) < begp .or. maxval(cptr%pfti) > endp) error=.true.
    if (minval(cptr%pftf) < begp .or. maxval(cptr%pftf) > endp) error=.true.
    if (error) then
      write(stderr,*) '   clm_ptrs_check: c index ranges - ERROR'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if
    if (myid==italk) write(stdout,*) '   clm_ptrs_check: c index ranges - OK'

    error = .false.
    if (minval(pptr%gridcell) < begg .or. &
        maxval(pptr%gridcell) > endg) error=.true.
    if (minval(pptr%landunit) < begl .or. &
        maxval(pptr%landunit) > endl) error=.true.
    if (minval(pptr%column) < begc .or. &
        maxval(pptr%column) > endc) error=.true.
    if (error) then
      write(stderr,*) '   clm_ptrs_check: p index ranges - ERROR'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if
    if (myid==italk) write(stdout,*) '   clm_ptrs_check: p index ranges - OK'

    !--- check that indices in arrays are monotonically increasing ---
    error = .false.
    do g = begg+1 , endg
      if (gptr%luni(g) < gptr%luni(g-1)) error = .true.
      if (gptr%lunf(g) < gptr%lunf(g-1)) error = .true.
      if (gptr%coli(g) < gptr%coli(g-1)) error = .true.
      if (gptr%colf(g) < gptr%colf(g-1)) error = .true.
      if (gptr%pfti(g) < gptr%pfti(g-1)) error = .true.
      if (gptr%pftf(g) < gptr%pftf(g-1)) error = .true.
      if (error) then
        write(stderr,*) '   clm_ptrs_check: g mono increasing - ERROR'
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if
    end do
    if (myid==italk) write(stdout,*) '   clm_ptrs_check: g mono increasing - OK'

    error = .false.
    do l=begl+1,endl
      if (lptr%gridcell(l) < lptr%gridcell(l-1)) error = .true.
      if (lptr%coli(l) < lptr%coli(l-1)) error = .true.
      if (lptr%colf(l) < lptr%colf(l-1)) error = .true.
      if (lptr%pfti(l) < lptr%pfti(l-1)) error = .true.
      if (lptr%pftf(l) < lptr%pftf(l-1)) error = .true.
      if (error) then
        write(stderr,*) '   clm_ptrs_check: l mono increasing - ERROR'
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if
    end do
    if (myid==italk) write(stdout,*) '   clm_ptrs_check: l mono increasing - OK'

    error = .false.
    do c = begc+1 , endc
      if (cptr%gridcell(c) < cptr%gridcell(c-1)) error = .true.
      if (cptr%landunit(c) < cptr%landunit(c-1)) error = .true.
      if (cptr%pfti(c) < cptr%pfti(c-1)) error = .true.
      if (cptr%pftf(c) < cptr%pftf(c-1)) error = .true.
      if (error) then
        write(stderr,*) '   clm_ptrs_check: c mono increasing - ERROR'
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if
    end do
    if (myid==italk) write(stdout,*) '   clm_ptrs_check: c mono increasing - OK'

    error = .false.
    do p = begp+1 , endp
      if (pptr%gridcell(p) < pptr%gridcell(p-1)) error = .true.
      if (pptr%landunit(p) < pptr%landunit(p-1)) error = .true.
      if (pptr%column  (p) < pptr%column  (p-1)) error = .true.
      if (error) then
        write(stderr,*) '   clm_ptrs_check: p mono increasing - ERROR'
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if
    end do
    if (myid==italk) write(stdout,*) '   clm_ptrs_check: p mono increasing - OK'

    !--- check that the tree is internally consistent ---
    error = .false.
    do g = begg , endg
      do l = gptr%luni(g) , gptr%lunf(g)
        if (lptr%gridcell(l) /= g) error = .true.
        do c = lptr%coli(l) , lptr%colf(l)
          if (cptr%gridcell(c) /= g) error = .true.
          if (cptr%landunit(c) /= l) error = .true.
          do p = cptr%pfti(c) , cptr%pftf(c)
            if (pptr%gridcell(p) /= g) error = .true.
            if (pptr%landunit(p) /= l) error = .true.
            if (pptr%column(p)   /= c) error = .true.
            if (error) then
              write(stderr,*) '   clm_ptrs_check: tree consistent - ERROR'
              call fatal(__FILE__,__LINE__,'clm now stopping')
            end if
          end do
        end do
      end do
    end do
    if (myid==italk) write(stdout,*) '   clm_ptrs_check: tree consistent - OK'
    if (myid==italk) write(stdout,*) ' '
  end subroutine clm_ptrs_check
  !
  ! Initialize vegetated landunit with competition
  !
  subroutine set_landunit_veg_compete(ltype,nw,gi,li,ci,pi,setdata)
    use mod_clm_type , only : clm3 , model_type , gridcell_type , &
            landunit_type , column_type , pft_type
    use mod_clm_subgrid , only : subgrid_get_gcellinfo
    use mod_clm_varpar , only : numpft , maxpatch_pft , numcft
    use mod_clm_varctl , only : allocate_all_vegpfts , create_crop_landunit
    implicit none
    integer(ik4) , intent(in) :: ltype        ! landunit type
    integer(ik4) , intent(in) :: nw           ! cell index
    integer(ik4) , intent(in) :: gi           ! gridcell index
    integer(ik4) , intent(inout) :: li        ! landunit index
    integer(ik4) , intent(inout) :: ci        ! column index
    integer(ik4) , intent(inout) :: pi        ! pft index
    logical , intent(in) :: setdata           ! set info or just compute
    integer(ik4) :: m            ! m index in wtxy(nw,m)
    integer(ik4) :: n            ! loop index
    integer(ik4) :: npfts        ! number of pfts in landunit
    integer(ik4) :: ncols        ! number of columns in landu
    integer(ik4) :: pitype       ! pft itype
    real(rk8) :: wtlunit2gcell   ! landunit weight in gridcell
    type(landunit_type) , pointer :: lptr  ! pointer to landunit
    type(column_type) , pointer :: cptr    ! pointer to column
    type(pft_type) , pointer :: pptr       ! pointer to pft

    ! Set decomposition properties

    call subgrid_get_gcellinfo(nw, nveg=npfts, wtveg=wtlunit2gcell)

    if (npfts > 0) then

      ! Set pointers into derived types for this module

      lptr => clm3%g%l
      cptr => clm3%g%l%c
      pptr => clm3%g%l%c%p

      ncols = 1

      li = li + 1
      ci = ci + 1

      if (setdata) then
        ! Set landunit properties
        lptr%ifspecial(li) = .false.
        lptr%lakpoi(li)    = .false.
        lptr%urbpoi(li)    = .false.
        lptr%itype(li)     = ltype

        lptr%gridcell (li) = gi
        lptr%wtgcell(li) = wtlunit2gcell

        ! Set column properties for this landunit (only one column on landunit)
        cptr%itype(ci)    = 1

        cptr%gridcell (ci) = gi
        cptr%wtgcell(ci) = wtlunit2gcell
        cptr%landunit (ci) = li
        cptr%wtlunit(ci) = 1.0D0
      end if ! setdata

      ! Set pft properties for this landunit

      if (create_crop_landunit) then
        do n = 1 , numpft+1-numcft
          pi = pi + 1
          pitype = n-1
          if (setdata) then
            pptr%mxy(pi)      = n
            pptr%itype(pi)    = pitype
            pptr%gridcell(pi) = gi
            pptr%landunit(pi) = li
            pptr%column (pi) = ci

            if (wtlunit2gcell > 0.D0) then
              pptr%wtgcell(pi) = 0.0D0
              pptr%wtlunit(pi) = 0.0D0
              pptr%wtcol(pi) = 0.0D0
              do m = 1 , maxpatch_pft
                if (vegxy(nw,m) == pitype) then
                  pptr%wtgcell(pi) = pptr%wtgcell(pi) + wtxy(nw,m)
                  pptr%wtlunit(pi) = pptr%wtlunit(pi) + &
                          wtxy(nw,m) / wtlunit2gcell
                  pptr%wtcol(pi) = pptr%wtcol(pi) + wtxy(nw,m) / wtlunit2gcell
                end if
              end do
            else  ! wtlunit2gcell == 0.D0
              ! TODO WJS: Temporarily setting this to equal weighting for all
              ! pfts. In the future, we could potentially get some info about
              ! this from the surface dataset, if it is changed to give
              ! pct_pft as % of the pft on the landunit
              pptr%wtgcell(pi) = 0.D0
              pptr%wtlunit(pi) = 1.D0 / (numpft+1-numcft)
              pptr%wtcol(pi)   = 1.D0 / (numpft+1-numcft)
            end if
          end if ! setdata
        end do
      else if (allocate_all_vegpfts) then
        do n = 1 , numpft+1
          pi = pi + 1
          pitype = n-1
          if (setdata) then
            pptr%mxy(pi)      = n
            pptr%itype(pi)    = pitype
            pptr%gridcell(pi) = gi
            pptr%landunit(pi) = li
            pptr%column (pi) = ci

            if (wtlunit2gcell > 0.D0) then
              pptr%wtgcell(pi) = 0.0D0
              pptr%wtlunit(pi) = 0.0D0
              pptr%wtcol(pi) = 0.0D0
              do m = 1 , maxpatch_pft
                if (vegxy(nw,m) == pitype) then
                  pptr%wtgcell(pi) = pptr%wtgcell(pi) + wtxy(nw,m)
                  pptr%wtlunit(pi) = pptr%wtlunit(pi) + &
                          wtxy(nw,m) / wtlunit2gcell
                  pptr%wtcol(pi) = pptr%wtcol(pi) + wtxy(nw,m) / wtlunit2gcell
                end if
              end do
            else  ! wtlunit2gcell == 0.D0
              ! TODO WJS: Temporarily setting this to equal weighting for all
              ! pfts. In the future, we could potentially get some info about
              ! this from the surface dataset, if it is changed to give
              ! pct_pft as % of the pft on the landunit
              pptr%wtgcell(pi) = 0.D0
              pptr%wtlunit(pi) = 1.D0 / (numpft+1)
              pptr%wtcol(pi)   = 1.D0 / (numpft+1)
            end if
          end if ! setdata
        end do
      else
        write(stderr,*) 'allocate_all_vegpfts=false is no longer supported'
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if
    end if
  end subroutine set_landunit_veg_compete
  !
  ! Initialize wet_ice_lake landunits that are non-urban
  ! (lake, wetland, glacier)
  !
  subroutine set_landunit_wet_ice_lake (ltype,nw,gi,li,ci,pi,setdata)
    use mod_clm_type , only : clm3 , model_type , gridcell_type , &
            landunit_type , column_type,pft_type
    use mod_clm_subgrid , only : subgrid_get_gcellinfo
    use mod_clm_varcon , only : istice , istwet , istdlak
    use mod_clm_varpar , only : npatch_lake , npatch_glacier , npatch_wet
    implicit none
    integer(ik4) , intent(in) :: ltype ! landunit type
    integer(ik4) , intent(in) :: nw    ! cell index
    integer(ik4) , intent(in) :: gi    ! gridcell index
    integer(ik4) , intent(inout) :: li ! landunit index
    integer(ik4) , intent(inout) :: ci ! column index
    integer(ik4) , intent(inout) :: pi ! pft index
    logical , intent(in) :: setdata    ! set info or just compute
    integer(ik4) :: m     ! m index in wtxy(nw,m)
    integer(ik4) :: ctype ! column type
    integer(ik4) :: npfts ! number of pfts in landunit
    integer(ik4) :: ncols ! number of columns in landu
    real(rk8) :: wtlunit2gcell  ! landunit weight in gridcell
    real(rk8) :: wtcol2lunit    ! col weight in landunit
    type(landunit_type) , pointer :: lptr ! pointer to landunit
    type(column_type) , pointer :: cptr   ! pointer to column
    type(pft_type) , pointer :: pptr      ! pointer to pft

    ! Set decomposition properties

    if (ltype == istwet) then
      call subgrid_get_gcellinfo(nw, nwetland=npfts, wtwetland=wtlunit2gcell)
      m = npatch_wet
    else if (ltype == istdlak) then
      call subgrid_get_gcellinfo(nw, nlake=npfts, wtlake=wtlunit2gcell)
      m = npatch_lake
    else if (ltype == istice) then
      call subgrid_get_gcellinfo(nw, nglacier=npfts, wtglacier=wtlunit2gcell)
      m = npatch_glacier
    else
      write(stderr,*)' set_landunit_wet_ice_lake: ltype of ',ltype,' not valid'
      write(stderr,*) &
              ' only istwet, istdlak, istice and ltypes are valid'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    if (npfts > 0) then
      ! Set pointers into derived types for this module
      lptr => clm3%g%l
      cptr => clm3%g%l%c
      pptr => clm3%g%l%c%p

      if ( npfts /=1 ) then
        write(stderr,*)' set_landunit_wet_ice_lake: compete landunit must'// &
               ' have one column and one pft '
        write(stderr,*)' current values of ncols, pfts=',ncols,npfts
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if

      ncols = 1

      ! Currently assume that each landunit only has only one column
      ! (of type 1) and that each column has its own pft

      wtcol2lunit = 1.0D0/ncols
      ctype = 1

      li = li + 1
      ci = ci + 1
      pi = pi + 1

      if (setdata) then

        ! Determine landunit properties

        lptr%itype    (li) = ltype
        lptr%ifspecial(li) = .true.
        lptr%urbpoi   (li) = .false.
        if (ltype == istdlak) then
          lptr%lakpoi(li) = .true.
        else
          lptr%lakpoi(li) = .false.
        end if

        lptr%gridcell (li) = gi
        lptr%wtgcell(li) = wtlunit2gcell

        ! Determine column and properties
        ! For the wet, ice or lake landunits it is assumed that each
        ! column has its own pft

        cptr%itype(ci)    = ctype

        cptr%gridcell (ci) = gi
        cptr%wtgcell(ci) = wtcol2lunit * wtlunit2gcell
        cptr%landunit (ci) = li
        cptr%wtlunit(ci) = wtcol2lunit

        ! Set pft properties

        pptr%mxy(pi)      = m
        pptr%itype(pi)    = vegxy(nw,m)

        pptr%gridcell (pi) = gi
        pptr%wtgcell(pi) = wtcol2lunit * wtlunit2gcell
        pptr%landunit (pi) = li
        pptr%wtlunit(pi) = wtcol2lunit
        pptr%column (pi) = ci
        pptr%wtcol(pi) = 1.0D0
      end if ! setdata
    end if     ! npfts > 0
  end subroutine set_landunit_wet_ice_lake
  !
  ! Initialize crop landunit without competition
  !
  subroutine set_landunit_crop_noncompete (ltype,nw,gi,li,ci,pi,setdata)
    use mod_clm_type , only : clm3 , model_type , gridcell_type , &
            landunit_type , column_type,pft_type
    use mod_clm_subgrid , only : subgrid_get_gcellinfo
    use mod_clm_varctl , only : create_crop_landunit
    use mod_clm_varpar , only : maxpatch_pft , numcft
    implicit none
    integer(ik4) , intent(in) :: ltype       ! landunit type
    integer(ik4) , intent(in) :: nw          ! cell index
    integer(ik4) , intent(in) :: gi          ! gridcell index
    integer(ik4) , intent(inout) :: li       ! landunit index
    integer(ik4) , intent(inout) :: ci       ! column index
    integer(ik4) , intent(inout) :: pi       ! pft index
    logical , intent(in) :: setdata          ! set info or just compute
    integer(ik4) :: m       ! m index in wtxy(nw,m)
    integer(ik4) :: npfts   ! number of pfts in landunit
    integer(ik4) :: ncols   ! number of columns in landu
    real(rk8) :: wtlunit2gcell ! landunit weight in gridcell
    type(landunit_type) , pointer :: lptr  ! pointer to landunit
    type(column_type) , pointer :: cptr    ! pointer to column
    type(pft_type) , pointer :: pptr       ! pointer to pft

    ! Set decomposition properties

    call subgrid_get_gcellinfo(nw, ncrop=npfts, wtcrop=wtlunit2gcell)

    if (npfts > 0) then

      ! Set pointers into derived types for this module

      lptr => clm3%g%l
      cptr => clm3%g%l%c
      pptr => clm3%g%l%c%p

      ! Set landunit properties - each column has its own pft

      ncols = npfts

      li = li + 1

      if (setdata) then
        lptr%itype(li)     = ltype
        lptr%ifspecial(li) = .false.
        lptr%lakpoi(li)    = .false.
        lptr%urbpoi(li)    = .false.
        lptr%gridcell (li) = gi
        lptr%wtgcell(li) = wtlunit2gcell
      end if ! setdata

      ! Set column and pft properties for this landunit
      ! (each column has its own pft)

      if (create_crop_landunit) then
        do m = maxpatch_pft-numcft+1 , maxpatch_pft
          ci = ci + 1
          pi = pi + 1

          if (setdata) then
            cptr%itype(ci)    = 1
            pptr%itype(pi)    = m - 1
            pptr%mxy(pi)      = m

            cptr%gridcell (ci) = gi
            cptr%wtgcell(ci) = wtxy(nw,m)
            cptr%landunit (ci) = li

            pptr%gridcell (pi) = gi
            pptr%wtgcell(pi) = wtxy(nw,m)
            pptr%landunit (pi) = li
            pptr%column (pi) = ci
            pptr%wtcol(pi) = 1.D0
            if (wtlunit2gcell > 0) then
              cptr%wtlunit(ci) = wtxy(nw,m) / wtlunit2gcell
              pptr%wtlunit(pi) = wtxy(nw,m) / wtlunit2gcell
            else
              ! TODO WJS: Temporarily setting this to equal weighting for
              ! all crop pfts. In the future, we could potentially get some
              ! info about this from the surface dataset, if it is changed
              ! to give pct_cft as % of the cft on the landunit
              cptr%wtlunit(ci) = 1.D0 / numcft
              pptr%wtlunit(pi) = 1.D0 / numcft
            end if
          end if ! setdata
        end do
      end if
    end if
  end subroutine set_landunit_crop_noncompete
  !
  ! Initialize urban landunits
  !
  subroutine set_landunit_urban (ltype,udenstype,nw,gi,li,ci,pi,setdata)
    use mod_clm_varcon , only : isturb , icol_roof , icol_sunwall , &
            icol_shadewall , icol_road_perv , icol_road_imperv ,    &
            udens_tbd , udens_hd , udens_md , udens_base
    use mod_clm_varpar , only : npatch_urban_tbd , npatch_urban_hd , &
            npatch_urban_md , maxpatch_urb
    use mod_clm_type , only : clm3 , model_type , gridcell_type , &
            landunit_type , column_type , pft_type
    use mod_clm_subgrid , only : subgrid_get_gcellinfo
    use mod_clm_urbaninput , only : urbinp
    implicit none
    integer(ik4) , intent(in) :: ltype        ! landunit type
    integer(ik4) , intent(in) :: udenstype    ! urban density type
    integer(ik4) , intent(in) :: nw           ! cell index
    integer(ik4) , intent(in) :: gi           ! gridcell index
    integer(ik4) , intent(inout) :: li        ! landunit index
    integer(ik4) , intent(inout) :: ci        ! column index
    integer(ik4) , intent(inout) :: pi        ! pft index
    logical , intent(in) :: setdata           ! set info or just compute
    integer(ik4) :: m          ! m index in wtxy(nw,m)
    integer(ik4) :: n          ! urban density type index
    integer(ik4) :: ctype      ! column type
    integer(ik4) :: npfts      ! number of pfts in landunit
    integer(ik4) :: ncols      ! number of columns in landunit
    integer(ik4) :: npatch     ! npatch for the given urban density class
    real(rk8) :: wtlunit2gcell ! weight relative to gridcell of landunit
    real(rk8) :: wtcol2lunit   ! weight of column with respect to landunit
    real(rk8) :: wtlunit_roof  ! weight of roof with respect to landunit
    ! weight of pervious road column with respect to total road
    real(rk8) :: wtroad_perv
    type(landunit_type) , pointer :: lptr ! pointer to landunit derived subtype
    type(column_type) , pointer :: cptr   ! pointer to column derived subtype
    type(pft_type) , pointer :: pptr      ! pointer to pft derived subtype

    ! Set decomposition properties, and set variables specific to urban
    ! density type

    select case (udenstype)
      case (udens_tbd)
        call subgrid_get_gcellinfo(nw, nurban_tbd=npfts, &
                wturban_tbd=wtlunit2gcell)
        npatch = npatch_urban_tbd
      case (udens_hd)
        call subgrid_get_gcellinfo(nw, nurban_hd=npfts, &
                 wturban_hd=wtlunit2gcell)
        npatch = npatch_urban_hd
      case (udens_md)
        call subgrid_get_gcellinfo(nw, nurban_md=npfts, &
                wturban_md=wtlunit2gcell)
        npatch = npatch_urban_md
      case default
        write(stderr,*)' set_landunit_urban: unknown udenstype: ', udenstype
        call fatal(__FILE__,__LINE__,'clm now stopping')
    end select

    n = udenstype - udens_base

    if (npfts > 0) then

      ! Set pointers into derived types for this module

      lptr => clm3%g%l
      cptr => clm3%g%l%c
      pptr => clm3%g%l%c%p

      ! Determine landunit properties - each columns has its own pft

      ncols = npfts

      li = li + 1
      if (setdata) then
        lptr%itype    (li) = ltype
        lptr%udenstype(li) = udenstype
        lptr%ifspecial(li) = .true.
        lptr%lakpoi   (li) = .false.
        lptr%urbpoi   (li) = .true.

        lptr%gridcell (li) = gi
        lptr%wtgcell  (li) = wtlunit2gcell
      end if

      ! Loop through columns for this landunit and set the column and pft
      ! properties For the urban landunits it is assumed that each column
      ! has its own pft

      do m = npatch , npatch + maxpatch_urb - 1
        if (wtxy(nw,m) > 0.D0) then

          wtlunit_roof = urbinp%wtlunit_roof(nw,n)
          wtroad_perv  = urbinp%wtroad_perv(nw,n)

          if (m == npatch  ) then
            ctype = icol_roof
            wtcol2lunit = wtlunit_roof
          else if (m == npatch+1) then
            ctype = icol_sunwall
            wtcol2lunit = (1. - wtlunit_roof)/3
          else if (m == npatch+2) then
            ctype = icol_shadewall
            wtcol2lunit = (1. - wtlunit_roof)/3
          else if (m == npatch+3) then
            ctype = icol_road_imperv
            wtcol2lunit = ((1. - wtlunit_roof)/3) * (1.-wtroad_perv)
          else if (m == npatch+4) then
            ctype = icol_road_perv
            wtcol2lunit = ((1. - wtlunit_roof)/3) * (wtroad_perv)
          end if

          ci = ci + 1
          pi = pi + 1

          if (setdata) then
            cptr%itype(ci)     = ctype

            cptr%gridcell (ci) = gi
            cptr%wtgcell  (ci) = wtcol2lunit * wtlunit2gcell
            cptr%landunit (ci) = li
            cptr%wtlunit  (ci) = wtcol2lunit

            pptr%mxy     (pi)  = m
            pptr%itype   (pi)  = vegxy(nw,m)

            pptr%gridcell(pi)  = gi
            pptr%wtgcell (pi)  = wtcol2lunit * wtlunit2gcell
            pptr%landunit(pi)  = li
            pptr%wtlunit (pi)  = wtcol2lunit
            pptr%column  (pi)  = ci
            pptr%wtcol   (pi)  = 1.0D0
          end if
        end if
      end do   ! end of loop through urban columns-pfts
    end if
  end subroutine set_landunit_urban

end module mod_clm_initgridcells
