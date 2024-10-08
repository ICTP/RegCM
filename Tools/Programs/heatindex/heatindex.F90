program calc_heatindex
  use iso_fortran_env
  use mod_heatindex
  use netcdf
  implicit none
  integer , parameter :: max_path = 1024
  integer :: i , j , ii , jj
  integer :: ntime , nx , ny
  integer :: itime , itime_bnds , itasmax , ihurs , ihtindx , icrs
  integer :: tm_ncid , rh_ncid , hi_ncid
  integer :: narg , nfstat , nfmode , nftype
  integer :: ndims , nvars , nglobalatts , numatts , unlimdimid
  character (len=max_path) :: progname
  character (len=max_path) :: tasmax_filename , hurs_filename
  character (len=max_path) :: hindex_filename
  character (len=nf90_max_name) :: item_name
  integer , allocatable , dimension(:) :: dimids , numdims , dimlens
  integer , allocatable , dimension(:,:) :: vdimids
  integer , allocatable , dimension(:) :: varids
  integer , dimension(3) :: istart , icount
  integer , dimension(1) :: its , itc
  integer , dimension(2) :: itbs , itbc

  real, dimension(1) :: var0d
  real, dimension(:) , allocatable :: var1d
  real, dimension(:,:) , allocatable :: var2d
  real(kind=real64), dimension(:,:) , allocatable :: tasmax
  real(kind=real64), dimension(:,:) , allocatable :: hurs
  real, dimension(:,:) , allocatable :: heat_index
  real, dimension(1) :: time
  real, dimension(2) :: time_bnds

  itime_bnds = -1

  CALL get_command_argument(0, progname)
  narg = command_argument_count()
  if ( narg < 2 ) then
    print *, 'Not enough arguments: need tasmax filename and hurs filename.'
    print *, 'Example: ', trim(progname),' tasmax.nc hurs.nc'
    stop 'ARGUMENTS'
  end if
  CALL get_command_argument(1, tasmax_filename)
  CALL get_command_argument(2, hurs_filename)
  nfstat = nf90_open(tasmax_filename,nf90_nowrite,tm_ncid)
  if ( nfstat /= nf90_noerr ) then
    print *, 'Cannot open ',trim(tasmax_filename)
    print *, nf90_strerror(nfstat)
    stop 'INPUT ERROR'
  end if
  nfstat = nf90_open(hurs_filename,nf90_nowrite,rh_ncid)
  if ( nfstat /= nf90_noerr ) then
    print *, 'Cannot open ',trim(hurs_filename)
    print *, nf90_strerror(nfstat)
    stop 'INPUT ERROR'
  end if

  hindex_filename = replacestr(basename(tasmax_filename),'tasmax','hindex')
  nfmode = ior(nf90_clobber,nf90_netcdf4)
  nfstat = nf90_create(hindex_filename, nfmode, hi_ncid)
  if ( nfstat /= nf90_noerr ) then
    print *, 'Cannot open ',trim(hindex_filename)
    print *, nf90_strerror(nfstat)
    stop 'OUTPUT ERROR'
  end if

  nfstat = nf90_inquire(tm_ncid, ndims, nvars, nglobalatts, unlimdimid)
  if ( nfstat /= nf90_noerr ) then
    print *, 'Cannot read specs from ',trim(tasmax_filename)
    print *, nf90_strerror(nfstat)
    stop 'READ ERROR'
  end if

  nfstat = nf90_inq_varid(rh_ncid, "hurs", ihurs)
  if ( nfstat /= nf90_noerr ) then
    print *, 'Cannot read variable hurs from ',trim(hurs_filename)
    print *, nf90_strerror(nfstat)
    stop 'READ ERROR'
  end if

  allocate(dimids(ndims))
  allocate(varids(nvars))
  allocate(numdims(nvars))
  allocate(vdimids(nvars,ndims))
  allocate(dimlens(ndims))
  do i = 1 , ndims
    nfstat = nf90_inquire_dimension(tm_ncid, i, item_name, len=dimlens(i))
    if ( nfstat /= nf90_noerr ) then
      print *, 'Cannot read specs from ',trim(tasmax_filename)
      print *, nf90_strerror(nfstat)
      stop 'READ ERROR'
    end if
    if ( i == unlimdimid .or. item_name == 'time' ) then
      nfstat = nf90_def_dim(hi_ncid, item_name, nf90_unlimited, dimids(i))
      if ( nfstat /= nf90_noerr ) then
        print *, 'Cannot define dimension',trim(item_name)
        print *, nf90_strerror(nfstat)
        stop 'READ ERROR'
      end if
    else
      nfstat = nf90_def_dim(hi_ncid, item_name, dimlens(i), dimids(i))
      if ( nfstat /= nf90_noerr ) then
        print *, 'Cannot define dimension',trim(item_name)
        print *, nf90_strerror(nfstat)
        stop 'WRITE ERROR'
      end if
    end if
  end do

  do i = 1 , nglobalatts
    nfstat = nf90_inq_attname(tm_ncid, nf90_global, i, item_name)
    if ( nfstat /= nf90_noerr ) then
      print *, 'Cannot inquire attribute from ',trim(tasmax_filename)
      print *, nf90_strerror(nfstat)
      stop 'READ ERROR'
    end if
    nfstat = nf90_copy_att(tm_ncid, nf90_global, item_name, &
      hi_ncid, nf90_global)
    if ( nfstat /= nf90_noerr ) then
      print *, 'Cannot copy attribute ',trim(item_name)
      print *, nf90_strerror(nfstat)
      stop 'WRITE ERROR'
    end if
  end do

  do i = 1 , nvars
    nfstat = nf90_inquire_variable(tm_ncid, i, item_name, nftype, &
                         numdims(i), vdimids(i,:), numatts)
    if ( nfstat /= nf90_noerr ) then
      print *, 'Cannot inquire variable from ',trim(tasmax_filename)
      print *, nf90_strerror(nfstat)
      stop 'READ ERROR'
    end if
    if ( item_name == 'time' ) then
      itime = i
    else if ( item_name == 'time_bnds' ) then
      itime_bnds = i
    else if ( item_name == 'rotated_pole' .or. item_name == 'crs' ) then
      icrs = i
    else if ( item_name == 'tasmax' ) then
      itasmax = i
      ihtindx = i
      item_name = 'heat_index'
    end if
    nfstat = nf90_def_var(hi_ncid, item_name, nftype, &
               vdimids(i,1:numdims(i)), varids(i))
    if ( nfstat /= nf90_noerr ) then
      print *, 'Cannot define variable ',trim(item_name)
      print *, nf90_strerror(nfstat)
      stop 'DEFINE ERROR'
    end if
    if ( i /= itasmax ) then
      do j = 1 , numatts
        nfstat = nf90_inq_attname(tm_ncid, i, j, item_name)
        if ( nfstat /= nf90_noerr ) then
          print *, 'Cannot inquire attribute from ',trim(tasmax_filename)
          print *, nf90_strerror(nfstat)
          stop 'READ ERROR'
        end if
        nfstat = nf90_copy_att(tm_ncid, i, item_name, &
                    hi_ncid, i)
        if ( nfstat /= nf90_noerr ) then
          print *, 'Cannot copy attribute ',trim(item_name)
          print *, nf90_strerror(nfstat)
          stop 'WRITE ERROR'
        end if
      end do
    else
      nfstat = nf90_put_att(hi_ncid, i, 'standard_name', 'heat_index')
      if ( nfstat /= nf90_noerr ) then
        print *, 'Cannot define attribute standard_name'
        print *, nf90_strerror(nfstat)
        stop 'WRITE ERROR'
      end if
      nfstat = nf90_put_att(hi_ncid, i, 'long_name', &
        'Heat Index as defined in Lu and Romps, '// &
        'Extending the Heat Index, JAMC, 2022')
      if ( nfstat /= nf90_noerr ) then
        print *, 'Cannot define attribute long_name'
        print *, nf90_strerror(nfstat)
        stop 'WRITE ERROR'
      end if
      nfstat = nf90_put_att(hi_ncid, i, 'units', 'K')
      if ( nfstat /= nf90_noerr ) then
        print *, 'Cannot define attribute units'
        print *, nf90_strerror(nfstat)
        stop 'WRITE ERROR'
      end if
      nfstat = nf90_put_att(hi_ncid, i, 'coordinates', 'height lat lon')
      if ( nfstat /= nf90_noerr ) then
        print *, 'Cannot define attribute coordinates'
        print *, nf90_strerror(nfstat)
        stop 'WRITE ERROR'
      end if
      nfstat = nf90_put_att(hi_ncid, i, 'grid_mapping', 'rotated_pole')
      if ( nfstat /= nf90_noerr ) then
        print *, 'Cannot define attribute grid_mapping'
        print *, nf90_strerror(nfstat)
        stop 'WRITE ERROR'
      end if
    end if
  end do

  nfstat = nf90_enddef(hi_ncid)
  if ( nfstat /= nf90_noerr ) then
    print *, 'Cannot end define in ',trim(hindex_filename)
    print *, nf90_strerror(nfstat)
    stop 'OUTPUT ERROR'
  end if

  do i = 1 , nvars
    if ( i == itime .or. i == itasmax .or. i == itime_bnds ) cycle
    if ( numdims(i) == 0 ) then
      if ( i == icrs ) cycle
      nfstat = nf90_get_var(tm_ncid,i,var0d)
      if ( nfstat /= nf90_noerr ) then
        print *, 'Cannot read variable from ',trim(tasmax_filename)
        print *, nf90_strerror(nfstat)
        stop 'READ ERROR'
      end if
      nfstat = nf90_put_var(hi_ncid,i,var0d)
      if ( nfstat /= nf90_noerr ) then
        print *, 'Cannot write variable to ',trim(hindex_filename)
        print *, nf90_strerror(nfstat)
        stop 'WRITE ERROR'
      end if
    else if ( numdims(i) == 1 ) then
      allocate(var1d(dimlens(vdimids(i,1))))
      nfstat = nf90_get_var(tm_ncid,i,var1d)
      if ( nfstat /= nf90_noerr ) then
        print *, 'Cannot read variable from ',trim(tasmax_filename)
        print *, nf90_strerror(nfstat)
        stop 'READ ERROR'
      end if
      nfstat = nf90_put_var(hi_ncid,i,var1d)
      if ( nfstat /= nf90_noerr ) then
        print *, 'Cannot write variable to ',trim(hindex_filename)
        print *, nf90_strerror(nfstat)
        stop 'WRITE ERROR'
      end if
      deallocate(var1d)
    elseif ( numdims(i) == 2 ) then
      allocate(var2d(dimlens(vdimids(i,1)),dimlens(vdimids(i,2))))
      nfstat = nf90_get_var(tm_ncid,i,var2d)
      if ( nfstat /= nf90_noerr ) then
        print *, 'Cannot read variable from ',trim(tasmax_filename)
        print *, nf90_strerror(nfstat)
        stop 'READ ERROR'
      end if
      nfstat = nf90_put_var(hi_ncid,i,var2d)
      if ( nfstat /= nf90_noerr ) then
        print *, 'Cannot write variable to ',trim(hindex_filename)
        print *, nf90_strerror(nfstat)
        stop 'WRITE ERROR'
      end if
      deallocate(var2d)
    end if
  end do

  ntime = dimlens(vdimids(itime,1))
  nx = dimlens(vdimids(itasmax,1))
  ny = dimlens(vdimids(itasmax,2))
  allocate(tasmax(nx,ny))
  allocate(hurs(nx,ny))
  allocate(heat_index(nx,ny))
  istart(1) = 1
  icount(1) = nx
  istart(2) = 1
  icount(2) = ny
  icount(3) = 1
  itc(1) = 1
  itbs(1) = 1
  itbc(1) = 2
  itbc(2) = 1
  do i = 1 , ntime
    print *, 'Time : ', i
    istart(3) = i
    its(1) = i
    itbs(2) = i
    nfstat = nf90_get_var(tm_ncid,itasmax,tasmax,istart,icount)
    if ( nfstat /= nf90_noerr ) then
      print *, 'Cannot read variable from ',trim(tasmax_filename)
      print *, nf90_strerror(nfstat)
      stop 'READ ERROR'
    end if
    nfstat = nf90_get_var(rh_ncid,ihurs,hurs,istart,icount)
    if ( nfstat /= nf90_noerr ) then
      print *, 'Cannot read variable from ',trim(hurs_filename)
      print *, nf90_strerror(nfstat)
      stop 'READ ERROR'
    end if
    nfstat = nf90_get_var(tm_ncid,itime,time,its,itc)
    if ( nfstat /= nf90_noerr ) then
      print *, 'Cannot read variable from ',trim(tasmax_filename)
      print *, nf90_strerror(nfstat)
      stop 'READ ERROR'
    end if
    if ( itime_bnds > 0 ) then
      nfstat = nf90_get_var(tm_ncid,itime_bnds,time_bnds,itbs,itbc)
      if ( nfstat /= nf90_noerr ) then
        print *, 'Cannot read variable from ',trim(tasmax_filename)
        print *, nf90_strerror(nfstat)
        stop 'READ ERROR'
      end if
    end if

!$OMP PARALLEL PRIVATE(ii)
!$OMP DO
    do jj = 1 , ny
      do ii = 1 , nx
        heat_index(ii,jj) = heatindex(tasmax(ii,jj),hurs(ii,jj)*0.01)
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    nfstat = nf90_put_var(hi_ncid,itime,time,its,itc)
    if ( nfstat /= nf90_noerr ) then
      print *, 'Cannot write variable time to ',trim(hindex_filename)
      print *, nf90_strerror(nfstat)
      stop 'READ ERROR'
    end if
    if ( itime_bnds > 0 ) then
      nfstat = nf90_put_var(hi_ncid,itime_bnds,time_bnds,itbs,itbc)
      if ( nfstat /= nf90_noerr ) then
        print *, 'Cannot write variable time_bounds to ',trim(hindex_filename)
        print *, nf90_strerror(nfstat)
        stop 'READ ERROR'
      end if
    end if
    nfstat = nf90_put_var(hi_ncid,ihtindx,heat_index,istart,icount)
    if ( nfstat /= nf90_noerr ) then
      print *, 'Cannot write variable heat_index to ',trim(hindex_filename)
      print *, nf90_strerror(nfstat)
      stop 'WRITE ERROR'
    end if

    nfstat = nf90_sync(hi_ncid)
    if ( nfstat /= nf90_noerr ) then
      print *, 'Cannot sync to ',trim(hindex_filename)
      print *, nf90_strerror(nfstat)
      stop 'WRITE ERROR'
    end if

  end do

  nfstat = nf90_close(hi_ncid)
  if ( nfstat /= nf90_noerr ) then
    print *, 'Cannot close ',trim(hindex_filename)
    print *, nf90_strerror(nfstat)
    stop 'OUTPUT ERROR'
  end if

  nfstat = nf90_close(tm_ncid)
  nfstat = nf90_close(rh_ncid)

  deallocate(dimids)
  deallocate(varids)
  deallocate(numdims)
  deallocate(vdimids)
  deallocate(dimlens)
  deallocate(tasmax)
  deallocate(hurs)
  deallocate(heat_index)

  contains

  pure recursive function replacestr(string,search,sub) result(mstring)
    implicit none
    character(len=*), intent(in) :: string, search, sub
    character(len=:), allocatable :: mstring
    integer :: i , stringlen , searchlen
    stringlen = len(string)
    searchlen = len(search)
    if ( stringlen == 0 .or. searchlen == 0 ) then
      mstring = ""
      return
    else if ( stringlen < searchlen ) then
      mstring = string
      return
    end if
    i = 1
    do
      if ( string(i:i+searchlen-1) == search ) then
        mstring = string(1:i-1)//sub// &
          replacestr(string(i+searchlen:stringlen),search,sub)
        exit
      end if
      if ( i+searchlen > stringlen ) then
        mstring = string
        exit
      end if
      i = i + 1
      cycle
    end do
  end function replacestr

  elemental pure function lower(str,istart,istop) result (string)
    implicit none
    character(*), intent(in) :: str
    character(len(str)) :: string
    integer , intent(in) , optional :: istart , istop
    integer :: i
    integer :: ibegin , iend
    string = str
    ibegin = 1
    if ( present(istart) ) then
      ibegin = max(ibegin,istart)
    end if
    iend = len_trim(str)
    if ( present(istop) ) then
      iend = min(iend,istop)
    end if
    do i = ibegin , iend
      select case (str(i:i))
        case ('A':'Z')
          string(i:i) = char(iachar(str(i:i))+32)
        case default
      end select
    end do
  end function lower

  elemental pure function upper(str,istart,istop) result (string)
    implicit none
    character(*), intent(in) :: str
    character(len(str)) :: string
    integer , intent(in) , optional :: istart , istop
    integer :: i
    integer :: ibegin , iend
    string = str
    ibegin = 1
    if ( present(istart) ) then
      ibegin = max(ibegin,istart)
    end if
    iend = len_trim(str)
    if ( present(istop) ) then
      iend = min(iend,istop)
    end if
    do i = ibegin , iend
      select case (str(i:i))
        case ('a':'z')
          string(i:i) = char(iachar(str(i:i))-32)
        case default
      end select
    end do
  end function upper

  subroutine split(input_line,array,delimiters,order,nulls)
    implicit none
    character(len=*) , intent(in) :: input_line
    character(len=*) , optional, intent(in) :: delimiters
    character(len=*) , optional , intent(in) :: order
    character(len=*) , optional , intent(in) :: nulls
    character(len=:) , allocatable , intent(out) :: array(:)

    integer :: n
    integer , allocatable , dimension(:) :: ibegin
    integer , allocatable , dimension(:) :: iterm
    character(len=:) , allocatable :: dlim
    character(len=:) , allocatable :: ordr
    character(len=:) , allocatable  :: nlls
    integer :: ii , iiii
    integer :: icount
    integer :: ilen
    integer :: i , j , k
    integer :: icol
    integer :: idlim
    integer :: ifound
    integer :: inotnull
    integer :: ireturn
    integer :: imax

    if ( present(delimiters) ) then
      if ( delimiters /= '' ) then
        dlim = delimiters
      else
        dlim = ' '//char(9)//char(10)//char(11)//char(12)//char(13)//char(0)
      end if
    else
      dlim = ' '//char(9)//char(10)//char(11)//char(12)//char(13)//char(0)
    end if
    idlim = len(dlim)

    if ( present(order) ) then
      ordr = lower(adjustl(order))
    else
      ordr='sequential'
    end if
    if ( present(nulls) ) then
      nlls = lower(adjustl(nulls))
    else
      nlls = 'ignore'
    end if

    n = len(input_line)+1
    allocate(ibegin(n))
    allocate(iterm(n))
    ibegin(:) = 1
    iterm(:) = 1

    ilen = len(input_line)
    icount = 0
    inotnull = 0
    imax = 0

    select case (ilen)
      case (0)
      case default
        icol = 1
        parseloop: do k = 1 , ilen , 1
          ibegin(k) = icol
          if ( index(dlim(1:idlim),input_line(icol:icol)) == 0 ) then
            iterm(k) = ilen
            do i = 1 , idlim
              ifound = index(input_line(ibegin(k):ilen),dlim(i:i))
              if ( ifound > 0 ) then
                iterm(k) = min(iterm(k),ifound+ibegin(k)-2)
              end if
            end do
            icol = iterm(k)+2
            inotnull = inotnull+1
          else
            iterm(k) = icol-1
            icol = icol+1
          end if
          imax = max(imax,iterm(k)-ibegin(k)+1)
          icount = k
          if ( icol > ilen ) then
            exit parseloop
          end if
        end do parseloop
    end select

    select case ( trim(adjustl(nlls)) )
      case ('ignore', '', 'ignoreend')
        ireturn = inotnull
      case default
        ireturn = icount
    end select
    allocate(character(len=imax) :: array(ireturn))
    select case ( trim(adjustl(ordr)) )
      case ('reverse', 'right')
        ii = ireturn
        iiii = -1
      case default
        ii = 1
        iiii = 1
    end select

    do j = 1 , icount
      if ( iterm(j) < ibegin(j) ) then
        select case ( trim(adjustl(nlls)) )
          case ('ignore', '', 'ignoreend')
          case default
            array(ii) = ' '
            ii = ii+iiii
        end select
      else
        array(ii) = input_line(ibegin(j):iterm(j))
        ii = ii+iiii
      end if
    end do
  end subroutine split

  function basename(path,suffix) result(base)
    implicit none
    character(*), intent(In) :: path
    logical, intent(in), optional :: suffix
    character(:), allocatable :: base
    character(:), allocatable :: file_parts(:)
    logical :: with_suffix

    if ( .not. present(suffix) ) then
      with_suffix = .true.
    else
      with_suffix = suffix
    end if
    call split(path,file_parts,delimiters='\/')
    if ( size(file_parts) > 0 ) then
      base = trim(file_parts(size(file_parts)))
    else
      base = ''
    end if
    if ( .not. with_suffix ) then
      call split(base,file_parts,delimiters='.')
      if ( size(file_parts) >= 2 ) then
        base = trim(file_parts(size(file_parts)-1))
      end if
    end if
  end function basename

end program calc_heatindex

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
