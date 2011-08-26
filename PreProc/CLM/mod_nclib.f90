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

module mod_nclib

  use mod_memutil
  use netcdf
  use m_stdio
  use m_die
  use m_realkinds

  contains

!-----------------------------------------------------------------------
!     Purpose:
!        This routine creates a netCDF file.
!     Aguments :
!        filnam  char  input   input file name
!        cdfid  int  input   the id of the file to be closed.
!        ierr  int  output  indicates possible ierrs found in this
!                            routine.
!                            ierr = 0   no ierrs detected.
!                            ierr = 1   ierr detected.
!     History:
!        Nov. 91  PPM  UW  Created.
!-----------------------------------------------------------------------
  subroutine crecdf (filnam, cdfid, phymin, phymax, ndim, ierr) 

  implicit none

  integer , intent (in) :: ndim
  character(256), intent(in) :: filnam
  real(sp) , dimension(ndim) :: phymin , phymax
  integer , intent (out) :: ierr , cdfid

  integer , parameter :: maxdim = 4
  character(64) :: attnam
  character(1)  :: chrid(maxdim)
  integer :: k
  data chrid /'x','y','z','a'/

!     create the netCDF file
  ierr = nf90_create(filnam,nf90_clobber, cdfid)
  if ( ierr/=nf90_noerr ) go to 920

!     define global attributes
  do k = 1 , ndim
    attnam(1:3) = 'dom'
    attnam(4:4) = chrid(k)
    attnam(5:7) = 'min'
    attnam = attnam(1:7)
    ierr = nf90_put_att(cdfid,nf90_global,attnam,phymin(k))
    if ( ierr/=nf90_noerr ) go to 920

    attnam(1:3) = 'dom'
    attnam(4:4) = chrid(k)
    attnam(5:7) = 'max'
    attnam = attnam(1:7)
    ierr = nf90_put_att(cdfid,nf90_global,attnam,phymax(k))
    if ( ierr/=nf90_noerr ) go to 920
  end do

!     End variable definitions.
  ierr = nf90_enddef(cdfid)
  if ( ierr/=nf90_noerr ) go to 920

!     normal exit
  return

!     ierr exit
 920  write(stderr,*) 'ERROR: An error occurred while attempting to ', &
               'create the data file in subroutine crecdf.'
  write(stderr,*) nf90_strerror(ierr)
  ierr = nf90_close(cdfid)
  ierr = 1

  end subroutine crecdf

  subroutine rcrecdf (filnam,cdfid,varmin,varmax,ndim,ierr)

  implicit none

  character(256) , intent(in) :: filnam
  integer , intent(in) :: ndim
  integer , intent(out) :: cdfid
  integer , intent(out) :: ierr
  real(sp) , dimension(ndim) , intent (in) :: varmin , varmax

  call crecdf (filnam,cdfid,varmin,varmax,3,ierr)

  end subroutine rcrecdf

!-----------------------------------------------------------------------
!     Purpose:
!        This routine closes an open netCDF file.
!     Aguments :
!        cdfid  int  input   the id of the file to be closed.
!        ierr  int  output  indicates possible ierrs found in this
!                            routine.
!                            ierr = 0   no ierrs detected.
!                            ierr = 1   ierr detected.
!     History:
!        Nov. 91  PPM  UW  Created.
!-----------------------------------------------------------------------

  subroutine clscdf (cdfid, ierr)

  implicit none

!     Argument declarations.
  integer , intent(in) :: cdfid
  integer , intent(out) :: ierr

!     Close requested file.
  ierr = nf90_close(cdfid)

  end subroutine clscdf

! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! write data in array arr to netcdf-file for timestep vtstep
! 
! arguments an its meaning
!    cdfid     integer            ID returned by crecdf      
!    varnam    character*10       name of the variable to be written
!    arr       real*8(ie,je,ke)   array containing data to be written
!    ie,je,ke  integer            dimensions of arr as used in declar
!    iadim     integer(3)         array with used dimensions
!    vtstep    real*8             timestep to be written
!    lname     character*20       long name of variable
!    vunit     character*13       units of variable
!    factor    real*8             scaling factor for output
!    offset    real*8             offset for output
!    vvarmin   real*8(3)          minimum values of dim in phys. space
!    vvarmax   real*8(3)          maximum values of dim in phys. space
!    izstag    integer            1 for full levels   0 otherwise
!    vmisdat   real*8             value for missing data

  subroutine writecdf(cdfid,varnam,arr,ie,je,ke,iadim,vtstep,       &
       lname,vunit,factor,offset,vvarmin,vvarmax,xlat1d,xlon1d,sigh,&
       izstag,vmisdat,iotype)

  implicit none

  integer , intent(in) :: cdfid , ie , je , ke
  character(64) :: varnam , lname , vunit
  real(sp) , dimension(ie,je,ke) , intent(in) :: arr
  real(sp) , dimension(ie) , intent(in) :: xlon1d
  real(sp) , dimension(je) , intent(in) :: xlat1d
  real(sp) , dimension(ke) , intent(in) :: sigh
  integer , dimension(3) , intent(in) :: iadim
  real(sp) , dimension(3) , intent (in) :: vvarmin , vvarmax
  real(dp) , intent(in) :: vtstep
  real(sp) , intent(in) :: offset , factor , vmisdat
  integer , intent(in) :: izstag , iotype

!     declarations needed for netcdf-stuff

  real(sp) , dimension(3) :: varmin , varmax , varstg
  integer , dimension(3) :: vardim
  real(sp) :: tstep
  integer :: idtest , it , ndims
  real(sp) :: rfac

! *********** declare some auxiliary variables **********

  integer :: i , k , iputlev
  integer :: ievar , jevar , kevar

!     convert real*8 input variables to real*4 variables
  do i = 1 , 3
    varmin(i) = vvarmin(i)
    varmax(i) = vvarmax(i)
  end do

  tstep  = real(vtstep)
  ievar  = iadim(1)
  jevar  = iadim(2)
  kevar  = iadim(3)

!     set dimensions used in call to putdef and putdat
  ndims = 4
  vardim(1) = ievar
  vardim(2) = jevar
  vardim(3) = kevar
  varstg(1) = 0.
  varstg(2) = 0.
  varstg(3) = -0.5
  if ( izstag/=0 ) varstg(3) = 0.0

  it = nf90_inq_varid(cdfid,varnam,idtest)

  if ( it/=nf90_noerr ) then
    call putdefcdf(cdfid,varnam,ndims,vmisdat,lname,vunit,          &
                   offset,factor,vardim,varmin,varmax,iotype,it)
    call putcoords(cdfid,ndims,vardim,ievar,jevar,kevar,            &
                   xlat1d,xlon1d,sigh,it)
  end if

!     loop over all vertical levels and write the data
  do k = 1, kevar
    if ( kevar==1 ) then
      iputlev = 0
    else
      iputlev = k
    end if

    if ( iotype==1 ) then        ! Integer*2 format
      rfac = 1./factor
      call putdatcdfi2(cdfid, varnam, vtstep, k,              &
                       iputlev, ievar, jevar, arr, ie, je, ke,      &
                       vmisdat, rfac, offset, it) 
    else if ( iotype==2 ) then   ! Real*4 format
      call putdatcdfr4(cdfid,varnam,vtstep,k,iputlev,         &
                       ievar,jevar,arr,ie,je,ke,it)
    end if
  
  end do   ! loop over vertical levels

  end subroutine writecdf

!------------------------------------------------------------------------
!     Purpose:
!        Get all times on the specified NetCDF file
!     Arguments: 
!        cdfid  int  input   identifier for NetCDF file
!        times  real input   array contains all time values on the file,
!                            dimensioned at least times(ntimes)
!        ntimes int  input   number of times on the file
!        ierr  int  output  ierrflag 
!     History:
!        Heini Wernli, ETHZ  
!        Christoph Schaer, ETHZ
!     Note:
!        This preliminary version does not define the times-array, but only
!        overwrites or extends an existing times-array.
!------------------------------------------------------------------------

  subroutine putcoords(cdfid,ndim,vardim,ie,je,ke,xlat1d,xlon1d,    &
                       sigh,ierr)

  implicit none

  integer , intent(in) :: cdfid , ie , je , ke , ndim
  integer , dimension(ndim) , intent(in) :: vardim
  real(sp) , dimension(ie) :: xlon1d
  real(sp) , dimension(je) :: xlat1d
  real(sp) , dimension(ke) :: sigh
  integer , intent(out) :: ierr

  integer :: ncvid

  ierr = 0

  ierr = nf90_inq_varid(cdfid,'lon',ncvid)
  if ( ierr/=nf90_noerr) go to 920
  ierr = nf90_put_var(cdfid,ncvid,xlon1d)
  if ( ierr/=nf90_noerr) go to 920
  ierr = nf90_inq_varid(cdfid,'lat',ncvid)
  if ( ierr/=nf90_noerr) go to 920
  ierr = nf90_put_var(cdfid,ncvid,xlat1d)
  if ( ierr/=nf90_noerr) go to 920
  if ( vardim(3) > 1) then
    ierr = nf90_inq_varid(cdfid,'level',ncvid)
    if ( ierr/=nf90_noerr) go to 920
    ierr = nf90_put_var(cdfid,ncvid,sigh)
    if ( ierr/=nf90_noerr) go to 920
  end if

  return

!     ierr exit
 920  write(stderr,*) 'ERROR: An error occurred while attempting to ',  &
               'write variable dimension values in putcoords.'
  write(stderr,*) nf90_strerror(ierr)
  ierr = nf90_close(cdfid)
  ierr = 1

  end subroutine putcoords

!-----------------------------------------------------------------------
!     Purpose:
!        This routine is called to write the data of a variable
!        to an IVE-NetCDF file for use with the IVE plotting package. 
!        Prior to calling this routine, the file must be opened with 
!        a call to opncdf (for extension) or crecdf (for creation), the 
!        variable must be defined with a call to putdef.
!     Arguments:
!        cdfid   int   input   file-identifier
!                              (must be obtained by calling routine 
!                              opncdf or crecdf)
!        varnam  char  input   the user-supplied variable name (must 
!                              previously be defined with a call to
!                              putdef)
!        time    real  input   the user-supplied time-level of the
!                              data to be written to the file (the time-
!                              levels stored in the file can be obtained
!                              with a call to gettimes). If 'time' is not
!                              yet known to the file, a knew time-level is
!                              allocated and appended to the times-array.
!        level   int input     the horizontal level(s) to be written 
!                              to the NetCDF file. Suppose that the
!                              variable is defined as (nx,ny,nz,nt).
!                              level>0: the call writes the subdomain
!                                       (1:nx,1:ny,level,itimes)
!                              level=0: the call writes the subdomain
!                                       (1:nx,1:ny,1:nz,itimes)
!                              Here itimes is the time-index corresponding
!                              to the value of 'time'. 
!        dat     real  output  data-array dimensioned sufficiently 
!                              large. The dimensions (nx,ny,nz)
!                              of the variable must previously be defined
!                              with a call to putdef. No previous 
!                              definition of the time-dimension is
!                              required.
!        ierr   int output    indicates possible ierrs found in this
!                              routine.
!                              ierr = 0   no ierrs detected.
!                              ierr = 1   the variable is not present on
!                                          the file.
!                              ierr = 2   the value of 'time' is new, but
!                                          appending it would yield a non
!                                          ascending times-array.
!                              ierr = 3   inconsistent value of level
!                              ierr =10   another ierr.
!     History:
!       March 93    Heini Wernli (ETHZ)      Created wr2cdf.
!       April 93    Bettina Messmer (ETHZ)   Created putdat.
!-----------------------------------------------------------------------

  subroutine putdatcdfi2(cdfid, varnam, time, k, level,             &
                         ievar, jevar, arr, ie, je, ke,             &
                         vmisdat, rfac, offset, ierr)

  implicit none

  integer , intent(in) :: cdfid , ievar , jevar , ie , je , ke , k ,&
                          level
  character(64) , intent(in) :: varnam
  real(sp) , dimension(ie,je,ke) , intent(in) :: arr
  real(sp) , intent(in) :: vmisdat , rfac , offset
  integer , intent(out) :: ierr
  real(dp) , intent(in) :: time

  integer(2) , dimension(ie*je*ke) :: dat
  real(sp) :: misdat
  real(dp) :: timeval
  real(dp) , dimension(2) :: dvrange
  integer , dimension(4) :: corner , edgeln , did , vardim
  integer :: ndims , ntime
  integer :: idtime , idvar , iflag
  integer :: i , j , ik , ij
  integer , dimension(1) :: istart
 
  integer(2) , parameter :: shfill = -32767_2
  ierr = 0

  ij = 0
  do j = 1 , jevar
    do i = 1 , ievar
      ij = ij + 1
      if ( arr(i,j,k)<vmisdat ) then 
        dat(ij) = shfill  ! the lowest integer value, should match min
      else 
        dat(ij) = nint((arr(i,j,k)-offset)*rfac)
      end if
    end do
  end do

!     get definitions of data
  call getdefi2 (cdfid,varnam,ndims,misdat,vardim,ierr)
  if ( ierr/=0 ) go to 920

!     get id of variable
  ierr = nf90_inq_varid(cdfid,varnam,idvar)
  if ( ierr/=nf90_noerr ) go to 920

!     get times-array 
  ierr = nf90_inq_dimid(cdfid,'time',did(4))
  if ( ierr/=nf90_noerr ) go to 920
  ierr = nf90_inquire_dimension(cdfid,did(4),len=ntime)
  if ( ierr/=nf90_noerr ) go to 920
  ierr = nf90_inq_varid(cdfid,'time',idtime)
  if ( ierr/=nf90_noerr ) go to 920
!     Check if a new time step is starting
  iflag = 0
  do  i = 1 , ntime
    istart(1) = i
    ierr = nf90_get_var(cdfid,idtime,timeval,istart)
    if ( ierr/=nf90_noerr ) go to 920
    if ( time==timeval ) iflag = i
  end do
  if ( iflag==0 ) then         ! new time step
    ntime = ntime + 1
    iflag = ntime
    istart(1) = ntime
    ierr = nf90_put_var(cdfid,idtime,time,istart)
    if ( ierr/=nf90_noerr ) go to 920
    ierr = nf90_get_att(cdfid,idtime,'actual_range',dvrange)
    if ( ierr/=nf90_noerr ) go to 920
    if ( ( dvrange(1)>time ) .or. ( dvrange(1)==0. ) ) &
            dvrange(1) = time 
    if ( ( dvrange(2)<time ) .or. ( dvrange(2)==0. ) ) &
            dvrange(2) = time
    ierr = nf90_put_att(cdfid,idtime,'actual_range',dvrange)
    if ( ierr/=nf90_noerr ) go to 920
  end if

!     Define data volume to write on the NetCDF file in index space
  corner(1) = 1               ! starting corner of data volume
  corner(2) = 1
  edgeln(1) = vardim(1)       ! edge lengthes of data volume
  edgeln(2) = vardim(2)
  if ( level==0 ) then
    ik = 3
  else
    ik = 4
    corner(3) = level
    edgeln(3) = 1
  end if
  corner(ik) = iflag
  edgeln(ik) = 1
  
!     Put data on NetCDF file
  ierr = nf90_put_var(cdfid,idvar,dat,corner(1:ik),edgeln(1:ik))
  if ( ierr/=nf90_noerr ) go to 920

!     Synchronize output to disk and close the files
  ierr = nf90_sync(cdfid)
  if ( ierr/=nf90_noerr ) go to 920

  return

!     ierr exit
 920  write(stderr,*) 'ERROR: An error occurred while attempting to ', &
               'write variable ', varnam, ' values ',               &
               'at time ', time, ' in putdatcdfi2.'
  write(stderr,*) nf90_strerror(ierr)
  ierr = nf90_close(cdfid)
  ierr = 1

  end subroutine putdatcdfi2

!-----------------------------------------------------------------------
!     Purpose:
!        This routine is called to define the dimensions and the
!        attributes of a variable on an IVE-NetCDF file for use with the
!        IVE plotting package. Prior to calling this routine, the file must
!        be opened with a call to opncdf (extend an existing file) or
!        crecdf (create a new file).
!     Arguments:
!        cdfid   int   input   file-identifier
!                              (can be obtained by calling routine 
!                              opncdf)
!        varnam  char  input   the user-supplied variable name.
!        ndim    int   input   the number of dimensions (ndim<=4). 
!                              Upon ndim=4, the fourth dimension of the
!                              variable is specified as 'unlimited'
!                              on the file (time-dimension). It can 
!                              later be extended to arbitrary length.
!        misdat  real  input   missing data value for the variable. 
!        vardim  int   input   the dimensions of the variable.
!                              Is dimensioned at least Min(3,ndim). 
!        varmin  real  input   the location in physical space of the
!                              origin of each variable.
!                              Is dimensioned at least Min(3,ndim). 
!        varmax  real  input   the extent of each variable in physical
!                              space.
!                              Is dimensioned at least Min(ndim). 
!        ierr   int   output  indicates possible ierrs found in this
!                              routine.
!                              ierr = 0   no ierrs detected.
!                              ierr =10   other ierrs detected.
!     History:
!       Apr. 93    Christoph Schaer (ETHZ)     Created.
!-----------------------------------------------------------------------

  subroutine putdefcdf(cdfid,varnam,ndim,misdat,clname,             &
                       clunits,offset,xscale,vardim,varmin,varmax,  &
                       iotype,ierr)

  implicit none

  integer , parameter :: maxdim = 4

  integer , intent(in) :: cdfid , ndim , iotype
  character(64) , intent(in) :: varnam , clname , clunits
  integer , dimension(ndim) :: vardim
  real(sp) , dimension(ndim) :: varmin , varmax
  real(sp) , intent(in) :: xscale , offset , misdat
  integer , intent(out) :: ierr

  character(64) :: dimnam , dimchk
  character(64) , dimension(10) :: dimnams
  character(5) , dimension(maxdim) :: rdim
  integer , dimension(10) :: dimvals
  integer , dimension(maxdim) :: did
  integer :: numdims , numvars , numgats , dimulim
  integer :: id , idtime , i , k , ik
  integer :: idcoor
  real(sp) , dimension(2) :: vrange
  real(dp) , dimension(2) :: dvrange
  character(64) , dimension(maxdim) :: long_name
  character(64) , dimension(maxdim) :: units
  integer(2) , parameter :: shfill = -32767

  data rdim /'lon','lat','level','time'/
  data long_name /'Longitude','Latitude','Height_Index','Time'/

  data units / 'degrees_east', 'degrees_north',                     &
               'level' , 'hours since 1900-1-1 00:00:0.0' /


!     initialize vardim to something negative
  vardim(4) = -99

!     Make sure ndim <= maxdim.
  if ( ndim>maxdim ) then
     ierr = 10
     go to 900
  end if

!     Read existing dimensions-declarations from the file
  ierr = nf90_inquire(cdfid,numdims,numvars,numgats,dimulim)
  if ( ierr/=nf90_noerr ) go to 920

  if ( numdims>0 ) then
    do i = 1 , numdims
      ierr = nf90_inquire_dimension(cdfid,i,dimnams(i),dimvals(i))
      if ( ierr/=nf90_noerr ) go to 920
    end do
  end if

!     put file into define mode
  ierr = nf90_redef(cdfid)
  if ( ierr/=nf90_noerr ) go to 920


!     define spatial dimensions
  ik = 0 
  do k = 1 , max0(ndim,3)
    if ( vardim(k)>1 ) then
      ik = ik + 1
      dimnam = rdim(k)
      did(ik) = -1
      if ( numdims>0 ) then
!           check if an existing dimension-declaration can be used
!           instead of defining a new dimension
        do i = 1 , numdims
          dimchk = dimnams(i)
          if ( ( dimnam(1:3) == dimchk(1:3) ) ) then 
            did(ik) = i
            exit
          end if
        end do
      end if
      if ( did(ik)<0 ) then
!         define the dimension and an array with name of dimension
        ierr = nf90_def_dim(cdfid,dimnam,vardim(k),did(ik))
        if ( ierr/=nf90_noerr ) go to 920
        ierr = nf90_def_var(cdfid,dimnam,nf90_float,did(ik),idcoor)
        if ( ierr/=nf90_noerr ) go to 920
        vrange(1) = varmin(k)
        vrange(2) = varmax(k)
        ierr = nf90_put_att(cdfid,idcoor,'long_name',long_name(k))
        if ( ierr/=nf90_noerr ) go to 920
        ierr = nf90_put_att(cdfid,idcoor,'units',units(k))
        if ( ierr/=nf90_noerr ) go to 920
        ierr = nf90_put_att(cdfid,idcoor,'actual_range',vrange)
        if ( ierr/=nf90_noerr ) go to 920
      end if
    end if
  end do


!     define the times-array
  if ( ndim==4 ) then
    ik = ik + 1
!       define dimension 'time'
    ierr = nf90_inq_dimid(cdfid,'time',did(ik))
    if ( ierr/=nf90_noerr ) then 
!         this dimension must first be defined
      ierr = nf90_def_dim(cdfid,'time',nf90_unlimited,did(ik))
      if ( ierr/=nf90_noerr ) go to 920
    end if
!       define array 'time'
    ierr = nf90_inq_varid(cdfid,'time',idtime)
    if ( ierr/=nf90_noerr ) then 
!         define the times-array
      ierr = nf90_def_var(cdfid,'time',nf90_double,did(ik),idtime)
      if ( ierr/=nf90_noerr ) go to 920
      ierr = nf90_put_att(cdfid,idtime,'long_name',long_name(4))
      if ( ierr/=nf90_noerr ) go to 920
      ierr = nf90_put_att(cdfid,idtime,'units',units(4))
      if ( ierr/=nf90_noerr ) go to 920
      dvrange(1)=0.0D0
      dvrange(2)=0.0D0
      ierr = nf90_put_att(cdfid,idtime,'actual_range',dvrange)
      if ( ierr/=nf90_noerr ) go to 920
    end if
  end if

!     define variable
  if ( iotype==1 ) then
    ierr = nf90_def_var(cdfid,varnam,nf90_short,did(1:ik),id)
  else if (iotype==2) then
    ierr = nf90_def_var(cdfid,varnam,nf90_float,did(1:ik),id)
  end if
  if ( ierr/=nf90_noerr ) go to 920

!     define long_name
  ierr = nf90_put_att(cdfid,id,'long_name',clname)
  if ( ierr/=nf90_noerr ) go to 920
!     define units
  ierr = nf90_put_att(cdfid,id,'units',clunits)
  if ( ierr/=nf90_noerr ) go to 920
!     define missing data value
  if ( iotype==1 ) then
    ierr = nf90_put_att(cdfid,id,'missing_value',shfill)
  else if ( iotype==2 ) then
    ierr = nf90_put_att(cdfid,id,'missing_value',misdat)
  end if
  if ( ierr/=nf90_noerr ) go to 920

  if ( iotype==1 ) then
!       define offset_value
    ierr = nf90_put_att(cdfid,id,'add_offset',offset)
    if ( ierr/=nf90_noerr ) go to 920

!       define scale_factor
    ierr = nf90_put_att(cdfid,id,'scale_factor',xscale)
    if ( ierr/=nf90_noerr ) go to 920
  end if

!     leave define mode
  ierr = nf90_enddef(cdfid)
  if ( ierr/=nf90_noerr ) go to 920

!     synchronise output to disk and exit
  ierr = nf90_sync(cdfid)
  if ( ierr/=nf90_noerr ) go to 920
  return

!     Error exits.
 900  write(stderr,*) '*ERROR*: When calling putdefcdf, the number of ',&
               'variable dimensions must be less or equal 4.'
  ierr = nf90_close(cdfid)
  ierr = 1
  return

 920  write(stderr,*) '*ERROR*: An error occurred while attempting to ',&
               'write the data file in subroutine putdefcdf.'
  write(stderr,*) nf90_strerror(ierr)
  ierr = nf90_close(cdfid)
  ierr = 1
  return
  end subroutine putdefcdf

!-----------------------------------------------------------------------
!     Purpose:
!        This routine is called to write the data of a variable
!        to an IVE-NetCDF file for use with the IVE plotting package. 
!        Prior to calling this routine, the file must be opened with 
!        a call to opncdf (for extension) or crecdf (for creation), the 
!        variable must be defined with a call to putdef.
!     Arguments:
!        cdfid   int   input   file-identifier
!                              (must be obtained by calling routine 
!                              opncdf or crecdf)
!        varnam  char  input   the user-supplied variable name (must 
!                              previously be defined with a call to
!                              putdef)
!        time    real  input   the user-supplied time-level of the
!                              data to be written to the file (the time-
!                              levels stored in the file can be obtained
!                              with a call to gettimes). If 'time' is not
!                              yet known to the file, a knew time-level is
!                              allocated and appended to the times-array.
!        level   int input     the horizontal level(s) to be written 
!                              to the NetCDF file. Suppose that the
!                              variable is defined as (nx,ny,nz,nt).
!                              level>0: the call writes the subdomain
!                                       (1:nx,1:ny,level,itimes)
!                              level=0: the call writes the subdomain
!                                       (1:nx,1:ny,1:nz,itimes)
!                              Here itimes is the time-index corresponding
!                              to the value of 'time'. 
!        dat     real  output  data-array dimensioned sufficiently 
!                              large. The dimensions (nx,ny,nz)
!                              of the variable must previously be defined
!                              with a call to putdef. No previous 
!                              definition of the time-dimension is
!                              required.
!        ierr   int output    indicates possible ierrs found in this
!                              routine.
!                              ierr = 0   no ierrs detected.
!                              ierr = 1   the variable is not present on
!                                          the file.
!                              ierr = 2   the value of 'time' is new, but
!                                          appending it would yield a non
!                                          ascending times-array.
!                              ierr = 3   inconsistent value of level
!                              ierr =10   another ierr.
!     History:
!       March 93    Heini Wernli (ETHZ)      Created wr2cdf.
!       April 93    Bettina Messmer (ETHZ)   Created putdat.
!-----------------------------------------------------------------------

  subroutine putdatcdfr4(cdfid, varnam, time, k, level, ievar,      &
                         jevar, arr, ie, je, ke, ierr)

  implicit none

  integer , intent(in) :: cdfid
  character(64) , intent(in) :: varnam
  real(dp) , intent(in) :: time
  integer , intent(in) :: k , level , ievar , jevar , ie ,  je , ke
  real(sp) , dimension(ie,je,ke) , intent(in) :: arr
  integer , intent(out) :: ierr
  integer , dimension(1) :: istart

!     Declaration of local variables

  real(sp), dimension(ie*je*ke) :: dat
  real(sp) :: misdat
  real(dp) :: timeval
  real(dp) , dimension(2) :: dvrange

  integer , dimension(4) :: corner , edgeln , did , vardim
  integer :: ndims , ntime
  integer :: idtime , idvar , iflag
  integer :: i , j , ik , ij

  ij = 0
  do j = 1 , jevar
    do i = 1 , ievar
      ij = ij + 1
      dat(ij) = arr(i,j,k)
    end do
  end do

!     get definitions of data
  call getdefcdfr4(cdfid,varnam,ndims,misdat,vardim,ierr)
  if ( ierr/=0 ) go to 920

!     get id of variable
  ierr = nf90_inq_varid(cdfid,varnam,idvar)
  if ( ierr/=nf90_noerr ) go to 920

!     get times-array 
  ierr = nf90_inq_dimid(cdfid,'time',did(4))
  if ( ierr/=nf90_noerr ) go to 920
  ierr = nf90_inquire_dimension(cdfid,did(4),len=ntime)
  if ( ierr/=nf90_noerr ) go to 920
  ierr = nf90_inq_varid(cdfid,'time',idtime)
  if ( ierr/=nf90_noerr ) go to 920

!     Check if a new time step is starting
  iflag = 0
  do i = 1 , ntime
    istart(1) = i
    ierr = nf90_get_var(cdfid,idtime,timeval,istart)
    if ( ierr/=nf90_noerr ) go to 920
    if ( time==timeval ) iflag = i
  end do

  if ( iflag==0 ) then              ! new time step
    ntime = ntime+1
    iflag = ntime
    istart(1) = ntime
    ierr = nf90_put_var(cdfid,idtime,time,istart)
    if ( ierr/=nf90_noerr ) go to 920
    ierr = nf90_get_att(cdfid,idtime,'actual_range',dvrange)
    if ( ierr/=nf90_noerr ) go to 920
    if ( ( dvrange(1)>time ) .or. ( dvrange(1)==0. ) )              &
             dvrange(1) = time 
    if ( ( dvrange(2)<time ) .or. ( dvrange(2)==0. ) )              &
             dvrange(2) = time
    ierr = nf90_put_att(cdfid,idtime,'actual_range',dvrange)
    if ( ierr/=nf90_noerr ) go to 920
  end if

!     Define data volume to write on the NetCDF file in index space
  corner(1) = 1               ! starting corner of data volume
  corner(2) = 1
  edgeln(1) = vardim(1)       ! edge lengthes of data volume
  edgeln(2) = vardim(2)
  if ( level==0 ) then
    ik = 3
  else
    ik = 4
    corner(3) = level
    edgeln(3) = 1
  end if
  corner(ik) = iflag
  edgeln(ik) = 1
  
!     Put data on NetCDF file

  ierr = nf90_put_var(cdfid,idvar,dat,corner(1:ik),edgeln(1:ik))
  if ( ierr/=nf90_noerr ) go to 920

!     Synchronize output to disk and close the files

  ierr = nf90_sync(cdfid)
  if ( ierr/=nf90_noerr ) go to 920

  return

 920  write(stderr,*) 'ERROR: An error occurred while attempting to ',  &
               'write variable ', varnam, ' values ',               &
               'at time ', time, ' in putdatcdfr4.'
  write(stderr,*) nf90_strerror(ierr)
  ierr = nf90_close(cdfid)
  ierr = 1

  end subroutine putdatcdfr4

!-----------------------------------------------------------------------
!     Purpose:
!        This routine is called to get the dimensions and attributes of 
!        a variable from an IVE-NetCDF file for use with the IVE plotting
!        package. Prior to calling this routine, the file must be opened
!        with a call to opncdf.
!     Arguments:
!        cdfid   int   input   file-identifier
!                              (can be obtained by calling routine 
!                              opncdf)
!        varnam  char  input   the user-supplied variable name.
!                              (can be obtained by calling routine 
!                              opncdf)
!        ndim    int   output  the number of dimensions (ndim<=4)
!        misdat  real  output  missing data value for the variable. 
!        vardim  int   output  the dimensions of the variable.
!                              Is dimensioned at least (ndim). 
!        ierr   int   output  indicates possible ierrs found in this
!                              routine.
!                              ierr = 0   no ierrs detected.
!                              ierr = 1   the variable is not on the file.
!                              ierr =10   other ierrs.
!     History:
!       Apr. 93    Christoph Schaer (ETHZ)     Created.
!-----------------------------------------------------------------------

  subroutine getdefi2 (cdfid, varnam, ndim, misdat, vardim, ierr)

  implicit none

  integer , parameter :: maxdim = 4

  integer , intent(in) :: cdfid
  integer , intent(out) :: ndim
  character(64), intent(in) :: varnam
  integer, dimension(4) , intent(out) :: vardim
  real(sp) , intent(out) :: misdat
  integer , intent(out) :: ierr

  character(64) , dimension(maxdim) :: dimnam
  character(64) :: vnam
  integer :: id , i , k
  integer :: ndims , nvars , ngatts , recdim
  integer , dimension(maxdim) :: dimsiz
  integer :: vartyp , nvatts

!     inquire for number of dimensions
  ierr = nf90_inquire(cdfid,ndims,nvars,ngatts,recdim)
  if ( ierr/=nf90_noerr ) go to 920

!     read dimension-table
  do  i = 1 , ndims 
    ierr = nf90_inquire_dimension(cdfid,i,dimnam(i),dimsiz(i))
    if ( ierr/=nf90_noerr ) go to 920
  end do

!     get id of the variable
  ierr = nf90_inq_varid(cdfid,varnam,id)
  if ( ierr/=nf90_noerr ) go to 920

!     inquire about variable
  ierr = nf90_inquire_variable(cdfid,id,vnam,vartyp,ndim,          &
                                vardim,nvatts)
  if ( ierr/=nf90_noerr ) go to 920
  if ( vartyp/=nf90_short ) then
    ierr = 1
    go to 910
  end if

!     Make sure ndim <= maxdim.
  if ( ndim>maxdim ) then
     ierr = 1
     go to 900
  end if

!     get dimensions from dimension-table
  do k = 1 , ndim 
    vardim(k) = dimsiz(vardim(k))
  end do

!     get missing data value
  ierr = nf90_get_att(cdfid,id,'missing_value',misdat)
  if ( ierr/=nf90_noerr ) go to 920

!     normal exit
  return

!     Error exits.
 900  write(stderr,*) '*ERROR*: When calling getdefi2, the number of ', &
               'variable dimensions must be less or equal 4.'
  ierr = nf90_close(cdfid)
  ierr = 1
  return

 910  write(stderr,*) '*ERROR*: The selected var could not be found ',&
               'or is of wrong type in the file by getdefi2.'
  ierr = nf90_close(cdfid)
  ierr = 1
  return

 920  write(stderr,*) '*ERROR*: An error occurred while attempting to ',&
               'read the data file in subroutine getdefi2.'
  write(stderr,*) nf90_strerror(ierr)
  ierr = nf90_close(cdfid)
  ierr = 1

  end subroutine getdefi2

!-----------------------------------------------------------------------
!     Purpose:
!        This routine is called to get the dimensions and attributes of 
!        a variable from an IVE-NetCDF file for use with the IVE plotting
!        package. Prior to calling this routine, the file must be opened
!        with a call to opncdf.
!     Arguments:
!        cdfid   int   input   file-identifier
!                              (can be obtained by calling routine 
!                              opncdf)
!        varnam  char  input   the user-supplied variable name.
!                              (can be obtained by calling routine 
!                              opncdf)
!        ndim    int   output  the number of dimensions (ndim<=4)
!        misdat  real  output  missing data value for the variable. 
!        vardim  int   output  the dimensions of the variable.
!                              Is dimensioned at least (ndim). 
!        ierr   int   output  indicates possible ierrs found in this
!                              routine.
!                              ierr = 0   no ierrs detected.
!                              ierr = 1   the variable is not on the file.
!                              ierr =10   other ierrs.
!     History:
!       Apr. 93    Christoph Schaer (ETHZ)     Created.
!-----------------------------------------------------------------------

  subroutine getdefcdfr4(cdfid, varnam, ndim, misdat, vardim, ierr)

  implicit none

  integer , intent(in) :: cdfid
  character(64) , intent(in) :: varnam
  integer , intent(out) :: ndim
  integer, dimension(4) , intent(out) :: vardim
  real(sp) , intent(out) :: misdat
  integer , intent(out) :: ierr

  character(64), dimension(10) :: dimnam
  character(64) :: vnam
  integer :: id , i , k
  integer :: ndims , nvars , ngatts , recdim
  integer , dimension(10) :: dimsiz
  integer :: vartyp , nvatts

  ierr = 0

!     inquire for number of dimensions
  ierr = nf90_inquire(cdfid,ndims,nvars,ngatts,recdim)
  if ( ierr/=nf90_noerr ) then
    write(stderr,*) 'Error in nf90_inquire'
    go to 920
  end if

!     read dimension-table
  do  i = 1 , ndims 
    ierr = nf90_inquire_dimension(cdfid,i,dimnam(i),dimsiz(i))
    if ( ierr/=nf90_noerr ) then
      write(stderr,*) 'Error nf90_inquire_dimension ', dimnam(i)
      go to 920
    end if
  end do

!     get id of the variable
  ierr = nf90_inq_varid(cdfid,varnam,id)
  if ( ierr/=nf90_noerr ) then
    write(stderr,*) 'Error nf90_inq_varid ', varnam
    go to 920
  end if

!     inquire about variable
  ierr = nf90_inquire_variable(cdfid,id,vnam,vartyp,ndim,           &
                                vardim,nvatts)
  if ( ierr/=nf90_noerr ) then
    write(stderr,*) 'Error nf90_inquire_variable ', varnam
    go to 920
  end if
  if ( vartyp/=nf90_float ) then
    ierr = 1
    go to 910
  end if

!     Make sure ndim <= maxdim.
  if ( ndim>4 ) then
     ierr = 1
     go to 900
  end if

!     get dimensions from dimension-table
  do k = 1 , ndims
    vardim(k) = dimsiz(vardim(k))
  end do

!     get missing data value
  ierr = nf90_get_att(cdfid,id,'missing_value',misdat)
  if ( ierr/=nf90_noerr ) then
    write(stderr,*) 'Error in nf90_get_att missing_value'
    go to 920
  end if

!     normal exit
  return

!     Error exits.
 900  write(stderr,*) '*ERROR*: When calling getdefcdfr4, the number ', &
                  'of variable dimensions must be less or equal 4.'
  ierr = nf90_close(cdfid)
  ierr = 1
  return

 910  write(stderr,*) '*ERROR*: The selected variable could not be ', &
             'found or is of wrong type in the file by getdefcdfr4.'
  ierr = nf90_close(cdfid)
  ierr = 1
  return

 920  write(stderr,*) '*ERROR*: An error occurred while attempting to ',&
               'read the data file in subroutine getdefcdfr4.'
  write(stderr,*) nf90_strerror(ierr)
  ierr = nf90_close(cdfid)
  ierr = 1

  end subroutine getdefcdfr4

  subroutine readcdfr4(idcdf,vnam,lnam,units,nlon1,nlon,nlat1,nlat, &
                       nlev1,nlev,ntim1,ntim,vals)
 
  implicit none
!
  integer :: idcdf , nlat , nlat1 , nlev , nlev1 , nlon , nlon1 ,   &
             ntim , ntim1
  character(64) :: lnam , units , vnam
  real(sp) , dimension(nlon,nlat,nlev,ntim) :: vals
  intent (in) nlat , nlat1 , nlev , nlev1 , nlon , nlon1 , ntim ,   &
              ntim1
!
  integer , dimension(4) :: icount , istart
  integer :: iflag , invarid
!
  istart(1) = nlon1
  icount(1) = nlon
  istart(2) = nlat1
  icount(2) = nlat
  istart(3) = nlev1
  icount(3) = nlev
  istart(4) = ntim1
  icount(4) = ntim
 
  iflag = nf90_inq_varid(idcdf,vnam,invarid)
  if (iflag /= nf90_noerr) go to 920

  iflag = nf90_get_var(idcdf,invarid,vals,istart,icount)
  if (iflag /= nf90_noerr) go to 920

  iflag = nf90_get_att(idcdf,invarid,'long_name',lnam)
  if (iflag /= nf90_noerr) go to 920

  iflag = nf90_get_att(idcdf,invarid,'units',units)
  if (iflag /= nf90_noerr) go to 920
 
  return

 920  write(stderr,*) 'ERROR: An error occurred while attempting to ', &
               'read variable ', vnam , ' at time ', ntim1
  write(stderr,*) nf90_strerror(iflag)

  call die('readcdfr4','READ ERROR',1)

  end subroutine readcdfr4

  subroutine readcdfr4_360(idcdf,vnam,lnam,units,nlon1,nlon,nlat1,  &
                           nlat,nlev1,nlev,ntim1,ntim,nglon,nglat,  &
                           nglev,ngtim,vals)
 
  implicit none
!
  integer :: idcdf , nglat , nglev , nglon , ngtim , nlat , nlat1 , &
             nlev , nlev1 , nlon , nlon1 , ntim , ntim1
  character(64) :: lnam , units , vnam
  real(sp) , dimension(nlon,nlat,nlev,ntim) :: vals
  intent (in) nglat , nglev , nglon , ngtim , nlat , nlat1 , nlev , &
              nlev1 , nlon , nlon1 , ntim , ntim1
  intent (out) vals
!
  integer :: i , iflag , ii , ilon5 , invarid , j , jj , k , kk ,   &
             l , ll , nlat2 , nlev2 , nlon2 , ntim2
  integer , dimension(4) :: icount , istart
  real(sp) , pointer , dimension(:,:,:,:) :: vals1 , vals2
! 
  istart(1) = 1
  icount(1) = nglon
  istart(2) = 1
  icount(2) = nglat
  istart(3) = 1
  icount(3) = nglev
  istart(4) = 1
  icount(4) = ngtim
!     /*get variable and attributes*/
  ilon5 = nglon/2
  call getmem4d(vals1,1,nglon,1,nglat,1,nglev,1,ngtim,'mod_nclib:vals1')
  call getmem4d(vals2,1,nglon,1,nglat,1,nglev,1,ngtim,'mod_nclib:vals2')

  iflag = nf90_inq_varid(idcdf,vnam,invarid)
  if (iflag /= nf90_noerr) go to 920

  iflag = nf90_get_var(idcdf,invarid,vals1,istart,icount)
  if (iflag /= nf90_noerr) go to 920

  do l = 1 , ngtim
    do k = 1 , nglev
      do j = 1 , nglat
        do i = 1 , nglon
          vals2(i,j,k,l) = vals1(i,j,k,l)
        end do
      end do
    end do
  end do
  do l = 1 , ngtim
    do k = 1 , nglev
      do j = 1 , nglat
        do i = 1 , nglon
          ii = ilon5 + i
          if ( ii>nglon ) ii = ii - nglon
          vals1(ii,j,k,l) = vals2(i,j,k,l)
        end do
      end do
    end do
  end do
 
  ntim2 = ntim1 + ntim - 1
  nlev2 = nlev1 + nlev - 1
  nlat2 = nlat1 + nlat - 1
  nlon2 = nlon1 + nlon - 1
  do l = ntim1 , ntim2
    do k = nlev1 , nlev2
      do j = nlat1 , nlat2
        do i = nlon1 , nlon2
          ll = l - ntim1 + 1
          kk = k - nlev1 + 1
          jj = j - nlat1 + 1
          ii = i - nlon1 + 1
          vals(ii,jj,kk,ll) = vals1(i,j,k,l)
        end do
      end do
    end do
  end do
 
  iflag = nf90_get_att(idcdf,invarid,'long_name',lnam)
  if (iflag /= nf90_noerr) go to 920

  iflag = nf90_get_att(idcdf,invarid,'units',units)
  if (iflag /= nf90_noerr) go to 920
 
  return

 920  write(stderr,*) 'ERROR: An error occurred while attempting to ', &
               'read variable ', vnam , ' at time ', ntim1
  write(stderr,*) nf90_strerror(iflag)

  call die('readcdfr4_360','READ ERROR',1)

  end subroutine readcdfr4_360

  subroutine readcdfr4_iso(idcdf,vnam,lnam,units,nlon1,nlon,nlat1,  &
                           nlat,nlev1,nlev,ntim1,ntim,vals)
 
  implicit none
!
  integer :: idcdf , nlat , nlat1 , nlev , nlev1 , nlon , nlon1 ,   &
             ntim , ntim1
  character(64) :: lnam
  character(64) :: units
  character(64) :: vnam
  real(sp) , dimension(nlon,nlat,nlev,ntim) :: vals
  intent (in) nlat , nlat1 , nlev , nlev1 , nlon , nlon1 , ntim ,   &
              ntim1
  intent (out) vals
!
  integer , dimension(4) :: icount , istart
  integer , dimension(2) :: icount1 , istart1
  integer :: iflag , invarid
!
  istart1(1) = nlon1
  icount1(1) = nlon
  istart1(2) = nlat1
  icount1(2) = nlat
 
  istart(1) = nlon1
  icount(1) = nlon
  istart(2) = nlat1
  icount(2) = nlat
  istart(3) = nlev1
  icount(3) = nlev
  istart(4) = ntim1
  icount(4) = ntim
 
  iflag = nf90_inq_varid(idcdf,vnam,invarid)
  if (iflag /= nf90_noerr) go to 920

  iflag = nf90_get_var(idcdf,invarid,vals,istart,icount)
  if (iflag /= nf90_noerr) go to 920

  iflag = nf90_get_att(idcdf,invarid,'long_name',lnam)
  if (iflag /= nf90_noerr) go to 920

  iflag = nf90_get_att(idcdf,invarid,'units',units)
  if (iflag /= nf90_noerr) go to 920
 
  return

 920  write(stderr,*) 'ERROR: An error occurred while attempting to ',  &
               'read variable ', vnam , ' at time ', ntim1
  write(stderr,*) nf90_strerror(iflag)

  call die('readcdfr4_iso','READ ERROR',1)

  end subroutine readcdfr4_iso

end module mod_nclib
