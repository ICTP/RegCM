c * compiled with: pgf90 ncepwin.f -L../env/liblinux/ -lnetcdf
      PROGRAM RWDATA
      implicit none
      include 'ncepwin.param'
 
      integer nxwin, nywin
      parameter (nxwin=(x2win-x1win)/ddeg+1
     &         , nywin=(y2win-y1win)/ddeg+1)

C ****** Working variables
      character infil*50, outfil*50, lnam*50, units*20, fulldir*64
     &        , ayr*4
      integer*8 jul1900, julncep, julday, idatex
      integer iyr, imo, idy, ihr, iwin, jwin, kwin, i, j
     &      , iyrold, ndim, ierr, nlevmax
     &      , idcdf, idcdfout, ncopn, ncdid, ncvid, ncnowrit
      integer*2 igrid2d(nx,ny), igrid3d(nx,ny,nlev), imisdat
     &        , igridwin2d(nxwin,nywin), igridwin3d(nxwin,nywin,nlevwin)
      parameter (ndim=3)
      integer idim(ndim)
      real plevwin(nlev), varmin(3), varmax(3), scale, offset
     &   , xwin, ywin
      real*8 xhro

      do klev=1,nlevwin
        plevwin(klev) = plev(klev)
        print*,plevwin(klev)
      end do
      idim(1) = nxwin
      idim(2) = nywin
      idim(3) = nlevwin
      varmin(1) = x1win
      varmin(2) = y1win
      varmin(3) = plev(1)
      varmax(1) = x2win
      varmax(2) = y2win
      varmax(3) = plev(nlevwin)

      print*,varmin
      print*,varmax

      do ifld=1,nfld
        idatex = idate1
        CALL JULIAN(idatex,dt,julday,julncep,jul1900
     &	   , iyr,imo,idy,ihr,ayr )
        iyrold = 0
        print*,vnam(ifld),idatex
        do while (idatex.le.idate2)
          if (iyrold.ne.iyr) then
            fulldir = trim(infildir)//ayr//'/'
            infil = trim(fulldir)//trim(vnamfil(ifld))//'.'//ayr//'.nc'
            print*,'INPUT FILE: ',infil
            idcdf = NCOPN(infil,ncnowrit,ierr) ! Open input file
            outfil = trim(outfildir)//trim(vnamfil(ifld))
     &               //'.'//trim(outfilext)//'.'//ayr//'.nc'
            print*,'OUTPUT FILE: ',outfil
            CALL RCRECDF(outfil,idcdfout,varmin,varmax,ndim,ierr)
          end if
          xhro = float(jul1900)
          print*,'READING/WRITING DATA: ',idatex,julncep
          if (numdim(ifld).eq.3) then
            idim(3) = nlevwin
            CALL READCDF4D(idcdf,infil,julncep,vnam(ifld)
     &           , igrid3d,scale,offset,imisdat,lnam,units)
            do klev=1,nlevwin
            do jwin=1,nywin
            do iwin=1,nxwin
              xwin = x1win+(iwin-1)*ddeg
              ywin = y1win+(jwin-1)*ddeg
              kwin = nlevwin - klev + 1
              i = nint((x1win-x1)/ddeg)+iwin
              j = ny - (nint((y1win-y1)/ddeg)+jwin) + 1
              if (i.le.0) i=nint((x1win+360.-x1)/ddeg)+iwin
              igridwin3d(iwin,jwin,kwin) = igrid3d(i,j,klev)
            end do
            end do
            end do
            CALL WRITECOARDSI2(idcdfout,vnam(ifld),igridwin3d
     &         , nxwin,nywin,nlevwin,idim,xhro
     &         , lnam,units
     &         , varmin,varmax,plevwin,0,imisdat,scale,offset)
          else if (numdim(ifld).eq.2) then
            CALL READCDF3D(idcdf,infil,julncep,vnam(ifld)
     &           , igrid2d,scale,offset,imisdat,lnam,units)
            do jwin=1,nywin
            do iwin=1,nxwin
              idim(3) = 1
              xwin = x1win+(iwin-1)*ddeg
              ywin = y1win+(jwin-1)*ddeg
              i = nint((x1win-x1)/ddeg)+iwin
              j = ny - (nint((y1win-y1)/ddeg)+jwin) + 1
              if (i.le.0) i=nint((x1win+360.-x1)/ddeg)+iwin
              igridwin2d(iwin,jwin) = igrid2d(i,j)
            end do
            end do
            CALL WRITECOARDSI2(idcdfout,vnam(ifld),igridwin2d
     &         , nxwin,nywin,1,idim,xhro
     &         , lnam,units
     &         , varmin,varmax,1.,0,imisdat,scale,offset)
          end if

          idatex = idatex + dt
          iyrold = iyr
          CALL JULIAN(idatex,dt,julday,julncep,jul1900
     &	     , iyr,imo,idy,ihr,ayr )

          if (iyrold.ne.iyr) then
            print*,'CLOSING FILE'
            call ncclos (idcdf, ierr)
          end if

        end do ! time while loop
        print*,'CLOSING FILE'
        call ncclos (idcdf, ierr)

      end do ! ifld loop

      END

      SUBROUTINE GETMINMAX(f,n1,n2,n3,vmin,vmax,vmisdat)

      implicit none
      integer i,j,k,n1,n2,n3
      real f(n1,n2,n3),vmin,vmax,vmisdat

      vmin=9.e30
      vmax=-9.e30
      do k=1,n3
      do j=1,n2
      do i=1,n1
        if (f(i,j,k).ne.vmisdat) then
          if (f(i,j,k).gt.vmax) vmax=f(i,j,k)
          if (f(i,j,k).lt.vmin) vmin=f(i,j,k)
        end if
      end do
      end do
      end do

      return
      end


      SUBROUTINE JULIAN( mdate, dt, julday, julncep, jul1900
     &	, iyr, imo, idy, ihr, ayr )

      implicit none
      integer lenmon(12), jprev(12), dt, j, ileap
     &      , iyrm1, iths, ihun, idec, ione
      integer*8 mdate, jul1900, julncep, julday, idate
      integer iyr, imo, idy, ihr
      character*4 ayr

      data lenmon / 31, 28, 31, 30, 31, 30
     A            , 31, 31, 30, 31, 30, 31 /

      idate = mdate / 100
      iyr = idate / 10000
      iyrm1 = iyr - 1
      imo = ( idate - iyr*10000 ) / 100
      idy = mod( idate, 100 )
      ihr = mdate - idate*100

      ileap = mod( iyr, 4 )
      if(ileap.eq.0) lenmon(2) = 29

      if (ihr.gt.23) then
        ihr = ihr - 24
        idy = idy + 1
      end if
      if (idy.gt.lenmon(imo)) then
        idy = 1
        imo = imo + 1
      end if
      if (imo.gt.12) then
        imo = 1
        iyr = iyr + 1
      end if
      mdate = iyr*1000000 + imo*10000 + idy*100 + ihr

      jprev(1) = 0
      do j=2,12
        jprev(j)  = jprev(j-1) + lenmon(j-1)
      end do

      julday = idy + jprev(imo) - 1
      julncep = (julday+1)*24/dt - (24/dt-1) + ihr/dt
      jul1900 = ((iyr-1900)*365 + julday + int((iyrm1-1900)/4))*24 + ihr

      iths = iyr/1000
      ihun = (iyr-iths*1000)/100
      idec = (iyr-iths*1000-ihun*100)/10
      ione = (iyr-iths*1000-ihun*100-idec*10)
      ayr = char(iths+48)//char(ihun+48)//char(idec+48)//char(ione+48)

      return
      end


      SUBROUTINE READCDF4D(idcdf,infil,itim,varnam
     &         , igrid3d,scale,offset,imisdat,lnam,units)

      implicit none
      include 'ncepwin.param'
      character infil*(*),varnam*(*),dumnam*31,lnam*(*),units*(*)
	
      real*4 scale,offset,error
      integer*2 igrid3d(nx,ny,nlev),imisdat
      integer*4 istart(4),icount(4)
     &	 ,idcdf,iflag,latid,lonid,levid,timid,invarid
     &	 ,latsiz,lonsiz,levsiz,timsiz
     &      , idcdfout, ncopn, ncdid, ncvid, ncnowrit
      integer itim

c *** Check to make sure dimensions match ***
      latid = NCDID(idcdf, 'lat', iflag)
      CALL NCDINQ(idcdf, latid, dumnam, latsiz, iflag)
      lonid = NCDID(idcdf, 'lon', iflag)
      CALL NCDINQ(idcdf, lonid, dumnam, lonsiz, iflag)
      levid = NCDID(idcdf, 'level', iflag)
      CALL NCDINQ(idcdf, levid, dumnam, levsiz, iflag)
      timid = NCDID(idcdf, 'time', iflag)
      CALL NCDINQ(idcdf, timid, dumnam, timsiz, iflag)
      if (latsiz.ne.ny .or. lonsiz.ne.nx .or. itim.gt.timsiz) then
        print*,'DECLARED DIMENSIONS AND DATA DEMENSIONS DO NOT MATCH'
        print*,'  LAT: latsiz=',latsiz,'ny=',ny
        print*,'  LON: lonsiz=',lonsiz,'nx=',nx
        print*,'  LEV: levsiz=',levsiz,'nlev=',nlev
        print*,' TIME: timsiz=',timsiz,'julncep=',itim
        print*,'latid,lonid,levid,timid'
        print*,latid,lonid,levid,timid
        print*,varnam
        stop 999
      end if
      istart(timid) = itim
      istart(levid) = 1
      istart(lonid) = 1
      istart(latid) = 1
      icount(timid) = 1
      icount(levid) = nlev
      icount(lonid) = nx
      icount(latid) = ny

c  /*get variable and attributes*/
      invarid = NCVID(idcdf, varnam, iflag)
      CALL NCVGT(idcdf, invarid, istart, icount, igrid3d, iflag)
      CALL NCAGT (idcdf, invarid, 'add_offset', offset, iflag)
      CALL NCAGT (idcdf, invarid, 'scale_factor', scale, iflag)
      CALL NCAGT (idcdf, invarid, 'missing_value', imisdat, iflag)
      CALL NCAGTC (idcdf, invarid, 'long_name', lnam, 50, iflag)
      CALL NCAGTC (idcdf, invarid, 'units', units, 50, iflag)

      return
      end

      SUBROUTINE READCDF3D(idcdf,infil,itim,varnam
     &         , igrid2d,scale,offset,imisdat,lnam,units)

      implicit none
      include 'ncepwin.param'
      character infil*(*),varnam*(*),dumnam*31,lnam*(*),units*(*)
	
      real*4 scale,offset,error
      integer*2 igrid2d(nx,ny),imisdat
      integer*4 istart(3),icount(3)
     &	 ,idcdf,iflag,latid,lonid,timid,invarid
     &	 ,latsiz,lonsiz,timsiz
     &      , idcdfout, ncopn, ncdid, ncvid, ncnowrit
      integer itim

c *** Check to make sure dimensions match ***
      latid = NCDID(idcdf, 'lat', iflag)
      CALL NCDINQ(idcdf, latid, dumnam, latsiz, iflag)
      lonid = NCDID(idcdf, 'lon', iflag)
      CALL NCDINQ(idcdf, lonid, dumnam, lonsiz, iflag)
      timid = NCDID(idcdf, 'time', iflag)
      CALL NCDINQ(idcdf, timid, dumnam, timsiz, iflag)
      if (latsiz.ne.ny .or. lonsiz.ne.nx .or. itim.gt.timsiz) then
        print*,'DECLARED DIMENSIONS AND DATA DEMENSIONS DO NOT MATCH'
        print*,'  LAT: latsiz=',latsiz,'ny=',ny
        print*,'  LON: lonsiz=',lonsiz,'nx=',nx
        print*,' TIME: timsiz=',timsiz,'julncep=',itim
        print*,'latid,lonid,timid'
        print*,latid,lonid,timid
        print*,varnam
        stop 999
      end if
      istart(timid) = itim
      istart(lonid) = 1
      istart(latid) = 1
      icount(timid) = 1
      icount(lonid) = nx
      icount(latid) = ny

c  /*get variable and attributes*/
      invarid = NCVID(idcdf, varnam, iflag)
      CALL NCVGT(idcdf, invarid, istart, icount, igrid2d, iflag)
      CALL NCAGT (idcdf, invarid, 'add_offset', offset, iflag)
      CALL NCAGT (idcdf, invarid, 'scale_factor', scale, iflag)
      CALL NCAGT (idcdf, invarid, 'missing_value', imisdat, iflag)
      CALL NCAGTC (idcdf, invarid, 'long_name', lnam, 50, iflag)
      CALL NCAGTC (idcdf, invarid, 'units', units, 50, iflag)

      return
      end


c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WRITECOARDSI2(cdfid,varnam,irr
     &         , ie,je,ke,idim,vtstep
     &         , lname,unit
     &         , varmin,varmax,sigh,izstag,imisdat,factor,offset)

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c write data in array irr to netcdf-file for timestep vtstep
c 
c arguments an its meaning
c    cdfid     integer            ID returned by crecdf      
c    varnam    character*10       name of the variable to be written
c    irr       integer*2(ie,je,ke) array containing data to be written
c    ie,je,ke  integer            dimensions of irr as used in declar
c    idim      integer(3)         array with used dimensions
c    vtstep    real*8             timestep to be written
c    lname     character*20       long name of variable
c    unit      character*13       units of variable
c    factor    real*8             scaling factor for output
c    offset    real*8             offset for output
c    vvarmin   real*8(3)          minimum values of dim in phys. space
c    vvarmax   real*8(3)          maximum values of dim in phys. space
c    izstag    integer            1 for full levels   0 otherwise
c    vmisdat   real*8             value for missing data


      integer    cdfid,idim(3),izstag,ie,je,ke
      integer*2  irr(ie,je,ke)
      real*8     vtstep
      character*(*)  varnam
      character*(*)  lname
      character*(*)  unit

      integer*2       ufeld(ie*je*ke)

      integer*2       imisdat

      include '../env/include/netcdf.inc'

c     declarations needed for netcdf-stuff
      real*4        varmin(3),varmax(3),varstg(3),tstep,sigh(ke)
      integer       strbeg,strend
      integer       idtest,it
      integer       vardim(3)
      real*4        offset,factor

c *********** declare some auxiliary variables **********

      integer i,j,k,iputlev
      integer ievar,jevar,kevar

      tstep  = vtstep
      ievar  = idim(1)
      jevar  = idim(2)
      kevar  = idim(3)

c     set dimensions used in call to putdef and putdat
      ndims=4
      vardim(1)=ievar
      vardim(2)=jevar
      vardim(3)=kevar
      varstg(1)=0.
      varstg(2)=0.
      varstg(3)=-0.5
      if (izstag.ne.0) varstg(3)=0.

c     test if variable varnam already exists
c     if not then define it
      call ncpopt(NCVERBOS)


      idtest = ncvid(cdfid,varnam(strbeg(varnam):strend(varnam)),it)
      if (it.ne.0) then

        call putdefCOARDS(cdfid, varnam, ndims, imisdat, lname, unit
     &            ,  offset, factor, vardim, varmin, varmax, ierr)

        if (ierr.ne.0) stop '*ERROR* in putdef in writecdf'
        call putcoords(cdfid,ndims,vardim,varmin,varmax,sigh,ierr)
        if (ierr.ne.0) stop '*ERROR* in putdef in grbtst'
      endif

c     loop over all vertical levels and write the data
      do k=1,kevar

c     write data from inputarray to outputarray
        ij = 0
        do j=1,jevar
        do i=1,ievar
          ij = ij + 1
          ufeld(ij)=irr(i,j,k)
        enddo
        enddo
      
        if (kevar.eq.1) then
          iputlev=0
        else
          iputlev=k
        endif

        call putdatCOARDS(cdfid,varnam(strbeg(varnam):strend(varnam)),
     &         vtstep, iputlev, ufeld, ierr) 
        if (ierr.ne.0) then
          print *,'*ERROR* in putdat: variable ',
     &          varnam(strbeg(varnam):strend(varnam))
        endif
      
      enddo   ! loop over vertical levels

      end



      subroutine getdefCOARDS (cdfid, varnam, ndim, misdat, 
     &                              vardim, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to get the dimensions and attributes of 
c        a variable from an IVE-NetCDF file for use with the IVE plotting
c        package. Prior to calling this routine, the file must be opened
c        with a call to opncdf.
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (can be obtained by calling routine 
c                              opncdf)
c        varnam  char  input   the user-supplied variable name.
c                              (can be obtained by calling routine 
c                              opncdf)
c        ndim    int   output  the number of dimensions (ndim<=4)
c        misdat  real  output  missing data value for the variable. 
c        vardim  int   output  the dimensions of the variable.
c                              Is dimensioned at least (ndim). 
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   the variable is not on the file.
c                              error =10   other errors.
c     History:
c       Apr. 93    Christoph Schaer (ETHZ)     Created.
c-----------------------------------------------------------------------

      include  '../env/include/netcdf.inc'

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      character *(*) varnam
      integer        vardim(*), ndim, error, cdfid
      real           misdat

c     Local variable declarations.
      character *(20) dimnam(MAXNCDIM),vnam
      integer         id,i,k
      integer         ndims,nvars,ngatts,recdim,dimsiz(MAXNCDIM)
      integer         vartyp,nvatts, ncopts

c     Get current value of error options.
      call ncgopt (ncopts)

c     make sure NetCDF-errors do not abort execution
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     inquire for number of dimensions
      call ncinq(cdfid,ndims,nvars,ngatts,recdim,error)
      if (error.eq.1) goto 920

c     read dimension-table
      do i=1,ndims 
        call ncdinq(cdfid,i,dimnam(i),dimsiz(i),error)
        if (error.gt.0) goto 920
      enddo

c     get id of the variable
      id=ncvid(cdfid,varnam,error)
      if (error.eq.1) goto 910

c     inquire about variable
      call ncvinq(cdfid,id,vnam,vartyp,ndim,vardim,nvatts,error)
      if (vartyp.ne.NCSHORT) error=1
      if (error.gt.0) goto 920

c     Make sure ndim <= MAXDIM.
      if (ndim.gt.MAXDIM) then
         error = 1
         go to 900
      endif

c     get dimensions from dimension-table
      do k=1,ndim 
        vardim(k)=dimsiz(vardim(k))
      enddo

c     get missing data value
      call ncagt(cdfid,id,'missing_value',misdat,error)
      if (error.gt.0) goto 920     

c     normal exit
      call ncpopt (ncopts)
      return

c     Error exits.
 900  write (6, *) '*ERROR*: When calling getcdf, the number of ',
     &             'variable dimensions must be less or equal 4.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 910  write (6, *) '*ERROR*: The selected variable could not be found ',       
     &             'in the file by getcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 920  write (6, *) '*ERROR*: An error occurred while attempting to ',
     &             'read the data file in subroutine getcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      end



      subroutine putdefCOARDS(cdfid,varnam,ndim,misdat,clname,clunits
     &           , offset, scale, vardim, varmin, varmax, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to define the dimensions and the
c        attributes of a variable on an IVE-NetCDF file for use with the
c        IVE plotting package. Prior to calling this routine, the file must
c        be opened with a call to opncdf (extend an existing file) or
c        crecdf (create a new file).
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (can be obtained by calling routine 
c                              opncdf)
c        varnam  char  input   the user-supplied variable name.
c        ndim    int   input   the number of dimensions (ndim<=4). 
c                              Upon ndim=4, the fourth dimension of the
c                              variable is specified as 'unlimited'
c                              on the file (time-dimension). It can 
c                              later be extended to arbitrary length.
c        misdat  real  input   missing data value for the variable. 
c        vardim  int   input   the dimensions of the variable.
c                              Is dimensioned at least Min(3,ndim). 
c        varmin  real  input   the location in physical space of the
c                              origin of each variable.
c                              Is dimensioned at least Min(3,ndim). 
c        varmax  real  input   the extent of each variable in physical
c                              space.
c                              Is dimensioned at least Min(ndim). 
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error =10   other errors detected.
c     History:
c       Apr. 93    Christoph Schaer (ETHZ)     Created.
c-----------------------------------------------------------------------

      include  '../env/include/netcdf.inc'

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      character *(*) varnam,clname,clunits
      integer        vardim(*), ndim, error, cdfid
      real           varmin(*), varmax(*)
      real           scale,offset
      integer*2      misdat

c     Local variable declarations.
      character *(20) dimnam,dimchk
      character *(20) dimnams(MAXNCDIM)
      character *(5)  rdima(MAXDIM)
      integer         dimvals(MAXNCDIM)
      integer         numdims,numvars,numgats,dimulim
      integer         id,did(MAXDIM),idtime,i,k,ik,ierr
      integer         ncopts

      integer         ibeg,iend
      integer         idcoor
      real            vrange(2)
      real*8          dvrange(2)
      character*(20)  long_name(MAXDIM)
      character*(30)  units(MAXDIM)
      data            rdima /'lon  ','lat  ','layer','time '/
      data            long_name /'Longitude', 'Latitude', 
     &                           'Height_Index', 'Time' /

      data            units / 'degrees_east', 'degrees_north',
     &                  'millibar' , 'hours since 1900-1-1 00:00:0.0' /


c     External function declarations.
      integer         strbeg, strend

C ********************************
CJSP FIX FOR 2D/3D
      dimnam = '                    '

c     Get current value of error options.
      call ncgopt (ncopts)

c     make sure NetCDF-errors do not abort execution
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     Make sure ndim <= MAXDIM.
      if (ndim.gt.MAXDIM) then
         error = 10
         go to 900
      endif

c     Read existing dimensions-declarations from the file
      call ncinq(cdfid,numdims,numvars,numgats,dimulim,error)
      if (numdims.gt.0) then
        do i=1,numdims
          call ncdinq(cdfid,i,dimnams(i),dimvals(i),error)
c         print *,dimnams(i),dimvals(i)
        enddo
      endif

c     put file into define mode
      call ncredf(cdfid,error)
      if (error.ne.0) goto 920

c     define spatial dimensions
      ik = 0 
      do k=1,min0(3,ndim)
       if (vardim(k).gt.1) then
        ik = ik + 1
c       define the default dimension-name
        ibeg=strbeg(rdima(k))
        iend=strend(rdima(k))
        dimnam(1:1+iend-ibeg)=rdima(k)(ibeg:iend)
        did(ik)=-1
        if (numdims.gt.0) then
c         check if an existing dimension-declaration can be used
c         instead of defining a new dimension
          do i=1,numdims
            dimchk=dimnams(i)
            if ((dimnam(1:3).eq.dimchk(1:3))) then 
              did(ik)=i
              goto 100
            endif
          enddo
 100      continue
        endif
        if ((did(ik).lt.0)) then
c         define the dimension and an array with name of dimension
          did(ik)=ncddef(cdfid,dimnam,vardim(k),error)
          idcoor = ncvdef(cdfid,dimnam,NCFLOAT,1,did(ik),ierr)
          if ((error.ne.0).or.(ierr.ne.0)) goto 920
          vrange(1) = varmin(k)
          vrange(2) = varmax(k)
          iend=strend(long_name(k))
          call ncaptc(cdfid,idcoor,'long_name',
     &                 NCCHAR,iend,long_name(k)(1:iend),ierr)
          iend=strend(units(k))
          call ncaptc(cdfid,idcoor,'units',
     &                 NCCHAR,iend,units(k)(1:iend),ierr)
          call ncapt (cdfid,idcoor,'actual_range',
     &                 NCFLOAT,2,vrange,ierr)
        endif
       endif
      enddo

c     define the times-array
      if (ndim.eq.4) then
        ik = ik +1
c       define dimension 'time'
        did(ik)=ncdid(cdfid,'time',ierr)
        if (ierr.ne.0) then 
c         this dimension must first be defined
          did(ik) = ncddef (cdfid,'time',NCUNLIM,ierr)
        endif
c       define array 'time'
        idtime=ncvid(cdfid,'time',ierr)
        if (ierr.ne.0) then 
c         define the times-array
          idtime = ncvdef (cdfid,'time',NCDOUBLE,1,did(ik),ierr)
          call ncaptc(cdfid,idtime,'long_name',NCCHAR,4,'Time',ierr)
          iend = strend(units(4))
          call ncaptc(cdfid,idtime,'units',
     &                 NCCHAR,iend,units(4)(1:iend),ierr)
          dvrange(1)=0.
          dvrange(2)=0.
          call ncapt (cdfid,idtime,'actual_range',
     &                 NCDOUBLE,2,dvrange,ierr)
        endif
      endif

c     define variable
      id=ncvdef(cdfid,varnam,NCSHORT,ik,did,error)
      if (error.ne.0) goto 920

c     define long_name
      iend=strend(clname)
      call ncaptc(cdfid,id,'long_name',NCCHAR,iend,clname,error)
      if (error.gt.0) goto 920     

c     define units
      iend = strend(clunits)
      call ncaptc(cdfid,id,'units',NCCHAR,iend,clunits,error)
      if (error.gt.0) goto 920     

c     define missing data value
      call ncapt(cdfid,id,'missing_value',NCSHORT,1,misdat,error)
      if (error.gt.0) goto 920     

c     define offset_value
      call ncapt(cdfid,id,'add_offset',NCFLOAT,1,offset,error)
      if (error.gt.0) goto 920     

c     define scale_factor
      call ncapt(cdfid,id,'scale_factor',NCFLOAT,1,scale,error)
      if (error.gt.0) goto 920     

c     leave define mode
      call ncendf(cdfid,error)
      if (error.gt.0) goto 920     

c     synchronise output to disk and exit
      call ncsnc (cdfid,error)
      call ncpopt (ncopts)
      return

c     Error exits.
 900  write (6, *) '*ERROR*: When calling putcdf, the number of ',
     &             'variable dimensions must be less or equal 4.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 920  write (6, *) '*ERROR*: An error occurred while attempting to ',
     &             'write the data file in subroutine putcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return
      end



      function strbeg (string)
c-----------------------------------------------------------------------
c     Purpose:
c        This function returns the position of the first nonblank
c           character in the input string.  Returns 0 if the entire
c           string is blank.
c     Arguments:
c        string  char  input  string to be examined.
c     History:
c-----------------------------------------------------------------------

c     Function declaration
      integer strbeg

c     Argument declarations
      character*(*) string

c     Local variable declarations.
      integer i

c     External function declarations.
      integer len

      do 10 i = 1, len(string)
         strbeg = i
         if (string(i:i).ne.' ' .and. string(i:i).ne.char(0)) then
            return
         endif
 10   continue
      strbeg = 0
      end


      function strend (string)
c-----------------------------------------------------------------------
c     Purpose:
c        This function returns the position of the last nonblank
c           character in the input string.  Returns 0 if the entire
c           string is blank.
c     Arguments:
c        string  char  input  string to be examined.
c     History:
c-----------------------------------------------------------------------

c     Function declaration
      integer strend

c     Argument declarations
      character*(*) string

c     Local variable declarations.
      integer i

c     External function declarations.
      integer len

      do 10 i = len(string), 1, -1
         strend = i
         if (string(i:i) .ne. ' ' .and. string(i:i) .ne. char(0)) then
            return
         endif
 10   continue
      strend = 0
      end

      subroutine clscdf (cdfid, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine closes an open netCDF file.
c     Aguments :
c        cdfid  int  input   the id of the file to be closed.
c        error  int  output  indicates possible errors found in this
c                            routine.
c                            error = 0   no errors detected.
c                            error = 1   error detected.
c     History:
c        Nov. 91  PPM  UW  Created.
c-----------------------------------------------------------------------

      include  '../env/include/netcdf.inc'

c     Argument declarations.
      integer      cdfid, error

c     Local variable declarations.
      integer      ncopts

c     Get current value of error options.
      call ncgopt (ncopts)

c     Make sure netCDF errors do not abort execution.
      call ncpopt (NCVERBOS)

c     Close requested file.
      call ncclos (cdfid, error)

c     Reset error options.
      call ncpopt (ncopts)

      end


      subroutine putdatCOARDS(cdfid, varnam, time, level, dat, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to write the data of a variable
c        to an IVE-NetCDF file for use with the IVE plotting package. 
c        Prior to calling this routine, the file must be opened with 
c        a call to opncdf (for extension) or crecdf (for creation), the 
c        variable must be defined with a call to putdef.
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (must be obtained by calling routine 
c                              opncdf or crecdf)
c        varnam  char  input   the user-supplied variable name (must 
c                              previously be defined with a call to
c                              putdef)
c        time    real  input   the user-supplied time-level of the
c                              data to be written to the file (the time-
c                              levels stored in the file can be obtained
c                              with a call to gettimes). If 'time' is not
c                              yet known to the file, a knew time-level is
c                              allocated and appended to the times-array.
c        level   int input     the horizontal level(s) to be written 
c                              to the NetCDF file. Suppose that the
c                              variable is defined as (nx,ny,nz,nt).
c                              level>0: the call writes the subdomain
c                                       (1:nx,1:ny,level,itimes)
c                              level=0: the call writes the subdomain
c                                       (1:nx,1:ny,1:nz,itimes)
c                              Here itimes is the time-index corresponding
c                              to the value of 'time'. 
c        dat     real  output  data-array dimensioned sufficiently 
c                              large. The dimensions (nx,ny,nz)
c                              of the variable must previously be defined
c                              with a call to putdef. No previous 
c                              definition of the time-dimension is
c                              required.
c        error   int output    indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   the variable is not present on
c                                          the file.
c                              error = 2   the value of 'time' is new, but
c                                          appending it would yield a non
c                                          ascending times-array.
c                              error = 3   inconsistent value of level
c                              error =10   another error.
c     History:
c       March 93    Heini Wernli (ETHZ)      Created wr2cdf.
c       April 93    Bettina Messmer (ETHZ)   Created putdat.
c-----------------------------------------------------------------------

      include  '../env/include/netcdf.inc'

C     Declaration of local variables

      character*(*) varnam
      character*(20) chars
      integer cdfid


      integer*2 	dat(*)
      real	misdat
      real*8    time, timeval
      real*8    dvrange(2)

      integer	corner(4),edgeln(4),did(4),vardim(4),ndims
      integer	error, ierr
      integer	level,ntime
      integer	idtime,idvar,iflag
      integer	i,ik

      integer   strend
 
      call ncpopt(NCVERBOS)

c     get definitions of data
      call getdefCOARDS (cdfid, varnam(1:strend(varnam)), ndims, misdat, 
     &                           vardim, ierr)
      if (ierr.ne.0)  print *,'*ERROR* in getdef in putdat'

c     get id of variable
      idvar=ncvid(cdfid,varnam(1:strend(varnam)),ierr)
      if (ierr.ne.0) print *,'*ERROR* in ncvid in putdat'

c     get times-array 
      did(4)=ncdid(cdfid,'time',ierr)
      if (ierr.ne.0) print *,'*ERROR* did(4) in putdat'
      call ncdinq(cdfid,did(4),chars,ntime,ierr)
      if (ierr.ne.0) print *,'*ERROR* in ncdinq in putdat'
      idtime=ncvid(cdfid,'time',ierr)
      if (ierr.ne.0) print *,'*ERROR* in ncvid in putdat'
C     Check if a new time step is starting
      iflag=0
      do i=1,ntime
        call ncvgt1(cdfid,idtime,i,timeval,ierr)
        if (ierr.ne.0) print *,'*ERROR* in ncvgt1 in putdat'
        if (time.eq.timeval) iflag=i
      enddo
      if (iflag.eq.0) then		! new time step
        ntime=ntime+1
        iflag=ntime
        idtime=ncvid(cdfid,'time',ierr)
        if (ierr.ne.0) print *, '*ERROR* in ncvid in putdat'
        call ncvpt1(cdfid,idtime,ntime,time,ierr)
        if (ierr.ne.0) print *, '*ERROR* in ncvpt1 in putdat'
        call ncagt (cdfid,idtime,'actual_range',dvrange,ierr)
        if ((dvrange(1).gt.time).or.(dvrange(1).eq.0.))
     &         dvrange(1) = time 
        if ((dvrange(2).lt.time).or.(dvrange(2).eq.0.))
     &         dvrange(2) = time
        call ncapt (cdfid,idtime,'actual_range',NCDOUBLE,
     &                2,dvrange,ierr)
      endif

C     Define data volume to write on the NetCDF file in index space
      corner(1)=1               ! starting corner of data volume
      corner(2)=1
      edgeln(1)=vardim(1)       ! edge lengthes of data volume
      edgeln(2)=vardim(2)
      if (level.eq.0) then
        ik = 3
      else
        ik = 4
        corner(3)=level
        edgeln(3)=1
      endif
      corner(ik)=iflag
      edgeln(ik)=1
      
C     Put data on NetCDF file

c      print *,'vor Aufruf ncvpt d.h. Daten schreiben in putdat '
c      print *,'cdfid ',cdfid
c      print *,'idvar ',idvar
c      print *,'corner ',corner
c      print *,'edgeln ',edgeln
     
c      print *,'dat(1)=', dat(1)

      call ncvpt(cdfid,idvar,corner,edgeln,dat,error)
      if (error.ne.0) then
        print *, '*ERROR* in ncvpt in putdat - Put data on NetCDF file'
      endif

C     Synchronize output to disk and close the files

      call ncsnc(cdfid,ierr)
      if (ierr.ne.0) print *, '*ERROR* in ncsnc in putdat'
      end


      subroutine putcoords(cdfid,ndim,vardim,varmin,varmax,sigh,ierr)
C------------------------------------------------------------------------
C     Purpose:
C        Get all times on the specified NetCDF file
C     Arguments: 
C        cdfid  int  input   identifier for NetCDF file
C        times	real input   array contains all time values on the file,
C                            dimensioned at least times(ntimes)
C        ntimes int  input   number of times on the file
C        error  int  output  errorflag 
C     History:
C        Heini Wernli, ETHZ  
C        Christoph Schaer, ETHZ
C     Note:
C        This preliminary version does not define the times-array, but only
C        overwrites or extends an existing times-array.
C------------------------------------------------------------------------

      integer   vardim(*)
      integer	ierr,i,k,ndim,ii
      real      varmin(*),varmax(*),sigh(*)

      integer	cdfid,idcoor
      integer	ncvid
      real      actval, delta
      character*(5)  rdima(3)

      data            rdima /'lon  ','lat  ','layer'/

      do k=1,max0(ndim,3)
        if (vardim(k).gt.1) then
          idcoor = ncvid(cdfid,rdima(k),ierr)
	  if (k.ne.3) then
            if (ierr.eq.0) then
              actval = varmin(k)
              delta = (varmax(k)-varmin(k))/float(vardim(k)-1)
              do i=1,vardim(k)
                call ncvpt1(cdfid,idcoor,i,actval,ierr)
                actval=actval+delta
              enddo
            endif
	  else
            if (ierr.eq.0) then
              do i=1,vardim(k)
		ii = vardim(k)-i+1
		actval=sigh(ii)
                call ncvpt1(cdfid,idcoor,i,actval,ierr)
              enddo
            endif
	  end if
        endif
      enddo

      return
      end



c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RCRECDF (filnam,cdfid,varmin,varmax,ndim,ierr)
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      character*(*) filnam
      integer       cdfid,ierr
      real*4        varmin(ndim),varmax(ndim)

      call crecdf (filnam,cdfid,varmin,
     &               varmax,3,ierr)

      end



      subroutine crecdf (filnam, cdfid, phymin, phymax, ndim,error) 
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to create a netCDF file for use with 
c        the UWGAP plotting package.
c           Any netCDF file written to must be closed with the call
c        'call clscdf(cdfid,error)', where cdfid and error are
c        as in the argumentlist below. 
c     Arguments:
c        filnam  char  input   the user-supplied netCDF file name.
c        cdfid   int   output  the file-identifier
c        phymin  real  input   the minimum physical dimension of the
c                              entire physical domain along each axis.
c                              phymin is dimensioned (ndim)
c        phymax  real  input   the maximum physical dimension of the
c                              entire physical domain along each axis.
c                              phymax is dimensioned (ndim)
c        ndim    int   input   the number of dimensions in the file
c                              (i.e. number of elements in phymin,
c                              phymax)
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   error detected.
c     History:
c        Nov. 91  PPM  UW  Created cr3df.
c        Jan. 92  CS   UW  Created crecdf.
c-----------------------------------------------------------------------

      include  '../env/include/netcdf.inc'

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      integer        ndim, error
      character *(*) filnam
      real           phymin(*), phymax(*)

c     Local variable declarations.
      character *(20) attnam
      character *(1)  chrid(MAXDIM)
      integer         cdfid, k, ibeg, iend, lenfil, ncopts
      data            chrid/'x','y','z','a'/

c     External function declarations.
      integer        strbeg, strend

c     Get current value of error options, and make sure netCDF-errors do 
c     not abort execution
      call ncgopt (ncopts)
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     create the netCDF file
      cdfid = nccre (filnam(1:strend(filnam)), NCCLOB, error)
      if (error.ne.0) go to 920

c     define global attributes
      do k=1,ndim
        attnam(1:3)='dom'
        attnam(4:4)=chrid(k)
        attnam(5:7)='min'
        attnam=attnam(1:7)
        call ncapt(cdfid,NCGLOBAL,attnam,NCFLOAT,1,phymin(k),error)
        if (error.gt.0) goto 920

        attnam(1:3)='dom'
        attnam(4:4)=chrid(k)
        attnam(5:7)='max'
        attnam=attnam(1:7)
        call ncapt(cdfid,NCGLOBAL,attnam,NCFLOAT,1,phymax(k),error)
        if (error.gt.0) goto 920
      enddo

c     End variable definitions.
      call ncendf (cdfid, error)
      if (error.gt.0) goto 920

c     normal exit
      call ncpopt (ncopts)
      return

c     error exit
 920  write (6, *) 'ERROR: An error occurred while attempting to ',
     &             'create the data file in subroutine crecdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      end
