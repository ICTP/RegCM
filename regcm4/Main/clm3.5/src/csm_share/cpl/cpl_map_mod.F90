!===============================================================================
! SVN $Id: cpl_map_mod.F90 3620 2007-03-21 17:46:49Z robj $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/cpl/cpl_map_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_map_mod -- mapping subsystem module
!
! !DESCRIPTION:
!    This module represents a major subsystem of cpl6. "Mapping" refers 
!    to the interpolation of 2d field data from one domain/grid to another.  It is 
!    often desirable that maps have the properties of being {\it smooth} and 
!    {\it conservative}.  Common mapping techniques are bilinear interpolation
!    and area-averaging.  Mapping in cpl6 is implemented by a sparse matrix multiply.
!    This module defines the {\tt cpl\_map} data type which hold all the information
!    needed to complete a mapping of a {\it bundle} from one {\it domain} to another.
!    It also includes routines to do the mapping and initialize the 
!    {\tt cpl\_map} data structure and check the properties of the weights from
!    the sparse matrix.
!
! !REVISION HISTORY:
!    2001-Aug-14 - B. Kauffman -- gathered all mapping routines into this module
!    2001-May-20 - T. Craig -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

module cpl_map_mod

! !USES:

   use mct_mod            ! mct interface
   use cpl_domain_mod     ! data type & methods
   use cpl_bundle_mod     ! data type & methods
   use cpl_comm_mod       ! global data
   use cpl_kind_mod       ! kinds
   use cpl_control_mod, only: dbug=>cpl_control_infoDBug
   use cpl_control_mod, only: bfbflag=>cpl_control_bfbflag
   use shr_sys_mod        ! flush
   use shr_mpi_mod        ! mpi layer

   implicit none

   private   ! except

! !PUBLIC TYPES:

   public :: cpl_map
   
   type cpl_map
     character(CL)            :: name    ! text ID of mapping data
     type(mct_sMatP)          :: sMatP   ! the mct sparse matrix Plus data type
     character(3)             :: newtype ! intermediate domain type: src or dst ?
     integer(IN)              :: IDtype  ! 0=normal, 1=identity(ID)
     type(mct_Avect)          :: areasrc ! area of src grid from mapping file
     type(mct_Avect)          :: areadst ! area of dst grid from mapping file
   end type cpl_map

! !PUBLIC MEMBER FUNCTIONS:

   public :: cpl_map_init      ! initialize a map
   public :: cpl_map_clean     ! clean/dealloc a map
   public :: cpl_map_info      ! obtain information about a map
   public :: cpl_map_bun       ! map from one bundle to another
   public :: cpl_map_npFix     ! fix NP values wrt mapping vector fields

!--- npFixNew3 and npFixNew3R are the same except New3R saves some information
!    between calls. This means it has to be always called for the same grids
!    (not fields).  New3 is totally general.  New3R improves performance
!    on the Cray X1E by reducing string work. ---
!--- npFixNew4R is a lower memory version of 3R

   interface cpl_map_npFix; module procedure cpl_map_npFixNew4R; end interface
!  interface cpl_map_npFix; module procedure cpl_map_npFixNew3R; end interface
!  interface cpl_map_npFix; module procedure cpl_map_npFixNew3; end interface
!  interface cpl_map_npFix; module procedure cpl_map_npFixNone; end interface

! !PUBLIC DATA MEMBERS:

   character(*),parameter,public :: cpl_map_areaAV_field = 'aream'

!EOP

   !--- module variables ---
   character(*),parameter :: modName = "cpl_map_mod"

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_init - Create a map between two domains
!
! !DESCRIPTION:
!    Initialize a {\tt map\_X} to interpolate data from {\it domain} {\tt dom\_src}
!    to {\it domain} {\tt dom\_dst}.  {\tt map\_X} is assigned the name {\tt mapName}.
!
!    Mapping weights are read from the file {\tt fileName}
! 
!    {\tt newdom} is either ``src'' or ``dst'' and specifies if the communication
!    needed to complete the mapping is done before, ``src''-based, or after,
!    ``dst''-based, the mapping.
!
!    If {\tt cpl\_control\_infoDBug} is greater than 1, then this routine
!    will perform consistency checks on the mapping weights.
!
! !REVISION HISTORY:
!    2006-Aug-06 - R. Jacob - change to use MCT SparseMatrixPlus
!    2001-Jun-14 - T. Craig - first functioning version
!    2001-Mar-20 - T. Craig, B. Kauffman, R. Jacob - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_map_init(map_X,dom_src,dom_dst,mapName,fileName,newdom)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_map)   ,intent(out)        :: map_X    ! map_X data
   type(cpl_domain),intent( in),target :: dom_src  ! map's source domain
   type(cpl_domain),intent( in),target :: dom_dst  ! map's destination domain
   character(*),intent( in)            :: mapName  ! map's ID string 
   character(*),intent( in)            :: fileName ! file containing map data
   character(*),intent( in)            :: newdom   ! which domain to alter (for mct)

!EOP

   !--- local ---
   type(mct_sMat ) :: sMati    ! initial sMat from read (either root or decomp)
   type(mct_Avect) :: areasrc0 ! area of src grid from mapping file
   type(mct_Avect) :: areadst0 ! area of dst grid from mapping file
   integer         :: rCode

   logical, parameter  :: MAP_READD = .true.   ! turns on low mem weights reader

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_init) '
   character(*),parameter :: F00 = "('(cpl_map_init) ',8a)"
   character(*),parameter :: F01 = "('(cpl_map_init) ',2(a,i8))"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

    write(6,F00) "initialize map: ",trim(mapName)

   !----------------------------------------------------------------------------
   ! read & test the SparseMatrix data on root processor only 
   !----------------------------------------------------------------------------
   if (MAP_READD) then
      call mct_sMat_readdnc(sMati, dom_src%gsMap, dom_dst%gsMap, newdom, &
           areasrc0,areadst0,fileName,cpl_comm_comp_pid, cpl_comm_comp)
   else
      if (cpl_comm_comp_pid == 0) then
         call cpl_map_read0(sMati,areasrc0,areadst0,fileName)
      endif
      call shr_mpi_barrier(cpl_comm_comp,subName//" barrier")
   endif

   !----------------------------------------------------------------------------
   ! scatter map_X data and create new/intermediate domain as required by MCT
   !----------------------------------------------------------------------------
   map_X%name    =  trim(mapName) // ", source file = " // trim(fileName)
   map_X%newtype =  newdom
   map_X%IDtype  =  0

   if (bfbflag) then
     map_X%newtype = 'src'
     write(6,*) subName,': bfbflag = ',bfbflag, &
       ': overwriting newtype from ',newdom,' to ',map_X%newtype
   endif

   call mct_aVect_clean(map_X%areasrc)
   call mct_aVect_scatter(areasrc0,map_X%areasrc, dom_src%gsMap,0,cpl_comm_comp,rCode)
   call mct_aVect_clean(map_X%areadst)
   call mct_aVect_scatter(areadst0,map_X%areadst, dom_dst%gsMap,0,cpl_comm_comp,rCode)
   if (dbug > 2) then
      if (cpl_comm_comp_pid == 0) then
         write(6,*) subName,'lSize of src & dest',&
                    mct_aVect_lSize(areasrc0), &
                    mct_aVect_lSize(areadst0)
         write(6,*) subName,'min/max src ',minval(areasrc0%rAttr(1,:)), &
                                           maxval(areasrc0%rAttr(1,:))
         write(6,*) subName,'min/max dst ',minval(areadst0%rAttr(1,:)), &
                                           maxval(areadst0%rAttr(1,:))
      endif
      call mct_aVect_info(4,map_X%areasrc,cpl_comm_comp,cpl_comm_comp_pid,istr='map areasrc')
      call mct_aVect_info(4,map_X%areadst,cpl_comm_comp,cpl_comm_comp_pid,istr='map areadst')
   endif


   if (MAP_READD) then
      if (map_X%newtype == "dst") then !--- create new/intermediate destination domain ---
         write(6,F00) 'build dst sMatPlus ...'  ; call shr_sys_flush(6)
         call mct_sMatP_init(map_X%sMatP, sMati, dom_src%gsMap, dom_dst%gsMap, &
                                 0,cpl_comm_comp,cpl_comm_mph_cid)
      elseif (map_X%newtype == "src") then !--- create new/intermediate source domain ---
         write(6,F00) 'build src sMatPlus...' ; call shr_sys_flush(6)
         call mct_sMatP_init(map_X%sMatP, sMati, dom_src%gsMap, dom_dst%gsMap, &
                                 0,cpl_comm_comp,cpl_comm_mph_cid)
      else !--- no other valid choices ---
         write(6,F00) 'ERROR: invalid newdom value = ',map_X%newtype
         call shr_sys_abort(trim(subName)//" invalid newdom value")
      endif

      if (dbug > 1 .or. rCode /= 0) then
         write(6,F01) 'scattered sMat rows x cols =',map_X%sMatP%Matrix%nrows,' x', &
         map_X%sMatP%Matrix%ncols
         write(6,F01) 'scattered sMat lSize  ',mct_sMat_lSize(map_X%sMatP%Matrix)
         write(6,F01) 'scattered sMat GnumEl ',mct_sMat_GNumEl(map_X%sMatP%Matrix, &
                                                                         cpl_comm_comp)
      end if
   else
      if (map_X%newtype == "dst") then !--- create new/intermediate destination domain ---
         write(6,F00) 'build dst sMatPlus ...'  ; call shr_sys_flush(6)
         call mct_sMatP_init(map_X%sMatP, sMati, dom_src%gsMap, dom_dst%gsMap, &
                                 'Yonly',0,cpl_comm_comp,cpl_comm_mph_cid)

         if (dbug > 1 .or. rCode /= 0) then
            write(6,F01) 'scattered sMat rows x cols =',map_X%sMatP%Matrix%nrows,' x', &
            map_X%sMatP%Matrix%ncols
            write(6,F01) 'scattered sMat lSize  ',mct_sMat_lSize(map_X%sMatP%Matrix)
            write(6,F01) 'scattered sMat GnumEl ',mct_sMat_GNumEl(map_X%sMatP%Matrix, &
                                                                            cpl_comm_comp)
         end if
      elseif (map_X%newtype == "src") then !--- create new/intermediate source domain ---
         write(6,F00) 'build src sMatPlus...' ; call shr_sys_flush(6)
         call mct_sMatP_init(map_X%sMatP, sMati, dom_src%gsMap, dom_dst%gsMap, &
                                 'Xonly',0,cpl_comm_comp,cpl_comm_mph_cid)
         if (dbug > 1 .or. rCode /= 0) then
            write(6,F01) 'scattered sMat rows x cols =',map_X%sMatP%Matrix%nrows,' x', &
            map_X%sMatP%Matrix%ncols
            write(6,F01) 'scattered sMat lSize  ',mct_sMat_lSize(map_X%sMatP%Matrix)
            write(6,F01) 'scattered sMat GnumEl ',mct_sMat_GNumEl(map_X%sMatP%Matrix, &
                                                                            cpl_comm_comp)
         end if
      else !--- no other valid choices ---
         write(6,F00) 'ERROR: invalid newdom value = ',map_X%newtype
         call shr_sys_abort(trim(subName)//" invalid newdom value")
      endif
   endif

#ifdef CPP_VECTOR
   !--- initialize the vector parts of the sMat ---
   call mct_sMatP_Vecinit(map_X%sMatP)
#endif

   call mct_sMat_clean(sMati)
   if (cpl_comm_comp_pid == 0) then
     call mct_aVect_clean(areasrc0)
     call mct_aVect_clean(areadst0)
   endif

   call cpl_map_MatPtest(map_X%sMatP)

   write(6,F00) "done initializing map: ",trim(mapName) ; call shr_sys_flush(6)
   call shr_sys_flush(6)

end subroutine cpl_map_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_clean - Deallocate a map data type
!
! !DESCRIPTION:
!    Deallocate all memory associated with the input map type {\tt mapping}.
! 
! !REVISION HISTORY:
!    2002-Jan-20 - T. Craig - first functioning version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_map_clean(mapping)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_map),intent(inout)    :: mapping ! mapping data

!EOP

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_clean) '
   character(*),parameter :: F00 = "('(cpl_map_clean) ',8a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   if (dbug>2) write(6,F00) 'cleaning map, name=',trim(mapping%name)

   mapping%name    = '<clean>'
   mapping%newtype = '<clean>'
   call mct_sMat_clean(mapping%sMatP%Matrix)
   call mct_aVect_clean(mapping%areasrc)
   call mct_aVect_clean(mapping%areadst)
   mapping%IDtype = 0

end subroutine cpl_map_clean

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_info - Print inforomation about the map
!
! !DESCRIPTION:
!    Print information about the map {\tt mapping} to stdout.
! 
! !REVISION HISTORY:
!    2002-Jan-14 - T. Craig - first functioning version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_map_info(mapping)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_map)   ,intent(in)                  :: mapping ! mapping data

!EOP

   !--- local ---
   integer(IN)         :: i,j,k,n      ! generic indicies
   real(R8)            :: mn,mx        ! min,max

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_info) '
   character(*),parameter :: F00 = '(a,a,2f12.3)'

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   write(6,*) ' '
   write(6,*) subName,' mapping%name:',trim(mapping%name)
   write(6,*) subName,' mapping%newtype:',trim(mapping%newtype)
   write(6,*) subName,' mapping%IDtype:',mapping%IDtype

   write(6,*) subName,' mapping%sMat: '
   write(6,*) subName,' mapping%sMatP%Matrix%nrows: ', &
                                  mapping%sMatP%Matrix%nrows
   write(6,*) subName,' mapping%sMatP%Matrix%ncols: ', &
                                  mapping%sMatP%Matrix%ncols

   call mct_aVect_info(2,mapping%sMatP%Matrix%data, &
                           cpl_comm_comp,cpl_comm_comp_pid)
   write(6,*) ' '
   call shr_sys_flush(6)

end subroutine cpl_map_info

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_bun - map a bundle from one domain to domain
!
! !DESCRIPTION:
!    Map input bundle {\tt buni} on one {\it domain} to output bundle {\tt buno}
!    on a different {\it domain} using the map {\tt mapx}.
!
!    All attributes in {\tt buni} and {\tt buno} which have the same name will
!    be mapped.  Data in {\tt buno} will be overwritten.
!
!    Both {\tt buni} and {\tt buno} must be initialized before calling this
!    routine
!
!    If the set of optional arguments {\tt bunfs}, {\tt fsname},
!    {\tt bunfd}, {\tt fsname} are all present then the data in {\tt buni} will be
!    multiplied by field {\tt fsname} from {\tt bunfs} before mapping
!    (Note: the data in {\tt buni} will NOT be altered) and the
!    data in {\tt buno} will be multipled by field {\tt fdname}
!    in bundle {\tt bunfs} after mapping.
!
!    If optional argument {\tt mvector} controls the use of the
!    vector-computer freindly versions of the mapping routine.
!    This can be used to override the default settings.
!    
!
! !REVISION HISTORY:
!    20May01 - T. Craig -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_bun(buni,buno,mapx,bunfs,fsname,bunfd,fdname,mvector)

! !USES:

   use shr_timer_mod       ! share timer routines

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout)       :: buni    ! input bundle
   type(cpl_bundle),intent(out)         :: buno    ! output bundle
   type(cpl_map)   ,intent(inout)       :: mapx    ! mapping between two domains
   type(cpl_bundle),intent(in),optional :: bunfs   ! src fraction input bundle
   character(*)    ,intent(in),optional :: fsname  ! name of field in bunfs
   type(cpl_bundle),intent(in),optional :: bunfd   ! dst fraction input bundle
   character(*)    ,intent(in),optional :: fdname  ! name of field in bunfd
   logical         ,intent(in),optional :: mvector  ! enable vector-friendly mapping

!EOP

   !--- local ---
   type(cpl_bundle):: buni_local  ! local copy of buni for optional arguments
   character(CL)   :: name        ! name string for temp bundle
   logical         :: normalize   ! true if optional arguments present
   logical         :: usevector   ! true if we want to use the vector-friendly code
   integer(IN)     :: n,m         ! generic indicies
   integer(IN)     :: npts        ! number of points (local) in an aVect field
   integer(IN)     :: nfld        ! number of fields (local) in an aVect field
   integer(IN)     :: kfld        ! field number of fld1 in bun1
   real(R8)        :: fdr         ! 1/bunfd

   logical,save :: first_call = .true.     ! first time in subroutine
   integer(IN),save :: tmap1,tmap2,tmap3,tmap4,tmap5,tmap6
   integer(IN),save :: tmap7,tmap8,tmap9,tmap10

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_bun) '
   character(*),parameter :: F00 = '("(cpl_map_bun) ",8a)'
   character(*),parameter :: F01 = '("(cpl_map_bun) ",3a,i5)'

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   if (first_call) then
     first_call = .false.
     call shr_timer_get(tmap1,'tmap1')
     call shr_timer_get(tmap2,'tmap2')
     call shr_timer_get(tmap3,'tmap3')
     call shr_timer_get(tmap4,'tmap4')
     call shr_timer_get(tmap5,'tmap5')
     call shr_timer_get(tmap6,'tmap6')
     call shr_timer_get(tmap7,'tmap7')
     call shr_timer_get(tmap8,'tmap8')
     call shr_timer_get(tmap9,'tmap9')
     call shr_timer_get(tmap10,'tmap10')
   endif

   !----------------------------------------------------------------------------
   ! By default, turn on vector parts if CPP_VECTOR is defined
   ! But behaviour can be changed by input argument.
   !----------------------------------------------------------------------------
#ifdef CPP_VECTOR
   usevector = .true.
#else
   usevector = .false.
#endif
   !--- override the defaults with the input argument ---
   if(present(mvector)) then
     usevector = mvector
   endif

   !----------------------------------------------------------------------------
   ! identity map => "mapping is a simply copy the bundle
   !----------------------------------------------------------------------------
   if (mapx%IDtype == 1) then
     call cpl_bundle_copy(buni,outbun=buno)
     return
   endif

   !----------------------------------------------------------------------------
   ! if normalization is requested, create a local copy of bunlde
   !----------------------------------------------------------------------------
   call shr_timer_start(tmap1)
   normalize = .false.
   if (present(bunfs) .and.present(bunfd).and.    &
       present(fsname).and.present(fdname)) then
     normalize = .true.
   endif
   if ((present(bunfs).and..not.present(bunfd)) .or.  &
       (present(bunfd).and..not.present(bunfs)) .or.  &
       (present(bunfs).and..not.present(fsname)) .or. &
       (present(bunfd).and..not.present(fdname))) then
     write(6,*) subName,' ERROR: optional arguments inconsistent '
     call shr_sys_abort(subName)
   endif
   call shr_timer_stop(tmap1)

   call shr_timer_start(tmap2)
   if (normalize) then
     name = trim(buni%name) // "_local"
     call cpl_bundle_initv(buni_local,name,buni,buni%dom)
!     call cpl_bundle_copy(buni,outbun=buni_local)
!     call cpl_bundle_mult(buni_local,bunfs,fsname)
      call cpl_bundle_mult(buni_local,bunfs,fsname,initbun=buni)
   endif
   call shr_timer_stop(tmap2)

   !----------------------------------------------------------------------------
   ! do the mapping 
   !----------------------------------------------------------------------------

   !--- map the data ---
   call shr_timer_start(tmap7)
   if (normalize) then
     call mct_Smat_AvMult(buni_local%data,mapx%sMatP,buno%data,vector=usevector)
   else
     call mct_Smat_AvMult(buni%data      ,mapx%sMatP,buno%data,vector=usevector)
   endif
   call shr_timer_stop(tmap7)


   !----------------------------------------------------------------------------
   ! finish normalization
   !----------------------------------------------------------------------------
   call shr_timer_start(tmap9)
   if (normalize) then
     call cpl_bundle_divide(buno,bunfd,fdname)
     call cpl_bundle_clean(buni_local)
   endif
   call shr_timer_stop(tmap9)

   if (buni%cnt /= 1) write(6,F01) &
       "WARNING:  bundle ",trim(buni%name)," has accum count =",buni%cnt
   buno%cnt = buni%cnt

end subroutine cpl_map_bun

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE:  cpl_map_read0 - Read in a SparseMatrix
!
! !DESCRIPTION: 
!     Read in mapping matrix data from a SCRIP netCDF data file.
!
! !REMARKS:
!   This routine must be generalized to
!   o deal with mapping data between arbitrary domains
!   o read in data files with arbitrary file names
!
! !SEE ALSO: 
!    mct/m_SparseMatrix.F90 (MCT source code)
!
! !REVISION HISTORY: 
!    2006-Nov-03 - R. Jacob - do allocates one at a time to lower memory footprint.
!    2002-Mar-19 - B. Kauffman -- first version
! 
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_read0(sMat,areasrc,areadst,fileName)

! !USES:

#include <netcdf.inc>

! !INPUT/OUTPUT PARAMETERS:

   type(mct_sMat)   ,intent(out) :: sMat   ! mapping data
   type(mct_Avect) ,intent(out):: areasrc ! area of src grid from mapping file
   type(mct_Avect) ,intent(out):: areadst ! area of dst grid from mapping file
   character(*),intent(in)  :: filename  ! netCDF file to read

! !EOP

   !--- local ---
   integer(IN)           :: n       ! generic loop indicies
   integer(IN)           :: na      ! size of source domain
   integer(IN)           :: nb      ! size of destination domain
   integer(IN)           :: ns      ! number of non-zero elements in matrix
   integer(IN)           :: ni,nj   ! number of row and col in the matrix
   integer(IN)           :: igrow   ! aVect index for matrix row
   integer(IN)           :: igcol   ! aVect index for matrix column
   integer(IN)           :: iwgt    ! aVect index for matrix element
   integer(IN)           :: iarea   ! aVect index for area

   real(R8)   ,allocatable :: rtemp(:)  ! reals
   integer(IN),allocatable :: itemp(:)  ! ints

   character,allocatable :: str(:)  ! variable length char string
   character(CL)         :: attstr  ! netCDF attribute name string
   integer(IN)           :: rcode   ! netCDF routine return code
   integer(IN)           :: fid     ! netCDF file      ID
   integer(IN)           :: vid     ! netCDF variable  ID
   integer(IN)           :: did     ! netCDF dimension ID

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_read0) '
   character(*),parameter :: F00 = '("(cpl_map_read0) ",4a)'
   character(*),parameter :: F01 = '("(cpl_map_read0) ",2(a,i7))'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   write(6,F00) "reading mapping matrix data..."

   !----------------------------------------------------------------------------
   ! open & read the file
   !----------------------------------------------------------------------------
   write(6,F00) "* file name                  : ",trim(fileName)
   rcode = nf_open(filename,NF_NOWRITE,fid)
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

   !--- get matrix dimensions ----------
   rcode = nf_inq_dimid (fid, 'n_s', did)  ! size of sparse matrix
   rcode = nf_inq_dimlen(fid, did  , ns)
   rcode = nf_inq_dimid (fid, 'n_a', did)  ! size of  input vector
   rcode = nf_inq_dimlen(fid, did  , na)
   rcode = nf_inq_dimid (fid, 'n_b', did)  ! size of output vector
   rcode = nf_inq_dimlen(fid, did  , nb)

   write(6,F01) "* matrix dimensions rows x cols :",na,' x',nb
   write(6,F01) "* number of non-zero elements: ",ns

   !----------------------------------------------------------------------------
   ! init the mct sMat data type
   !----------------------------------------------------------------------------
   ! mct_sMat_init must be given the number of rows and columns that
   ! would be in the full matrix.  Nrows= size of output vector=nb.
   ! Ncols = size of input vector = na.
   call mct_sMat_init(sMat, nb, na, ns)

   igrow = mct_sMat_indexIA(sMat,'grow')
   igcol = mct_sMat_indexIA(sMat,'gcol')
   iwgt  = mct_sMat_indexRA(sMat,'weight')

   !--- read and load matrix weights ---
   allocate(rtemp(ns),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate weights',rcode)

   rcode = nf_inq_varid     (fid,'S'  ,vid)
   rcode = nf_get_var_double(fid,vid  ,rtemp  )
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

   sMat%data%rAttr(iwgt ,:) =   rtemp(:)

   deallocate(rtemp, stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: deallocate weights',rcode)

   !--- read and load rows ---
   allocate(itemp(ns),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate rows',rcode)

   rcode = nf_inq_varid     (fid,'row',vid)
   rcode = nf_get_var_int   (fid,vid  ,itemp)
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

   sMat%data%iAttr(igrow,:) = itemp(:)

   !--- read and load columns ---
   itemp(:) = 0

   rcode = nf_inq_varid     (fid,'col',vid)
   rcode = nf_get_var_int   (fid,vid  ,itemp)
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

   sMat%data%iAttr(igcol,:) = itemp(:)

   deallocate(itemp, stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: deallocate cols',rcode)

   !--- read and load area_a ---
   allocate(rtemp(na),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate area_a',rcode)

   rcode = nf_inq_varid     (fid,'area_a',vid)
   rcode = nf_get_var_double(fid,vid     ,rtemp)
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

   call mct_aVect_init(areasrc,' ',cpl_map_areaAV_field,na)

   iarea = mct_aVect_indexRA(areasrc,cpl_map_areaAV_field)
   areasrc%rAttr(iarea,1:na) = rtemp(1:na)

   deallocate(rtemp, stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: deallocate area_a',rcode)

   !--- read and load area_b ---
   allocate(rtemp(nb),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate area_b',rcode)

   rcode = nf_inq_varid     (fid,'area_b',vid)
   rcode = nf_get_var_double(fid,vid     ,rtemp)
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
   rcode = nf_close(fid)

   call mct_aVect_init(areadst,' ',cpl_map_areaAV_field,nb)

   iarea = mct_aVect_indexRA(areadst,cpl_map_areaAV_field)
   areadst%rAttr(iarea,1:nb) = rtemp(1:nb)

   deallocate(rtemp, stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: deallocate area_b',rcode)

   write(6,F00) "... done reading file"
   call shr_sys_flush(6)

end subroutine cpl_map_read0

!===============================================================================
! !BOP ===========================================================================
!
! !IROUTINE:  cpl_map_Mattest - consistancy test for a uniprocessor SparseMatrix
!
! !DESCRIPTION: 
!     Perform a variety of tests on a *uniprocessor* mct sparse-matrix data type
!
! !SEE ALSO: 
!    mct/m_SparseMatrix.F90 (MCT source code)
!
! !REVISION HISTORY: 
!    2002-Mar-19 - B. Kauffman -- first version
! 
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_Mattest(sMat)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(mct_sMat),intent(in) :: sMat  ! mapping data

! !EOP

   !--- local ---
   integer(IN)    :: ns           ! size of aVect
   integer(IN)    :: igrow        ! aVect index for matrix row
   integer(IN)    :: igcol        ! aVect index for matrix column
   integer(IN)    :: iwgt         ! aVect index for matrix element

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_Mattest) '
   character(*),parameter :: F00 = '("(cpl_map_Mattest) ",4a)'
   character(*),parameter :: F01 = '("(cpl_map_Mattest) ",2(a,i9))'
   character(*),parameter :: F02 = '("(cpl_map_Mattest) ",2(a,es11.3))'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   write(6,F00) "consistancy check on mapping matrix data..."


   ns    = mct_sMat_lSize  (sMat)
   igRow = mct_sMat_indexIA(sMat,'grow')
   igCol = mct_sMat_indexIA(sMat,'gcol')
   iwgt  = mct_sMat_indexRA(sMat,'weight')

   write(6,F01) '* number of rows x cols = ',sMat%nRows, &
                                       ' x ',sMat%nCols

   write(6,F01) "* number non-zero elements =",ns

   if (ns > 0) then
      write(6,F02) "* min/max   element value = ",                  &
                      minval(sMat%data%rAttr(iwgt ,:)),", ", &
                      maxval(sMat%data%rAttr(iwgt ,:))

      write(6,F01) "* first/last row    value = ",          &
                      sMat%data%iAttr(igRow, 1),", ", &
                      sMat%data%iAttr(igRow,ns)
      write(6,F01) "* min/max    row    value = ",                  &
                      minval(sMat%data%iAttr(igRow,:)),", ", &
                      maxval(sMat%data%iAttr(igRow,:))

      write(6,F01) "* first/last column value = ",          &
                      sMat%data%iAttr(igCol, 1),", ", &
                      sMat%data%iAttr(igCol,ns)
      write(6,F01) "* min/max    column value = ",                  &
                      minval(sMat%data%iAttr(igCol,:)),", ", &
                      maxval(sMat%data%iAttr(igCol,:))
   endif

   call shr_sys_flush(6)

end subroutine cpl_map_Mattest

!===============================================================================
! !BOP ===========================================================================
!
! !IROUTINE:  cpl_map_MatPtest - consistancy test for a SparseMatrixPlus
!
! !DESCRIPTION: 
!     Perform a variety of tests on an mct sparse-matrix plus data type
!
! !SEE ALSO: 
!    mct/m_SparseMatrixPlus.F90 (MCT source code)
!
! !REVISION HISTORY: 
!    2007-Jan-17 - T. Craig -- first version
! 
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_MatPtest(sMatP)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(mct_sMatP),intent(in) :: sMatP  ! mapping data

! !EOP

   type(mct_string) :: mctStr  ! temporary mct/mpeu string var
   character(CL)        :: strategy

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_MatPtest) '
   character(*),parameter :: F00 = '("(cpl_map_MatPtest) ",4a)'
   character(*),parameter :: F01 = '("(cpl_map_MatPtest) ",2(a,i9))'
   character(*),parameter :: F02 = '("(cpl_map_MatPtest) ",2(a,es11.3))'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   write(6,F00) "consistancy check on sMatP mapping matrix data..."

!---this results in internal compiler error on lightning, tcraig, jan, 2007.
!   write(6,F00) "strategy = ",mct_string_tochar(sMatP%Strategy)

   mctStr = sMatP%Strategy
   strategy = mct_string_tochar(mctStr)
   write(6,F00) "strategy = ",trim(strategy)

   write(6,F01) "tag      = ",sMatP%Tag

   call cpl_map_mattest(sMatP%Matrix)

   call shr_sys_flush(6)

end subroutine cpl_map_MatPtest

!===============================================================================
! !BOP ===========================================================================
!
! !IROUTINE: cpl_map_npFixNone - just return, do nothing.
!
! !DESCRIPTION:
!    Stub for north pole correction of vector fields.  This one just returns
!    and does no adjustments.
!
! !REVISION HISTORY:
!    1Feb03 - T. Craig -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_npFixNone(buni,buno,fld1,fld2)

! !USES:

   use cpl_const_mod

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout):: buni    ! input bundle
   type(cpl_bundle),intent(out)  :: buno    ! output bundle
   character(*)    ,intent(in)   :: fld1    ! name of first input field
   character(*)    ,intent(in)   :: fld2    ! name of second input field

! !EOP

   !--- local ---
   logical,save :: first_call = .true. ! flags 1st invocation of routine

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_npFixNone) '
   character(*),parameter :: F00 = "('(cpl_map_npFixNone) ',8a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   if (first_call) then
     write(6,F00) ' WARNING, no npFix adjustments '
     first_call = .false.
   endif

end subroutine cpl_map_npFixNone

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_npFixNew3 - correct the north pole mapping of velocity fields
!
! !DESCRIPTION:
!    Correct the north pole mapping of velocity fields from the atm to ocn
!    grids.  This assumes the input grid is a regular lat/lon with the north
!    pole surrounded by the last latitude line in the input array.  The
!    longitudes in the last latitude must be ordered and equally spaced.
!
! !REVISION HISTORY:
!    29Aug03 - T. Craig -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_npFixNew3(buni,buno,fld1,fld2)

! !USES:

   use cpl_const_mod
   use shr_timer_mod       ! share timer routines

#include <mpif.h>         ! mpi library include file

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout):: buni    ! input bundle
   type(cpl_bundle),intent(out)  :: buno    ! output bundle
   character(*)    ,intent(in)   :: fld1    ! name of first input field
   character(*)    ,intent(in)   :: fld2    ! name of second input field

!EOP
!===============================================================================
!    This version (New3) is the same as New except it saves the gGrid data
!    type from the first call.  This assumes the gGrid used in all calls
!    to npfix is the same.  It is different from New2 in that it doesn't
!    use gather to compute npu and npv.  This is bfb with New2 and New on 
!    9/1/2003.
!===============================================================================

   !--- local ---
   integer(IN)  :: n,m                   ! generic indices
   integer(IN)  :: n1,n2,n3              ! generic indices
   integer(IN)  :: kui,kvi               ! field indices
   integer(IN)  :: kuo,kvo               ! field indices
   integer(IN)  :: kin                   ! index index
   integer(IN)  :: nmin,nmax             ! indices of highest latitude in input
   integer(IN)  :: npts                  ! local number of points in an aV
   integer(IN)  :: num                   ! number of points at highest latitude
   integer(IN)  :: kloni                 ! longitude index on input domain
   integer(IN)  :: klati                 ! latitude index on input domain
   integer(IN)  :: klono                 ! longitude index on output domain
   integer(IN)  :: klato                 ! latitude index on output domain
   integer(IN)  :: index                 ! index value
   integer(IN)  :: nn1,nn2               ! local n to global n values
   real(R8)     :: rindex                ! index value
   real(R8)     :: latmax                ! value of highest latitude
   real(R8)     :: olon,olat             ! output bundle lon/lat
   real(R8)     :: ilon,ilat             ! input bundle lon/lat
   real(R8)     :: npu,npv               ! np velocity fields relative to lon
   real(R8)     :: theta1,theta2         ! angles for trig functions
   real(R8),allocatable :: ilon1(:)      ! lon of input grid at highest latitude
   real(R8),allocatable :: ilat1(:)      ! lat of input grid at highest latitude
   real(R8),allocatable :: ilon2(:)      ! lon of input grid at highest latitude
   real(R8),allocatable :: ilat2(:)      ! lat of input grid at highest latitude
   real(R8)     :: w1,w2,w3,w4           ! weights
   real(R8)     :: f1,f2,f3,f4           ! function values
   real(R8)     :: alpha,beta            ! used to generate weights
   real(R8)     :: rtmp                  ! real temporary
   real(R8),allocatable :: lData(:,:)    ! last lat local input bundle data
                                         ! also compressed global data
   real(R8),allocatable :: gData(:,:,:)  ! last lat gathered input bundle data
   type(mct_aVect),save :: gGrid     ! global/gathered input bundle data
   integer(IN),parameter :: pid0 = 0     ! root process pid = zero
   logical      :: found                 ! search for new interpolation
   integer(IN)  :: rcode                 ! error code
   integer(IN)  :: np1                   ! n+1 or tmp
   real(R8)     :: ilon1x                ! tmp
   logical,save :: first_call = .true.   ! flags 1st invocation of routine
   integer(IN),save :: tnpf1,tnpf2,tnpf3,tnpf4,tnpf5,tnpf6,tnpf7,tnpf8,tnpf9

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_npFixNew3) '
   character(*),parameter :: F00 = "('(cpl_map_npFixNew3) ',8a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   if (first_call) then
     write(6,F00) " compute bilinear weights & indicies for NP region."
     call shr_timer_get(tnpf1,'tnpf1')
     call shr_timer_get(tnpf2,'tnpf2')
     call shr_timer_get(tnpf3,'tnpf3')
     call shr_timer_get(tnpf4,'tnpf4')
     call shr_timer_get(tnpf5,'tnpf5')
     call shr_timer_get(tnpf6,'tnpf6')
     call shr_timer_get(tnpf7,'tnpf7')
     call shr_timer_get(tnpf8,'tnpf8')
     call shr_timer_get(tnpf9,'tnpf9')
     call mct_aVect_gather(buni%dom%lGrid,gGrid,buni%dom%gsMap,pid0,cpl_comm_comp,rcode)
     if (cpl_comm_comp_pid /= pid0) call mct_aVect_clean(gGrid)
     call mct_aVect_bcast(gGrid,pid0,cpl_comm_comp,rcode)
     first_call = .false.
   endif

   call shr_timer_start(tnpf1)

   kui = mct_aVect_indexRA(buni%data,fld1,perrWith=subName)
   kvi = mct_aVect_indexRA(buni%data,fld2,perrWith=subName)
   kuo = mct_aVect_indexRA(buno%data,fld1,perrWith=subName)
   kvo = mct_aVect_indexRA(buno%data,fld2,perrWith=subName)

   kin   = mct_aVect_indexRA(buni%dom%lGrid,"index",perrWith=subName)
   klati = mct_aVect_indexRA(buni%dom%lGrid,"lat",perrWith=subName)
   kloni = mct_aVect_indexRA(buni%dom%lGrid,"lon",perrWith=subName)
   klato = mct_aVect_indexRA(buno%dom%lGrid,"lat",perrWith=subName)
   klono = mct_aVect_indexRA(buno%dom%lGrid,"lon",perrWith=subName)

   nmin = (buni%dom%ni)*(buni%dom%nj-1) + 1
   nmax = buni%dom%n
   num  = buni%dom%ni

   call shr_timer_stop(tnpf1)
   call shr_timer_start(tnpf2)

   !--- barrier not required but interesting for timing. ---
!  call shr_mpi_barrier(cpl_comm_comp,subName//" barrier")

   call shr_timer_stop(tnpf2)

   allocate(lData(3,num))
   allocate(gData(3,num,cpl_comm_comp_npe))
   lData = 0._R8
   gData = 0._R8
   npts = mct_aVect_lSize(buni%data)
   m = 0   
   do n=1,npts
     rindex = buni%dom%lGrid%rAttr(kin,n)
     if (rindex.ge.nmin) then
       m=m+1
       lData(1,m) = rindex
       lData(2,m) = buni%data%rAttr(kui,n)
       lData(3,m) = buni%data%rAttr(kvi,n)
     endif
   enddo

   call MPI_ALLGATHER(lData,3*num,MPI_REAL8, &
     gData,3*num,MPI_REAL8,cpl_comm_comp,rcode)

   if (rcode.ne.0) then
     write(6,*) trim(subName),' rcode error ',rcode
     call shr_sys_abort()
   endif

   call shr_timer_start(tnpf3)

   m = 0
   lData = 0._R8
   do n2=1,num
   do n3=1,cpl_comm_comp_npe
     if (gData(1,n2,n3).gt.0.1_R8) then
       m = m+1
       index = nint(gData(1,n2,n3)) - nmin + 1
       lData(1:3,index) = gData(1:3,n2,n3)
     endif
   enddo
   enddo
   if (m.ne.num) write(6,*) trim(subName),' error allgather ',m,num
   do n2=1,num
     if (lData(1,n2).lt.0.1_R8) then
       write(6,*) trim(subName),' error allgather2 ',n2
     endif
   enddo

   call shr_timer_stop(tnpf3)
   call shr_timer_start(tnpf4)

   allocate(ilon1(num))
   allocate(ilon2(num))
   allocate(ilat1(num))
   allocate(ilat2(num))

   call shr_timer_stop(tnpf4)
   call shr_timer_start(tnpf5)

   latmax = gGrid%rAttr(klati,nmin)
   npu = 0._R8
   npv = 0._R8
   do n = 1,num
     np1 = mod(n,num)+1
     nn1 = nmin + n - 1
     nn2 = nmin + np1 - 1 
     rtmp = gGrid%rAttr(kloni,nn1)
     ilon1(n) = mod(rtmp+360._R8,360._R8)
     rtmp = gGrid%rAttr(kloni,nn2)
     ilon2(n) = mod(rtmp+360._R8,360._R8)
     ilat1(n) =     gGrid%rAttr(klati,nn1)
     ilat2(n) =     gGrid%rAttr(klati,nn2)
     if (ilon2(n) < ilon1(n)) ilon2(n) = ilon2(n) + 360._R8

     latmax = max(latmax,ilat1(n))

     theta1 = ilon1(n)*cpl_const_deg2rad
     npu = npu + cos(theta1)*lData(2,n) &
               - sin(theta1)*lData(3,n)
     npv = npv + sin(theta1)*lData(2,n) &
               + cos(theta1)*lData(3,n)
   enddo
   npu = npu / real(num,R8)
   npv = npv / real(num,R8)

   call shr_timer_stop(tnpf5)
   call shr_timer_start(tnpf6)

   npts = mct_aVect_lSize(buno%data)
   do m = 1,npts
     olat = buno%dom%lGrid%rAttr(klato,m)
     if (olat >= latmax) then
       rtmp = buno%dom%lGrid%rAttr(klono,m)
       olon = mod(rtmp,360._R8)
       n = 1
       found = .false.
       do while (n <= num .and. .not.found )
         if (    olon >= ilon1(n) .and. olon < ilon2(n) .or.   &
            olon+360._R8 >= ilon1(n) .and. olon < ilon2(n)) then
           np1 = mod(n,num)+1
           ilat = (ilat1(n) + ilat2(n)) * 0.5_R8
           if (ilon2(n) == ilon1(n)) then
             alpha = 0.5_R8
           else if (    olon >= ilon1(n) .and. olon < ilon2(n)) then
             alpha = (olon - ilon1(n)) / (ilon2(n) - ilon1(n))
           else if (olon+360._R8>= ilon1(n) .and. olon < ilon2(n)) then
             alpha = (olon+360._R8 - ilon1(n)) / (ilon2(n) - ilon1(n))
           else
             write(6,*) subName,' ERROR: olon ',olon,ilon1(n),ilon2(n)
           endif
           if (ilat >= 90._R8) then
             beta  = 1.0_R8
           else
             beta  = (olat - ilat) / (90._R8 - ilat)
           endif
           w1 = (1.0_R8-alpha)*(1.0_R8-beta)
           w2 = (    alpha)*(1.0_R8-beta)
           w3 = (    alpha)*(    beta)
           w4 = (1.0_R8-alpha)*(    beta)

           theta1 = ilon1(n)*cpl_const_deg2rad
           theta2 = ilon2(n)*cpl_const_deg2rad

           f1 = lData(2,n)
           f2 = lData(2,np1)
           f3 =  cos(theta1)*npu + sin(theta1)*npv
           f4 =  cos(theta2)*npu + sin(theta2)*npv
           rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
           buno%data%rAttr(kuo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4

           f1 = lData(3,n)
           f2 = lData(3,np1)
           f3 = -sin(theta1)*npu + cos(theta1)*npv
           f4 = -sin(theta2)*npu + cos(theta2)*npv
           rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
           buno%data%rAttr(kvo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4
           found = .true.
         endif
         n = n + 1     ! normal increment
       enddo
       if ( .not.found ) then
         write(6,*) subName,' ERROR: found = false ',found,m,olon,olat
       endif
     endif
   end do

   call shr_timer_stop(tnpf6)
   call shr_timer_start(tnpf7)

   deallocate(gData)
   deallocate(lData)
   deallocate(ilon1)
   deallocate(ilon2)
   deallocate(ilat1)
   deallocate(ilat2)

   call shr_timer_stop(tnpf7)

end subroutine cpl_map_npFixNew3

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_npFixNew3R - correct the north pole mapping of velocity fields
!
! !DESCRIPTION:
!    Correct the north pole mapping of velocity fields from the atm to ocn
!    grids.  This assumes the input grid is a regular lat/lon with the north
!    pole surrounded by the last latitude line in the input array.  The
!    longitudes in the last latitude must be ordered and equally spaced.
!
!    npFixNew3 and npFixNew3R are the same except New3R saves some information
!    between calls. This means it has to be always called for the same grids
!    (not fields).  New3 is totally general.  New3R improves performance
!    on the Cray X1E by reducing string work.  The R is "reuse"

! !REVISION HISTORY:
!    15Apr06 - T. Craig -- modified New3 for X1E performance, reuse
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_npFixNew3R(buni,buno,fld1,fld2)

! !USES:

   use cpl_const_mod
   use shr_timer_mod       ! share timer routines

#include <mpif.h>         ! mpi library include file

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout):: buni    ! input bundle
   type(cpl_bundle),intent(out)  :: buno    ! output bundle
   character(*)    ,intent(in)   :: fld1    ! name of first input field
   character(*)    ,intent(in)   :: fld2    ! name of second input field

!EOP
!===============================================================================
!    This version (New3R) is the same as New except it saves the gGrid data
!    type from the first call.  This assumes the gGrid used in all calls
!    to npfix is the same.  It is different from New2 in that it doesn't
!    use gather to compute npu and npv.  This is bfb with New2 and New on 
!    9/1/2003.
!===============================================================================

   !--- local ---
   integer(IN)  :: n,m                   ! generic indices
   integer(IN)  :: n1,n2,n3              ! generic indices
   integer(IN)  :: kui,kvi               ! field indices
   integer(IN)  :: kuo,kvo               ! field indices
   integer(IN)  :: kin                   ! index index
   integer(IN)  :: nmin,nmax             ! indices of highest latitude in input
   integer(IN)  :: npts                  ! local number of points in an aV
   integer(IN)  :: num                   ! number of points at highest latitude
   integer(IN)  :: kloni                 ! longitude index on input domain
   integer(IN)  :: klati                 ! latitude index on input domain
   integer(IN)  :: klono                 ! longitude index on output domain
   integer(IN)  :: klato                 ! latitude index on output domain
   integer(IN)  :: index                 ! index value
   integer(IN)  :: nn1,nn2               ! local n to global n values
   real(R8)     :: rindex                ! index value
   real(R8)     :: latmax                ! value of highest latitude
   real(R8)     :: olon,olat             ! output bundle lon/lat
   real(R8)     :: ilon,ilat             ! input bundle lon/lat
   real(R8)     :: npu,npv               ! np velocity fields relative to lon
   real(R8)     :: theta1,theta2         ! angles for trig functions
   real(R8),allocatable :: ilon1(:)      ! lon of input grid at highest latitude
   real(R8),allocatable :: ilat1(:)      ! lat of input grid at highest latitude
   real(R8),allocatable :: ilon2(:)      ! lon of input grid at highest latitude
   real(R8),allocatable :: ilat2(:)      ! lat of input grid at highest latitude
   real(R8)     :: w1,w2,w3,w4           ! weights
   real(R8)     :: f1,f2,f3,f4           ! function values
   real(R8)     :: alpha,beta            ! used to generate weights
   real(R8)     :: rtmp                  ! real temporary
   real(R8),allocatable :: lData(:,:)    ! last lat local input bundle data
                                         ! also compressed global data
   real(R8),allocatable :: gData(:,:,:)  ! last lat gathered input bundle data
   type(mct_aVect),save :: gGrid     ! global/gathered input bundle data
   real(R8),allocatable,save :: alphafound(:) ! list of found alphas
   real(R8),allocatable,save :: betafound(:)  ! list of found betas
   integer(IN),allocatable,save :: mfound(:)  ! list of found ms
   integer(IN),allocatable,save :: nfound(:)  ! list of found ns
   integer(IN),save  :: cntfound              ! number found
   integer(IN)  :: cnt                   ! loop counter
   integer(IN),parameter :: pid0 = 0     ! root process pid = zero
   logical      :: found                 ! search for new interpolation
   integer(IN)  :: rcode                 ! error code
   integer(IN)  :: np1                   ! n+1 or tmp
   real(R8)     :: ilon1x                ! tmp
   logical,save :: first_call = .true.   ! flags 1st invocation of routine
   integer(IN),save :: tnpf1,tnpf2,tnpf3,tnpf4,tnpf5,tnpf6,tnpf7,tnpf8,tnpf9

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_npFixNew3R) '
   character(*),parameter :: F00 = "('(cpl_map_npFixNew3R) ',8a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   if (first_call) then
     write(6,F00) " compute bilinear weights & indicies for NP region."
     call shr_timer_get(tnpf1,'tnpf1')
     call shr_timer_get(tnpf2,'tnpf2')
     call shr_timer_get(tnpf3,'tnpf3')
     call shr_timer_get(tnpf4,'tnpf4')
     call shr_timer_get(tnpf5,'tnpf5')
     call shr_timer_get(tnpf6,'tnpf6')
     call shr_timer_get(tnpf7,'tnpf7')
     call shr_timer_get(tnpf8,'tnpf8')
     call shr_timer_get(tnpf9,'tnpf9')
     call shr_timer_start(tnpf9)
     call mct_aVect_gather(buni%dom%lGrid,gGrid,buni%dom%gsMap,pid0,cpl_comm_comp,rcode)
     if (cpl_comm_comp_pid /= pid0) call mct_aVect_clean(gGrid)
     call mct_aVect_bcast(gGrid,pid0,cpl_comm_comp,rcode)
     call shr_timer_stop(tnpf9)
   endif

   call shr_timer_start(tnpf1)

   kui = mct_aVect_indexRA(buni%data,fld1,perrWith=subName)
   kvi = mct_aVect_indexRA(buni%data,fld2,perrWith=subName)
   kuo = mct_aVect_indexRA(buno%data,fld1,perrWith=subName)
   kvo = mct_aVect_indexRA(buno%data,fld2,perrWith=subName)

   kin   = mct_aVect_indexRA(buni%dom%lGrid,"index",perrWith=subName)
   klati = mct_aVect_indexRA(buni%dom%lGrid,"lat",perrWith=subName)
   kloni = mct_aVect_indexRA(buni%dom%lGrid,"lon",perrWith=subName)
   klato = mct_aVect_indexRA(buno%dom%lGrid,"lat",perrWith=subName)
   klono = mct_aVect_indexRA(buno%dom%lGrid,"lon",perrWith=subName)

   nmin = (buni%dom%ni)*(buni%dom%nj-1) + 1
   nmax = buni%dom%n
   num  = buni%dom%ni

   call shr_timer_stop(tnpf1)
   call shr_timer_start(tnpf2)

   !--- barrier not required but interesting for timing. ---
!  call shr_mpi_barrier(cpl_comm_comp,subName//" barrier")

   allocate(lData(3,num))
   allocate(gData(3,num,cpl_comm_comp_npe))
   lData = 0._R8
   gData = 0._R8
   npts = mct_aVect_lSize(buni%data)
   m = 0   
   do n=1,npts
     rindex = buni%dom%lGrid%rAttr(kin,n)
     if (rindex.ge.nmin) then
       m=m+1
       lData(1,m) = rindex
       lData(2,m) = buni%data%rAttr(kui,n)
       lData(3,m) = buni%data%rAttr(kvi,n)
     endif
   enddo

   call MPI_ALLGATHER(lData,3*num,MPI_REAL8, &
     gData,3*num,MPI_REAL8,cpl_comm_comp,rcode)

   if (rcode.ne.0) then
     write(6,*) trim(subName),' rcode error ',rcode
     call shr_sys_abort()
   endif

   call shr_timer_stop(tnpf2)
   call shr_timer_start(tnpf3)

   m = 0
   lData = 0._R8
   do n2=1,num
   do n3=1,cpl_comm_comp_npe
     if (gData(1,n2,n3).gt.0.1_R8) then
       m = m+1
       index = nint(gData(1,n2,n3)) - nmin + 1
       lData(1:3,index) = gData(1:3,n2,n3)
     endif
   enddo
   enddo
   if (m.ne.num) write(6,*) trim(subName),' error allgather ',m,num
   do n2=1,num
     if (lData(1,n2).lt.0.1_R8) then
       write(6,*) trim(subName),' error allgather2 ',n2
     endif
   enddo

   call shr_timer_stop(tnpf3)
   call shr_timer_start(tnpf4)

   allocate(ilon1(num))
   allocate(ilon2(num))
   allocate(ilat1(num))
   allocate(ilat2(num))

   call shr_timer_stop(tnpf4)
   call shr_timer_start(tnpf5)

   latmax = gGrid%rAttr(klati,nmin)
   npu = 0._R8
   npv = 0._R8
   do n = 1,num
     np1 = mod(n,num)+1
     nn1 = nmin + n - 1
     nn2 = nmin + np1 - 1 
     rtmp = gGrid%rAttr(kloni,nn1)
     ilon1(n) = mod(rtmp+360._R8,360._R8)
     rtmp = gGrid%rAttr(kloni,nn2)
     ilon2(n) = mod(rtmp+360._R8,360._R8)
     ilat1(n) =     gGrid%rAttr(klati,nn1)
     ilat2(n) =     gGrid%rAttr(klati,nn2)
     if (ilon2(n) < ilon1(n)) ilon2(n) = ilon2(n) + 360._R8

     latmax = max(latmax,ilat1(n))

     theta1 = ilon1(n)*cpl_const_deg2rad
     npu = npu + cos(theta1)*lData(2,n) &
               - sin(theta1)*lData(3,n)
     npv = npv + sin(theta1)*lData(2,n) &
               + cos(theta1)*lData(3,n)
   enddo
   npu = npu / real(num,R8)
   npv = npv / real(num,R8)

   call shr_timer_stop(tnpf5)
   call shr_timer_start(tnpf6)

   if (first_call) then
   allocate(mfound(npts),nfound(npts),alphafound(npts),betafound(npts))
   cntfound = 0
   npts = mct_aVect_lSize(buno%data)
   do m = 1,npts
     olat = buno%dom%lGrid%rAttr(klato,m)
     if (olat >= latmax) then
       rtmp = buno%dom%lGrid%rAttr(klono,m)
       olon = mod(rtmp,360._R8)
       n = 1
       found = .false.
       do while (n <= num .and. .not.found )
         if (    olon >= ilon1(n) .and. olon < ilon2(n) .or.   &
            olon+360._R8 >= ilon1(n) .and. olon < ilon2(n)) then
           ilat = (ilat1(n) + ilat2(n)) * 0.5_R8
           if (ilon2(n) == ilon1(n)) then
             alpha = 0.5_R8
           else if (    olon >= ilon1(n) .and. olon < ilon2(n)) then
             alpha = (olon - ilon1(n)) / (ilon2(n) - ilon1(n))
           else if (olon+360._R8>= ilon1(n) .and. olon < ilon2(n)) then
             alpha = (olon+360._R8 - ilon1(n)) / (ilon2(n) - ilon1(n))
           else
             write(6,*) subName,' ERROR: olon ',olon,ilon1(n),ilon2(n)
           endif
           if (ilat >= 90._R8) then
             beta  = 1.0_R8
           else
             beta  = (olat - ilat) / (90._R8 - ilat)
           endif

           cntfound = cntfound + 1
           mfound(cntfound) = m
           nfound(cntfound) = n
           alphafound(cntfound) = alpha
           betafound(cntfound) = beta
           found = .true.

         endif
         n = n + 1     ! normal increment
       enddo
       if ( .not.found ) then
         write(6,*) subName,' ERROR: found = false ',found,m,olon,olat
       endif
     endif
   end do
   endif

   call shr_timer_stop(tnpf6)
   call shr_timer_start(tnpf8)

!DIR$ CONCURRENT
   do cnt = 1,cntfound
           m = mfound(cnt)
           n = nfound(cnt)
           np1 = mod(n,num)+1
           alpha = alphafound(cnt)
           beta = betafound(cnt)

           w1 = (1.0_R8-alpha)*(1.0_R8-beta)
           w2 = (    alpha)*(1.0_R8-beta)
           w3 = (    alpha)*(    beta)
           w4 = (1.0_R8-alpha)*(    beta)

           theta1 = ilon1(n)*cpl_const_deg2rad
           theta2 = ilon2(n)*cpl_const_deg2rad

           f1 = lData(2,n)
           f2 = lData(2,np1)
           f3 =  cos(theta1)*npu + sin(theta1)*npv
           f4 =  cos(theta2)*npu + sin(theta2)*npv
           rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
           buno%data%rAttr(kuo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4

           f1 = lData(3,n)
           f2 = lData(3,np1)
           f3 = -sin(theta1)*npu + cos(theta1)*npv
           f4 = -sin(theta2)*npu + cos(theta2)*npv
           rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
           buno%data%rAttr(kvo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4
   enddo

   call shr_timer_stop(tnpf8)

   call shr_timer_start(tnpf7)

   deallocate(gData)
   deallocate(lData)
   deallocate(ilon1)
   deallocate(ilon2)
   deallocate(ilat1)
   deallocate(ilat2)
   first_call = .false.

   call shr_timer_stop(tnpf7)

end subroutine cpl_map_npFixNew3R

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_npFixNew4R - correct the north pole mapping of velocity fields
!
! !DESCRIPTION:
!    Correct the north pole mapping of velocity fields from the atm to ocn
!    grids.  This assumes the input grid is a regular lat/lon with the north
!    pole surrounded by the last latitude line in the input array.  The
!    longitudes in the last latitude must be ordered and equally spaced.
!
!    4R is a low memory version of 3R.

! !REVISION HISTORY:
!    12Feb07 - T. Craig -- modified New3R to reduce memory
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_npFixNew4R(buni,buno,fld1,fld2)

! !USES:

   use cpl_const_mod
   use shr_timer_mod       ! share timer routines

#include <mpif.h>         ! mpi library include file

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout):: buni    ! input bundle
   type(cpl_bundle),intent(out)  :: buno    ! output bundle
   character(*)    ,intent(in)   :: fld1    ! name of first input field
   character(*)    ,intent(in)   :: fld2    ! name of second input field

!EOP
!===============================================================================
!    This version (New4R) is the same as 3R except it uses a lot less memory
!    and is a bit faster.  Like 3R, it saves data between calls and so
!    assumes the input grid remains constant for all calls.  This is bfb
!    with 3R on bluevista as of 2/15/07.
!===============================================================================

   !--- local ---
   integer(IN)  :: n,m                   ! generic indices
   integer(IN)  :: n1,n2,n3              ! generic indices
   integer(IN)  :: kui,kvi               ! field indices
   integer(IN)  :: kuo,kvo               ! field indices
   integer(IN)  :: kin                   ! index index
   integer(IN)  :: nmin,nmax             ! indices of highest latitude in input
   integer(IN)  :: npts                  ! local number of points in an aV
   integer(IN)  :: num                   ! number of points at highest latitude
   integer(IN)  :: kloni                 ! longitude index on input domain
   integer(IN)  :: klati                 ! latitude index on input domain
   integer(IN)  :: klono                 ! longitude index on output domain
   integer(IN)  :: klato                 ! latitude index on output domain
   integer(IN)  :: index                 ! index value
   real(R8)     :: rindex                ! index value
   real(R8)     :: latmax                ! value of highest latitude
   real(R8)     :: olon,olat             ! output bundle lon/lat
   real(R8)     :: ilon,ilat             ! input bundle lon/lat
   real(R8)     :: npu,npv               ! np velocity fields relative to lon
   real(R8)     :: theta1,theta2         ! angles for trig functions
   real(R8),allocatable,save :: ilon1(:)      ! lon of input grid at highest latitude
   real(R8),allocatable,save :: ilat1(:)      ! lat of input grid at highest latitude
   real(R8),allocatable,save :: ilon2(:)      ! lon of input grid at highest latitude
   real(R8),allocatable,save :: ilat2(:)      ! lat of input grid at highest latitude
   real(R8),allocatable  :: rarray(:)    ! temporary array
   real(R8),allocatable  :: rarray2(:,:)    ! temporary array
   real(R8)     :: w1,w2,w3,w4           ! weights
   real(R8)     :: f1,f2,f3,f4           ! function values
   real(R8)     :: alpha,beta            ! used to generate weights
   real(R8)     :: rtmp                  ! real temporary
   real(R8),allocatable :: lData(:,:)    ! last lat local input bundle data
                                         ! also compressed global data
   real(R8),allocatable,save :: alphafound(:) ! list of found alphas
   real(R8),allocatable,save :: betafound(:)  ! list of found betas
   integer(IN),allocatable,save :: mfound(:)  ! list of found ms
   integer(IN),allocatable,save :: nfound(:)  ! list of found ns
   real(R8),allocatable      :: rfound(:)   ! temporary for copy
   integer(IN),allocatable   :: ifound(:)   ! temporary for copy
   integer(IN),save  :: cntfound              ! number found
   integer(IN)  :: cnt                   ! loop counter
   integer(IN),parameter :: pid0 = 0     ! root process pid = zero
   logical      :: found                 ! search for new interpolation
   integer(IN)  :: rcode                 ! error code
   integer(IN)  :: np1                   ! n+1 or tmp
   real(R8)     :: ilon1x                ! tmp
   logical,save :: first_call = .true.   ! flags 1st invocation of routine
   integer(IN),save :: tnpf1,tnpf2,tnpf3,tnpf4,tnpf5,tnpf6,tnpf7,tnpf8,tnpf9

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_npFixNew4R) '
   character(*),parameter :: F00 = "('(cpl_map_npFixNew4R) ',8a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   if (first_call) then
     write(6,F00) " compute bilinear weights & indicies for NP region."
     call shr_timer_get(tnpf1,'tnpf1')
     call shr_timer_get(tnpf2,'tnpf2')
     call shr_timer_get(tnpf3,'tnpf3')
     call shr_timer_get(tnpf4,'tnpf4')
     call shr_timer_get(tnpf5,'tnpf5')
     call shr_timer_get(tnpf6,'tnpf6')
     call shr_timer_get(tnpf7,'tnpf7')
     call shr_timer_get(tnpf8,'tnpf8')
     call shr_timer_get(tnpf9,'tnpf9')
   endif

   call shr_timer_start(tnpf1)

   kui = mct_aVect_indexRA(buni%data,fld1,perrWith=subName)
   kvi = mct_aVect_indexRA(buni%data,fld2,perrWith=subName)
   kuo = mct_aVect_indexRA(buno%data,fld1,perrWith=subName)
   kvo = mct_aVect_indexRA(buno%data,fld2,perrWith=subName)

   kin   = mct_aVect_indexRA(buni%dom%lGrid,"index",perrWith=subName)
   klati = mct_aVect_indexRA(buni%dom%lGrid,"lat",perrWith=subName)
   kloni = mct_aVect_indexRA(buni%dom%lGrid,"lon",perrWith=subName)
   klato = mct_aVect_indexRA(buno%dom%lGrid,"lat",perrWith=subName)
   klono = mct_aVect_indexRA(buno%dom%lGrid,"lon",perrWith=subName)

   nmin = (buni%dom%ni)*(buni%dom%nj-1) + 1
   nmax = buni%dom%n
   num  = buni%dom%ni

   call shr_timer_stop(tnpf1)

   if (first_call) then
     call shr_timer_start(tnpf4)
     allocate(ilon1(num))
     allocate(ilon2(num))
     allocate(ilat1(num))
     allocate(ilat2(num))

     ilon1 = 0._r8
     ilon2 = 0._r8
     ilat1 = 0._r8
     ilat2 = 0._r8

     npts = mct_aVect_lSize(buni%dom%lGrid)
     do m=1,npts
       rindex = buni%dom%lGrid%rAttr(kin,m)
       if (rindex.ge.nmin) then
         n = nint(rindex) - nmin + 1
         rtmp = buni%dom%lGrid%rAttr(kloni,m)
         ilon1(n) = mod(rtmp+360._R8,360._R8)
         ilat1(n) = buni%dom%lGrid%rAttr(klati,m)
       endif
     enddo

     do n = 1,num
       np1 = mod(n,num)+1
       ilat2(n) = ilat1(np1)
       ilon2(n) = ilon1(np1)
       if (ilon2(n) < ilon1(n)) ilon2(n) = ilon2(n) + 360._R8
     enddo

     !--- all gather local data, MPI_SUM is low memory and simple
     !--- but is a performance penalty compared to gatherv and copy
     !--- or a fancy send/recv 

     allocate(rarray(num))
     rarray = ilat1
     call MPI_ALLREDUCE(rarray,ilat1,num,MPI_REAL8,MPI_SUM,cpl_comm_comp,rcode)
     if (rcode.ne.0) then
       write(6,*) trim(subName),' ilat1 rcode error ',rcode
       call shr_sys_abort()
     endif
     rarray = ilon1
     call MPI_ALLREDUCE(rarray,ilon1,num,MPI_REAL8,MPI_SUM,cpl_comm_comp,rcode)
     if (rcode.ne.0) then
       write(6,*) trim(subName),' ilon1 rcode error ',rcode
       call shr_sys_abort()
     endif
     rarray = ilat2
     call MPI_ALLREDUCE(rarray,ilat2,num,MPI_REAL8,MPI_SUM,cpl_comm_comp,rcode)
     if (rcode.ne.0) then
       write(6,*) trim(subName),' ilat2 rcode error ',rcode
       call shr_sys_abort()
     endif
     rarray = ilon2
     call MPI_ALLREDUCE(rarray,ilon2,num,MPI_REAL8,MPI_SUM,cpl_comm_comp,rcode)
     if (rcode.ne.0) then
       write(6,*) trim(subName),' ilon2 rcode error ',rcode
       call shr_sys_abort()
     endif
     deallocate(rarray)

     latmax = maxval(ilat1)

     call shr_timer_stop(tnpf4)

     !--- compute weights and save them ---

     call shr_timer_start(tnpf6)

     npts = mct_aVect_lSize(buno%data)
     allocate(mfound(npts),nfound(npts),alphafound(npts),betafound(npts))
     cntfound = 0
     do m = 1,npts
       olat = buno%dom%lGrid%rAttr(klato,m)
       if (olat >= latmax) then
         rtmp = buno%dom%lGrid%rAttr(klono,m)
         olon = mod(rtmp,360._R8)
         n = 1
         found = .false.
         do while (n <= num .and. .not.found )
           if (    olon >= ilon1(n) .and. olon < ilon2(n) .or.   &
              olon+360._R8 >= ilon1(n) .and. olon < ilon2(n)) then
             ilat = (ilat1(n) + ilat2(n)) * 0.5_R8
             if (ilon2(n) == ilon1(n)) then
               alpha = 0.5_R8
             else if (    olon >= ilon1(n) .and. olon < ilon2(n)) then
               alpha = (olon - ilon1(n)) / (ilon2(n) - ilon1(n))
             else if (olon+360._R8>= ilon1(n) .and. olon < ilon2(n)) then
               alpha = (olon+360._R8 - ilon1(n)) / (ilon2(n) - ilon1(n))
             else
               write(6,*) subName,' ERROR: olon ',olon,ilon1(n),ilon2(n)
             endif
             if (ilat >= 90._R8) then
               beta  = 1.0_R8
             else
               beta  = (olat - ilat) / (90._R8 - ilat)
             endif

             cntfound = cntfound + 1
             mfound(cntfound) = m
             nfound(cntfound) = n
             alphafound(cntfound) = alpha
             betafound(cntfound) = beta
             found = .true.

           endif
           n = n + 1     ! normal increment
         enddo
         if ( .not.found ) then
           write(6,*) subName,' ERROR: found = false ',found,m,olon,olat
         endif
       endif
     end do

     allocate(ifound(npts))
     ifound = mfound
     deallocate(mfound)
     if (cntfound > 0) then
        allocate(mfound(cntfound))
        mfound(1:cntfound) = ifound(1:cntfound)
     endif

     ifound = nfound
     deallocate(nfound)
     if (cntfound > 0) then
        allocate(nfound(cntfound))
        nfound(1:cntfound) = ifound(1:cntfound)
     endif
     deallocate(ifound)

     allocate(rfound(npts))
     rfound = alphafound
     deallocate(alphafound)
     if (cntfound > 0) then
        allocate(alphafound(cntfound))
        alphafound(1:cntfound) = rfound(1:cntfound)
     endif

     rfound = betafound
     deallocate(betafound)
     if (cntfound > 0) then
        allocate(betafound(cntfound))
        betafound(1:cntfound) = rfound(1:cntfound)
     endif
     deallocate(rfound)
   
     call shr_timer_stop(tnpf6)

   endif

   call shr_timer_start(tnpf2)

   !--- barrier not required but interesting for timing. ---
!  call shr_mpi_barrier(cpl_comm_comp,subName//" barrier")

   !--- extract index, u, v from buni ---

   allocate(lData(3,num))
   lData = 0._R8
   npts = mct_aVect_lSize(buni%data)
   do n=1,npts
     rindex = buni%dom%lGrid%rAttr(kin,n)
     if (rindex.ge.nmin) then
       m = nint(rindex) - nmin + 1
       lData(1,m) = rindex
       lData(2,m) = buni%data%rAttr(kui,n)
       lData(3,m) = buni%data%rAttr(kvi,n)
     endif
   enddo

   !--- all gather local data, MPI_SUM is low memory and simple
   !--- but is a performance penalty compared to gatherv and copy

   allocate(rarray2(3,num))
   rarray2=lData
   call MPI_ALLREDUCE(rarray2,lData,3*num,MPI_REAL8,MPI_SUM,cpl_comm_comp,rcode)
   deallocate(rarray2)

   if (rcode.ne.0) then
     write(6,*) trim(subName),' rcode error ',rcode
     call shr_sys_abort()
   endif

   do n2=1,num
     if (lData(1,n2).lt.0.1_R8) then
       write(6,*) trim(subName),' error allreduce ',n2
     endif
   enddo

   call shr_timer_stop(tnpf2)

!   call shr_timer_start(tnpf3)
!   call shr_timer_stop(tnpf3)

   !--- compute npu, npv (pole data) and initialize ilon,ilat arrays ---

   call shr_timer_start(tnpf5)

   npu = 0._R8
   npv = 0._R8
   do n = 1,num
     theta1 = ilon1(n)*cpl_const_deg2rad
     npu = npu + cos(theta1)*lData(2,n) &
               - sin(theta1)*lData(3,n)
     npv = npv + sin(theta1)*lData(2,n) &
               + cos(theta1)*lData(3,n)
   enddo
   npu = npu / real(num,R8)
   npv = npv / real(num,R8)

   call shr_timer_stop(tnpf5)

   !--- compute updated pole vectors ---

   call shr_timer_start(tnpf8)

!DIR$ CONCURRENT
   do cnt = 1,cntfound
      m = mfound(cnt)
      n = nfound(cnt)
      np1 = mod(n,num)+1
      alpha = alphafound(cnt)
      beta = betafound(cnt)

      w1 = (1.0_R8-alpha)*(1.0_R8-beta)
      w2 = (    alpha)*(1.0_R8-beta)
      w3 = (    alpha)*(    beta)
      w4 = (1.0_R8-alpha)*(    beta)

      theta1 = ilon1(n)*cpl_const_deg2rad
      theta2 = ilon2(n)*cpl_const_deg2rad

      f1 = lData(2,n)
      f2 = lData(2,np1)
      f3 =  cos(theta1)*npu + sin(theta1)*npv
      f4 =  cos(theta2)*npu + sin(theta2)*npv
      rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
      buno%data%rAttr(kuo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4

      f1 = lData(3,n)
      f2 = lData(3,np1)
      f3 = -sin(theta1)*npu + cos(theta1)*npv
      f4 = -sin(theta2)*npu + cos(theta2)*npv
      rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
      buno%data%rAttr(kvo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4
   enddo

   call shr_timer_stop(tnpf8)

   call shr_timer_start(tnpf7)

   deallocate(lData)
   first_call = .false.

   call shr_timer_stop(tnpf7)

end subroutine cpl_map_npFixNew4R

!===============================================================================
! !BOP ===========================================================================
!
! !IROUTINE: cpl_map_npFixOld - correct the north pole mapping of velocity fields
!
! !DESCRIPTION:
!    Correct the north pole mapping of velocity fields from the atm to ocn
!    grids.  This assumes the input grid is a regular lat/lon with the north
!    pole surrounded by the last latitude line in the input array.  The
!    longitudes in the last latitude must be ordered.
!
! !REVISION HISTORY:
!    20Sep02 - T. Craig -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_npFixOld(buni,buno,fld1,fld2)

! !USES:

   use cpl_const_mod

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout):: buni    ! input bundle
   type(cpl_bundle),intent(inout):: buno    ! output bundle
   character(*)    ,intent(in)   :: fld1    ! name of first input field
   character(*)    ,intent(in)   :: fld2    ! name of second input field

! !EOP

   !----- local -----
   integer :: n           ! index relative to output array (subroutine arg)
   integer :: npts        ! number of output array points that need fixing
   integer :: nn          ! index that goes from 1 to npts
   integer :: m           ! index relative to input array (subroutine arg)
   integer :: mpts        ! number of grid points in top row of input grid
   integer :: mm          ! index that goes from 1 to mpts
   integer :: m1st        ! 1st m st y(m)=max input latitude

   real(R8)    :: npu         ! north-pole u-velocity
   real(R8)    :: npv         ! north-pole v-velocity

   real(R8)   ,allocatable ::  wght(:,:) ! 4 bilin map weights
   integer    ,allocatable ::    mw(:,:) ! index on input grid
   integer    ,allocatable ::    nw  (:) ! index on output grid
   real(R8)   ,allocatable :: Sa_in_u(:) ! local NP u-velocity "on" input grid
   real(R8)   ,allocatable :: Sa_in_v(:) ! local NP v-velocity "on" input grid
   real(R8)   ,allocatable ::  theta (:) ! rotation angle wrt 0 deg east (radians)

   real(R8)    :: alpha,beta    ! used to compute bilin weights
   real(R8)    :: dx,dy         ! used to compute biline weights
   real(R8)    :: xm,ym         ! x/y-coords of output data for bilin weights
   real(R8)    :: xl,xh,yl,yh   ! x/y-coords of  input data for bilin weights
   real(R8)    :: w1,w2,w3,w4   ! 4 bilin weights
   integer     :: m1,m2,m3,m4   ! 4 bilin indicies wrt input grid
   integer     :: mm3,mm4       ! bilin indicies wrt fabricated NP input grid row
   real(R8)    :: f1,f2,f3,f4   ! 4 bilin function values on input grid

   logical :: first_call = .true. ! flags 1st invocation of routine
   integer :: timer_n             ! system-clock-timer number

   type(mct_aVect) :: gDatai      ! global/gathered input bundle data
   type(mct_aVect) :: gGridi      ! global/gathered input bundle data
   type(mct_aVect) :: gDatao      ! global/gathered output bundle data
   type(mct_aVect) :: gGrido      ! global/gathered output bundle data

   integer(IN),parameter :: pid0 = 0     ! root process pid = zero
   integer(IN)  :: n_a,n_o               ! a/o input/output sizes
   integer(IN)  :: kui,kvi               ! field indices
   integer(IN)  :: kuo,kvo               ! field indices
   integer(IN)  :: kloni                 ! longitude index on input domain
   integer(IN)  :: klati                 ! latitude index on input domain
   integer(IN)  :: klono                 ! longitude index on output domain
   integer(IN)  :: klato                 ! latitude index on output domain
   integer(IN)  :: rcode                 ! error code

   save  

   !----- formats -----
   character(*),parameter :: subName = '(cpl_map_npFixOld) '
   character(*),parameter :: F00 = "('(cpl_map_npFixOld) ',8a)"

!-------------------------------------------------------------------------------
! PURPOSE:
!    Compute/correct atm wind velocties on ocn (or ice) grid.
!
! ASSUMPTIONS:
!  o the input (atm) grid has a "top row" with a constant latitude value
!  o the ocn & ice have identical grids, this routine works for both ocn & ice
!
! NOTATION:
!
!           ^  f4      f3
!       b   |
!       e  dy      f                interpolate to f from f1,f2,f3,f4
!       t   |                     
!       a   v  f1      f2
!              <--dx--->
!               alpha
!
!   f  is located at (xm,ym)
!   f1 is located at (xl,yl)
!   f2 is located at (xh,yl)
!   f3 is located at (xh,yh)
!   f4 is located at (xl,yh)
!
!-------------------------------------------------------------------------------

   call mct_aVect_gather(buni%data,gDatai,buni%dom%gsMap,pid0,cpl_comm_comp,rcode)
   call mct_aVect_gather(buni%dom%lGrid,gGridi,buni%dom%gsMap,pid0,cpl_comm_comp,rcode)
   call mct_aVect_gather(buno%data,gDatao,buno%dom%gsMap,pid0,cpl_comm_comp,rcode)
   call mct_aVect_gather(buno%dom%lGrid,gGrido,buno%dom%gsMap,pid0,cpl_comm_comp,rcode)


if (cpl_comm_comp_pid == pid0 ) then

   n_a = buni%dom%n
   n_o = buno%dom%n

   kui = mct_aVect_indexRA(buni%data,fld1,perrWith=subName)
   kvi = mct_aVect_indexRA(buni%data,fld2,perrWith=subName)
   kuo = mct_aVect_indexRA(buno%data,fld1,perrWith=subName)
   kvo = mct_aVect_indexRA(buno%data,fld2,perrWith=subName)

   klati = mct_aVect_indexRA(buni%dom%lGrid,"lat",perrWith=subName)
   kloni = mct_aVect_indexRA(buni%dom%lGrid,"lon",perrWith=subName)
   klato = mct_aVect_indexRA(buno%dom%lGrid,"lat",perrWith=subName)
   klono = mct_aVect_indexRA(buno%dom%lGrid,"lon",perrWith=subName)

   !----------------------------------------------------------------------------
   ! one-time setup computations
   !----------------------------------------------------------------------------
   if (first_call) then
      write(6,F00) " compute bilinear weights & indicies for NP region."

      !--- determine number of points in top row of atm (input) data ---
      mpts = 0
      m1st = 0
      do m=1,n_a
        if (gGridi%rAttr(klati,m) >= gGridi%rAttr(klati,n_a)) then
           m1st = m               ! 1st point in top row
           mpts = n_a - m1st + 1  ! number of points in top row
           exit
        end if 
      end do

      !--- allocate memory for a new, north-pole row ---
      allocate(Sa_in_u(mpts))  ! u-velocity
      allocate(Sa_in_v(mpts))  ! v-velocity
      allocate( theta (mpts))  ! rotation angle wrt prime-meridian

      !--- compute (u,v) rotation angles relative to prime-meridian ---
      do m=m1st,n_a
        mm = m - m1st + 1                   ! index into local NP row
        theta(mm) = gGridi%rAttr(kloni,m)*cpl_const_deg2rad ! units are radians
      end do

      !--- determine number of points that need correcting (output) data ---
      npts = 0
      do n=1,n_o
        if (gGrido%rAttr(klato,n) > gGridi%rAttr(klati,n_a)) npts = npts + 1
      end do

      !--- allocate bilin-map memory given number of output points involved ---
      allocate(   nw(  npts) ) ! indices to output points
      allocate(   mw(4,npts) ) ! for each output point, 4 input  indicies
      allocate( wght(4,npts) ) ! for each output point, 4 bilinear weights

      !--- compute bilin-map data that utilizes new/local NP row ---
      nn = 0
      do n=1,n_o
        if (gGrido%rAttr(klato,n) > gGridi%rAttr(klati,n_a)) then
           nn = nn + 1
           nw(nn) = n              ! output index associated with weights
           xm = gGrido%rAttr(klono,n)     ! x-coord of output
           ym = gGrido%rAttr(klato,n)     ! y-coord of output
           yl = gGridi%rAttr(klati,n_a)   ! constant y-coord of top atm row
           yh = 90.0_R8               ! constant y-coord of north pole
           if ( xm < gGridi%rAttr(kloni,m1st)) then
             m1 = n_a              ! input indicies associated with weights
             m2 = m1st
             m3 = m2   + mpts      ! note: this index goes beyond actual data
             m4 = m1   + mpts      ! note: this index goes beyond actual data
             xl = gGridi%rAttr(kloni,m1) - 360.0_R8
             xh = gGridi%rAttr(kloni,m2)
           else if ( xm < gGridi%rAttr(kloni,n_a) ) then
             do m=m1st,n_a - 1
               if ( xm < gGridi%rAttr(kloni,m+1) ) then
                 m1 = m            ! input indicies associated with weights
                 m2 = m+1
                 m3 = m2  + mpts   ! note: this index goes beyond actual data
                 m4 = m1  + mpts   ! note: this index goes beyond actual data
                 xl = gGridi%rAttr(kloni,m1)
                 xh = gGridi%rAttr(kloni,m2)
                 exit
               end if
             end do
           else
             m1 = n_a              ! input indicies associated with weights
             m2 = m1st
             m3 = m2   + mpts      ! note: this index goes beyond actual data
             m4 = m1   + mpts      ! note: this index goes beyond actual data
             xl = gGridi%rAttr(kloni,m1)
             xh = gGridi%rAttr(kloni,m2) + 360.0_R8
           end if
           dx = xh - xl 
           dy = yh - yl 
           alpha = (xm - xl)/dx
           beta  = (ym - yl)/dy
           wght(1,nn) = (1.0_R8-alpha)*(1.0_R8-beta) ! 4 input bilinear weights
           wght(2,nn) = (    alpha)*(1.0_R8-beta)
           wght(3,nn) = (    alpha)*(    beta)
           wght(4,nn) = (1.0_R8-alpha)*(    beta) 
           mw  (1,nn) = m1                     ! input indicies wrt with weights
           mw  (2,nn) = m2   
           mw  (3,nn) = m3  
           mw  (4,nn) = m4 
        end if
      end do

      first_call = .false.
   end if

   !----------------------------------------------------------------------------
   ! fabricate multiple instances of atm north-pole data w/ various rotations
   !----------------------------------------------------------------------------

   npu = 0  ! avg NP value with y-basis vector pointing north along
   npv = 0  ! prime-meridian and x-basis vector pointing east at prime-meridian
   do mm=1,mpts
      m = n_a - mpts + mm
      npu = npu + cos(theta(mm))*gDatai%rAttr(kui,m) - sin(theta(mm))*gDatai%rAttr(kvi,m)
      npv = npv + sin(theta(mm))*gDatai%rAttr(kui,m) + cos(theta(mm))*gDatai%rAttr(kvi,m)
   end do
   npu = npu/mpts
   npv = npv/mpts

   !--- copies of NP value with basis-vectors aligned with cell's longitude
   do mm=1,mpts
      Sa_in_u(mm) =  cos(theta(mm))*npu + sin(theta(mm))*npv
      Sa_in_v(mm) = -sin(theta(mm))*npu + cos(theta(mm))*npv
   end do

   !----------------------------------------------------------------------------
   ! bilinear interpolation from atm to select set of ocn/ice grid points
   !----------------------------------------------------------------------------
   do nn=1,npts
      n   = nw(nn)            ! output grid index

      m1  = mw(1,nn)          ! input grid indicies on atm grid
      m2  = mw(2,nn)          ! input grid indicies on atm grid
      mm3 = mw(3,nn) - n_a    ! input grid indicies on local "NP grid" array
      mm4 = mw(4,nn) - n_a    ! input grid indicies on local "NP grid" array

      w1 = wght(1,nn)         ! 4 mapping weights
      w2 = wght(2,nn)
      w3 = wght(3,nn)
      w4 = wght(4,nn)

      f1 = gDatai%rAttr(kui,m1) ! 4 input velocities
      f2 = gDatai%rAttr(kui,m2)
      f3 = Sa_in_u(mm3)        ! local NP pole value
      f4 = Sa_in_u(mm4)        ! local NP pole value
      gDatao%rAttr(kuo,n) = w1*f1 + w2*f2 + w3*f3 + w4*f4

      f1 = gDatai%rAttr(kvi,m1)
      f2 = gDatai%rAttr(kvi,m2)
      f3 = Sa_in_v(mm3)        ! local NP pole value
      f4 = Sa_in_v(mm4)        ! local NP pole value
      gDatao%rAttr(kvo,n) = w1*f1 + w2*f2 + w3*f3 + w4*f4

   end do

endif

   call mct_aVect_clean(buno%data) 
   call mct_aVect_scatter(gDatao,buno%data,buno%dom%gsMap,pid0,cpl_comm_comp,rcode)

   if (cpl_comm_comp_pid == pid0 ) then
      call mct_aVect_clean(gDatai)
      call mct_aVect_clean(gGridi)
      call mct_aVect_clean(gDatao)
      call mct_aVect_clean(gGrido)
   endif

end subroutine cpl_map_npFixOld

!===============================================================================
!===============================================================================

end module cpl_map_mod

