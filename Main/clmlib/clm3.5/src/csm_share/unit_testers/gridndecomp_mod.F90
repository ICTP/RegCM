module GridnDecomp_mod
use shr_kind_mod,  only: SHR_KIND_R8
use shr_sys_mod,   only: shr_sys_abort
implicit none
private

public GridnDecomp_CenterToVertices
public GridnDecomp_CalcDecomp
public GridnDecomp_GetLocalGIndices
public setRandomSeed

  save
  integer :: nlats, nlons, npts
  integer, allocatable :: nproc(:)
  integer, allocatable :: indx(:,:,:)

contains

subroutine GridnDecomp_CenterToVertices( centerLon, centerLat, coordLon, coordLat, tweak )
  implicit none

  real(SHR_KIND_R8), intent(in) :: centerLat(:)
  real(SHR_KIND_R8), intent(in) :: centerLon(:)
  real(SHR_KIND_R8), intent(out) :: coordLon(:)
  real(SHR_KIND_R8), intent(out) :: coordLat(:)
  logical, optional, intent(in)  :: tweak

  integer :: nd
  integer :: lon, lat, irow
  real(SHR_KIND_R8) :: diff1, diff2
  logical :: spectral, tweaked

  if ( present(tweak) )then
     tweaked = tweak
  else
     tweaked = .false.
  end if
  !
  ! Make sure dimensions are correct
  !
  nlats = size(centerLat)
  nlons = size(centerLon)
  nd    = size(coordLon)
  if ( nd /= (nlons+1) )then
     call shr_sys_abort( "dimension of coordlon /= dimension of centerLon+1" )
  end if
  nd    = size(coordLat)
  if ( nd /= (nlats+1) ) &
     call shr_sys_abort( "dimension of coordlat /= dimension of centerLat+1" )

  do lon = 1, nlons-1
       coordLon(lon+1) = (centerLon(lon+1) + centerLon(lon))*0.5_SHR_KIND_R8
  end do
  coordLon(1) =  coordLon(2)-(centerLon(2)-centerLon(1))
  coordLon(nlons+1)   =  (360.0_SHR_KIND_R8 + centerLon(nlons))*0.5_SHR_KIND_R8

  spectral = .false.
  diff1 = abs(centerLat(2) - centerLat(1))
  do lat = 1, nlats-1
      coordLat(lat+1) = (centerLat(lat+1) + centerLat(lat))*0.5_SHR_KIND_R8
      diff2 = abs(centerLat(lat+1) - centerLat(lat))
      if ( mod(nlats,2) == 0 .and. mod(nlons,2) == 0 .and. &
           abs(diff1 - diff2) > 1.0e-5 ) spectral = .true.
  end do
  !  coordLat(nlats)   = 2.0_SHR_KIND_R8*centerLat(nlats) - coordLat(nlats-1) 
  coordLat(1)   = coordLat(2)-(centerLat(2) - centerLat(1)) 
  coordLat(nlats+1) = coordLat(nlats) + (coordLat(nlats) - coordLat(nlats-1))
  !
  ! If this is a spectral grid -- figure it out differently
  ! In this case make sure the vertices average to the centered values
  !
  if ( spectral )then
     coordLat(nlats/2+1) = 0.0_SHR_KIND_R8
     do irow = 1, nlats/2
        coordLat(nlats/2+irow+1) = centerLat(nlats/2+irow)*2.0_SHR_KIND_R8 &
                               - coordLat(nlats/2+irow)
     end do
     do irow = 1, nlats/2
        coordLat(nlats/2+1-irow) = -coordLat(nlats/2+1+irow)
     end do

     write(6,*) 'Spectral: CoordLat = ', CoordLat(:)
  else
     write(6,*) 'FV:       CoordLat = ', CoordLat(:)
  end if
  if ( tweaked )then
     do lon = 3, nlons-2
        coordLon(lon) = coordLon(lon) + 0.5
     end do
     do lat = 3, nlats-2
        coordLat(lat) = coordLat(lat) + 0.5
     end do
  end if

  npts = nlons*nlats
 
end subroutine GridnDecomp_CenterToVertices

subroutine GridnDecomp_CalcDecomp( MasterRank, npes, rank, mpicom, MyCount, landGrid, &
                                   GlobalMask )
  use shr_mpi_mod,             only: shr_mpi_recv, shr_mpi_send, shr_mpi_bcast
  implicit none

  integer, intent(in)  :: MasterRank
  integer, intent(in)  :: npes
  integer, intent(in)  :: rank
  integer, intent(in)  :: mpicom
  integer, intent(out) :: MyCount
  logical, intent(in),  optional :: landGrid
  integer,              optional, pointer :: GlobalMask(:,:)

  integer :: i, j, mxpes, proc
  intrinsic :: random_number
  real(shr_kind_r8) :: rand
  integer, allocatable :: buffer(:)
  logical              :: landGr

  if ( present(landGrid) )then
     landGr = landGrid
  else
     landGr = .false.
  end if
  if ( allocated(nproc) ) deallocate( nproc )
  allocate( nproc(0:npes-1) )
  mxpes = (npts / npes) + 1
  if ( allocated(indx) ) deallocate( indx )
  allocate( indx(0:npes-1,mxpes,2) )
  indx(:,:,:) = -1
  nproc(:) = 0
  if ( present(GlobalMask) ) allocate( GlobalMask(nlons,nlats) )
  if ( MasterRank == Rank)then
    do j = 1, nlats

       do i = 1, nlons
          call random_number(rand) 
          proc = nint(rand*(npes-1))
          call random_number(rand) 
          if ( landGr .and. (rand < 0.3) )then
              if ( present(GlobalMask) ) GlobalMask(i,j) = 0
              cycle   ! Drop random points for land
          else
              if ( present(GlobalMask) ) GlobalMask(i,j) = 1
          end if
          if (nproc(proc) >= mxpes)then
             do while( nproc(proc) >= mxpes )
               call random_number(rand) 
               proc = nint(rand*(npes-1))
             end do
          end if
          nproc(proc) = nproc(proc) + 1
          indx(proc,nproc(proc),1) = i
          indx(proc,nproc(proc),2) = j
       end do

    end do
    if ( present(GlobalMask) ) call shr_mpi_bcast( GlobalMask(:,:), mpicom )
    do proc = 0, npes-1
       if ( proc /= MasterRank )then
          call shr_mpi_send(nproc,proc,211,mpicom)
          allocate(buffer(nproc(proc)))
          buffer(:) = indx(proc,1:nproc(proc),1)
          call shr_mpi_send(buffer,proc,411,mpicom)
          buffer(:) = indx(proc,1:nproc(proc),2)
          call shr_mpi_send(buffer,proc,611,mpicom)
          deallocate(buffer)
       end if
    end do
  else
     if ( present(GlobalMask) ) call shr_mpi_bcast( GlobalMask(:,:), mpicom )
     call shr_mpi_recv(nproc,MasterRank,211,mpicom)
     allocate(buffer(nproc(rank)))
     call shr_mpi_recv(buffer,MasterRank,411,mpicom)
     indx(rank,1:nproc(rank),1) = buffer(:)
     call shr_mpi_recv(buffer,MasterRank,611,mpicom)
     indx(rank,1:nproc(rank),2) = buffer(:)
     deallocate(buffer)
  end if
  MyCount = nproc(rank)
end subroutine GridnDecomp_CalcDecomp

subroutine GridnDecomp_GetLocalGIndices( MyCount, rank, MyIndices )
  implicit none

  integer, intent(in) :: MyCount
  integer, intent(in) :: rank
  integer, intent(out) :: MyIndices(MyCount,2)

  MyIndices(1:MyCount,1) = indx(rank,1:MyCount,1)
  MyIndices(1:MyCount,2) = indx(rank,1:MyCount,2)
  if ( any(MyIndices(1:MyCount,1:2) == -1) ) &
     call shr_sys_abort( "Index == -1 -- somethings wrong" )
end subroutine GridnDecomp_GetLocalGIndices

subroutine SetRandomSeed( value )
  integer, intent(in) :: value

  integer, allocatable :: seed(:)
  integer              :: seedsize

  call random_seed( size=seedsize )
  allocate( seed(seedsize) )
  seed(:) = value
  call random_seed( put=seed )
  deallocate( seed )

end subroutine SetRandomSeed

end module GridnDecomp_mod
