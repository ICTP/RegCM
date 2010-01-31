      module parrad

      use regcm_param

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: plon = ix - 1
      integer , parameter :: plev = kx
      integer , parameter :: plat = 1
      integer , parameter :: pcnst = 1
      integer , parameter :: plevmx = 4
      integer , parameter :: plevp = plev + 1
      integer , parameter :: nxpt = 1
      integer , parameter :: jintmx = 1
      integer , parameter :: plond = plon
      integer , parameter :: platd = plat
      integer , parameter :: plevd = plev*(3+pcnst)
      integer , parameter :: plevr = kx
      integer , parameter :: plevrp = plevr + 1
      integer , parameter :: plngbuf = 512*                             &
                           & ((plond*plevrp*plevrp+plond*plevr*4+       &
                           & plond*plevrp)/512+1)
      integer , parameter :: plnlv = plon*plev
      integer , parameter :: plndlv = plond*plev
      integer , parameter :: pbflnb = (7+2*pcnst)*plndlv +              &
                           &          (15+plevmx+pcnst)*plond
      integer , parameter :: pbflna = (3+3*plev)*plond
      integer , parameter :: pbflnm1 = (1+2*plev)*plond
      integer , parameter :: pflenb = ((pbflnb+pbflnm1)/512+1)*512
      integer , parameter :: pflena = (pbflna/512+1)*512
      integer , parameter :: plenalcl = ((pbflna+3*plndlv+plond)/512+1) &
                           & *512
      integer , parameter :: plexbuf = (((1+7*plev)*plond)/512+1)*512
      integer , parameter :: ptapes = 6
      integer , parameter :: pflds = 92 + 8*pcnst + 2*(pcnst-1) + plevmx
      integer , parameter :: ptifld = 11
      integer , parameter :: ptvsfld = 1
      integer , parameter :: ptvofld = 2
      integer , parameter :: plenhis = 37
      integer , parameter :: plenhcs = 89
      integer , parameter :: plenhi = plenhis + 3*pflds
      integer , parameter :: plenhc = plenhcs + 2*pflds
      integer , parameter :: plenhr = 3*(2*plev+1) + 2*plat
      integer , parameter :: ptilenis = plenhis
      integer , parameter :: ptilencs = plenhcs
      integer , parameter :: ptileni = ptilenis + 3*ptifld
      integer , parameter :: ptilenc = ptilencs + 2*ptifld
      integer , parameter :: ptolenis = plenhis
      integer , parameter :: ptolencs = plenhcs
      integer , parameter :: ptvoleni = ptolenis + 3*ptvofld
      integer , parameter :: ptvolenc = ptolencs + 2*ptvofld
      integer , parameter :: ptslenis = plenhis
      integer , parameter :: ptslencs = plenhcs
      integer , parameter :: ptvsleni = ptslenis + 3*ptvsfld
      integer , parameter :: ptvslenc = ptslencs + 2*ptvsfld

      end module parrad
