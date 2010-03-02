      module mod_domain
      implicit none
!
      character(6) , parameter :: iproj = 'LAMCON'
      integer , parameter :: iy = 34 , jx = 51 , kz = 18 , nsg = 1
      real , parameter :: ds = 60.0 , ptop = 5.0 , clat = 45.39 ,       &
                        & clon = 13.48 , plat = clat , plon = clon ,    &
                        & truelatl = 30. , truelath = 60.
      integer , parameter :: ntypec = 10 , ntypec_s = 10
      real , parameter :: h2opct = 75.
      logical , parameter :: ifanal = .true. , smthbdy = .false. ,      &
                           & lakadj = .false.
      integer , parameter :: igrads = 1 , ibigend = 1 , ibyte = 4
      logical , parameter :: fudge_lnd = .false. ,                      &
                           & fudge_lnd_s = .false. ,                    &
                           & fudge_tex = .false. , fudge_tex_s = .false.
      character(50) , parameter :: filout = '../../Input/DOMAIN.INFO' , &
                                 & filctl = '../../Input/DOMAIN.CTL'
      integer , parameter :: idate1 = 1989010100 , idate2 = 1989020100
      character(5) , parameter :: dattyp = 'ERAIN' , ssttyp = 'ERSST'
      logical , parameter :: ehso4 = .false.
      character(4) , parameter :: lsmtyp = 'BATS'
      integer , parameter :: nveg = 20
      character(7) , parameter :: aertyp = 'AER00D0'
      integer , parameter :: ntex = 17 , nproc = 0
      end module mod_domain
