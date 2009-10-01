
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine radtap
c
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/param1.cb'
      include 'Commons/param2.cb'
      include 'Commons/parrad.cb'
      include 'Commons/bats.cb'
#ifdef MPP1
      include 'Common2/bats.cb2'
      include 'Common2/rad.cb'
#else
      include 'Commons/bats.cb2'
      include 'Commons/rad.cb'
#endif
      include 'Commons/date.cb'
      include 'Commons/iunits.cb'
      integer nrad2d,nrad3d
      parameter(nrad2d=21,nrad3d=5)
#ifdef MPP1
      include 'mpiregcm.h'
      real*4  frad2d_io(mjx-2,ix-2,nrad2d)
      real*4  frad3d_io(mjx-2,ix-2,kx,nrad3d)
      COMMON /OUTRADio/ frad2d_io,frad3d_io
#else
      real*4  frad2d(jx-2,ix-2,nrad2d),frad3d(jx-2,ix-2,kx,nrad3d)
      COMMON /OUTRAD/ frad2d,frad3d
#endif
      integer nrcrad
      COMMON /RADREC/ nrcrad
C
C Radiation resolution and I/O parameters
C
c
      integer       i    ! latitude index
      integer       j    ! longitude index
      integer       k    ! level index
      integer n
c
      real*8  cgsmks     ! conversion factor for fluxes from cgs to mks
c
      data cgsmks / 1.e-3 /
c
c compute total radiative heating flux for the surface,
c converting units from cgs to mks:
c
      print *, 'Writing rad fields at ktau = ', ktau,idatex
      if (iotyp.eq.1) then
        do n=1,nrad3d
          if (n.ne.1) then ! skip h20mmr
            do k=kx,1,-1
              nrcrad=nrcrad+1
#ifdef MPP1
              write(iutrad,rec=nrcrad)
     &                     ((frad3d_io(j,i,k,n),j=1,mjx-2),i=1,ix-2)
#else
              write(iutrad,rec=nrcrad)
     &                     ((frad3d(j,i,k,n),j=1,jx-2),i=1,ix-2)
#endif
           enddo
          end if
        enddo
        do n=1,nrad2d
          if (n.lt.10) then ! skip everything from alb (10 to 21)
            nrcrad=nrcrad+1
#ifdef MPP1
            write(iutrad,rec=nrcrad)
     &                   ((frad2d_io(j,i,n),j=1,mjx-2),i=1,ix-2)
#else
            write(iutrad,rec=nrcrad)
     &                   ((frad2d(j,i,n),j=1,jx-2),i=1,ix-2)
#endif
          end if
        enddo
      else if (iotyp.eq.2) then
        write(iutrad) idatex
        do n=1,nrad3d
          if (n.ne.1) then ! skip h20mmr
            do k=kx,1,-1
#ifdef MPP1
              write(iutrad) ((frad3d_io(j,i,k,n),j=1,mjx-2),i=1,ix-2)
#else
              write(iutrad) ((frad3d(j,i,k,n),j=1,jx-2),i=1,ix-2)
#endif
            end do
          end if
        end do
        do n=1,nrad2d
          if (n.lt.10) then ! skip everything from alb (10 to 21)
#ifdef MPP1
            write(iutrad) ((frad2d_io(j,i,n),j=1,mjx-2),i=1,ix-2)
#else
            write(iutrad) ((frad2d(j,i,n),j=1,jx-2),i=1,ix-2)
#endif
          end if
        end do
      end if

      return
      end
