      subroutine calcrh2d(f2d,nx,ny,nfld,nta,nqa,npsrf,nrh,vmisdat)
 
      use mod_constants , only : svp1 , svp2 , svp3 , ep2
      implicit none
!
! Dummy arguments
!
      integer :: nfld , npsrf , nqa , nrh , nta , nx , ny
      real(4) :: vmisdat
      real(4) , dimension(nx,ny,nfld) :: f2d
      intent (in) nfld , npsrf , nqa , nrh , nta , nx , ny , vmisdat
      intent (inout) f2d
!
! Local variables
!
      integer :: i , j
      real(4) :: pres , qa , qs , satvp , ta
!
      do j = 1 , ny
        do i = 1 , nx
          pres = f2d(i,j,npsrf)
          if ( pres>0. ) then
            ta = f2d(i,j,nta)
            qa = f2d(i,j,nqa)
            if ( ta>273.15 ) then
              satvp = svp1*exp(svp2*(ta-273.15)/(ta-svp3))
                                                         ! SAT'N VAP PRES
            else
              satvp = svp1*exp(22.514-6.15E3/ta)
            end if
            qs = ep2*satvp/(pres-satvp)              ! SAT. MIXING RATIO
            f2d(i,j,nrh) = qa/qs
          else
            f2d(i,j,nrh) = vmisdat
          end if
        end do
      end do
 
      end subroutine calcrh2d
