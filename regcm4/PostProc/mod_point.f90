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

      module mod_point
      implicit none
      integer :: ndrag , net , nlwd , nlwn , nprc , npsrf , npt ,       &
               & nqanm , nrha , nrnfs , nsh , nsmr , nsmu , nsnow ,     &
               & nswi , nswn , ntamax , ntamin , ntanm , ntf , ntg ,    &
               & ntgmax , ntgmin , nux , nvx , kxpbl , psmin , w10max
      integer :: npsi , nslpi , ntgi
      integer :: ndivi , nhgti , nqvi , nrhi , ntdi , nthi , nti , nui ,&
               & nvi , nvori
      integer :: nbf , npsa , nrt , nslp , nsmt , ntgb
      integer :: ndiva , nhgt , nomega , nqca , nqva , nrh , nta ,      &
               & ntda , ntha , nua , nva , nvora
      integer :: nclrls , nclrlt , nclrss , nclrst , nfirtp , nflw ,    &
               & nfsw , nsabtp , nsolin
      integer :: ncld , nclwp , nqrl , nqrs
      integer :: nsdrag , nset , nsprc , nspsrf , nspt , nsqanm ,       &
               & nsrha , nsrnfs , nssh , nssmr , nssmu , nssnow ,       &
               & nstanm , nstf , nstg , nsux , nsvx

      data nui , nvi , nti , nqvi , nrhi , ntdi , nthi/1 , 2 , 3 , 4 ,  &
         & 5 , 6 , 7/
      data nvori , ndivi , nhgti/8 , 9 , 10/

      data npsi , ntgi , nslpi/1 , 2 , 3/

      data nua , nva , nomega , nta , nqva , nqca , nrh , nhgt/1 , 2 ,  &
         & 3 , 4 , 5 , 6 , 7 , 8/
      data ntha , ntda , nvora , ndiva/9 , 10 , 11 , 12/

      data npsa , nrt , ntgb , nsmt , nbf , nslp/1 , 2 , 3 , 4 , 5 , 6/

      data nux , nvx , ndrag , ntg , ntf , ntanm , nqanm/1 , 2 , 3 , 4 ,&
         & 5 , 6 , 7/
      data nsmu , nsmr , npt , net , nrnfs , nsnow , nsh/8 , 9 , 10 ,   &
         & 11 , 12 , 13 , 14/
      data nlwn , nswn , nlwd , nswi , nprc , npsrf , kxpbl/15 , 16 ,   &
         & 17 , 18 , 19 , 20 , 21/
      data ntgmax , ntgmin , ntamax , ntamin , w10max , psmin/22 , 23 , &
         & 24 , 25 , 26 , 27/
      data nrha/28/

      data nsux , nsvx , nsdrag , nstg , nstf , nstanm/1 , 2 , 3 , 4 ,  &
         & 5 , 6/
      data nsqanm , nssmu , nssmr , nspt , nset , nsrnfs/7 , 8 , 9 ,    &
         & 10 , 11 , 12/
      data nssnow , nssh , nsprc , nspsrf , nsrha/13 , 14 , 15 , 16 ,   &
         & 17/

      data ncld , nclwp , nqrs , nqrl/1 , 2 , 3 , 4/
      data nfsw , nflw , nclrst , nclrss , nclrlt/1 , 2 , 3 , 4 , 5/
      data nclrls , nsolin , nsabtp , nfirtp/6 , 7 , 8 , 9/

      end module mod_point
