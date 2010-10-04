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

      module mod_fudge

      implicit none

      contains

      subroutine lndfudge(fudge,ch,lndout,htgrid,iy,jx,char_lnd)

      implicit none
!
! Dummy arguments
!
      character(*) :: char_lnd
      logical :: fudge,there
      integer :: iy , jx
      character(1) , dimension(iy,jx) :: ch
      real(4) , dimension(iy,jx) :: htgrid , lndout
      intent (in) char_lnd , fudge , iy , jx
      intent (inout) ch , htgrid , lndout
!
! Local variables
!
      integer :: i , j
!
      if ( fudge ) then
        inquire (file=char_lnd,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' , char_lnd ,                          &
             &' FILE:  FILE DOES NOT EXIST'
           stop ' IN SUBROUTINE lndfudge'
        endif 
        open (13,file=char_lnd,form='formatted')
        do i = iy , 1 , -1
          read (13,99001) (ch(i,j),j=1,jx)
        end do
        close (13)
        do j = 1 , jx
          do i = 1 , iy
            if ( ch(i,j)==' ' ) then
              lndout(i,j) = 15.
            else if ( ch(i,j)=='1' ) then
              lndout(i,j) = 1.
            else if ( ch(i,j)=='2' ) then
              lndout(i,j) = 2.
            else if ( ch(i,j)=='3' ) then
              lndout(i,j) = 3.
            else if ( ch(i,j)=='4' ) then
              lndout(i,j) = 4.
            else if ( ch(i,j)=='5' ) then
              lndout(i,j) = 5.
            else if ( ch(i,j)=='6' ) then
              lndout(i,j) = 6.
            else if ( ch(i,j)=='7' ) then
              lndout(i,j) = 7.
            else if ( ch(i,j)=='8' ) then
              lndout(i,j) = 8.
            else if ( ch(i,j)=='9' ) then
              lndout(i,j) = 9.
            else if ( ch(i,j)=='A' ) then
              lndout(i,j) = 10.
            else if ( ch(i,j)=='B' ) then
              lndout(i,j) = 11.
            else if ( ch(i,j)=='C' ) then
              lndout(i,j) = 12.
            else if ( ch(i,j)=='D' ) then
              lndout(i,j) = 13.
            else if ( ch(i,j)=='E' ) then
              lndout(i,j) = 14.
            else if ( ch(i,j)=='F' ) then
              lndout(i,j) = 15.
            else if ( ch(i,j)=='G' ) then
              lndout(i,j) = 16.
            else if ( ch(i,j)=='H' ) then
              lndout(i,j) = 17.
            else if ( ch(i,j)=='I' ) then
              lndout(i,j) = 18.
            else if ( ch(i,j)=='J' ) then
              lndout(i,j) = 19.
            else if ( ch(i,j)=='K' ) then
              lndout(i,j) = 20.
            else if ( nint(lndout(i,j))==0 ) then
!               ch(i,j) = 'X'
              ch(i,j) = ' '
            else
              write (*,*) 'LANDUSE MASK exceed the limit'
              stop
            end if
!_fix         if(nint(lndout(i,j)).eq.15) htgrid(i,j) = 0.0
            if ( htgrid(i,j)<0.1 .and. nint(lndout(i,j))==15 )        &
               & htgrid(i,j) = 0.0
          end do
        end do
      else
        do j = 1 , jx
          do i = 1 , iy
            if ( nint(lndout(i,j))==15 .or. nint(lndout(i,j))==0 ) then
              ch(i,j) = ' '
            else if ( nint(lndout(i,j))==1 ) then
              ch(i,j) = '1'
            else if ( nint(lndout(i,j))==2 ) then
              ch(i,j) = '2'
            else if ( nint(lndout(i,j))==3 ) then
              ch(i,j) = '3'
            else if ( nint(lndout(i,j))==4 ) then
              ch(i,j) = '4'
            else if ( nint(lndout(i,j))==5 ) then
              ch(i,j) = '5'
            else if ( nint(lndout(i,j))==6 ) then
              ch(i,j) = '6'
            else if ( nint(lndout(i,j))==7 ) then
              ch(i,j) = '7'
            else if ( nint(lndout(i,j))==8 ) then
              ch(i,j) = '8'
            else if ( nint(lndout(i,j))==9 ) then
              ch(i,j) = '9'
            else if ( nint(lndout(i,j))==10 ) then
              ch(i,j) = 'A'
            else if ( nint(lndout(i,j))==11 ) then
              ch(i,j) = 'B'
            else if ( nint(lndout(i,j))==12 ) then
              ch(i,j) = 'C'
            else if ( nint(lndout(i,j))==13 ) then
              ch(i,j) = 'D'
            else if ( nint(lndout(i,j))==14 ) then
              ch(i,j) = 'E'
            else if ( nint(lndout(i,j))==16 ) then
              ch(i,j) = 'G'
            else if ( nint(lndout(i,j))==17 ) then
              ch(i,j) = 'H'
            else if ( nint(lndout(i,j))==18 ) then
              ch(i,j) = 'I'
            else if ( nint(lndout(i,j))==19 ) then
              ch(i,j) = 'J'
            else if ( nint(lndout(i,j))==20 ) then
              ch(i,j) = 'K'
            else
              write (*,*) 'LANDUSE MASK' , nint(lndout(i,j)) ,        &
                         &'exceed the limit'
              stop
            end if
          end do
        end do
        open (13,file=char_lnd,form='formatted')
        do i = iy , 1 , -1
          write (13,99001) (ch(i,j),j=1,jx)
        end do
        close (13)
      end if
99001 format (132A1)
      end subroutine lndfudge

      subroutine texfudge(fudge,ch,texout,htgrid,iy,jx,char_tex)
      implicit none
!
! Dummy arguments
!
      character(*) :: char_tex
      logical :: fudge, there
      integer :: iy , jx
      character(1) , dimension(iy,jx) :: ch
      real(4) , dimension(iy,jx) :: htgrid , texout
      intent (in) char_tex , fudge , iy , jx
      intent (out) htgrid
      intent (inout) ch , texout
!
! Local variables
!
      integer :: i , j
!
      if ( fudge ) then
         inquire (file=char_tex,exist=there)
             if ( .not.there ) then
               print * , 'ERROR OPENING ' , char_tex ,                          &
               &' FILE:  FILE DOES NOT EXIST'
               stop ' IN SUBROUTINE texfudge'
             endif 
        open (13,file=char_tex,form='formatted')
        do i = iy , 1 , -1
          read (13,99001) (ch(i,j),j=1,jx)
        end do
        close (13)
        do j = 1 , jx
          do i = 1 , iy
            if ( ch(i,j)==' ' ) then
              texout(i,j) = 14.
            else if ( ch(i,j)=='1' ) then
              texout(i,j) = 1.
            else if ( ch(i,j)=='2' ) then
              texout(i,j) = 2.
            else if ( ch(i,j)=='3' ) then
              texout(i,j) = 3.
            else if ( ch(i,j)=='4' ) then
              texout(i,j) = 4.
            else if ( ch(i,j)=='5' ) then
              texout(i,j) = 5.
            else if ( ch(i,j)=='6' ) then
              texout(i,j) = 6.
            else if ( ch(i,j)=='7' ) then
              texout(i,j) = 7.
            else if ( ch(i,j)=='8' ) then
              texout(i,j) = 8.
            else if ( ch(i,j)=='9' ) then
              texout(i,j) = 9.
            else if ( ch(i,j)=='A' ) then
              texout(i,j) = 10.
            else if ( ch(i,j)=='B' ) then
              texout(i,j) = 11.
            else if ( ch(i,j)=='C' ) then
              texout(i,j) = 12.
            else if ( ch(i,j)=='D' ) then
              texout(i,j) = 13.
            else if ( ch(i,j)=='E' ) then
              texout(i,j) = 14.
            else if ( ch(i,j)=='F' ) then
              texout(i,j) = 15.
            else if ( ch(i,j)=='G' ) then
              texout(i,j) = 16.
            else if ( ch(i,j)=='H' ) then
              texout(i,j) = 17.
            else if ( nint(texout(i,j))==0 ) then
!             ch(i,j) = 'X'
              ch(i,j) = ' '
            else
              write (*,*) 'TEXTURE TYPE exceed the limit'
              stop
            end if
            if ( nint(texout(i,j))==14 ) htgrid(i,j) = 0.0
          end do
        end do
      else
        do j = 1 , jx
          do i = 1 , iy
            if ( nint(texout(i,j))==14 ) then
              ch(i,j) = ' '
            else if ( nint(texout(i,j))==1 ) then
              ch(i,j) = '1'
            else if ( nint(texout(i,j))==2 ) then
              ch(i,j) = '2'
            else if ( nint(texout(i,j))==3 ) then
              ch(i,j) = '3'
            else if ( nint(texout(i,j))==4 ) then
              ch(i,j) = '4'
            else if ( nint(texout(i,j))==5 ) then
              ch(i,j) = '5'
            else if ( nint(texout(i,j))==6 ) then
              ch(i,j) = '6'
            else if ( nint(texout(i,j))==7 ) then
              ch(i,j) = '7'
            else if ( nint(texout(i,j))==8 ) then
              ch(i,j) = '8'
            else if ( nint(texout(i,j))==9 ) then
              ch(i,j) = '9'
            else if ( nint(texout(i,j))==10 ) then
              ch(i,j) = 'A'
            else if ( nint(texout(i,j))==11 ) then
              ch(i,j) = 'B'
            else if ( nint(texout(i,j))==12 ) then
              ch(i,j) = 'C'
            else if ( nint(texout(i,j))==13 ) then
              ch(i,j) = 'D'
            else if ( nint(texout(i,j))==15 ) then
              ch(i,j) = 'F'
            else if ( nint(texout(i,j))==16 ) then
              ch(i,j) = 'G'
            else if ( nint(texout(i,j))==17 ) then
              ch(i,j) = 'H'
            else
              write (*,*) 'TEXTURE TYPE' , nint(texout(i,j)) ,          &
                         &'exceed the limit'
              stop
            end if
          end do
        end do
        open (13,file=char_tex,form='formatted')
        do i = iy , 1 , -1
          write (13,99001) (ch(i,j),j=1,jx)
        end do
        close (13)
      end if
99001 format (132A1)
      end subroutine texfudge

      subroutine lakeadj(lnduse,htgrid,dhlake,xlat,xlon,imx,jmx,i_lake)
 
      implicit none
!
! PARAMETER definitions
!
      real(4) , parameter :: zsuperior   = 183.
      real(4) , parameter :: zhuron      = 176.
      real(4) , parameter :: zmichigan   = 176.
      real(4) , parameter :: zerie       = 174.
      real(4) , parameter :: zontario    =  75.

      real(4) , parameter :: zcaspian    = -28.
      real(4) , parameter :: zmatano     = 382.

      real(4) , parameter :: zvictoria   =1134.
      real(4) , parameter :: zaral       =  53.
      real(4) , parameter :: ztanganyima = 773.
      real(4) , parameter :: zbaikal     = 456.
      real(4) , parameter :: zgreatbear  = 186.
      real(4) , parameter :: zwinnipeg   = 217.
      real(4) , parameter :: zbalkhash   = 341.4
      real(4) , parameter :: zladoga     =   4.8
      real(4) , parameter :: zsap        =   5.
      real(4) , parameter :: zmaracaibo  =  -1.
      real(4) , parameter :: ztungting   =  11.
      real(4) , parameter :: zpatos      =  -1.
      
      real(4) , parameter :: zonega      =  35.
      real(4) , parameter :: zeyre       =  -9.5
      real(4) , parameter :: zvolta      =  85.
      real(4) , parameter :: ztiticaca   =3812.
      real(4) , parameter :: zleopold2   = 340.
      real(4) , parameter :: znicaragua  =  32.
      real(4) , parameter :: zathabasca  = 213.
      real(4) , parameter :: zturkana    = 360.4
      real(4) , parameter :: zsmallwood  = 471.
      real(4) , parameter :: zissyk_kool =1606.
      real(4) , parameter :: zaswan_dam  = 183.
      real(4) , parameter :: zurmia      =1270.
      real(4) , parameter :: zkujbyshevsk=  53.
      real(4) , parameter :: ztorrens    =  30.
      real(4) , parameter :: zreindeer   = 337.
      real(4) , parameter :: zvanern     =  44.
      real(4) , parameter :: zbukhtarmins= 387.
      real(4) , parameter :: zbratskoye  = 402.
      real(4) , parameter :: zkariba     = 485.
      real(4) , parameter :: zalbert     = 615.
      real(4) , parameter :: znasser     = 183.
      real(4) , parameter :: zwinnipegosi= 253.
      real(4) , parameter :: znetilling  =  30.
      real(4) , parameter :: zgreatsalt  =1283.
      real(4) , parameter :: zmanitoba   = 248.
      real(4) , parameter :: zchienghai  =3196.
      real(4) , parameter :: ztaimyr     =   3.7
      real(4) , parameter :: zrybinsk    =102.4
      real(4) , parameter :: znipigon    = 320.
      real(4) , parameter :: zgairdner   =  34.
      real(4) , parameter :: zmweru      = 922.
      real(4) , parameter :: zwoods      = 323.
      real(4) , parameter :: zpeipus     =  30.
      real(4) , parameter :: zsobradinho = 392.5
      real(4) , parameter :: zkhanka     =  68.9
      real(4) , parameter :: zvan        =1646.
      real(4) , parameter :: zdubawnt    = 236.
      real(4) , parameter :: ztana       =1788.
      real(4) , parameter :: zchudsko    =  30.
      real(4) , parameter :: zuvs        = 759.
      real(4) , parameter :: zpoyang     =  16.5
      real(4) , parameter :: zvolgogradsk=  15.
      real(4) , parameter :: zmelville   =  -1.
      real(4) , parameter :: zamadjuak   = 113.
      real(4) , parameter :: zlop_nor    = 780.
      real(4) , parameter :: zhelmand    = 488.
      real(4) , parameter :: zrukwa      = 793.
      real(4) , parameter :: zmirim      =  -1.
      real(4) , parameter :: zcaniapiscau= 535.
      real(4) , parameter :: zhovsgol    =1645.
      real(4) , parameter :: zdongting   =  33.5
      real(4) , parameter :: zcaborabassa= 314.
      real(4) , parameter :: ztsimlyansko=  36.
      real(4) , parameter :: zwollaston  = 398.
      real(4) , parameter :: zalakol     = 347.
      real(4) , parameter :: ziliamna    =  15.
      real(4) , parameter :: zdensu      =  14.1
      real(4) , parameter :: znam        =4627.
      real(4) , parameter :: zlagrande2  = 173.5
      real(4) , parameter :: ztucurui    =  72.
      real(4) , parameter :: ztai        =   3.1
      real(4) , parameter :: zlagrande3  = 256.
      real(4) , parameter :: zsouthernind= 258.
      real(4) , parameter :: zbalbina    =  50.
      real(4) , parameter :: zedward     = 912.
      real(4) , parameter :: zkrementchug=  81.
      real(4) , parameter :: zprimavera  = 259.
      real(4) , parameter :: zbuenosaires= 217.
      real(4) , parameter :: zkivu       =1460.
      real(4) , parameter :: zviluyskoe  = 244.
      real(4) , parameter :: zkakhovskoye=  16.
      real(4) , parameter :: zzeiskoye   = 315.
      real(4) , parameter :: zmistassini = 375.
      real(4) , parameter :: zmichikamau = 460.
      real(4) , parameter :: znueltin    = 278.
      real(4) , parameter :: zchany      = 106.
      real(4) , parameter :: zkrasnoyarsk= 243.
      real(4) , parameter :: zmarchiquita=  69.
      real(4) , parameter :: zhungtze    =  12.3
      real(4) , parameter :: zmanicouagan= 366.
      real(4) , parameter :: zhammer     =  10.
      real(4) , parameter :: znamu       =4718.
      real(4) , parameter :: zust_ilimsko= 296.
      real(4) , parameter :: zkamskoye   = 108.5
      real(4) , parameter :: zokeechobee =   4.
      real(4) , parameter :: zvattern    =  89.
      real(4) , parameter :: zcarrera    = 215.
      real(4) , parameter :: zkaptchagays= 485.
      real(4) , parameter :: zsaratovskoy=  28.
      real(4) , parameter :: zbaker      =   2.
      real(4) , parameter :: zwilliston  = 672.1
      real(4) , parameter :: zmartre     = 265.
      real(4) , parameter :: zhar_us     =1153.
      real(4) , parameter :: zsaimaa     =  76.
      real(4) , parameter :: zchilwa     = 622.
      real(4) , parameter :: zhulun      = 543.
      real(4) , parameter :: zkyoga      = 914.
      real(4) , parameter :: zterminos   =  -1.
      real(4) , parameter :: zqilin      =4530.
      real(4) , parameter :: zkurisches  =  -1.
      real(4) , parameter :: zpontchartr =   1.
      real(4) , parameter :: zyacyreta   =  82.
      real(4) , parameter :: zgorkovskoye=  84.
      real(4) , parameter :: ztengiz     = 304.
      real(4) , parameter :: zchad       = 280.
      real(4) , parameter :: zbangweulu  =1140.
      real(4) , parameter :: ztuz        = 925.
      real(4) , parameter :: zargentino  = 187.
      real(4) , parameter :: zseul       = 357.
      real(4) , parameter :: zclaire     = 213.
      real(4) , parameter :: zhyargas    =1029.
      real(4) , parameter :: zselawik    =   1.
      real(4) , parameter :: ztangra     =4724.
      real(4) , parameter :: zbaghrash   =1038.
      real(4) , parameter :: zmanzala    =  -1.
      real(4) , parameter :: zsevan      =1905.
      real(4) , parameter :: zitaipu     = 220.
      real(4) , parameter :: zmoose      = 255.
      real(4) , parameter :: zronge      = 366.
      real(4) , parameter :: zyathkyed   = 141.
      real(4) , parameter :: zkasba      = 336.
      real(4) , parameter :: zcedar      = 256.
      real(4) , parameter :: zgouin      = 404.
      real(4) , parameter :: zluang      =  -1.
      real(4) , parameter :: zsheksninsko= 113.
      real(4) , parameter :: zeau_claire = 238.
      real(4) , parameter :: zvygozersko =  89.3
      real(4) , parameter :: zilha_soltei= 328.
      real(4) , parameter :: zbecharof   =   4.
      real(4) , parameter :: zlesser     = 577.
      real(4) , parameter :: zred        = 358.
      real(4) , parameter :: zabaya      =1285.
      real(4) , parameter :: zcree       = 487.
      real(4) , parameter :: zmalaren    =   0.3
      real(4) , parameter :: zchamplain  =  29.38
      real(4) , parameter :: zvotkinskoye=  89.
      real(4) , parameter :: zstclain    = 175.
      real(4) , parameter :: zchapala    =1524.
      real(4) , parameter :: zcaratasca  =  -1.
      real(4) , parameter :: zaberdeen   =  80.
      real(4) , parameter :: zchongon    =  15.
      real(4) , parameter :: zpaijanne   =  78.
      real(4) , parameter :: ztoba       = 905.
      real(4) , parameter :: zbas_dor    =  -1.
      real(4) , parameter :: zweishan    =  32.
      real(4) , parameter :: znovosibirsk= 113.7
      real(4) , parameter :: zviedma     = 250.
      real(4) , parameter :: zsongkhla   =   0.
      real(4) , parameter :: zebi        = 213.
      real(4) , parameter :: zgods       = 178.
      real(4) , parameter :: zsaint_john =  99.6
      real(4) , parameter :: zinari      = 119.
      real(4) , parameter :: zbienville  = 391.
      real(4) , parameter :: zisland     =  20.1
      real(4) , parameter :: zmanagua    =  37.4
      real(4) , parameter :: zopinaca    = 216.
      real(4) , parameter :: zohiggins   = 250.
      real(4) , parameter :: ztakiyuak   = 381.
      real(4) , parameter :: zbositeng   =1048.
!
! Dummy arguments
!
      integer :: imx , jmx , i_lake
      real(4) , dimension(imx,jmx) :: htgrid , dhlake , xlat , xlon
      integer , dimension(imx,jmx) :: lnduse
      intent (in) imx , jmx , i_lake , xlat , xlon
      intent (inout) htgrid, dhlake, lnduse
!
! Local variables
!
      integer :: i , j
      real(4) :: xx , yy , ds1deg , conv05
!
!     ****  ADJUST GREAT LAKE ELEVATION **** C
!
      ds1deg = 6371.*3.1415926/180.
      conv05 = 3.1415926/360.

      do i = 1 , imx
        do j = 1 , jmx
          if ( lnduse(i,j)==14 ) then
            xx = xlon(i,j)
            if(xx.gt.180.)  xx = xx - 360.
            if(xx.le.-180.) xx = xx + 360.
            yy = xlat(i,j)

            if ( yy<=47.0 .and. yy>=41.0 .and. xx<= 42.0 .and.          &
               & xx>=27.0 ) then       ! Black Sea should be ocean
               htgrid(i,j) = 0.0
               lnduse(i,j) = 15
               write(*,*) 'Black Sea :', 0.0 , i , j

            else if ( yy<=43.2 .and. yy>=41.0 .and. xx<=-78.0 .and.     &
               & xx>=-84.0 ) then                       ! LAKE ERIE
              print * , '**** ADUJUSTING LAKE ERIE LEVEL ****'
              print * , '     NEW:' , zerie , '    OLD:' , htgrid(i,j) ,&
                  & i , j
              htgrid(i,j) = zerie
              if(i_lake==1) dhlake(i,j) = 17.7 ! 64.0 
            else if ( yy<=46.4 .and. yy>=43.0 .and. xx<=-79.9 .and.     &
                    & yy>=-85.0 ) then                  ! LAKE HURON
              print * , '**** ADUJUSTING LAKE HURON LEVEL ****'
              print * , '     NEW:' , zhuron , '    OLD:' , htgrid(i,j) &
                  & , i , j
              htgrid(i,j) = zhuron
              if(i_lake==1) dhlake(i,j) = 53.0 ! 228.0 
            else if ( yy<=44.5 .and. yy>=43.2 .and. xx<=-75.0 .and.     &
                    & yy>=-79.9 ) then                  ! LAKE ONTARIO
              print * , '**** ADUJUSTING LAKE ONTARIO LEVEL ****'
              print * , '     NEW:' , zontario , '    OLD:' , htgrid(i,j) &
                  & , i , j
              htgrid(i,j) = zontario
              if(i_lake==1) dhlake(i,j) = 86.0 ! 224.0 
            else if ( yy<=49.4 .and. yy>=46.2 .and. xx<=-84.2 .and.     &
                    & xx>=-93.0 ) then                  ! LAKE SUPERIOR
              print * , '**** ADUJUSTING LAKE SUPERIOR LEVEL ****'
              print * , '     NEW:' , zsuperior , '    OLD:' , htgrid(i,j) , &
                  & i , j
              htgrid(i,j) = zsuperior
              if(i_lake==1) dhlake(i,j) = 148.0 ! 406.0 
            else if ( yy<=46.2 .and. yy>=41.0 .and. xx<=-84.8 .and.     &
                    & xx>=-89.0 ) then                  ! LAKE MICHIGAN
              print * , '**** ADUJUSTING LAKE MICHIGAN LEVEL ****'
              print * , '     NEW:' , zmichigan , '    OLD:' , htgrid(i,j) ,&
                  & i , j
              htgrid(i,j) = zmichigan
              if(i_lake==1) dhlake(i,j) = 84.0 ! 281.0 
!           else if(ds1deg*sqrt( ((xx+88.28)*cos((yy+47.53)*conv05))**2  &
!                             +(yy-47.53)**2).lt.2.*sqrt(82367.)) then
!             print * , '**** ADUJUSTING LAKE superior LEVEL ****'
!             print * , '     NEW:' , zsuperior , '    OLD:' , htgrid(i,j) ,&
!                & i , j
!             htgrid(i,j) = zsuperior
!             if(i_lake==1) dhlake(i,j) = 148.0 ! 406.0 
!           else if(ds1deg*sqrt( ((xx+82.28)*cos((yy+44.67)*conv05))**2  &
!                             +(yy-44.67)**2).lt.2.*sqrt(59570.)) then
!             print * , '**** ADUJUSTING LAKE huron LEVEL ****'
!             print * , '     NEW:' , zhuron , '    OLD:' , htgrid(i,j) ,&
!                & i , j
!             htgrid(i,j) = zhuron
!             if(i_lake==1) dhlake(i,j) = 53.0 ! 228.0 
!           else if(ds1deg*sqrt( ((xx+86.23)*cos((yy+43.90)*conv05))**2  &
!                             +(yy-43.90)**2).lt.2.*sqrt(58016.)) then
!             print * , '**** ADUJUSTING LAKE michigan LEVEL ****'
!             print * , '     NEW:' , zmichigan , '    OLD:' , htgrid(i,j) ,&
!                & i , j
!             htgrid(i,j) = zmichigan
!             if(i_lake==1) dhlake(i,j) = 84.0 ! 281.0 
!           else if(ds1deg*sqrt( ((xx+81.13)*cos((yy+41.73)*conv05))**2  &
!                             +(yy-41.73)**2).lt.2.*sqrt(25821.)) then
!             print * , '**** ADUJUSTING LAKE erie LEVEL ****'
!             print * , '     NEW:' , zerie , '    OLD:' , htgrid(i,j) ,&
!                & i , j
!             htgrid(i,j) = zerie
!             if(i_lake==1) dhlake(i,j) = 17.7 ! 64.0 
!           else if(ds1deg*sqrt( ((xx+78.03)*cos((yy+43.70)*conv05))**2  &
!                             +(yy-43.70)**2).lt.2.*sqrt(19009.)) then
!             print * , '**** ADUJUSTING LAKE ontario LEVEL ****'
!             print * , '     NEW:' , zontario , '    OLD:' , htgrid(i,j) ,&
!                & i , j
!             htgrid(i,j) = zontario
!             if(i_lake==1) dhlake(i,j) = 86.0 ! 224.0 

            else if(ds1deg*sqrt( ((xx-50.00)*cos((yy+42.00)*conv05))**2  &
                              +(yy-42.00)**2).lt.2.*sqrt(374000.)) then
              print * , '**** ADUJUSTING LAKE caspian LEVEL ****'
              print * , '     NEW:' , zcaspian , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zcaspian
              if(i_lake==1) dhlake(i,j) = 182.0 ! 1025.0 
            else if(ds1deg*sqrt( ((xx-121.48)*cos((yy-2.43)*conv05))**2  &
                              +(yy+2.43)**2).lt.2.*sqrt(164080.)) then
              print * , '**** ADUJUSTING LAKE matano LEVEL ****'
              print * , '     NEW:' , zmatano , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zmatano
              if(i_lake==1) dhlake(i,j) = 300.0 ! 590.0 

            else if(ds1deg*sqrt( ((xx-33.27)*cos((yy-1.67)*conv05))**2  &
                              +(yy+1.67)**2).lt.2.*sqrt(68800.)) then
              print * , '**** ADUJUSTING LAKE victoria LEVEL ****'
              print * , '     NEW:' , zvictoria , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zvictoria
              if(i_lake==1) dhlake(i,j) = 40.0 ! 84.0 
            else if(ds1deg*sqrt( ((xx-60.00)*cos((yy+45.00)*conv05))**2  &
                              +(yy-45.00)**2).lt.2.*sqrt(64500.)) then
              print * , '**** ADUJUSTING LAKE aral LEVEL ****'
              print * , '     NEW:' , zaral , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zaral
              if(i_lake==1) dhlake(i,j) = 16.0 ! 67.0 
            else if(ds1deg*sqrt( ((xx-30.17)*cos((yy-6.08)*conv05))**2  &
                              +(yy+6.08)**2).lt.2.*sqrt(32000.)) then
              print * , '**** ADUJUSTING LAKE tanganyima LEVEL ****'
              print * , '     NEW:' , ztanganyima , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = ztanganyima
              if(i_lake==1) dhlake(i,j) = 572.0 ! 1471.0 
            else if(ds1deg*sqrt( ((xx-106.78)*cos((yy+53.62)*conv05))**2  &
                              +(yy-53.62)**2).lt.2.*sqrt(31500.)) then
              print * , '**** ADUJUSTING LAKE baikal LEVEL ****'
              print * , '     NEW:' , zbaikal , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zbaikal
              if(i_lake==1) dhlake(i,j) = 740.0 ! 1741.0 
            else if(ds1deg*sqrt( ((xx+120.50)*cos((yy+66.00)*conv05))**2  &
                              +(yy-66.00)**2).lt.2.*sqrt(31153.)) then
              print * , '**** ADUJUSTING LAKE greatbear LEVEL ****'
              print * , '     NEW:' , zgreatbear , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zgreatbear
              if(i_lake==1) dhlake(i,j) = 71.7 ! 446.0 
            else if(ds1deg*sqrt( ((xx+97.77)*cos((yy+52.10)*conv05))**2  &
                              +(yy-52.10)**2).lt.2.*sqrt(23750.)) then
              print * , '**** ADUJUSTING LAKE winnipeg LEVEL ****'
              print * , '     NEW:' , zwinnipeg , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zwinnipeg
              if(i_lake==1) dhlake(i,j) = 12.0 ! 36.0 
            else if(ds1deg*sqrt( ((xx-76.42)*cos((yy+45.73)*conv05))**2  &
                              +(yy-45.73)**2).lt.2.*sqrt(18200.)) then
              print * , '**** ADUJUSTING LAKE balkhash LEVEL ****'
              print * , '     NEW:' , zbalkhash , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zbalkhash
              if(i_lake==1) dhlake(i,j) = 5.8 ! 25.6 
            else if(ds1deg*sqrt( ((xx-31.37)*cos((yy+60.83)*conv05))**2  &
                              +(yy-60.83)**2).lt.2.*sqrt(18135.)) then
              print * , '**** ADUJUSTING LAKE ladoga LEVEL ****'
              print * , '     NEW:' , zladoga , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zladoga
              if(i_lake==1) dhlake(i,j) = 51.0 ! 230.0 
            else if(ds1deg*sqrt( ((xx-104.00)*cos((yy+13.00)*conv05))**2  &
                              +(yy-13.00)**2).lt.2.*sqrt(16000.)) then
              print * , '**** ADUJUSTING LAKE sap LEVEL ****'
              print * , '     NEW:' , zsap , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zsap
              if(i_lake==1) dhlake(i,j) = 6.0 ! 12.0 
            else if(ds1deg*sqrt( ((xx-71.50)*cos((yy+9.67)*conv05))**2  &
                              +(yy-9.67)**2).lt.2.*sqrt(13010.)) then
              print * , '**** ADUJUSTING LAKE maracaibo LEVEL ****'
              print * , '     NEW:' , zmaracaibo , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zmaracaibo
              if(i_lake==1) dhlake(i,j) = 30.0 ! 60.0 
            else if(ds1deg*sqrt( ((xx-112.75)*cos((yy+29.30)*conv05))**2  &
                              +(yy-29.30)**2).lt.2.*sqrt(12000.)) then
              print * , '**** ADUJUSTING LAKE tungting LEVEL ****'
              print * , '     NEW:' , ztungting , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = ztungting
              if(i_lake==1) dhlake(i,j) = 5.0 ! 10.0 
            else if(ds1deg*sqrt( ((xx+51.25)*cos((yy-31.10)*conv05))**2  &
                              +(yy+31.10)**2).lt.2.*sqrt(10140.)) then
              print * , '**** ADUJUSTING LAKE patos LEVEL ****'
              print * , '     NEW:' , zpatos , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zpatos
              if(i_lake==1) dhlake(i,j) = 5.0 ! ????? 
      
            else if(ds1deg*sqrt( ((xx-35.37)*cos((yy+61.92)*conv05))**2  &
                              +(yy-61.92)**2).lt.2.*sqrt(9890.)) then
              print * , '**** ADUJUSTING LAKE onega LEVEL ****'
              print * , '     NEW:' , zonega , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zonega
              if(i_lake==1) dhlake(i,j) = 30.0 ! 120.0 
            else if(ds1deg*sqrt( ((xx-137.33)*cos((yy-28.50)*conv05))**2  &
                              +(yy+28.50)**2).lt.2.*sqrt(9690.)) then
              print * , '**** ADUJUSTING LAKE eyre LEVEL ****'
              print * , '     NEW:' , zeyre , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zeyre
              if(i_lake==1) dhlake(i,j) = 3.1 ! ????? 
            else if(ds1deg*sqrt( ((xx-1.00)*cos((yy+7.70)*conv05))**2  &
                              +(yy-7.70)**2).lt.2.*sqrt(8502.)) then
              print * , '**** ADUJUSTING LAKE volta LEVEL ****'
              print * , '     NEW:' , zvolta , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zvolta
              if(i_lake==1) dhlake(i,j) = 18.8 ! 75.0 
            else if(ds1deg*sqrt( ((xx+69.57)*cos((yy-15.62)*conv05))**2  &
                              +(yy+15.62)**2).lt.2.*sqrt(8372.)) then
              print * , '**** ADUJUSTING LAKE titicaca LEVEL ****'
              print * , '     NEW:' , ztiticaca , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = ztiticaca
              if(i_lake==1) dhlake(i,j) = 107.0 ! 281.0 
            else if(ds1deg*sqrt( ((xx-18.33)*cos((yy-2.00)*conv05))**2  &
                              +(yy+2.00)**2).lt.2.*sqrt(8210.)) then
              print * , '**** ADUJUSTING LAKE leopold2 LEVEL ****'
              print * , '     NEW:' , zleopold2 , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zleopold2
              if(i_lake==1) dhlake(i,j) = 5.0 ! 12.0 
            else if(ds1deg*sqrt( ((xx+85.50)*cos((yy+11.50)*conv05))**2  &
                              +(yy-11.50)**2).lt.2.*sqrt(8150.)) then
              print * , '**** ADUJUSTING LAKE nicaragua LEVEL ****'
              print * , '     NEW:' , znicaragua , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = znicaragua
              if(i_lake==1) dhlake(i,j) = 35.0 ! 70.0 
            else if(ds1deg*sqrt( ((xx+108.67)*cos((yy+59.12)*conv05))**2  &
                              +(yy-59.12)**2).lt.2.*sqrt(7900.)) then
              print * , '**** ADUJUSTING LAKE athabasca LEVEL ****'
              print * , '     NEW:' , zathabasca , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zathabasca
              if(i_lake==1) dhlake(i,j) = 26.0 ! 120.0 
            else if(ds1deg*sqrt( ((xx-36.28)*cos((yy+3.48)*conv05))**2  &
                              +(yy-3.48)**2).lt.2.*sqrt(6750.)) then
              print * , '**** ADUJUSTING LAKE turkana LEVEL ****'
              print * , '     NEW:' , zturkana , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zturkana
              if(i_lake==1) dhlake(i,j) = 30.2 ! 109.0 
            else if(ds1deg*sqrt( ((xx+64.00)*cos((yy+54.00)*conv05))**2  &
                              +(yy-54.00)**2).lt.2.*sqrt(6475.)) then
              print * , '**** ADUJUSTING LAKE smallwood LEVEL ****'
              print * , '     NEW:' , zsmallwood , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zsmallwood
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx-77.25)*cos((yy+42.92)*conv05))**2  &
                              +(yy-42.92)**2).lt.2.*sqrt(6236.)) then
              print * , '**** ADUJUSTING LAKE issyk_kool LEVEL ****'
              print * , '     NEW:' , zissyk_kool , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zissyk_kool
              if(i_lake==1) dhlake(i,j) = 270.0 ! 668.0 
            else if(ds1deg*sqrt( ((xx-31.68)*cos((yy+22.20)*conv05))**2  &
                              +(yy-22.20)**2).lt.2.*sqrt(6000.)) then
              print * , '**** ADUJUSTING LAKE aswan_dam LEVEL ****'
              print * , '     NEW:' , zaswan_dam , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zaswan_dam
              if(i_lake==1) dhlake(i,j) = 110.0 ! 162.0 
            else if(ds1deg*sqrt( ((xx-45.50)*cos((yy+37.67)*conv05))**2  &
                              +(yy-37.67)**2).lt.2.*sqrt(5960.)) then
              print * , '**** ADUJUSTING LAKE urmia LEVEL ****'
              print * , '     NEW:' , zurmia , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zurmia
              if(i_lake==1) dhlake(i,j) = 5.0 ! 16.0 
            else if(ds1deg*sqrt( ((xx-49.50)*cos((yy+54.75)*conv05))**2  &
                              +(yy-54.75)**2).lt.2.*sqrt(5900.)) then
              print * , '**** ADUJUSTING LAKE kujbyshevsk LEVEL ****'
              print * , '     NEW:' , zkujbyshevsk , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zkujbyshevsk
              if(i_lake==1) dhlake(i,j) = 9.8 ! 41.0 
            else if(ds1deg*sqrt( ((xx-137.83)*cos((yy-31.00)*conv05))**2  &
                              +(yy+31.00)**2).lt.2.*sqrt(5780.)) then
              print * , '**** ADUJUSTING LAKE torrens LEVEL ****'
              print * , '     NEW:' , ztorrens , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = ztorrens
              if(i_lake==1) dhlake(i,j) = 0.5 ! 1.0 
            else if(ds1deg*sqrt( ((xx+102.62)*cos((yy+57.33)*conv05))**2  &
                              +(yy-57.33)**2).lt.2.*sqrt(5650.)) then
              print * , '**** ADUJUSTING LAKE reindeer LEVEL ****'
              print * , '     NEW:' , zreindeer , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zreindeer
              if(i_lake==1) dhlake(i,j) = 17.0 ! 219.0 
            else if(ds1deg*sqrt( ((xx-13.23)*cos((yy+58.88)*conv05))**2  &
                              +(yy-58.88)**2).lt.2.*sqrt(5648.)) then
              print * , '**** ADUJUSTING LAKE vanern LEVEL ****'
              print * , '     NEW:' , zvanern , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zvanern
              if(i_lake==1) dhlake(i,j) = 27.0 ! 106.0 
            else if(ds1deg*sqrt( ((xx-83.83)*cos((yy+48.00)*conv05))**2  &
                              +(yy-48.00)**2).lt.2.*sqrt(5510.)) then
              print * , '**** ADUJUSTING LAKE bukhtarmins LEVEL ****'
              print * , '     NEW:' , zbukhtarmins , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zbukhtarmins
              if(i_lake==1) dhlake(i,j) = 5.0 ! 10.0 
            else if(ds1deg*sqrt( ((xx-102.20)*cos((yy+54.50)*conv05))**2  &
                              +(yy-54.50)**2).lt.2.*sqrt(5478.)) then
              print * , '**** ADUJUSTING LAKE bratskoye LEVEL ****'
              print * , '     NEW:' , zbratskoye , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zbratskoye
              if(i_lake==1) dhlake(i,j) = 31.1 ! 150.0 
            else if(ds1deg*sqrt( ((xx-27.87)*cos((yy-17.27)*conv05))**2  &
                              +(yy+17.27)**2).lt.2.*sqrt(5400.)) then
              print * , '**** ADUJUSTING LAKE kariba LEVEL ****'
              print * , '     NEW:' , zkariba , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zkariba
              if(i_lake==1) dhlake(i,j) = 31.0 ! 78.0 
            else if(ds1deg*sqrt( ((xx-30.92)*cos((yy+1.67)*conv05))**2  &
                              +(yy-1.67)**2).lt.2.*sqrt(5300.)) then
              print * , '**** ADUJUSTING LAKE albert LEVEL ****'
              print * , '     NEW:' , zalbert , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zalbert
              if(i_lake==1) dhlake(i,j) = 30.0 ! 58.0 
            else if(ds1deg*sqrt( ((xx-32.12)*cos((yy+22.98)*conv05))**2  &
                              +(yy-22.98)**2).lt.2.*sqrt(5248.)) then
              print * , '**** ADUJUSTING LAKE nasser LEVEL ****'
              print * , '     NEW:' , znasser , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = znasser
              if(i_lake==1) dhlake(i,j) = 25.2 ! 130.0 
            else if(ds1deg*sqrt( ((xx+100.58)*cos((yy+52.60)*conv05))**2  &
                              +(yy-52.60)**2).lt.2.*sqrt(5150.)) then
              print * , '**** ADUJUSTING LAKE winnipegosi LEVEL ****'
              print * , '     NEW:' , zwinnipegosi , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zwinnipegosi
              if(i_lake==1) dhlake(i,j) = 4.24 ! 18.3 
            else if(ds1deg*sqrt( ((xx+70.17)*cos((yy+66.53)*conv05))**2  &
                              +(yy-66.53)**2).lt.2.*sqrt(5050.)) then
              print * , '**** ADUJUSTING LAKE netilling LEVEL ****'
              print * , '     NEW:' , znetilling , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = znetilling
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx+112.67)*cos((yy+41.17)*conv05))**2  &
                              +(yy-41.17)**2).lt.2.*sqrt(5000.)) then
              print * , '**** ADUJUSTING LAKE greatsalt LEVEL ****'
              print * , '     NEW:' , zgreatsalt , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zgreatsalt
              if(i_lake==1) dhlake(i,j) = 6.0 ! 12.0 
            else if(ds1deg*sqrt( ((xx+98.77)*cos((yy+50.98)*conv05))**2  &
                              +(yy-50.98)**2).lt.2.*sqrt(4610.)) then
              print * , '**** ADUJUSTING LAKE manitoba LEVEL ****'
              print * , '     NEW:' , zmanitoba , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zmanitoba
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx-100.38)*cos((yy+36.67)*conv05))**2  &
                              +(yy-36.67)**2).lt.2.*sqrt(4583.)) then
              print * , '**** ADUJUSTING LAKE chienghai LEVEL ****'
              print * , '     NEW:' , zchienghai , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zchienghai
              if(i_lake==1) dhlake(i,j) = 18.6 ! ????? 
            else if(ds1deg*sqrt( ((xx-103.00)*cos((yy+74.50)*conv05))**2  &
                              +(yy-74.50)**2).lt.2.*sqrt(4560.)) then
              print * , '**** ADUJUSTING LAKE taimyr LEVEL ****'
              print * , '     NEW:' , ztaimyr , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = ztaimyr
              if(i_lake==1) dhlake(i,j) = 2.8 ! 26.0 
            else if(ds1deg*sqrt( ((xx-38.83)*cos((yy+58.05)*conv05))**2  &
                              +(yy-58.05)**2).lt.2.*sqrt(4550.)) then
              print * , '**** ADUJUSTING LAKE rybinsk LEVEL ****'
              print * , '     NEW:' , zrybinsk , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zrybinsk
              if(i_lake==1) dhlake(i,j) = 5.6 ! 28.0 
            else if(ds1deg*sqrt( ((xx+88.62)*cos((yy+49.92)*conv05))**2  &
                              +(yy-49.92)**2).lt.2.*sqrt(4510.)) then
              print * , '**** ADUJUSTING LAKE nipigon LEVEL ****'
              print * , '     NEW:' , znipigon , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = znipigon
              if(i_lake==1) dhlake(i,j) = 54.9 ! 165.0 
            else if(ds1deg*sqrt( ((xx-136.00)*cos((yy-31.58)*conv05))**2  &
                              +(yy+31.58)**2).lt.2.*sqrt(4470.)) then
              print * , '**** ADUJUSTING LAKE gairdner LEVEL ****'
              print * , '     NEW:' , zgairdner , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zgairdner
              if(i_lake==1) dhlake(i,j) = 0.5 ! 1.0 
            else if(ds1deg*sqrt( ((xx-28.75)*cos((yy-9.00)*conv05))**2  &
                              +(yy+9.00)**2).lt.2.*sqrt(4350.)) then
              print * , '**** ADUJUSTING LAKE mweru LEVEL ****'
              print * , '     NEW:' , zmweru , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zmweru
              if(i_lake==1) dhlake(i,j) = 7.0 ! 37.0 
            else if(ds1deg*sqrt( ((xx+94.65)*cos((yy+49.25)*conv05))**2  &
                              +(yy-49.25)**2).lt.2.*sqrt(4350.)) then
              print * , '**** ADUJUSTING LAKE woods LEVEL ****'
              print * , '     NEW:' , zwoods , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zwoods
              if(i_lake==1) dhlake(i,j) = 11.0 ! 21.0 
            else if(ds1deg*sqrt( ((xx-30.87)*cos((yy+57.32)*conv05))**2  &
                              +(yy-57.32)**2).lt.2.*sqrt(4300.)) then
              print * , '**** ADUJUSTING LAKE peipus LEVEL ****'
              print * , '     NEW:' , zpeipus , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zpeipus
              if(i_lake==1) dhlake(i,j) = 15.0 ! 25.0 
            else if(ds1deg*sqrt( ((xx+40.83)*cos((yy-9.58)*conv05))**2  &
                              +(yy+9.58)**2).lt.2.*sqrt(4220.)) then
              print * , '**** ADUJUSTING LAKE sobradinho LEVEL ****'
              print * , '     NEW:' , zsobradinho , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zsobradinho
              if(i_lake==1) dhlake(i,j) = 8.6 ! 30.0 
            else if(ds1deg*sqrt( ((xx-132.42)*cos((yy+44.92)*conv05))**2  &
                              +(yy-44.92)**2).lt.2.*sqrt(4190.)) then
              print * , '**** ADUJUSTING LAKE khanka LEVEL ****'
              print * , '     NEW:' , zkhanka , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zkhanka
              if(i_lake==1) dhlake(i,j) = 4.5 ! 6.5 
            else if(ds1deg*sqrt( ((xx-43.40)*cos((yy+38.65)*conv05))**2  &
                              +(yy-38.65)**2).lt.2.*sqrt(3713.)) then
              print * , '**** ADUJUSTING LAKE van LEVEL ****'
              print * , '     NEW:' , zvan , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zvan
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx+101.53)*cos((yy+63.10)*conv05))**2  &
                              +(yy-63.10)**2).lt.2.*sqrt(3630.)) then
              print * , '**** ADUJUSTING LAKE dubawnt LEVEL ****'
              print * , '     NEW:' , zdubawnt , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zdubawnt
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx-37.38)*cos((yy+11.60)*conv05))**2  &
                              +(yy-11.60)**2).lt.2.*sqrt(3600.)) then
              print * , '**** ADUJUSTING LAKE tana LEVEL ****'
              print * , '     NEW:' , ztana , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = ztana
              if(i_lake==1) dhlake(i,j) = 9.0 ! 14.0 
            else if(ds1deg*sqrt( ((xx-27.62)*cos((yy+58.50)*conv05))**2  &
                              +(yy-58.50)**2).lt.2.*sqrt(3558.)) then
              print * , '**** ADUJUSTING LAKE chudsko LEVEL ****'
              print * , '     NEW:' , zchudsko , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zchudsko
              if(i_lake==1) dhlake(i,j) = 7.1 ! 15.3 
            else if(ds1deg*sqrt( ((xx-92.82)*cos((yy+50.33)*conv05))**2  &
                              +(yy-50.33)**2).lt.2.*sqrt(3350.)) then
              print * , '**** ADUJUSTING LAKE uvs LEVEL ****'
              print * , '     NEW:' , zuvs , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zuvs
              if(i_lake==1) dhlake(i,j) = 6.0 ! 10.0 
            else if(ds1deg*sqrt( ((xx-116.28)*cos((yy+29.08)*conv05))**2  &
                              +(yy-29.08)**2).lt.2.*sqrt(3210.)) then
              print * , '**** ADUJUSTING LAKE poyang LEVEL ****'
              print * , '     NEW:' , zpoyang , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zpoyang
              if(i_lake==1) dhlake(i,j) = 8.4 ! 25.0 
            else if(ds1deg*sqrt( ((xx-45.50)*cos((yy+50.33)*conv05))**2  &
                              +(yy-50.33)**2).lt.2.*sqrt(3120.)) then
              print * , '**** ADUJUSTING LAKE volgogradsk LEVEL ****'
              print * , '     NEW:' , zvolgogradsk , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zvolgogradsk
              if(i_lake==1) dhlake(i,j) = 10.0 ! 41.0 
            else if(ds1deg*sqrt( ((xx+59.43)*cos((yy+53.75)*conv05))**2  &
                              +(yy-53.75)**2).lt.2.*sqrt(3069.)) then
              print * , '**** ADUJUSTING LAKE melville LEVEL ****'
              print * , '     NEW:' , zmelville , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zmelville
              if(i_lake==1) dhlake(i,j) = 130.0 ! 256.0 
            else if(ds1deg*sqrt( ((xx+71.18)*cos((yy+64.98)*conv05))**2  &
                              +(yy-64.98)**2).lt.2.*sqrt(3060.)) then
              print * , '**** ADUJUSTING LAKE amadjuak LEVEL ****'
              print * , '     NEW:' , zamadjuak , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zamadjuak
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx-90.50)*cos((yy+40.50)*conv05))**2  &
                              +(yy-40.50)**2).lt.2.*sqrt(3010.)) then
              print * , '**** ADUJUSTING LAKE lop_nor LEVEL ****'
              print * , '     NEW:' , zlop_nor , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zlop_nor
              if(i_lake==1) dhlake(i,j) = 1.0 ! 2.0 
            else if(ds1deg*sqrt( ((xx-61.17)*cos((yy+31.00)*conv05))**2  &
                              +(yy-31.00)**2).lt.2.*sqrt(3000.)) then
              print * , '**** ADUJUSTING LAKE helmand LEVEL ****'
              print * , '     NEW:' , zhelmand , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zhelmand
              if(i_lake==1) dhlake(i,j) = 6.0 ! 11.0 
            else if(ds1deg*sqrt( ((xx-32.42)*cos((yy-8.00)*conv05))**2  &
                              +(yy+8.00)**2).lt.2.*sqrt(3000.)) then
              print * , '**** ADUJUSTING LAKE rukwa LEVEL ****'
              print * , '     NEW:' , zrukwa , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zrukwa
              if(i_lake==1) dhlake(i,j) = 0.5 ! 1.0 
            else if(ds1deg*sqrt( ((xx+52.83)*cos((yy-32.75)*conv05))**2  &
                              +(yy+32.75)**2).lt.2.*sqrt(2970.)) then
              print * , '**** ADUJUSTING LAKE mirim LEVEL ****'
              print * , '     NEW:' , zmirim , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zmirim
              if(i_lake==1) dhlake(i,j) = 5.0 ! 10.0 
            else if(ds1deg*sqrt( ((xx+69.37)*cos((yy+54.28)*conv05))**2  &
                              +(yy-54.28)**2).lt.2.*sqrt(2892.5)) then
              print * , '**** ADUJUSTING LAKE caniapiscau LEVEL ****'
              print * , '     NEW:' , zcaniapiscau , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zcaniapiscau
              if(i_lake==1) dhlake(i,j) = 11.84 ! 49.0 
            else if(ds1deg*sqrt( ((xx-101.32)*cos((yy+51.03)*conv05))**2  &
                              +(yy-51.03)**2).lt.2.*sqrt(2770.)) then
              print * , '**** ADUJUSTING LAKE hovsgol LEVEL ****'
              print * , '     NEW:' , zhovsgol , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zhovsgol
              if(i_lake==1) dhlake(i,j) = 138.0 ! 267.0 
            else if(ds1deg*sqrt( ((xx-112.42)*cos((yy+29.42)*conv05))**2  &
                              +(yy-29.42)**2).lt.2.*sqrt(2740.)) then
              print * , '**** ADUJUSTING LAKE dongting LEVEL ****'
              print * , '     NEW:' , zdongting , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zdongting
              if(i_lake==1) dhlake(i,j) = 6.7 ! 30.8 
            else if(ds1deg*sqrt( ((xx-31.57)*cos((yy-15.73)*conv05))**2  &
                              +(yy+15.73)**2).lt.2.*sqrt(2739.)) then
              print * , '**** ADUJUSTING LAKE caborabassa LEVEL ****'
              print * , '     NEW:' , zcaborabassa , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zcaborabassa
              if(i_lake==1) dhlake(i,j) = 20.9 ! 157.0 
            else if(ds1deg*sqrt( ((xx+43.08)*cos((yy-0.80)*conv05))**2  &
                              +(yy+0.80)**2).lt.2.*sqrt(2703.)) then
              print * , '**** ADUJUSTING LAKE tsimlyansko LEVEL ****'
              print * , '     NEW:' , ztsimlyansko , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = ztsimlyansko
              if(i_lake==1) dhlake(i,j) = 8.8 ! 35.0 
            else if(ds1deg*sqrt( ((xx+103.47)*cos((yy+58.28)*conv05))**2  &
                              +(yy-58.28)**2).lt.2.*sqrt(2681.)) then
              print * , '**** ADUJUSTING LAKE wollaston LEVEL ****'
              print * , '     NEW:' , zwollaston , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zwollaston
              if(i_lake==1) dhlake(i,j) = 20.6 ! 97.0 
            else if(ds1deg*sqrt( ((xx-81.75)*cos((yy+46.08)*conv05))**2  &
                              +(yy-46.08)**2).lt.2.*sqrt(2650.)) then
              print * , '**** ADUJUSTING LAKE alakol LEVEL ****'
              print * , '     NEW:' , zalakol , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zalakol
              if(i_lake==1) dhlake(i,j) = 22.0 ! 54.0 
            else if(ds1deg*sqrt( ((xx+155.00)*cos((yy+59.50)*conv05))**2  &
                              +(yy-59.50)**2).lt.2.*sqrt(2590.)) then
              print * , '**** ADUJUSTING LAKE iliamna LEVEL ****'
              print * , '     NEW:' , ziliamna , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = ziliamna
              if(i_lake==1) dhlake(i,j) = 150.0 ! 299.0 
            else if(ds1deg*sqrt( ((xx-0.35)*cos((yy+5.55)*conv05))**2  &
                              +(yy-5.55)**2).lt.2.*sqrt(2564.)) then
              print * , '**** ADUJUSTING LAKE densu LEVEL ****'
              print * , '     NEW:' , zdensu , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zdensu
              if(i_lake==1) dhlake(i,j) = 10.5 ! 15.6 
            else if(ds1deg*sqrt( ((xx-90.50)*cos((yy+30.75)*conv05))**2  &
                              +(yy-30.75)**2).lt.2.*sqrt(2500.)) then
              print * , '**** ADUJUSTING LAKE nam LEVEL ****'
              print * , '     NEW:' , znam , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = znam
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx+76.95)*cos((yy+53.75)*conv05))**2  &
                              +(yy-53.75)**2).lt.2.*sqrt(2485.5)) then
              print * , '**** ADUJUSTING LAKE lagrande2 LEVEL ****'
              print * , '     NEW:' , zlagrande2 , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zlagrande2
              if(i_lake==1) dhlake(i,j) = 21.8 ! 137.0 
            else if(ds1deg*sqrt( ((xx+49.75)*cos((yy-4.75)*conv05))**2  &
                              +(yy+4.75)**2).lt.2.*sqrt(2430.)) then
              print * , '**** ADUJUSTING LAKE tucurui LEVEL ****'
              print * , '     NEW:' , ztucurui , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = ztucurui
              if(i_lake==1) dhlake(i,j) = 18.9 ! 75.0 
            else if(ds1deg*sqrt( ((xx-120.25)*cos((yy+31.25)*conv05))**2  &
                              +(yy-31.25)**2).lt.2.*sqrt(2427.8)) then
              print * , '**** ADUJUSTING LAKE tai LEVEL ****'
              print * , '     NEW:' , ztai , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = ztai
              if(i_lake==1) dhlake(i,j) = 1.9 ! 2.6 
            else if(ds1deg*sqrt( ((xx+75.00)*cos((yy+53.83)*conv05))**2  &
                              +(yy-53.83)**2).lt.2.*sqrt(2420.)) then
              print * , '**** ADUJUSTING LAKE lagrande3 LEVEL ****'
              print * , '     NEW:' , zlagrande3 , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zlagrande3
              if(i_lake==1) dhlake(i,j) = 24.4 ! 90.0 
            else if(ds1deg*sqrt( ((xx+98.82)*cos((yy+57.13)*conv05))**2  &
                              +(yy-57.13)**2).lt.2.*sqrt(2391.)) then
              print * , '**** ADUJUSTING LAKE southernind LEVEL ****'
              print * , '     NEW:' , zsouthernind , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zsouthernind
              if(i_lake==1) dhlake(i,j) = 9.8 ! 30.0 
            else if(ds1deg*sqrt( ((xx+59.47)*cos((yy-1.92)*conv05))**2  &
                              +(yy+1.92)**2).lt.2.*sqrt(2360.)) then
              print * , '**** ADUJUSTING LAKE balbina LEVEL ****'
              print * , '     NEW:' , zbalbina , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zbalbina
              if(i_lake==1) dhlake(i,j) = 7.4 ! 30.0 
            else if(ds1deg*sqrt( ((xx-29.57)*cos((yy-0.88)*conv05))**2  &
                              +(yy+0.88)**2).lt.2.*sqrt(2325.)) then
              print * , '**** ADUJUSTING LAKE edward LEVEL ****'
              print * , '     NEW:' , zedward , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zedward
              if(i_lake==1) dhlake(i,j) = 17.0 ! 112.0 
            else if(ds1deg*sqrt( ((xx-32.00)*cos((yy+49.33)*conv05))**2  &
                              +(yy-49.33)**2).lt.2.*sqrt(2250.)) then
              print * , '**** ADUJUSTING LAKE krementchug LEVEL ****'
              print * , '     NEW:' , zkrementchug , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zkrementchug
              if(i_lake==1) dhlake(i,j) = 6.0 ! 20.0 
            else if(ds1deg*sqrt( ((xx+52.97)*cos((yy-22.47)*conv05))**2  &
                              +(yy+22.47)**2).lt.2.*sqrt(2250.)) then
              print * , '**** ADUJUSTING LAKE primavera LEVEL ****'
              print * , '     NEW:' , zprimavera , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zprimavera
              if(i_lake==1) dhlake(i,j) = 10.0 ! 18.3 
            else if(ds1deg*sqrt( ((xx+72.00)*cos((yy-46.50)*conv05))**2  &
                              +(yy+46.50)**2).lt.2.*sqrt(2240.)) then
              print * , '**** ADUJUSTING LAKE buenosaires LEVEL ****'
              print * , '     NEW:' , zbuenosaires , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zbuenosaires
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx-29.43)*cos((yy-2.00)*conv05))**2  &
                              +(yy+2.00)**2).lt.2.*sqrt(2220.)) then
              print * , '**** ADUJUSTING LAKE kivu LEVEL ****'
              print * , '     NEW:' , zkivu , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zkivu
              if(i_lake==1) dhlake(i,j) = 240.0 ! 480.0 
            else if(ds1deg*sqrt( ((xx-111.00)*cos((yy+62.78)*conv05))**2  &
                              +(yy-62.78)**2).lt.2.*sqrt(2170.)) then
              print * , '**** ADUJUSTING LAKE viluyskoe LEVEL ****'
              print * , '     NEW:' , zviluyskoe , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zviluyskoe
              if(i_lake==1) dhlake(i,j) = 30.0 ! 60.0 
            else if(ds1deg*sqrt( ((xx-34.33)*cos((yy+47.20)*conv05))**2  &
                              +(yy-47.20)**2).lt.2.*sqrt(2150.)) then
              print * , '**** ADUJUSTING LAKE kakhovskoye LEVEL ****'
              print * , '     NEW:' , zkakhovskoye , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zkakhovskoye
              if(i_lake==1) dhlake(i,j) = 8.5 ! 24.0 
            else if(ds1deg*sqrt( ((xx-128.12)*cos((yy+54.28)*conv05))**2  &
                              +(yy-54.28)**2).lt.2.*sqrt(2119.)) then
              print * , '**** ADUJUSTING LAKE zeiskoye LEVEL ****'
              print * , '     NEW:' , zzeiskoye , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zzeiskoye
              if(i_lake==1) dhlake(i,j) = 50.0 ! 93.0 
            else if(ds1deg*sqrt( ((xx+73.42)*cos((yy+50.88)*conv05))**2  &
                              +(yy-50.88)**2).lt.2.*sqrt(2115.)) then
              print * , '**** ADUJUSTING LAKE mistassini LEVEL ****'
              print * , '     NEW:' , zmistassini , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zmistassini
              if(i_lake==1) dhlake(i,j) = 75.0 ! 183.0 
            else if(ds1deg*sqrt( ((xx+64.10)*cos((yy+54.12)*conv05))**2  &
                              +(yy-54.12)**2).lt.2.*sqrt(2030.)) then
              print * , '**** ADUJUSTING LAKE michikamau LEVEL ****'
              print * , '     NEW:' , zmichikamau , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zmichikamau
              if(i_lake==1) dhlake(i,j) = 40.0 ! 80.0 
            else if(ds1deg*sqrt( ((xx+99.97)*cos((yy+60.20)*conv05))**2  &
                              +(yy-60.20)**2).lt.2.*sqrt(2030.)) then
              print * , '**** ADUJUSTING LAKE nueltin LEVEL ****'
              print * , '     NEW:' , znueltin , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = znueltin
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx-77.50)*cos((yy+54.82)*conv05))**2  &
                              +(yy-54.82)**2).lt.2.*sqrt(2010.)) then
              print * , '**** ADUJUSTING LAKE chany LEVEL ****'
              print * , '     NEW:' , zchany , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zchany
              if(i_lake==1) dhlake(i,j) = 2.2 ! 5.5 
            else if(ds1deg*sqrt( ((xx-91.70)*cos((yy+54.83)*conv05))**2  &
                              +(yy-54.83)**2).lt.2.*sqrt(2000.)) then
              print * , '**** ADUJUSTING LAKE krasnoyarsk LEVEL ****'
              print * , '     NEW:' , zkrasnoyarsk , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zkrasnoyarsk
              if(i_lake==1) dhlake(i,j) = 37.0 ! 105.0 
            else if(ds1deg*sqrt( ((xx+62.00)*cos((yy-30.70)*conv05))**2  &
                              +(yy+30.70)**2).lt.2.*sqrt(1984.)) then
              print * , '**** ADUJUSTING LAKE marchiquita LEVEL ****'
              print * , '     NEW:' , zmarchiquita , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zmarchiquita
              if(i_lake==1) dhlake(i,j) = 7.3 ! ????? 
            else if(ds1deg*sqrt( ((xx-118.67)*cos((yy+33.33)*conv05))**2  &
                              +(yy-33.33)**2).lt.2.*sqrt(1960.)) then
              print * , '**** ADUJUSTING LAKE hungtze LEVEL ****'
              print * , '     NEW:' , zhungtze , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zhungtze
              if(i_lake==1) dhlake(i,j) = 1.4 ! ????? 
            else if(ds1deg*sqrt( ((xx+68.72)*cos((yy+50.63)*conv05))**2  &
                              +(yy-50.63)**2).lt.2.*sqrt(1950.)) then
              print * , '**** ADUJUSTING LAKE manicouagan LEVEL ****'
              print * , '     NEW:' , zmanicouagan , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zmanicouagan
              if(i_lake==1) dhlake(i,j) = 85.0 ! 350.0 
            else if(ds1deg*sqrt( ((xx-47.17)*cos((yy+30.83)*conv05))**2  &
                              +(yy-30.83)**2).lt.2.*sqrt(1940.)) then
              print * , '**** ADUJUSTING LAKE hammer LEVEL ****'
              print * , '     NEW:' , zhammer , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zhammer
              if(i_lake==1) dhlake(i,j) = 1.0 ! 2.0 
            else if(ds1deg*sqrt( ((xx-90.50)*cos((yy+30.67)*conv05))**2  &
                              +(yy-30.67)**2).lt.2.*sqrt(1920.)) then
              print * , '**** ADUJUSTING LAKE namu LEVEL ****'
              print * , '     NEW:' , znamu , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = znamu
              if(i_lake==1) dhlake(i,j) = 26.0 ! ????? 
            else if(ds1deg*sqrt( ((xx-102.70)*cos((yy+57.20)*conv05))**2  &
                              +(yy-57.20)**2).lt.2.*sqrt(1920.)) then
              print * , '**** ADUJUSTING LAKE ust_ilimsko LEVEL ****'
              print * , '     NEW:' , zust_ilimsko , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zust_ilimsko
              if(i_lake==1) dhlake(i,j) = 30.7 ! 97.0 
            else if(ds1deg*sqrt( ((xx-56.25)*cos((yy+58.53)*conv05))**2  &
                              +(yy-58.53)**2).lt.2.*sqrt(1910.)) then
              print * , '**** ADUJUSTING LAKE kamskoye LEVEL ****'
              print * , '     NEW:' , zkamskoye , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zkamskoye
              if(i_lake==1) dhlake(i,j) = 6.4 ! 28.6 
            else if(ds1deg*sqrt( ((xx+80.87)*cos((yy+26.93)*conv05))**2  &
                              +(yy-26.93)**2).lt.2.*sqrt(1894.)) then
              print * , '**** ADUJUSTING LAKE okeechobee LEVEL ****'
              print * , '     NEW:' , zokeechobee , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zokeechobee
              if(i_lake==1) dhlake(i,j) = 2.0 ! 4.5 
            else if(ds1deg*sqrt( ((xx-14.55)*cos((yy+58.32)*conv05))**2  &
                              +(yy-58.32)**2).lt.2.*sqrt(1856.)) then
              print * , '**** ADUJUSTING LAKE vattern LEVEL ****'
              print * , '     NEW:' , zvattern , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zvattern
              if(i_lake==1) dhlake(i,j) = 39.9 ! 128.0 
            else if(ds1deg*sqrt( ((xx+72.00)*cos((yy-46.58)*conv05))**2  &
                              +(yy+46.58)**2).lt.2.*sqrt(1850.)) then
              print * , '**** ADUJUSTING LAKE carrera LEVEL ****'
              print * , '     NEW:' , zcarrera , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zcarrera
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx-77.95)*cos((yy+43.83)*conv05))**2  &
                              +(yy-43.83)**2).lt.2.*sqrt(1850.)) then
              print * , '**** ADUJUSTING LAKE kaptchagays LEVEL ****'
              print * , '     NEW:' , zkaptchagays , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zkaptchagays
              if(i_lake==1) dhlake(i,j) = 15.2 ! 40.0 
            else if(ds1deg*sqrt( ((xx-48.95)*cos((yy+52.70)*conv05))**2  &
                              +(yy-52.70)**2).lt.2.*sqrt(1830.)) then
              print * , '**** ADUJUSTING LAKE saratovskoy LEVEL ****'
              print * , '     NEW:' , zsaratovskoy , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zsaratovskoy
              if(i_lake==1) dhlake(i,j) = 7.0 ! 28.0 
            else if(ds1deg*sqrt( ((xx+95.27)*cos((yy+64.15)*conv05))**2  &
                              +(yy-64.15)**2).lt.2.*sqrt(1780.)) then
              print * , '**** ADUJUSTING LAKE baker LEVEL ****'
              print * , '     NEW:' , zbaker , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zbaker
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx+123.58)*cos((yy+56.02)*conv05))**2  &
                              +(yy-56.02)**2).lt.2.*sqrt(1779.)) then
              print * , '**** ADUJUSTING LAKE williston LEVEL ****'
              print * , '     NEW:' , zwilliston , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zwilliston
              if(i_lake==1) dhlake(i,j) = 43.3 ! 166.0 
            else if(ds1deg*sqrt( ((xx+117.95)*cos((yy+63.32)*conv05))**2  &
                              +(yy-63.32)**2).lt.2.*sqrt(1776.)) then
              print * , '**** ADUJUSTING LAKE martre LEVEL ****'
              print * , '     NEW:' , zmartre , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zmartre
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx-92.17)*cos((yy+48.00)*conv05))**2  &
                              +(yy-48.00)**2).lt.2.*sqrt(1760.)) then
              print * , '**** ADUJUSTING LAKE har_us LEVEL ****'
              print * , '     NEW:' , zhar_us , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zhar_us
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx-28.25)*cos((yy+61.25)*conv05))**2  &
                              +(yy-61.25)**2).lt.2.*sqrt(1760.)) then
              print * , '**** ADUJUSTING LAKE saimaa LEVEL ****'
              print * , '     NEW:' , zsaimaa , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zsaimaa
              if(i_lake==1) dhlake(i,j) = 40.0 ! 82.0 
            else if(ds1deg*sqrt( ((xx-35.67)*cos((yy-15.33)*conv05))**2  &
                              +(yy+15.33)**2).lt.2.*sqrt(1750.)) then
              print * , '**** ADUJUSTING LAKE chilwa LEVEL ****'
              print * , '     NEW:' , zchilwa , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zchilwa
              if(i_lake==1) dhlake(i,j) = 1.0 ! 2.7 
            else if(ds1deg*sqrt( ((xx-117.37)*cos((yy+48.92)*conv05))**2  &
                              +(yy-48.92)**2).lt.2.*sqrt(1730.5)) then
              print * , '**** ADUJUSTING LAKE hulun LEVEL ****'
              print * , '     NEW:' , zhulun , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zhulun
              if(i_lake==1) dhlake(i,j) = 5.0 ! 7.0 
            else if(ds1deg*sqrt( ((xx-33.25)*cos((yy+1.50)*conv05))**2  &
                              +(yy-1.50)**2).lt.2.*sqrt(1720.)) then
              print * , '**** ADUJUSTING LAKE kyoga LEVEL ****'
              print * , '     NEW:' , zkyoga , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zkyoga
              if(i_lake==1) dhlake(i,j) = 3.0 ! 5.7 
            else if(ds1deg*sqrt( ((xx+91.55)*cos((yy+18.62)*conv05))**2  &
                              +(yy-18.62)**2).lt.2.*sqrt(1660.)) then
              print * , '**** ADUJUSTING LAKE terminos LEVEL ****'
              print * , '     NEW:' , zterminos , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zterminos
              if(i_lake==1) dhlake(i,j) = 0.4 ! 1.0 
            else if(ds1deg*sqrt( ((xx-89.00)*cos((yy+31.83)*conv05))**2  &
                              +(yy-31.83)**2).lt.2.*sqrt(1640.)) then
              print * , '**** ADUJUSTING LAKE qilin LEVEL ****'
              print * , '     NEW:' , zqilin , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zqilin
              if(i_lake==1) dhlake(i,j) = 33.0 ! ????? 
            else if(ds1deg*sqrt( ((xx-21.00)*cos((yy+55.00)*conv05))**2  &
                              +(yy-55.00)**2).lt.2.*sqrt(1620.)) then
              print * , '**** ADUJUSTING LAKE kurisches LEVEL ****'
              print * , '     NEW:' , zkurisches , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zkurisches
              if(i_lake==1) dhlake(i,j) = 4.0 ! 10.0 
            else if(ds1deg*sqrt( ((xx+90.12)*cos((yy+30.22)*conv05))**2  &
                              +(yy-30.22)**2).lt.2.*sqrt(1620.)) then
              print * , '**** ADUJUSTING LAKE pontchartr LEVEL ****'
              print * , '     NEW:' , zpontchartr , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zpontchartr
              if(i_lake==1) dhlake(i,j) = 2.0 ! 5.0 
            else if(ds1deg*sqrt( ((xx+56.80)*cos((yy-27.50)*conv05))**2  &
                              +(yy+27.50)**2).lt.2.*sqrt(1600.)) then
              print * , '**** ADUJUSTING LAKE yacyreta LEVEL ****'
              print * , '     NEW:' , zyacyreta , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zyacyreta
              if(i_lake==1) dhlake(i,j) = 13.1 ! 26.0 
            else if(ds1deg*sqrt( ((xx-41.50)*cos((yy+57.12)*conv05))**2  &
                              +(yy-57.12)**2).lt.2.*sqrt(1590.)) then
              print * , '**** ADUJUSTING LAKE gorkovskoye LEVEL ****'
              print * , '     NEW:' , zgorkovskoye , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zgorkovskoye
              if(i_lake==1) dhlake(i,j) = 5.5 ! 22.0 
            else if(ds1deg*sqrt( ((xx-68.95)*cos((yy+50.40)*conv05))**2  &
                              +(yy-50.40)**2).lt.2.*sqrt(1590.)) then
              print * , '**** ADUJUSTING LAKE tengiz LEVEL ****'
              print * , '     NEW:' , ztengiz , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = ztengiz
              if(i_lake==1) dhlake(i,j) = 7.0 ! 8.0 
            else if(ds1deg*sqrt( ((xx-14.17)*cos((yy+13.33)*conv05))**2  &
                              +(yy-13.33)**2).lt.2.*sqrt(1540.)) then
              print * , '**** ADUJUSTING LAKE chad LEVEL ****'
              print * , '     NEW:' , zchad , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zchad
              if(i_lake==1) dhlake(i,j) = 4.11 ! 10.5 
            else if(ds1deg*sqrt( ((xx-29.75)*cos((yy-11.08)*conv05))**2  &
                              +(yy+11.08)**2).lt.2.*sqrt(1510.)) then
              print * , '**** ADUJUSTING LAKE bangweulu LEVEL ****'
              print * , '     NEW:' , zbangweulu , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zbangweulu
              if(i_lake==1) dhlake(i,j) = 4.0 ! 10.0 
            else if(ds1deg*sqrt( ((xx-35.42)*cos((yy+38.95)*conv05))**2  &
                              +(yy-38.95)**2).lt.2.*sqrt(1500.)) then
              print * , '**** ADUJUSTING LAKE tuz LEVEL ****'
              print * , '     NEW:' , ztuz , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = ztuz
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx+72.75)*cos((yy-50.33)*conv05))**2  &
                              +(yy+50.33)**2).lt.2.*sqrt(1466.)) then
              print * , '**** ADUJUSTING LAKE argentino LEVEL ****'
              print * , '     NEW:' , zargentino , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zargentino
              if(i_lake==1) dhlake(i,j) = 150.0 ! 500.0 
            else if(ds1deg*sqrt( ((xx+92.40)*cos((yy+50.43)*conv05))**2  &
                              +(yy-50.43)**2).lt.2.*sqrt(1450.)) then
              print * , '**** ADUJUSTING LAKE seul LEVEL ****'
              print * , '     NEW:' , zseul , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zseul
              if(i_lake==1) dhlake(i,j) = 7.1 ! 47.2 
            else if(ds1deg*sqrt( ((xx+112.07)*cos((yy+58.58)*conv05))**2  &
                              +(yy-58.58)**2).lt.2.*sqrt(1410.)) then
              print * , '**** ADUJUSTING LAKE claire LEVEL ****'
              print * , '     NEW:' , zclaire , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zclaire
              if(i_lake==1) dhlake(i,j) = 1.2 ! 2.0 
            else if(ds1deg*sqrt( ((xx-93.50)*cos((yy+49.00)*conv05))**2  &
                              +(yy-49.00)**2).lt.2.*sqrt(1407.)) then
              print * , '**** ADUJUSTING LAKE hyargas LEVEL ****'
              print * , '     NEW:' , zhyargas , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zhyargas
              if(i_lake==1) dhlake(i,j) = 20.0 ! 80.0 
            else if(ds1deg*sqrt( ((xx+160.15)*cos((yy+66.50)*conv05))**2  &
                              +(yy-66.50)**2).lt.2.*sqrt(1400.)) then
              print * , '**** ADUJUSTING LAKE selawik LEVEL ****'
              print * , '     NEW:' , zselawik , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zselawik
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx-86.37)*cos((yy+31.00)*conv05))**2  &
                              +(yy-31.00)**2).lt.2.*sqrt(1400.)) then
              print * , '**** ADUJUSTING LAKE tangra LEVEL ****'
              print * , '     NEW:' , ztangra , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = ztangra
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx-87.00)*cos((yy+42.00)*conv05))**2  &
                              +(yy-42.00)**2).lt.2.*sqrt(1380.)) then
              print * , '**** ADUJUSTING LAKE baghrash LEVEL ****'
              print * , '     NEW:' , zbaghrash , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zbaghrash
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx-32.00)*cos((yy+31.25)*conv05))**2  &
                              +(yy-31.25)**2).lt.2.*sqrt(1360.)) then
              print * , '**** ADUJUSTING LAKE manzala LEVEL ****'
              print * , '     NEW:' , zmanzala , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zmanzala
              if(i_lake==1) dhlake(i,j) = 0.5 ! 1.0 
            else if(ds1deg*sqrt( ((xx-45.33)*cos((yy+40.42)*conv05))**2  &
                              +(yy-40.42)**2).lt.2.*sqrt(1360.)) then
              print * , '**** ADUJUSTING LAKE sevan LEVEL ****'
              print * , '     NEW:' , zsevan , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zsevan
              if(i_lake==1) dhlake(i,j) = 41.0 ! 86.0 
            else if(ds1deg*sqrt( ((xx+51.00)*cos((yy-23.00)*conv05))**2  &
                              +(yy+23.00)**2).lt.2.*sqrt(1350.)) then
              print * , '**** ADUJUSTING LAKE itaipu LEVEL ****'
              print * , '     NEW:' , zitaipu , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zitaipu
              if(i_lake==1) dhlake(i,j) = 9.0 ! 18.3 
            else if(ds1deg*sqrt( ((xx+100.15)*cos((yy+53.95)*conv05))**2  &
                              +(yy-53.95)**2).lt.2.*sqrt(1340.)) then
              print * , '**** ADUJUSTING LAKE moose LEVEL ****'
              print * , '     NEW:' , zmoose , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zmoose
              if(i_lake==1) dhlake(i,j) = 3.0 ! 7.0 
            else if(ds1deg*sqrt( ((xx+104.83)*cos((yy+50.58)*conv05))**2  &
                              +(yy-50.58)**2).lt.2.*sqrt(1330.)) then
              print * , '**** ADUJUSTING LAKE ronge LEVEL ****'
              print * , '     NEW:' , zronge , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zronge
              if(i_lake==1) dhlake(i,j) = 14.6 ! 42.1 
            else if(ds1deg*sqrt( ((xx+97.82)*cos((yy+62.58)*conv05))**2  &
                              +(yy-62.58)**2).lt.2.*sqrt(1330.)) then
              print * , '**** ADUJUSTING LAKE yathkyed LEVEL ****'
              print * , '     NEW:' , zyathkyed , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zyathkyed
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx+102.25)*cos((yy+60.55)*conv05))**2  &
                              +(yy-60.55)**2).lt.2.*sqrt(1320.)) then
              print * , '**** ADUJUSTING LAKE kasba LEVEL ****'
              print * , '     NEW:' , zkasba , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zkasba
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx+99.83)*cos((yy+53.35)*conv05))**2  &
                              +(yy-53.35)**2).lt.2.*sqrt(1320.)) then
              print * , '**** ADUJUSTING LAKE cedar LEVEL ****'
              print * , '     NEW:' , zcedar , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zcedar
              if(i_lake==1) dhlake(i,j) = 4.18 ! 10.0 
            else if(ds1deg*sqrt( ((xx+75.00)*cos((yy+48.57)*conv05))**2  &
                              +(yy-48.57)**2).lt.2.*sqrt(1302.)) then
              print * , '**** ADUJUSTING LAKE gouin LEVEL ****'
              print * , '     NEW:' , zgouin , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zgouin
              if(i_lake==1) dhlake(i,j) = 69.14 ! ????? 
            else if(ds1deg*sqrt( ((xx-100.25)*cos((yy+7.50)*conv05))**2  &
                              +(yy-7.50)**2).lt.2.*sqrt(1290.)) then
              print * , '**** ADUJUSTING LAKE luang LEVEL ****'
              print * , '     NEW:' , zluang , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zluang
              if(i_lake==1) dhlake(i,j) = 0.5 ! 1.0 
            else if(ds1deg*sqrt( ((xx-37.58)*cos((yy+60.23)*conv05))**2  &
                              +(yy-60.23)**2).lt.2.*sqrt(1290.)) then
              print * , '**** ADUJUSTING LAKE sheksninsko LEVEL ****'
              print * , '     NEW:' , zsheksninsko , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zsheksninsko
              if(i_lake==1) dhlake(i,j) = 4.5 ! 20.0 
            else if(ds1deg*sqrt( ((xx+74.37)*cos((yy+56.13)*conv05))**2  &
                              +(yy-56.13)**2).lt.2.*sqrt(1250.)) then
              print * , '**** ADUJUSTING LAKE eau_claire LEVEL ****'
              print * , '     NEW:' , zeau_claire , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zeau_claire
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx-34.92)*cos((yy+63.43)*conv05))**2  &
                              +(yy-63.43)**2).lt.2.*sqrt(1250.)) then
              print * , '**** ADUJUSTING LAKE vygozersko LEVEL ****'
              print * , '     NEW:' , zvygozersko , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zvygozersko
              if(i_lake==1) dhlake(i,j) = 7.5 ! 20.0 
            else if(ds1deg*sqrt( ((xx+51.37)*cos((yy-20.37)*conv05))**2  &
                              +(yy+20.37)**2).lt.2.*sqrt(1195.)) then
              print * , '**** ADUJUSTING LAKE ilha_soltei LEVEL ****'
              print * , '     NEW:' , zilha_soltei , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zilha_soltei
              if(i_lake==1) dhlake(i,j) = 25.0 ! 46.0 
            else if(ds1deg*sqrt( ((xx+156.38)*cos((yy+57.93)*conv05))**2  &
                              +(yy-57.93)**2).lt.2.*sqrt(1190.)) then
              print * , '**** ADUJUSTING LAKE becharof LEVEL ****'
              print * , '     NEW:' , zbecharof , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zbecharof
              if(i_lake==1) dhlake(i,j) = 40.0 ! 92.0 
            else if(ds1deg*sqrt( ((xx+115.45)*cos((yy+55.42)*conv05))**2  &
                              +(yy-55.42)**2).lt.2.*sqrt(1170.)) then
              print * , '**** ADUJUSTING LAKE lesser LEVEL ****'
              print * , '     NEW:' , zlesser , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zlesser
              if(i_lake==1) dhlake(i,j) = 11.7 ! 21.0 
            else if(ds1deg*sqrt( ((xx+94.92)*cos((yy+48.03)*conv05))**2  &
                              +(yy-48.03)**2).lt.2.*sqrt(1170.)) then
              print * , '**** ADUJUSTING LAKE red LEVEL ****'
              print * , '     NEW:' , zred , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zred
              if(i_lake==1) dhlake(i,j) = 3.0 ! 9.0 
            else if(ds1deg*sqrt( ((xx-37.63)*cos((yy+6.12)*conv05))**2  &
                              +(yy-6.12)**2).lt.2.*sqrt(1160.)) then
              print * , '**** ADUJUSTING LAKE abaya LEVEL ****'
              print * , '     NEW:' , zabaya , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zabaya
              if(i_lake==1) dhlake(i,j) = 7.0 ! 13.0 
            else if(ds1deg*sqrt( ((xx+106.63)*cos((yy+57.48)*conv05))**2  &
                              +(yy-57.48)**2).lt.2.*sqrt(1152.)) then
              print * , '**** ADUJUSTING LAKE cree LEVEL ****'
              print * , '     NEW:' , zcree , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zcree
              if(i_lake==1) dhlake(i,j) = 8.2 ! 16.8 
            else if(ds1deg*sqrt( ((xx-17.03)*cos((yy+59.52)*conv05))**2  &
                              +(yy-59.52)**2).lt.2.*sqrt(1140.)) then
              print * , '**** ADUJUSTING LAKE malaren LEVEL ****'
              print * , '     NEW:' , zmalaren , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zmalaren
              if(i_lake==1) dhlake(i,j) = 11.9 ! 61.0 
            else if(ds1deg*sqrt( ((xx-73.25)*cos((yy+43.65)*conv05))**2  &
                              +(yy-43.65)**2).lt.2.*sqrt(1130.)) then
              print * , '**** ADUJUSTING LAKE champlain LEVEL ****'
              print * , '     NEW:' , zchamplain , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zchamplain
              if(i_lake==1) dhlake(i,j) = 22.8 ! 123.0 
            else if(ds1deg*sqrt( ((xx-54.78)*cos((yy+57.28)*conv05))**2  &
                              +(yy-57.28)**2).lt.2.*sqrt(1120.)) then
              print * , '**** ADUJUSTING LAKE votkinskoye LEVEL ****'
              print * , '     NEW:' , zvotkinskoye , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zvotkinskoye
              if(i_lake==1) dhlake(i,j) = 8.4 ! 28.0 
            else if(ds1deg*sqrt( ((xx+82.63)*cos((yy+42.98)*conv05))**2  &
                              +(yy-42.98)**2).lt.2.*sqrt(1114.)) then
              print * , '**** ADUJUSTING LAKE stclain LEVEL ****'
              print * , '     NEW:' , zstclain , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zstclain
              if(i_lake==1) dhlake(i,j) = 3.6 ! 9.0 
            else if(ds1deg*sqrt( ((xx+103.05)*cos((yy+20.20)*conv05))**2  &
                              +(yy-20.20)**2).lt.2.*sqrt(1112.)) then
              print * , '**** ADUJUSTING LAKE chapala LEVEL ****'
              print * , '     NEW:' , zchapala , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zchapala
              if(i_lake==1) dhlake(i,j) = 7.2 ! 10.5 
            else if(ds1deg*sqrt( ((xx+83.92)*cos((yy+15.38)*conv05))**2  &
                              +(yy-15.38)**2).lt.2.*sqrt(1110.)) then
              print * , '**** ADUJUSTING LAKE caratasca LEVEL ****'
              print * , '     NEW:' , zcaratasca , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zcaratasca
              if(i_lake==1) dhlake(i,j) = 1.5 ! 5.0 
            else if(ds1deg*sqrt( ((xx+98.95)*cos((yy+64.57)*conv05))**2  &
                              +(yy-64.57)**2).lt.2.*sqrt(1100.)) then
              print * , '**** ADUJUSTING LAKE aberdeen LEVEL ****'
              print * , '     NEW:' , zaberdeen , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zaberdeen
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx+80.10)*cos((yy-2.22)*conv05))**2  &
                              +(yy+2.22)**2).lt.2.*sqrt(1100.)) then
              print * , '**** ADUJUSTING LAKE chongon LEVEL ****'
              print * , '     NEW:' , zchongon , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zchongon
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx-25.58)*cos((yy+61.73)*conv05))**2  &
                              +(yy-61.73)**2).lt.2.*sqrt(1100.)) then
              print * , '**** ADUJUSTING LAKE paijanne LEVEL ****'
              print * , '     NEW:' , zpaijanne , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zpaijanne
              if(i_lake==1) dhlake(i,j) = 17.0 ! 98.0 
            else if(ds1deg*sqrt( ((xx-98.85)*cos((yy+2.65)*conv05))**2  &
                              +(yy-2.65)**2).lt.2.*sqrt(1100.)) then
              print * , '**** ADUJUSTING LAKE toba LEVEL ****'
              print * , '     NEW:' , ztoba , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = ztoba
              if(i_lake==1) dhlake(i,j) = 100.0 ! 529.0 
            else if(ds1deg*sqrt( ((xx+60.78)*cos((yy+44.92)*conv05))**2  &
                              +(yy-44.92)**2).lt.2.*sqrt(1099.)) then
              print * , '**** ADUJUSTING LAKE bas_dor LEVEL ****'
              print * , '     NEW:' , zbas_dor , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zbas_dor
              if(i_lake==1) dhlake(i,j) = 15.0 ! 70.0 
            else if(ds1deg*sqrt( ((xx-116.95)*cos((yy+34.88)*conv05))**2  &
                              +(yy-34.88)**2).lt.2.*sqrt(1097.)) then
              print * , '**** ADUJUSTING LAKE weishan LEVEL ****'
              print * , '     NEW:' , zweishan , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zweishan
              if(i_lake==1) dhlake(i,j) = 1.0 ! 2.76 
            else if(ds1deg*sqrt( ((xx-82.20)*cos((yy+54.28)*conv05))**2  &
                              +(yy-54.28)**2).lt.2.*sqrt(1090.)) then
              print * , '**** ADUJUSTING LAKE novosibirsk LEVEL ****'
              print * , '     NEW:' , znovosibirsk , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = znovosibirsk
              if(i_lake==1) dhlake(i,j) = 9.0 ! 25.0 
            else if(ds1deg*sqrt( ((xx+72.58)*cos((yy-49.58)*conv05))**2  &
                              +(yy+49.58)**2).lt.2.*sqrt(1090.)) then
              print * , '**** ADUJUSTING LAKE viedma LEVEL ****'
              print * , '     NEW:' , zviedma , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zviedma
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx-100.37)*cos((yy+7.48)*conv05))**2  &
                              +(yy-7.48)**2).lt.2.*sqrt(1082.)) then
              print * , '**** ADUJUSTING LAKE songkhla LEVEL ****'
              print * , '     NEW:' , zsongkhla , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zsongkhla
              if(i_lake==1) dhlake(i,j) = 1.4 ! 2.0 
            else if(ds1deg*sqrt( ((xx-82.92)*cos((yy+44.92)*conv05))**2  &
                              +(yy-44.92)**2).lt.2.*sqrt(1070.)) then
              print * , '**** ADUJUSTING LAKE ebi LEVEL ****'
              print * , '     NEW:' , zebi , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zebi
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx+94.18)*cos((yy+54.57)*conv05))**2  &
                              +(yy-54.57)**2).lt.2.*sqrt(1060.)) then
              print * , '**** ADUJUSTING LAKE gods LEVEL ****'
              print * , '     NEW:' , zgods , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zgods
              if(i_lake==1) dhlake(i,j) = 13.2 ! 75.3 
            else if(ds1deg*sqrt( ((xx+72.12)*cos((yy+48.58)*conv05))**2  &
                              +(yy-48.58)**2).lt.2.*sqrt(1053.)) then
              print * , '**** ADUJUSTING LAKE saint_john LEVEL ****'
              print * , '     NEW:' , zsaint_john , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zsaint_john
              if(i_lake==1) dhlake(i,j) = 11.4 ! 63.1 
            else if(ds1deg*sqrt( ((xx-27.75)*cos((yy+69.00)*conv05))**2  &
                              +(yy-69.00)**2).lt.2.*sqrt(1050.)) then
              print * , '**** ADUJUSTING LAKE inari LEVEL ****'
              print * , '     NEW:' , zinari , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zinari
              if(i_lake==1) dhlake(i,j) = 14.4 ! 96.0 
            else if(ds1deg*sqrt( ((xx+72.88)*cos((yy+55.08)*conv05))**2  &
                              +(yy-55.08)**2).lt.2.*sqrt(1045.)) then
              print * , '**** ADUJUSTING LAKE bienville LEVEL ****'
              print * , '     NEW:' , zbienville , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zbienville
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx+94.67)*cos((yy+53.78)*conv05))**2  &
                              +(yy-53.78)**2).lt.2.*sqrt(1040.)) then
              print * , '**** ADUJUSTING LAKE island LEVEL ****'
              print * , '     NEW:' , zisland , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zisland
              if(i_lake==1) dhlake(i,j) = 20.1 ! 59.4
            else if(ds1deg*sqrt( ((xx+86.35)*cos((yy+12.35)*conv05))**2  &
                              +(yy-12.35)**2).lt.2.*sqrt(1040.)) then
              print * , '**** ADUJUSTING LAKE managua LEVEL ****'
              print * , '     NEW:' , zmanagua , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zmanagua
              if(i_lake==1) dhlake(i,j) = 7.8 ! 80.0 
            else if(ds1deg*sqrt( ((xx+76.62)*cos((yy+52.43)*conv05))**2  &
                              +(yy-52.43)**2).lt.2.*sqrt(1040.)) then
              print * , '**** ADUJUSTING LAKE opinaca LEVEL ****'
              print * , '     NEW:' , zopinaca , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zopinaca
              if(i_lake==1) dhlake(i,j) = 8.2 ! 51.0 
            else if(ds1deg*sqrt( ((xx+72.25)*cos((yy-48.75)*conv05))**2  &
                              +(yy+48.75)**2).lt.2.*sqrt(1038.)) then
              print * , '**** ADUJUSTING LAKE ohiggins LEVEL ****'
              print * , '     NEW:' , zohiggins , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zohiggins
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ?????
            else if(ds1deg*sqrt( ((xx+113.17)*cos((yy+66.57)*conv05))**2  &
                              +(yy-66.57)**2).lt.2.*sqrt(1030.)) then
              print * , '**** ADUJUSTING LAKE takiyuak LEVEL ****'
              print * , '     NEW:' , ztakiyuak , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = ztakiyuak
!             if(i_lake==1) dhlake(i,j) = 9.9999 ! ????? 
            else if(ds1deg*sqrt( ((xx-87.05)*cos((yy+42.08)*conv05))**2  &
                              +(yy-42.08)**2).lt.2.*sqrt(1010.)) then
              print * , '**** ADUJUSTING LAKE bositeng LEVEL ****'
              print * , '     NEW:' , zbositeng , '    OLD:' , htgrid(i,j) ,&
                 & i , j
              htgrid(i,j) = zbositeng
              if(i_lake==1) dhlake(i,j) = 7.66 ! 16.0 
            else
            end if
          end if
        end do
      end do
 
      end subroutine lakeadj

      end module mod_fudge
