!
!  Copyright 2023 RSPACE developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!     **********  radialhartree8e.f90 02/24/2011-01  **********

      module mod_radialhartree
      implicit none

      contains

      subroutine radialhartree_01(nradmx,lrhomx,nradct,pl,rho,vh,radial,dradial)
      implicit none
      integer nradmx,lrhomx,nradct,pl
      real*8 rho(nradmx,lrhomx),vh(nradmx,lrhomx)
      real*8 radial(nradmx),dradial(nradmx)
      integer ir,la,lla
      real*8 pi
      real*8 vh10,vh11,vh20,vh21,vh30
      pi=dacos(-1.0d0)

      do la=1,(2*pl+1)**2
!$omp do
        do ir=1,nradct
          vh(ir,la)=0.0d0
        end do
      end do
!$omp barrier
!$omp do schedule(static,1)
      do la=1,(2*pl+1)**2
        if (la .eq. 1) lla=0
        if ((la .ge. 2) .and. (la .le. 4)) lla=1
        if ((la .ge. 5) .and. (la .le. 9)) lla=2
        if ((la .ge. 10) .and. (la .le. 16)) lla=3
        if ((la .ge. 17) .and. (la .le. 25)) lla=4
        vh20=0.0d0
        vh30=0.0d0
        do ir=2,nradct-1
          vh20=vh20+rho(ir,la)/radial(ir)**(lla-1)*dradial(ir)
          vh30=vh30+rho(ir,la)*radial(ir)**(lla+2)*dradial(ir)
        end do
        vh10=0.0d0
        do ir=2,nradct-1
          vh11=rho(ir,la)*radial(ir)**(lla+2)*dradial(ir)
          vh21=rho(ir,la)/radial(ir)**(lla-1)*dradial(ir)
          vh10=vh10+vh11
          vh20=vh20-vh21
          vh(ir,la)=vh(ir,la)+4.0d0*pi*(vh10/radial(ir)**(lla+1) &
              +radial(ir)**lla*vh20-radial(ir)**lla &
              /radial(nradct)**(2*lla+1)*vh30)/(2.0d0*dble(lla)+1.0d0)
        end do
      end do
      return
      end subroutine


      subroutine radialhartree_02(nradmx,lrhomx,nrmax,npoint,pl,vhr,vhm,yylm)
      implicit none
      integer nradmx,nrmax,lrhomx,npoint,pl
      real*8 vhr(nradmx,npoint),vhm(nradmx,lrhomx),yylm(npoint,lrhomx)
      integer il,ir

      select case (pl)
      case (0)
!$omp do
        do il=1,npoint
          do ir=2,nrmax-1
            vhr(ir,il)=  vhm(ir, 1)*yylm(il, 1)
          end do
        end do
      case (1)
!$omp do
        do il=1,npoint
          do ir=2,nrmax-1
            vhr(ir,il)=  vhm(ir, 1)*yylm(il, 1) &
                       + vhm(ir, 2)*yylm(il, 2) &
                       + vhm(ir, 3)*yylm(il, 3) &
                       + vhm(ir, 4)*yylm(il, 4) &
                       + vhm(ir, 5)*yylm(il, 5) &
                       + vhm(ir, 6)*yylm(il, 6) &
                       + vhm(ir, 7)*yylm(il, 7) &
                       + vhm(ir, 8)*yylm(il, 8) &
                       + vhm(ir, 9)*yylm(il, 9)
          end do
        end do
      case (2)
!$omp do
        do il=1,npoint
          do ir=2,nrmax-1
            vhr(ir,il)=  vhm(ir, 1)*yylm(il, 1) &
                       + vhm(ir, 2)*yylm(il, 2) &
                       + vhm(ir, 3)*yylm(il, 3) &
                       + vhm(ir, 4)*yylm(il, 4) &
                       + vhm(ir, 5)*yylm(il, 5) &
                       + vhm(ir, 6)*yylm(il, 6) &
                       + vhm(ir, 7)*yylm(il, 7) &
                       + vhm(ir, 8)*yylm(il, 8) &
                       + vhm(ir, 9)*yylm(il, 9) &
                       + vhm(ir,10)*yylm(il,10) &
                       + vhm(ir,11)*yylm(il,11) &
                       + vhm(ir,12)*yylm(il,12) &
                       + vhm(ir,13)*yylm(il,13) &
                       + vhm(ir,14)*yylm(il,14) &
                       + vhm(ir,15)*yylm(il,15) &
                       + vhm(ir,16)*yylm(il,16) &
                       + vhm(ir,17)*yylm(il,17) &
                       + vhm(ir,18)*yylm(il,18) &
                       + vhm(ir,19)*yylm(il,19) &
                       + vhm(ir,20)*yylm(il,20) &
                       + vhm(ir,21)*yylm(il,21) &
                       + vhm(ir,22)*yylm(il,22) &
                       + vhm(ir,23)*yylm(il,23) &
                       + vhm(ir,24)*yylm(il,24) &
                       + vhm(ir,25)*yylm(il,25)
          end do
        end do
      end select
      return
      end subroutine


      end module
