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
!     **********  ksprecondition8f.f90 03/04/2012-01  **********

      module mod_ksprecondition
      implicit none

      contains

      subroutine ksprecondition(ncpx,ncpy,ncpz,nf,ncol,vre,avre)
      implicit none
      integer, intent(in)::ncpx,ncpy,ncpz,nf,ncol
      real*8, intent(in)::vre(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf,ncol)
      real*8, intent(out):: avre(ncpx,ncpy,ncpz,ncol)

      integer ix,iy,iz,is
      real*8 xx,yy,zz

      xx=1.0d0/12.0d0
      yy=1.0d0/12.0d0
      zz=1.0d0/12.0d0

      do is=1,ncol
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,is)=xx*vre(ix-1,iy,iz,is) &
                       +yy*vre(ix,iy-1,iz,is) &
                       +zz*vre(ix,iy,iz-1,is) &
                    +0.5d0*vre(ix,iy,iz,is) &
                       +zz*vre(ix,iy,iz+1,is) &
                       +yy*vre(ix,iy+1,iz,is) &
                       +xx*vre(ix+1,iy,iz,is)
      end do
      end do
      end do
!$omp end do nowait
      end do
      return
      end subroutine


      subroutine ksprecondition_c(ncpx,ncpy,ncpz,nf,ncol,vcm,avcm)
      implicit none
      integer, intent(in)::ncpx,ncpy,ncpz,nf,ncol
      complex*16,intent(in) ::vcm(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf,ncol)
      complex*16,intent(out)::avcm(ncpx,ncpy,ncpz,ncol)

      integer ix,iy,iz,is
      real*8 xx,yy,zz

      xx=1.0d0/12.0d0
      yy=1.0d0/12.0d0
      zz=1.0d0/12.0d0

      do is=1,ncol
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avcm(ix,iy,iz,is)=xx*vcm(ix-1,iy,iz,is) &
                       +yy*vcm(ix,iy-1,iz,is) &
                       +zz*vcm(ix,iy,iz-1,is) &
                    +0.5d0*vcm(ix,iy,iz,is) &
                       +zz*vcm(ix,iy,iz+1,is) &
                       +yy*vcm(ix,iy+1,iz,is) &
                       +xx*vcm(ix+1,iy,iz,is)
      end do
      end do
      end do
!$omp end do nowait
      end do
      return
      end subroutine


      end module
