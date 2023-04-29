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
!     **********  kslaplacian8f.F90 02/05/2014-01  **********

      module mod_kslaplacian
      implicit none
      real*8,allocatable::acfd(:)

      contains

      subroutine kslaplacian_initialize(nf)
      use mod_tools, only:tools_definefdcoef
      implicit none
      integer nf
      allocate(acfd(0:8))
      call tools_definefdcoef(nf,acfd(0),acfd(1),acfd(2),acfd(3),acfd(4),acfd(5),acfd(6),acfd(7),acfd(8))
      end subroutine


      subroutine kslaplacian_finalize
      deallocate(acfd)
      end subroutine


      subroutine kslaplacian_r(ncpx,ncpy,ncpz,neigmx,nf,dx,dy,dz,vre,veff,l,avre)
      use mod_stopp
      implicit none
      integer, intent(in)::ncpx,ncpy,ncpz,neigmx,l
      integer, intent(in)::nf
      real*8, intent(in)::dx,dy,dz
      real*8, intent(in)::vre(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf)
      real*8, intent(in)::veff(ncpx,ncpy,ncpz)
      real*8, intent(inout):: avre(ncpx,ncpy,ncpz,neigmx)
      integer ix,iy,iz
      real*8, allocatable::xx(:),yy(:),zz(:)

      allocate(xx(0:8),yy(0:8),zz(0:8))

      if ((nf<1) .and. (nf>8)) &
       call stopp('error in kslaplacian_r! This subroutine is for nf in 1-8.')

      xx(0)=acfd(0)/dx/dx*(-0.5d0)
      yy(0)=acfd(0)/dy/dy*(-0.5d0)
      zz(0)=acfd(0)/dz/dz*(-0.5d0)
      xx(1)=acfd(1)/dx/dx*(-0.5d0)
      yy(1)=acfd(1)/dy/dy*(-0.5d0)
      zz(1)=acfd(1)/dz/dz*(-0.5d0)
      xx(2)=acfd(2)/dx/dx*(-0.5d0)
      yy(2)=acfd(2)/dy/dy*(-0.5d0)
      zz(2)=acfd(2)/dz/dz*(-0.5d0)
      xx(3)=acfd(3)/dx/dx*(-0.5d0)
      yy(3)=acfd(3)/dy/dy*(-0.5d0)
      zz(3)=acfd(3)/dz/dz*(-0.5d0)
      xx(4)=acfd(4)/dx/dx*(-0.5d0)
      yy(4)=acfd(4)/dy/dy*(-0.5d0)
      zz(4)=acfd(4)/dz/dz*(-0.5d0)
      xx(5)=acfd(5)/dx/dx*(-0.5d0)
      yy(5)=acfd(5)/dy/dy*(-0.5d0)
      zz(5)=acfd(5)/dz/dz*(-0.5d0)
      xx(6)=acfd(6)/dx/dx*(-0.5d0)
      yy(6)=acfd(6)/dy/dy*(-0.5d0)
      zz(6)=acfd(6)/dz/dz*(-0.5d0)
      xx(7)=acfd(7)/dx/dx*(-0.5d0)
      yy(7)=acfd(7)/dy/dy*(-0.5d0)
      zz(7)=acfd(7)/dz/dz*(-0.5d0)
      xx(8)=acfd(8)/dx/dx*(-0.5d0)
      yy(8)=acfd(8)/dy/dy*(-0.5d0)
      zz(8)=acfd(8)/dz/dz*(-0.5d0)

!     ==========  the second derivative for nf=1  ==========
      if (nf .eq. 1) then
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +xx(1)*vre(ix-1,iy,iz) &
            +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vre(ix,iy,iz) &
            +xx(1)*vre(ix+1,iy,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +yy(1)*vre(ix,iy-1,iz) &
            +yy(1)*vre(ix,iy+1,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +zz(1)*vre(ix,iy,iz-1) &
            +zz(1)*vre(ix,iy,iz+1)
      end do
      end do
      end do
      end if
!     ======================================================

!     ==========  the second derivative for nf=2  ==========
      if (nf .eq. 2) then
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +xx(2)*vre(ix-2,iy,iz) &
            +xx(1)*vre(ix-1,iy,iz) &
            +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vre(ix,iy,iz) &
            +xx(1)*vre(ix+1,iy,iz) &
            +xx(2)*vre(ix+2,iy,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +yy(2)*vre(ix,iy-2,iz) &
            +yy(1)*vre(ix,iy-1,iz) &
            +yy(1)*vre(ix,iy+1,iz) &
            +yy(2)*vre(ix,iy+2,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +zz(2)*vre(ix,iy,iz-2) &
            +zz(1)*vre(ix,iy,iz-1) &
            +zz(1)*vre(ix,iy,iz+1) &
            +zz(2)*vre(ix,iy,iz+2)
      end do
      end do
      end do
      end if
!     ======================================================

!     ==========  the second derivative for nf=3  ==========
      if (nf .eq. 3) then
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +xx(3)*vre(ix-3,iy,iz) &
            +xx(2)*vre(ix-2,iy,iz) &
            +xx(1)*vre(ix-1,iy,iz) &
            +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vre(ix,iy,iz) &
            +xx(1)*vre(ix+1,iy,iz) &
            +xx(2)*vre(ix+2,iy,iz) &
            +xx(3)*vre(ix+3,iy,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +yy(3)*vre(ix,iy-3,iz) &
            +yy(2)*vre(ix,iy-2,iz) &
            +yy(1)*vre(ix,iy-1,iz) &
            +yy(1)*vre(ix,iy+1,iz) &
            +yy(2)*vre(ix,iy+2,iz) &
            +yy(3)*vre(ix,iy+3,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +zz(3)*vre(ix,iy,iz-3) &
            +zz(2)*vre(ix,iy,iz-2) &
            +zz(1)*vre(ix,iy,iz-1) &
            +zz(1)*vre(ix,iy,iz+1) &
            +zz(2)*vre(ix,iy,iz+2) &
            +zz(3)*vre(ix,iy,iz+3)
      end do
      end do
      end do
      end if
!     ======================================================

!     ==========  the second derivative for nf=4  ==========
      if (nf .eq. 4) then
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +xx(4)*vre(ix-4,iy,iz) &
            +xx(3)*vre(ix-3,iy,iz) &
            +xx(2)*vre(ix-2,iy,iz) &
            +xx(1)*vre(ix-1,iy,iz) &
            +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vre(ix,iy,iz) &
            +xx(1)*vre(ix+1,iy,iz) &
            +xx(2)*vre(ix+2,iy,iz) &
            +xx(3)*vre(ix+3,iy,iz) &
            +xx(4)*vre(ix+4,iy,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +yy(4)*vre(ix,iy-4,iz) &
            +yy(3)*vre(ix,iy-3,iz) &
            +yy(2)*vre(ix,iy-2,iz) &
            +yy(1)*vre(ix,iy-1,iz) &
            +yy(1)*vre(ix,iy+1,iz) &
            +yy(2)*vre(ix,iy+2,iz) &
            +yy(3)*vre(ix,iy+3,iz) &
            +yy(4)*vre(ix,iy+4,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +zz(4)*vre(ix,iy,iz-4) &
            +zz(3)*vre(ix,iy,iz-3) &
            +zz(2)*vre(ix,iy,iz-2) &
            +zz(1)*vre(ix,iy,iz-1) &
            +zz(1)*vre(ix,iy,iz+1) &
            +zz(2)*vre(ix,iy,iz+2) &
            +zz(3)*vre(ix,iy,iz+3) &
            +zz(4)*vre(ix,iy,iz+4)
      end do
      end do
      end do
      end if
!     ======================================================

!     ==========  the second derivative for nf=5  ==========
      if (nf .eq. 5) then
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +xx(5)*vre(ix-5,iy,iz) &
            +xx(4)*vre(ix-4,iy,iz) &
            +xx(3)*vre(ix-3,iy,iz) &
            +xx(2)*vre(ix-2,iy,iz) &
            +xx(1)*vre(ix-1,iy,iz) &
            +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vre(ix,iy,iz) &
            +xx(1)*vre(ix+1,iy,iz) &
            +xx(2)*vre(ix+2,iy,iz) &
            +xx(3)*vre(ix+3,iy,iz) &
            +xx(4)*vre(ix+4,iy,iz) &
            +xx(5)*vre(ix+5,iy,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +yy(5)*vre(ix,iy-5,iz) &
            +yy(4)*vre(ix,iy-4,iz) &
            +yy(3)*vre(ix,iy-3,iz) &
            +yy(2)*vre(ix,iy-2,iz) &
            +yy(1)*vre(ix,iy-1,iz) &
            +yy(1)*vre(ix,iy+1,iz) &
            +yy(2)*vre(ix,iy+2,iz) &
            +yy(3)*vre(ix,iy+3,iz) &
            +yy(4)*vre(ix,iy+4,iz) &
            +yy(5)*vre(ix,iy+5,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +zz(5)*vre(ix,iy,iz-5) &
            +zz(4)*vre(ix,iy,iz-4) &
            +zz(3)*vre(ix,iy,iz-3) &
            +zz(2)*vre(ix,iy,iz-2) &
            +zz(1)*vre(ix,iy,iz-1) &
            +zz(1)*vre(ix,iy,iz+1) &
            +zz(2)*vre(ix,iy,iz+2) &
            +zz(3)*vre(ix,iy,iz+3) &
            +zz(4)*vre(ix,iy,iz+4) &
            +zz(5)*vre(ix,iy,iz+5)
      end do
      end do
      end do
      end if
!     ======================================================

!     ==========  the second derivative for nf=6  ==========
      if (nf .eq. 6) then
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +xx(6)*vre(ix-6,iy,iz) &
            +xx(5)*vre(ix-5,iy,iz) &
            +xx(4)*vre(ix-4,iy,iz) &
            +xx(3)*vre(ix-3,iy,iz) &
            +xx(2)*vre(ix-2,iy,iz) &
            +xx(1)*vre(ix-1,iy,iz) &
            +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vre(ix,iy,iz) &
            +xx(1)*vre(ix+1,iy,iz) &
            +xx(2)*vre(ix+2,iy,iz) &
            +xx(3)*vre(ix+3,iy,iz) &
            +xx(4)*vre(ix+4,iy,iz) &
            +xx(5)*vre(ix+5,iy,iz) &
            +xx(6)*vre(ix+6,iy,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +yy(6)*vre(ix,iy-6,iz) &
            +yy(5)*vre(ix,iy-5,iz) &
            +yy(4)*vre(ix,iy-4,iz) &
            +yy(3)*vre(ix,iy-3,iz) &
            +yy(2)*vre(ix,iy-2,iz) &
            +yy(1)*vre(ix,iy-1,iz) &
            +yy(1)*vre(ix,iy+1,iz) &
            +yy(2)*vre(ix,iy+2,iz) &
            +yy(3)*vre(ix,iy+3,iz) &
            +yy(4)*vre(ix,iy+4,iz) &
            +yy(5)*vre(ix,iy+5,iz) &
            +yy(6)*vre(ix,iy+6,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +zz(6)*vre(ix,iy,iz-6) &
            +zz(5)*vre(ix,iy,iz-5) &
            +zz(4)*vre(ix,iy,iz-4) &
            +zz(3)*vre(ix,iy,iz-3) &
            +zz(2)*vre(ix,iy,iz-2) &
            +zz(1)*vre(ix,iy,iz-1) &
            +zz(1)*vre(ix,iy,iz+1) &
            +zz(2)*vre(ix,iy,iz+2) &
            +zz(3)*vre(ix,iy,iz+3) &
            +zz(4)*vre(ix,iy,iz+4) &
            +zz(5)*vre(ix,iy,iz+5) &
            +zz(6)*vre(ix,iy,iz+6)
      end do
      end do
      end do
      end if
!     ======================================================

!     ==========  the second derivative for nf=7  ==========
      if (nf .eq. 7) then
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +xx(7)*vre(ix-7,iy,iz) &
            +xx(6)*vre(ix-6,iy,iz) &
            +xx(5)*vre(ix-5,iy,iz) &
            +xx(4)*vre(ix-4,iy,iz) &
            +xx(3)*vre(ix-3,iy,iz) &
            +xx(2)*vre(ix-2,iy,iz) &
            +xx(1)*vre(ix-1,iy,iz) &
            +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vre(ix,iy,iz) &
            +xx(1)*vre(ix+1,iy,iz) &
            +xx(2)*vre(ix+2,iy,iz) &
            +xx(3)*vre(ix+3,iy,iz) &
            +xx(4)*vre(ix+4,iy,iz) &
            +xx(5)*vre(ix+5,iy,iz) &
            +xx(6)*vre(ix+6,iy,iz) &
            +xx(7)*vre(ix+7,iy,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +yy(7)*vre(ix,iy-7,iz) &
            +yy(6)*vre(ix,iy-6,iz) &
            +yy(5)*vre(ix,iy-5,iz) &
            +yy(4)*vre(ix,iy-4,iz) &
            +yy(3)*vre(ix,iy-3,iz) &
            +yy(2)*vre(ix,iy-2,iz) &
            +yy(1)*vre(ix,iy-1,iz) &
            +yy(1)*vre(ix,iy+1,iz) &
            +yy(2)*vre(ix,iy+2,iz) &
            +yy(3)*vre(ix,iy+3,iz) &
            +yy(4)*vre(ix,iy+4,iz) &
            +yy(5)*vre(ix,iy+5,iz) &
            +yy(6)*vre(ix,iy+6,iz) &
            +yy(7)*vre(ix,iy+7,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +zz(7)*vre(ix,iy,iz-7) &
            +zz(6)*vre(ix,iy,iz-6) &
            +zz(5)*vre(ix,iy,iz-5) &
            +zz(4)*vre(ix,iy,iz-4) &
            +zz(3)*vre(ix,iy,iz-3) &
            +zz(2)*vre(ix,iy,iz-2) &
            +zz(1)*vre(ix,iy,iz-1) &
            +zz(1)*vre(ix,iy,iz+1) &
            +zz(2)*vre(ix,iy,iz+2) &
            +zz(3)*vre(ix,iy,iz+3) &
            +zz(4)*vre(ix,iy,iz+4) &
            +zz(5)*vre(ix,iy,iz+5) &
            +zz(6)*vre(ix,iy,iz+6) &
            +zz(7)*vre(ix,iy,iz+7)
      end do
      end do
      end do
      end if
!     ======================================================

!     ==========  the second derivative for nf=8  ==========
      if (nf .eq. 8) then
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +xx(8)*vre(ix-8,iy,iz) &
            +xx(7)*vre(ix-7,iy,iz) &
            +xx(6)*vre(ix-6,iy,iz) &
            +xx(5)*vre(ix-5,iy,iz) &
            +xx(4)*vre(ix-4,iy,iz) &
            +xx(3)*vre(ix-3,iy,iz) &
            +xx(2)*vre(ix-2,iy,iz) &
            +xx(1)*vre(ix-1,iy,iz) &
            +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vre(ix,iy,iz) &
            +xx(1)*vre(ix+1,iy,iz) &
            +xx(2)*vre(ix+2,iy,iz) &
            +xx(3)*vre(ix+3,iy,iz) &
            +xx(4)*vre(ix+4,iy,iz) &
            +xx(5)*vre(ix+5,iy,iz) &
            +xx(6)*vre(ix+6,iy,iz) &
            +xx(7)*vre(ix+7,iy,iz) &
            +xx(8)*vre(ix+8,iy,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +yy(8)*vre(ix,iy-8,iz) &
            +yy(7)*vre(ix,iy-7,iz) &
            +yy(6)*vre(ix,iy-6,iz) &
            +yy(5)*vre(ix,iy-5,iz) &
            +yy(4)*vre(ix,iy-4,iz) &
            +yy(3)*vre(ix,iy-3,iz) &
            +yy(2)*vre(ix,iy-2,iz) &
            +yy(1)*vre(ix,iy-1,iz) &
            +yy(1)*vre(ix,iy+1,iz) &
            +yy(2)*vre(ix,iy+2,iz) &
            +yy(3)*vre(ix,iy+3,iz) &
            +yy(4)*vre(ix,iy+4,iz) &
            +yy(5)*vre(ix,iy+5,iz) &
            +yy(6)*vre(ix,iy+6,iz) &
            +yy(7)*vre(ix,iy+7,iz) &
            +yy(8)*vre(ix,iy+8,iz)
      end do
      end do
      end do
!cdir nodep
!$omp do
!ocl norecurrence(avre)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
      avre(ix,iy,iz,l)=avre(ix,iy,iz,l) &
            +zz(8)*vre(ix,iy,iz-8) &
            +zz(7)*vre(ix,iy,iz-7) &
            +zz(6)*vre(ix,iy,iz-6) &
            +zz(5)*vre(ix,iy,iz-5) &
            +zz(4)*vre(ix,iy,iz-4) &
            +zz(3)*vre(ix,iy,iz-3) &
            +zz(2)*vre(ix,iy,iz-2) &
            +zz(1)*vre(ix,iy,iz-1) &
            +zz(1)*vre(ix,iy,iz+1) &
            +zz(2)*vre(ix,iy,iz+2) &
            +zz(3)*vre(ix,iy,iz+3) &
            +zz(4)*vre(ix,iy,iz+4) &
            +zz(5)*vre(ix,iy,iz+5) &
            +zz(6)*vre(ix,iy,iz+6) &
            +zz(7)*vre(ix,iy,iz+7) &
            +zz(8)*vre(ix,iy,iz+8)
      end do
      end do
      end do
      end if
!     ======================================================
      deallocate(xx,yy,zz)
      return
      end subroutine


      subroutine kslaplacian_c(ncpx,ncpy,ncpz,neigmx,ncol,nvef,nf,dx,dy,dz,xmax,ymax,zmax,skpx,skpy,skpz &
                               ,vcm,veff,l,avcm,workc)
      use mod_mpi, only: myrx,myry,myrz
      implicit none
      integer,   intent(in)::ncpx,ncpy,ncpz,neigmx,ncol,nvef,l
      integer,   intent(in)::nf
      real*8,    intent(in)::dx,dy,dz
      real*8,    intent(in)::xmax,ymax,zmax
      real*8,    intent(in)::skpx,skpy,skpz
      complex*16,intent(inout)::vcm(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf,ncol)
      real*8,    intent(in)::veff(ncpx,ncpy,ncpz,nvef)
      complex*16,intent(inout)::avcm(ncpx,ncpy,ncpz,neigmx,ncol)
      complex*16,intent(inout)::workc(ncpx,ncpy,ncpz,ncol)
      integer ix,iy,iz,jx,jy,jz,ns,nsv
      real*8 x,y,z
      complex*16 cuniti
      cuniti=dcmplx(0.0d0,1.0d0)

      if (ncol == 1) then
!$omp do
        do iz=-(nf-1),ncpz+nf
        do iy=-(nf-1),ncpy+nf
        do ix=-(nf-1),ncpx+nf
          jz=myrz*ncpz+iz
          jy=myry*ncpy+iy
          jx=myrx*ncpx+ix
          x=jx*dx-xmax-0.5d0*dx
          y=jy*dy-ymax-0.5d0*dy
          z=jz*dz-zmax-0.5d0*dz
          vcm(ix,iy,iz,1)=exp(cuniti*(skpx*x+skpy*y+skpz*z))*vcm(ix,iy,iz,1)
        end do
        end do
        end do
      else
!$omp do
        do iz=-(nf-1),ncpz+nf
        do iy=-(nf-1),ncpy+nf
        do ix=-(nf-1),ncpx+nf
          jz=myrz*ncpz+iz
          jy=myry*ncpy+iy
          jx=myrx*ncpx+ix
          x=jx*dx-xmax-0.5d0*dx
          y=jy*dy-ymax-0.5d0*dy
          z=jz*dz-zmax-0.5d0*dz
          vcm(ix,iy,iz,   1)=exp(cuniti*(skpx*x+skpy*y+skpz*z))*vcm(ix,iy,iz,   1)
          vcm(ix,iy,iz,ncol)=exp(cuniti*(skpx*x+skpy*y+skpz*z))*vcm(ix,iy,iz,ncol)
        end do
        end do
        end do
      end if

!$omp single
      workc=dcmplx(0.0d0,0.0d0)
!$omp end single
      do ns= 1,ncol
        nsv= 1
        if (nvef>1) nsv= ns
        call kslaplacian_c_01(  &
         ncpx,ncpy,ncpz,nf,dx,dy,dz,  &
         vcm(1-nf,1-nf,1-nf,ns),veff(1,1,1,nsv),workc(1,1,1,ns))
!$op barrier
      enddo
      if (nvef==4) call kslaplacian_c_02(ncpx,ncpy,ncpz,nf,vcm,veff,workc)

      if (ncol == 1) then
!$omp do
        do iz=1,ncpz
        do iy=1,ncpy
        do ix=1,ncpx
          jz=myrz*ncpz+iz
          jy=myry*ncpy+iy
          jx=myrx*ncpx+ix
          x=jx*dx-xmax-0.5d0*dx
          y=jy*dy-ymax-0.5d0*dy
          z=jz*dz-zmax-0.5d0*dz
          avcm(ix,iy,iz,l,1)=avcm(ix,iy,iz,l,1)+exp(-cuniti*(skpx*x+skpy*y+skpz*z))*workc(ix,iy,iz,1)
        end do
        end do
        end do
      else
!$omp do
        do iz=1,ncpz
        do iy=1,ncpy
        do ix=1,ncpx
          jz=myrz*ncpz+iz
          jy=myry*ncpy+iy
          jx=myrx*ncpx+ix
          x=jx*dx-xmax-0.5d0*dx
          y=jy*dy-ymax-0.5d0*dy
          z=jz*dz-zmax-0.5d0*dz
          avcm(ix,iy,iz,l,   1)=avcm(ix,iy,iz,l,   1)+exp(-cuniti*(skpx*x+skpy*y+skpz*z))*workc(ix,iy,iz,   1)
          avcm(ix,iy,iz,l,ncol)=avcm(ix,iy,iz,l,ncol)+exp(-cuniti*(skpx*x+skpy*y+skpz*z))*workc(ix,iy,iz,ncol)
        end do
        end do
        end do
      end if
      return
      end subroutine kslaplacian_c


      subroutine kslaplacian_c_01(ncpx,ncpy,ncpz,nf,dx,dy,dz,vcm,veff,avcm)
      use mod_stopp
      implicit none
      integer,   intent(in)::ncpx,ncpy,ncpz
      integer,   intent(in)::nf
      real*8,    intent(in)::dx,dy,dz
      complex*16,intent(in)::vcm(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf)
      real*8,    intent(in)::veff(ncpx,ncpy,ncpz)
      complex*16,intent(inout)::avcm(ncpx,ncpy,ncpz)
      integer ix,iy,iz
      real*8, allocatable::xx(:),yy(:),zz(:)

      allocate(xx(0:8),yy(0:8),zz(0:8))

      if ((nf<1) .and. (nf>8)) &
       call stopp('error in kslaplacian_c! This subroutine is for nf in 1-8.')

      xx(0)=acfd(0)/dx/dx*(-0.5d0)
      yy(0)=acfd(0)/dy/dy*(-0.5d0)
      zz(0)=acfd(0)/dz/dz*(-0.5d0)
      xx(1)=acfd(1)/dx/dx*(-0.5d0)
      yy(1)=acfd(1)/dy/dy*(-0.5d0)
      zz(1)=acfd(1)/dz/dz*(-0.5d0)
      xx(2)=acfd(2)/dx/dx*(-0.5d0)
      yy(2)=acfd(2)/dy/dy*(-0.5d0)
      zz(2)=acfd(2)/dz/dz*(-0.5d0)
      xx(3)=acfd(3)/dx/dx*(-0.5d0)
      yy(3)=acfd(3)/dy/dy*(-0.5d0)
      zz(3)=acfd(3)/dz/dz*(-0.5d0)
      xx(4)=acfd(4)/dx/dx*(-0.5d0)
      yy(4)=acfd(4)/dy/dy*(-0.5d0)
      zz(4)=acfd(4)/dz/dz*(-0.5d0)
      xx(5)=acfd(5)/dx/dx*(-0.5d0)
      yy(5)=acfd(5)/dy/dy*(-0.5d0)
      zz(5)=acfd(5)/dz/dz*(-0.5d0)
      xx(6)=acfd(6)/dx/dx*(-0.5d0)
      yy(6)=acfd(6)/dy/dy*(-0.5d0)
      zz(6)=acfd(6)/dz/dz*(-0.5d0)
      xx(7)=acfd(7)/dx/dx*(-0.5d0)
      yy(7)=acfd(7)/dy/dy*(-0.5d0)
      zz(7)=acfd(7)/dz/dz*(-0.5d0)
      xx(8)=acfd(8)/dx/dx*(-0.5d0)
      yy(8)=acfd(8)/dy/dy*(-0.5d0)
      zz(8)=acfd(8)/dz/dz*(-0.5d0)

!     ==========  the second derivative for nf=1  ==========
      if (nf .eq. 1) then
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +xx(1)*vcm(ix-1,iy,iz) &
          +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vcm(ix,iy,iz) &
          +xx(1)*vcm(ix+1,iy,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +yy(1)*vcm(ix,iy-1,iz) &
          +yy(1)*vcm(ix,iy+1,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +zz(1)*vcm(ix,iy,iz-1) &
          +zz(1)*vcm(ix,iy,iz+1)
      end do
      end do
      end do
!$omp end do nowait
      end if
!     ======================================================

!     ==========  the second derivative for nf=2  ==========
      if (nf .eq. 2) then
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +xx(2)*vcm(ix-2,iy,iz) &
          +xx(1)*vcm(ix-1,iy,iz) &
          +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vcm(ix,iy,iz) &
          +xx(1)*vcm(ix+1,iy,iz) &
          +xx(2)*vcm(ix+2,iy,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +yy(2)*vcm(ix,iy-2,iz) &
          +yy(1)*vcm(ix,iy-1,iz) &
          +yy(1)*vcm(ix,iy+1,iz) &
          +yy(2)*vcm(ix,iy+2,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +zz(2)*vcm(ix,iy,iz-2) &
          +zz(1)*vcm(ix,iy,iz-1) &
          +zz(1)*vcm(ix,iy,iz+1) &
          +zz(2)*vcm(ix,iy,iz+2)
      end do
      end do
      end do
!$omp end do nowait
      end if
!     ======================================================

!     ==========  the second derivative for nf=3  ==========
      if (nf .eq. 3) then
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +xx(3)*vcm(ix-3,iy,iz) &
          +xx(2)*vcm(ix-2,iy,iz) &
          +xx(1)*vcm(ix-1,iy,iz) &
          +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vcm(ix,iy,iz) &
          +xx(1)*vcm(ix+1,iy,iz) &
          +xx(2)*vcm(ix+2,iy,iz) &
          +xx(3)*vcm(ix+3,iy,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +yy(3)*vcm(ix,iy-3,iz) &
          +yy(2)*vcm(ix,iy-2,iz) &
          +yy(1)*vcm(ix,iy-1,iz) &
          +yy(1)*vcm(ix,iy+1,iz) &
          +yy(2)*vcm(ix,iy+2,iz) &
          +yy(3)*vcm(ix,iy+3,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +zz(3)*vcm(ix,iy,iz-3) &
          +zz(2)*vcm(ix,iy,iz-2) &
          +zz(1)*vcm(ix,iy,iz-1) &
          +zz(1)*vcm(ix,iy,iz+1) &
          +zz(2)*vcm(ix,iy,iz+2) &
          +zz(3)*vcm(ix,iy,iz+3)
      end do
      end do
      end do
!$omp end do nowait
      end if
!     ======================================================

!     ==========  the second derivative for nf=4  ==========
      if (nf .eq. 4) then
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +xx(4)*vcm(ix-4,iy,iz) &
          +xx(3)*vcm(ix-3,iy,iz) &
          +xx(2)*vcm(ix-2,iy,iz) &
          +xx(1)*vcm(ix-1,iy,iz) &
          +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vcm(ix,iy,iz) &
          +xx(1)*vcm(ix+1,iy,iz) &
          +xx(2)*vcm(ix+2,iy,iz) &
          +xx(3)*vcm(ix+3,iy,iz) &
          +xx(4)*vcm(ix+4,iy,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +yy(4)*vcm(ix,iy-4,iz) &
          +yy(3)*vcm(ix,iy-3,iz) &
          +yy(2)*vcm(ix,iy-2,iz) &
          +yy(1)*vcm(ix,iy-1,iz) &
          +yy(1)*vcm(ix,iy+1,iz) &
          +yy(2)*vcm(ix,iy+2,iz) &
          +yy(3)*vcm(ix,iy+3,iz) &
          +yy(4)*vcm(ix,iy+4,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +zz(4)*vcm(ix,iy,iz-4) &
          +zz(3)*vcm(ix,iy,iz-3) &
          +zz(2)*vcm(ix,iy,iz-2) &
          +zz(1)*vcm(ix,iy,iz-1) &
          +zz(1)*vcm(ix,iy,iz+1) &
          +zz(2)*vcm(ix,iy,iz+2) &
          +zz(3)*vcm(ix,iy,iz+3) &
          +zz(4)*vcm(ix,iy,iz+4)
      end do
      end do
      end do
!$omp end do nowait
      end if
!     ======================================================

!     ==========  the second derivative for nf=5  ==========
      if (nf .eq. 5) then
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +xx(5)*vcm(ix-5,iy,iz) &
          +xx(4)*vcm(ix-4,iy,iz) &
          +xx(3)*vcm(ix-3,iy,iz) &
          +xx(2)*vcm(ix-2,iy,iz) &
          +xx(1)*vcm(ix-1,iy,iz) &
          +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vcm(ix,iy,iz) &
          +xx(1)*vcm(ix+1,iy,iz) &
          +xx(2)*vcm(ix+2,iy,iz) &
          +xx(3)*vcm(ix+3,iy,iz) &
          +xx(4)*vcm(ix+4,iy,iz) &
          +xx(5)*vcm(ix+5,iy,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +yy(5)*vcm(ix,iy-5,iz) &
          +yy(4)*vcm(ix,iy-4,iz) &
          +yy(3)*vcm(ix,iy-3,iz) &
          +yy(2)*vcm(ix,iy-2,iz) &
          +yy(1)*vcm(ix,iy-1,iz) &
          +yy(1)*vcm(ix,iy+1,iz) &
          +yy(2)*vcm(ix,iy+2,iz) &
          +yy(3)*vcm(ix,iy+3,iz) &
          +yy(4)*vcm(ix,iy+4,iz) &
          +yy(5)*vcm(ix,iy+5,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +zz(5)*vcm(ix,iy,iz-5) &
          +zz(4)*vcm(ix,iy,iz-4) &
          +zz(3)*vcm(ix,iy,iz-3) &
          +zz(2)*vcm(ix,iy,iz-2) &
          +zz(1)*vcm(ix,iy,iz-1) &
          +zz(1)*vcm(ix,iy,iz+1) &
          +zz(2)*vcm(ix,iy,iz+2) &
          +zz(3)*vcm(ix,iy,iz+3) &
          +zz(4)*vcm(ix,iy,iz+4) &
          +zz(5)*vcm(ix,iy,iz+5)
      end do
      end do
      end do
!$omp end do nowait
      end if
!     ======================================================

!     ==========  the second derivative for nf=6  ==========
      if (nf .eq. 6) then
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +xx(6)*vcm(ix-6,iy,iz) &
          +xx(5)*vcm(ix-5,iy,iz) &
          +xx(4)*vcm(ix-4,iy,iz) &
          +xx(3)*vcm(ix-3,iy,iz) &
          +xx(2)*vcm(ix-2,iy,iz) &
          +xx(1)*vcm(ix-1,iy,iz) &
          +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vcm(ix,iy,iz) &
          +xx(1)*vcm(ix+1,iy,iz) &
          +xx(2)*vcm(ix+2,iy,iz) &
          +xx(3)*vcm(ix+3,iy,iz) &
          +xx(4)*vcm(ix+4,iy,iz) &
          +xx(5)*vcm(ix+5,iy,iz) &
          +xx(6)*vcm(ix+6,iy,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +yy(6)*vcm(ix,iy-6,iz) &
          +yy(5)*vcm(ix,iy-5,iz) &
          +yy(4)*vcm(ix,iy-4,iz) &
          +yy(3)*vcm(ix,iy-3,iz) &
          +yy(2)*vcm(ix,iy-2,iz) &
          +yy(1)*vcm(ix,iy-1,iz) &
          +yy(1)*vcm(ix,iy+1,iz) &
          +yy(2)*vcm(ix,iy+2,iz) &
          +yy(3)*vcm(ix,iy+3,iz) &
          +yy(4)*vcm(ix,iy+4,iz) &
          +yy(5)*vcm(ix,iy+5,iz) &
          +yy(6)*vcm(ix,iy+6,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +zz(6)*vcm(ix,iy,iz-6) &
          +zz(5)*vcm(ix,iy,iz-5) &
          +zz(4)*vcm(ix,iy,iz-4) &
          +zz(3)*vcm(ix,iy,iz-3) &
          +zz(2)*vcm(ix,iy,iz-2) &
          +zz(1)*vcm(ix,iy,iz-1) &
          +zz(1)*vcm(ix,iy,iz+1) &
          +zz(2)*vcm(ix,iy,iz+2) &
          +zz(3)*vcm(ix,iy,iz+3) &
          +zz(4)*vcm(ix,iy,iz+4) &
          +zz(5)*vcm(ix,iy,iz+5) &
          +zz(6)*vcm(ix,iy,iz+6)
      end do
      end do
      end do
!$omp end do nowait
      end if
!     ======================================================

!     ==========  the second derivative for nf=7  ==========
      if (nf .eq. 7) then
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +xx(7)*vcm(ix-7,iy,iz) &
          +xx(6)*vcm(ix-6,iy,iz) &
          +xx(5)*vcm(ix-5,iy,iz) &
          +xx(4)*vcm(ix-4,iy,iz) &
          +xx(3)*vcm(ix-3,iy,iz) &
          +xx(2)*vcm(ix-2,iy,iz) &
          +xx(1)*vcm(ix-1,iy,iz) &
          +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vcm(ix,iy,iz) &
          +xx(1)*vcm(ix+1,iy,iz) &
          +xx(2)*vcm(ix+2,iy,iz) &
          +xx(3)*vcm(ix+3,iy,iz) &
          +xx(4)*vcm(ix+4,iy,iz) &
          +xx(5)*vcm(ix+5,iy,iz) &
          +xx(6)*vcm(ix+6,iy,iz) &
          +xx(7)*vcm(ix+7,iy,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +yy(7)*vcm(ix,iy-7,iz) &
          +yy(6)*vcm(ix,iy-6,iz) &
          +yy(5)*vcm(ix,iy-5,iz) &
          +yy(4)*vcm(ix,iy-4,iz) &
          +yy(3)*vcm(ix,iy-3,iz) &
          +yy(2)*vcm(ix,iy-2,iz) &
          +yy(1)*vcm(ix,iy-1,iz) &
          +yy(1)*vcm(ix,iy+1,iz) &
          +yy(2)*vcm(ix,iy+2,iz) &
          +yy(3)*vcm(ix,iy+3,iz) &
          +yy(4)*vcm(ix,iy+4,iz) &
          +yy(5)*vcm(ix,iy+5,iz) &
          +yy(6)*vcm(ix,iy+6,iz) &
          +yy(7)*vcm(ix,iy+7,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +zz(7)*vcm(ix,iy,iz-7) &
          +zz(6)*vcm(ix,iy,iz-6) &
          +zz(5)*vcm(ix,iy,iz-5) &
          +zz(4)*vcm(ix,iy,iz-4) &
          +zz(3)*vcm(ix,iy,iz-3) &
          +zz(2)*vcm(ix,iy,iz-2) &
          +zz(1)*vcm(ix,iy,iz-1) &
          +zz(1)*vcm(ix,iy,iz+1) &
          +zz(2)*vcm(ix,iy,iz+2) &
          +zz(3)*vcm(ix,iy,iz+3) &
          +zz(4)*vcm(ix,iy,iz+4) &
          +zz(5)*vcm(ix,iy,iz+5) &
          +zz(6)*vcm(ix,iy,iz+6) &
          +zz(7)*vcm(ix,iy,iz+7)
      end do
      end do
      end do
!$omp end do nowait
      end if
!     ======================================================

!     ==========  the second derivative for nf=8  ==========
      if (nf .eq. 8) then
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +xx(8)*vcm(ix-8,iy,iz) &
          +xx(7)*vcm(ix-7,iy,iz) &
          +xx(6)*vcm(ix-6,iy,iz) &
          +xx(5)*vcm(ix-5,iy,iz) &
          +xx(4)*vcm(ix-4,iy,iz) &
          +xx(3)*vcm(ix-3,iy,iz) &
          +xx(2)*vcm(ix-2,iy,iz) &
          +xx(1)*vcm(ix-1,iy,iz) &
          +(xx(0)+yy(0)+zz(0)+veff(ix,iy,iz))*vcm(ix,iy,iz) &
          +xx(1)*vcm(ix+1,iy,iz) &
          +xx(2)*vcm(ix+2,iy,iz) &
          +xx(3)*vcm(ix+3,iy,iz) &
          +xx(4)*vcm(ix+4,iy,iz) &
          +xx(5)*vcm(ix+5,iy,iz) &
          +xx(6)*vcm(ix+6,iy,iz) &
          +xx(7)*vcm(ix+7,iy,iz) &
          +xx(8)*vcm(ix+8,iy,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +yy(8)*vcm(ix,iy-8,iz) &
          +yy(7)*vcm(ix,iy-7,iz) &
          +yy(6)*vcm(ix,iy-6,iz) &
          +yy(5)*vcm(ix,iy-5,iz) &
          +yy(4)*vcm(ix,iy-4,iz) &
          +yy(3)*vcm(ix,iy-3,iz) &
          +yy(2)*vcm(ix,iy-2,iz) &
          +yy(1)*vcm(ix,iy-1,iz) &
          +yy(1)*vcm(ix,iy+1,iz) &
          +yy(2)*vcm(ix,iy+2,iz) &
          +yy(3)*vcm(ix,iy+3,iz) &
          +yy(4)*vcm(ix,iy+4,iz) &
          +yy(5)*vcm(ix,iy+5,iz) &
          +yy(6)*vcm(ix,iy+6,iz) &
          +yy(7)*vcm(ix,iy+7,iz) &
          +yy(8)*vcm(ix,iy+8,iz)
      end do
      end do
      end do
!$omp end do nowait
!cdir nodep
!$omp do
!ocl norecurrence(avcm)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz)=avcm(ix,iy,iz) &
          +zz(8)*vcm(ix,iy,iz-8) &
          +zz(7)*vcm(ix,iy,iz-7) &
          +zz(6)*vcm(ix,iy,iz-6) &
          +zz(5)*vcm(ix,iy,iz-5) &
          +zz(4)*vcm(ix,iy,iz-4) &
          +zz(3)*vcm(ix,iy,iz-3) &
          +zz(2)*vcm(ix,iy,iz-2) &
          +zz(1)*vcm(ix,iy,iz-1) &
          +zz(1)*vcm(ix,iy,iz+1) &
          +zz(2)*vcm(ix,iy,iz+2) &
          +zz(3)*vcm(ix,iy,iz+3) &
          +zz(4)*vcm(ix,iy,iz+4) &
          +zz(5)*vcm(ix,iy,iz+5) &
          +zz(6)*vcm(ix,iy,iz+6) &
          +zz(7)*vcm(ix,iy,iz+7) &
          +zz(8)*vcm(ix,iy,iz+8)
      end do
      end do
      end do
!$omp end do nowait
      end if
!     ======================================================
      deallocate(xx,yy,zz)
      return
      end subroutine kslaplacian_c_01


      subroutine kslaplacian_c_02(ncpx,ncpy,ncpz,nf,vcm,veff,avcm)
      implicit none
      integer,   intent(in)::ncpx,ncpy,ncpz,nf
      complex*16,intent(in)::vcm(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf,2)
      real*8,    intent(in)::veff(ncpx,ncpy,ncpz,4)
      complex*16,intent(inout)::avcm(ncpx,ncpy,ncpz,2)
      integer ix,iy,iz

!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        avcm(ix,iy,iz,1)= avcm(ix,iy,iz,1) +dcmplx(veff(ix,iy,iz,3),-veff(ix,iy,iz,4))*vcm(ix,iy,iz,2)
        avcm(ix,iy,iz,2)= avcm(ix,iy,iz,2) +dcmplx(veff(ix,iy,iz,3), veff(ix,iy,iz,4))*vcm(ix,iy,iz,1)
      end do
      end do
      end do

      end subroutine kslaplacian_c_02


      end module

