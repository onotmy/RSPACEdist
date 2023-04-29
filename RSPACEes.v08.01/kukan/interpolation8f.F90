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
!     **********  interpolation8f.f90 07/01/2011-01  **********

      module mod_interpolation
      implicit none

      contains

      subroutine interpolation(ncpx,ncpy,ncpz,nmesh,v_dd,v_cc,nperi,nf,is,ndisp)
      use mod_mpi, only: myrx,myry,myrz,nprocx,nprocy,nprocz
      implicit none
      integer nmesh,ncpx,ncpy,ncpz,nperi,nf,ndisp,is
      integer ncpx_d,ncpy_d,ncpz_d
      integer nf1,nf2,nf3
      integer i,j,k,ix,iy,iz,jx,jy,jz,kx,ky,kz,lx,ly,lz
      real*8 tmp,tmp1,r
      real*8 v_cc(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf)
      real*8 v_dd(-(nf*nmesh-1):ncpx*nmesh+nf*nmesh,-(nf*nmesh-1):ncpy*nmesh+nf*nmesh,-(nf*nmesh-1):ncpz*nmesh+nf*nmesh)
      real*8,allocatable::aweight(:,:)
      real*8,allocatable::vtmp1(:,:,:)
      real*8,allocatable::vtmp2(:,:,:)
      ncpx_d=ncpx*nmesh
      ncpy_d=ncpy*nmesh
      ncpz_d=ncpz*nmesh
      nf1=nf-1
      nf2=nf*nmesh
      nf3=nf2-1
      allocate(aweight(-nf+1:nf,0:nmesh-1))
      allocate(vtmp1(-nf3:ncpx_d+nf2,-nf1:ncpy+nf,-nf1:ncpz+nf))
      allocate(vtmp2(-nf3:ncpx_d+nf2,-nf3:ncpy_d+nf2,-nf1:ncpz+nf))

      if ((is .ne. 1) .and. (is .ne. 2)) then
        write(ndisp,*) 'error in interpolation! IS should be 1 or 2!'
        stop
      end if

      aweight=1.0d0
      tmp=1.0d0/nmesh
      tmp1=0.5d0
      if (mod(nmesh,2) .eq. 1) tmp1=0.0d0
      do k=0,nmesh-1
        r=(dfloat(k)+tmp1)*tmp
        do j=-nf1,nf
          do i=-nf1,nf
            if (i .ne. j) aweight(j,k)=aweight(j,k)*(r-dfloat(i))/(dfloat(j)-dfloat(i))
          end do
        end do
      end do

!     **********  forward interpolation  **********
      if (is .eq. 1) then
      v_dd=0.0d0
      vtmp1=0.0d0
      vtmp2=0.0d0
!$omp parallel default(shared) private(jx,jy,jz,ix,iy,iz,kx,ky,kz,lx,ly,lz)
!cdir nodep
!$omp do
      do jz=-nf1,ncpz+nf
      do jy=-nf1,ncpy+nf
      do jx=1,ncpx_d
        do ix=-nf1,nf
          kx=mod(jx+(nmesh-1)/2,nmesh)
          lx=(jx+(nmesh-1)/2)/nmesh
          vtmp1(jx,jy,jz)=vtmp1(jx,jy,jz)+aweight(ix,kx)*v_cc(lx+ix,jy,jz)
        end do
      end do
      end do
      end do
!cdir nodep
!$omp do
      do jz=-nf1,ncpz+nf
      do jy=1,ncpy_d
      do jx=1,ncpx_d
        do iy=-nf1,nf
          ky=mod(jy+(nmesh-1)/2,nmesh)
          ly=(jy+(nmesh-1)/2)/nmesh
          vtmp2(jx,jy,jz)=vtmp2(jx,jy,jz)+aweight(iy,ky)*vtmp1(jx,ly+iy,jz)
        end do
      end do
      end do
      end do
!cdir nodep
!$omp do
      do jz=1,ncpz_d
      do jy=1,ncpy_d
      do jx=1,ncpx_d
        do iz=-nf1,nf
          kz=mod(jz+(nmesh-1)/2,nmesh)
          lz=(jz+(nmesh-1)/2)/nmesh
          v_dd(jx,jy,jz)=v_dd(jx,jy,jz)+aweight(iz,kz)*vtmp2(jx,jy,lz+iz)
        end do
      end do
      end do
      end do
!$omp end parallel
      end if ! (is .eq. 1)
!     *********************************************

!     **********  backward interpolation  **********
      if (is .eq. 2) then
      v_cc=0.0d0
      vtmp1=0.0d0
      vtmp2=0.0d0
!$omp parallel default(shared) private(jx,jy,jz,ix,iy,iz,kx,ky,kz,lx,ly,lz)
      if ((myrx .ne. 0) .or. (nperi .gt. 0)) then
!$omp do
        do jz=-nf3,ncpz_d+nf2
        do jy=-nf3,ncpy_d+nf2
        do jx=-nf3,0
          v_dd(jx,jy,jz)=0.0d0
        end do
        end do
        end do
      end if
      if ((myrx .ne. nprocx-1) .or. (nperi .gt. 0)) then
!$omp do
        do jz=-nf3,ncpz_d+nf2
        do jy=-nf3,ncpy_d+nf2
        do jx=ncpx_d+1,ncpx_d+nf2
          v_dd(jx,jy,jz)=0.0d0
        end do
        end do
        end do
      end if
      if ((myry .ne. 0) .or. (nperi .gt. 1)) then
!$omp do
        do jz=-nf3,ncpz_d+nf2
        do jy=-nf3,0
        do jx=-nf3,ncpx_d+nf2
          v_dd(jx,jy,jz)=0.0d0
        end do
        end do
        end do
      end if
      if ((myry .ne. nprocy-1) .or. (nperi .gt. 1)) then
!$omp do
        do jz=-nf3,ncpz_d+nf2
        do jy=ncpy_d+1,ncpy_d+nf2
        do jx=-nf3,ncpx_d+nf2
          v_dd(jx,jy,jz)=0.0d0
        end do
        end do
        end do
      end if
      if ((myrz .ne. 0) .or. (nperi .gt. 2)) then
!$omp do
        do jz=-nf3,0
        do jy=-nf3,ncpy_d+nf2
        do jx=-nf3,ncpx_d+nf2
          v_dd(jx,jy,jz)=0.0d0
        end do
        end do
        end do
      end if
      if ((myrz .ne. nprocz-1) .or. (nperi .gt. 2)) then
!$omp do
        do jz=ncpz_d+1,ncpz_d+nf2
        do jy=-nf3,ncpy_d+nf2
        do jx=-nf3,ncpx_d+nf2
          v_dd(jx,jy,jz)=0.0d0
        end do
        end do
        end do
      end if

      do iz=-nf1,nf
        do jz=-nf3,ncpz_d+nf2
        lz=(jz+nf*nmesh+(nmesh-1)/2)/nmesh-nf
!         check whether the computed coarse grid point is in the subdomain
          if ((lz+iz+nf1)*(ncpz+nf-(lz+iz)) .ge. 0) then
            kz=mod(jz+nf*nmesh+(nmesh-1)/2,nmesh)
!cdir nodep
!$omp do
            do jy=-nf3,ncpy_d+nf2
              do jx=-nf3,ncpx_d+nf2
                vtmp2(jx,jy,lz+iz)=vtmp2(jx,jy,lz+iz)+aweight(iz,kz)*v_dd(jx,jy,jz)
              end do
            end do
          end if
        end do
      end do

      do iy=-nf1,nf
!$omp do
        do jz=-nf1,ncpz+nf
          do jy=-nf3,ncpy_d+nf2
            ly=(jy+nf*nmesh+(nmesh-1)/2)/nmesh-nf
!           check whether the computed coarse grid point is in the subdomain
            if ((ly+iy+nf1)*(ncpy+nf-(ly+iy)) .ge. 0) then
              ky=mod(jy+nf*nmesh+(nmesh-1)/2,nmesh)
!cdir nodep
              do jx=-nf3,ncpx_d+nf2
                vtmp1(jx,ly+iy,jz)=vtmp1(jx,ly+iy,jz)+aweight(iy,ky)*vtmp2(jx,jy,jz)
              end do
            end if
          end do
        end do
      end do

      do ix=-nf1,nf
!$omp do
        do jz=-nf1,ncpz+nf
          do jy=-nf1,ncpy+nf
            do jx=-nf3,ncpx_d+nf2
              lx=(jx+nf*nmesh+(nmesh-1)/2)/nmesh-nf
!             check whether the computed coarse grid point is in the subdomain
              if ((lx+ix+nf1)*(ncpx+nf-(lx+ix)) .ge. 0) then
                kx=mod(jx+nf*nmesh+(nmesh-1)/2,nmesh)
                v_cc(lx+ix,jy,jz)=v_cc(lx+ix,jy,jz)+aweight(ix,kx)*vtmp1(jx,jy,jz)
              end if
            end do
          end do
        end do
      end do
!$omp end parallel

      v_cc=v_cc/nmesh**3
      end if ! (is .eq. 2)
!     **********************************************

      deallocate(aweight,vtmp1,vtmp2)
      return
      end subroutine


      end module
