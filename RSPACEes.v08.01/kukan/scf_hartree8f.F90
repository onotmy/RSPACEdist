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
! **********  scf_hartree8f.F90 06/12/2022-01  **********

module mod_scf_hartree
use mod_mpi
use mod_overlap_finitedifference, only:overlap_finitedifference_init,overlap_finitedifference_final, &
                                       overlap_finitedifference_r,overlap_fdcheck_r
implicit none

contains

subroutine scf_hartree( &
 nd,ncpx_d,ncpy_d,ncpz_d,nfh,ndisp,nperi,ncgres,ncgmin,ncgmax,xmax,ymax,zmax,epsvh, & ! <
 rho_aug_dense,boundx,boundy,boundz,                                                & ! <
 vh_dense)                                                                            ! X
use mod_tools, only:tools_definefdcoef
implicit none
integer,intent(in)   ::nd
integer,intent(in)   ::ncpx_d,ncpy_d,ncpz_d
integer,intent(in)   ::nfh,ndisp,nperi
integer,intent(in)   ::ncgres,ncgmin,ncgmax
real*8, intent(in)   ::xmax,ymax,zmax
real*8, intent(in)   ::epsvh
real*8, intent(in)   ::rho_aug_dense(ncpx_d,ncpy_d,ncpz_d)
real*8, intent(in)   ::boundx(-(nfh-1)*(1-nd)+nd:nfh*(1-nd)+nd,ncpy_d*(1-nd)+nd,ncpz_d*(1-nd)+nd)
real*8, intent(in)   ::boundy(ncpx_d*(1-nd)+nd,-(nfh-1)*(1-nd)+nd:nfh*(1-nd)+nd,ncpz_d*(1-nd)+nd)
real*8, intent(in)   ::boundz(ncpx_d*(1-nd)+nd,ncpy_d*(1-nd)+nd,-(nfh-1)*(1-nd)+nd:nfh*(1-nd)+nd)
real*8, intent(inout)::vh_dense(ncpx_d,ncpy_d,ncpz_d)

integer ix,iy,iz
real*8 pi
real*8 ddx,ddy,ddz
real*8 drtx1,drtx2,drtx3,drtx4,drtx5,drtx6
real*8 drty1,drty2,drty3,drty4,drty5,drty6
real*8 drtz1,drtz2,drtz3,drtz4,drtz5,drtz6

real*8,allocatable:: acfdh(:)
real*8,allocatable:: v_dense(:,:,:) &
                   ,av_dense(:,:,:) &
                 ,bvec_dense(:,:,:) &
                 ,rvec_dense(:,:,:) &
                 ,pvec_dense(:,:,:) &
                 ,avhs_dense(:,:,:)
  allocate(acfdh(0:8))
  allocate(v_dense(-(nfh-1):ncpx_d+nfh,-(nfh-1):ncpy_d+nfh,-(nfh-1):ncpz_d+nfh))
  allocate(  av_dense(ncpx_d,ncpy_d,ncpz_d) &
          ,bvec_dense(ncpx_d,ncpy_d,ncpz_d) &
          ,rvec_dense(ncpx_d,ncpy_d,ncpz_d) &
          ,pvec_dense(ncpx_d,ncpy_d,ncpz_d) &
          ,avhs_dense(ncpx_d,ncpy_d,ncpz_d))

  pi=dacos(-1.0d0)
  ddx=2.0d0*xmax/nprocx/ncpx_d
  ddy=2.0d0*ymax/nprocy/ncpy_d
  ddz=2.0d0*zmax/nprocz/ncpz_d

  call overlap_finitedifference_init(ncpx_d,ncpy_d,ncpz_d,nfh,1,0)
  call tools_definefdcoef(nfh,acfdh(0),acfdh(1),acfdh(2),acfdh(3),acfdh(4),acfdh(5),acfdh(6),acfdh(7),acfdh(8))

  if ((nfh .ne. 4) .and. (nfh .ne. 5) .and. (nfh .ne. 6)) then
    if (myrank_glbl .eq. 0) write(ndisp,*) 'error in scf_hartree! This subroutine is for nfh ne 4, 5 or 6.'
      call mpi_abort(mpi_comm_world,mpij)
    stop
  end if

  drtx1=acfdh(1)/ddx/ddx*(-0.5d0)
  drty1=acfdh(1)/ddy/ddy*(-0.5d0)
  drtz1=acfdh(1)/ddz/ddz*(-0.5d0)
  drtx2=acfdh(2)/ddx/ddx*(-0.5d0)
  drty2=acfdh(2)/ddy/ddy*(-0.5d0)
  drtz2=acfdh(2)/ddz/ddz*(-0.5d0)
  drtx3=acfdh(3)/ddx/ddx*(-0.5d0)
  drty3=acfdh(3)/ddy/ddy*(-0.5d0)
  drtz3=acfdh(3)/ddz/ddz*(-0.5d0)
  drtx4=acfdh(4)/ddx/ddx*(-0.5d0)
  drty4=acfdh(4)/ddy/ddy*(-0.5d0)
  drtz4=acfdh(4)/ddz/ddz*(-0.5d0)
  drtx5=acfdh(5)/ddx/ddx*(-0.5d0)
  drty5=acfdh(5)/ddy/ddy*(-0.5d0)
  drtz5=acfdh(5)/ddz/ddz*(-0.5d0)
  drtx6=acfdh(6)/ddx/ddx*(-0.5d0)
  drty6=acfdh(6)/ddy/ddy*(-0.5d0)
  drtz6=acfdh(6)/ddz/ddz*(-0.5d0)

!$omp parallel default(shared) private(ix)
!$omp do
  do ix=1,ncpx_d*ncpy_d*ncpz_d
    bvec_dense(ix,1,1)=2.0d0*pi*rho_aug_dense(ix,1,1)
  end do
!$omp end parallel

! ==========  set the boundary for nfh=4  ==========
  if (nfh .eq. 4) then
  if (nperi .le. 2) then
  if (myrz .eq. 0) then
  do iy=1,ncpy_d
  do ix=1,ncpx_d
        bvec_dense(ix,iy,1)=bvec_dense(ix,iy,1)   &
             -drtz1*(boundz(ix,iy, 0)) &
             -drtz2*(boundz(ix,iy,-1)) &
             -drtz3*(boundz(ix,iy,-2)) &
             -drtz4*(boundz(ix,iy,-3))
        bvec_dense(ix,iy,2)=bvec_dense(ix,iy,2)   &
             -drtz2*(boundz(ix,iy, 0)) &
             -drtz3*(boundz(ix,iy,-1)) &
             -drtz4*(boundz(ix,iy,-2))
        bvec_dense(ix,iy,3)=bvec_dense(ix,iy,3)   &
             -drtz3*(boundz(ix,iy, 0)) &
             -drtz4*(boundz(ix,iy,-1))
        bvec_dense(ix,iy,4)=bvec_dense(ix,iy,4)   &
             -drtz4*(boundz(ix,iy, 0))
  end do
  end do
  end if
  if (myrz .eq. nprocz-1) then
  do iy=1,ncpy_d
  do ix=1,ncpx_d
        bvec_dense(ix,iy,ncpz_d  )=bvec_dense(ix,iy,ncpz_d  ) &
             -drtz1*(boundz(ix,iy,1))     &
             -drtz2*(boundz(ix,iy,2))     &
             -drtz3*(boundz(ix,iy,3))     &
             -drtz4*(boundz(ix,iy,4))
        bvec_dense(ix,iy,ncpz_d-1)=bvec_dense(ix,iy,ncpz_d-1) &
             -drtz2*(boundz(ix,iy,1))     &
             -drtz3*(boundz(ix,iy,2))     &
             -drtz4*(boundz(ix,iy,3))
        bvec_dense(ix,iy,ncpz_d-2)=bvec_dense(ix,iy,ncpz_d-2) &
             -drtz3*(boundz(ix,iy,1))     &
             -drtz4*(boundz(ix,iy,2))
        bvec_dense(ix,iy,ncpz_d-3)=bvec_dense(ix,iy,ncpz_d-3) &
             -drtz4*(boundz(ix,iy,1))
  end do
  end do
  end if
  end if
  if (nperi .le. 1) then
  if (myry .eq. 0) then
  do iz=1,ncpz_d
  do ix=1,ncpx_d
        bvec_dense(ix,1,iz)=bvec_dense(ix,1,iz)   &
             -drty1*(boundy(ix, 0,iz)) &
             -drty2*(boundy(ix,-1,iz)) &
             -drty3*(boundy(ix,-2,iz)) &
             -drty4*(boundy(ix,-3,iz))
        bvec_dense(ix,2,iz)=bvec_dense(ix,2,iz)   &
             -drty2*(boundy(ix, 0,iz)) &
             -drty3*(boundy(ix,-1,iz)) &
             -drty4*(boundy(ix,-2,iz))
        bvec_dense(ix,3,iz)=bvec_dense(ix,3,iz)   &
             -drty3*(boundy(ix, 0,iz)) &
             -drty4*(boundy(ix,-1,iz))
        bvec_dense(ix,4,iz)=bvec_dense(ix,4,iz)   &
             -drty4*(boundy(ix, 0,iz))
  end do
  end do
  end if
  if (myry .eq. nprocy-1) then
  do iz=1,ncpz_d
  do ix=1,ncpx_d
        bvec_dense(ix,ncpy_d  ,iz)=bvec_dense(ix,ncpy_d  ,iz) &
             -drty1*(boundy(ix,1,iz))     &
             -drty2*(boundy(ix,2,iz))     &
             -drty3*(boundy(ix,3,iz))     &
             -drty4*(boundy(ix,4,iz))
        bvec_dense(ix,ncpy_d-1,iz)=bvec_dense(ix,ncpy_d-1,iz) &
             -drty2*(boundy(ix,1,iz))     &
             -drty3*(boundy(ix,2,iz))     &
             -drty4*(boundy(ix,3,iz))
        bvec_dense(ix,ncpy_d-2,iz)=bvec_dense(ix,ncpy_d-2,iz) &
             -drty3*(boundy(ix,1,iz))     &
             -drty4*(boundy(ix,2,iz))
        bvec_dense(ix,ncpy_d-3,iz)=bvec_dense(ix,ncpy_d-3,iz) &
             -drty4*(boundy(ix,1,iz))
  end do
  end do
  end if
  end if
  if (nperi .eq. 0) then
  if (myrx .eq. 0) then
  do iz=1,ncpz_d
  do iy=1,ncpy_d
        bvec_dense(1,iy,iz)=bvec_dense(1,iy,iz)   &
             -drtx1*(boundx( 0,iy,iz)) &
             -drtx2*(boundx(-1,iy,iz)) &
             -drtx3*(boundx(-2,iy,iz)) &
             -drtx4*(boundx(-3,iy,iz))
        bvec_dense(2,iy,iz)=bvec_dense(2,iy,iz)   &
             -drtx2*(boundx( 0,iy,iz)) &
             -drtx3*(boundx(-1,iy,iz)) &
             -drtx4*(boundx(-2,iy,iz))
        bvec_dense(3,iy,iz)=bvec_dense(3,iy,iz)   &
             -drtx3*(boundx( 0,iy,iz)) &
             -drtx4*(boundx(-1,iy,iz))
        bvec_dense(4,iy,iz)=bvec_dense(4,iy,iz)   &
             -drtx4*(boundx( 0,iy,iz))
  end do
  end do
  end if
  if (myrx .eq. nprocx-1) then
  do iz=1,ncpz_d
  do iy=1,ncpy_d
        bvec_dense(ncpx_d  ,iy,iz)=bvec_dense(ncpx_d  ,iy,iz) &
             -drtx1*(boundx(1,iy,iz))     &
             -drtx2*(boundx(2,iy,iz))     &
             -drtx3*(boundx(3,iy,iz))     &
             -drtx4*(boundx(4,iy,iz))
        bvec_dense(ncpx_d-1,iy,iz)=bvec_dense(ncpx_d-1,iy,iz) &
             -drtx2*(boundx(1,iy,iz))     &
             -drtx3*(boundx(2,iy,iz))     &
             -drtx4*(boundx(3,iy,iz))
        bvec_dense(ncpx_d-2,iy,iz)=bvec_dense(ncpx_d-2,iy,iz) &
             -drtx3*(boundx(1,iy,iz))     &
             -drtx4*(boundx(2,iy,iz))
        bvec_dense(ncpx_d-3,iy,iz)=bvec_dense(ncpx_d-3,iy,iz) &
             -drtx4*(boundx(1,iy,iz))
  end do
  end do
  end if
  end if
  end if
! ==================================================

! ==========  set the boundary for nfh=5  ==========
  if (nfh .eq. 5) then
  if (nperi .le. 2) then
  if (myrz .eq. 0) then
  do iy=1,ncpy_d
  do ix=1,ncpx_d
        bvec_dense(ix,iy,1)=bvec_dense(ix,iy,1)   &
             -drtz1*(boundz(ix,iy, 0)) &
             -drtz2*(boundz(ix,iy,-1)) &
             -drtz3*(boundz(ix,iy,-2)) &
             -drtz4*(boundz(ix,iy,-3)) &
             -drtz5*(boundz(ix,iy,-4))
        bvec_dense(ix,iy,2)=bvec_dense(ix,iy,2)   &
             -drtz2*(boundz(ix,iy, 0)) &
             -drtz3*(boundz(ix,iy,-1)) &
             -drtz4*(boundz(ix,iy,-2)) &
             -drtz5*(boundz(ix,iy,-3))
        bvec_dense(ix,iy,3)=bvec_dense(ix,iy,3)   &
             -drtz3*(boundz(ix,iy, 0)) &
             -drtz4*(boundz(ix,iy,-1)) &
             -drtz5*(boundz(ix,iy,-2))
        bvec_dense(ix,iy,4)=bvec_dense(ix,iy,4)   &
             -drtz4*(boundz(ix,iy, 0)) &
             -drtz5*(boundz(ix,iy,-1))
        bvec_dense(ix,iy,5)=bvec_dense(ix,iy,5)   &
             -drtz5*(boundz(ix,iy, 0))
  end do
  end do
  end if
  if (myrz .eq. nprocz-1) then
  do iy=1,ncpy_d
  do ix=1,ncpx_d
        bvec_dense(ix,iy,ncpz_d  )=bvec_dense(ix,iy,ncpz_d  ) &
             -drtz1*(boundz(ix,iy,1)) &
             -drtz2*(boundz(ix,iy,2)) &
             -drtz3*(boundz(ix,iy,3)) &
             -drtz4*(boundz(ix,iy,4)) &
             -drtz5*(boundz(ix,iy,5))
        bvec_dense(ix,iy,ncpz_d-1)=bvec_dense(ix,iy,ncpz_d-1) &
             -drtz2*(boundz(ix,iy,1)) &
             -drtz3*(boundz(ix,iy,2)) &
             -drtz4*(boundz(ix,iy,3)) &
             -drtz5*(boundz(ix,iy,4))
        bvec_dense(ix,iy,ncpz_d-2)=bvec_dense(ix,iy,ncpz_d-2) &
             -drtz3*(boundz(ix,iy,1)) &
             -drtz4*(boundz(ix,iy,2)) &
             -drtz5*(boundz(ix,iy,3))
        bvec_dense(ix,iy,ncpz_d-3)=bvec_dense(ix,iy,ncpz_d-3) &
             -drtz4*(boundz(ix,iy,1)) &
             -drtz5*(boundz(ix,iy,2))
        bvec_dense(ix,iy,ncpz_d-4)=bvec_dense(ix,iy,ncpz_d-4) &
             -drtz5*(boundz(ix,iy,1))
  end do
  end do
  end if
  end if
  if (nperi .le. 1) then
  if (myry .eq. 0) then
  do iz=1,ncpz_d
  do ix=1,ncpx_d
        bvec_dense(ix,1,iz)=bvec_dense(ix,1,iz)   &
             -drty1*(boundy(ix, 0,iz)) &
             -drty2*(boundy(ix,-1,iz)) &
             -drty3*(boundy(ix,-2,iz)) &
             -drty4*(boundy(ix,-3,iz)) &
             -drty5*(boundy(ix,-4,iz))
        bvec_dense(ix,2,iz)=bvec_dense(ix,2,iz)   &
             -drty2*(boundy(ix, 0,iz)) &
             -drty3*(boundy(ix,-1,iz)) &
             -drty4*(boundy(ix,-2,iz)) &
             -drty5*(boundy(ix,-3,iz))
        bvec_dense(ix,3,iz)=bvec_dense(ix,3,iz)   &
             -drty3*(boundy(ix, 0,iz)) &
             -drty4*(boundy(ix,-1,iz)) &
             -drty5*(boundy(ix,-2,iz))
        bvec_dense(ix,4,iz)=bvec_dense(ix,4,iz)   &
             -drty4*(boundy(ix, 0,iz)) &
             -drty5*(boundy(ix,-1,iz))
        bvec_dense(ix,5,iz)=bvec_dense(ix,5,iz)   &
             -drty5*(boundy(ix, 0,iz))
  end do
  end do
  end if
  if (myry .eq. nprocy-1) then
  do iz=1,ncpz_d
  do ix=1,ncpx_d
        bvec_dense(ix,ncpy_d  ,iz)=bvec_dense(ix,ncpy_d  ,iz) &
             -drty1*(boundy(ix,1,iz)) &
             -drty2*(boundy(ix,2,iz)) &
             -drty3*(boundy(ix,3,iz)) &
             -drty4*(boundy(ix,4,iz)) &
             -drty5*(boundy(ix,5,iz))
        bvec_dense(ix,ncpy_d-1,iz)=bvec_dense(ix,ncpy_d-1,iz) &
             -drty2*(boundy(ix,1,iz)) &
             -drty3*(boundy(ix,2,iz)) &
             -drty4*(boundy(ix,3,iz)) &
             -drty5*(boundy(ix,4,iz))
        bvec_dense(ix,ncpy_d-2,iz)=bvec_dense(ix,ncpy_d-2,iz) &
             -drty3*(boundy(ix,1,iz)) &
             -drty4*(boundy(ix,2,iz)) &
             -drty5*(boundy(ix,3,iz))
        bvec_dense(ix,ncpy_d-3,iz)=bvec_dense(ix,ncpy_d-3,iz) &
             -drty4*(boundy(ix,1,iz)) &
             -drty5*(boundy(ix,2,iz))
        bvec_dense(ix,ncpy_d-4,iz)=bvec_dense(ix,ncpy_d-4,iz) &
             -drty5*(boundy(ix,1,iz))
  end do
  end do
  end if
  end if
  if (nperi .eq. 0) then
  if (myrx .eq. 0) then
  do iz=1,ncpz_d
  do iy=1,ncpy_d
        bvec_dense(1,iy,iz)=bvec_dense(1,iy,iz)   &
             -drtx1*(boundx( 0,iy,iz)) &
             -drtx2*(boundx(-1,iy,iz)) &
             -drtx3*(boundx(-2,iy,iz)) &
             -drtx4*(boundx(-3,iy,iz)) &
             -drtx5*(boundx(-4,iy,iz))
        bvec_dense(2,iy,iz)=bvec_dense(2,iy,iz)   &
             -drtx2*(boundx( 0,iy,iz)) &
             -drtx3*(boundx(-1,iy,iz)) &
             -drtx4*(boundx(-2,iy,iz)) &
             -drtx5*(boundx(-3,iy,iz))
        bvec_dense(3,iy,iz)=bvec_dense(3,iy,iz)   &
             -drtx3*(boundx( 0,iy,iz)) &
             -drtx4*(boundx(-1,iy,iz)) &
             -drtx5*(boundx(-2,iy,iz))
        bvec_dense(4,iy,iz)=bvec_dense(4,iy,iz)   &
             -drtx4*(boundx( 0,iy,iz)) &
             -drtx5*(boundx(-1,iy,iz))
        bvec_dense(5,iy,iz)=bvec_dense(5,iy,iz)   &
             -drtx5*(boundx( 0,iy,iz))
  end do
  end do
  end if
  if (myrx .eq. nprocx-1) then
  do iz=1,ncpz_d
  do iy=1,ncpy_d
        bvec_dense(ncpx_d  ,iy,iz)=bvec_dense(ncpx_d  ,iy,iz) &
             -drtx1*(boundx(1,iy,iz)) &
             -drtx2*(boundx(2,iy,iz)) &
             -drtx3*(boundx(3,iy,iz)) &
             -drtx4*(boundx(4,iy,iz)) &
             -drtx5*(boundx(5,iy,iz))
        bvec_dense(ncpx_d-1,iy,iz)=bvec_dense(ncpx_d-1,iy,iz) &
             -drtx2*(boundx(1,iy,iz)) &
             -drtx3*(boundx(2,iy,iz)) &
             -drtx4*(boundx(3,iy,iz)) &
             -drtx5*(boundx(4,iy,iz))
        bvec_dense(ncpx_d-2,iy,iz)=bvec_dense(ncpx_d-2,iy,iz) &
             -drtx3*(boundx(1,iy,iz)) &
             -drtx4*(boundx(2,iy,iz)) &
             -drtx5*(boundx(3,iy,iz))
        bvec_dense(ncpx_d-3,iy,iz)=bvec_dense(ncpx_d-3,iy,iz) &
             -drtx4*(boundx(1,iy,iz)) &
             -drtx5*(boundx(2,iy,iz))
        bvec_dense(ncpx_d-4,iy,iz)=bvec_dense(ncpx_d-4,iy,iz) &
             -drtx5*(boundx(1,iy,iz))
  end do
  end do
  end if
  end if
  end if
! =================================================

! ==========  set the boundary for nfh=6  ==========
  if (nfh .eq. 6) then
  if (nperi .le. 2) then
  if (myrz .eq. 0) then
  do iy=1,ncpy_d
  do ix=1,ncpx_d
        bvec_dense(ix,iy,1)=bvec_dense(ix,iy,1)   &
             -drtz1*(boundz(ix,iy, 0)) &
             -drtz2*(boundz(ix,iy,-1)) &
             -drtz3*(boundz(ix,iy,-2)) &
             -drtz4*(boundz(ix,iy,-3)) &
             -drtz5*(boundz(ix,iy,-4)) &
             -drtz6*(boundz(ix,iy,-5))
        bvec_dense(ix,iy,2)=bvec_dense(ix,iy,2)   &
             -drtz2*(boundz(ix,iy, 0)) &
             -drtz3*(boundz(ix,iy,-1)) &
             -drtz4*(boundz(ix,iy,-2)) &
             -drtz5*(boundz(ix,iy,-3)) &
             -drtz6*(boundz(ix,iy,-4))
        bvec_dense(ix,iy,3)=bvec_dense(ix,iy,3)   &
             -drtz3*(boundz(ix,iy, 0)) &
             -drtz4*(boundz(ix,iy,-1)) &
             -drtz5*(boundz(ix,iy,-2)) &
             -drtz6*(boundz(ix,iy,-3))
        bvec_dense(ix,iy,4)=bvec_dense(ix,iy,4)   &
             -drtz4*(boundz(ix,iy, 0)) &
             -drtz5*(boundz(ix,iy,-1)) &
             -drtz6*(boundz(ix,iy,-2))
        bvec_dense(ix,iy,5)=bvec_dense(ix,iy,5)   &
             -drtz5*(boundz(ix,iy, 0)) &
             -drtz6*(boundz(ix,iy,-1))
        bvec_dense(ix,iy,6)=bvec_dense(ix,iy,6)   &
             -drtz6*(boundz(ix,iy, 0))
  end do
  end do
  end if
  if (myrz .eq. nprocz-1) then
  do iy=1,ncpy_d
  do ix=1,ncpx_d
        bvec_dense(ix,iy,ncpz_d  )=bvec_dense(ix,iy,ncpz_d  ) &
             -drtz1*(boundz(ix,iy,1)) &
             -drtz2*(boundz(ix,iy,2)) &
             -drtz3*(boundz(ix,iy,3)) &
             -drtz4*(boundz(ix,iy,4)) &
             -drtz5*(boundz(ix,iy,5)) &
             -drtz6*(boundz(ix,iy,6))
        bvec_dense(ix,iy,ncpz_d-1)=bvec_dense(ix,iy,ncpz_d-1) &
             -drtz2*(boundz(ix,iy,1)) &
             -drtz3*(boundz(ix,iy,2)) &
             -drtz4*(boundz(ix,iy,3)) &
             -drtz5*(boundz(ix,iy,4)) &
             -drtz6*(boundz(ix,iy,5))
        bvec_dense(ix,iy,ncpz_d-2)=bvec_dense(ix,iy,ncpz_d-2) &
             -drtz3*(boundz(ix,iy,1)) &
             -drtz4*(boundz(ix,iy,2)) &
             -drtz5*(boundz(ix,iy,3)) &
             -drtz6*(boundz(ix,iy,4))
        bvec_dense(ix,iy,ncpz_d-3)=bvec_dense(ix,iy,ncpz_d-3) &
             -drtz4*(boundz(ix,iy,1)) &
             -drtz5*(boundz(ix,iy,2)) &
             -drtz6*(boundz(ix,iy,3))
        bvec_dense(ix,iy,ncpz_d-4)=bvec_dense(ix,iy,ncpz_d-4) &
             -drtz5*(boundz(ix,iy,1)) &
             -drtz6*(boundz(ix,iy,2))
        bvec_dense(ix,iy,ncpz_d-5)=bvec_dense(ix,iy,ncpz_d-5) &
             -drtz6*(boundz(ix,iy,1))
  end do
  end do
  end if
  end if
  if (nperi .le. 1) then
  if (myry .eq. 0) then
  do iz=1,ncpz_d
  do ix=1,ncpx_d
        bvec_dense(ix,1,iz)=bvec_dense(ix,1,iz)   &
             -drty1*(boundy(ix, 0,iz)) &
             -drty2*(boundy(ix,-1,iz)) &
             -drty3*(boundy(ix,-2,iz)) &
             -drty4*(boundy(ix,-3,iz)) &
             -drty5*(boundy(ix,-4,iz)) &
             -drty6*(boundy(ix,-5,iz))
        bvec_dense(ix,2,iz)=bvec_dense(ix,2,iz)   &
             -drty2*(boundy(ix, 0,iz)) &
             -drty3*(boundy(ix,-1,iz)) &
             -drty4*(boundy(ix,-2,iz)) &
             -drty5*(boundy(ix,-3,iz)) &
             -drty6*(boundy(ix,-4,iz))
        bvec_dense(ix,3,iz)=bvec_dense(ix,3,iz)   &
             -drty3*(boundy(ix, 0,iz)) &
             -drty4*(boundy(ix,-1,iz)) &
             -drty5*(boundy(ix,-2,iz)) &
             -drty6*(boundy(ix,-3,iz))
        bvec_dense(ix,4,iz)=bvec_dense(ix,4,iz)   &
             -drty4*(boundy(ix, 0,iz)) &
             -drty5*(boundy(ix,-1,iz)) &
             -drty6*(boundy(ix,-2,iz))
        bvec_dense(ix,5,iz)=bvec_dense(ix,5,iz)   &
             -drty5*(boundy(ix, 0,iz)) &
             -drty6*(boundy(ix,-1,iz))
        bvec_dense(ix,6,iz)=bvec_dense(ix,6,iz)   &
             -drty6*(boundy(ix, 0,iz))
  end do
  end do
  end if
  if (myry .eq. nprocy-1) then
  do iz=1,ncpz_d
  do ix=1,ncpx_d
        bvec_dense(ix,ncpy_d  ,iz)=bvec_dense(ix,ncpy_d  ,iz) &
             -drty1*(boundy(ix,1,iz)) &
             -drty2*(boundy(ix,2,iz)) &
             -drty3*(boundy(ix,3,iz)) &
             -drty4*(boundy(ix,4,iz)) &
             -drty5*(boundy(ix,5,iz)) &
             -drty6*(boundy(ix,6,iz))
        bvec_dense(ix,ncpy_d-1,iz)=bvec_dense(ix,ncpy_d-1,iz) &
             -drty2*(boundy(ix,1,iz)) &
             -drty3*(boundy(ix,2,iz)) &
             -drty4*(boundy(ix,3,iz)) &
             -drty5*(boundy(ix,4,iz)) &
             -drty6*(boundy(ix,5,iz))
        bvec_dense(ix,ncpy_d-2,iz)=bvec_dense(ix,ncpy_d-2,iz) &
             -drty3*(boundy(ix,1,iz)) &
             -drty4*(boundy(ix,2,iz)) &
             -drty5*(boundy(ix,3,iz)) &
             -drty6*(boundy(ix,4,iz))
        bvec_dense(ix,ncpy_d-3,iz)=bvec_dense(ix,ncpy_d-3,iz) &
             -drty4*(boundy(ix,1,iz)) &
             -drty5*(boundy(ix,2,iz)) &
             -drty6*(boundy(ix,3,iz))
        bvec_dense(ix,ncpy_d-4,iz)=bvec_dense(ix,ncpy_d-4,iz) &
             -drty5*(boundy(ix,1,iz)) &
             -drty6*(boundy(ix,2,iz))
        bvec_dense(ix,ncpy_d-5,iz)=bvec_dense(ix,ncpy_d-5,iz) &
             -drty6*(boundy(ix,1,iz))
  end do
  end do
  end if
  end if
  if (nperi .eq. 0) then
  if (myrx .eq. 0) then
  do iz=1,ncpz_d
  do iy=1,ncpy_d
        bvec_dense(1,iy,iz)=bvec_dense(1,iy,iz)   &
             -drtx1*(boundx( 0,iy,iz)) &
             -drtx2*(boundx(-1,iy,iz)) &
             -drtx3*(boundx(-2,iy,iz)) &
             -drtx4*(boundx(-3,iy,iz)) &
             -drtx5*(boundx(-4,iy,iz)) &
             -drtx6*(boundx(-5,iy,iz))
        bvec_dense(2,iy,iz)=bvec_dense(2,iy,iz)   &
             -drtx2*(boundx( 0,iy,iz)) &
             -drtx3*(boundx(-1,iy,iz)) &
             -drtx4*(boundx(-2,iy,iz)) &
             -drtx5*(boundx(-3,iy,iz)) &
             -drtx6*(boundx(-4,iy,iz))
        bvec_dense(3,iy,iz)=bvec_dense(3,iy,iz)   &
             -drtx3*(boundx( 0,iy,iz)) &
             -drtx4*(boundx(-1,iy,iz)) &
             -drtx5*(boundx(-2,iy,iz)) &
             -drtx6*(boundx(-3,iy,iz))
        bvec_dense(4,iy,iz)=bvec_dense(4,iy,iz)   &
             -drtx4*(boundx( 0,iy,iz)) &
             -drtx5*(boundx(-1,iy,iz)) &
             -drtx6*(boundx(-2,iy,iz))
        bvec_dense(5,iy,iz)=bvec_dense(5,iy,iz)   &
             -drtx5*(boundx( 0,iy,iz)) &
             -drtx6*(boundx(-1,iy,iz))
        bvec_dense(6,iy,iz)=bvec_dense(6,iy,iz)   &
             -drtx6*(boundx( 0,iy,iz))
  end do
  end do
  end if
  if (myrx .eq. nprocx-1) then
  do iz=1,ncpz_d
  do iy=1,ncpy_d
        bvec_dense(ncpx_d  ,iy,iz)=bvec_dense(ncpx_d  ,iy,iz) &
             -drtx1*(boundx(1,iy,iz)) &
             -drtx2*(boundx(2,iy,iz)) &
             -drtx3*(boundx(3,iy,iz)) &
             -drtx4*(boundx(4,iy,iz)) &
             -drtx5*(boundx(5,iy,iz)) &
             -drtx6*(boundx(6,iy,iz))
        bvec_dense(ncpx_d-1,iy,iz)=bvec_dense(ncpx_d-1,iy,iz) &
             -drtx2*(boundx(1,iy,iz)) &
             -drtx3*(boundx(2,iy,iz)) &
             -drtx4*(boundx(3,iy,iz)) &
             -drtx5*(boundx(4,iy,iz)) &
             -drtx6*(boundx(5,iy,iz))
        bvec_dense(ncpx_d-2,iy,iz)=bvec_dense(ncpx_d-2,iy,iz) &
             -drtx3*(boundx(1,iy,iz)) &
             -drtx4*(boundx(2,iy,iz)) &
             -drtx5*(boundx(3,iy,iz)) &
             -drtx6*(boundx(4,iy,iz))
        bvec_dense(ncpx_d-3,iy,iz)=bvec_dense(ncpx_d-3,iy,iz) &
             -drtx4*(boundx(1,iy,iz)) &
             -drtx5*(boundx(2,iy,iz)) &
             -drtx6*(boundx(3,iy,iz))
        bvec_dense(ncpx_d-4,iy,iz)=bvec_dense(ncpx_d-4,iy,iz) &
             -drtx5*(boundx(1,iy,iz)) &
             -drtx6*(boundx(2,iy,iz))
        bvec_dense(ncpx_d-5,iy,iz)=bvec_dense(ncpx_d-5,iy,iz) &
             -drtx6*(boundx(1,iy,iz))
  end do
  end do
  end if
  end if
  end if
! ==================================================

  call scf_hartree_cg(ncpx_d,ncpy_d,ncpz_d &
         ,nperi,nfh,ncgres,ncgmin,ncgmax &
         ,ndisp &
         ,xmax,ymax,zmax &
         ,epsvh &
         ,vh_dense,bvec_dense,v_dense,pvec_dense &
         ,rvec_dense,av_dense,avhs_dense &
         ,acfdh)

  call overlap_finitedifference_final

  deallocate(acfdh,v_dense,av_dense,bvec_dense,rvec_dense,pvec_dense,avhs_dense)
  return
end subroutine


!     **********  cg7g.f90 05/10/2008-01  **********
!     05/10/2008 routine using kake3d was omitted.
!     04/15/2008 sentence for display was changed.
!     04/06/2008 module vhcalc_arrays was added.
!     04/04/2008 criteria for the conversence of hartree potential was changed.
!     03/01/2008 barrier was added after critical
!     01/20/2008 Comments were added.

!     This file contains
!      + cg
!      + cgpoisson_01-05

subroutine scf_hartree_cg(ncpx_d,ncpy_d,ncpz_d &
             ,nperi,nfh,ncgres,ncgmin,ncgmax &
             ,ndisp &
             ,xmax,ymax,zmax &
             ,epsvh &
             ,vh_dense,bvec_dense,v_dense,pvec_dense &
             ,rvec_dense,av_dense,avhs_dense &
             ,acfdh)

implicit none
integer, intent(in) :: ncpx_d,ncpy_d,ncpz_d
integer, intent(in) :: nperi,nfh,ncgres,ncgmin,ncgmax
integer, intent(in) :: ndisp
real*8, intent(in) :: xmax,ymax,zmax
real*8, intent(in) :: epsvh
real*8 vh_dense(ncpx_d,ncpy_d,ncpz_d)
real*8 v_dense(-(nfh-1):ncpx_d+nfh,-(nfh-1):ncpy_d+nfh,-(nfh-1):ncpz_d+nfh)
real*8 av_dense(ncpx_d,ncpy_d,ncpz_d) &
    ,bvec_dense(ncpx_d,ncpy_d,ncpz_d) &
    ,rvec_dense(ncpx_d,ncpy_d,ncpz_d) &
    ,pvec_dense(ncpx_d,ncpy_d,ncpz_d) &
    ,avhs_dense(ncpx_d,ncpy_d,ncpz_d)
real*8, intent(in) :: acfdh(0:8)

integer ix,iy,iz
integer loopcg
real*8 bve,bveall,bveloc
real*8 r1r,r1rall,r2r,r2rall,dp
real*8 pap,papall
real*8 cgtmp1,cgtmp2
real*8 alph,beta
real*8 ddx,ddy,ddz
logical isopen

  ddx=2.0d0*xmax/nprocx/ncpx_d
  ddy=2.0d0*ymax/nprocy/ncpy_d
  ddz=2.0d0*zmax/nprocz/ncpz_d

  loopcg=0
  bve=0.0d0
!$omp parallel default(shared) private(ix,iy,iz,bveloc)
  if (nperi .eq. 3) then
    bveloc=0.0d0
!$omp do
    do iz=1,ncpz_d
    do iy=1,ncpy_d
    do ix=1,ncpx_d
       bveloc=bveloc+bvec_dense(ix,iy,iz)*ddx*ddy*ddz
    end do
    end do
    end do
!$omp critical
    bve=bve+bveloc
!$omp end critical
!$omp barrier
!$omp single
    call mpi_allreduce(bve,bveall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
    bve=bveall
    bve=bve/(8.0d0*xmax*ymax*zmax)
!$omp end single
!$omp do
    do ix=1,ncpx_d*ncpy_d*ncpz_d
      bvec_dense(ix,1,1)=bvec_dense(ix,1,1)-bve
    end do
!$omp end do nowait
  end if

!$omp do
  do ix=1,(ncpx_d+2*nfh)*(ncpy_d+2*nfh)*(ncpz_d+2*nfh)
    v_dense(-nfh+ix,-nfh+1,-nfh+1)=0.0d0
  end do
!$omp do
  do iz=1,ncpz_d
  do iy=1,ncpy_d
  do ix=1,ncpx_d
     v_dense(ix,iy,iz)=vh_dense(ix,iy,iz)
  end do
  end do
  end do
!$omp single
  call overlap_finitedifference_r(nperi,ncpx_d,ncpy_d,ncpz_d,nfh,nfh-1,nfh,v_dense)
  call overlap_fdcheck_r(nperi,ncpx_d,ncpy_d,ncpz_d,nfh,nfh-1,nfh,v_dense)
!$omp end single
  call scf_hartee_laplacian(ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,v_dense,av_dense,acfdh,nfh,ndisp)
!$omp barrier
!$omp do
  do ix=1,ncpx_d*ncpy_d*ncpz_d
    avhs_dense(ix,1,1)=av_dense(ix,1,1)
  end do
!$omp end do nowait

  call scf_hartee_cg_01(ncpx_d*ncpy_d*ncpz_d,rvec_dense,bvec_dense,av_dense)
!$omp barrier

!$omp do
  do ix=1,ncpx_d*ncpy_d*ncpz_d
    pvec_dense(ix,1,1)=rvec_dense(ix,1,1)
  end do
!$omp end do nowait
!$omp end parallel

 1000 loopcg=loopcg+1
!$omp parallel default(shared) private(ix,iy,iz)
!$omp do
  do iz=1,ncpz_d
  do iy=1,ncpy_d
  do ix=1,ncpx_d
     v_dense(ix,iy,iz)=pvec_dense(ix,iy,iz)
  end do
  end do
  end do

!$omp single
  call overlap_finitedifference_r(nperi,ncpx_d,ncpy_d,ncpz_d,nfh,nfh-1,nfh,v_dense)
  call overlap_fdcheck_r(nperi,ncpx_d,ncpy_d,ncpz_d,nfh,nfh-1,nfh,v_dense)
!$omp end single
  call scf_hartee_laplacian(ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,v_dense,av_dense,acfdh,nfh,ndisp)
!$omp barrier

!$omp single
  cgtmp1=0.0d0
  cgtmp2=0.0d0
!$omp end single
  call scf_hartee_cg_02(ncpx_d*ncpy_d*ncpz_d,rvec_dense,pvec_dense,av_dense,cgtmp1,cgtmp2)
!$omp barrier
!$omp single
  r1r=cgtmp1
  pap=cgtmp2

  call mpi_allreduce(pap,papall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  call mpi_allreduce(r1r,r1rall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  pap=papall
  r1r=r1rall

  alph=r1r/pap
!$omp end single
  call scf_hartee_cg_03(ncpx_d*ncpy_d*ncpz_d,vh_dense,pvec_dense,av_dense,avhs_dense,alph)
!$omp barrier

!!$omp single
!  call overlap_finitedifference_r(nperi,ncpx_d,ncpy_d,ncpz_d,nfh,nfh-1,nfh,v_dense)
!  call overlap_fdcheck_r(nperi,ncpx_d,ncpy_d,ncpz_d,nfh,nfh-1,nfh,v_dense)
!!$omp end single
!  call scf_hartee_laplacian(ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,v_dense,av_dense,acfdh,nfh,ndisp)
!!$omp barrier

!$omp single
  cgtmp1=0.0d0
!$omp end single
  call scf_hartee_cg_04(ncpx_d*ncpy_d*ncpz_d,rvec_dense,bvec_dense,avhs_dense,cgtmp1)
!$omp barrier
!$omp single
  r2r=cgtmp1

  call mpi_allreduce(r2r,r2rall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  r2r=r2rall

  beta=r2r/r1r
  if (mod(loopcg,ncgres) .eq. 0) beta=0.0d0

!$omp end single
  call scf_hartee_cg_05(ncpx_d*ncpy_d*ncpz_d,pvec_dense,rvec_dense,beta)
!$omp end parallel

!  if (myrank_glbl .eq. 0) then
!  write ( 6,*,err=9999) 'cgloop=',loopcg,'   ','cgdp=',r2r
!  write ( 6,*,err=9999) r1r,beta,pap
!  end if

  dp=dsqrt(r2r*ddx*ddy*ddz/(8.0d0*xmax*ymax*zmax))
  if (loopcg .lt. ncgmin) go to 1000
  if (loopcg .ge. ncgmax) go to 1010
  if (dp .gt. epsvh) go to 1000
1010 if (myrank_glbl .eq. 0) then
    write (ndisp,*,err=9999) 'convergence of Poisson equation'
    write (ndisp,*,err=9999) 'cgloop=',loopcg,'   ','cgdp=',dp
    inquire(unit=11,opened=isopen)
    if (isopen) then
      write (11   ,*,err=9999) 'convergence of Poisson equation'
      write (11   ,*,err=9999) 'cgloop=',loopcg,'   ','cgdp=',dp
    end if
  end if
  return
9999 continue
  call mpi_abort(mpi_comm_world,mpij)
  stop
end subroutine


subroutine scf_hartee_cg_01(n,r,b,av)
implicit none
integer i,n
real*8 av(n),b(n),r(n)

!$omp do
  do i=1,n
     r(i)=b(i)-av(i)
  end do
  return
  end subroutine


  subroutine scf_hartee_cg_02(n,r,p,av,cgtmp1,cgtmp2)
  implicit none
  integer i,n
  real*8 cgtmp1,cgtmp2,r1r,pap
  real*8 av(n),r(n),p(n)

  r1r=0.0d0
  pap=0.0d0
!$omp do
  do i=1,n
     r1r=r1r+r(i)*r(i)
     pap=pap+p(i)*av(i)
  end do
!$omp end do nowait
!$omp critical
  cgtmp1=cgtmp1+r1r
  cgtmp2=cgtmp2+pap
!$omp end critical
  return
end subroutine


subroutine scf_hartee_cg_03(n,vh,p,av,avhs,alph)
implicit none
integer i,n
real*8 alph
real*8 vh(n),av(n),p(n),avhs(n)

!$omp do
  do i=1,n
     vh(i)=vh(i)+alph*p(i)
     avhs(i)=avhs(i)+alph*av(i)
  end do
!$omp end do nowait
  return
end subroutine


subroutine scf_hartee_cg_04(n,r,b,avhs,cgtmp1)
implicit none
integer i,n
real*8 cgtmp1,r2rloc
real*8 b(n) &
      ,r(n) &
      ,avhs(n)
  r2rloc=0.0d0
!$omp do
  do i=1,n
     r(i)=b(i)-avhs(i)
     r2rloc=r2rloc+r(i)*r(i)
  end do
!$omp end do nowait
!$omp critical
  cgtmp1=cgtmp1+r2rloc
!$omp end critical
  return
end subroutine


subroutine scf_hartee_cg_05(n,p,r,beta)
implicit none
integer i,n
real*8 beta
real*8 p(n),r(n)
!$omp do
  do i=1,n
    p(i)=r(i)+beta*p(i)
  end do
!$omp end do nowait
  return
end subroutine


subroutine scf_hartee_laplacian(ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,v_dense,av_dense,acfdh,nfh,ndisp)
implicit none

integer, intent(in) :: nfh,ndisp
integer, intent(in) :: ncpx_d,ncpy_d,ncpz_d
real*8, intent(in) :: ddx,ddy,ddz
real*8, intent(in) :: acfdh(0:8)
real*8, intent(in) :: v_dense(-(nfh-1):ncpx_d+nfh,-(nfh-1):ncpy_d+nfh,-(nfh-1):ncpz_d+nfh)
real*8, intent(out) :: av_dense(ncpx_d,ncpy_d,ncpz_d)

integer ix,iy,iz
real*8 hadx,hady,hadz
real*8 xx1,xx2,xx3,xx4,xx5,xx6
real*8 yy1,yy2,yy3,yy4,yy5,yy6
real*8 zz1,zz2,zz3,zz4,zz5,zz6
real*8 diagonal

  if ((nfh .ne. 4) .and. (nfh .ne. 5) .and. (nfh .ne. 6)) then
    if (myrank_glbl .eq. 0) write(ndisp,*) 'error in Poisson_laplacian! This subroutine is for nfh ne 4, 5 or 6.'
    call mpi_abort(mpi_comm_world,mpij)
    stop
  end if

  hadx=acfdh(0)/(ddx*ddx)*(-0.5d0)
  hady=acfdh(0)/(ddy*ddy)*(-0.5d0)
  hadz=acfdh(0)/(ddz*ddz)*(-0.5d0)
  xx1=acfdh(1)/ddx/ddx*(-0.5d0)
  yy1=acfdh(1)/ddy/ddy*(-0.5d0)
  zz1=acfdh(1)/ddz/ddz*(-0.5d0)
  xx2=acfdh(2)/ddx/ddx*(-0.5d0)
  yy2=acfdh(2)/ddy/ddy*(-0.5d0)
  zz2=acfdh(2)/ddz/ddz*(-0.5d0)
  xx3=acfdh(3)/ddx/ddx*(-0.5d0)
  yy3=acfdh(3)/ddy/ddy*(-0.5d0)
  zz3=acfdh(3)/ddz/ddz*(-0.5d0)
  xx4=acfdh(4)/ddx/ddx*(-0.5d0)
  yy4=acfdh(4)/ddy/ddy*(-0.5d0)
  zz4=acfdh(4)/ddz/ddz*(-0.5d0)
  xx5=acfdh(5)/ddx/ddx*(-0.5d0)
  yy5=acfdh(5)/ddy/ddy*(-0.5d0)
  zz5=acfdh(5)/ddz/ddz*(-0.5d0)
  xx6=acfdh(6)/ddx/ddx*(-0.5d0)
  yy6=acfdh(6)/ddy/ddy*(-0.5d0)
  zz6=acfdh(6)/ddz/ddz*(-0.5d0)

  diagonal=hadx+hady+hadz

! ==========  the second derivative for nfh=4  ==========
  if (nfh .eq. 4) then
!$omp do
  do iz=1,ncpz_d
  do iy=1,ncpy_d
  do ix=1,ncpx_d
  av_dense(ix,iy,iz)= &
        +xx4*v_dense(ix-4,iy,iz) &
        +xx3*v_dense(ix-3,iy,iz) &
        +xx2*v_dense(ix-2,iy,iz) &
        +xx1*v_dense(ix-1,iy,iz) &
        +diagonal*v_dense(ix,iy,iz) &
        +xx1*v_dense(ix+1,iy,iz) &
        +xx2*v_dense(ix+2,iy,iz) &
        +xx3*v_dense(ix+3,iy,iz) &
        +xx4*v_dense(ix+4,iy,iz)
  end do
  end do
  end do
!$omp end do nowait
!$omp do
!ocl norecurrence(av_dense)
  do iz=1,ncpz_d
  do iy=1,ncpy_d
  do ix=1,ncpx_d
  av_dense(ix,iy,iz)=av_dense(ix,iy,iz) &
        +yy4*v_dense(ix,iy-4,iz) &
        +yy3*v_dense(ix,iy-3,iz) &
        +yy2*v_dense(ix,iy-2,iz) &
        +yy1*v_dense(ix,iy-1,iz) &
        +yy1*v_dense(ix,iy+1,iz) &
        +yy2*v_dense(ix,iy+2,iz) &
        +yy3*v_dense(ix,iy+3,iz) &
        +yy4*v_dense(ix,iy+4,iz)
  end do
  end do
  end do
!$omp end do nowait
!$omp do
!ocl norecurrence(av_dense)
  do iz=1,ncpz_d
  do iy=1,ncpy_d
  do ix=1,ncpx_d
  av_dense(ix,iy,iz)=av_dense(ix,iy,iz) &
        +zz4*v_dense(ix,iy,iz-4) &
        +zz3*v_dense(ix,iy,iz-3) &
        +zz2*v_dense(ix,iy,iz-2) &
        +zz1*v_dense(ix,iy,iz-1) &
        +zz1*v_dense(ix,iy,iz+1) &
        +zz2*v_dense(ix,iy,iz+2) &
        +zz3*v_dense(ix,iy,iz+3) &
        +zz4*v_dense(ix,iy,iz+4)
  end do
  end do
  end do
!$omp end do nowait
  end if
! =======================================================

! ==========  the second derivative for nfh=5  ==========
  if (nfh .eq. 5) then
!$omp do
  do iz=1,ncpz_d
  do iy=1,ncpy_d
  do ix=1,ncpx_d
  av_dense(ix,iy,iz)= &
        +xx5*v_dense(ix-5,iy,iz) &
        +xx4*v_dense(ix-4,iy,iz) &
        +xx3*v_dense(ix-3,iy,iz) &
        +xx2*v_dense(ix-2,iy,iz) &
        +xx1*v_dense(ix-1,iy,iz) &
        +diagonal*v_dense(ix,iy,iz) &
        +xx1*v_dense(ix+1,iy,iz) &
        +xx2*v_dense(ix+2,iy,iz) &
        +xx3*v_dense(ix+3,iy,iz) &
        +xx4*v_dense(ix+4,iy,iz) &
        +xx5*v_dense(ix+5,iy,iz)
  end do
  end do
  end do
!$omp end do nowait
!$omp do
!ocl norecurrence(av_dense)
  do iz=1,ncpz_d
  do iy=1,ncpy_d
  do ix=1,ncpx_d
  av_dense(ix,iy,iz)=av_dense(ix,iy,iz) &
        +yy5*v_dense(ix,iy-5,iz) &
        +yy4*v_dense(ix,iy-4,iz) &
        +yy3*v_dense(ix,iy-3,iz) &
        +yy2*v_dense(ix,iy-2,iz) &
        +yy1*v_dense(ix,iy-1,iz) &
        +yy1*v_dense(ix,iy+1,iz) &
        +yy2*v_dense(ix,iy+2,iz) &
        +yy3*v_dense(ix,iy+3,iz) &
        +yy4*v_dense(ix,iy+4,iz) &
        +yy5*v_dense(ix,iy+5,iz)
  end do
  end do
  end do
!$omp end do nowait
!$omp do
!ocl norecurrence(av_dense)
  do iz=1,ncpz_d
  do iy=1,ncpy_d
  do ix=1,ncpx_d
  av_dense(ix,iy,iz)=av_dense(ix,iy,iz) &
        +zz5*v_dense(ix,iy,iz-5) &
        +zz4*v_dense(ix,iy,iz-4) &
        +zz3*v_dense(ix,iy,iz-3) &
        +zz2*v_dense(ix,iy,iz-2) &
        +zz1*v_dense(ix,iy,iz-1) &
        +zz1*v_dense(ix,iy,iz+1) &
        +zz2*v_dense(ix,iy,iz+2) &
        +zz3*v_dense(ix,iy,iz+3) &
        +zz4*v_dense(ix,iy,iz+4) &
        +zz5*v_dense(ix,iy,iz+5)
  end do
  end do
  end do
!$omp end do nowait
  end if
! =======================================================

! ==========  the second derivative for nfh=6  ==========
  if (nfh .eq. 6) then
!$omp do
  do iz=1,ncpz_d
  do iy=1,ncpy_d
  do ix=1,ncpx_d
  av_dense(ix,iy,iz)= &
        +xx6*v_dense(ix-6,iy,iz) &
        +xx5*v_dense(ix-5,iy,iz) &
        +xx4*v_dense(ix-4,iy,iz) &
        +xx3*v_dense(ix-3,iy,iz) &
        +xx2*v_dense(ix-2,iy,iz) &
        +xx1*v_dense(ix-1,iy,iz) &
        +diagonal*v_dense(ix,iy,iz) &
        +xx1*v_dense(ix+1,iy,iz) &
        +xx2*v_dense(ix+2,iy,iz) &
        +xx3*v_dense(ix+3,iy,iz) &
        +xx4*v_dense(ix+4,iy,iz) &
        +xx5*v_dense(ix+5,iy,iz) &
        +xx6*v_dense(ix+6,iy,iz)
  end do
  end do
  end do
!$omp end do nowait
!$omp do
!ocl norecurrence(av_dense)
  do iz=1,ncpz_d
  do iy=1,ncpy_d
  do ix=1,ncpx_d
  av_dense(ix,iy,iz)=av_dense(ix,iy,iz) &
        +yy6*v_dense(ix,iy-6,iz) &
        +yy5*v_dense(ix,iy-5,iz) &
        +yy4*v_dense(ix,iy-4,iz) &
        +yy3*v_dense(ix,iy-3,iz) &
        +yy2*v_dense(ix,iy-2,iz) &
        +yy1*v_dense(ix,iy-1,iz) &
        +yy1*v_dense(ix,iy+1,iz) &
        +yy2*v_dense(ix,iy+2,iz) &
        +yy3*v_dense(ix,iy+3,iz) &
        +yy4*v_dense(ix,iy+4,iz) &
        +yy5*v_dense(ix,iy+5,iz) &
        +yy6*v_dense(ix,iy+6,iz)
  end do
  end do
  end do
!$omp end do nowait
!$omp do
!ocl norecurrence(av_dense)
  do iz=1,ncpz_d
  do iy=1,ncpy_d
  do ix=1,ncpx_d
  av_dense(ix,iy,iz)=av_dense(ix,iy,iz) &
        +zz6*v_dense(ix,iy,iz-6) &
        +zz5*v_dense(ix,iy,iz-5) &
        +zz4*v_dense(ix,iy,iz-4) &
        +zz3*v_dense(ix,iy,iz-3) &
        +zz2*v_dense(ix,iy,iz-2) &
        +zz1*v_dense(ix,iy,iz-1) &
        +zz1*v_dense(ix,iy,iz+1) &
        +zz2*v_dense(ix,iy,iz+2) &
        +zz3*v_dense(ix,iy,iz+3) &
        +zz4*v_dense(ix,iy,iz+4) &
        +zz5*v_dense(ix,iy,iz+5) &
        +zz6*v_dense(ix,iy,iz+6)
  end do
  end do
  end do
!$omp end do nowait
  end if
! =======================================================
  return
end subroutine


end module
