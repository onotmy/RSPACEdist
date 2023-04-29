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
! **********  ggaxyz8f.F90 06/18/2014-01  **********

module mod_vxpot_ggaxyz
implicit none
contains

subroutine vxpot_ggaxyz(nperi,nfh,ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,rho,abgr,ggabg,ggr)
use mod_overlap_finitedifference, only:overlap_finitedifference_init,&
                                       overlap_finitedifference_r,overlap_fdcheck_r,&
                                       overlap_finitedifference_final
integer,intent(in) :: nperi,nfh,ncpx_d,ncpy_d,ncpz_d
real*8 ,intent(in) :: ddx,ddy,ddz
real*8 ,intent(in) :: rho(ncpx_d,ncpy_d,ncpz_d)
real*8 ,intent(out):: abgr(ncpx_d,ncpy_d,ncpz_d),ggabg(ncpx_d,ncpy_d,ncpz_d),ggr(ncpx_d,ncpy_d,ncpz_d)

! .. Local Arrays ..
real*8,allocatable :: ggrho(:,:,:),ggatmp(:,:,:),drhox(:,:,:),drhoy(:,:,:),drhoz(:,:,:)

! .. Local Scholars ..
integer :: ix,iy,iz,nfh1
real*8  :: xxx0,yyy0,zzz0,xxx1,yyy1,zzz1,xxx2,yyy2,zzz2,xxx3,yyy3,zzz3,xxx4,yyy4,zzz4
real*8  :: xxx5,yyy5,zzz5,xxx6,yyy6,zzz6,xxx7,yyy7,zzz7,xxx8,yyy8,zzz8
real*8  :: aa1,aa2,aa3,aa4,aa5,aa6,aa7,aa8
real*8  :: acfd0,acfd1,acfd2,acfd3,acfd4,acfd5,acfd6,acfd7,acfd8
real*8  :: udx,udy,udz,udx2,udy2,udz2
real*8  :: udxdy,udydz,udzdx,udydx,udzdy,udxdz
real*8  :: tdx,tdy,tdz,tdx2,tdy2,tdz2
real*8  :: tdxdy,tdydz,tdzdx
real*8  :: rhot,rhou,d,gxzet,gyzet,gzzet,zet,grhtgz
real*8  :: uabgr,uggabg,uggr,tabgr,tggabg,tggr
real*8  :: rs,vxu,exu

  nfh1=nfh-1
  allocate(ggrho (-nfh1:ncpx_d+nfh,-nfh1:ncpy_d+nfh,-nfh1:ncpz_d+nfh))
  allocate(ggatmp(-nfh1:ncpx_d+nfh,-nfh1:ncpy_d+nfh,-nfh1:ncpz_d+nfh))
  allocate(drhox (-nfh1:ncpx_d+nfh,-nfh1:ncpy_d+nfh,-nfh1:ncpz_d+nfh))
  allocate(drhoy (-nfh1:ncpx_d+nfh,-nfh1:ncpy_d+nfh,-nfh1:ncpz_d+nfh))
  allocate(drhoz (-nfh1:ncpx_d+nfh,-nfh1:ncpy_d+nfh,-nfh1:ncpz_d+nfh))

  call overlap_finitedifference_init(ncpx_d,ncpy_d,ncpz_d,nfh,1,0)

  ggrho(:,:,:)=0.d0
  drhox(:,:,:)=0.d0
  drhoy(:,:,:)=0.d0
  drhoz(:,:,:)=0.d0

  do iz=1,ncpz_d
  do iy=1,ncpy_d
  do ix=1,ncpx_d
    ggrho(ix,iy,iz)=rho(ix,iy,iz)
  end do
  end do
  end do

  ggatmp(:,:,:)=ggrho(:,:,:)
  call overlap_finitedifference_r(nperi,ncpx_d,ncpy_d,ncpz_d,nfh,nfh1,nfh,ggatmp)
  call overlap_fdcheck_r(nperi,ncpx_d,ncpy_d,ncpz_d,nfh,nfh1,nfh,ggatmp)
  ggrho(:,:,:)=ggatmp(:,:,:)

  select case (nfh)
  case (4)
    aa1=-672.0d0/840.0d0
    aa2=168.0d0/840.0d0
    aa3=-32.0d0/840.0d0
    aa4=3.0d0/840.0d0
    aa5=0.0d0
    aa6=0.0d0
    aa7=0.0d0
    aa8=0.0d0
    acfd0=-205.0d0/72.0d0
    acfd1=8.0d0/5.0d0
    acfd2=-1.0d0/5.0d0
    acfd3=8.0d0/315.0d0
    acfd4=-1.0d0/560.0d0
    acfd5=0.0d0
    acfd6=0.0d0
    acfd7=0.0d0
    acfd8=0.0d0
  case (5)
    aa1=77.0d0/24.0d0
    aa2=-107.0d0/12.0d0
    aa3=39.0d0/4.0d0
    aa4=-61.0d0/12.0d0
    aa5=25.0d0/24.0d0
    aa6=0.0d0
    aa7=0.0d0
    aa8=0.0d0
    acfd0=-36883.0d0/12600.0d0
    acfd1=5.0d0/3.0d0
    acfd2=-5.0d0/21.0d0
    acfd3=5.0d0/126.0d0
    acfd4=-5.0d0/1008.0d0
    acfd5=1.0d0/3150.0d0
    acfd6=0.0d0
    acfd7=0.0d0
    acfd8=0.0d0
  case(6)
    aa1= 4.35d0
    aa2=-14.625d0
    aa3= 21.166666666666667d0
    aa4=-16.5d0
    aa5= 6.75d0
    aa6=-1.1416666666666667d0
    aa7=0.0d0
    aa8=0.0d0
    acfd0=-5369.0d0/1800.0d0
    acfd1=12.0d0/7.0d0
    acfd2=-15.0d0/56.0d0
    acfd3=10.0d0/189.0d0
    acfd4=-1.0d0/112.0d0
    acfd5=2.0d0/1925.0d0
    acfd6=-1.0d0/16632.0d0
    acfd7=0.0d0
    acfd8=0.0d0
  end select
  xxx0=acfd0/ddx/ddx
  yyy0=acfd0/ddy/ddy
  zzz0=acfd0/ddz/ddz
  xxx1=acfd1/ddx/ddx
  yyy1=acfd1/ddy/ddy
  zzz1=acfd1/ddz/ddz
  xxx2=acfd2/ddx/ddx
  yyy2=acfd2/ddy/ddy
  zzz2=acfd2/ddz/ddz
  xxx3=acfd3/ddx/ddx
  yyy3=acfd3/ddy/ddy
  zzz3=acfd3/ddz/ddz
  xxx4=acfd4/ddx/ddx
  yyy4=acfd4/ddy/ddy
  zzz4=acfd4/ddz/ddz
  xxx5=acfd5/ddx/ddx
  yyy5=acfd5/ddy/ddy
  zzz5=acfd5/ddz/ddz
  xxx6=acfd6/ddx/ddx
  yyy6=acfd6/ddy/ddy
  zzz6=acfd6/ddz/ddz
  xxx7=acfd7/ddx/ddx
  yyy7=acfd7/ddy/ddy
  zzz7=acfd7/ddz/ddz
  xxx8=acfd8/ddx/ddx
  yyy8=acfd8/ddy/ddy
  zzz8=acfd8/ddz/ddz

  select case (nfh)
  case (4)
!$omp parallel default(shared),private(ix,iy,iz)
!$omp do
    do iz=1,ncpz_d
    do iy=1,ncpy_d
    do ix=1,ncpx_d
       drhox(ix,iy,iz) = (aa4*ggrho(ix-4,iy,iz)+aa3*ggrho(ix-3,iy,iz) &
                         +aa2*ggrho(ix-2,iy,iz)+aa1*ggrho(ix-1,iy,iz) &
                         -aa1*ggrho(ix+1,iy,iz)-aa2*ggrho(ix+2,iy,iz) &
                         -aa3*ggrho(ix+3,iy,iz)-aa4*ggrho(ix+4,iy,iz))/ddx
       drhoy(ix,iy,iz) = (aa4*ggrho(ix,iy-4,iz)+aa3*ggrho(ix,iy-3,iz) &
                         +aa2*ggrho(ix,iy-2,iz)+aa1*ggrho(ix,iy-1,iz) &
                         -aa1*ggrho(ix,iy+1,iz)-aa2*ggrho(ix,iy+2,iz) &
                         -aa3*ggrho(ix,iy+3,iz)-aa4*ggrho(ix,iy+4,iz))/ddy
       drhoz(ix,iy,iz) = (aa4*ggrho(ix,iy,iz-4)+aa3*ggrho(ix,iy,iz-3) &
                         +aa2*ggrho(ix,iy,iz-2)+aa1*ggrho(ix,iy,iz-1) &
                         -aa1*ggrho(ix,iy,iz+1)-aa2*ggrho(ix,iy,iz+2) &
                         -aa3*ggrho(ix,iy,iz+3)-aa4*ggrho(ix,iy,iz+4))/ddz
    end do
    end do
    end do
!$omp end parallel
  case (5)
!$omp parallel default(shared),private(ix,iy,iz)
!$omp do
    do iz=1,ncpz_d
    do iy=1,ncpy_d
    do ix=1,ncpx_d
       drhox(ix,iy,iz) = (aa5*ggrho(ix-5,iy,iz)+aa4*ggrho(ix-4,iy,iz) &
                         +aa3*ggrho(ix-3,iy,iz)+aa2*ggrho(ix-2,iy,iz) &
                         +aa1*ggrho(ix-1,iy,iz)-aa1*ggrho(ix+1,iy,iz) &
                         -aa2*ggrho(ix+2,iy,iz)-aa3*ggrho(ix+3,iy,iz) &
                         -aa4*ggrho(ix+4,iy,iz)-aa5*ggrho(ix+5,iy,iz))/ddx
       drhoy(ix,iy,iz) = (aa5*ggrho(ix,iy-5,iz)+aa4*ggrho(ix,iy-4,iz) &
                         +aa3*ggrho(ix,iy-3,iz)+aa2*ggrho(ix,iy-2,iz) &
                         +aa1*ggrho(ix,iy-1,iz)-aa1*ggrho(ix,iy+1,iz) &
                         -aa2*ggrho(ix,iy+2,iz)-aa3*ggrho(ix,iy+3,iz) &
                         -aa4*ggrho(ix,iy+4,iz)-aa5*ggrho(ix,iy+5,iz))/ddy
       drhoz(ix,iy,iz) = (aa5*ggrho(ix,iy,iz-5)+aa4*ggrho(ix,iy,iz-4) &
                         +aa3*ggrho(ix,iy,iz-3)+aa2*ggrho(ix,iy,iz-2) &
                         +aa1*ggrho(ix,iy,iz-1)-aa1*ggrho(ix,iy,iz+1) &
                         -aa2*ggrho(ix,iy,iz+2)-aa3*ggrho(ix,iy,iz+3) &
                         -aa4*ggrho(ix,iy,iz+4)-aa5*ggrho(ix,iy,iz+5))/ddz
    end do
    end do
    end do
!$omp end parallel
  case (6)
!$omp parallel default(shared),private(ix,iy,iz)
!$omp do
    do iz=1,ncpz_d
    do iy=1,ncpy_d
    do ix=1,ncpx_d
       drhox(ix,iy,iz) = (aa6*ggrho(ix-6,iy,iz)+aa5*ggrho(ix-5,iy,iz) &
                         +aa4*ggrho(ix-4,iy,iz)+aa3*ggrho(ix-3,iy,iz) &
                         +aa2*ggrho(ix-2,iy,iz)+aa1*ggrho(ix-1,iy,iz) &
                         -aa1*ggrho(ix+1,iy,iz)-aa2*ggrho(ix+2,iy,iz) &
                         -aa3*ggrho(ix+3,iy,iz)-aa4*ggrho(ix+4,iy,iz) &
                         -aa5*ggrho(ix+5,iy,iz)-aa6*ggrho(ix+6,iy,iz))/ddx
       drhoy(ix,iy,iz) = (aa6*ggrho(ix,iy-6,iz)+aa5*ggrho(ix,iy-5,iz) &
                         +aa4*ggrho(ix,iy-4,iz)+aa3*ggrho(ix,iy-3,iz) &
                         +aa2*ggrho(ix,iy-2,iz)+aa1*ggrho(ix,iy-1,iz) &
                         -aa1*ggrho(ix,iy+1,iz)-aa2*ggrho(ix,iy+2,iz) &
                         -aa3*ggrho(ix,iy+3,iz)-aa4*ggrho(ix,iy+4,iz) &
                         -aa5*ggrho(ix,iy+5,iz)-aa6*ggrho(ix,iy+6,iz))/ddy
       drhoz(ix,iy,iz) = (aa6*ggrho(ix,iy,iz-6)+aa5*ggrho(ix,iy,iz-5) &
                         +aa4*ggrho(ix,iy,iz-4)+aa3*ggrho(ix,iy,iz-3) &
                         +aa2*ggrho(ix,iy,iz-2)+aa1*ggrho(ix,iy,iz-1) &
                         -aa1*ggrho(ix,iy,iz+1)-aa2*ggrho(ix,iy,iz+2) &
                         -aa3*ggrho(ix,iy,iz+3)-aa4*ggrho(ix,iy,iz+4) &
                         -aa5*ggrho(ix,iy,iz+5)-aa6*ggrho(ix,iy,iz+6))/ddz
    end do
    end do
    end do
!$omp end parallel
  end select

  ggatmp(:,:,:)=drhox(:,:,:)
  call overlap_finitedifference_r(nperi,ncpx_d,ncpy_d,ncpz_d,nfh,nfh1,nfh,ggatmp)
  call overlap_fdcheck_r(nperi,ncpx_d,ncpy_d,ncpz_d,nfh,nfh1,nfh,ggatmp)
  drhox(:,:,:)=ggatmp(:,:,:)

  ggatmp(:,:,:)=drhoy(:,:,:)
  call overlap_finitedifference_r(nperi,ncpx_d,ncpy_d,ncpz_d,nfh,nfh1,nfh,ggatmp)
  call overlap_fdcheck_r(nperi,ncpx_d,ncpy_d,ncpz_d,nfh,nfh1,nfh,ggatmp)
  drhoy(:,:,:)=ggatmp(:,:,:)

  ggatmp(:,:,:)=drhoz(:,:,:)
  call overlap_finitedifference_r(nperi,ncpx_d,ncpy_d,ncpz_d,nfh,nfh1,nfh,ggatmp)
  call overlap_fdcheck_r(nperi,ncpx_d,ncpy_d,ncpz_d,nfh,nfh1,nfh,ggatmp)
  drhoz(:,:,:)=ggatmp(:,:,:)

  select case (nfh)
  case (4)
!$omp parallel default(firstprivate),shared(ggrho,drhox,drhoy,drhoz,abgr,ggabg,ggr)
!$omp do
    do iz=1,ncpz_d
    do iy=1,ncpy_d
    do ix=1,ncpx_d
       udx   =drhox(ix,iy,iz)
       udy   =drhoy(ix,iy,iz)
       udz   =drhoz(ix,iy,iz)
       udx2  = xxx4*ggrho(ix-4,iy,iz)+xxx3*ggrho(ix-3,iy,iz) &
              +xxx2*ggrho(ix-2,iy,iz)+xxx1*ggrho(ix-1,iy,iz) &
              +xxx0*ggrho(ix  ,iy,iz)+xxx1*ggrho(ix+1,iy,iz) &
              +xxx2*ggrho(ix+2,iy,iz)+xxx3*ggrho(ix+3,iy,iz) &
              +xxx4*ggrho(ix+4,iy,iz)
       udy2  = yyy4*ggrho(ix,iy-4,iz)+yyy3*ggrho(ix,iy-3,iz) &
              +yyy2*ggrho(ix,iy-2,iz)+yyy1*ggrho(ix,iy-1,iz) &
              +yyy0*ggrho(ix,iy  ,iz)+yyy1*ggrho(ix,iy+1,iz) &
              +yyy2*ggrho(ix,iy+2,iz)+yyy3*ggrho(ix,iy+3,iz) &
              +yyy4*ggrho(ix,iy+4,iz)
       udz2  = zzz4*ggrho(ix,iy,iz-4)+zzz3*ggrho(ix,iy,iz-3) &
              +zzz2*ggrho(ix,iy,iz-2)+zzz1*ggrho(ix,iy,iz-1) &
              +zzz0*ggrho(ix,iy,iz  )+zzz1*ggrho(ix,iy,iz+1) &
              +zzz2*ggrho(ix,iy,iz+2)+zzz3*ggrho(ix,iy,iz+3) &
              +zzz4*ggrho(ix,iy,iz+4)
       udxdy = (aa4*drhoy(ix-4,iy,iz)+aa3*drhoy(ix-3,iy,iz) &
               +aa2*drhoy(ix-2,iy,iz)+aa1*drhoy(ix-1,iy,iz) &
               -aa1*drhoy(ix+1,iy,iz)-aa2*drhoy(ix+2,iy,iz) &
               -aa3*drhoy(ix+3,iy,iz)-aa4*drhoy(ix+4,iy,iz))/ddx
       udydz = (aa4*drhoz(ix,iy-4,iz)+aa3*drhoz(ix,iy-3,iz) &
               +aa2*drhoz(ix,iy-2,iz)+aa1*drhoz(ix,iy-1,iz) &
               -aa1*drhoz(ix,iy+1,iz)-aa2*drhoz(ix,iy+2,iz) &
               -aa3*drhoz(ix,iy+3,iz)-aa4*drhoz(ix,iy+4,iz))/ddy
       udzdx = (aa4*drhox(ix,iy,iz-4)+aa3*drhox(ix,iy,iz-3) &
               +aa2*drhox(ix,iy,iz-2)+aa1*drhox(ix,iy,iz-1) &
               -aa1*drhox(ix,iy,iz+1)-aa2*drhox(ix,iy,iz+2) &
               -aa3*drhox(ix,iy,iz+3)-aa4*drhox(ix,iy,iz+4))/ddz
       udydx = (aa4*drhox(ix,iy-4,iz)+aa3*drhox(ix,iy-3,iz) &
               +aa2*drhox(ix,iy-2,iz)+aa1*drhox(ix,iy-1,iz) &
               -aa1*drhox(ix,iy+1,iz)-aa2*drhox(ix,iy+2,iz) &
               -aa3*drhox(ix,iy+3,iz)-aa4*drhox(ix,iy+4,iz))/ddy
       udzdy = (aa4*drhoy(ix,iy,iz-4)+aa3*drhoy(ix,iy,iz-3) &
               +aa2*drhoy(ix,iy,iz-2)+aa1*drhoy(ix,iy,iz-1) &
               -aa1*drhoy(ix,iy,iz+1)-aa2*drhoy(ix,iy,iz+2) &
               -aa3*drhoy(ix,iy,iz+3)-aa4*drhoy(ix,iy,iz+4))/ddz
       udxdz = (aa4*drhoz(ix-4,iy,iz)+aa3*drhoz(ix-3,iy,iz) &
               +aa2*drhoz(ix-2,iy,iz)+aa1*drhoz(ix-1,iy,iz) &
               -aa1*drhoz(ix+1,iy,iz)-aa2*drhoz(ix+2,iy,iz) &
               -aa3*drhoz(ix+3,iy,iz)-aa4*drhoz(ix+4,iy,iz))/ddx
       udxdy =0.5d0*(udxdy+udydx)
       udydz =0.5d0*(udydz+udzdy)
       udzdx =0.5d0*(udzdx+udxdz)
       uabgr =sqrt(udx**2+udy**2+udz**2)
       uggabg=(udx**2*udx2+udy**2*udy2+udz**2*udz2 &
              +2.d0*(udx*udy*udxdy+udy*udz*udydz+udz*udx*udzdx))/uabgr
       uggr  =udx2+udy2+udz2
       abgr(ix,iy,iz)=uabgr
       ggabg(ix,iy,iz)=uggabg
       ggr(ix,iy,iz)=uggr
    end do
    end do
    end do
!$omp end parallel
  case (5)
!$omp parallel default(firstprivate),shared(ggrho,drhox,drhoy,drhoz,abgr,ggabg,ggr)
!$omp do
    do iz=1,ncpz_d
    do iy=1,ncpy_d
    do ix=1,ncpx_d
       udx   =drhox(ix,iy,iz)
       udy   =drhoy(ix,iy,iz)
       udz   =drhoz(ix,iy,iz)
       udx2  = xxx5*ggrho(ix-5,iy,iz)+xxx4*ggrho(ix-4,iy,iz) &
              +xxx3*ggrho(ix-3,iy,iz)+xxx2*ggrho(ix-2,iy,iz) &
              +xxx1*ggrho(ix-1,iy,iz)+xxx0*ggrho(ix  ,iy,iz) &
              +xxx1*ggrho(ix+1,iy,iz)+xxx2*ggrho(ix+2,iy,iz) &
              +xxx3*ggrho(ix+3,iy,iz)+xxx4*ggrho(ix+4,iy,iz) &
              +xxx5*ggrho(ix+5,iy,iz)
       udy2  = yyy5*ggrho(ix,iy-5,iz)+yyy4*ggrho(ix,iy-4,iz) &
              +yyy3*ggrho(ix,iy-3,iz)+yyy2*ggrho(ix,iy-2,iz) &
              +yyy1*ggrho(ix,iy-1,iz)+yyy0*ggrho(ix,iy  ,iz) &
              +yyy1*ggrho(ix,iy+1,iz)+yyy2*ggrho(ix,iy+2,iz) &
              +yyy3*ggrho(ix,iy+3,iz)+yyy4*ggrho(ix,iy+4,iz) &
              +yyy5*ggrho(ix,iy+5,iz)
       udz2  = zzz5*ggrho(ix,iy,iz-5)+zzz4*ggrho(ix,iy,iz-4) &
              +zzz3*ggrho(ix,iy,iz-3)+zzz2*ggrho(ix,iy,iz-2) &
              +zzz1*ggrho(ix,iy,iz-1)+zzz0*ggrho(ix,iy,iz  ) &
              +zzz1*ggrho(ix,iy,iz+1)+zzz2*ggrho(ix,iy,iz+2) &
              +zzz3*ggrho(ix,iy,iz+3)+zzz4*ggrho(ix,iy,iz+4) &
              +zzz5*ggrho(ix,iy,iz+5)
       udxdy = (aa5*drhoy(ix-5,iy,iz)+aa4*drhoy(ix-4,iy,iz) &
               +aa3*drhoy(ix-3,iy,iz)+aa2*drhoy(ix-2,iy,iz) &
               +aa1*drhoy(ix-1,iy,iz)-aa1*drhoy(ix+1,iy,iz) &
               -aa2*drhoy(ix+2,iy,iz)-aa3*drhoy(ix+3,iy,iz) &
               -aa4*drhoy(ix+4,iy,iz)-aa5*drhoy(ix+5,iy,iz))/ddx
       udydz = (aa5*drhoz(ix,iy-5,iz)+aa4*drhoz(ix,iy-4,iz) &
               +aa3*drhoz(ix,iy-3,iz)+aa2*drhoz(ix,iy-2,iz) &
               +aa1*drhoz(ix,iy-1,iz)-aa1*drhoz(ix,iy+1,iz) &
               -aa2*drhoz(ix,iy+2,iz)-aa3*drhoz(ix,iy+3,iz) &
               -aa4*drhoz(ix,iy+4,iz)-aa5*drhoz(ix,iy+5,iz))/ddy
       udzdx = (aa5*drhox(ix,iy,iz-5)+aa4*drhox(ix,iy,iz-4) &
               +aa3*drhox(ix,iy,iz-3)+aa2*drhox(ix,iy,iz-2) &
               +aa1*drhox(ix,iy,iz-1)-aa1*drhox(ix,iy,iz+1) &
               -aa2*drhox(ix,iy,iz+2)-aa3*drhox(ix,iy,iz+3) &
               -aa4*drhox(ix,iy,iz+4)-aa5*drhox(ix,iy,iz+5))/ddz
       udydx = (aa5*drhox(ix,iy-5,iz)+aa4*drhox(ix,iy-4,iz) &
               +aa3*drhox(ix,iy-3,iz)+aa2*drhox(ix,iy-2,iz) &
               +aa1*drhox(ix,iy-1,iz)-aa1*drhox(ix,iy+1,iz) &
               -aa2*drhox(ix,iy+2,iz)-aa3*drhox(ix,iy+3,iz) &
               -aa4*drhox(ix,iy+4,iz)-aa5*drhox(ix,iy+5,iz))/ddy
       udzdy = (aa5*drhoy(ix,iy,iz-5)+aa4*drhoy(ix,iy,iz-4) &
               +aa3*drhoy(ix,iy,iz-3)+aa2*drhoy(ix,iy,iz-2) &
               +aa1*drhoy(ix,iy,iz-1)-aa1*drhoy(ix,iy,iz+1) &
               -aa2*drhoy(ix,iy,iz+2)-aa3*drhoy(ix,iy,iz+3) &
               -aa4*drhoy(ix,iy,iz+4)-aa5*drhoy(ix,iy,iz+5))/ddz
       udxdz = (aa5*drhoz(ix-5,iy,iz)+aa4*drhoz(ix-4,iy,iz) &
               +aa3*drhoz(ix-3,iy,iz)+aa2*drhoz(ix-2,iy,iz) &
               +aa1*drhoz(ix-1,iy,iz)-aa1*drhoz(ix+1,iy,iz) &
               -aa2*drhoz(ix+2,iy,iz)-aa3*drhoz(ix+3,iy,iz) &
               -aa4*drhoz(ix+4,iy,iz)-aa5*drhoz(ix+5,iy,iz))/ddx
       udxdy =0.5d0*(udxdy+udydx)
       udydz =0.5d0*(udydz+udzdy)
       udzdx =0.5d0*(udzdx+udxdz)
       uabgr =sqrt(udx**2+udy**2+udz**2)
       uggabg=(udx**2*udx2+udy**2*udy2+udz**2*udz2 &
              +2.d0*(udx*udy*udxdy+udy*udz*udydz+udz*udx*udzdx))/uabgr
       uggr  =udx2+udy2+udz2
       abgr(ix,iy,iz)=uabgr
       ggabg(ix,iy,iz)=uggabg
       ggr(ix,iy,iz)=uggr
    end do
    end do
    end do
!$omp end parallel
  case (6)
!$omp parallel default(firstprivate),shared(ggrho,drhox,drhoy,drhoz,abgr,ggabg,ggr)
!$omp do
    do iz=1,ncpz_d
    do iy=1,ncpy_d
    do ix=1,ncpx_d
       udx   =drhox(ix,iy,iz)
       udy   =drhoy(ix,iy,iz)
       udz   =drhoz(ix,iy,iz)
       udx2  = xxx6*ggrho(ix-6,iy,iz)+xxx5*ggrho(ix-5,iy,iz) &
              +xxx4*ggrho(ix-4,iy,iz)+xxx3*ggrho(ix-3,iy,iz) &
              +xxx2*ggrho(ix-2,iy,iz)+xxx1*ggrho(ix-1,iy,iz) &
              +xxx0*ggrho(ix  ,iy,iz)+xxx1*ggrho(ix+1,iy,iz) &
              +xxx2*ggrho(ix+2,iy,iz)+xxx3*ggrho(ix+3,iy,iz) &
              +xxx4*ggrho(ix+4,iy,iz)+xxx5*ggrho(ix+5,iy,iz) &
              +xxx6*ggrho(ix+6,iy,iz)
       udy2  = yyy6*ggrho(ix,iy-6,iz)+yyy5*ggrho(ix,iy-5,iz) &
              +yyy4*ggrho(ix,iy-4,iz)+yyy3*ggrho(ix,iy-3,iz) &
              +yyy2*ggrho(ix,iy-2,iz)+yyy1*ggrho(ix,iy-1,iz) &
              +yyy0*ggrho(ix,iy  ,iz)+yyy1*ggrho(ix,iy+1,iz) &
              +yyy2*ggrho(ix,iy+2,iz)+yyy3*ggrho(ix,iy+3,iz) &
              +yyy4*ggrho(ix,iy+4,iz)+yyy5*ggrho(ix,iy+5,iz) &
              +yyy6*ggrho(ix,iy+6,iz)
       udz2  = zzz6*ggrho(ix,iy,iz-6)+zzz5*ggrho(ix,iy,iz-5) &
              +zzz4*ggrho(ix,iy,iz-4)+zzz3*ggrho(ix,iy,iz-3) &
              +zzz2*ggrho(ix,iy,iz-2)+zzz1*ggrho(ix,iy,iz-1) &
              +zzz0*ggrho(ix,iy,iz  )+zzz1*ggrho(ix,iy,iz+1) &
              +zzz2*ggrho(ix,iy,iz+2)+zzz3*ggrho(ix,iy,iz+3) &
              +zzz4*ggrho(ix,iy,iz+4)+zzz5*ggrho(ix,iy,iz+5) &
              +zzz6*ggrho(ix,iy,iz+6)
       udxdy = (aa6*drhoy(ix-6,iy,iz)+aa5*drhoy(ix-5,iy,iz) &
               +aa4*drhoy(ix-4,iy,iz)+aa3*drhoy(ix-3,iy,iz) &
               +aa2*drhoy(ix-2,iy,iz)+aa1*drhoy(ix-1,iy,iz) &
               -aa1*drhoy(ix+1,iy,iz)-aa2*drhoy(ix+2,iy,iz) &
               -aa3*drhoy(ix+3,iy,iz)-aa4*drhoy(ix+4,iy,iz) &
               -aa5*drhoy(ix+5,iy,iz)-aa6*drhoy(ix+6,iy,iz))/ddx
       udydz = (aa6*drhoz(ix,iy-6,iz)+aa5*drhoz(ix,iy-5,iz) &
               +aa4*drhoz(ix,iy-4,iz)+aa3*drhoz(ix,iy-3,iz) &
               +aa2*drhoz(ix,iy-2,iz)+aa1*drhoz(ix,iy-1,iz) &
               -aa1*drhoz(ix,iy+1,iz)-aa2*drhoz(ix,iy+2,iz) &
               -aa3*drhoz(ix,iy+3,iz)-aa4*drhoz(ix,iy+4,iz) &
               -aa5*drhoz(ix,iy+5,iz)-aa6*drhoz(ix,iy+6,iz))/ddy
       udzdx = (aa6*drhox(ix,iy,iz-6)+aa5*drhox(ix,iy,iz-5) &
               +aa4*drhox(ix,iy,iz-4)+aa3*drhox(ix,iy,iz-3) &
               +aa2*drhox(ix,iy,iz-2)+aa1*drhox(ix,iy,iz-1) &
               -aa1*drhox(ix,iy,iz+1)-aa2*drhox(ix,iy,iz+2) &
               -aa3*drhox(ix,iy,iz+3)-aa4*drhox(ix,iy,iz+4) &
               -aa5*drhox(ix,iy,iz+5)-aa6*drhox(ix,iy,iz+6))/ddz
       udydx = (aa6*drhox(ix,iy-6,iz)+aa5*drhox(ix,iy-5,iz) &
               +aa4*drhox(ix,iy-4,iz)+aa3*drhox(ix,iy-3,iz) &
               +aa2*drhox(ix,iy-2,iz)+aa1*drhox(ix,iy-1,iz) &
               -aa1*drhox(ix,iy+1,iz)-aa2*drhox(ix,iy+2,iz) &
               -aa3*drhox(ix,iy+3,iz)-aa4*drhox(ix,iy+4,iz) &
               -aa5*drhox(ix,iy+5,iz)-aa6*drhox(ix,iy+6,iz))/ddy
       udzdy = (aa6*drhoy(ix,iy,iz-6)+aa5*drhoy(ix,iy,iz-5) &
               +aa4*drhoy(ix,iy,iz-4)+aa3*drhoy(ix,iy,iz-3) &
               +aa2*drhoy(ix,iy,iz-2)+aa1*drhoy(ix,iy,iz-1) &
               -aa1*drhoy(ix,iy,iz+1)-aa2*drhoy(ix,iy,iz+2) &
               -aa3*drhoy(ix,iy,iz+3)-aa4*drhoy(ix,iy,iz+4) &
               -aa5*drhoy(ix,iy,iz+5)-aa6*drhoy(ix,iy,iz+6))/ddz
       udxdz = (aa6*drhoz(ix-6,iy,iz)+aa5*drhoz(ix-5,iy,iz) &
               +aa4*drhoz(ix-4,iy,iz)+aa3*drhoz(ix-3,iy,iz) &
               +aa2*drhoz(ix-2,iy,iz)+aa1*drhoz(ix-1,iy,iz) &
               -aa1*drhoz(ix+1,iy,iz)-aa2*drhoz(ix+2,iy,iz) &
               -aa3*drhoz(ix+3,iy,iz)-aa4*drhoz(ix+4,iy,iz) &
               -aa5*drhoz(ix+5,iy,iz)-aa6*drhoz(ix+6,iy,iz))/ddx
       udxdy =0.5d0*(udxdy+udydx)
       udydz =0.5d0*(udydz+udzdy)
       udzdx =0.5d0*(udzdx+udxdz)
       uabgr =sqrt(udx**2+udy**2+udz**2)
       uggabg=(udx**2*udx2+udy**2*udy2+udz**2*udz2 &
              +2.d0*(udx*udy*udxdy+udy*udz*udydz+udz*udx*udzdx))/uabgr
       uggr  =udx2+udy2+udz2
       abgr(ix,iy,iz)=uabgr
       ggabg(ix,iy,iz)=uggabg
       ggr(ix,iy,iz)=uggr
    end do
    end do
    end do
!$omp end parallel
  end select

  deallocate(ggrho,ggatmp,drhox,drhoy,drhoz)
  call overlap_finitedifference_final

return
end subroutine vxpot_ggaxyz


end module
