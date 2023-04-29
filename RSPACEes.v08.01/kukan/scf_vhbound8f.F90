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
! **********  scf_vhbound8f.f90 06/15/2014-01  **********

module mod_scf_vhbound
implicit none

contains

subroutine scf_vhbound( &
 natom,num_spe,nperi,indspe,nfh,ncpx_d,ncpy_d,ncpz_d, & ! <
 atmpole,vboundx,vboundy,vboundz,                     & ! <
 boundx,boundy,boundz)                                  ! >
use mod_mpi,       only: myrx,myry,myrz,nprocx,nprocy,nprocz
implicit none
integer, intent(in)::natom,num_spe,nperi,nfh
integer, intent(in)::ncpx_d,ncpy_d,ncpz_d
integer, intent(in)::indspe(natom)
real*8, intent(in)::atmpole(10,natom)
real*8, intent(in)::vboundx(-(nfh-1):nfh,ncpy_d,ncpz_d,9*(3-nperi)/3+1,(natom-1)*(3-nperi)/3+1)
real*8, intent(in)::vboundy(ncpx_d,-(nfh-1):nfh,ncpz_d,9*(3-nperi)/2+1,(natom-1)*(3-nperi)/2+1)
real*8, intent(in)::vboundz(ncpx_d,ncpy_d,-(nfh-1):nfh,9*(1-nperi/3)+1,(natom-1)*(1-nperi/3)+1)
real*8, intent(out)::boundx(-(nfh-1):nfh,ncpy_d,ncpz_d),boundy(ncpx_d,-(nfh-1):nfh,ncpz_d),boundz(ncpx_d,ncpy_d,-(nfh-1):nfh)
integer na,ix,iy,iz,jx,jy,jz
real*8 cp1,x,y,z,r,rin,r2in,r3in,r5in,derfdrr
real*8 atptmp,edpolx,edpoly,edpolz,eqpoxx,eqpoyy,eqpozz,eqpoxy,eqpoyz,eqpozx,eqpoxy2,eqpoyz2,eqpozx2
real*8 vmopm,vdpmx,vdpmy,vdpmz,vqmxx,vqmyy,vqmzz,vqmxy,vqmyz,vqmzx

!$omp single
  boundx(:,:,:)=0.0d0
  boundy(:,:,:)=0.0d0
  boundz(:,:,:)=0.0d0
!$omp end single
!$omp barrier

  do na=1,natom
  atptmp=atmpole( 1,na)
  edpolx=atmpole( 2,na)
  edpoly=atmpole( 3,na)
  edpolz=atmpole( 4,na)
  eqpoxx=atmpole( 5,na)
  eqpoyy=atmpole( 6,na)
  eqpozz=atmpole( 7,na)
  eqpoxy=atmpole( 8,na)
  eqpoyz=atmpole( 9,na)
  eqpozx=atmpole(10,na)
  eqpoxy2=2.0d0*eqpoxy
  eqpoyz2=2.0d0*eqpoyz
  eqpozx2=2.0d0*eqpozx
  if ((nperi<1) .and. (myrx==0)) then
!$omp do
    do iz=1,ncpz_d
    do iy=1,ncpy_d
    do ix=-(nfh-1),0
      boundx(ix,iy,iz)=boundx(ix,iy,iz)+atptmp*vboundx(ix,iy,iz,1,na) &
         +edpolx*vboundx(ix,iy,iz,2,na)+edpoly*vboundx(ix,iy,iz,3,na)+edpolz*vboundx(ix,iy,iz,4,na) &
         +eqpoxx*vboundx(ix,iy,iz,5,na)+eqpoyy*vboundx(ix,iy,iz,6,na)+eqpozz*vboundx(ix,iy,iz,7,na) &
         +eqpoxy2*vboundx(ix,iy,iz,8,na)+eqpoyz2*vboundx(ix,iy,iz,9,na)+eqpozx2*vboundx(ix,iy,iz,10,na)
    end do
    end do
    end do
  end if
!$omp barrier
  if ((nperi<1) .and. (myrx==nprocx-1)) then
!$omp do
    do iz=1,ncpz_d
    do iy=1,ncpy_d
    do ix=1,nfh
      boundx(ix,iy,iz)=boundx(ix,iy,iz)+atptmp*vboundx(ix,iy,iz,1,na) &
         +edpolx*vboundx(ix,iy,iz,2,na)+edpoly*vboundx(ix,iy,iz,3,na)+edpolz*vboundx(ix,iy,iz,4,na) &
         +eqpoxx*vboundx(ix,iy,iz,5,na)+eqpoyy*vboundx(ix,iy,iz,6,na)+eqpozz*vboundx(ix,iy,iz,7,na) &
         +eqpoxy2*vboundx(ix,iy,iz,8,na)+eqpoyz2*vboundx(ix,iy,iz,9,na)+eqpozx2*vboundx(ix,iy,iz,10,na)
    end do
    end do
    end do
  end if
!$omp barrier
  if ((nperi<2) .and. (myry==0)) then
!$omp do
    do iz=1,ncpz_d
    do iy=-(nfh-1),0
    do ix=1,ncpx_d
      boundy(ix,iy,iz)=boundy(ix,iy,iz)+atptmp*vboundy(ix,iy,iz,1,na) &
         +edpolx*vboundy(ix,iy,iz,2,na)+edpoly*vboundy(ix,iy,iz,3,na)+edpolz*vboundy(ix,iy,iz,4,na) &
         +eqpoxx*vboundy(ix,iy,iz,5,na)+eqpoyy*vboundy(ix,iy,iz,6,na)+eqpozz*vboundy(ix,iy,iz,7,na) &
         +eqpoxy2*vboundy(ix,iy,iz,8,na)+eqpoyz2*vboundy(ix,iy,iz,9,na)+eqpozx2*vboundy(ix,iy,iz,10,na)
    end do
    end do
    end do
  end if
!$omp barrier
  if ((nperi<2) .and. (myry==nprocy-1)) then
!$omp do
    do iz=1,ncpz_d
    do iy=1,nfh
    do ix=1,ncpx_d
      boundy(ix,iy,iz)=boundy(ix,iy,iz)+atptmp*vboundy(ix,iy,iz,1,na) &
         +edpolx*vboundy(ix,iy,iz,2,na)+edpoly*vboundy(ix,iy,iz,3,na)+edpolz*vboundy(ix,iy,iz,4,na) &
         +eqpoxx*vboundy(ix,iy,iz,5,na)+eqpoyy*vboundy(ix,iy,iz,6,na)+eqpozz*vboundy(ix,iy,iz,7,na) &
         +eqpoxy2*vboundy(ix,iy,iz,8,na)+eqpoyz2*vboundy(ix,iy,iz,9,na)+eqpozx2*vboundy(ix,iy,iz,10,na)
    end do
    end do
    end do
  end if
!$omp barrier
  if ((nperi<3) .and. (myrz==0)) then
!$omp do
    do iz=-(nfh-1),0
    do iy=1,ncpy_d
    do ix=1,ncpx_d
      boundz(ix,iy,iz)=boundz(ix,iy,iz)+atptmp*vboundz(ix,iy,iz,1,na) &
         +edpolx*vboundz(ix,iy,iz,2,na)+edpoly*vboundz(ix,iy,iz,3,na)+edpolz*vboundz(ix,iy,iz,4,na) &
         +eqpoxx*vboundz(ix,iy,iz,5,na)+eqpoyy*vboundz(ix,iy,iz,6,na)+eqpozz*vboundz(ix,iy,iz,7,na) &
         +eqpoxy2*vboundz(ix,iy,iz,8,na)+eqpoyz2*vboundz(ix,iy,iz,9,na)+eqpozx2*vboundz(ix,iy,iz,10,na)
    end do
    end do
    end do
  end if
!$omp barrier
  if ((nperi<3) .and. (myrz==nprocz-1)) then
!$omp do
    do iz=1,nfh
    do iy=1,ncpy_d
    do ix=1,ncpx_d
      boundz(ix,iy,iz)=boundz(ix,iy,iz)+atptmp*vboundz(ix,iy,iz,1,na) &
         +edpolx*vboundz(ix,iy,iz,2,na)+edpoly*vboundz(ix,iy,iz,3,na)+edpolz*vboundz(ix,iy,iz,4,na) &
         +eqpoxx*vboundz(ix,iy,iz,5,na)+eqpoyy*vboundz(ix,iy,iz,6,na)+eqpozz*vboundz(ix,iy,iz,7,na) &
         +eqpoxy2*vboundz(ix,iy,iz,8,na)+eqpoyz2*vboundz(ix,iy,iz,9,na)+eqpozx2*vboundz(ix,iy,iz,10,na)
    end do
    end do
    end do
  end if
  end do

return
end subroutine


end module
