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
! **********  scf_fuzzycellmoment8e.f90 01/06/2020-01  **********

module mod_scf_fuzzycellmoment
implicit none
contains


subroutine scf_fuzzycellmoment( &
 key_natpri_inps,key_pp_paw,nperi,nspv,natom,num_spe,num_ppcell_d,        & ! <
 num_list_d,nmesh,ncpx,ncpy,ncpz,                                         & ! <
 ntyppp,indspe,natprid,natinfd,napsd,ndatx,ndaty,ndatz,lstdx,lstdy,lstdz, & ! <
 dx,dy,dz,xmax,ymax,zmax,                                                 & ! <
 atx,aty,atz,rhosmt,rhoaug3d,pwei,                                        & ! <
 atmpole)                                                                   ! >
use mod_mpi
implicit none
integer,intent(in) ::key_natpri_inps,key_pp_paw
integer,intent(in) ::nperi,nspv,natom,num_spe,num_ppcell_d,num_list_d,nmesh,ncpx,ncpy,ncpz
integer,intent(in) ::natprid(natom),natinfd(natom),napsd(natom),ndatx(natom),ndaty(natom),ndatz(natom)
integer,intent(in) ::lstdx(num_list_d,num_ppcell_d),lstdy(num_list_d,num_ppcell_d),lstdz(num_list_d,num_ppcell_d)
integer,intent(in) ::ntyppp(num_spe),indspe(natom)
real*8, intent(in) ::dx,dy,dz,xmax,ymax,zmax
real*8, intent(in) ::atx(natom),aty(natom),atz(natom)
real*8, intent(in) ::rhoaug3d(num_list_d,num_ppcell_d)
real*8, intent(in) ::rhosmt(ncpx,ncpy,ncpz,nspv)
real*8, intent(in) ::pwei(ncpx,ncpy,ncpz,natom)
real*8, intent(out)::atmpole(10,natom)
integer l
integer na,ixyz,ix,iy,iz,jx,jy,jz,is,ns
real*8  ddx,ddy,ddz,dxyz,ddxyz,atomx,atomy,atomz,tptmp,dpolx,dpoly,dpolz,qpoxx,qpoyy,qpozz,qpoxy,qpoyz,qpozx,x,y,z,rhowei

!$omp do
  do ix=1,10*natom
    atmpole(ix,1)=0.0d0
  end do

  ddx=dx/nmesh
  ddy=dy/nmesh
  ddz=dz/nmesh
  dxyz=dx*dy*dz
  ddxyz=ddx*ddy*ddz

  ns=1
  if (nspv>1) ns=2

  do na=1,natom
  atomx=atx(na)
  atomy=aty(na)
  atomz=atz(na)
  tptmp=0.0d0
  dpolx=0.0d0
  dpoly=0.0d0
  dpolz=0.0d0
  qpoxx=0.0d0
  qpoyy=0.0d0
  qpozz=0.0d0
  qpoxy=0.0d0
  qpoyz=0.0d0
  qpozx=0.0d0
  if (nperi<3) then
    if ((ntyppp(indspe(na))==key_pp_paw) .and. (natprid(na)==key_natpri_inps)) then
!$omp do
      do ixyz=1,natinfd(na)
        ix=lstdx(ixyz,napsd(na))
        iy=lstdy(ixyz,napsd(na))
        iz=lstdz(ixyz,napsd(na))
        x=ix*ddx-(atx(na)-(ndatx(na)*ddx-0.5d0*ddx))
        y=iy*ddy-(aty(na)-(ndaty(na)*ddy-0.5d0*ddy))
        z=iz*ddz-(atz(na)-(ndatz(na)*ddz-0.5d0*ddz))
        rhowei=rhoaug3d(ixyz,napsd(na))*ddxyz
        tptmp=tptmp+rhowei
        dpolx=dpolx+x*rhowei
        dpoly=dpoly+y*rhowei
        dpolz=dpolz+z*rhowei
        qpoxx=qpoxx+0.5d0*x*x*rhowei
        qpoyy=qpoyy+0.5d0*y*y*rhowei
        qpozz=qpozz+0.5d0*z*z*rhowei
        qpoxy=qpoxy+0.5d0*x*y*rhowei
        qpoyz=qpoyz+0.5d0*y*z*rhowei
        qpozx=qpozx+0.5d0*z*x*rhowei
      end do
    end if
    do is=1,ns
      select case (nperi)
      case(0)
!$omp do
        do iz=1,ncpz
        do iy=1,ncpy
        do ix=1,ncpx
          jx=myrx*ncpx+ix
          jy=myry*ncpy+iy
          jz=myrz*ncpz+iz
          x=(jx*dx-xmax-0.5d0*dx)-atomx
          y=(jy*dy-ymax-0.5d0*dy)-atomy
          z=(jz*dz-zmax-0.5d0*dz)-atomz
          rhowei=rhosmt(ix,iy,iz,is)*pwei(ix,iy,iz,na)*dxyz
          tptmp=tptmp+rhowei
          dpolx=dpolx+x*rhowei
          dpoly=dpoly+y*rhowei
          dpolz=dpolz+z*rhowei
          qpoxx=qpoxx+0.5d0*x*x*rhowei
          qpoyy=qpoyy+0.5d0*y*y*rhowei
          qpozz=qpozz+0.5d0*z*z*rhowei
          qpoxy=qpoxy+0.5d0*x*y*rhowei
          qpoyz=qpoyz+0.5d0*y*z*rhowei
          qpozx=qpozx+0.5d0*z*x*rhowei
        end do
        end do
        end do
      case(1)
!$omp do
        do iz=1,ncpz
        do iy=1,ncpy
        do ix=1,ncpx
          jx=myrx*ncpx+ix
          jy=myry*ncpy+iy
          jz=myrz*ncpz+iz
          x=(jx*dx-xmax-0.5d0*dx)-atomx
          y=(jy*dy-ymax-0.5d0*dy)-atomy
          z=(jz*dz-zmax-0.5d0*dz)-atomz
          if (x .gt. xmax) x=x-2.0d0*xmax
          if (x .lt.-xmax) x=x+2.0d0*xmax
          rhowei=rhosmt(ix,iy,iz,is)*pwei(ix,iy,iz,na)*dxyz
          tptmp=tptmp+rhowei
          dpolx=dpolx+x*rhowei
          dpoly=dpoly+y*rhowei
          dpolz=dpolz+z*rhowei
          qpoxx=qpoxx+0.5d0*x*x*rhowei
          qpoyy=qpoyy+0.5d0*y*y*rhowei
          qpozz=qpozz+0.5d0*z*z*rhowei
          qpoxy=qpoxy+0.5d0*x*y*rhowei
          qpoyz=qpoyz+0.5d0*y*z*rhowei
          qpozx=qpozx+0.5d0*z*x*rhowei
        end do
        end do
        end do
      case(2)
!$omp do
        do iz=1,ncpz
        do iy=1,ncpy
        do ix=1,ncpx
          jx=myrx*ncpx+ix
          jy=myry*ncpy+iy
          jz=myrz*ncpz+iz
          x=(jx*dx-xmax-0.5d0*dx)-atomx
          y=(jy*dy-ymax-0.5d0*dy)-atomy
          z=(jz*dz-zmax-0.5d0*dz)-atomz
          if (x .gt. xmax) x=x-2.0d0*xmax
          if (x .lt.-xmax) x=x+2.0d0*xmax
          if (y .gt. ymax) y=y-2.0d0*ymax
          if (y .lt.-ymax) y=y+2.0d0*ymax
          rhowei=rhosmt(ix,iy,iz,is)*pwei(ix,iy,iz,na)*dxyz
          tptmp=tptmp+rhowei
          dpolx=dpolx+x*rhowei
          dpoly=dpoly+y*rhowei
          dpolz=dpolz+z*rhowei
          qpoxx=qpoxx+0.5d0*x*x*rhowei
          qpoyy=qpoyy+0.5d0*y*y*rhowei
          qpozz=qpozz+0.5d0*z*z*rhowei
          qpoxy=qpoxy+0.5d0*x*y*rhowei
          qpoyz=qpoyz+0.5d0*y*z*rhowei
          qpozx=qpozx+0.5d0*z*x*rhowei
        end do
        end do
        end do
      end select
    end do
!$omp critical
    atmpole( 1,na)=atmpole( 1,na)+tptmp
    atmpole( 2,na)=atmpole( 2,na)+dpolx
    atmpole( 3,na)=atmpole( 3,na)+dpoly
    atmpole( 4,na)=atmpole( 4,na)+dpolz
    atmpole( 5,na)=atmpole( 5,na)+qpoxx
    atmpole( 6,na)=atmpole( 6,na)+qpoyy
    atmpole( 7,na)=atmpole( 7,na)+qpozz
    atmpole( 8,na)=atmpole( 8,na)+qpoxy
    atmpole( 9,na)=atmpole( 9,na)+qpoyz
    atmpole(10,na)=atmpole(10,na)+qpozx
!$omp end critical
!$omp barrier
  end if
  if (nperi==3) then
    if (natprid(na)==1) then
!$omp do
      do ixyz=1,natinfd(na)
        rhowei=rhoaug3d(ixyz,napsd(na))*ddxyz
        tptmp=tptmp+rhowei
      end do
    end if
    do is=1,ns
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        rhowei=rhosmt(ix,iy,iz,is)*pwei(ix,iy,iz,na)*dxyz
        tptmp=tptmp+rhowei
      end do
      end do
      end do
    end do
!$omp critical
    atmpole( 1,na)=atmpole( 1,na)+tptmp
!$omp end critical
!$omp barrier
  end if
  end do ! na

!$omp single
    do na=1,natom
      do l=1,10
        call mpi_allreduce(atmpole(l,na),tptmp,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
        atmpole(l,na)=tptmp
      end do
    end do
!$omp end single
!$omp barrier

  return
end subroutine


end module
