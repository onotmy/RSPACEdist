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
! **********  trans8e.f90 11/19/2013-01  **********

      module mod_trans
      implicit none
      contains


      subroutine trans_c2d_smtcharge(nmesh,nspv,nperi,ndisp,nf,ncpx,ncpy,ncpz,rhosmt,rho_dense)
      use mod_overlap_interpolation, only:overlap_interpolation
      use mod_interpolation, only:interpolation
      implicit none
      integer, intent(in)::nmesh,nspv,nperi,ndisp,nf
      integer, intent(in)::ncpx,ncpy,ncpz
      real*8, intent(in)::rhosmt(ncpx,ncpy,ncpz,nspv)
      real*8, intent(out)::rho_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh,nspv)
      integer nf1,ns,ix,iy,iz,ncpx_d,ncpy_d,ncpz_d
      real*8,allocatable::vtmp_cc(:,:,:),vtmp_dd(:,:,:)
      ncpx_d=ncpx*nmesh
      ncpy_d=ncpy*nmesh
      ncpz_d=ncpz*nmesh
      nf1=nf-1

      allocate(vtmp_cc(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf))
      allocate(vtmp_dd(-(nf*nmesh-1):ncpx_d+nf*nmesh,-(nf*nmesh-1):ncpy_d+nf*nmesh,-(nf*nmesh-1):ncpz_d+nf*nmesh))
      vtmp_cc=0.0d0
      do ns=1,nspv
        do iz=1,ncpz
        do iy=1,ncpy
        do ix=1,ncpx
          vtmp_cc(ix,iy,iz)=rhosmt(ix,iy,iz,ns)
        end do
        end do
        end do
        call overlap_interpolation(ndisp,nperi,ncpx,ncpy,ncpz,nf,nf1,vtmp_cc,1)
        call interpolation(ncpx,ncpy,ncpz,nmesh,vtmp_dd,vtmp_cc,nperi,nf,1,ndisp)
        do iz=1,ncpz_d
        do iy=1,ncpy_d
        do ix=1,ncpx_d
          rho_dense(ix,iy,iz,ns)=vtmp_dd(ix,iy,iz)
        end do
        end do
        end do
      end do
      deallocate(vtmp_cc,vtmp_dd)

      return
      end subroutine trans_c2d_smtcharge


      subroutine trans_d2c_smtcharge( &
       nmesh,nperi,ndisp,nf,ncpx,ncpy,ncpz, & ! <
       rho_aug_dense,                       & ! <
       rho_coarse)                            ! >
      use mod_overlap_interpolation, only:overlap_interpolation
      use mod_interpolation, only:interpolation
      implicit none
      integer, intent(in)::nmesh,nperi,ndisp,nf
      integer, intent(in)::ncpx,ncpy,ncpz
      real*8, intent(in)::rho_aug_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
      real*8, intent(out)::rho_coarse(ncpx,ncpy,ncpz)
      integer nf1,ix,iy,iz,ncpx_d,ncpy_d,ncpz_d
      real*8,allocatable::vtmp_cc(:,:,:),vtmp_dd(:,:,:)
      ncpx_d=ncpx*nmesh
      ncpy_d=ncpy*nmesh
      ncpz_d=ncpz*nmesh
      nf1=nf-1

      allocate(vtmp_cc(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf))
      allocate(vtmp_dd(-(nf*nmesh-1):ncpx_d+nf*nmesh,-(nf*nmesh-1):ncpy_d+nf*nmesh,-(nf*nmesh-1):ncpz_d+nf*nmesh))
      vtmp_dd=0.0d0
      do iz=1,ncpz_d
      do iy=1,ncpy_d
      do ix=1,ncpx_d
        vtmp_dd(ix,iy,iz)=rho_aug_dense(ix,iy,iz)
      end do
      end do
      end do
      call interpolation(ncpx,ncpy,ncpz,nmesh,vtmp_dd,vtmp_cc,nperi,nf,2,ndisp)
      call overlap_interpolation(ndisp,nperi,ncpx,ncpy,ncpz,nf,nf1,vtmp_cc,2)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        rho_coarse(ix,iy,iz)=vtmp_cc(ix,iy,iz)
      end do
      end do
      end do
      deallocate(vtmp_cc,vtmp_dd)

      return
      end subroutine trans_d2c_smtcharge


      subroutine trans_d2c_vh(natom,num_spe,nmesh,nperi,ndisp,nf, & ! <
                              ncpx,ncpy,ncpz,                     & ! <
                              xmax,ymax,zmax,                     & ! <
                              indspe,                             & ! <
                              atmpole,cp,atx,aty,atz,             & ! <
                              vh_coarse,                          & ! >
                              vh_dense)                             ! <
      use mod_mpi
      use mod_overlap_interpolation, only:overlap_interpolation
      use mod_interpolation, only:interpolation
      implicit none
      integer, intent(in)::natom,num_spe,nmesh,nperi,ndisp,nf
      integer, intent(in)::ncpx,ncpy,ncpz
      real*8, intent(in)::xmax,ymax,zmax
      integer,intent(in)::indspe(natom)
      real*8, intent(in)::atmpole(10,natom)
      real*8, intent(in)::cp(8,num_spe)
      real*8, intent(in)::atx(natom),aty(natom),atz(natom)
      real*8, intent(out)::vh_coarse(ncpx,ncpy,ncpz)
      real*8, intent(in)::vh_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
      integer nf1,ix,iy,iz,ncpx_d,ncpy_d,ncpz_d
      real*8 ddx,ddy,ddz
      real*8,allocatable::vtmp_cc(:,:,:),vtmp_dd(:,:,:)
      ncpx_d=ncpx*nmesh
      ncpy_d=ncpy*nmesh
      ncpz_d=ncpz*nmesh
      nf1=nf-1
      ddx=2.0d0*xmax/(ncpx_d*nprocx)
      ddy=2.0d0*ymax/(ncpy_d*nprocy)
      ddz=2.0d0*zmax/(ncpz_d*nprocz)

      allocate(vtmp_cc(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf))
      allocate(vtmp_dd(-(nf*nmesh-1):ncpx_d+nf*nmesh,-(nf*nmesh-1):ncpy_d+nf*nmesh,-(nf*nmesh-1):ncpz_d+nf*nmesh))
      vtmp_dd=0.0d0
      do iz=1,ncpz_d
      do iy=1,ncpy_d
      do ix=1,ncpx_d
        vtmp_dd(ix,iy,iz)=vh_dense(ix,iy,iz)
      end do
      end do
      end do
      if (nperi<3) then 
        call trans_d2c_vh_01(natom,num_spe,nmesh,nperi,nf,  & ! <
                             ncpx_d,ncpy_d,ncpz_d,          & ! <
                             ddx,ddy,ddz,xmax,ymax,zmax,    & ! <
                             indspe,atmpole,cp,atx,aty,atz, & ! <
                             vtmp_dd)                         ! X
      endif 
      call interpolation(ncpx,ncpy,ncpz,nmesh,vtmp_dd,vtmp_cc,nperi,nf,2,ndisp)
      call overlap_interpolation(ndisp,nperi,ncpx,ncpy,ncpz,nf,nf1,vtmp_cc,2)
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        vh_coarse(ix,iy,iz)=vtmp_cc(ix,iy,iz)
      end do
      end do
      end do
      deallocate(vtmp_cc,vtmp_dd)
      return
      end subroutine trans_d2c_vh


      subroutine trans_d2c_vh_01(natom,num_spe,nmesh,nperi,nf,  & ! <
                                 ncpx_d,ncpy_d,ncpz_d,          & ! <
                                 ddx,ddy,ddz,xmax,ymax,zmax,    & ! <
                                 indspe,atmpole,cp,atx,aty,atz, & ! <
                                 vtmp_dd)                         ! X
      use mod_mpi
      implicit none
      integer,intent(in)::natom,num_spe,nmesh,nperi,nf,ncpx_d,ncpy_d,ncpz_d
      real*8, intent(in)::ddx,ddy,ddz,xmax,ymax,zmax
      integer,intent(in)::indspe(natom)
      real*8, intent(in)::atmpole(10,natom)
      real*8, intent(in)::cp(8,num_spe)
      real*8, intent(in)::atx(natom),aty(natom),atz(natom)
      real*8, intent(inout):: vtmp_dd(-(nf*nmesh-1):ncpx_d+nf*nmesh,-(nf*nmesh-1):ncpy_d+nf*nmesh,-(nf*nmesh-1):ncpz_d+nf*nmesh)
      integer nf2,nf3
      integer na,ix,iy,iz,jx,jy,jz,iix,iiy,iiz
      real*8 cp1,x,y,z,r,rin,r2in,r3in,r5in,derfdrr
      real*8 atptmp,edpolx,edpoly,edpolz,eqpoxx,eqpoyy,eqpozz,eqpoxy,eqpoyz,eqpozx,eqpoxy2,eqpoyz2,eqpozx2
      real*8 vmopm,vdpmx,vdpmy,vdpmz,vqmxx,vqmyy,vqmzz,vqmxy,vqmyz,vqmzx
      nf2=nf*nmesh
      nf3=nf2-1

      do na=1,natom
      cp1=cp(1,indspe(na))
      atptmp=atmpole( 1,na)/cp1
      edpolx=atmpole( 2,na)/cp1
      edpoly=atmpole( 3,na)/cp1
      edpolz=atmpole( 4,na)/cp1
      eqpoxx=atmpole( 5,na)/cp1
      eqpoyy=atmpole( 6,na)/cp1
      eqpozz=atmpole( 7,na)/cp1
      eqpoxy=atmpole( 8,na)/cp1
      eqpoyz=atmpole( 9,na)/cp1
      eqpozx=atmpole(10,na)/cp1
      eqpoxy2=2.0d0*eqpoxy
      eqpoyz2=2.0d0*eqpoyz
      eqpozx2=2.0d0*eqpozx

      if (nperi .eq. 0) then
      if (myrx .eq. 0) then
        do iz=-nf3,ncpz_d+nf2
        do iy=-nf3,ncpy_d+nf2
        do ix=-nf3,0
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=atx(na)-(jx*ddx-xmax-0.5d0*ddx)
          y=aty(na)-(jy*ddy-ymax-0.5d0*ddy)
          z=atz(na)-(jz*ddz-zmax-0.5d0*ddz)
          r=dsqrt(x*x+y*y+z*z)
          rin=1.0d0/r
          r2in=rin*rin
          r3in=r2in*rin
          r5in=r3in*r2in
          vmopm=-cp1*rin
          derfdrr=-r3in
          vdpmx= cp1*derfdrr*x
          vdpmy= cp1*derfdrr*y
          vdpmz= cp1*derfdrr*z
          vqmxx=-cp1*(3.0d0*x*x*r5in-r3in)
          vqmyy=-cp1*(3.0d0*y*y*r5in-r3in)
          vqmzz=-cp1*(3.0d0*z*z*r5in-r3in)
          vqmxy=-cp1*3.0d0*x*y*r5in
          vqmyz=-cp1*3.0d0*y*z*r5in
          vqmzx=-cp1*3.0d0*z*x*r5in
          vtmp_dd(ix,iy,iz)=vtmp_dd(ix,iy,iz)-atptmp*vmopm &
         -edpolx* vdpmx-edpoly* vdpmy-edpolz* vdpmz &
         -eqpoxx* vqmxx-eqpoyy* vqmyy-eqpozz* vqmzz &
         -eqpoxy2*vqmxy-eqpoyz2*vqmyz-eqpozx2*vqmzx
        end do
        end do
        end do
      end if
      if (myrx .eq. nprocx-1) then
        do iz=-nf3,ncpz_d+nf2
        do iy=-nf3,ncpy_d+nf2
        do ix=1,nf2
          iix=ix+ncpx_d
          jx=myrx*ncpx_d+ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=atx(na)-(jx*ddx-xmax-0.5d0*ddx)
          y=aty(na)-(jy*ddy-ymax-0.5d0*ddy)
          z=atz(na)-(jz*ddz-zmax-0.5d0*ddz)
          r=dsqrt(x*x+y*y+z*z)
          rin=1.0d0/r
          r2in=rin*rin
          r3in=r2in*rin
          r5in=r3in*r2in
          vmopm=-cp1*rin
          derfdrr=-r3in
          vdpmx= cp1*derfdrr*x
          vdpmy= cp1*derfdrr*y
          vdpmz= cp1*derfdrr*z
          vqmxx=-cp1*(3.0d0*x*x*r5in-r3in)
          vqmyy=-cp1*(3.0d0*y*y*r5in-r3in)
          vqmzz=-cp1*(3.0d0*z*z*r5in-r3in)
          vqmxy=-cp1*3.0d0*x*y*r5in
          vqmyz=-cp1*3.0d0*y*z*r5in
          vqmzx=-cp1*3.0d0*z*x*r5in
          vtmp_dd(iix,iy,iz)=vtmp_dd(iix,iy,iz)-atptmp*vmopm &
           -edpolx* vdpmx-edpoly* vdpmy-edpolz* vdpmz &
           -eqpoxx* vqmxx-eqpoyy* vqmyy-eqpozz* vqmzz &
           -eqpoxy2*vqmxy-eqpoyz2*vqmyz-eqpozx2*vqmzx
        end do
        end do
        end do
      end if
      end if

      if (nperi .le. 1) then
      if (myry .eq. 0) then
        do iz=-nf3,ncpz_d+nf2
        do iy=-nf3,0
        do ix=-nf3,ncpx_d+nf2
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=atx(na)-(jx*ddx-xmax-0.5d0*ddx)
          y=aty(na)-(jy*ddy-ymax-0.5d0*ddy)
          z=atz(na)-(jz*ddz-zmax-0.5d0*ddz)
          r=dsqrt(x*x+y*y+z*z)
          rin=1.0d0/r
          r2in=rin*rin
          r3in=r2in*rin
          r5in=r3in*r2in
          vmopm=-cp1*rin
          derfdrr=-r3in
          vdpmx= cp1*derfdrr*x
          vdpmy= cp1*derfdrr*y
          vdpmz= cp1*derfdrr*z
          vqmxx=-cp1*(3.0d0*x*x*r5in-r3in)
          vqmyy=-cp1*(3.0d0*y*y*r5in-r3in)
          vqmzz=-cp1*(3.0d0*z*z*r5in-r3in)
          vqmxy=-cp1*3.0d0*x*y*r5in
          vqmyz=-cp1*3.0d0*y*z*r5in
          vqmzx=-cp1*3.0d0*z*x*r5in
          vtmp_dd(ix,iy,iz)=vtmp_dd(ix,iy,iz)-atptmp*vmopm &
         -edpolx* vdpmx-edpoly* vdpmy-edpolz* vdpmz &
         -eqpoxx* vqmxx-eqpoyy* vqmyy-eqpozz* vqmzz &
         -eqpoxy2*vqmxy-eqpoyz2*vqmyz-eqpozx2*vqmzx
        end do
        end do
        end do
      end if
      if (myry .eq. nprocy-1) then
        do iz=-nf3,ncpz_d+nf2
        do iy=1,nf2
        do ix=-nf3,ncpx_d+nf2
          iiy=iy+ncpy_d
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=atx(na)-(jx*ddx-xmax-0.5d0*ddx)
          y=aty(na)-(jy*ddy-ymax-0.5d0*ddy)
          z=atz(na)-(jz*ddz-zmax-0.5d0*ddz)
          r=dsqrt(x*x+y*y+z*z)
          rin=1.0d0/r
          r2in=rin*rin
          r3in=r2in*rin
          r5in=r3in*r2in
          vmopm=-cp1*rin
          derfdrr=-r3in
          vdpmx= cp1*derfdrr*x
          vdpmy= cp1*derfdrr*y
          vdpmz= cp1*derfdrr*z
          vqmxx=-cp1*(3.0d0*x*x*r5in-r3in)
          vqmyy=-cp1*(3.0d0*y*y*r5in-r3in)
          vqmzz=-cp1*(3.0d0*z*z*r5in-r3in)
          vqmxy=-cp1*3.0d0*x*y*r5in
          vqmyz=-cp1*3.0d0*y*z*r5in
          vqmzx=-cp1*3.0d0*z*x*r5in
          vtmp_dd(ix,iiy,iz)=vtmp_dd(ix,iiy,iz)-atptmp*vmopm &
           -edpolx* vdpmx-edpoly* vdpmy-edpolz* vdpmz &
           -eqpoxx* vqmxx-eqpoyy* vqmyy-eqpozz* vqmzz &
           -eqpoxy2*vqmxy-eqpoyz2*vqmyz-eqpozx2*vqmzx
        end do
        end do
        end do
      end if
      end if

      if (nperi .le. 2) then
      if (myrz .eq. 0) then
        do iz=-nf3,0
        do iy=-nf3,ncpy_d+nf2
        do ix=-nf3,ncpx_d+nf2
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=atx(na)-(jx*ddx-xmax-0.5d0*ddx)
          y=aty(na)-(jy*ddy-ymax-0.5d0*ddy)
          z=atz(na)-(jz*ddz-zmax-0.5d0*ddz)
          r=dsqrt(x*x+y*y+z*z)
          rin=1.0d0/r
          r2in=rin*rin
          r3in=r2in*rin
          r5in=r3in*r2in
          vmopm=-cp1*rin
          derfdrr=-r3in
          vdpmx= cp1*derfdrr*x
          vdpmy= cp1*derfdrr*y
          vdpmz= cp1*derfdrr*z
          vqmxx=-cp1*(3.0d0*x*x*r5in-r3in)
          vqmyy=-cp1*(3.0d0*y*y*r5in-r3in)
          vqmzz=-cp1*(3.0d0*z*z*r5in-r3in)
          vqmxy=-cp1*3.0d0*x*y*r5in
          vqmyz=-cp1*3.0d0*y*z*r5in
          vqmzx=-cp1*3.0d0*z*x*r5in
          vtmp_dd(ix,iy,iz)=vtmp_dd(ix,iy,iz)-atptmp*vmopm &
         -edpolx* vdpmx-edpoly* vdpmy-edpolz* vdpmz &
         -eqpoxx* vqmxx-eqpoyy* vqmyy-eqpozz* vqmzz &
         -eqpoxy2*vqmxy-eqpoyz2*vqmyz-eqpozx2*vqmzx
        end do
        end do
        end do
      end if
      if(myrz .eq. nprocz-1) then
        do iz=1,nf2
        do iy=-nf3,ncpy_d+nf2
        do ix=-nf3,ncpx_d+nf2
          iiz=iz+ncpz_d
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+ncpz_d+iz
          x=atx(na)-(jx*ddx-xmax-0.5d0*ddx)
          y=aty(na)-(jy*ddy-ymax-0.5d0*ddy)
          z=atz(na)-(jz*ddz-zmax-0.5d0*ddz)
          r=dsqrt(x*x+y*y+z*z)
          rin=1.0d0/r
          r2in=rin*rin
          r3in=r2in*rin
          r5in=r3in*r2in
          vmopm=-cp1*rin
          derfdrr=-r3in
          vdpmx= cp1*derfdrr*x
          vdpmy= cp1*derfdrr*y
          vdpmz= cp1*derfdrr*z
          vqmxx=-cp1*(3.0d0*x*x*r5in-r3in)
          vqmyy=-cp1*(3.0d0*y*y*r5in-r3in)
          vqmzz=-cp1*(3.0d0*z*z*r5in-r3in)
          vqmxy=-cp1*3.0d0*x*y*r5in
          vqmyz=-cp1*3.0d0*y*z*r5in
          vqmzx=-cp1*3.0d0*z*x*r5in
          vtmp_dd(ix,iy,iiz)=vtmp_dd(ix,iy,iiz)-atptmp*vmopm &
           -edpolx* vdpmx-edpoly* vdpmy-edpolz* vdpmz &
           -eqpoxx* vqmxx-eqpoyy* vqmyy-eqpozz* vqmzz &
           -eqpoxy2*vqmxy-eqpoyz2*vqmyz-eqpozx2*vqmzx
        end do
        end do
        end do
      end if
      end if

      end do

      return
      end subroutine trans_d2c_vh_01


      subroutine trans_d2c_vxc(nmesh,nspv,nperi,ndisp,nf,ncpx,ncpy,ncpz,vx,vx_dense)
      use mod_overlap_interpolation, only:overlap_interpolation
      use mod_interpolation, only:interpolation
      implicit none
      integer, intent(in)::nmesh,nspv,nperi,ndisp,nf,ncpx,ncpy,ncpz
      real*8, intent(out)::vx(ncpx,ncpy,ncpz,nspv)
      real*8, intent(in)::vx_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh,nspv)
      integer nf1,ns,ix,iy,iz,ncpx_d,ncpy_d,ncpz_d
      real*8 tmp1
      real*8,allocatable::vtmp_cc(:,:,:),vtmp_dd(:,:,:)
      ncpx_d=ncpx*nmesh
      ncpy_d=ncpy*nmesh
      ncpz_d=ncpz*nmesh
      nf1=nf-1
      allocate(vtmp_cc(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf))
      allocate(vtmp_dd(-(nf*nmesh-1):ncpx_d+nf*nmesh,-(nf*nmesh-1):ncpy_d+nf*nmesh,-(nf*nmesh-1):ncpz_d+nf*nmesh))
      vtmp_dd=0.0d0
      tmp1=-1.865282921333129D-005 ! This value is obtained by ex. co. pot. of Vosko, Wilk, and Nusair.
      call trans_d2c_vxc_01(nperi,nmesh,nf,ncpx_d,ncpy_d,ncpz_d,tmp1,vtmp_dd)
      do ns=1,nspv
        do iz=1,ncpz_d
        do iy=1,ncpy_d
        do ix=1,ncpx_d
          vtmp_dd(ix,iy,iz)=vx_dense(ix,iy,iz,ns)
        end do
        end do
        end do
        call interpolation(ncpx,ncpy,ncpz,nmesh,vtmp_dd,vtmp_cc,nperi,nf,2,ndisp)
        call overlap_interpolation(ndisp,nperi,ncpx,ncpy,ncpz,nf,nf1,vtmp_cc,2)
        do iz=1,ncpz
        do iy=1,ncpy
        do ix=1,ncpx
          vx(ix,iy,iz,ns)=vtmp_cc(ix,iy,iz)
        end do
        end do
        end do
      end do
      deallocate(vtmp_cc,vtmp_dd)
      return
      end subroutine trans_d2c_vxc


      subroutine trans_d2c_vxc_01(nperi,nmesh,nf,ncpx_d,ncpy_d,ncpz_d,tmp1,vtmp_dd)
      use mod_mpi
      implicit none
      integer ix,iy,iz,iix,iiy,iiz,nf2,nf3
      integer, intent(in)::nperi,nmesh,nf,ncpx_d,ncpy_d,ncpz_d
      real*8, intent(in):: tmp1
      real*8, intent(inout)::vtmp_dd(-(nf*nmesh-1):ncpx_d+nf*nmesh,-(nf*nmesh-1):ncpy_d+nf*nmesh,-(nf*nmesh-1):ncpz_d+nf*nmesh)
      nf2=nf*nmesh
      nf3=nf2-1

      if (nperi .eq. 0) then
        if(myrx .eq. 0) then
          do ix=-nf3,0
            vtmp_dd(ix,:,:)=tmp1
          end do
        end if
        if(myrx .eq. nprocx-1) then
          do ix=1,nf2
            iix=ix+ncpx_d
            vtmp_dd(iix,:,:)=tmp1
          end do
        end if
      end if
      if (nperi .le. 1) then
        if(myry .eq. 0) then
          do iy=-nf3,0
            vtmp_dd(:,iy,:)=tmp1
          end do
        end if
        if(myry .eq. nprocy-1) then
          do iy=1,nf2
            iiy=iy+ncpy_d
            vtmp_dd(:,iiy,:)=tmp1
          end do
        end if
      end if
      if (nperi .le. 2) then
        if(myrz .eq. 0) then
          do iz=-nf3,0
            vtmp_dd(:,:,iz)=tmp1
          end do
        end if
        if(myrz .eq. nprocz-1) then
          do iz=1,nf2
            iiz=iz+ncpz_d
            vtmp_dd(:,:,iiz)=tmp1
          end do
        end if
      end if

      return
      end subroutine trans_d2c_vxc_01


      end module
