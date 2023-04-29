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
! **********  fuzzycell8f.f90 06/15/2014-01  **********

module mod_fuzzycell
implicit none

contains

subroutine fuzzycell(natom,nperi,ncpx,ncpy,ncpz,xmax,ymax,zmax,dx,dy,dz,atx,aty,atz,pwei)
use mod_mpi
implicit none
integer, intent(in)::natom,nperi
integer, intent(in)::ncpx,ncpy,ncpz
real*8, intent(in)::xmax,ymax,zmax
real*8, intent(in)::dx,dy,dz
real*8, intent(in)::atx(natom),aty(natom),atz(natom)
real*8, intent(out)::pwei(ncpx,ncpy,ncpz,natom)
real*8 x,y,z,rlen,x0,y0,z0,x1,y1,z1,r0,r1,amyu,fun1,fun2,fun3,scal
integer na,ix,iy,iz,na0,na1,jx0,jy0,jz0,jx,jy,jz,nx,ny,nz
real*8,allocatable::vtmp(:,:,:)
  allocate(vtmp(ncpx,ncpy,ncpz))

! ----------  clear  ----------
  pwei=1.0d0
! -----------------------------

  select case (nperi)
  case (0)
    nx=0
    ny=0
    nz=0
    do na1=1,natom
    do na0=1,natom
      do jz0=-nz,nz
      do jy0=-ny,ny
      do jx0=-nx,nx
        if (na0 .ne. na1) then
          x=atx(na0)+jx0*2.0d0*xmax-atx(na1)
          y=aty(na0)+jy0*2.0d0*ymax-aty(na1)
          z=atz(na0)+jz0*2.0d0*zmax-atz(na1)
          if ((dabs(x)-3.0d0*xmax.le.1.0d-8).and.(dabs(y)-3.0d0*ymax.le.1.0d-8).and.(dabs(z)-3.0d0*zmax.le.1.0d-8)) then
            rlen=dsqrt(x*x+y*y+z*z)
            do iz=1,ncpz
            do iy=1,ncpy
            do ix=1,ncpx
              jx=myrx*ncpx+ix
              jy=myry*ncpy+iy
              jz=myrz*ncpz+iz
              x0=jx*dx-xmax-0.5d0*dx-atx(na0)-jx0*xmax*2.0d0
              y0=jy*dy-ymax-0.5d0*dy-aty(na0)-jy0*ymax*2.0d0
              z0=jz*dz-zmax-0.5d0*dz-atz(na0)-jz0*zmax*2.0d0
              x1=jx*dx-xmax-0.5d0*dx-atx(na1)
              y1=jy*dy-ymax-0.5d0*dy-aty(na1)
              z1=jz*dz-zmax-0.5d0*dz-atz(na1)
              r0=dsqrt(x0*x0+y0*y0+z0*z0)
              r1=dsqrt(x1*x1+y1*y1+z1*z1)
              amyu=(r1-r0)/rlen
              if (amyu .gt.  1.0d0) amyu= 1.0d0
              if (amyu .lt. -1.0d0) amyu=-1.0d0
              fun1=1.5d0*amyu-0.5d0*amyu*amyu*amyu
              fun2=1.5d0*fun1-0.5d0*fun1*fun1*fun1
              fun3=1.5d0*fun2-0.5d0*fun2*fun2*fun2
              scal=0.5d0-0.5d0*fun3
              pwei(ix,iy,iz,na1)=pwei(ix,iy,iz,na1)*scal
            end do
            end do
            end do
          end if
        end if
      end do
      end do
      end do
    end do
    end do
  case (1)
    nx=3
    ny=0
    nz=0
    do na1=1,natom
    do na0=1,natom
      do jz0=-nz,nz
      do jy0=-ny,ny
      do jx0=-nx,nx
        if (na0 .ne. na1) then
          x=atx(na0)+jx0*2.0d0*xmax-atx(na1)
          y=aty(na0)+jy0*2.0d0*ymax-aty(na1)
          z=atz(na0)+jz0*2.0d0*zmax-atz(na1)
          if ((dabs(x)-3.0d0*xmax.le.1.0d-8).and.(dabs(y)-3.0d0*ymax.le.1.0d-8).and.(dabs(z)-3.0d0*zmax.le.1.0d-8)) then
            rlen=dsqrt(x*x+y*y+z*z)
            do iz=1,ncpz
            do iy=1,ncpy
            do ix=1,ncpx
              jx=myrx*ncpx+ix
              jy=myry*ncpy+iy
              jz=myrz*ncpz+iz
              x0=jx*dx-xmax-0.5d0*dx-atx(na0)-jx0*xmax*2.0d0
              y0=jy*dy-ymax-0.5d0*dy-aty(na0)-jy0*ymax*2.0d0
              z0=jz*dz-zmax-0.5d0*dz-atz(na0)-jz0*zmax*2.0d0
              x1=jx*dx-xmax-0.5d0*dx-atx(na1)
              y1=jy*dy-ymax-0.5d0*dy-aty(na1)
              z1=jz*dz-zmax-0.5d0*dz-atz(na1)
              if (x1 .gt. xmax) then
                x1=x1-2.0d0*xmax
                x0=x0-2.0d0*xmax
              end if
              if (x1 .lt.-xmax) then
                x1=x1+2.0d0*xmax
                x0=x0+2.0d0*xmax
              end if
              r0=dsqrt(x0*x0+y0*y0+z0*z0)
              r1=dsqrt(x1*x1+y1*y1+z1*z1)
              amyu=(r1-r0)/rlen
              if (amyu .gt.  1.0d0) amyu= 1.0d0
              if (amyu .lt. -1.0d0) amyu=-1.0d0
              fun1=1.5d0*amyu-0.5d0*amyu*amyu*amyu
              fun2=1.5d0*fun1-0.5d0*fun1*fun1*fun1
              fun3=1.5d0*fun2-0.5d0*fun2*fun2*fun2
              scal=0.5d0-0.5d0*fun3
              pwei(ix,iy,iz,na1)=pwei(ix,iy,iz,na1)*scal
            end do
            end do
            end do
          end if
        end if
      end do
      end do
      end do
    end do
    end do
  case (2)
    nx=3
    ny=3
    nz=0
    do na1=1,natom
    do na0=1,natom
      do jz0=-nz,nz
      do jy0=-ny,ny
      do jx0=-nx,nx
        if (na0 .ne. na1) then
          x=atx(na0)+jx0*2.0d0*xmax-atx(na1)
          y=aty(na0)+jy0*2.0d0*ymax-aty(na1)
          z=atz(na0)+jz0*2.0d0*zmax-atz(na1)
          if ((dabs(x)-3.0d0*xmax.le.1.0d-8).and.(dabs(y)-3.0d0*ymax.le.1.0d-8).and.(dabs(z)-3.0d0*zmax.le.1.0d-8)) then
            rlen=dsqrt(x*x+y*y+z*z)
            do iz=1,ncpz
            do iy=1,ncpy
            do ix=1,ncpx
              jx=myrx*ncpx+ix
              jy=myry*ncpy+iy
              jz=myrz*ncpz+iz
              x0=jx*dx-xmax-0.5d0*dx-atx(na0)-jx0*xmax*2.0d0
              y0=jy*dy-ymax-0.5d0*dy-aty(na0)-jy0*ymax*2.0d0
              z0=jz*dz-zmax-0.5d0*dz-atz(na0)-jz0*zmax*2.0d0
              x1=jx*dx-xmax-0.5d0*dx-atx(na1)
              y1=jy*dy-ymax-0.5d0*dy-aty(na1)
              z1=jz*dz-zmax-0.5d0*dz-atz(na1)
              if (x1 .gt. xmax) then
                x1=x1-2.0d0*xmax
                x0=x0-2.0d0*xmax
              end if
              if (x1 .lt.-xmax) then
                x1=x1+2.0d0*xmax
                x0=x0+2.0d0*xmax
              end if
              if (y1 .gt. ymax) then
                y1=y1-2.0d0*ymax
                y0=y0-2.0d0*ymax
              end if
              if (y1 .lt.-ymax) then
                y1=y1+2.0d0*ymax
                y0=y0+2.0d0*ymax
              end if
              r0=dsqrt(x0*x0+y0*y0+z0*z0)
              r1=dsqrt(x1*x1+y1*y1+z1*z1)
              amyu=(r1-r0)/rlen
              if (amyu .gt.  1.0d0) amyu= 1.0d0
              if (amyu .lt. -1.0d0) amyu=-1.0d0
              fun1=1.5d0*amyu-0.5d0*amyu*amyu*amyu
              fun2=1.5d0*fun1-0.5d0*fun1*fun1*fun1
              fun3=1.5d0*fun2-0.5d0*fun2*fun2*fun2
              scal=0.5d0-0.5d0*fun3
              pwei(ix,iy,iz,na1)=pwei(ix,iy,iz,na1)*scal
            end do
            end do
            end do
          end if
        end if
      end do
      end do
      end do
    end do
    end do
  case (3)
    nx=3
    ny=3
    nz=3
    do na1=1,natom
    do na0=1,natom
      do jz0=-nz,nz
      do jy0=-ny,ny
      do jx0=-nx,nx
        if (na0 .ne. na1) then
          x=atx(na0)+jx0*2.0d0*xmax-atx(na1)
          y=aty(na0)+jy0*2.0d0*ymax-aty(na1)
          z=atz(na0)+jz0*2.0d0*zmax-atz(na1)
          if ((dabs(x)-3.0d0*xmax.le.1.0d-8).and.(dabs(y)-3.0d0*ymax.le.1.0d-8).and.(dabs(z)-3.0d0*zmax.le.1.0d-8)) then
            rlen=dsqrt(x*x+y*y+z*z)
            do iz=1,ncpz
            do iy=1,ncpy
            do ix=1,ncpx
              jx=myrx*ncpx+ix
              jy=myry*ncpy+iy
              jz=myrz*ncpz+iz
              x0=jx*dx-xmax-0.5d0*dx-atx(na0)-jx0*xmax*2.0d0
              y0=jy*dy-ymax-0.5d0*dy-aty(na0)-jy0*ymax*2.0d0
              z0=jz*dz-zmax-0.5d0*dz-atz(na0)-jz0*zmax*2.0d0
              x1=jx*dx-xmax-0.5d0*dx-atx(na1)
              y1=jy*dy-ymax-0.5d0*dy-aty(na1)
              z1=jz*dz-zmax-0.5d0*dz-atz(na1)
              if (x1 .gt. xmax) then
                x1=x1-2.0d0*xmax
                x0=x0-2.0d0*xmax
              end if
              if (x1 .lt.-xmax) then
                x1=x1+2.0d0*xmax
                x0=x0+2.0d0*xmax
              end if
              if (y1 .gt. ymax) then
                y1=y1-2.0d0*ymax
                y0=y0-2.0d0*ymax
              end if
              if (y1 .lt.-ymax) then
                y1=y1+2.0d0*ymax
                y0=y0+2.0d0*ymax
              end if
              if (z1 .gt. zmax) then
                z1=z1-2.0d0*zmax
                z0=z0-2.0d0*zmax
              end if
              if (z1 .lt.-zmax) then
                z1=z1+2.0d0*zmax
                z0=z0+2.0d0*zmax
              end if
              r0=dsqrt(x0*x0+y0*y0+z0*z0)
              r1=dsqrt(x1*x1+y1*y1+z1*z1)
              amyu=(r1-r0)/rlen
              if (amyu .gt.  1.0d0) amyu= 1.0d0
              if (amyu .lt. -1.0d0) amyu=-1.0d0
              fun1=1.5d0*amyu-0.5d0*amyu*amyu*amyu
              fun2=1.5d0*fun1-0.5d0*fun1*fun1*fun1
              fun3=1.5d0*fun2-0.5d0*fun2*fun2*fun2
              scal=0.5d0-0.5d0*fun3
              pwei(ix,iy,iz,na1)=pwei(ix,iy,iz,na1)*scal
            end do
            end do
            end do
          end if
        end if
      end do
      end do
      end do
    end do
    end do
  end select

  do ix=1,ncpx*ncpy*ncpz
    vtmp(ix,1,1)=0.0d0
  end do
  do na=1,natom
    do ix=1,ncpx*ncpy*ncpz
       vtmp(ix,1,1)=vtmp(ix,1,1)+pwei(ix,1,1,na)
    end do
  end do
  do ix=1,ncpx*ncpy*ncpz
     vtmp(ix,1,1)=1.0d0/vtmp(ix,1,1)
  end do
  do na=1,natom
  do ix=1,ncpx*ncpy*ncpz
     pwei(ix,1,1,na)=pwei(ix,1,1,na)*vtmp(ix,1,1)
  end do
  end do

  deallocate(vtmp)
  return
end subroutine


end module
