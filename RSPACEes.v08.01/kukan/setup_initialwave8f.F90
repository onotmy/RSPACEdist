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
! **********  setup_initialwave8f.F90 11/19/2013-01  **********


module mod_setup_initialwave
implicit none

contains

subroutine setup_initialwave(nrc,natom,num_spe,nums,ncol,numk,neigmx,ncpx,ncpy,ncpz,zs_pre,pol_pre, & ! <
                            xmax,ymax,zmax,dx,dy,dz,                                                & ! <
                            mspe,indspe,                                                            & ! <
                            cp,atx,aty,atz,                                                         & ! X
                            svecre,sveccm,sval)                                                       ! >
use mod_mpi
use mod_stopp
implicit none
integer,   intent(in)   ::nrc,natom,num_spe,nums,ncol,numk,neigmx,ncpx,ncpy,ncpz
real*8,    intent(in)   ::zs_pre,pol_pre
real*8,    intent(in)   ::xmax,ymax,zmax,dx,dy,dz
integer,   intent(in)   ::mspe(num_spe),indspe(natom)
real*8,    intent(inout)::cp(8,num_spe)
real*8,    intent(inout)::atx(natom),aty(natom),atz(natom)
real*8,    intent(out)  ::svecre(ncpx*(1-nrc)+nrc,ncpy*(1-nrc)+nrc,ncpz*(1-nrc)+nrc &
                               ,neigmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk*(1-nrc)+nrc)
complex*16,intent(out)  ::sveccm(ncpx*nrc-nrc+1,ncpy*nrc-nrc+1,ncpz*nrc-nrc+1,neigmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1)
real*8,    intent(out)  ::sval(neigmx,nums+1-ncol,numk)

integer na,nk,ns,ix,iy,iz,jx,jy,jz,ntran,l
real*8 x,y,z,r,zs,yy1,yy2,yy3,yy4,yy5,yy6,yy7,yy8,yy9,tmp,snorm,snoall,tmpw
real*8 pi,sqfourpi,sq3fourpi,sq1516pi,sq0516pi,sq1504pi

  if (neigmx>9*natom) call stopp( &
   'setup_initialwave: Only 9 intitial wf per (spin and atom) are generated (neigmx is too large) !')

  pi=dacos(-1.0d0)
  sqfourpi=dsqrt(1.0d0/4.0d0/pi)
  sq3fourpi=dsqrt(3.0d0/4.0d0/pi)
  sq1516pi=dsqrt(15.0d0/16.0d0/pi)
  sq0516pi=dsqrt( 5.0d0/16.0d0/pi)
  sq1504pi=dsqrt(15.0d0/4.0d0/pi)

  call mpi_bcast(atx,natom,mpi_double_precision,0,mpicom_kpt,mpij)
  call mpi_bcast(aty,natom,mpi_double_precision,0,mpicom_kpt,mpij)
  call mpi_bcast(atz,natom,mpi_double_precision,0,mpicom_kpt,mpij)
  call mpi_bcast(cp,8*num_spe,mpi_double_precision,0,mpicom_kpt,mpij) ! only cp(1,:) needed

! **********  configration of scf data  **********
! ----------  zero clear  ----------
  svecre=0.0d0
  sveccm=dcmplx(0.0d0,0.0d0)
  sval=0.0d0
! ----------------------------------

  if (nrc==0) then
! ++++++++++  calculation of eigenvector  ++++++++++
    ntran=0
    do na=1,natom
      if ((mspe(indspe(na)) .ge. 21) .and. (mspe(indspe(na)) .le. 30)) ntran=1
      if ((mspe(indspe(na)) .ge. 39) .and. (mspe(indspe(na)) .le. 48)) ntran=1
      if ((mspe(indspe(na)) .ge. 57) .and. (mspe(indspe(na)) .le. 80)) ntran=1
    end do
    do nk=1,numk
    do ns=1,nums
    tmpw=1.0d0
    if (ns .eq. 2) tmpw=pol_pre
    do na=1,natom
       do iz=1,ncpz
       jz=myrz*ncpz+iz
       do iy=1,ncpy
       jy=myry*ncpy+iy
       do ix=1,ncpx
       jx=myrx*ncpx+ix
       x=jx*dx-xmax-0.5d0*dx-atx(na)
       y=jy*dy-ymax-0.5d0*dy-aty(na)
       z=jz*dz-zmax-0.5d0*dz-atz(na)
       if (x .gt. xmax) x=x-2.0d0*xmax
       if (x .lt.-xmax) x=x+2.0d0*xmax
       if (y .gt. ymax) y=y-2.0d0*ymax
       if (y .lt.-ymax) y=y+2.0d0*ymax
       if (z .gt. zmax) z=z-2.0d0*zmax
       if (z .lt.-zmax) z=z+2.0d0*zmax
       r=dsqrt(x*x+y*y+z*z)
       zs=cp(1,indspe(na))*zs_pre*tmpw
       yy1= sqfourpi
       yy2= sq3fourpi*x
       yy3= sq3fourpi*y
       yy4= sq3fourpi*z
       yy5= sq1516pi*(x*x-y*y)
       yy6=-sq1504pi*z*x
       yy7= sq0516pi*(3.0d0*z*z-r*r)
       yy8=-sq1504pi*y*z
       yy9= sq1504pi*x*y
       if ((neigmx .le. natom*4) .and. (ntran .eq. 0)) then
          if (na                 .le. neigmx) svecre(ix,iy,iz,na                ,ns,nk)=dexp(-zs*r)*yy1
          if (natom+(na-1)*3+1   .le. neigmx) svecre(ix,iy,iz,natom+(na-1)*3+1  ,ns,nk)=dexp(-zs*r)*yy2
          if (natom+(na-1)*3+2   .le. neigmx) svecre(ix,iy,iz,natom+(na-1)*3+2  ,ns,nk)=dexp(-zs*r)*yy3
          if (natom+(na-1)*3+3   .le. neigmx) svecre(ix,iy,iz,natom+(na-1)*3+3  ,ns,nk)=dexp(-zs*r)*yy4
       else
          if (na                 .le. neigmx) svecre(ix,iy,iz,na                ,ns,nk)=dexp(-zs*r)*yy1
          if (natom+(na-1)*3+1   .le. neigmx) svecre(ix,iy,iz,natom+(na-1)*3+1  ,ns,nk)=dexp(-zs*r)*yy2
          if (natom+(na-1)*3+2   .le. neigmx) svecre(ix,iy,iz,natom+(na-1)*3+2  ,ns,nk)=dexp(-zs*r)*yy3
          if (natom+(na-1)*3+3   .le. neigmx) svecre(ix,iy,iz,natom+(na-1)*3+3  ,ns,nk)=dexp(-zs*r)*yy4
          if (natom*4+(na-1)*5+1 .le. neigmx) svecre(ix,iy,iz,natom*4+(na-1)*5+1,ns,nk)=dexp(-zs*r)*yy5
          if (natom*4+(na-1)*5+2 .le. neigmx) svecre(ix,iy,iz,natom*4+(na-1)*5+2,ns,nk)=dexp(-zs*r)*yy6
          if (natom*4+(na-1)*5+3 .le. neigmx) svecre(ix,iy,iz,natom*4+(na-1)*5+3,ns,nk)=dexp(-zs*r)*yy7
          if (natom*4+(na-1)*5+4 .le. neigmx) svecre(ix,iy,iz,natom*4+(na-1)*5+4,ns,nk)=dexp(-zs*r)*yy8
          if (natom*4+(na-1)*5+5 .le. neigmx) svecre(ix,iy,iz,natom*4+(na-1)*5+5,ns,nk)=dexp(-zs*r)*yy9
       end if
       end do
       end do
       end do
    end do
    end do
    end do

!   __________  normalization __________
    do nk=1,numk
      do ns=1,nums
        do l=1,neigmx
          tmp=0.0d0
          snorm=0.0d0
          do iz=1,ncpz
          do iy=1,ncpy
          do ix=1,ncpx
            snorm=snorm+svecre(ix,iy,iz,l,ns,nk)**2
          end do
          end do
          end do
          tmp=tmp+snorm
          snorm=tmp
          call mpi_allreduce(snorm,snoall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
          snorm=snoall*dx*dy*dz
          snorm=1.0d0/dsqrt(snorm)
          do iz=1,ncpz
          do iy=1,ncpy
          do ix=1,ncpx
            svecre(ix,iy,iz,l,ns,nk)=svecre(ix,iy,iz,l,ns,nk)*snorm
          end do
          end do
          end do
        end do
      end do
    end do
! ************************************************
  else
! ++++++++++  calculation of eigenvector  ++++++++++
    ntran=0
    do na=1,natom
      if ((mspe(indspe(na)) .ge. 21) .and. (mspe(indspe(na)) .le. 30)) ntran=1
      if ((mspe(indspe(na)) .ge. 39) .and. (mspe(indspe(na)) .le. 48)) ntran=1
      if ((mspe(indspe(na)) .ge. 57) .and. (mspe(indspe(na)) .le. 80)) ntran=1
    end do
    do nk=1,numk
    do ns=1,nums
    tmpw=1.0d0
    if (ns .eq. 2) tmpw=pol_pre
    do na=1,natom
       do iz=1,ncpz
       jz=myrz*ncpz+iz
       do iy=1,ncpy
       jy=myry*ncpy+iy
       do ix=1,ncpx
       jx=myrx*ncpx+ix
       x=jx*dx-xmax-0.5d0*dx-atx(na)
       y=jy*dy-ymax-0.5d0*dy-aty(na)
       z=jz*dz-zmax-0.5d0*dz-atz(na)
       if (x .gt. xmax) x=x-2.0d0*xmax
       if (x .lt.-xmax) x=x+2.0d0*xmax
       if (y .gt. ymax) y=y-2.0d0*ymax
       if (y .lt.-ymax) y=y+2.0d0*ymax
       if (z .gt. zmax) z=z-2.0d0*zmax
       if (z .lt.-zmax) z=z+2.0d0*zmax
       r=dsqrt(x*x+y*y+z*z)
       zs=cp(1,indspe(na))*zs_pre*tmpw
       yy1= sqfourpi
       yy2= sq3fourpi*x
       yy3= sq3fourpi*y
       yy4= sq3fourpi*z
       yy5= sq1516pi*(x*x-y*y)
       yy6=-sq1504pi*z*x
       yy7= sq0516pi*(3.0d0*z*z-r*r)
       yy8=-sq1504pi*y*z
       yy9= sq1504pi*x*y
       if ((neigmx .le. natom*4) .and. (ntran .eq. 0)) then
          if (na                 .le. neigmx) sveccm(ix,iy,iz,na                ,ns,nk)=dcmplx(dexp(-zs*r)*yy1,0.0d0)
          if (natom+(na-1)*3+1   .le. neigmx) sveccm(ix,iy,iz,natom+(na-1)*3+1  ,ns,nk)=dcmplx(dexp(-zs*r)*yy2,0.0d0)
          if (natom+(na-1)*3+2   .le. neigmx) sveccm(ix,iy,iz,natom+(na-1)*3+2  ,ns,nk)=dcmplx(dexp(-zs*r)*yy3,0.0d0)
          if (natom+(na-1)*3+3   .le. neigmx) sveccm(ix,iy,iz,natom+(na-1)*3+3  ,ns,nk)=dcmplx(dexp(-zs*r)*yy4,0.0d0)
       else
          if (na                 .le. neigmx) sveccm(ix,iy,iz,na                ,ns,nk)=dcmplx(dexp(-zs*r)*yy1,0.0d0)
          if (natom+(na-1)*3+1   .le. neigmx) sveccm(ix,iy,iz,natom+(na-1)*3+1  ,ns,nk)=dcmplx(dexp(-zs*r)*yy2,0.0d0)
          if (natom+(na-1)*3+2   .le. neigmx) sveccm(ix,iy,iz,natom+(na-1)*3+2  ,ns,nk)=dcmplx(dexp(-zs*r)*yy3,0.0d0)
          if (natom+(na-1)*3+3   .le. neigmx) sveccm(ix,iy,iz,natom+(na-1)*3+3  ,ns,nk)=dcmplx(dexp(-zs*r)*yy4,0.0d0)
          if (natom*4+(na-1)*5+1 .le. neigmx) sveccm(ix,iy,iz,natom*4+(na-1)*5+1,ns,nk)=dcmplx(dexp(-zs*r)*yy5,0.0d0)
          if (natom*4+(na-1)*5+2 .le. neigmx) sveccm(ix,iy,iz,natom*4+(na-1)*5+2,ns,nk)=dcmplx(dexp(-zs*r)*yy6,0.0d0)
          if (natom*4+(na-1)*5+3 .le. neigmx) sveccm(ix,iy,iz,natom*4+(na-1)*5+3,ns,nk)=dcmplx(dexp(-zs*r)*yy7,0.0d0)
          if (natom*4+(na-1)*5+4 .le. neigmx) sveccm(ix,iy,iz,natom*4+(na-1)*5+4,ns,nk)=dcmplx(dexp(-zs*r)*yy8,0.0d0)
          if (natom*4+(na-1)*5+5 .le. neigmx) sveccm(ix,iy,iz,natom*4+(na-1)*5+5,ns,nk)=dcmplx(dexp(-zs*r)*yy9,0.0d0)
       end if
       end do
       end do
       end do
    end do
    end do
    end do

!   __________  normalization __________
    do nk=1,numk
      do ns=1,nums
        do l=1,neigmx
          tmp=0.0d0
          snorm=0.0d0
          do iz=1,ncpz
          do iy=1,ncpy
          do ix=1,ncpx
            snorm=snorm+abs(sveccm(ix,iy,iz,l,ns,nk))**2
          end do
          end do
          end do
          tmp=tmp+snorm
          snorm=tmp
          call mpi_allreduce(snorm,snoall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
          snorm=snoall*dx*dy*dz
          snorm=1.0d0/dsqrt(snorm)
          do iz=1,ncpz
          do iy=1,ncpy
          do ix=1,ncpx
            sveccm(ix,iy,iz,l,ns,nk)=sveccm(ix,iy,iz,l,ns,nk)*snorm
          end do
          end do
          end do
        end do
      end do
    end do
! ************************************************
  end if

  return
  end subroutine setup_initialwave


end module
