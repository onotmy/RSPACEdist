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
! **********  setup_pseudopotentials8f.f90 06/20/2018-01  **********

module mod_setup_pseudopot
implicit none
contains


subroutine setup_pseudopot( &
 num_spe,nradmx,nprmx,nprjmx,lmx,lpmx,npr,nrprj,nradct, & ! <
 radial,dradial,grwein,gmaxqp,gmax,sss0,akv0,awf,pwf,   & ! <
 nwexp,nprj,nlind,noind,sss,akv,rfac,wail,wpil)           ! >
use mod_mpi
implicit none
integer,intent(in) ::num_spe,nradmx,nprmx,nprjmx,lmx
integer,intent(in) ::lpmx(num_spe),npr(0:lmx-1,num_spe),nrprj(num_spe),nradct(num_spe)
real*8, intent(in) ::radial(nradmx,num_spe),dradial(nradmx,num_spe),grwein(15)
real*8, intent(in) ::gmaxqp,gmax(num_spe)
real*8, intent(in) ::sss0(nprmx*lmx,nprmx*lmx,num_spe),akv0(nprmx*lmx,nprmx*lmx,num_spe)
real*8, intent(in) ::awf(nradmx,nprmx*lmx,num_spe),pwf(nradmx,nprmx*lmx,num_spe)
integer,intent(out)::nwexp(num_spe),nprj(num_spe),nlind(nprjmx,num_spe),noind(nprjmx,num_spe)
real*8, intent(out)::sss(nprjmx,nprjmx,num_spe),akv(nprjmx,nprjmx,num_spe)
real*8, intent(out)::rfac(0:2*(lmx-1),num_spe)
real*8, intent(out)::wail((nprmx**2+nprmx)/2,lmx,num_spe),wpil((nprmx**2+nprmx)/2,lmx,num_spe)
integer:: ispe,ir,i,j,i1,i2,ie,lla,n
real*8 :: grmax,rcut,rfac0,r,dr,augbase,tmp

  if (myrank_glbl==0) then

  do ispe= 1,num_spe
! ================================================================
! ----------  assign matrix elements  ----------
    sss(:,:,ispe)=0.0d0
    akv(:,:,ispe)=0.0d0
    i=1
    j=1
    if (nrprj(ispe) .ge. 2) then
      sss(i,i,ispe)=sss0(j  ,j  ,ispe)
      akv(i,i,ispe)=akv0(j  ,j  ,ispe)
      i=i+1
      j=j+1
      if (npr(0,ispe) .eq. 2) then
        sss(i-1,i  ,ispe)=sss0(j-1,j  ,ispe)
        sss(i  ,i  ,ispe)=sss0(j  ,j  ,ispe)
        akv(i-1,i  ,ispe)=akv0(j-1,j  ,ispe)
        akv(i  ,i  ,ispe)=akv0(j  ,j  ,ispe)
        i=i+1
      end if
      j=j+1
    end if
    if (nrprj(ispe) .ge. 4) then
      sss(i  ,i  ,ispe)=sss0(j  ,j  ,ispe)
      sss(i+1,i+1,ispe)=sss0(j  ,j  ,ispe)
      sss(i+2,i+2,ispe)=sss0(j  ,j  ,ispe)
      akv(i  ,i  ,ispe)=akv0(j  ,j  ,ispe)
      akv(i+1,i+1,ispe)=akv0(j  ,j  ,ispe)
      akv(i+2,i+2,ispe)=akv0(j  ,j  ,ispe)
      i=i+3
      j=j+1
      if (npr(1,ispe) .eq. 2) then
        sss(i  ,i  ,ispe)=sss0(j  ,j  ,ispe)
        sss(i+1,i+1,ispe)=sss0(j  ,j  ,ispe)
        sss(i+2,i+2,ispe)=sss0(j  ,j  ,ispe)
        sss(i-3,i  ,ispe)=sss0(j-1,j  ,ispe)
        sss(i-2,i+1,ispe)=sss0(j-1,j  ,ispe)
        sss(i-1,i+2,ispe)=sss0(j-1,j  ,ispe)
        akv(i  ,i  ,ispe)=akv0(j  ,j  ,ispe)
        akv(i+1,i+1,ispe)=akv0(j  ,j  ,ispe)
        akv(i+2,i+2,ispe)=akv0(j  ,j  ,ispe)
        akv(i-3,i  ,ispe)=akv0(j-1,j  ,ispe)
        akv(i-2,i+1,ispe)=akv0(j-1,j  ,ispe)
        akv(i-1,i+2,ispe)=akv0(j-1,j  ,ispe)
        i=i+3
      end if
      j=j+1
    end if
    if (nrprj(ispe) .ge. 6) then
      sss(i  ,i  ,ispe)=sss0(j  ,j  ,ispe)
      sss(i+1,i+1,ispe)=sss0(j  ,j  ,ispe)
      sss(i+2,i+2,ispe)=sss0(j  ,j  ,ispe)
      sss(i+3,i+3,ispe)=sss0(j  ,j  ,ispe)
      sss(i+4,i+4,ispe)=sss0(j  ,j  ,ispe)
      akv(i  ,i  ,ispe)=akv0(j  ,j  ,ispe)
      akv(i+1,i+1,ispe)=akv0(j  ,j  ,ispe)
      akv(i+2,i+2,ispe)=akv0(j  ,j  ,ispe)
      akv(i+3,i+3,ispe)=akv0(j  ,j  ,ispe)
      akv(i+4,i+4,ispe)=akv0(j  ,j  ,ispe)
      i=i+5
      j=j+1
      if (npr(2,ispe) .eq. 2) then
        sss(i  ,i  ,ispe)=sss0(j  ,j  ,ispe)
        sss(i+1,i+1,ispe)=sss0(j  ,j  ,ispe)
        sss(i+2,i+2,ispe)=sss0(j  ,j  ,ispe)
        sss(i+3,i+3,ispe)=sss0(j  ,j  ,ispe)
        sss(i+4,i+4,ispe)=sss0(j  ,j  ,ispe)
        sss(i-5,i  ,ispe)=sss0(j-1,j  ,ispe)
        sss(i-4,i+1,ispe)=sss0(j-1,j  ,ispe)
        sss(i-3,i+2,ispe)=sss0(j-1,j  ,ispe)
        sss(i-2,i+3,ispe)=sss0(j-1,j  ,ispe)
        sss(i-1,i+4,ispe)=sss0(j-1,j  ,ispe)
        akv(i  ,i  ,ispe)=akv0(j  ,j  ,ispe)
        akv(i+1,i+1,ispe)=akv0(j  ,j  ,ispe)
        akv(i+2,i+2,ispe)=akv0(j  ,j  ,ispe)
        akv(i+3,i+3,ispe)=akv0(j  ,j  ,ispe)
        akv(i+4,i+4,ispe)=akv0(j  ,j  ,ispe)
        akv(i-5,i  ,ispe)=akv0(j-1,j  ,ispe)
        akv(i-4,i+1,ispe)=akv0(j-1,j  ,ispe)
        akv(i-3,i+2,ispe)=akv0(j-1,j  ,ispe)
        akv(i-2,i+3,ispe)=akv0(j-1,j  ,ispe)
        akv(i-1,i+4,ispe)=akv0(j-1,j  ,ispe)
      end if
    end if

    do j=1,nprjmx
      do i=1,j-1
        sss(j,i,ispe)=sss(i,j,ispe)
        akv(j,i,ispe)=akv(i,j,ispe)
      end do
    end do
! ------------------------------------------
    ! construct pseudo density
    ! --> estimate exponent for weinert-construction
    grmax=gmax(ispe)*gmaxqp*radial(nradct(ispe),ispe)
!    ie= transfer( minloc( abs(grwein(:)-grmax)) ,1)
    tmp=10.0d5
    do i=1,15
      if (dabs(grwein(i)-grmax) .lt. tmp) then
        ie=i
        tmp=dabs(grwein(i)-grmax)
      end if
    end do
    nwexp(ispe)= max(ie -1 ,0)
! ================================================================

! ==========  compute bases of compensation charge  ==========
! Although the coefficient is discussed in J. Math. Phys. 22, 2433 (1981),
! we do not have to use it because it is numerically normalized.
    rcut=radial(nradct(ispe),ispe)
    rfac0=1.0d0
    do lla=0,2*(lmx-1)
      n=nwexp(ispe)-lla
      tmp=0.0d0
      do ir=2,nradct(ispe)-1
        r=radial(ir,ispe)
        dr=dradial(ir,ispe)
        augbase=rfac0*(r/rcut)**lla*(1.0d0-(r/rcut)*(r/rcut))**n
        tmp=tmp+r*r*augbase*dr*r**lla
      end do
      rfac(lla,ispe)=rfac0/tmp
    end do
! ============================================================

! ==========  compute spherical momentum and orbital indexes  ==========
    i=0
    do lla=0,lpmx(ispe)
      do i2=1,npr(lla,ispe)
        do i1=1,2*lla+1
          i= i+1
          nlind(i,ispe)= lla**2 +i1
          noind(i,ispe)= lla*nprmx +i2
        enddo
      enddo
    enddo
    nprj(ispe)= i
! ======================================================================

! ==========  compute overlap of partial waves ===================
  do lla=0,lpmx(ispe)
    i= 0
    do i1=1,npr(lla,ispe)
      do i2=1,i1
        i= i+1
        wail(i,lla+1,ispe)= 0.0d0
        wpil(i,lla+1,ispe)= 0.0d0
        do ir= 2,nradct(ispe)
          wail(i,lla+1,ispe)= wail(i,lla+1,ispe)  &
           + dradial(ir,ispe)*awf(ir,lla*nprmx+i1,ispe)*awf(ir,lla*nprmx+i2,ispe)
          wpil(i,lla+1,ispe)= wpil(i,lla+1,ispe)  &
           + dradial(ir,ispe)*pwf(ir,lla*nprmx+i1,ispe)*pwf(ir,lla*nprmx+i2,ispe)
        enddo
      enddo
    enddo
  enddo
! ================================================================
  enddo ! ispe= 1,num_spe

  end if ! (myrank_glbl==0)

  call mpi_bcast(sss,nprjmx*nprjmx*num_spe,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(akv,nprjmx*nprjmx*num_spe,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(nprj,num_spe,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(nwexp,num_spe,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(rfac,(2*lmx-1)*num_spe,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(nlind,nprjmx*num_spe,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(noind,nprjmx*num_spe,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(wail,((nprmx**2+nprmx)/2)*lmx*num_spe,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(wpil,((nprmx**2+nprmx)/2)*lmx*num_spe,mpi_double_precision,0,mpicom_space,mpij)

end subroutine setup_pseudopot


subroutine setup_pseudopot_broadcast( &
 num_spe,nprmx,nprjmx,lmx,              & ! <
 ntyppp,nprj,nlind,noind,sss,wail,wpil)   ! X
use mod_mpi
implicit none
integer,intent(in)   ::num_spe,nprmx,nprjmx,lmx
integer,intent(inout)::ntyppp(num_spe) ! read in read_ppfile, not modified in setup_pseudopot
integer,intent(inout)::nprj(num_spe),nlind(nprjmx,num_spe),noind(nprjmx,num_spe)
real*8, intent(inout)::sss(nprjmx,nprjmx,num_spe)
real*8, intent(inout)::wail((nprmx**2+nprmx)/2,lmx,num_spe),wpil((nprmx**2+nprmx)/2,lmx,num_spe)

  call mpi_bcast(ntyppp,num_spe,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(nprj,num_spe,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(nlind,nprjmx*num_spe,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(noind,nprjmx*num_spe,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(sss,nprjmx*nprjmx*num_spe,mpi_double_precision,0,mpicom_kpt,mpij)
  call mpi_bcast(wail,((nprmx**2+nprmx)/2)*lmx*num_spe,mpi_double_precision,0,mpicom_kpt,mpij)
  call mpi_bcast(wpil,((nprmx**2+nprmx)/2)*lmx*num_spe,mpi_double_precision,0,mpicom_kpt,mpij)

end subroutine setup_pseudopot_broadcast

end module
