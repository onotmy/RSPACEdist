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
! **********  kscg8f.f90 06/20/2018-01  **********

module mod_kscg
implicit none
contains

subroutine kscg_r(natom,num_spe,neigmx,nums,nspv,numk,nperi,northo,ns,nk,nf,nkscg,nsdmax,    & ! <
                 nprecon_cg,nprjmx,num_list,num_ppcell,ndisp,                                & ! <
                 ncpx,ncpy,ncpz,                                                             & ! <
                 key_ortho_nocmpt_innerproduct,key_ortho_cmpt_innerproduct,                  & ! <
                 key_natpri_in,key_natpri_inps,key_natpri_out,                               & ! <
                 dx,dy,dz,epssd,                                                             & ! <
                 indspe,natpri,naps,nprj,natinf,lstvec2,latom,ntyppp,                        & ! <
                 dij,veff,vnlocp,sss,                                                        & ! <
                 ksconv,ksitmax,                                                             & ! X
                 sval,rspsep,svecre,ssvre,                                                   & ! X
                 residual_states,                                                            & ! >
                 hsvre)                                                                        ! >
use mod_mpi
use mod_stopp 
use mod_overlap_finitedifference, only:overlap_finitedifference_init,overlap_finitedifference_final, &
                                       overlap_finitedifference_r,overlap_fdcheck_r
use mod_nonlocaloperation, only:nonlocaloperation_r_01,nonlocaloperation_r_02
use mod_ksprecondition, only:ksprecondition
use mod_kslaplacian, only:kslaplacian_r,kslaplacian_initialize,kslaplacian_finalize
use mod_orthogonalization, only: &
 orthogonalization_r_01,orthogonalization_r_03,orthogonalization_r_04,orthogonalization_r_05,orthogonalization_r_06,multr
implicit none
integer, intent(in)   ::natom,num_spe,neigmx,nums,nspv,numk,nperi,northo,ns,nk,nf,nkscg,nsdmax
integer, intent(in)   ::nprecon_cg,nprjmx,num_list,num_ppcell,ndisp
integer, intent(in)   ::ncpx,ncpy,ncpz
integer, intent(in)   ::key_ortho_nocmpt_innerproduct,key_ortho_cmpt_innerproduct
integer, intent(in)   ::key_natpri_in,key_natpri_inps,key_natpri_out
real*8,  intent(in)   ::dx,dy,dz
real*8,  intent(in)   ::epssd
integer, intent(in)   ::indspe(natom),natpri(natom),naps(natom),nprj(num_spe)
integer, intent(in)   ::natinf(natom),lstvec2(num_list,num_ppcell),latom(natom),ntyppp(num_spe)
real*8,  intent(in)   ::dij(nprjmx,nprjmx,nums,natom)
real*8,  intent(in)   ::veff(ncpx,ncpy,ncpz,nspv)
real*8,  intent(in)   ::vnlocp(num_list,nprjmx,num_ppcell),sss(nprjmx,nprjmx,num_spe)
integer, intent(inout)::ksitmax
logical, intent(inout)::ksconv
real*8,  intent(inout)::sval(neigmx,nums,numk)
real*8,  intent(inout)::rspsep(nprjmx,num_ppcell,neigmx,nums,numk)
real*8,  intent(inout)::svecre(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8,  intent(inout)::ssvre(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8,  intent(out)  ::residual_states(neigmx)
real*8,  intent(out)  ::hsvre(ncpx,ncpy,ncpz,neigmx)
real*8 r2rr,valre,r1rr,valall,tmpall,vecnor,tmp,pap,sapre,papall,sapreall
integer na,loop,iiter,i,l,ix,iy,iz,iaps,key,nsv
real*8,allocatable::&
       rvecre(:,:,:) &
      ,pvecre(:,:,:) &
        ,pvre(:,:,:) &
        ,rvre(:,:,:) &
       ,spvre(:,:,:)
real*8,allocatable::vre(:,:,:)
real*8,allocatable::avre(:,:,:)
real*8,allocatable::asvre(:,:,:)
real*8,allocatable::workr(:)
real*8,allocatable::avreal(:),avrall(:)
real*8,allocatable::rrpsep(:,:),rppsep(:,:)     ! for isolated system
real*8,allocatable::rpsep(:,:)     ! for isolated system
real*8,allocatable::avr(:)
integer::nsdmax1,nsdmax2
logical::doreset,l_converged
integer nblas03,nblas04
  nblas03=500000
  nblas04=500000

  nsv= min(ns,nspv)

  allocate(avreal(neigmx),avrall(neigmx))
  allocate(rrpsep(nprjmx,num_ppcell),rppsep(nprjmx,num_ppcell))
  allocate(rpsep(nprjmx,num_ppcell))
  allocate(rvecre(ncpx,ncpy,ncpz) &
          ,pvecre(ncpx,ncpy,ncpz) &
            ,pvre(ncpx,ncpy,ncpz) &
            ,rvre(ncpx,ncpy,ncpz) &
           ,spvre(ncpx,ncpy,ncpz))
  allocate(vre(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf))
  allocate(avre(ncpx,ncpy,ncpz))
  allocate(asvre(ncpx,ncpy,ncpz))
  allocate(workr(1))
  allocate(avr(num_list))
  rrpsep=0.0d0

  call kslaplacian_initialize(nf)
  call overlap_finitedifference_init(ncpx,ncpy,ncpz,nf,1,0)

  nsdmax1= 10
  nsdmax2=200

  if (myrank_glbl==0) write(ndisp,fmt='(3(A,i5))') 'nsdmax=',nsdmax,' , nsdmax1=',nsdmax1,' , nsdmax2=',nsdmax2

  if (myrank_glbl==0) write(ndisp,*,err=9999) '  k ,spin,band,     eigen value    ,   residual norm   ,# of its.'

  do l=1,neigmx
! **********  orthonormalize wave function (svec)  **********
  if (northo .eq. 1) key=key_ortho_nocmpt_innerproduct
  if (northo .eq. 2) key=key_ortho_cmpt_innerproduct
  do loop=1,northo
    call orthogonalization_r_01(natom,neigmx,nprjmx,nums,numk,nk,ns,l,num_list,num_ppcell,num_spe,key,ndisp, & ! <
                                ncpx,ncpy,ncpz,                                                              & ! <
                                key_ortho_cmpt_innerproduct,                                                 & ! <
                                key_natpri_in,key_natpri_inps,                                               & ! <
                                dx,dy,dz,                                                                    & ! <
                                nprj,indspe,natpri,naps,natinf,lstvec2,latom,                                & ! <
                                vnlocp,sss,                                                                  & ! <
                                svecre,ssvre,rspsep)                                                           ! >
  end do
  !$omp parallel default(shared) private(i,ix,iy,iz)
  !$omp do
  do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)
    vre(-nf+i,-nf+1,-nf+1)=0.0d0
  end do
  !$omp do
  do i=1,ncpx*ncpy*ncpz
    hsvre(i,1,1,l)=0.0d0
  end do
  !$omp end do nowait
  !$omp do
  do iz=1,ncpz
  do iy=1,ncpy
  do ix=1,ncpx
    vre(ix,iy,iz)=svecre(ix,iy,iz,l,ns,nk)
  end do
  end do
  end do
  !$omp end do nowait
  !$omp barrier
  !$omp single
  call overlap_finitedifference_r(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,vre)
  !$omp end single
  !$omp do
  do i=1,nprjmx*num_ppcell
    rpsep(i,1)=rspsep(i,1,l,ns,nk)
  end do
  !$omp end do nowait
  !$omp barrier
  call nonlocaloperation_r_02(1,natom,num_spe,nums,nprjmx,ns,num_list,num_ppcell, & ! <
                              ncpx,ncpy,ncpz,                                     & ! <
                              key_natpri_in,key_natpri_inps,                      & ! <
                              dx,dy,dz,                                           & ! <
                              nprj,indspe,natinf,lstvec2,natpri,naps,             & ! <
                              vnlocp,rpsep,dij,                                   & ! <
                              hsvre(1,1,1,l),                                     & ! X
                              avr)                                                  ! W
  !$omp barrier
  !$omp single
  call overlap_fdcheck_r(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,vre)
  !$omp end single
  call kslaplacian_r(ncpx,ncpy,ncpz,neigmx,nf,dx,dy,dz,vre,veff(1,1,1,nsv),l,hsvre)
  !$omp end parallel
!     ***********************************************************

!     **********  conjugate gradient loop  **********
  iiter=0
  l_converged= .false.
  do while ((iiter < nsdmax).and.(.not. l_converged))
  iiter=iiter+1

  doreset= (iiter==1).or.(iiter==nsdmax1)
  if (nsdmax2>0)  doreset= doreset .or. (mod(iiter-1,nsdmax2)==0)
  if (doreset) then
    r2rr=1.0d0
  !$omp parallel default(shared) private(i)
  !$omp do
  do i=1,ncpx*ncpy*ncpz
    pvecre(i,1,1)=0.0d0
    rvecre(i,1,1)=0.0d0
    pvre(i,1,1)=0.0d0
    spvre(i,1,1)=0.0d0
  end do
  !$omp end do nowait
  !$omp do
  do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)
    vre(-nf+i,-nf+1,-nf+1)=0.0d0
  end do
  !$omp end do nowait
  !$omp do
  do i=1,nprjmx*num_ppcell
    rppsep(i,1)=0.0d0
  end do
  !$omp end do nowait
  !$omp end parallel
  end if

! ==========  compute eigenvalue and residual vector (rvec)  ==========
  valre=0.0d0
  r1rr=r2rr
  !$omp parallel default(shared)
  call kscg_r_01(ncpx*ncpy*ncpz,svecre(1,1,1,l,ns,nk),hsvre(1,1,1,l),valre)
  !$omp barrier
  !$omp single
  call mpi_allreduce(valre,valall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  valre= valall
  valre=valre*dx*dy*dz
  sval(l,ns,nk)=valre
  !$omp end single
  call kscg_r_02(ncpx*ncpy*ncpz,hsvre(1,1,1,l),ssvre(1,1,1,l,ns,nk),rvecre,sval(l,ns,nk))
  !$omp end parallel
! =====================================================================

! ==========  precondition for CG  ==========
  !$omp parallel default(shared) private(i,ix,iy,iz)
  if (nprecon_cg .eq. 1) then
    !$omp do
    do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)
      vre(-nf+i,-nf+1,-nf+1)=0.0d0
    end do
    !$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
      vre(ix,iy,iz)=rvecre(ix,iy,iz)
    end do
    end do
    end do
    !$omp end do nowait
    !$omp barrier
    !$omp single
    call overlap_finitedifference_r(nperi,ncpx,ncpy,ncpz,nf,nf-1,1,vre)
    call overlap_fdcheck_r(nperi,ncpx,ncpy,ncpz,nf,nf-1,1,vre)
    !$omp end single
    call ksprecondition(ncpx,ncpy,ncpz,nf,1,vre,rvre)
  else
    !$omp do
    do i=1,ncpx*ncpy*ncpz
      rvre(i,1,1)=rvecre(i,1,1)
    end do
    !$omp end do nowait
  end if
  !$omp end parallel
! ===========================================

! ==========  orthonormalize preconditioned residual vector  ==========
  do loop=1,northo
    !$omp parallel default(shared) private(i,ix,iy,iz)
    !$omp do
    do i=1,neigmx
      avreal(i)=0.0d0
    end do
    !$omp end do nowait
    !$omp barrier
    call orthogonalization_r_04(ncpx*ncpy*ncpz,1,neigmx,1,rvre,1,l,ssvre(1,1,1,1,ns,nk),avreal)
    !$omp barrier
    !$omp single
    call mpi_allreduce(avreal,avrall,l,mpi_double_precision,mpi_sum,mpicom_space,mpij)
    !$omp end single
    !$omp do
    do i=1,l
      avreal(i)=avrall(i)*dx*dy*dz
    end do
    !$omp end parallel
    if (l > nblas03) then
      call dgemv('n',ncpx*ncpy*ncpz,l,-1.0d0,svecre(1,1,1,1,ns,nk),ncpx*ncpy*ncpz &
          ,avreal,1,1.0d0,rvre,1)
    else
      !$omp parallel default(shared)
    call orthogonalization_r_03(ncpx*ncpy*ncpz,neigmx,l,rvre,svecre(1,1,1,1,ns,nk),avreal)
      !$omp end parallel
    end if
  end do
  !$omp parallel default(shared) private(i,ix,iy,iz)
  !$omp do
  do i=1,nprjmx*num_ppcell
    rpsep(i,1)=0.0d0
  end do
  call nonlocaloperation_r_01(natom,nprjmx,num_spe,num_ppcell,num_list,     & ! <
                              ncpx,ncpy,ncpz,                               & ! <
                              key_natpri_in,key_natpri_inps,                & ! <
                              nprj,indspe,natinf,natpri,naps,lstvec2,       & ! <
                              vnlocp,rvre,                                  & ! <
                              rpsep)                                          ! X
  !$omp barrier
  !$omp single
  do i=1,natom
    na=latom(i)
    if (natpri(na) /= key_natpri_out) then
      iaps=naps(na)
      call mpi_allreduce(rpsep(1,iaps),rrpsep(1,iaps),nprj(indspe(na)),mpi_double_precision,mpi_sum,mpicom_atom(na),mpij)
    end if
  end do
  !$omp end single
  !$omp do
  do i=1,ncpx*ncpy*ncpz
    asvre(i,1,1)=rvre(i,1,1)
  end do
  call nonlocaloperation_r_02(0,natom,num_spe,nums,nprjmx,1,num_list,num_ppcell, & ! <
                              ncpx,ncpy,ncpz,                                    & ! <
                              key_natpri_in,key_natpri_inps,                     & ! <
                              dx,dy,dz,                                          & ! <
                              nprj,indspe,natinf,lstvec2,natpri,naps,            & ! <
                              vnlocp,rrpsep,sss,                                 & ! <
                              asvre,                                             & ! X
                              avr)                                                 ! W
  !$omp barrier
  !$omp single
  workr(1)=0.0d0
  !$omp end single
  call orthogonalization_r_04(ncpx*ncpy*ncpz,1,1,1,rvecre(1,1,1),1,1,rvre(1,1,1),workr(1))
  !$omp end parallel
  call mpi_allreduce(workr(1),tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  r2rr= tmpall
  r2rr=r2rr*dx*dy*dz
  residual_states(l)=r2rr
! =====================================================================

! ==========  compute conjugate vector (pvec)  ==========
  !$omp parallel default(shared) private(i)
  if (nkscg .eq. 1) then
    call kscg_r_03(ncpx*ncpy*ncpz,pvecre,rvre,pvre,spvre,asvre,r1rr,r2rr)
    !$omp barrier
    !$omp do
    do i=1,nprjmx*num_ppcell
      rppsep(i,1)=rrpsep(i,1)+r2rr/r1rr*rppsep(i,1)
    end do
    !$omp end do nowait
  else
    !$omp do
    do i=1,ncpx*ncpy*ncpz
      pvecre(i,1,1)=rvre(i,1,1)
      spvre(i,1,1)=asvre(i,1,1)
    end do
    !$omp end do nowait
    !$omp do
    do i=1,nprjmx*num_ppcell
      rppsep(i,1)=rrpsep(i,1)
    end do
    !$omp end do nowait
  end if
  !$omp end parallel
! =======================================================

! ==========  orthonormalize conjugate vector (pvec)  ==========
  !$omp parallel default(shared) private(i)
  !$omp do
  do i=1,nprjmx*num_ppcell
    rpsep(i,1)=rppsep(i,1)
  end do
  !$omp end do nowait
  !$omp do
  do i=1,ncpx*ncpy*ncpz
    asvre(i,1,1)=spvre(i,1,1)
  end do
  !$omp end do nowait
  !$omp end parallel
  do loop=1,northo
    !$omp parallel default(shared) private(i)
    !$omp do
    do i=1,neigmx
      avreal(i)=0.0d0
    end do
    !$omp end parallel
    if (l > nblas04) then
      call dgemv('t',ncpx*ncpy*ncpz,l,1.0d0,ssvre(1,1,1,1,ns,nk),ncpx*ncpy*ncpz &
              ,pvecre,1,0.0d0,avreal,1)
    else
      !$omp parallel default(shared)
      call orthogonalization_r_04(ncpx*ncpy*ncpz,1,neigmx,1,pvecre,1,l,ssvre(1,1,1,1,ns,nk),avreal)
      !$omp end parallel
    end if
    call mpi_allreduce(avreal,avrall,l,mpi_double_precision,mpi_sum,mpicom_space,mpij)
    !$omp parallel default(shared) private(i)
    !$omp do
    do i=1,l
      avreal(i)=avrall(i)*dx*dy*dz
    end do
    !$omp end parallel
    if (l > nblas03) then
      call dgemv('n',ncpx*ncpy*ncpz,l,-1.0d0,svecre(1,1,1,1,ns,nk),ncpx*ncpy*ncpz &
          ,avreal,1,1.0d0,pvecre,1)
      call dgemv('n',ncpx*ncpy*ncpz,l,-1.0d0,ssvre(1,1,1,1,ns,nk),ncpx*ncpy*ncpz &
          ,avreal,1,1.0d0,asvre,1)
      call dgemv('n',nprjmx*num_ppcell,l,-1.0d0,rspsep(1,1,1,ns,nk),nprjmx*num_ppcell &
          ,avreal,1,1.0d0,rpsep,1)
    else
      !$omp parallel default(shared)
      call orthogonalization_r_03(ncpx*ncpy*ncpz,neigmx,l,pvecre,svecre(1,1,1,1,ns,nk),avreal)
      call orthogonalization_r_03(ncpx*ncpy*ncpz,neigmx,l,asvre,ssvre(1,1,1,1,ns,nk),avreal)
      call multr(nprjmx*num_ppcell,1,l,rpsep,rspsep(1,1,1,ns,nk),avreal)
      !$omp end parallel
    end if
  end do
  vecnor=0.0d0
  !$omp parallel default(shared) private(i)
  call orthogonalization_r_05(ncpx*ncpy*ncpz,pvecre,asvre,vecnor)
  !$omp barrier
  !$omp single
  call mpi_allreduce(vecnor,tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  vecnor= tmpall
  vecnor=vecnor*dx*dy*dz
  if (vecnor<=0.0d0) call stopp('kscg_r: pvec has negative norm')
  tmp=1.0d0/dsqrt(vecnor)
  !$omp end single
  call orthogonalization_r_06(ncpx*ncpy*ncpz,pvecre,tmp)
  call orthogonalization_r_06(ncpx*ncpy*ncpz,asvre,tmp)
  !$omp do
  do i=1,nprjmx*num_ppcell
    rpsep(i,1)=rpsep(i,1)*tmp
  end do
  !$omp end do nowait
  !$omp end parallel
! ==============================================================

! ==========  update eigenvector (svec)  ==========
! --  rpsep is computed by the preceding routine
  !$omp parallel default(shared) private(i,ix,iy,iz)
  !$omp do
  do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)
    vre(-nf+i,-nf+1,-nf+1)=0.0d0
  end do
  !$omp do
  do i=1,ncpx*ncpy*ncpz
    avre(i,1,1)=0.0d0
  end do
  !$omp end do nowait
  !$omp do
  do iz=1,ncpz
  do iy=1,ncpy
  do ix=1,ncpx
    vre(ix,iy,iz)=pvecre(ix,iy,iz)
  end do
  end do
  end do
  !$omp single
  call overlap_finitedifference_r(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,vre)
  !$omp end single
  call nonlocaloperation_r_02(1,natom,num_spe,nums,nprjmx,ns,num_list,num_ppcell, & ! <
                              ncpx,ncpy,ncpz,                                     & ! <
                              key_natpri_in,key_natpri_inps,                      & ! <
                              dx,dy,dz,                                           & ! <
                              nprj,indspe,natinf,lstvec2,natpri,naps,             & ! <
                              vnlocp,rpsep,dij,                                   & ! <
                              avre,                                               & ! X
                              avr)                                                  ! W
  !$omp barrier
  !$omp single
  call overlap_fdcheck_r(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,vre)
  !$omp end single
  call kslaplacian_r(ncpx,ncpy,ncpz,1,nf,dx,dy,dz,vre,veff(1,1,1,nsv),1,avre)
  !$omp barrier
  !$omp single
  pap=0.0d0
  sapre=0.0d0
  !$omp end single
  call kscg_r_04(ncpx*ncpy*ncpz,pvecre,svecre(1,1,1,l,ns,nk),avre,pap,sapre)
  !$omp barrier
  !$omp single
  call mpi_allreduce(pap,papall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  call mpi_allreduce(sapre,sapreall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  pap=papall*dx*dy*dz
  sapre=sapreall*dx*dy*dz
  !$omp end single
  call kscg_r_05(ncpx*ncpy*ncpz,nprjmx*num_ppcell,pvecre,hsvre(1,1,1,l),avre,asvre  &
                 ,svecre(1,1,1,l,ns,nk),ssvre(1,1,1,l,ns,nk),rpsep,rspsep(1,1,l,ns,nk) &
                 ,pap,sapre,sval(l,ns,nk),l_converged)
  !$omp end parallel
! =================================================

  l_converged= (dsqrt(residual_states(l)) < epssd)
  ksconv= ksconv .and. (l_converged .or. (iiter < nsdmax))

  end do ! iiter

  if (myrank_glbl==0) write (ndisp,'(i4,1x,i4,1x,i4,1x,2e20.10,1x,i3)',err=9999) &
    nk,ns,l,sval(l,ns,nk),dsqrt(residual_states(l)),iiter
  ksitmax= max(iiter,ksitmax)
 end do ! l
! ***********************************************

  call kslaplacian_finalize
  call overlap_finitedifference_final

  deallocate(vre)
  deallocate(avre)
  deallocate(asvre)
  deallocate(workr)
  deallocate(avreal,avrall)
  deallocate(rrpsep,rppsep)
  deallocate(rpsep)
  deallocate(rvecre,pvecre,pvre,rvre,spvre)
  deallocate(avr)

  return
9999 continue
  call mpi_abort(mpi_comm_world)
end subroutine kscg_r


subroutine kscg_c(nso,natom,num_spe,neigmx,nums,ncol,nspv,numk,nperi,northo,ns1,nk,nf,nkscg,nsdmax, & ! <
                 nprecon_cg,nprjmx,num_list,num_ppcell,ndisp,                                       & ! <
                 ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                                               & ! <
                 key_ortho_nocmpt_innerproduct,key_ortho_cmpt_innerproduct,                         & ! <
                 key_natpri_in,key_natpri_inps,key_natpri_out,key_soc_calc,                         & ! <
                 dx,dy,dz,skpx,skpy,skpz,epssd,                                                     & ! <
                 indspe,natpri,naps,nprj,natinf,lstvec2,latom,natsoc,ntyppp,                        & ! <
                 lstx,lsty,lstz,natx,naty,natz,                                                     & ! <
                 dij,dijsoc,veff,vnlocp,sss,                                                        & ! <
                 ksconv,ksitmax,                                                                    & ! X
                 sval,cspsep,sveccm,ssvcm,                                                          & ! X
                 residual_states,                                                                   & ! >
                 hsvcm)                                                                               ! >
use mod_mpi
use mod_stopp 
use mod_overlap_finitedifference, only:overlap_finitedifference_init,overlap_finitedifference_final, &
                                       overlap_finitedifference_c,overlap_fdcheck_c
use mod_nonlocaloperation, only:nonlocaloperation_c_01,nonlocaloperation_c_02
use mod_ksprecondition, only:ksprecondition_c
use mod_kslaplacian, only:kslaplacian_c,kslaplacian_initialize,kslaplacian_finalize
use mod_orthogonalization, only: &
 orthogonalization_c_01,orthogonalization_c_03,orthogonalization_c_04, &
 orthogonalization_c_05,orthogonalization_c_06,multc
implicit none
integer,    intent(in)   ::nso,natom,num_spe,neigmx,nums,ncol,nspv,numk,nperi,northo,ns1,nk,nf,nkscg,nsdmax
integer,    intent(in)   ::nprecon_cg,nprjmx,num_list,num_ppcell,ndisp
integer,    intent(in)   ::ncpx,ncpy,ncpz,npxmax,npymax,npzmax
integer,    intent(in)   ::key_ortho_nocmpt_innerproduct,key_ortho_cmpt_innerproduct
integer,    intent(in)   ::key_natpri_in,key_natpri_inps,key_natpri_out,key_soc_calc
real*8,     intent(in)   ::dx,dy,dz
real*8,     intent(in)   ::skpx,skpy,skpz
real*8,     intent(in)   ::epssd
integer,    intent(in)   ::indspe(natom),natpri(natom),naps(natom),nprj(num_spe)
integer,    intent(in)   ::natinf(natom),lstvec2(num_list,num_ppcell),latom(natom),natsoc(natom),ntyppp(num_spe)
integer,    intent(in)   ::lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell)
integer,    intent(in)   ::natx(natom),naty(natom),natz(natom)
real*8,     intent(in)   ::dij(nprjmx,nprjmx,nums*ncol,natom)
real*8,     intent(in)   ::dijsoc(nprjmx*nso-nso+1,nprjmx*nso-nso+1,3*nso-nso+1,natom*nso-nso+1)
real*8,     intent(in)   ::veff(ncpx,ncpy,ncpz,nspv)
real*8,     intent(in)   ::vnlocp(num_list,nprjmx,num_ppcell),sss(nprjmx,nprjmx,num_spe)
integer,    intent(inout)::ksitmax
logical,    intent(inout)::ksconv
real*8,     intent(inout)::sval(neigmx,nums+1-ncol,numk)
complex*16, intent(inout)::cspsep(nprjmx,num_ppcell,neigmx,nums,numk)
complex*16, intent(inout)::sveccm(ncpx,ncpy,ncpz,neigmx,nums,numk)
complex*16, intent(inout)::ssvcm(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8,     intent(out)  ::residual_states(neigmx)
complex*16, intent(out)  ::hsvcm(ncpx,ncpy,ncpz,neigmx,ncol)
real*8 valre,valall,tmpall,vecnor,tmp,pap,papall
real*8 xmax,ymax,zmax
integer na,loop,iiter,l,ix,iy,iz,i,iaps,key,ns2,nsv,nvef
complex*16 sapcm,sapcmall
complex*16 r1rcm,r2rcm,r2rcmall
complex*16,allocatable::&
       rveccm(:,:,:,:) &
      ,pveccm(:,:,:,:) &
        ,pvcm(:,:,:,:) &
        ,rvcm(:,:,:,:) &
       ,spvcm(:,:,:,:)
complex*16,allocatable::vcm(:,:,:,:)
complex*16,allocatable::avcm(:,:,:,:)
complex*16,allocatable::asvcm(:,:,:,:)
complex*16,allocatable::workc(:,:,:,:)
complex*16,allocatable::avcmplx(:),avcmplxall(:)
complex*16,allocatable::vcccm(:,:)
complex*16,allocatable::avc(:,:)
complex*16,allocatable::crpsep(:,:,:),cppsep(:,:,:) ! for periodic or noncollinear system
complex*16,allocatable::cpsep(:,:,:) ! for periodic or noncollinear system
complex*16,allocatable::workcc(:)
logical   :: doreset,l_converged
integer   :: nsdmax1,nsdmax2
integer nblas03,nblas04
complex*16 calph,cbeta
  nblas03=500000
  nblas04=500000
  calph=dcmplx(1.0d0,0.0d0)
  cbeta=dcmplx(0.0d0,0.0d0)

  nvef=1
  if (ncol==2) nvef=nspv
  nsv= min(ns1,nspv)

  xmax=(ncpx*nprocx)*dx*0.5d0
  ymax=(ncpy*nprocy)*dy*0.5d0
  zmax=(ncpz*nprocz)*dz*0.5d0
  allocate(avcmplx(neigmx),avcmplxall(neigmx))
  allocate(crpsep(nprjmx,num_ppcell,ncol),cppsep(nprjmx,num_ppcell,ncol))
  allocate(cpsep(nprjmx,num_ppcell,ncol))
  allocate(rveccm(ncpx,ncpy,ncpz,ncol) &
          ,pveccm(ncpx,ncpy,ncpz,ncol) &
            ,pvcm(ncpx,ncpy,ncpz,ncol) &
            ,rvcm(ncpx,ncpy,ncpz,ncol) &
           ,spvcm(ncpx,ncpy,ncpz,ncol))
  allocate(vcm(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf,ncol))
  allocate(avcm(ncpx,ncpy,ncpz,ncol))
  allocate(asvcm(ncpx,ncpy,ncpz,ncol))
  allocate(workc(ncpx,ncpy,ncpz,ncol))
  allocate(vcccm(num_list,ncol),avc(num_list,ncol))
  allocate(workcc(1))
  crpsep=dcmplx(0.0d0,0.0d0)

  call kslaplacian_initialize(nf)
  call overlap_finitedifference_init(ncpx,ncpy,ncpz,nf,ncol,1)

  nsdmax1= 10
  nsdmax2= 200

  if (myrank_glbl==0) write(ndisp,fmt='(3(A,i5))') 'nsdmax=',nsdmax,' , nsdmax1=',nsdmax1,' , nsdmax2=',nsdmax2

  if (myrank_glbl==0) write(ndisp,*,err=9999) '  k ,spin,band,     eigen value    ,   residual norm   ,# of its.'

  do l=1,neigmx

! **********  orthonormalize wave function (svec)  **********
  if (northo .eq. 1) key=key_ortho_nocmpt_innerproduct
  if (northo .eq. 2) key=key_ortho_cmpt_innerproduct
  do loop=1,northo
    call orthogonalization_c_01(natom,neigmx,nprjmx,nums,ncol,numk,nk,ns1,l,num_list,num_ppcell,num_spe,key,ndisp, & ! <
                                ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                                               & ! <
                                key_ortho_cmpt_innerproduct,                                                       & ! <
                                key_natpri_in,key_natpri_inps,                                                     & ! <
                                dx,dy,dz,skpx,skpy,skpz,                                                           & ! <
                                nprj,indspe,natpri,naps,natinf,lstvec2,latom,                                      & ! <
                                lstx,lsty,lstz,natx,naty,natz,                                                     & ! <
                                vnlocp,sss,                                                                        & ! <
                                sveccm,ssvcm,cspsep)                                                                 ! X
  end do
  !$omp parallel default(shared) private(i,ix,iy,iz)
  !$omp do
  do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)*ncol
    vcm(-nf+i,-nf+1,-nf+1,1)=dcmplx(0.0d0,0.0d0)
  end do
  do ns2= 1,ncol
    !$omp do
    do i=1,ncpx*ncpy*ncpz
      hsvcm(i,1,1,l,ns2)=dcmplx(0.0d0,0.0d0)
    end do
    !$omp end do nowait
    !$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
      vcm(ix,iy,iz,ns2)=sveccm(ix,iy,iz,l,ns1+ns2-1,nk)
    end do
    end do
    end do
    !$omp end do nowait
  end do ! ns2
  !$omp barrier
  !$omp single
  call overlap_finitedifference_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,ncol,vcm)
  !$omp end single
  do ns2= 1,ncol
    !$omp do
    do i=1,nprjmx*num_ppcell
      cpsep(i,1,ns2)= cspsep(i,1,l,ns1+ns2-1,nk)
    end do
    !$omp end do nowait
  end do
  !$omp barrier
  call nonlocaloperation_c_02(1,nso,natom,num_spe,neigmx,nums,nprjmx,ncol,ns1,l,num_list,num_ppcell, & ! <
                              ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                                   & ! <
                              key_natpri_in,key_natpri_inps,key_soc_calc,                            & ! <
                              dx,dy,dz,skpx,skpy,skpz,                                               & ! <
                              nprj,indspe,natinf,lstvec2,natpri,naps,natsoc,                         & ! <
                              lstx,lsty,lstz,natx,naty,natz,                                         & ! <
                              vnlocp,cpsep,dij,dijsoc,                                               & ! <
                              hsvcm,                                                                 & ! X
                              vcccm,avc)                                                               ! W
  !$omp barrier
  !$omp single
  call overlap_fdcheck_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,ncol,vcm)
  !$omp end single
  call kslaplacian_c(ncpx,ncpy,ncpz,neigmx,ncol,nvef,nf,dx,dy,dz,xmax,ymax,zmax,skpx,skpy,skpz &
                     ,vcm,veff(1,1,1,nsv),l,hsvcm,workc)
  !$omp end parallel
!     ***********************************************************

!     **********  conjugate gradient loop  **********
  iiter=0
  l_converged= .false.
  do while ((iiter < nsdmax).and.(.not. l_converged))
  iiter=iiter+1

  doreset= (iiter==1).or.(iiter==nsdmax1)
  if (nsdmax2>0)  doreset= doreset .or. (mod(iiter-1,nsdmax2)==0)
  if (doreset) then
    r2rcm=dcmplx(1.0d0,0.0d0)
    !$omp parallel default(shared) private(i)
    !$omp do
    do i=1,ncpx*ncpy*ncpz*ncol
      pveccm(i,1,1,1)=dcmplx(0.0d0,0.0d0)
      rveccm(i,1,1,1)=dcmplx(0.0d0,0.0d0)
      pvcm(i,1,1,1)=dcmplx(0.0d0,0.0d0)
      spvcm(i,1,1,1)=dcmplx(0.0d0,0.0d0)
    end do
    !$omp end do nowait
    !$omp do
    do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)*ncol
      vcm(-nf+i,-nf+1,-nf+1,1)=dcmplx(0.0d0,0.0d0)
    end do
    !$omp end do nowait
    !$omp do
    do i=1,nprjmx*num_ppcell*ncol
      cppsep(i,1,1)=dcmplx(0.0d0,0.0d0)
    end do
    !$omp end do nowait
    !$omp end parallel
  end if

! ==========  compute eigenvalue and residual vector (rvec)  ==========
  valre=0.0d0
  r1rcm=r2rcm
  !$omp parallel default(shared)
  call kscg_c_01(ncpx*ncpy*ncpz,neigmx,ncol,l,sveccm(1,1,1,1,ns1,nk),hsvcm,valre)
  !$omp barrier
  !$omp single
  call mpi_allreduce(valre,valall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  valre= valall
  valre=valre*dx*dy*dz
  sval(l,ns1,nk)=valre
  !$omp end single
  do ns2= 1,ncol
    call kscg_c_02(ncpx*ncpy*ncpz,hsvcm(1,1,1,l,ns2),ssvcm(1,1,1,l,ns1+ns2-1,nk),rveccm(1,1,1,ns2),sval(l,ns1,nk))
  end do
  !$omp end parallel
! =====================================================================

! ==========  precondition for CG  ==========
  !$omp parallel default(shared) private(i,ix,iy,iz)
  if (nprecon_cg .eq. 1) then
    !$omp do
    do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)*ncol
      vcm(-nf+i,-nf+1,-nf+1,1)=dcmplx(0.0d0,0.0d0)
    end do
    do ns2= 1,ncol
      !$omp do
      do i=1,ncpx*ncpy*ncpz
        avcm(i,1,1,ns2)=dcmplx(0.0d0,0.0d0)
      end do
      !$omp end do nowait
      !$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        vcm(ix,iy,iz,ns2)=rveccm(ix,iy,iz,ns2)
      end do
      end do
      end do
      !$omp end do nowait
    end do ! ns2
    !$omp barrier
    !$omp single
    call overlap_finitedifference_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,1,ncol,vcm)
    call overlap_fdcheck_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,1,ncol,vcm)
    !$omp end single
    call ksprecondition_c(ncpx,ncpy,ncpz,nf,ncol,vcm,rvcm)
  else
    !$omp do
    do i=1,ncpx*ncpy*ncpz*ncol
      rvcm(i,1,1,1)=rveccm(i,1,1,1)
    end do
    !$omp end do nowait
  end if
  !$omp end parallel
! ===========================================

! ==========  orthonormalize preconditioned residual vector  ==========
  do loop=1,northo
    !$omp parallel default(shared) private(i)
    !$omp do
    do i=1,neigmx
      avcmplx(i)=dcmplx(0.0d0,0.0d0)
    end do
    !$omp end do nowait
    !$omp barrier
    call orthogonalization_c_04(ncpx*ncpy*ncpz,1,neigmx,ncol,1,rvcm,1,l,ssvcm(1,1,1,1,ns1,nk),avcmplx)
    !$omp barrier
    !$omp single
    call mpi_allreduce(avcmplx,avcmplxall,l,mpi_double_complex,mpi_sum,mpicom_space,mpij)
    !$omp end single
    !$omp do
    do i=1,l
      avcmplx(i)=avcmplxall(i)*dx*dy*dz
    end do
    !$omp end parallel
    if (l > nblas03) then
      do ns2=1,ncol
        call zgemv('n',ncpx*ncpy*ncpz,l,-calph,sveccm(1,1,1,1,ns1+ns2-1,nk),ncpx*ncpy*ncpz &
            ,avcmplx,1,calph,rvcm(1,1,1,ns2),1)
      end do
    else
      !$omp parallel default(shared)
      call orthogonalization_c_03(ncpx*ncpy*ncpz,neigmx,l,ncol,rvcm,sveccm(1,1,1,1,ns1,nk),avcmplx)
      !$omp end parallel
    end if
  end do
  !$omp parallel default(shared) private(i)
  !$omp do
  do i=1,nprjmx*num_ppcell*ncol
    cpsep(i,1,1)=dcmplx(0.0d0,0.0d0)
  end do
  call nonlocaloperation_c_01(natom,nprjmx,num_spe,num_ppcell,num_list,1,ncol,1, & ! <
                              ncpx,ncpy,ncpz,                                    & ! <
                              key_natpri_in,key_natpri_inps,                     & ! <
                              dx,dy,dz,skpx,skpy,skpz,                           & ! <
                              nprj,indspe,natinf,natpri,naps,lstvec2,            & ! <
                              lstx,lsty,lstz,natx,naty,natz,                     & ! <
                              vnlocp,rvcm,                                       & ! <
                              cpsep,                                             & ! X
                              vcccm)                                               ! W
  !$omp barrier
  do ns2= 1,ncol
    !$omp single
    do i=1,natom
      na=latom(i)
      if (natpri(na) /= key_natpri_out) then
        iaps=naps(na)
        call mpi_allreduce(cpsep(1,iaps,ns2),crpsep(1,iaps,ns2),nprj(indspe(na)),mpi_double_complex,mpi_sum,mpicom_atom(na),mpij)
      end if
    end do
    !$omp end single
  end do ! ns2
  !$omp do
  do i=1,ncpx*ncpy*ncpz*ncol
    asvcm(i,1,1,1)=rvcm(i,1,1,1)
  end do
  call nonlocaloperation_c_02(0,0,natom,num_spe,1,nums,nprjmx,ncol,1,1,num_list,num_ppcell, & ! <
                              ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                          & ! <
                              key_natpri_in,key_natpri_inps,key_soc_calc,                   & ! <
                              dx,dy,dz,skpx,skpy,skpz,                                      & ! <
                              nprj,indspe,natinf,lstvec2,natpri,naps,natsoc(1),             & ! <
                              lstx,lsty,lstz,natx,naty,natz,                                & ! <
                              vnlocp,crpsep,sss,dijsoc,                                     & ! <
                              asvcm,                                                        & ! X
                              vcccm,avc)                                                      ! W
  !$omp barrier
  !$omp single
  workcc(1)=dcmplx(0.0d0,0.0d0)
  !$omp end single
  do ns2= 1,ncol
    call orthogonalization_c_04(ncpx*ncpy*ncpz,1,1,1,1,rveccm(1,1,1,ns2),1,1,rvcm(1,1,1,ns2),workcc(1))     ! ono10-06 new
  end do
  !$omp end parallel
  call mpi_allreduce(workcc(1),r2rcmall,1,mpi_double_complex,mpi_sum,mpicom_space,mpij)
  r2rcm=r2rcmall*dx*dy*dz
  residual_states(l)=abs(r2rcm)
! =====================================================================

! ==========  compute conjugate vector (pvec)  ==========
  !$omp parallel default(shared) private(i)
  if (nkscg .eq. 1) then
    call kscg_c_03(ncpx*ncpy*ncpz,ncol,pveccm,rvcm,pvcm,spvcm,asvcm,r1rcm,r2rcm)
    !$omp barrier
    !$omp do
    do i=1,nprjmx*num_ppcell*ncol
      cppsep(i,1,1)=crpsep(i,1,1)+r2rcm/r1rcm*cppsep(i,1,1)
    end do
    !$omp end do nowait
  else
    !$omp do
    do i=1,ncpx*ncpy*ncpz*ncol
      pveccm(i,1,1,1)=rvcm(i,1,1,1)
      spvcm(i,1,1,1)=asvcm(i,1,1,1)
    end do
    !$omp end do nowait
    !$omp do
    do i=1,nprjmx*num_ppcell*ncol
      cppsep(i,1,1)=crpsep(i,1,1)
    end do
    !$omp end do nowait
  end if
  !$omp end parallel
! =======================================================

! ==========  orthonormalize conjugate vector (pvec)  ==========
  !$omp parallel default(shared) private(i)
  !$omp do
  do i=1,nprjmx*num_ppcell*ncol
    cpsep(i,1,1)=cppsep(i,1,1)
  end do
  !$omp end do nowait
  !$omp do
  do i=1,ncpx*ncpy*ncpz*ncol
    asvcm(i,1,1,1)=spvcm(i,1,1,1)
  end do
  !$omp end parallel
  do loop=1,northo
    !$omp parallel default(shared) private(i)
    !$omp do
    do i=1,neigmx
      avcmplx(i)=dcmplx(0.0d0,0.0d0)
    end do
    !$omp end parallel
    if (l > nblas04) then
      do ns2=1,ncol
        call zgemv('c',ncpx*ncpy*ncpz,l,calph,ssvcm(1,1,1,1,ns1+ns2-1,nk),ncpx*ncpy*ncpz &
                ,pveccm(1,1,1,ns1+ns2-1),1,cbeta,avcmplx,1)
      end do
    else
      !$omp parallel default(shared)
      call orthogonalization_c_04(ncpx*ncpy*ncpz,1,neigmx,ncol,1,pveccm,1,l,ssvcm(1,1,1,1,ns1,nk),avcmplx)
      !$omp end parallel
    end if
    call mpi_allreduce(avcmplx,avcmplxall,l,mpi_double_complex,mpi_sum,mpicom_space,mpij)
    !$omp parallel default(shared) private(i)
    !$omp do
    do i=1,l
      avcmplx(i)=avcmplxall(i)*dx*dy*dz
    end do
    !$omp end parallel
    if (l > nblas03) then
      do ns2= 1,ncol
        call zgemv('n',ncpx*ncpy*ncpz,l,-calph,sveccm(1,1,1,1,ns1+ns2-1,nk),ncpx*ncpy*ncpz &
            ,avcmplx,1,calph,pveccm(1,1,1,ns2),1)
        call zgemv('n',ncpx*ncpy*ncpz,l,-calph,ssvcm(1,1,1,1,ns1+ns2-1,nk),ncpx*ncpy*ncpz &
            ,avcmplx,1,calph,asvcm(1,1,1,ns2),1)
        call zgemv('n',nprjmx*num_ppcell,l,-calph,cspsep(1,1,1,ns1+ns2-1,nk),nprjmx*num_ppcell &
            ,avcmplx,1,calph,cpsep(1,1,ns2),1)
      end do
    else
      !$omp parallel default(shared) private(i)
      call orthogonalization_c_03(ncpx*ncpy*ncpz,neigmx,l,ncol,pveccm,sveccm(1,1,1,1,ns1,nk),avcmplx)
      call orthogonalization_c_03(ncpx*ncpy*ncpz,neigmx,l,ncol,asvcm,ssvcm(1,1,1,1,ns1,nk),avcmplx)
      do ns2= 1,ncol
        call multc(nprjmx*num_ppcell,1,l,cpsep(1,1,ns2),cspsep(1,1,1,ns1+ns2-1,nk),avcmplx)
      end do
      !$omp end parallel
    end if
  end do ! loop=1,northo
  vecnor=0.0d0
  !$omp parallel default(shared) private(i)
  call orthogonalization_c_05(ncpx*ncpy*ncpz,1,1,ncol,1,pveccm,1,asvcm,vecnor)
  !$omp barrier
  !$omp single
  call mpi_allreduce(vecnor,tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  vecnor= tmpall
  vecnor=vecnor*dx*dy*dz
  if (vecnor<=0.0d0) call stopp('kscg_c: pvec has negative norm')
  tmp=1.0d0/dsqrt(vecnor)
  !$omp end single
  call orthogonalization_c_06(ncpx*ncpy*ncpz*ncol,pveccm,tmp)
  call orthogonalization_c_06(ncpx*ncpy*ncpz*ncol,asvcm,tmp)
  !$omp do
  do i=1,nprjmx*num_ppcell*ncol
    cpsep(i,1,1)=cpsep(i,1,1)*tmp
  end do
  !$omp end do nowait
  !$omp end parallel
! ==============================================================

! ==========  update eigenvector (svec)  ==========
! --  cpsep is computed by the preceding routine
  !$omp parallel default(shared) private(i,ix,iy,iz)
  !$omp do
  do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)*ncol
    vcm(-nf+i,-nf+1,-nf+1,1)=dcmplx(0.0d0,0.0d0)
  end do
  do ns2= 1,ncol
    !$omp do
    do i=1,ncpx*ncpy*ncpz
      avcm(i,1,1,ns2)=dcmplx(0.0d0,0.0d0)
    end do
    !$omp end do nowait
    !$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       vcm(ix,iy,iz,ns2)=pveccm(ix,iy,iz,ns2)
    end do
    end do
    end do
  end do ! ns2
  !$omp single
  call overlap_finitedifference_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,ncol,vcm)
  !$omp end single
  call nonlocaloperation_c_02(1,nso,natom,num_spe,1,nums,nprjmx,ncol,ns1,1,num_list,num_ppcell, & ! <
                              ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                              & ! <
                              key_natpri_in,key_natpri_inps,key_soc_calc,                       & ! <
                              dx,dy,dz,skpx,skpy,skpz,                                          & ! <
                              nprj,indspe,natinf,lstvec2,natpri,naps,natsoc,                    & ! <
                              lstx,lsty,lstz,natx,naty,natz,                                    & ! <
                              vnlocp,cpsep,dij,dijsoc,                                          & ! <
                              avcm,                                                             & ! X
                              vcccm,avc)                                                          ! W
  !$omp barrier
  !$omp single
  call overlap_fdcheck_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,ncol,vcm)
  !$omp end single
  call kslaplacian_c(ncpx,ncpy,ncpz,1,ncol,nvef,nf,dx,dy,dz,xmax,ymax,zmax,skpx,skpy,skpz &
                     ,vcm,veff(1,1,1,nsv),1,avcm,workc)
  !$omp barrier
  !$omp single
  pap=0.0d0
  sapcm=dcmplx(0.0d0,0.0d0)
  !$omp end single
  call kscg_c_04(ncpx*ncpy*ncpz,neigmx,ncol,pveccm,l,sveccm(1,1,1,1,ns1,nk),avcm,pap,sapcm)
  !$omp barrier
  !$omp single
  call mpi_allreduce(pap,papall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  call mpi_allreduce(sapcm,sapcmall,1,mpi_double_complex,mpi_sum,mpicom_space,mpij)
  pap  =papall*dx*dy*dz
  sapcm=sapcmall*dx*dy*dz
  !$omp end single
  call kscg_c_05(ncpx*ncpy*ncpz,nprjmx*num_ppcell,neigmx,ncol,pap,sapcm,sval(l,ns1,nk) &
                ,pveccm,avcm,asvcm,cpsep &
                ,l,hsvcm,sveccm(1,1,1,1,ns1,nk) &
                ,ssvcm(1,1,1,1,ns1,nk) &
                ,cspsep(1,1,1,ns1,nk),l_converged)
  !$omp end parallel
! =================================================

  l_converged= (dsqrt(residual_states(l)) < epssd)
  ksconv= ksconv .and. (l_converged .or. (iiter < nsdmax))

  end do ! iiter

  if (myrank_glbl==0) write (ndisp,'(i4,1x,i4,1x,i4,1x,2e20.10,1x,i5)',err=9999) &
    nk,ns1,l,sval(l,ns1,nk),dsqrt(residual_states(l)),iiter
  ksitmax= max(iiter,ksitmax)
  end do ! l
! ***********************************************

  call kslaplacian_finalize
  call overlap_finitedifference_final

  deallocate(vcm)
  deallocate(avcm)
  deallocate(asvcm)
  deallocate(workc)
  deallocate(avcmplx,avcmplxall)
  deallocate(crpsep,cppsep)
  deallocate(cpsep)
  deallocate(rveccm,pveccm,pvcm,rvcm,spvcm)
  deallocate(vcccm,avc)
  deallocate(workcc)

  return
9999 continue
  call mpi_abort(mpi_comm_world)
end subroutine kscg_c


subroutine kscg_r_01(n,a,b,t)
implicit none
real*8 a(n),b(n)
real*8 t,t_loc
integer n,i
  t_loc=0.0d0
!$omp do
  do i=1,n
     t_loc=t_loc+a(i)*b(i)
  end do
!$omp end do nowait
!$omp critical
  t=t+t_loc
!$omp end critical
  return
end subroutine kscg_r_01


subroutine kscg_c_01(n,neigmx,ncol,l,a,b,t)
implicit none
integer,   intent(in)    :: n,neigmx,ncol,l
complex*16,intent(in)    :: a(n,neigmx,ncol),b(n,neigmx,ncol)
real*8,    intent(inout) :: t
real*8 t_loc
integer i,ns
  t_loc=0.0d0
  do ns= 1,ncol
!$omp do
    do i=1,n
      t_loc=t_loc+dreal(a(i,l,ns))*dreal(b(i,l,ns))+dimag(a(i,l,ns))*dimag(b(i,l,ns))
    end do
!$omp end do nowait
  end do
!$omp critical
  t=t+t_loc
!$omp end critical
  return
end subroutine kscg_c_01


subroutine kscg_r_02(n,a,b,c,d)
implicit none
integer, intent(in)  :: n
real*8,  intent(in)  :: a(n),b(n),d
real*8,  intent(out) :: c(n)
integer i
!$omp do
  do i=1,n
     c(i)=-(a(i)-d*b(i))
  end do
!$omp end do nowait
  return
end subroutine kscg_r_02


subroutine kscg_c_02(n,a,b,c,d)
implicit none
integer, intent(in)    :: n
real*8,intent(in)      :: d
complex*16,intent(in)  :: a(n),b(n)
complex*16,intent(out) :: c(n)
integer i
!$omp do
  do i=1,n
     c(i)=-(a(i)-d*b(i))
  end do
!$omp end do nowait
  return
end subroutine kscg_c_02


subroutine kscg_r_03(n,pvecre,rvre,pvre,spvre,asvre,r1rr,r2rr)
implicit none
real*8 r1rr,r2rr,tmp
integer i,n
real*8 pvecre(n),rvre(n),pvre(n),spvre(n),asvre(n)
  tmp=r2rr/r1rr
!ocl norecurrence(pvecre,pvre,spvre,asvre)
!$omp do
  do i=1,n
     pvre(i)=rvre(i)+tmp*pvre(i)
     pvecre(i)=pvre(i)
     spvre(i)=asvre(i)+tmp*spvre(i)
  end do
!$omp end do nowait
  return
end subroutine kscg_r_03


subroutine kscg_c_03(n,ncol,pveccm,rvcm,pvcm,spvcm,asvcm,r1rcm,r2rcm)
implicit none
integer, intent(in)     :: n,ncol
complex*16,intent(in)   :: r1rcm,r2rcm
complex*16,intent(in)   :: rvcm(n,ncol),asvcm(n,ncol)
complex*16,intent(inout):: pvcm(n,ncol),spvcm(n,ncol)
complex*16,intent(out)  :: pveccm(n,ncol)
complex*16 ctmp0
integer i,ns

  ctmp0=dconjg(r1rcm)*r2rcm/(abs(r1rcm)**2)
  do ns= 1,ncol
!ocl norecurrence(pveccm,pvcm,spvcm,asvcm)
!$omp do
    do i=1,n
       pveccm(i,ns)=rvcm(i,ns)+ctmp0*pvcm(i,ns)
       pvcm(i,ns)=pveccm(i,ns)
       spvcm(i,ns)=asvcm(i,ns)+ctmp0*spvcm(i,ns)
    end do
!$omp end do nowait
  end do
  return
end subroutine kscg_c_03


subroutine kscg_r_04(n,pvecre,svecre,avre,pap,sapre)
implicit none
real*8 pap,sapre,pappap,sapresap
integer i,n
real*8 pvecre(n),svecre(n),avre(n)
  pappap=0.0d0
  sapresap=0.0d0
!$omp do
  do i=1,n
    pappap=pappap+pvecre(i)*avre(i)
    sapresap=sapresap+svecre(i)*avre(i)
  end do
!$omp end do nowait
!$omp critical
  pap=pap+pappap
  sapre=sapre+sapresap
!$omp end critical
  return
end subroutine kscg_r_04


subroutine kscg_c_04(n,neigmx,ncol,pveccm,l,sveccm,avcm,pap,sapcm)
implicit none
integer, intent(in)     :: n,neigmx,ncol,l
complex*16,intent(in)   :: pveccm(n,ncol),sveccm(n,neigmx,ncol),avcm(n,ncol)
real*8,  intent(inout)  :: pap
complex*16,intent(inout):: sapcm
integer i,ns
real*8 pappap
complex*16 sapcmsap
  pappap=0.0d0
  sapcmsap=dcmplx(0.0d0,0.0d0)
  do ns= 1,ncol
!$omp do
    do i=1,n
       pappap=pappap+dreal(pveccm(i,ns))*dreal(avcm(i,ns))+dimag(pveccm(i,ns))*dimag(avcm(i,ns))
       sapcmsap=sapcmsap+dconjg(sveccm(i,l,ns))*avcm(i,ns)
    end do
!$omp end do nowait
  end do
!$omp critical
  pap=pap+pappap
  sapcm=sapcm+sapcmsap
!$omp end critical
  return
end subroutine kscg_c_04


subroutine kscg_r_05(n,n1,pvecre,hsvre,avre,asvre,svecre,ssvre,rpsep,rspsep,pap,sapre,sval,l_converged)
implicit none
real*8 pap,sapre,alph,beta,gamma,sval,eps
integer n,n1,i
real*8 pvecre(n),hsvre(n),avre(n),asvre(n),svecre(n),ssvre(n)
real*8 rpsep(n1),rspsep(n1)
logical l_converged
  eps=1.0d-14
  alph=1.0d0
  beta=0.0d0
  gamma=(sval+pap)*0.5d0-dsqrt((sval-pap)*(sval-pap)*0.25d0+sapre*sapre)
  if (dabs((sval-gamma)/gamma) .gt. eps) then
     alph=sapre/dsqrt(sapre*sapre+(sval-gamma)*(sval-gamma))
     beta=-(sval-gamma)/dsqrt(sapre*sapre+(sval-gamma)*(sval-gamma))
  end if
  if (alph<eps) then
    alph=-alph
    beta=-beta
    l_converged =.true.
  end if
!cdir nodep
!$omp do
  do i=1,n
     svecre(i)=alph*svecre(i)+beta*pvecre(i)
     hsvre(i)=alph*hsvre(i)+beta*avre(i)
     ssvre(i)=alph*ssvre(i)+beta*asvre(i)
  end do
!cdir nodep
!$omp do
  do i=1,n1
    rspsep(i)=alph*rspsep(i)+beta*rpsep(i)
  end do
  return
end subroutine kscg_r_05


subroutine kscg_c_05(n,n1,neigmx,ncol,pap,sapcm,sval,pveccm,avcm,asvcm,cpsep  &
                    ,l,hsvcm,sveccm,ssvcm,cspsep,l_converged)
implicit none
integer,    intent(in)    :: n,n1,neigmx,ncol, l
real*8,     intent(in)    :: pap,sval
complex*16, intent(in)    :: sapcm
complex*16, intent(in)    :: pveccm(n,ncol),avcm(n,ncol),asvcm(n,ncol)
complex*16, intent(in)    :: cpsep(n1,ncol)
complex*16, intent(inout) :: hsvcm(n,neigmx,ncol),sveccm(n,neigmx,ncol),ssvcm(n,neigmx,ncol)
complex*16, intent(inout) :: cspsep(n1,neigmx,ncol)
logical,    intent(inout) :: l_converged
integer i,ns2
complex*16 calph,cbeta
real*8 alph,gamma,eps
  eps=1.0d-14
  gamma=(sval+pap)*0.5d0-dsqrt((sval-pap)*(sval-pap)*0.25d0+dreal(dconjg(sapcm)*sapcm))
  if (dabs((sval-gamma)/gamma) .gt. eps) then
    calph=sapcm/dsqrt(dreal(dconjg(sapcm)*sapcm)+(sval-gamma)*(sval-gamma))
    cbeta=-(sval-gamma)/dsqrt(dreal(dconjg(sapcm)*sapcm)+(sval-gamma)*(sval-gamma))
    alph=abs(calph)
    if (alph>eps) cbeta=cbeta*alph/calph
  else
    alph  = 1.0d0
    cbeta = dcmplx(0.0d0,0.0d0)
    l_converged = .true.
  end if
  do ns2=1,ncol
!cdir nodep
!$omp do
    do i=1,n
      sveccm(i,l,ns2)=alph*sveccm(i,l,ns2)+cbeta*pveccm(i,ns2)
      hsvcm(i,l,ns2) =alph*hsvcm(i,l,ns2)+cbeta*avcm(i,ns2)
      ssvcm(i,l,ns2) =alph*ssvcm(i,l,ns2)+cbeta*asvcm(i,ns2)
    end do
!cdir nodep
!$omp do
    do i=1,n1
      cspsep(i,l,ns2)=alph*cspsep(i,l,ns2)+cbeta*cpsep(i,ns2)
    end do
  end do ! ns2
  return
end subroutine kscg_c_05


end module
