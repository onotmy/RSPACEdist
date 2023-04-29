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
! **********  ksdiis8f.F90 06/20/2018-01  **********

module mod_ksdiis
implicit none

contains

subroutine ksdiis_r(l_hsv,natom,num_spe,neigmx,nums,nspv,numk,nperi,ndiismax,     & ! <
                   nprecon_diis,nk,ns,nf,nprjmx,num_list,num_ppcell,ndisp,northo, & ! <
                   ncpx,ncpy,ncpz,                                                & ! <
                   key_ortho_cmpt_innerproduct,                                   & ! <
                   key_natpri_in,key_natpri_inps,                                 & ! <
                   eps_eig_diis,alambda_diis,                                     & ! <
                   ratio_diis,alambda_min,alambda_max,epssd,                      & ! <
                   xmax,ymax,zmax,                                                & ! <
                   indspe,natpri,naps,nprj,natinf,lstvec2,latom,ntyppp,           & ! <
                   dij,veff,vnlocp,sss,                                           & ! <
                   ksconv,ksitmax,                                                & ! X
                   rspsep,sval,svecre,ssvre,                                      & ! X
                   residual_states,                                               & ! X
                   hsvre)                                                           ! >
use mod_mpi
use mod_overlap_finitedifference, only:overlap_finitedifference_init,overlap_finitedifference_final, &
                                       overlap_finitedifference_r,overlap_fdcheck_r
use mod_nonlocaloperation, only:nonlocaloperation_r_01,nonlocaloperation_r_02
use mod_orthogonalization, only:orthogonalization_r_01,orthogonalization_r_04,orthogonalization_r_05, &
                                orthogonalization_r_06
use mod_ksprecondition, only:ksprecondition
use mod_kslaplacian, only:kslaplacian_r,kslaplacian_initialize,kslaplacian_finalize
implicit none
integer,parameter    ::nitermax=2
logical,intent(in)   ::l_hsv   
integer,intent(in)   ::natom,num_spe,neigmx,nums,nspv,numk,nperi,ndiismax
integer,intent(in)   ::nprecon_diis,nk,ns,nf,nprjmx,num_list,num_ppcell,ndisp,northo
integer,intent(in)   ::ncpx,ncpy,ncpz
integer,intent(in)   ::key_ortho_cmpt_innerproduct
integer,intent(in)   ::key_natpri_in,key_natpri_inps
real*8, intent(in)   ::eps_eig_diis,alambda_diis,ratio_diis,alambda_min,alambda_max,epssd
real*8, intent(in)   ::xmax,ymax,zmax
integer,intent(in)   ::indspe(natom),natpri(natom),naps(natom),nprj(num_spe)
integer,intent(in)   ::natinf(natom),lstvec2(num_list,num_ppcell),latom(natom),ntyppp(num_spe)
real*8, intent(in)   ::dij(nprjmx,nprjmx,nums,natom),veff(ncpx,ncpy,ncpz,nspv)
real*8, intent(in)   ::vnlocp(num_list,nprjmx,num_ppcell),sss(nprjmx,nprjmx,num_spe)
logical,intent(inout)::ksconv
integer,intent(inout)::ksitmax
real*8, intent(inout)::rspsep(nprjmx,num_ppcell,neigmx,nums,numk)
real*8, intent(inout)::sval(neigmx,nums,numk)
real*8, intent(inout)::svecre(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8, intent(inout)::ssvre(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8, intent(inout)::residual_states(neigmx)
real*8, intent(out)  ::hsvre(ncpx,ncpy,ncpz,neigmx)
real*8  dx,dy,dz
real*8  valre,valall,vecnor,tmp,residual_statesold,residual_statesfirst
logical l_converged
integer na,iiter,jiter,i,j,k,l,ix,iy,iz,ierr,iaps
integer nmat
real*8  aaa,bbb,alambda,alambda1
real*8,allocatable::rvecre(:,:,:,:),rvre(:,:,:)
real*8,allocatable::phire(:,:,:,:),hphire(:,:,:,:),sphire(:,:,:,:)
real*8,allocatable::vre(:,:,:)
real*8,allocatable::avre(:,:,:)
real*8,allocatable::asvre(:,:,:)
real*8,allocatable::tpsplr(:,:)
real*8,allocatable::rpsep(:,:),rphisep(:,:,:)     ! for isolated system
real*8,allocatable::eigen(:),work(:),work2(:)
real*8,allocatable::rrmat(:,:),rsmat(:,:),rsmat0(:,:),rvect(:,:)
real*8,allocatable::rrmat2(:,:),rvect2(:,:)
real*8,allocatable::sumtmp(:),sumtmpall(:)
real*8,allocatable::avr(:)
real*8,allocatable::rtmp(:),rtmpall(:)

  dx=2.0d0*xmax/(ncpx*nprocx)
  dy=2.0d0*ymax/(ncpy*nprocy)
  dz=2.0d0*zmax/(ncpz*nprocz)

  if (nspv/=nums) then
    write(ndisp,*) 'WARNING!! Real diis routine is skipped for nspv/=nums !'
    return
  endif

  allocate(vre(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf))
  allocate(avre(ncpx,ncpy,ncpz))
  allocate(asvre(ncpx,ncpy,ncpz))
  allocate(tpsplr(nprjmx,num_ppcell))
  allocate(rpsep(nprjmx,num_ppcell))
  allocate(rphisep(nprjmx,num_ppcell,0:nitermax))
  allocate(rvecre(ncpx,ncpy,ncpz,0:nitermax),rvre(ncpx,ncpy,ncpz))
  allocate(phire(ncpx,ncpy,ncpz,0:nitermax),sphire(ncpx,ncpy,ncpz,0:nitermax),hphire(ncpx,ncpy,ncpz,0:nitermax))
  allocate(sumtmp(4),sumtmpall(4))
  allocate(avr(num_list))
  allocate(rtmp(nitermax*(nitermax+1)),rtmpall(nitermax*(nitermax+1)))

  call kslaplacian_initialize(nf)
  call overlap_finitedifference_init(ncpx,ncpy,ncpz,nf,1,0)

  jiter=0
  do while (jiter+1 < ndiismax)
  jiter=jiter+1

  if (myrank_glbl .eq. 0) then
    write(ndisp,*,err=9999) 'DIIS iter. #',jiter
    write(ndisp,*,err=9999) '  k ,spin,band,     eigen value    ,   residual norm   ,matrix dim.'
  end if

  do l=1,neigmx
!$omp parallel default(shared) private(i)
!$omp do
  do i=1,ncpx*ncpy*ncpz
    phire(i,1,1,0)=svecre(i,1,1,l,ns,nk)
    sphire(i,1,1,0)=ssvre(i,1,1,l,ns,nk)
  end do
!$omp do
  do i=1,nprjmx*num_ppcell
    rphisep(i,1,0)=rspsep(i,1,l,ns,nk)
  end do
!$omp end parallel

  iiter=0
  l_converged= .false.
  do while ((iiter+1 < nitermax).and.(.not. l_converged))
  iiter=iiter+1
! ==========  compute eigenvalue and residual vector (rvec)  ==========
!$omp parallel default(shared) private(i,ix,iy,iz)
!$omp do
  do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)
    vre(-nf+i,-nf+1,-nf+1)=0.0d0
  end do
!$omp do
  do i=1,ncpx*ncpy*ncpz
    avre(i,1,1)=0.0d0
  end do
!$omp do
  do iz=1,ncpz
  do iy=1,ncpy
  do ix=1,ncpx
    vre(ix,iy,iz)=phire(ix,iy,iz,iiter-1)
  end do
  end do
  end do
!$omp end parallel

  if (iiter .eq. 1) then
!$omp parallel default(shared) private(i,ix,iy,iz)
!$omp do
    do i=1,nprjmx*num_ppcell
      rpsep(i,1)=rphisep(i,1,iiter-1)
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
    call kslaplacian_r(ncpx,ncpy,ncpz,1,nf,dx,dy,dz,vre,veff(1,1,1,ns),1,avre)
!$omp barrier
!$omp do
    do i=1,ncpx*ncpy*ncpz
      hphire(i,1,1,iiter-1)=avre(i,1,1)
    end do
!$omp end parallel
  end if
  valre=0.0d0
!$omp parallel default(shared)
  call orthogonalization_r_05(ncpx*ncpy*ncpz,phire(1,1,1,iiter-1),hphire(1,1,1,iiter-1),valre)
!$omp barrier
!$omp single
  call mpi_allreduce(valre,valall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  valre=valall*dx*dy*dz
  sval(l,ns,nk)=valre
!$omp end single
  call ksdiis_r_01( &
   ncpx*ncpy*ncpz,-valre,hphire(1,1,1,iiter-1),sphire(1,1,1,iiter-1), & ! <
   rvecre(1,1,1,iiter-1))                                                ! >
!$omp end parallel
! =====================================================================

! ==========  precondition ==========
!$omp parallel default(shared) private(i,ix,iy,iz)
  if (nprecon_diis .eq. 1) then
!$omp do
    do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)
      vre(-nf+i,-nf+1,-nf+1)=0.0d0
    end do
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
      vre(ix,iy,iz)=rvecre(ix,iy,iz,iiter-1)
    end do
    end do
    end do
!$omp single
    call overlap_finitedifference_r(nperi,ncpx,ncpy,ncpz,nf,nf-1,1,vre)
    call overlap_fdcheck_r(nperi,ncpx,ncpy,ncpz,nf,nf-1,1,vre)
!$omp end single
    call ksprecondition(ncpx,ncpy,ncpz,nf,1,vre,rvre)
  else
!$omp do
    do i=1,ncpx*ncpy*ncpz
      rvre(i,1,1)=rvecre(i,1,1,iiter-1)
    end do
  end if
!$omp end parallel
! ===================================

!$omp parallel default(shared) private(i,ix,iy,iz)
!$omp do
  do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)
    vre(-nf+i,-nf+1,-nf+1)=0.0d0
  end do
!$omp do
  do i=1,ncpx*ncpy*ncpz
    avre(i,1,1)=0.0d0
  end do
!$omp do
  do iz=1,ncpz
  do iy=1,ncpy
  do ix=1,ncpx
    vre(ix,iy,iz)=rvre(ix,iy,iz)
  end do
  end do
  end do
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
!$omp do
  do i=1,nprjmx*num_ppcell
    tpsplr(i,1)=rpsep(i,1)
  end do
!$omp single
  do i=1,natom
    na=latom(i)
    if (natpri(na) .ge. key_natpri_inps) then
      iaps=naps(na)
      call mpi_allreduce(tpsplr(1,iaps),rpsep(1,iaps),nprj(indspe(na)),mpi_double_precision,mpi_sum,mpicom_atom(na),mpij)
    end if
  end do
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
  call kslaplacian_r(ncpx,ncpy,ncpz,1,nf,dx,dy,dz,vre,veff(1,1,1,ns),1,avre)
!$omp barrier
!$omp do
  do i=1,ncpx*ncpy*ncpz
    asvre(i,1,1)=rvre(i,1,1)
  end do
  call nonlocaloperation_r_02(0,natom,num_spe,nums,nprjmx,1,num_list,num_ppcell, & ! <
                              ncpx,ncpy,ncpz,                                    & ! <
                              key_natpri_in,key_natpri_inps,                     & ! <
                              dx,dy,dz,                                          & ! <
                              nprj,indspe,natinf,lstvec2,natpri,naps,            & ! <
                              vnlocp,rpsep,sss,                                  & ! <
                              asvre,                                             & ! X
                              avr)                                                 ! W
!$omp end parallel

  if (iiter .eq. 1) then
! ==========  compute \alambda ==========
  alambda1=alambda_diis
  sumtmp(1:4)=0.0d0
!$omp parallel default(shared) private(i,ix,iy,iz)
!$omp do
  do i=1,nprjmx*num_ppcell
    rphisep(i,1,1)=rphisep(i,1,0)+alambda1*rpsep(i,1)
  end do
  call ksdiis_r_01( &
   ncpx*ncpy*ncpz,alambda1,phire(1,1,1,0),rvre,   & ! <
   phire(1,1,1,1))                                  ! >
  call ksdiis_r_01( &
   ncpx*ncpy*ncpz,alambda1,hphire(1,1,1,0),avre,  & ! <
   hphire(1,1,1,1))                                 ! >
  call ksdiis_r_01( &
   ncpx*ncpy*ncpz,alambda1,sphire(1,1,1,0),asvre, & ! <
   sphire(1,1,1,1))                                 ! >
  call orthogonalization_r_05(ncpx*ncpy*ncpz,rvre,hphire(1,1,1,0),sumtmp(1))
  call orthogonalization_r_05(ncpx*ncpy*ncpz,rvre,sphire(1,1,1,0),sumtmp(2))
  call orthogonalization_r_05(ncpx*ncpy*ncpz,phire(1,1,1,1),hphire(1,1,1,1),sumtmp(3))
  call orthogonalization_r_05(ncpx*ncpy*ncpz,phire(1,1,1,1),sphire(1,1,1,1),sumtmp(4))
!$omp end parallel
  call mpi_allreduce(sumtmp,sumtmpall,4,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  aaa=2.0d0*(sumtmpall(1)-sval(l,ns,nk)*sumtmpall(2))*dx*dy*dz
  valre=sumtmpall(3)/sumtmpall(4)
  bbb=(valre-sval(l,ns,nk)-aaa*alambda1)/(alambda1*alambda1)
! I think alambda is in opposite direction to that in PRB54 11169 (1996)
! because the residual vector is in opposite sign to the conjugate vector of CG.
! 04/19/2006 T.Ono.
!  alambda=-0.5d0*aaa/bbb
  alambda=0.5d0*aaa/bbb
  if (alambda .gt. alambda_max) alambda=alambda_max
  if (alambda .lt. alambda_min) alambda=alambda_min
! ======================================
  end if !iiter

! ==========  compute trial phi ==========
!$omp parallel default(shared) private(i,ix,iy,iz)
  call ksdiis_r_01( &
   ncpx*ncpy*ncpz,alambda,phire(1,1,1,iiter-1),rvre,   & ! <
   phire(1,1,1,iiter))                                   ! >
  call ksdiis_r_01( &
   ncpx*ncpy*ncpz,alambda,hphire(1,1,1,iiter-1),avre,  & ! <
   hphire(1,1,1,iiter))                                  ! >
  call ksdiis_r_01( &
   ncpx*ncpy*ncpz,alambda,sphire(1,1,1,iiter-1),asvre, & ! <
   sphire(1,1,1,iiter))                                  ! >
!$omp do
  do i=1,nprjmx*num_ppcell
    rphisep(i,1,iiter)=rphisep(i,1,iiter-1)+alambda*rpsep(i,1)
  end do
!$omp end parallel
! ========================================

! ==========  normalize trial phi  ==========
!  vecnor=0.0d0
!!$omp parallel default(shared) private(i)
!  call orthogonalization_r_05(ncpx*ncpy*ncpz,phire(1,1,1,iiter),sphire(1,1,1,iiter),vecnor)
!!$omp barrier
!!$omp single
!  call mpi_allreduce(vecnor,tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
!  vecnor=tmpall*dx*dy*dz
!  tmp=1.0d0/dsqrt(vecnor)
!!$omp end single
!  call orthogonalization_r_06(ncpx*ncpy*ncpz,phire(1,1,1,iiter),tmp)
!  call orthogonalization_r_06(ncpx*ncpy*ncpz,hphire(1,1,1,iiter),tmp)
!  call orthogonalization_r_06(ncpx*ncpy*ncpz,sphire(1,1,1,iiter),tmp)
!!$omp barrier
!!$omp do
!  do i=1,nprjmx*num_ppcell
!    rphisep(i,1,iiter)=rphisep(i,1,iiter)*tmp
!  end do
!!$omp end parallel
! ===========================================

! ==========  compute eigenvalue and residual vector (rvec)  ==========
!  valre=0.0d0
!!$omp parallel shared(phire,hphire,valre)
!  call orthogonalization_r_05(ncpx*ncpy*ncpz,phire(1,1,1,iiter),hphire(1,1,1,iiter),valre)
!!$omp end parallel
!  call mpi_allreduce(valre,valall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
!  valre=valall*dx*dy*dz
!$omp parallel default(shared) private(j)
  do j=0,iiter
    call ksdiis_r_01( &
     ncpx*ncpy*ncpz,-sval(l,ns,nk),hphire(1,1,1,j),sphire(1,1,1,j), & ! <
     rvecre(1,1,1,j))                                                 ! >
  end do
!$omp end parallel
! =====================================================================

! ==========  diis step  ==========
! ----------  set up matrix  ----------
  allocate(eigen(iiter+1),work((iiter+1)*3))
  allocate(rrmat(iiter+1,iiter+1),rsmat(iiter+1,iiter+1),rsmat0(iiter+1,iiter+1),rvect(iiter+1,iiter+1))


  k=0
  do j=0,iiter
    do i=0,j
      rtmp(k+1:k+2)=dcmplx(0.0d0,0.0d0)
!$omp parallel default(shared)
      call orthogonalization_r_04(ncpx*ncpy*ncpz,1,1,1,sphire(1,1,1,j),1,1,phire(1,1,1,i),rtmp(k+1))
      call orthogonalization_r_04(ncpx*ncpy*ncpz,1,1,1,rvecre(1,1,1,j),1,1,rvecre(1,1,1,i),rtmp(k+2))
!$omp end parallel
      k=k+2
    end do
  end do
  call mpi_allreduce(rtmp,rtmpall,(iiter+1)*(iiter+2),mpi_double_precision,mpi_sum,mpicom_space,mpij)
  k=0
  do j=0,iiter
    do i=0,j
      rsmat(i+1,j+1)=rtmpall(k+1)*dx*dy*dz
      rsmat(j+1,i+1)=rtmpall(k+1)*dx*dy*dz
      rrmat(i+1,j+1)=rtmpall(k+2)*dx*dy*dz
      rrmat(j+1,i+1)=rtmpall(k+2)*dx*dy*dz
      k=k+2
    end do
  end do
  vecnor=rsmat(iiter+1,iiter+1)
  tmp=1.0d0/dsqrt(vecnor)
!$omp parallel default(shared)
  call orthogonalization_r_06(ncpx*ncpy*ncpz,phire(1,1,1,iiter),tmp)
  call orthogonalization_r_06(ncpx*ncpy*ncpz,hphire(1,1,1,iiter),tmp)
  call orthogonalization_r_06(ncpx*ncpy*ncpz,sphire(1,1,1,iiter),tmp)
  call orthogonalization_r_06(ncpx*ncpy*ncpz,rvecre(1,1,1,iiter),tmp)
!$omp do
  do i=1,nprjmx*num_ppcell
    rphisep(i,1,iiter)=rphisep(i,1,iiter)*tmp
  end do
!$omp end parallel
  do i=1,iiter+1
    rsmat(iiter+1,i)=rsmat(iiter+1,i)*tmp
    rrmat(iiter+1,i)=rrmat(iiter+1,i)*tmp
  end do
  do i=1,iiter+1
    rsmat(i,iiter+1)=rsmat(i,iiter+1)*tmp
    rrmat(i,iiter+1)=rrmat(i,iiter+1)*tmp
  end do
  rsmat0=rsmat
! -------------------------------------
  call dsyev('v','u',iiter+1,rsmat,iiter+1,eigen,work,(iiter+1)*3,ierr)
  do i=1,iiter+1
    work(i)=eigen(i)
  end do
  do j=1,iiter+1
    do i=1,iiter+1
      rvect(i,j)=rsmat(i,iiter+2-j)
    end do
    eigen(j)=work(iiter+2-j)
  end do
! ----------  eliminate non positive eigenvalues  ----------
  nmat=0
  k=0
  do while (k < iiter+1)
    nmat=nmat+1
    k=k+1
    if (eigen(nmat) .lt. eps_eig_diis) then
      do i=nmat,iiter-1
        rvect(:,i)=rvect(:,i+1)
        eigen(i)=eigen(i+1)
      end do
      nmat=nmat-1
!      write(6,*) 'rsmat is not positive definite!!!!'
    end if
  end do
! ----------------------------------------------------------
  if (nmat .gt. 1) then
! ----------  set up new matrix  ----------
    do i=1,nmat
      rvect(:,i)=rvect(:,i)/dsqrt(eigen(i))
    end do
    allocate(rrmat2(nmat,nmat),rvect2(nmat,nmat),work2(nmat*3))
    rsmat=0.0d0
    do k=1,iiter+1
      do j=1,nmat
        do i=1,iiter+1
          rsmat(k,j)=rsmat(k,j)+rrmat(k,i)*rvect(i,j)
        end do
      end do
    end do
    rrmat2=0.0d0
    do k=1,nmat
      do j=1,nmat
        do i=1,iiter+1
          rrmat2(k,j)=rrmat2(k,j)+rvect(i,k)*rsmat(i,j)
        end do
      end do
    end do
! -----------------------------------------
    ierr=0
    call dsyev('v','u',nmat,rrmat2,nmat,eigen,work2,nmat*3,ierr)
    do i=1,nmat
      work2(i)=eigen(i)
    end do
    do j=1,nmat
      do i=1,nmat
        rvect2(i,j)=rrmat2(i,nmat+1-j)
      end do
      eigen(j)=work2(nmat+1-j)
    end do
    if ((ierr/=0) .and. (myrx**2+myry**2+myrz**2==0)) &
       write(ndisp,*) 'error occurs in dsyev of diis!!','myr_kpt=',myr_kpt,'nk=',nk,'ns=',ns,'l=',l
!    write(ndisp,*) 'eigen',(eigen(i),i=1,iiter+1)
    if (iiter .gt. 1) residual_statesold=residual_states(l)
    if (iiter .eq. 2) residual_statesfirst=residual_states(l)/residual_statesold
    residual_states(l)=eigen(nmat)
    if (iiter .gt. 2) then
!   Be careful!! If iiter=1, residual_statesfirst=0.
      if (residual_states(l)/residual_statesold/residual_statesfirst .lt. ratio_diis) then
        iiter=iiter-1
        nmat=nmat-1
        residual_states(l)=residual_statesold
        deallocate(rrmat2,rvect2,work2)
        l_converged= .true.
        goto 1010
      end if
    end if
! ----------  replicate eigenvectors  ----------
    rrmat=0.0d0
    do k=1,nmat
      do j=1,iiter+1
        do i=1,nmat
          rrmat(j,k)=rrmat(j,k)+rvect(j,i)*rvect2(i,k)
        end do
      end do
    end do
    tmp=0.0d0
    do j=1,iiter+1
      do i=1,iiter+1
        tmp=tmp+rrmat(i,nmat)*rsmat0(i,j)*rrmat(j,nmat)
      end do
    end do
    rrmat(:,nmat)=rrmat(:,nmat)/dsqrt(tmp)
    deallocate(rrmat2,rvect2,work2)
! ----------------------------------------------
  else
    rrmat(:,nmat)=0.0d0
    rrmat(iiter+1,nmat)=1.0d0
  end if
! ----------  compute new wave function  ----------
!$omp parallel default(shared) private(i,j)
!$omp do
  do i=1,ncpx*ncpy*ncpz
    avre(i,1,1)=0.0d0
  end do
  call ksdiis_r_03( &
   ncpx*ncpy*ncpz,iiter+1,phire(1,1,1,0),rrmat(1,nmat), & ! <
   avre)                                                   ! X
!$omp barrier
!$omp do
  do i=1,ncpx*ncpy*ncpz
    phire(i,1,1,iiter)=avre(i,1,1)
    avre(i,1,1)=0.0d0
  end do
  call ksdiis_r_03( &
   ncpx*ncpy*ncpz,iiter+1,sphire(1,1,1,0),rrmat(1,nmat), & ! <
   avre)                                                    ! X
!$omp barrier
!$omp do
  do i=1,ncpx*ncpy*ncpz
    sphire(i,1,1,iiter)=avre(i,1,1)
    avre(i,1,1)=0.0d0
  end do
  call ksdiis_r_03( &
   ncpx*ncpy*ncpz,iiter+1,hphire(1,1,1,0),rrmat(1,nmat), & ! <
   avre)                                                    ! X
!$omp barrier
!$omp do
  do i=1,ncpx*ncpy*ncpz
    hphire(i,1,1,iiter)=avre(i,1,1)
  end do
!$omp do
  do i=1,nprjmx*num_ppcell
    rpsep(i,1)=0.0d0
  end do
  do j=1,iiter+1
!$omp do
    do i=1,nprjmx*num_ppcell
      rpsep(i,1)=rpsep(i,1)+rrmat(j,nmat)*rphisep(i,1,j-1)
    end do
  end do
!$omp do
  do i=1,nprjmx*num_ppcell
    rphisep(i,1,iiter)=rpsep(i,1)
  end do
!$omp end parallel
! -------------------------------------------------
! ----------  check for DIIS  ----------
!  avre=0.0d0
!!$omp parallel default(shared)
!  call ksdiis_r_03( &
!   ncpx*ncpy*ncpz,iiter+1,rvecre(1,1,1,0),rrmat(1,nmat), & ! <
!   avre)                                                    ! X
!!$omp end parallel
!  rvecre(:,:,:,iiter)=avre
!  vecnor=0.0d0
!!$omp parallel shared(rvecre,rvecre,vecnor)
!  call orthogonalization_r_05(ncpx*ncpy*ncpz,rvecre(1,1,1,iiter),rvecre(1,1,1,iiter),vecnor)
!!$omp end parallel
!  call mpi_allreduce(vecnor,tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
!  vecnor=tmpall*dx*dy*dz
!  write(ndisp,*) 'norm of residual vector=',vecnor
! --------------------------------------

  1010 continue

  deallocate(eigen,work)
  deallocate(rrmat,rsmat,rsmat0,rvect)
! ======================================

!  if (myrank_glbl .eq. 0) then
!  write (ndisp,9002,err=9999) nk,ns,l,sval(l,ns,nk),dsqrt(residual_states(l)),nmat
!  end if

  end do ! iiter

!$omp parallel default(shared) private(i,ix,iy,iz)
!$omp do
  do i=1,ncpx*ncpy*ncpz
    svecre(i,1,1,l,ns,nk)=phire(i,1,1,iiter)
    ssvre(i,1,1,l,ns,nk)=sphire(i,1,1,iiter)
    hsvre(i,1,1,l)=hphire(i,1,1,iiter)
  end do
!$omp do
  do i=1,nprjmx*num_ppcell
    rspsep(i,1,l,ns,nk)=rphisep(i,1,iiter)
  end do
!$omp end parallel

  if (myrank_glbl .eq. 0) then
  write (ndisp,'(i4,1x,i4,1x,i6,1x,2e20.10,1x,i3)',err=9999) &
   nk,ns,l,sval(l,ns,nk),dsqrt(residual_states(l)),nmat
  end if

  ksitmax= max(iiter,ksitmax)
  ksconv = ksconv .and. (dsqrt(residual_states(l)) <= epssd)

  end do ! l
! ***********************************************

  if (l_hsv) then

  if (myrank_glbl .eq. 0) then
    write (ndisp,*) 'orthogonalization was implemented.'
  end if

! **********  orthonormalize wave function (svec)  **********
    do i=1,northo
      call orthogonalization_r_01(natom,neigmx,nprjmx,nums,numk,nk,ns,0,num_list,num_ppcell,num_spe, & ! <
                                  key_ortho_cmpt_innerproduct,ndisp,                                 & ! <
                                  ncpx,ncpy,ncpz,                                                    & ! <
                                  key_ortho_cmpt_innerproduct,                                       & ! <
                                  key_natpri_in,key_natpri_inps,                                     & ! <
                                  dx,dy,dz,                                                          & ! <
                                  nprj,indspe,natpri,naps,natinf,lstvec2,latom,                      & ! <
                                  vnlocp,sss,                                                        & ! <
                                  svecre,ssvre,rspsep)                                                 ! X
    end do
! ***********************************************************

! **********  compute innerproduct between H and \psi  **********
!$omp parallel default(shared) private(i,ix,iy,iz)
   do l=1,neigmx
!$omp do
      do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)
        vre(-nf+i,-nf+1,-nf+1)=0.0d0
      end do
!$omp do
      do i=1,ncpx*ncpy*ncpz
        hsvre(i,1,1,l)=0.0d0
      end do
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        vre(ix,iy,iz)=svecre(ix,iy,iz,l,ns,nk)
      end do
      end do
      end do
!$omp do
      do i=1,nprjmx*num_ppcell
        rpsep(i,1)=rspsep(i,1,l,ns,nk)
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
                                  hsvre(1,1,1,l),                                     & ! X
                                  avr)                                                  ! W
!$omp barrier
!$omp single
      call overlap_fdcheck_r(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,vre)
!$omp end single
      call kslaplacian_r(ncpx,ncpy,ncpz,neigmx,nf,dx,dy,dz,vre,veff(1,1,1,ns),l,hsvre)
!$omp barrier
    end do
!$omp end parallel
! ***************************************************************

  end if

  end do !jiter

  call kslaplacian_finalize
  call overlap_finitedifference_final

  deallocate(vre)
  deallocate(avre)
  deallocate(asvre)
  deallocate(tpsplr)
  deallocate(rpsep)
  deallocate(rphisep)
  deallocate(rvecre,rvre)
  deallocate(phire,sphire,hphire)
  deallocate(sumtmp,sumtmpall)
  deallocate(avr)
  deallocate(rtmp,rtmpall)
  return
9999 continue
  call mpi_abort(mpi_comm_world)
end subroutine


subroutine ksdiis_c(l_hsv,nso,natom,num_spe,neigmx,nums,ncol,nspv,numk,nperi,ndiismax, & ! <
                   nprecon_diis,nk,ns1,nf,nprjmx,num_list,num_ppcell,ndisp,northo,     & ! <
                   ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                                & ! <
                   key_ortho_cmpt_innerproduct,                                        & ! <
                   key_natpri_in,key_natpri_inps,key_soc_calc,                         & ! <
                   eps_eig_diis,alambda_diis,                                          & ! <
                   ratio_diis,alambda_min,alambda_max,epssd,                           & ! <
                   xmax,ymax,zmax,skpx,skpy,skpz,                                      & ! <
                   indspe,natpri,naps,nprj,natinf,lstvec2,latom,natsoc,ntyppp,         & ! <
                   lstx,lsty,lstz,natx,naty,natz,                                      & ! <
                   dij,dijsoc,veff,vnlocp,sss,                                         & ! <
                   ksconv,ksitmax,                                                     & ! X
                   cspsep,sval,sveccm,ssvcm,                                           & ! X
                   residual_states,                                                    & ! X
                   hsvcm)                                                                ! >
use mod_mpi
use mod_overlap_finitedifference, only:overlap_finitedifference_init,overlap_finitedifference_final, &
                                       overlap_finitedifference_c,overlap_fdcheck_c
use mod_nonlocaloperation, only:nonlocaloperation_c_01,nonlocaloperation_c_02
use mod_orthogonalization, only:orthogonalization_c_01,orthogonalization_c_04,orthogonalization_c_05, &
                                orthogonalization_c_06
use mod_ksprecondition, only:ksprecondition_c
use mod_kslaplacian, only:kslaplacian_c,kslaplacian_initialize,kslaplacian_finalize
implicit none
integer,parameter ::nb0=64
integer,parameter    ::nitermax=2
logical,   intent(in)   ::l_hsv   
integer,   intent(in)   ::nso,natom,num_spe,neigmx,nums,ncol,nspv,numk,nperi,ndiismax
integer,   intent(in)   ::nprecon_diis,nk,ns1,nf,nprjmx,num_list,num_ppcell,ndisp,northo
integer,   intent(in)   ::ncpx,ncpy,ncpz,npxmax,npymax,npzmax
integer,   intent(in)   ::key_ortho_cmpt_innerproduct
integer,   intent(in)   ::key_natpri_in,key_natpri_inps,key_soc_calc
real*8,    intent(in)   ::eps_eig_diis,alambda_diis,ratio_diis,alambda_min,alambda_max,epssd
real*8,    intent(in)   ::xmax,ymax,zmax,skpx,skpy,skpz
integer,   intent(in)   ::indspe(natom),natpri(natom),naps(natom),nprj(num_spe)
integer,   intent(in)   ::natinf(natom),lstvec2(num_list,num_ppcell),latom(natom),natsoc(natom),ntyppp(num_spe)
integer,   intent(in)   ::lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell)
integer,   intent(in)   ::natx(natom),naty(natom),natz(natom)
real*8,    intent(in)   ::dij(nprjmx,nprjmx,nums*ncol,natom)
real*8,    intent(in)   ::dijsoc(nprjmx*nso-nso+1,nprjmx*nso-nso+1,3*nso-nso+1,natom*nso-nso+1)
real*8,    intent(in)   ::veff(ncpx,ncpy,ncpz,nspv)
real*8,    intent(in)   ::vnlocp(num_list,nprjmx,num_ppcell),sss(nprjmx,nprjmx,num_spe)
logical,   intent(inout)::ksconv
integer,   intent(inout)::ksitmax
complex*16,intent(inout)::cspsep(nprjmx,num_ppcell,neigmx,nums,numk)
real*8,    intent(inout)::sval(neigmx,nums+1-ncol,numk)
complex*16,intent(inout)::sveccm(ncpx,ncpy,ncpz,neigmx,nums,numk)
complex*16,intent(inout)::ssvcm(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8,    intent(inout)::residual_states(neigmx)
complex*16,intent(out)  ::hsvcm(ncpx,ncpy,ncpz,neigmx,ncol)
real*8  dx,dy,dz
real*8  vecnor,tmp
logical l_convergeda
integer na,iiter,jiter,i,j,k,l,ix,iy,iz,ierr,iaps,ns2,nsv,nvef,ib,nb
real*8  aaa,bbb,alambda1
real*8  valre(nb0),valall(nb0),residual_statesold(nb0),residual_statesfirst(nb0)
logical l_converged(nb0)
integer nmat(nb0)
real*8  alambda(nb0)

complex*16,allocatable::rveccm(:,:,:,:,:,:),rvcm(:,:,:,:,:)
complex*16,allocatable::phicm(:,:,:,:,:,:),hphicm(:,:,:,:,:,:),sphicm(:,:,:,:,:,:)
complex*16,allocatable::vcm(:,:,:,:,:)
complex*16,allocatable::avcm(:,:,:,:,:)
complex*16,allocatable::asvcm(:,:,:,:,:)
complex*16,allocatable::workc(:,:,:,:)
complex*16,allocatable::cpsep(:,:,:,:),cphisep(:,:,:,:,:)
real*8,    allocatable::eigen(:,:),work(:),work2(:)
complex*16,allocatable::crmat(:,:,:),csmat(:,:,:),csmat0(:,:,:),cvect(:,:,:)
complex*16,allocatable::crmat2(:,:),cvect2(:,:),cwork(:),cwork2(:)
real*8,    allocatable::sumtmp(:,:),sumtmpall(:,:)
complex*16,allocatable::vcccm(:,:)
complex*16,allocatable::avc(:,:)
complex*16,allocatable::ctmp(:),ctmpall(:)

  nvef=1
  if (ncol==2) nvef=nspv
  nsv= min(ns1,nspv)

  dx=2.0d0*xmax/(ncpx*nprocx)
  dy=2.0d0*ymax/(ncpy*nprocy)
  dz=2.0d0*zmax/(ncpz*nprocz)

  allocate(vcm(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf,ncol,nb0))
  allocate(avcm(ncpx,ncpy,ncpz,ncol,nb0))
  allocate(asvcm(ncpx,ncpy,ncpz,ncol,nb0))
  allocate(workc(ncpx,ncpy,ncpz,ncol))
  allocate(cpsep(nprjmx,num_ppcell,ncol,nb0))
  allocate(cphisep(nprjmx,num_ppcell,ncol,0:nitermax,nb0))
  allocate(rveccm(ncpx,ncpy,ncpz,ncol,0:nitermax,nb0),rvcm(ncpx,ncpy,ncpz,ncol,nb0))
  allocate(phicm(ncpx,ncpy,ncpz,ncol,0:nitermax,nb0),sphicm(ncpx,ncpy,ncpz,ncol,0:nitermax,nb0) &
                                                    ,hphicm(ncpx,ncpy,ncpz,ncol,0:nitermax,nb0))
  allocate(sumtmp(4,nb0),sumtmpall(4,nb0))
  allocate(vcccm(num_list,ncol))
  allocate(avc(num_list,ncol))
  i=max(nitermax*(nitermax+1)*nb0,nprjmx*ncol*nb0)
  allocate(ctmp(i),ctmpall(i))

  call kslaplacian_initialize(nf)
  call overlap_finitedifference_init(ncpx,ncpy,ncpz,nf,ncol,1)

  jiter=0
  do while (jiter+1 < ndiismax)
  jiter=jiter+1

  if (myrank_glbl .eq. 0) then
    write(ndisp,*,err=9999) 'DIIS iter. #',jiter
    write(ndisp,*,err=9999) '  k ,spin,band,     eigen value    ,   residual norm   ,matrix dim.'
  end if

  do l=1,neigmx,nb0
  nb=nb0
  if (neigmx-l+1 < nb0) nb=neigmx-l+1

!$omp parallel default(shared) private(i)
    do ns2= 1,ncol
      do ib=1,nb
!$omp do
        do i=1,ncpx*ncpy*ncpz
          phicm( i,1,1,ns2,0,ib)=sveccm(i,1,1,l+ib-1,ns1+ns2-1,nk)
          sphicm(i,1,1,ns2,0,ib)=ssvcm(i,1,1,l+ib-1,ns1+ns2-1,nk)
        end do
!$omp do
        do i=1,nprjmx*num_ppcell
          cphisep(i,1,ns2,0,ib)=cspsep(i,1,l+ib-1,ns1+ns2-1,nk)
        end do
      end do
    end do
!$omp end parallel

  iiter=0
  l_converged(1:nb)= .false.
  l_convergeda= .false.
  do while ((iiter+1 < nitermax).and.(.not. l_convergeda))
  iiter=iiter+1

! ==========  compute eigenvalue and residual vector (rvec)  ==========
  do ib=1,nb
!$omp parallel default(shared) private(i,ix,iy,iz)
!$omp do
    do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)*ncol
      vcm(-nf+i,-nf+1,-nf+1,1,ib)=dcmplx(0.0d0,0.0d0)
    end do
!$omp do
    do i=1,ncpx*ncpy*ncpz*ncol
      avcm(i,1,1,1,ib)=dcmplx(0.0d0,0.0d0)
    end do
    do ns2= 1,ncol
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        vcm(ix,iy,iz,ns2,ib)=phicm(ix,iy,iz,ns2,iiter-1,ib)
      end do
      end do
      end do
    end do
!$omp end parallel
    if (iiter .eq. 1) then
!$omp parallel default(shared) private(i,ix,iy,iz)
!$omp do
      do i=1,nprjmx*num_ppcell*ncol
        cpsep(i,1,1,ib)=cphisep(i,1,1,iiter-1,ib)
      end do
!$omp single
      call overlap_finitedifference_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,ncol,vcm(-nf+1,-nf+1,-nf+1,1,ib))
!$omp end single
      call nonlocaloperation_c_02(1,nso,natom,num_spe,1,nums,nprjmx,ncol,ns1,1,num_list,num_ppcell, & ! <
                                  ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                              & ! <
                                  key_natpri_in,key_natpri_inps,key_soc_calc,                       & ! <
                                  dx,dy,dz,skpx,skpy,skpz,                                          & ! <
                                  nprj,indspe,natinf,lstvec2,natpri,naps,natsoc,                    & ! <
                                  lstx,lsty,lstz,natx,naty,natz,                                    & ! <
                                  vnlocp,cpsep(1,1,1,ib),dij,dijsoc,                                & ! <
                                  avcm(1,1,1,1,ib),                                                 & ! X
                                  vcccm,avc)                                                          ! W
!$omp barrier
!$omp single
      call overlap_fdcheck_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,ncol,vcm(-nf+1,-nf+1,-nf+1,1,ib))
!$omp end single
      call kslaplacian_c(ncpx,ncpy,ncpz,1,ncol,nvef,nf,dx,dy,dz,xmax,ymax,zmax,skpx,skpy,skpz &
                        ,vcm(-nf+1,-nf+1,-nf+1,1,ib),veff(1,1,1,nsv),1,avcm(1,1,1,1,ib),workc)
!$omp barrier
!$omp do
      do i=1,ncpx*ncpy*ncpz*ncol
        hphicm(i,1,1,1,iiter-1,ib)=avcm(i,1,1,1,ib)
      end do
!$omp end parallel
    end if
    valre(ib)=0.0d0
!$omp parallel default(shared)
    call orthogonalization_c_05(ncpx*ncpy*ncpz,1,1,ncol,1,phicm(1,1,1,1,iiter-1,ib),1,hphicm(1,1,1,1,iiter-1,ib),valre(ib))
!$omp end parallel
  end do !ib
  call mpi_allreduce(valre,valall,nb,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  do ib=1,nb
    valre(ib)=valall(ib)*dx*dy*dz
    sval(l+ib-1,ns1,nk)=valre(ib)
!$omp parallel default(shared)
    call ksdiis_c_01( &
     ncpx*ncpy*ncpz*ncol,-valre(ib),hphicm(1,1,1,1,iiter-1,ib),sphicm(1,1,1,1,iiter-1,ib), & ! <
     rveccm(1,1,1,1,iiter-1,ib))                                                       ! >
!$omp end parallel
  end do !ib
! =====================================================================

! ==========  precondition  ==========
!$omp parallel default(shared) private(i,ix,iy,iz)
  if (nprecon_diis .eq. 1) then
    do ib=1,nb
!$omp do
      do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)*ncol
        vcm(-nf+i,-nf+1,-nf+1,1,ib)=dcmplx(0.0d0,0.0d0)
      end do
      do ns2=1,ncol
!$omp do
        do iz=1,ncpz
        do iy=1,ncpy
        do ix=1,ncpx
          vcm(ix,iy,iz,ns2,ib)=rveccm(ix,iy,iz,ns2,iiter-1,ib)
        end do
        end do
        end do
      end do
!$omp single
      call overlap_finitedifference_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,1,ncol,vcm(-nf+1,-nf+1,-nf+1,1,ib))
      call overlap_fdcheck_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,1,ncol,vcm(-nf+1,-nf+1,-nf+1,1,ib))
!$omp end single
      call ksprecondition_c(ncpx,ncpy,ncpz,nf,ncol,vcm(-nf+1,-nf+1,-nf+1,1,ib),rvcm(1,1,1,1,ib))
    end do ! ib
  else
    do ib=1,nb
!$omp do
      do i=1,ncpx*ncpy*ncpz*ncol
        rvcm(i,1,1,1,ib)=rveccm(i,1,1,1,iiter-1,ib)
      end do
    end do
  end if
!$omp end parallel
! ====================================

  do ib=1,nb
!$omp parallel default(shared) private(i,ix,iy,iz)
!$omp do
    do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)*ncol
      vcm(-nf+i,-nf+1,-nf+1,1,ib)=dcmplx(0.0d0,0.0d0)
    end do
!$omp do
    do i=1,ncpx*ncpy*ncpz*ncol
      avcm(i,1,1,1,ib)=dcmplx(0.0d0,0.0d0)
    end do
    do ns2=1,ncol
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        vcm(ix,iy,iz,ns2,ib)=rvcm(ix,iy,iz,ns2,ib)
      end do
      end do
      end do
    end do
!$omp do
    do i=1,nprjmx*num_ppcell*ncol
      cpsep(i,1,1,ib)=dcmplx(0.0d0,0.0d0)
    end do
    call nonlocaloperation_c_01(natom,nprjmx,num_spe,num_ppcell,num_list,1,ncol,1, & ! <
                                ncpx,ncpy,ncpz,                                    & ! <
                                key_natpri_in,key_natpri_inps,                     & ! <
                                dx,dy,dz,skpx,skpy,skpz,                           & ! <
                                nprj,indspe,natinf,natpri,naps,lstvec2,            & ! <
                                lstx,lsty,lstz,natx,naty,natz,                     & ! <
                                vnlocp,rvcm(1,1,1,1,ib),                           & ! <
                                cpsep(1,1,1,ib),                                   & ! X
                                vcccm)                                               ! W
!$omp end parallel
  end do !ib

  do k=1,natom
    na=latom(k)
    if (natpri(na) .ge. key_natpri_inps) then
      iaps=naps(na)
      j=0
      do ns2=1,ncol
        do ib=1,nb
          do i=1,nprj(indspe(na))
            j=j+1
            ctmp(j)=cpsep(i,iaps,ns2,ib)
          end do
        end do
      end do
      call mpi_allreduce(ctmp,ctmpall,nprj(indspe(na))*nb*ncol,mpi_double_complex,mpi_sum,mpicom_atom(na),mpij)
      j=0
      do ns2=1,ncol
        do ib=1,nb
          do i=1,nprj(indspe(na))
            j=j+1
            cpsep(i,iaps,ns2,ib)=ctmpall(j)
          end do
        end do
      end do
    end if
  end do
  do ib=1,nb
    call overlap_finitedifference_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,ncol,vcm(-nf+1,-nf+1,-nf+1,1,ib))
!$omp parallel default(shared) private(i,ix,iy,iz)
    call nonlocaloperation_c_02(1,nso,natom,num_spe,1,nums,nprjmx,ncol,ns1,1,num_list,num_ppcell, & ! <
                                ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                              & ! <
                                key_natpri_in,key_natpri_inps,key_soc_calc,                       & ! <
                                dx,dy,dz,skpx,skpy,skpz,                                          & ! <
                                nprj,indspe,natinf,lstvec2,natpri,naps,natsoc,                    & ! <
                                lstx,lsty,lstz,natx,naty,natz,                                    & ! <
                                vnlocp,cpsep(1,1,1,ib),dij,dijsoc,                                & ! <
                                avcm(1,1,1,1,ib),                                                 & ! X
                                vcccm,avc)                                                          ! W
!$omp barrier
!$omp single
    call overlap_fdcheck_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,ncol,vcm(-nf+1,-nf+1,-nf+1,1,ib))
!$omp end single
    call kslaplacian_c(ncpx,ncpy,ncpz,1,ncol,nvef,nf,dx,dy,dz,xmax,ymax,zmax,skpx,skpy,skpz,vcm(-nf+1,-nf+1,-nf+1,1,ib) &
                      ,veff(1,1,1,nsv),1,avcm(1,1,1,1,ib),workc)
!$omp barrier
!$omp do
    do i=1,ncpx*ncpy*ncpz*ncol
      asvcm(i,1,1,1,ib)=rvcm(i,1,1,1,ib)
    end do
    call nonlocaloperation_c_02(0,0,natom,num_spe,1,nums,nprjmx,ncol,1,1,num_list,num_ppcell, & ! <
                                ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                          & ! <
                                key_natpri_in,key_natpri_inps,key_soc_calc,                   & ! <
                                dx,dy,dz,skpx,skpy,skpz,                                      & ! <
                                nprj,indspe,natinf,lstvec2,natpri,naps,natsoc(1),             & ! <
                                lstx,lsty,lstz,natx,naty,natz,                                & ! <
                                vnlocp,cpsep(1,1,1,ib),sss,dijsoc(1,1,1,1),                   & ! <
                                asvcm(1,1,1,1,ib),                                            & ! X
                                vcccm,avc)                                                      ! W
!$omp end parallel
  end do ! ib

  if (iiter .eq. 1) then
! ==========  compute \alambda ==========
  do ib=1,nb
    alambda1=alambda_diis
    sumtmp(1:4,ib)=0.0d0
!$omp parallel default(shared) private(i,ix,iy,iz)
!$omp do
    do i=1,nprjmx*num_ppcell*ncol
      cphisep(i,1,1,1,ib)=cphisep(i,1,1,0,ib)+alambda1*cpsep(i,1,1,ib)
    end do
    call ksdiis_c_01( &
     ncpx*ncpy*ncpz*ncol,alambda1,phicm(1,1,1,1,0,ib),rvcm(1,1,1,1,ib),   & ! <
     phicm(1,1,1,1,1,ib))                                       ! >
    call ksdiis_c_01( &
     ncpx*ncpy*ncpz*ncol,alambda1,hphicm(1,1,1,1,0,ib),avcm(1,1,1,1,ib),  & ! <
     hphicm(1,1,1,1,1,ib))                                      ! >
    call ksdiis_c_01( &
     ncpx*ncpy*ncpz*ncol,alambda1,sphicm(1,1,1,1,0,ib),asvcm(1,1,1,1,ib), & ! <
     sphicm(1,1,1,1,1,ib))                                      ! >
    call orthogonalization_c_05(ncpx*ncpy*ncpz,1,1,ncol,1,rvcm(1,1,1,1,ib),1,hphicm(1,1,1,1,0,ib),sumtmp(1,ib))
    call orthogonalization_c_05(ncpx*ncpy*ncpz,1,1,ncol,1,rvcm(1,1,1,1,ib),1,sphicm(1,1,1,1,0,ib),sumtmp(2,ib))
    call orthogonalization_c_05(ncpx*ncpy*ncpz,1,1,ncol,1,phicm(1,1,1,1,1,ib),1,hphicm(1,1,1,1,1,ib),sumtmp(3,ib))
    call orthogonalization_c_05(ncpx*ncpy*ncpz,1,1,ncol,1,phicm(1,1,1,1,1,ib),1,sphicm(1,1,1,1,1,ib),sumtmp(4,ib))
!$omp end parallel
  end do ! ib
  call mpi_allreduce(sumtmp,sumtmpall,4*nb,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  do ib=1,nb
    aaa=2.0d0*(sumtmpall(1,ib)-sval(l+ib-1,ns1,nk)*sumtmpall(2,ib))*dx*dy*dz
    bbb=(sumtmpall(3,ib)/sumtmpall(4,ib)-sval(l+ib-1,ns1,nk)-aaa*alambda1)/(alambda1*alambda1)
! I think alambda is in opposite direction to that in PRB54 11169 (1996)
! because the residual vector is in opposite sign to the conjugate vector of CG.
! 04/19/2006 T.Ono.
!  alambda=-0.5d0*aaa/bbb
    alambda(ib)=0.5d0*aaa/bbb
    if (alambda(ib) .gt. alambda_max) alambda(ib)=alambda_max
    if (alambda(ib) .lt. alambda_min) alambda(ib)=alambda_min
  end do ! ib
! ======================================
  end if !iiter

! ==========  compute trial phi ==========
!$omp parallel default(shared) private(i,ix,iy,iz)
  do ib=1,nb
    call ksdiis_c_01( &
     ncpx*ncpy*ncpz*ncol,alambda(ib),phicm(1,1,1,1,iiter-1,ib),rvcm(1,1,1,1,ib),   & ! <
     phicm(1,1,1,1,iiter,ib))                                        ! >
    call ksdiis_c_01( &
     ncpx*ncpy*ncpz*ncol,alambda(ib),hphicm(1,1,1,1,iiter-1,ib),avcm(1,1,1,1,ib),  & ! <
     hphicm(1,1,1,1,iiter,ib))                                       ! >
    call ksdiis_c_01( &
     ncpx*ncpy*ncpz*ncol,alambda(ib),sphicm(1,1,1,1,iiter-1,ib),asvcm(1,1,1,1,ib), & ! <
     sphicm(1,1,1,1,iiter,ib))                                       ! >
!$omp do
    do i=1,nprjmx*num_ppcell*ncol
      cphisep(i,1,1,iiter,ib)=cphisep(i,1,1,iiter-1,ib)+alambda(ib)*cpsep(i,1,1,ib)
    end do
  end do ! ib
!$omp end parallel
! ========================================

! ==========  normalize trial phi  ==========
!  vecnor=0.0d0
!!$omp parallel default(shared) private(i,ix,iy,iz)
!  call orthogonalization_c_05(ncpx*ncpy*ncpz,1,1,ncol,1,phicm(1,1,1,1,iiter),1,sphicm(1,1,1,1,iiter),vecnor)
!!$omp barrier
!!$omp single
!  call mpi_allreduce(vecnor,tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
!  vecnor=tmpall*dx*dy*dz
!  tmp=1.0d0/dsqrt(vecnor)
!!$omp end single
!  call orthogonalization_c_06(ncpx*ncpy*ncpz*ncol,phicm(1,1,1,1,iiter),tmp)
!  call orthogonalization_c_06(ncpx*ncpy*ncpz*ncol,hphicm(1,1,1,1,iiter),tmp)
!  call orthogonalization_c_06(ncpx*ncpy*ncpz*ncol,sphicm(1,1,1,1,iiter),tmp)
!!$omp barrier
!!$omp do
!  do i=1,nprjmx*num_ppcell*ncol
!    cphisep(i,1,1,iiter)=cphisep(i,1,1,iiter)*tmp
!  end do
!!$omp end parallel
! ===========================================

! ==========  compute eigenvalue and residual vector (rvec)  ==========
!  valre=0.0d0
!!$omp parallel default(shared)
!  call orthogonalization_c_05(ncpx*ncpy*ncpz,1,1,ncol,1,phicm(1,1,1,1,iiter) &
!                 ,1,hphicm(1,1,1,1,iiter),valre)
!!$omp end parallel
!  call mpi_allreduce(valre,valall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
!  valre=valall*dx*dy*dz
!$omp parallel default(shared) private(j)
  do ib=1,nb
    do j=0,iiter
      call ksdiis_c_01( &
       ncpx*ncpy*ncpz*ncol,-sval(l+ib-1,ns1,nk),hphicm(1,1,1,1,j,ib),sphicm(1,1,1,1,j,ib), & ! <
       rveccm(1,1,1,1,j,ib))                                                         ! >
    end do
  end do ! ib
!$omp end parallel
! =====================================================================

! ==========  diis step  ==========
! ----------  set up matrix  ----------
  allocate(eigen(iiter+1,nb),work((iiter+1)*3),cwork((iiter+1)*3))
  allocate(crmat(iiter+1,iiter+1,nb),csmat(iiter+1,iiter+1,nb),csmat0(iiter+1,iiter+1,nb),cvect(iiter+1,iiter+1,nb))
  do ib=1,nb
    k=0
    do j=0,iiter
      do i=0,j
        ctmp(k+1+(iiter+1)*(iiter+2)*(ib-1))=dcmplx(0.0d0,0.0d0)
        ctmp(k+2+(iiter+1)*(iiter+2)*(ib-1))=dcmplx(0.0d0,0.0d0)
!$omp parallel default(shared)
        call orthogonalization_c_04(ncpx*ncpy*ncpz,1,1,ncol,1,sphicm(1,1,1,1,j,ib),1,1,phicm(1,1,1,1,i,ib) &
          ,ctmp(k+1+(iiter+1)*(iiter+2)*(ib-1)))
        call orthogonalization_c_04(ncpx*ncpy*ncpz,1,1,ncol,1,rveccm(1,1,1,1,j,ib),1,1,rveccm(1,1,1,1,i,ib) &
          ,ctmp(k+2+(iiter+1)*(iiter+2)*(ib-1)))
!$omp end parallel
        k=k+2
      end do
    end do
  end do ! ib
  call mpi_allreduce(ctmp,ctmpall,(iiter+1)*(iiter+2)*nb,mpi_double_complex,mpi_sum,mpicom_space,mpij)
  do ib=1,nb
    k=0
    do j=0,iiter
      do i=0,j
        csmat(i+1,j+1,ib)=ctmpall(k+1+(iiter+1)*(iiter+2)*(ib-1))*dx*dy*dz
        csmat(j+1,i+1,ib)=dconjg(csmat(i+1,j+1,ib))
        crmat(i+1,j+1,ib)=ctmpall(k+2+(iiter+1)*(iiter+2)*(ib-1))*dx*dy*dz
        crmat(j+1,i+1,ib)=dconjg(crmat(i+1,j+1,ib))
        k=k+2
      end do
    end do
    vecnor=dreal(csmat(iiter+1,iiter+1,ib))
    tmp=1.0d0/dsqrt(vecnor)
!$omp parallel default(shared) private(i)
    call orthogonalization_c_06(ncpx*ncpy*ncpz*ncol,phicm(1,1,1,1,iiter,ib),tmp)
    call orthogonalization_c_06(ncpx*ncpy*ncpz*ncol,hphicm(1,1,1,1,iiter,ib),tmp)
    call orthogonalization_c_06(ncpx*ncpy*ncpz*ncol,sphicm(1,1,1,1,iiter,ib),tmp)
    call orthogonalization_c_06(ncpx*ncpy*ncpz*ncol,rveccm(1,1,1,1,iiter,ib),tmp)
!$omp do
    do i=1,nprjmx*num_ppcell*ncol
      cphisep(i,1,1,iiter,ib)=cphisep(i,1,1,iiter,ib)*tmp
    end do
!$omp end parallel
    do i=1,iiter+1
      csmat(iiter+1,i,ib)=csmat(iiter+1,i,ib)*tmp
      crmat(iiter+1,i,ib)=crmat(iiter+1,i,ib)*tmp
    end do
    do i=1,iiter+1
      csmat(i,iiter+1,ib)=csmat(i,iiter+1,ib)*tmp
      crmat(i,iiter+1,ib)=crmat(i,iiter+1,ib)*tmp
    end do
    csmat0(:,:,ib)=csmat(:,:,ib)
! -------------------------------------
    call zheev('v','u',iiter+1,csmat(1,1,ib),iiter+1,eigen(1,ib),cwork,(iiter+1)*3,work,ierr)
    do i=1,iiter+1
      work(i)=eigen(i,ib)
    end do
    do j=1,iiter+1
      do i=1,iiter+1
        cvect(i,j,ib)=csmat(i,iiter+2-j,ib)
      end do
      eigen(j,ib)=work(iiter+2-j)
    end do
! ----------  eliminate non positive eigenvalues  ----------
    nmat(ib)=0
    k=0
    do while (k < iiter+1)
      nmat(ib)=nmat(ib)+1
      k=k+1
      if (eigen(nmat(ib),ib) .lt. eps_eig_diis) then
        do i=nmat(ib),iiter-1
          cvect(:,i,ib)=cvect(:,i+1,ib)
          eigen(i,ib)=eigen(i+1,ib)
        end do
        nmat(ib)=nmat(ib)-1
!        write(6,*) 'smat is not positive definite!!!!'
      end if
    end do
! ----------------------------------------------------------
    if (nmat(ib) .gt. 1) then
! ----------  set up new matrix  ----------
      do j=1,nmat(ib)
        do i=1,nmat(ib)
          cvect(i,j,ib)=cvect(i,j,ib)/dsqrt(eigen(j,ib))
        end do
      end do
      allocate(crmat2(nmat(ib),nmat(ib)),cvect2(nmat(ib),nmat(ib)),work2(nmat(ib)*3),cwork2(nmat(ib)*3))
      csmat(:,:,ib)=dcmplx(0.0d0,0.0d0)
      do k=1,iiter+1
        do j=1,nmat(ib)
          do i=1,iiter+1
            csmat(k,j,ib)=csmat(k,j,ib)+crmat(k,i,ib)*cvect(i,j,ib)
          end do
        end do
      end do
      crmat2=dcmplx(0.0d0,0.0d0)
      do k=1,nmat(ib)
        do j=1,nmat(ib)
          do i=1,iiter+1
            crmat2(k,j)=crmat2(k,j)+dconjg(cvect(i,k,ib))*csmat(i,j,ib)
          end do
        end do
      end do
! -----------------------------------------
      ierr=0
      call zheev('v','u',nmat(ib),crmat2,nmat(ib),eigen(1,ib),cwork2,nmat(ib)*3,work2,ierr)
      do i=1,nmat(ib)
        work2(i)=eigen(i,ib)
      end do
      do j=1,nmat(ib)
        do i=1,nmat(ib)
          cvect2(i,j)=crmat2(i,nmat(ib)+1-j)
        end do
        eigen(j,ib)=work2(nmat(ib)+1-j)
      end do
      if ((ierr/=0) .and. (myrx**2+myry**2+myrz**2==0)) &
         write(ndisp,*) 'error occurs in dsyev of diis!!','myr_kpt=',myr_kpt,'nk=',nk,'ns=',ns1,'l=',l
!    write(ndisp,*) 'eigen',(eigen(i),i=1,iiter+1)
      if (iiter .gt. 1) residual_statesold(ib)=residual_states(l+ib-1)
      if (iiter .eq. 2) residual_statesfirst(ib)=residual_states(l+ib-1)/residual_statesold(ib)
      residual_states(l+ib-1)=eigen(nmat(ib),ib)
      if (iiter .gt. 2) then
!   Be careful!! If iiter=1, residual_statesfirst=0.
        if (residual_states(l+ib-1)/residual_statesold(ib)/residual_statesfirst(ib) .lt. ratio_diis) then
          iiter=iiter-1
          nmat(ib)=nmat(ib)-1
          residual_states(l+ib-1)=residual_statesold(ib)
          deallocate(crmat2,cvect2,work2,cwork2)
          l_converged(ib)= .true.
          goto 1010
        end if
      end if
! ----------  replicate eigenvectors  ----------
      crmat(:,:,ib)=dcmplx(0.0d0,0.0d0)
      do k=1,nmat(ib)
        do j=1,iiter+1
          do i=1,nmat(ib)
            crmat(j,k,ib)=crmat(j,k,ib)+cvect(j,i,ib)*cvect2(i,k)
          end do
        end do
      end do
      tmp=0.0d0
      do j=1,iiter+1
        do i=1,iiter+1
          tmp=tmp+dreal(dconjg(crmat(i,nmat(ib),ib))*csmat0(i,j,ib)*crmat(j,nmat(ib),ib))
        end do
      end do
      crmat(:,nmat(ib),ib)=crmat(:,nmat(ib),ib)/dsqrt(tmp)
      deallocate(crmat2,cvect2,work2,cwork2)
! ----------------------------------------------
    else
      crmat(:,nmat(ib),ib)=dcmplx(0.0d0,0.0d0)
      crmat(iiter+1,nmat(ib),ib)=dcmplx(1.0d0,0.0d0)
    end if !(nmat(ib) .gt. 1)
! ----------  compute new wave function  ----------
!$omp parallel default(shared) private(i)
!$omp do
    do i=1,ncpx*ncpy*ncpz*ncol
      avcm(i,1,1,1,ib)=dcmplx(0.0d0,0.0d0)
    end do
    call ksdiis_c_03( &
     ncpx*ncpy*ncpz*ncol,iiter+1,phicm(1,1,1,1,0,ib),crmat(1,nmat(ib),ib), & ! <
     avcm(1,1,1,1,ib))                                                          ! X
!$omp barrier
!$omp do
    do i=1,ncpx*ncpy*ncpz*ncol
      phicm(i,1,1,1,iiter,ib)=avcm(i,1,1,1,ib)
      avcm(i,1,1,1,ib)=dcmplx(0.0d0,0.0d0)
    end do
    call ksdiis_c_03( &
     ncpx*ncpy*ncpz*ncol,iiter+1,sphicm(1,1,1,1,0,ib),crmat(1,nmat(ib),ib), & ! <
     avcm(1,1,1,1,ib))                                                           ! X
!$omp barrier
!$omp do
    do i=1,ncpx*ncpy*ncpz*ncol
      sphicm(i,1,1,1,iiter,ib)=avcm(i,1,1,1,ib)
      avcm(i,1,1,1,ib)=dcmplx(0.0d0,0.0d0)
    end do
    call ksdiis_c_03( &
     ncpx*ncpy*ncpz*ncol,iiter+1,hphicm(1,1,1,1,0,ib),crmat(1,nmat(ib),ib), & ! <
     avcm(1,1,1,1,ib))                                                           ! X
!$omp barrier
!$omp do
    do i=1,ncpx*ncpy*ncpz*ncol
      hphicm(i,1,1,1,iiter,ib)=avcm(i,1,1,1,ib)
    end do
!$omp do
    do i=1,nprjmx*num_ppcell*ncol
      cpsep(i,1,1,ib)=dcmplx(0.0d0,0.0d0)
    end do
    do j=1,iiter+1
!$omp do
      do i=1,nprjmx*num_ppcell*ncol
        cpsep(i,1,1,ib)=cpsep(i,1,1,ib)+crmat(j,nmat(ib),ib)*cphisep(i,1,1,j-1,ib)
      end do
    end do
!$omp do
    do i=1,nprjmx*num_ppcell*ncol
      cphisep(i,1,1,iiter,ib)=cpsep(i,1,1,ib)
    end do
!$omp end parallel
! -------------------------------------------------
! ----------  check for DIIS  ----------
!  avre=0.0d0
!  avim=0.0d0
!!$omp parallel default(shared)
!  call ksdiis_c_03(ncpx*ncpy*ncpz,iiter+1,rveccm(1,1,1,0),crmat(1,nmat(ib)),avcm)
!!$omp end parallel
!  rveccm(:,:,:,:,iiter)=avcm
!  vecnor=0.0d0
!!$omp parallel default(shared)
!  call orthogonalization_c_05(ncpx*ncpy*ncpz,1,1,ncol,1,rveccm(1,1,1,1,iiter),1,rveccm(1,1,1,1,iiter),vecnor)
!!$omp end parallel
!  call mpi_allreduce(vecnor,tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
!  vecnor=tmpall*dx*dy*dz
!  write(ndisp,*) 'norm of residual vector=',vecnor
! --------------------------------------

  1010 continue
  end do ! ib

  deallocate(eigen,work,cwork)
  deallocate(crmat,csmat,csmat0,cvect)
! ======================================

!  if (myrank_glbl .eq. 0) then
!  write (ndisp,9002,err=9999) nk,ns1,l,sval(l,ns1,nk),dsqrt(residual_states(l)),nmat(ib)
!  end if

  l_convergeda=.true.
  do ib=1,nb
    if (.not. (l_converged(ib))) l_convergeda=.false.
  end do

  end do ! iiter

  do ib=1,nb
!$omp parallel default(shared) private(i,ix,iy,iz)
    do ns2= 1,ncol
!$omp do
      do i=1,ncpx*ncpy*ncpz
        sveccm(i,1,1,l+ib-1,ns1+ns2-1,nk)=phicm(i,1,1,ns2,iiter,ib)
        ssvcm(i,1,1,l+ib-1,ns1+ns2-1,nk)=sphicm(i,1,1,ns2,iiter,ib)
        hsvcm(i,1,1,l+ib-1,ns2)=hphicm(i,1,1,ns2,iiter,ib)
      end do
!$omp do
      do i=1,nprjmx*num_ppcell
        cspsep(i,1,l+ib-1,ns1+ns2-1,nk)=cphisep(i,1,ns2,iiter,ib)
      end do
    end do
!$omp end parallel

    if (myrank_glbl .eq. 0) then
      write (ndisp,'(i4,1x,i4,1x,i6,1x,2e20.10,1x,i3)',err=9999) &
       nk,ns1,l+ib-1,sval(l+ib-1,ns1,nk),dsqrt(residual_states(l+ib-1)),nmat(ib)
    end if

    ksitmax= max(iiter,ksitmax)
    ksconv = ksconv .and. (dsqrt(residual_states(l+ib-1)) <= epssd)

    end do ! ib
  end do ! l
! ***********************************************

  if (l_hsv) then

    if (myrank_glbl .eq. 0) then
      write (ndisp,*) 'orthogonalization was implemented.'
    end if

! **********  orthonormalize wave function (svec)  **********
    do i=1,northo
      call orthogonalization_c_01(natom,neigmx,nprjmx,nums,ncol,numk,nk,ns1,0,num_list,num_ppcell,  & ! <
                                  num_spe,key_ortho_cmpt_innerproduct,ndisp,                        & ! <
                                  ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                              & ! <
                                  key_ortho_cmpt_innerproduct,                                      & ! <
                                  key_natpri_in,key_natpri_inps,                                    & ! <
                                  dx,dy,dz,skpx,skpy,skpz,                                          & ! <
                                  nprj,indspe,natpri,naps,natinf,lstvec2,latom,                     & ! <
                                  lstx,lsty,lstz,natx,naty,natz,                                    & ! <
                                  vnlocp,sss,                                                       & ! <
                                  sveccm,ssvcm,cspsep)                                                ! X
    end do
! ***********************************************************

! **********  compute innerproduct between H and \psi  **********
!$omp parallel default(shared) private(i,ix,iy,iz)
    do l=1,neigmx
      do ns2=1,ncol
!$omp do
        do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)
          vcm(-nf+i,-nf+1,-nf+1,ns2,1)=dcmplx(0.0d0,0.0d0)
        end do
!$omp do
        do i=1,ncpx*ncpy*ncpz
          hsvcm(i,1,1,l,ns2)=dcmplx(0.0d0,0.0d0)
        end do
!$omp do
        do iz=1,ncpz
        do iy=1,ncpy
        do ix=1,ncpx
          vcm(ix,iy,iz,ns2,1)=sveccm(ix,iy,iz,l,ns1+ns2-1,nk)
        end do
        end do
        end do
!$omp do
        do i=1,nprjmx*num_ppcell
          cpsep(i,1,ns2,1)=cspsep(i,1,l,ns1+ns2-1,nk)
        end do
      end do ! ns2
!$omp single
      call overlap_finitedifference_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,ncol,vcm(-nf+1,-nf+1,-nf+1,1,1))
!$omp end single
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
      call overlap_fdcheck_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,ncol,vcm(-nf+1,-nf+1,-nf+1,1,1))
!$omp end single
      call kslaplacian_c(ncpx,ncpy,ncpz,neigmx,ncol,nvef,nf,dx,dy,dz,xmax,ymax,zmax,skpx,skpy,skpz, &
       vcm(-nf+1,-nf+1,-nf+1,1,1),veff(1,1,1,nsv),l,hsvcm,workc)
!$omp barrier
    end do ! l
!$omp end parallel
! ***************************************************************

  end if

  end do !jiter

  call kslaplacian_finalize
  call overlap_finitedifference_final

  deallocate(vcm)
  deallocate(avcm)
  deallocate(asvcm)
  deallocate(workc)
  deallocate(cpsep)
  deallocate(cphisep)
  deallocate(rveccm,rvcm)
  deallocate(phicm,sphicm,hphicm)
  deallocate(sumtmp,sumtmpall)
  deallocate(vcccm)
  deallocate(avc)
  deallocate(ctmp,ctmpall)

  return
9999 continue
  call mpi_abort(mpi_comm_world)
end subroutine


subroutine ksdiis_r_01( &
 n,d,a,b, & ! <
 c)         ! >
implicit none
integer,intent(in) ::n
real*8, intent(in) ::d
real*8, intent(in) ::a(n),b(n)
real*8, intent(out)::c(n)
integer i
!$omp do
  do i=1,n
    c(i)=a(i)+d*b(i)
  end do
  return
end subroutine ksdiis_r_01


subroutine ksdiis_c_01( &
 n,d,a,b, & ! <
 c)         ! >
implicit none
integer,   intent(in) ::n
real*8,    intent(in) ::d
complex*16,intent(in) ::a(n),b(n)
complex*16,intent(out)::c(n)
integer i
!$omp do
  do i=1,n
    c(i)=a(i)+d*b(i)
  end do
  return
end subroutine


subroutine ksdiis_r_03( &
 n,l,b,r, & ! <
 a)         ! X
implicit none
integer,intent(in)   ::n,l
real*8, intent(in)   ::b(n,l),r(l)
real*8, intent(inout)::a(n)
integer i,j
  j=0
  do while (j<l)
    select case(l-j)
    case (2:)
!$omp do
      do i=1,n
        a(i)=a(i)+r(j+1)*b(i,j+1)+r(j+2)*b(i,j+2)
      end do
!$omp end do nowait
      j=j+2
    case (1)
!$omp do
      do i=1,n
        a(i)=a(i)+r(j+1)*b(i,j+1)
      end do
!$omp end do nowait
      j=j+1
    end select
  end do
  return
end subroutine


subroutine ksdiis_c_03( &
 n,l,cb,c, & ! <
 ca)         ! X
implicit none
integer,   intent(in)   ::n,l
complex*16,intent(in)   ::cb(n,l),c(l)
complex*16,intent(inout)::ca(n)
integer i,j
  j=0
  do while (j<l)
    select case(l-j)
    case (2:)
!$omp do
      do i=1,n
        ca(i)=ca(i)+c(j+1)*cb(i,j+1)+c(j+2)*cb(i,j+2)
      end do
!$omp end do nowait
      j=j+2
    case (1)
!$omp do
      do i=1,n
        ca(i)=ca(i)+c(j+1)*cb(i,j+1)
      end do
!$omp end do nowait
      j=j+1
    end select
  end do
  return
end subroutine


end module
