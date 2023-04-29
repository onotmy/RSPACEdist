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
! **********  orthogonalization8f.F90 11/21/2013-01  **********

module mod_orthogonalization
implicit none

contains

subroutine orthogonalization_r(natom,neigmx,nprjmx,nums,numk,num_ppcell,num_spe,northo,num_list,nflag,ndisp, & ! <
                       ncpx,ncpy,ncpz,                                                                       & ! <
                       key_ortho_cmpt_innerproduct,                                                          & ! <
                       key_natpri_in,key_natpri_inps,                                                        & ! <
                       dx,dy,dz,                                                                             & ! <
                       nprj,indspe,natpri,naps,natinf,lstvec2,latom,ntyppp,                                  & ! <
                       vnlocp,sss,                                                                           & ! <
                       svecre,ssvre,rspsep)                                                                    ! X
implicit none
integer,intent(in)::natom,neigmx,nprjmx,nums,numk,num_ppcell,num_spe,northo,num_list,nflag,ndisp
integer,intent(in)::ncpx,ncpy,ncpz
integer,intent(in)::key_ortho_cmpt_innerproduct
integer,intent(in)::key_natpri_in,key_natpri_inps
real*8, intent(in)::dx,dy,dz
integer,intent(in)::indspe(natom),natpri(natom),naps(natom),nprj(num_spe)
integer,intent(in)::natinf(natom),lstvec2(num_list,num_ppcell),latom(natom),ntyppp(num_spe)
real*8, intent(in)::vnlocp(num_list,nprjmx,num_ppcell),sss(nprjmx,nprjmx,num_spe)
real*8, intent(inout)::svecre(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8, intent(inout)::ssvre(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8, intent(inout)::rspsep(nprjmx,num_ppcell,neigmx,nums,numk)
integer nk,ns,loop
  do nk=1,numk
    do ns=1,nums
      do loop=1,northo
        call orthogonalization_r_01(natom,neigmx,nprjmx,nums,numk,nk,ns,0,num_list,num_ppcell,num_spe,nflag,ndisp, & ! <
                                    ncpx,ncpy,ncpz,                                                                & ! <
                                    key_ortho_cmpt_innerproduct,                                                   & ! <
                                    key_natpri_in,key_natpri_inps,                                                 & ! <
                                    dx,dy,dz,                                                                      & ! <
                                    nprj,indspe,natpri,naps,natinf,lstvec2,latom,                                  & ! <
                                    vnlocp,sss,                                                                    & ! <
                                    svecre,ssvre,rspsep)                                                             ! X
      end do
    end do
  end do
  return
end subroutine


subroutine orthogonalization_c(natom,neigmx,nprjmx,nums,ncol,numk,num_ppcell,num_spe,northo,num_list,nflag,ndisp, & ! <
                               ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                                               & ! <
                               key_ortho_cmpt_innerproduct,                                                       & ! <
                               key_natpri_in,key_natpri_inps,                                                     & ! <
                               dx,dy,dz,                                                                          & ! <
                               nprj,indspe,natpri,naps,natinf,lstvec2,latom,ntyppp,                               & ! <
                               lstx,lsty,lstz,natx,naty,natz,                                                     & ! <
                               skpxx,skpyy,skpzz,                                                                 & ! <
                               vnlocp,sss,                                                                        & ! <
                               sveccm,ssvcm,cspsep)                                                                 ! X
implicit none
integer,   intent(in)::natom,neigmx,nprjmx,nums,ncol,numk,num_ppcell,num_spe,northo,num_list,nflag,ndisp
integer,   intent(in)::ncpx,ncpy,ncpz,npxmax,npymax,npzmax
integer,   intent(in)::key_ortho_cmpt_innerproduct
integer,   intent(in)::key_natpri_in,key_natpri_inps
real*8,    intent(in)::dx,dy,dz
real*8,    intent(in)::skpxx(numk),skpyy(numk),skpzz(numk)
integer,   intent(in)::indspe(natom),natpri(natom),naps(natom),nprj(num_spe)
integer,   intent(in)::natinf(natom),lstvec2(num_list,num_ppcell),latom(natom),ntyppp(num_spe)
integer,   intent(in)::lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell)
integer,   intent(in)::natx(natom),naty(natom),natz(natom)
real*8,    intent(in)::vnlocp(num_list,nprjmx,num_ppcell),sss(nprjmx,nprjmx,num_spe)
complex*16,intent(inout)::sveccm(ncpx,ncpy,ncpz,neigmx,nums,numk)
complex*16,intent(inout)::ssvcm(ncpx,ncpy,ncpz,neigmx,nums,numk)
complex*16,intent(inout)::cspsep(nprjmx,num_ppcell,neigmx,nums,numk)
integer nk,ns,loop

  do nk=1,numk
    do ns=1,nums+1-ncol
      do loop=1,northo
        call orthogonalization_c_01(natom,neigmx,nprjmx,nums,ncol,numk,nk,ns,0,num_list,num_ppcell,num_spe,nflag,ndisp, & ! <
                                    ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                                                & ! <
                                    key_ortho_cmpt_innerproduct,                                                        & ! <
                                    key_natpri_in,key_natpri_inps,                                                      & ! <
                                    dx,dy,dz,skpxx(nk),skpyy(nk),skpzz(nk),                                             & ! <
                                    nprj,indspe,natpri,naps,natinf,lstvec2,latom,                                       & ! <
                                    lstx,lsty,lstz,natx,naty,natz,                                                      & ! <
                                    vnlocp,sss,                                                                         & ! <
                                    sveccm,ssvcm,cspsep)                                                                  ! X
      end do
    end do
  end do
  return
end subroutine


subroutine orthogonalization_r_01(natom,neigmx,nprjmx,nums,numk,nk,ns,l_inp,num_list,num_ppcell,num_spe,nflag,ndisp, & ! <
                                  ncpx,ncpy,ncpz,                                                                    & ! <
                                  key_ortho_cmpt_innerproduct,                                                       & ! <
                                  key_natpri_in,key_natpri_inps,                                                     & ! <
                                  dx,dy,dz,                                                                          & ! <
                                  nprj,indspe,natpri,naps,natinf,lstvec2,latom,                                      & ! <
                                  vnlocp,sss,                                                                        & ! <
                                  svecre,ssvre,rspsep)                                                                 ! X
use mod_stopp
use mod_nonlocaloperation, only:nonlocaloperation_r_01,nonlocaloperation_r_02
use mod_disptime
implicit none
integer,intent(in)::natom,neigmx,nprjmx,nums,numk,nk,ns,l_inp,num_list,num_ppcell,num_spe,nflag,ndisp
integer,intent(in)::ncpx,ncpy,ncpz
integer,intent(in)::key_ortho_cmpt_innerproduct
integer,intent(in)::key_natpri_in,key_natpri_inps
integer,intent(in)::indspe(natom),natpri(natom),naps(natom),nprj(num_spe)
integer,intent(in)::natinf(natom),lstvec2(num_list,num_ppcell),latom(natom)
real*8, intent(in)::dx,dy,dz
real*8, intent(in)::vnlocp(num_list,nprjmx,num_ppcell),sss(nprjmx,nprjmx,num_spe)
real*8, intent(inout)::svecre(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8, intent(inout)::ssvre(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8, intent(inout)::rspsep(nprjmx,num_ppcell,neigmx,nums,numk)
integer na,iaps,i,j,l,l0,l_stt,l_end
real*8 tmp,tmpall,vecnor
real*8,allocatable::avreal(:),avrall(:),amatr(:,:),amatr_a(:,:)
real*8,allocatable::rtmp(:,:),rtmpall(:,:)
real*8,allocatable::rpsep(:,:)
real*8,allocatable::avr(:)
real*8 stime
integer mstep,nblas03,nblas04

  nblas03=500000
  nblas04=500000
  mstep=32
  if (neigmx .gt. 1024) mstep=64
  if (neigmx .gt. 2048) mstep=128
  if (neigmx .gt. 4096) mstep=256

  allocate(avreal(neigmx),avrall(neigmx),amatr(mstep,mstep),amatr_a(mstep,mstep))
  allocate(avr(num_list))

  if (l_inp .eq. 0) then
    l_stt=1
    l_end=neigmx
  else
    l_stt=l_inp
    l_end=l_inp
  end if

  call starttime(stime)
  if (nflag .eq. key_ortho_cmpt_innerproduct) then
    allocate(rpsep(nprjmx,num_ppcell),rtmp(nprjmx*neigmx,num_ppcell),rtmpall(nprjmx*neigmx,num_ppcell))
    !$omp parallel default(shared) private(i)
    !$omp do
    do i=1,nprjmx*neigmx*num_ppcell
      rtmpall(i,1)=0.0d0
    end do
    !$omp end parallel
    do l=l_stt,l_end
      !$omp parallel default(shared) private(i)
      !$omp do
      do i=1,nprjmx*num_ppcell
        rpsep(i,1)=0.0d0
      end do
      call nonlocaloperation_r_01(natom,nprjmx,num_spe,num_ppcell,num_list,     & ! <
                                  ncpx,ncpy,ncpz,                               & ! <
                                  key_natpri_in,key_natpri_inps,                & ! <
                                  nprj,indspe,natinf,natpri,naps,lstvec2,       & ! <
                                  vnlocp,svecre(1,1,1,l,ns,nk),                 & ! <
                                  rpsep)                                          ! X
      !$omp end parallel
      do i=1,natom
        na=latom(i)
        if ((natpri(na) == key_natpri_in) .or. (natpri(na) == key_natpri_inps)) then
          iaps=naps(na)
          do mpij=1,nprj(indspe(na))
            rtmp(mpij+nprj(indspe(na))*(l-l_stt),iaps)=rpsep(mpij,iaps)
          end do
        end if
      end do
    end do !l
    do i=1,natom
      na=latom(i)
      if ((natpri(na) == key_natpri_in) .or. (natpri(na) == key_natpri_inps)) then
        iaps=naps(na)
        call mpi_allreduce(rtmp(1,iaps),rtmpall(1,iaps),nprj(indspe(na))*(l_end-l_stt+1),mpi_double_precision,mpi_sum,mpicom_atom(na),mpij)
      end if
    end do
    do l=l_stt,l_end
      do i=1,natom
        na=latom(i)
        if ((natpri(na) == key_natpri_in) .or. (natpri(na) == key_natpri_inps)) then
          iaps=naps(na)
          do mpij=1,nprj(indspe(na))
            rspsep(mpij,iaps,l,ns,nk)=rtmpall(mpij+nprj(indspe(na))*(l-l_stt),iaps)
          end do
        end if
      end do ! i
      !$omp parallel default(shared) private(i)
      !$omp do
      do i=1,nprjmx*num_ppcell
        rpsep(i,1)=rspsep(i,1,l,ns,nk)
      end do
      !$omp do
      do i=1,ncpx*ncpy*ncpz
        ssvre(i,1,1,l,ns,nk)=svecre(i,1,1,l,ns,nk)
      end do
      call nonlocaloperation_r_02(0,natom,num_spe,nums,nprjmx,1,num_list,num_ppcell, & ! <
                                  ncpx,ncpy,ncpz,                                    & ! <
                                  key_natpri_in,key_natpri_inps,                     & ! <
                                  dx,dy,dz,                                          & ! <
                                  nprj,indspe,natinf,lstvec2,natpri,naps,            & ! <
                                  vnlocp,rpsep,sss,                                  & ! <
                                  ssvre(1,1,1,l,ns,nk),                              & ! X
                                  avr)                                                 ! W
      !$omp end parallel
    end do !l
    deallocate(rpsep,rtmp,rtmpall)
  end if ! (nflag .eq. key_ortho_cmpt_innerproduct)
  if (nflag .eq. key_ortho_cmpt_innerproduct) call endtime(ndisp,stime,'[TI_sub] setup <p|\psi> :')

  l=l_stt
  do while (l < l_end+1)
    if ((l>mstep).and.(l_end-l+1>=mstep)) then
      do j=1,l-1,mstep
        !$omp parallel default(shared) private(i)
        !$omp do
        do i=1,mstep*mstep
          amatr(i,1)=0.0d0
        end do
        !$omp end parallel
        call dgemm('t','n',mstep,mstep,ncpx*ncpy*ncpz,1.0d0,svecre(1,1,1,l,ns,nk),ncpx*ncpy*ncpz, &
                   ssvre(1,1,1,j,ns,nk),ncpx*ncpy*ncpz,0.0d0,amatr,mstep)
        call mpi_allreduce(amatr,amatr_a,mstep*mstep,mpi_double_precision,mpi_sum,mpicom_space,mpij)
        !$omp parallel default(shared) private(i)
        !$omp do
        do i=1,mstep*mstep
          amatr(i,1)=amatr_a(i,1)*dx*dy*dz
        end do
        !$omp end parallel
        call dgemm('n','t',ncpx*ncpy*ncpz,mstep,mstep,1.0d0,svecre(1,1,1,j,ns,nk),ncpx*ncpy*ncpz, &
                   -amatr,mstep,1.0d0,svecre(1,1,1,l,ns,nk),ncpx*ncpy*ncpz)
        call dgemm('n','t',ncpx*ncpy*ncpz,mstep,mstep,1.0d0,ssvre(1,1,1,j,ns,nk),ncpx*ncpy*ncpz, &
                   -amatr,mstep,1.0d0,ssvre(1,1,1,l,ns,nk),ncpx*ncpy*ncpz)
        !$omp parallel default(shared)
        call multr(nprjmx*num_ppcell,mstep,mstep,rspsep(1,1,l,ns,nk),rspsep(1,1,j,ns,nk),amatr)
        !$omp end parallel
      end do
      vecnor=0.0d0
      !$omp parallel default(shared) private(i)
      call orthogonalization_r_05(ncpx*ncpy*ncpz,svecre(1,1,1,l,ns,nk),ssvre(1,1,1,l,ns,nk),vecnor)
      !$omp barrier
      !$omp single
      call mpi_allreduce(vecnor,tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
      vecnor=tmpall*dx*dy*dz
      if (vecnor<=0.0d0) call stopp('orthogonalization_r_01 (1): wf has negative norm')
      tmp=1.0d0/dsqrt(vecnor)
      !$omp end single
      call orthogonalization_r_06(ncpx*ncpy*ncpz,svecre(1,1,1,l,ns,nk),tmp)
      call orthogonalization_r_06(ncpx*ncpy*ncpz,ssvre(1,1,1,l,ns,nk),tmp)
      !$omp do
      do i=1,nprjmx*num_ppcell
        rspsep(i,1,l,ns,nk)=rspsep(i,1,l,ns,nk)*tmp
      end do
      !$omp end parallel
      do l0=l+1,l+mstep-1
        !$omp parallel default(shared) private(i)
        !$omp do
        do i=1,neigmx
          avreal(i)=0.0d0
        end do
        !$omp end parallel
        if (l0-1 > nblas04) then
          call dgemv('t',ncpx*ncpy*ncpz,l0-l,1.0d0,ssvre(1,1,1,l,ns,nk),ncpx*ncpy*ncpz &
              ,svecre(1,1,1,l0,ns,nk),1,0.0d0,avreal,1)
        else
        !$omp parallel default(shared)
          call orthogonalization_r_04(ncpx*ncpy*ncpz,neigmx,neigmx,l0,svecre(1,1,1,1,ns,nk) &
                          ,l,l0-1,ssvre(1,1,1,1,ns,nk),avreal)
        !$omp end parallel
        end if
        call mpi_allreduce(avreal,avrall,l0-l,mpi_double_precision,mpi_sum,mpicom_space,mpij)
        !$omp parallel default(shared) private(i)
        !$omp do
        do i=1,neigmx
          avreal(i)=avrall(i)*dx*dy*dz
        end do
        !$omp end parallel
        if (l0-l > nblas03) then
          call dgemv('n',ncpx*ncpy*ncpz,l0-l,-1.0d0,svecre(1,1,1,l,ns,nk),ncpx*ncpy*ncpz &
              ,avreal,1,1.0d0,svecre(1,1,1,l0,ns,nk),1)
          call dgemv('n',ncpx*ncpy*ncpz,l0-l,-1.0d0,ssvre(1,1,1,l,ns,nk),ncpx*ncpy*ncpz &
              ,avreal,1,1.0d0,ssvre(1,1,1,l0,ns,nk),1)
          call dgemv('n',nprjmx*num_ppcell,l0-l,-1.0d0,rspsep(1,1,l,ns,nk),nprjmx*num_ppcell &
              ,avreal,1,1.0d0,rspsep(1,1,l0,ns,nk),1)
        else
          !$omp parallel default(shared)
          call orthogonalization_r_03(ncpx*ncpy*ncpz,neigmx,l0-l,svecre(1,1,1,l0,ns,nk),svecre(1,1,1,l,ns,nk),avreal)
          call orthogonalization_r_03(ncpx*ncpy*ncpz,neigmx,l0-l,ssvre(1,1,1,l0,ns,nk),ssvre(1,1,1,l,ns,nk),avreal)
          call multr(nprjmx*num_ppcell,1,l0-l,rspsep(1,1,l0,ns,nk),rspsep(1,1,l,ns,nk),avreal)
          !$omp end parallel
        end if
        vecnor=0.0d0
        !$omp parallel default(shared) private(i)
        call orthogonalization_r_05(ncpx*ncpy*ncpz,svecre(1,1,1,l0,ns,nk),ssvre(1,1,1,l0,ns,nk),vecnor)
        !$omp barrier
        !$omp single
        tmpall=0.0d0
        call mpi_allreduce(vecnor,tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
        vecnor=tmpall*dx*dy*dz
        if (vecnor<=0.0d0) call stopp('orthogonalization_r_01 (2): wf has negative norm')
        tmp=1.0d0/dsqrt(vecnor)
        !$omp end single
        call orthogonalization_r_06(ncpx*ncpy*ncpz,svecre(1,1,1,l0,ns,nk),tmp)
        call orthogonalization_r_06(ncpx*ncpy*ncpz,ssvre(1,1,1,l0,ns,nk),tmp)
        !$omp do
        do i=1,nprjmx*num_ppcell
          rspsep(i,1,l0,ns,nk)=rspsep(i,1,l0,ns,nk)*tmp
        end do
        !$omp end parallel
      end do
      l=l+mstep
    else
      !$omp parallel default(shared) private(i)
      !$omp do
      do i=1,neigmx
        avreal(i)=0.0d0
      end do
      !$omp end parallel
      if (l > nblas04) then
        call dgemv('t',ncpx*ncpy*ncpz,l-1,1.0d0,ssvre(1,1,1,1,ns,nk),ncpx*ncpy*ncpz &
            ,svecre(1,1,1,l,ns,nk),1,0.0d0,avreal,1)
      else
      !$omp parallel default(shared)
        call orthogonalization_r_04(ncpx*ncpy*ncpz,neigmx,neigmx,l,svecre(1,1,1,1,ns,nk),1,l-1,ssvre(1,1,1,1,ns,nk),avreal)
      !$omp end parallel
      end if
      call mpi_allreduce(avreal,avrall,l-1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
      !$omp parallel default(shared) private(i)
      do i=1,neigmx
        avreal(i)=avrall(i)*dx*dy*dz
      end do
      !$omp end parallel
      if (l-1 > nblas03) then
        call dgemv('n',ncpx*ncpy*ncpz,l-1,-1.0d0,svecre(1,1,1,1,ns,nk),ncpx*ncpy*ncpz &
            ,avreal,1,1.0d0,svecre(1,1,1,l,ns,nk),1)
        call dgemv('n',ncpx*ncpy*ncpz,l-1,-1.0d0,ssvre(1,1,1,1,ns,nk),ncpx*ncpy*ncpz &
            ,avreal,1,1.0d0,ssvre(1,1,1,l,ns,nk),1)
        call dgemv('n',nprjmx*num_ppcell,l-1,-1.0d0,rspsep(1,1,1,ns,nk),nprjmx*num_ppcell &
            ,avreal,1,1.0d0,rspsep(1,1,l,ns,nk),1)
      else
        !$omp parallel default(shared)
        call orthogonalization_r_03(ncpx*ncpy*ncpz,neigmx,l-1,svecre(1,1,1,l,ns,nk),svecre(1,1,1,1,ns,nk),avreal)
        call orthogonalization_r_03(ncpx*ncpy*ncpz,neigmx,l-1,ssvre(1,1,1,l,ns,nk),ssvre(1,1,1,1,ns,nk),avreal)
        call multr(nprjmx*num_ppcell,1,l-1,rspsep(1,1,l,ns,nk),rspsep(1,1,1,ns,nk),avreal)
        !$omp end parallel
      end if
      vecnor=0.0d0
      !$omp parallel default(shared) private(i)
      call orthogonalization_r_05(ncpx*ncpy*ncpz,svecre(1,1,1,l,ns,nk),ssvre(1,1,1,l,ns,nk),vecnor)
      !$omp barrier
      !$omp single
      call mpi_allreduce(vecnor,tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
      vecnor=tmpall*dx*dy*dz
      if (vecnor<=0.0d0) call stopp('orthogonalization_r_01 (3): wf has negative norm')
      tmp=1.0d0/dsqrt(vecnor)
      !$omp end single
      call orthogonalization_r_06(ncpx*ncpy*ncpz,svecre(1,1,1,l,ns,nk),tmp)
      call orthogonalization_r_06(ncpx*ncpy*ncpz,ssvre(1,1,1,l,ns,nk),tmp)
      !$omp do
      do i=1,nprjmx*num_ppcell
        rspsep(i,1,l,ns,nk)=rspsep(i,1,l,ns,nk)*tmp
      end do
      !$omp end parallel
      l=l+1
    end if

  end do !l

  deallocate(avreal,avrall,amatr,amatr_a)
  deallocate(avr)
  return
end subroutine


subroutine orthogonalization_c_01(natom,neigmx,nprjmx,nums,ncol,numk,nk,ns1,l_inp,num_list,num_ppcell,num_spe,nflag,ndisp, & ! <
                                  ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                                                     & ! <
                                  key_ortho_cmpt_innerproduct,                                                             & ! <
                                  key_natpri_in,key_natpri_inps,                                                           & ! <
                                  dx,dy,dz,skpx,skpy,skpz,                                                                 & ! <
                                  nprj,indspe,natpri,naps,natinf,lstvec2,latom,                                            & ! <
                                  lstx,lsty,lstz,natx,naty,natz,                                                           & ! <
                                  vnlocp,sss,                                                                              & ! <
                                  sveccm,ssvcm,cspsep)                                                                       ! X
use mod_stopp
use mod_nonlocaloperation, only:nonlocaloperation_c_01,nonlocaloperation_c_02
use mod_disptime
implicit none
integer, intent(in)::natom,neigmx,nprjmx,nums,ncol,numk,nk,ns1,l_inp,num_list,num_ppcell,num_spe,nflag,ndisp
integer, intent(in)::ncpx,ncpy,ncpz,npxmax,npymax,npzmax
integer, intent(in)::key_ortho_cmpt_innerproduct
integer, intent(in)::key_natpri_in,key_natpri_inps
real*8,  intent(in)::dx,dy,dz,skpx,skpy,skpz
integer, intent(in)::indspe(natom),natpri(natom),naps(natom),nprj(num_spe)
integer, intent(in)::natinf(natom),lstvec2(num_list,num_ppcell),latom(natom)
integer, intent(in)::lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell)
integer, intent(in)::natx(natom),naty(natom),natz(natom)
real*8,  intent(in)::vnlocp(num_list,nprjmx,num_ppcell),sss(nprjmx,nprjmx,num_spe)
complex*16, intent(inout)::sveccm(ncpx,ncpy,ncpz,neigmx,nums,numk)
complex*16, intent(inout)::ssvcm(ncpx,ncpy,ncpz,neigmx,nums,numk)
complex*16, intent(inout)::cspsep(nprjmx,num_ppcell,neigmx,nums,numk)
real*8 tmp,tmpall,vecnor
integer na,iaps,ns2,i,j,l,l0,l_stt,l_end
integer,   allocatable::natdummy(:)
real*8,    allocatable::dijdummy(:,:,:,:) 
complex*16,allocatable::cmat(:,:),cmat_a(:,:)
complex*16,allocatable::avcmplx(:),avcmplxall(:)
complex*16,allocatable::cpsep(:,:,:),ctmp(:,:),ctmpall(:,:)
complex*16,allocatable::vcccm(:,:)
complex*16,allocatable::avc(:,:)
real*8 stime
integer mstep,nblas03,nblas04
complex*16 calph,cbeta

  nblas03=500000
  nblas04=500000
  mstep=32
  if (neigmx .gt. 1024) mstep=64
  if (neigmx .gt. 2048) mstep=128
  if (neigmx .gt. 4096) mstep=256
  calph=dcmplx(1.0d0,0.0d0)
  cbeta=dcmplx(0.0d0,0.0d0)

  allocate(avcmplx(neigmx),avcmplxall(neigmx))
  allocate(cmat(mstep,mstep),cmat_a(mstep,mstep))
  allocate(vcccm(num_list,ncol))
  allocate(avc(num_list,ncol))

  if (l_inp .eq. 0) then
    l_stt=1
    l_end=neigmx
  else
    l_stt=l_inp
    l_end=l_inp
  end if

  call starttime(stime)
  if (nflag .eq. key_ortho_cmpt_innerproduct) then
    allocate(cpsep(nprjmx,num_ppcell,ncol),ctmp(nprjmx*neigmx*ncol,num_ppcell),ctmpall(nprjmx*neigmx*ncol,num_ppcell))
    allocate(natdummy(1),dijdummy(1,1,1,1))
    natdummy=0
    dijdummy=0.0d0
    !$omp parallel default(shared) private(i)
    !$omp do
    do i=1,nprjmx*neigmx*num_ppcell
      ctmpall(i,1)=dcmplx(0.0d0,0.0d0)
    end do
    !$omp end parallel
    do l=l_stt,l_end
      !$omp parallel default(shared) private(i)
      !$omp do
      do i=1,nprjmx*num_ppcell*ncol
        cpsep(i,1,1)=dcmplx(0.0d0,0.0d0)
      end do
      call nonlocaloperation_c_01(natom,nprjmx,num_spe,num_ppcell,num_list,neigmx,ncol,l, & ! <
                                  ncpx,ncpy,ncpz,                                         & ! <
                                  key_natpri_in,key_natpri_inps,                          & ! <
                                  dx,dy,dz,skpx,skpy,skpz,                                & ! <
                                  nprj,indspe,natinf,natpri,naps,lstvec2,                 & ! <
                                  lstx,lsty,lstz,natx,naty,natz,                          & ! <
                                  vnlocp,sveccm(1,1,1,1,ns1,nk),                          & ! <
                                  cpsep,                                                  & ! X
                                  vcccm)                                                    ! W
      !$omp end parallel
      do ns2= 1,ncol
        do i=1,natom
          na=latom(i)
          if ((natpri(na) == key_natpri_in) .or. (natpri(na) == key_natpri_inps)) then
            iaps=naps(na)
            do mpij=1,nprj(indspe(na))
              ctmp(mpij+nprj(indspe(na))*(l-l_stt)+nprj(indspe(na))*(l_end-l_stt+1)*(ns2-1),iaps)=cpsep(mpij,iaps,ns2)
            end do
          end if
        end do
      end do ! ns2
    end do !l
    do i=1,natom
      na=latom(i)
      if ((natpri(na) == key_natpri_in) .or. (natpri(na) == key_natpri_inps)) then
        iaps=naps(na)
        call mpi_allreduce(ctmp(1,iaps),ctmpall(1,iaps),nprj(indspe(na))*(l_end-l_stt+1)*ncol,mpi_double_complex,mpi_sum,mpicom_atom(na),mpij)
      end if
    end do
    do l=l_stt,l_end
      do ns2= 1,ncol
        do i=1,natom
          na=latom(i)
          if ((natpri(na) == key_natpri_in) .or. (natpri(na) == key_natpri_inps)) then
            iaps=naps(na)
            do mpij=1,nprj(indspe(na))
              cspsep(mpij,iaps,l,ns1+ns2-1,nk)=ctmpall(mpij+nprj(indspe(na))*(l-l_stt)+nprj(indspe(na))*neigmx*(ns2-1),iaps)
            end do
          end if
        end do ! i
      end do ! ns2
      !$omp parallel default(shared) private(i)
      do ns2= 1,ncol
        !$omp do
        do i=1,nprjmx*num_ppcell
          cpsep(i,1,ns2)=cspsep(i,1,l,ns1+ns2-1,nk)
        end do
        !$omp do
        do i=1,ncpx*ncpy*ncpz
          ssvcm(i,1,1,l,ns1+ns2-1,nk)=sveccm(i,1,1,l,ns1+ns2-1,nk)
        end do
      end do ! ns2
      call nonlocaloperation_c_02(0,0,natom,num_spe,neigmx,nums,nprjmx,ncol,1,l,num_list,num_ppcell, & ! <
                                  ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                               & ! <
                                  key_natpri_in,key_natpri_inps,0,                                   & ! <
                                  dx,dy,dz,skpx,skpy,skpz,                                           & ! <
                                  nprj,indspe,natinf,lstvec2,natpri,naps,natdummy,                   & ! <
                                  lstx,lsty,lstz,natx,naty,natz,                                     & ! <
                                  vnlocp,cpsep,sss,dijdummy,                                         & ! <
                                  ssvcm(1,1,1,1,ns1,nk),                                             & ! X
                                  vcccm,avc)                                                           ! W
      !$omp end parallel
    end do !l
    deallocate(cpsep,ctmp,ctmpall)
    deallocate(natdummy,dijdummy)
  end if ! (nflag .eq. key_ortho_cmpt_innerproduct)
  if (nflag .eq. key_ortho_cmpt_innerproduct) call endtime(ndisp,stime,'[TI_sub] setup <p|\psi> :')

  l=l_stt
  do while (l < l_end+1)
    if ((l>mstep).and.(l_end-l+1>=mstep)) then
      do j=1,l-1,mstep
        !$omp parallel default(shared) private(i)
        !$omp do
        do i=1,mstep*mstep
          cmat(i,1)=dcmplx(0.0d0,0.0d0)
        end do
        !$omp end parallel
        do ns2= 1,ncol
          call zgemm('c','n',mstep,mstep,ncpx*ncpy*ncpz,calph,sveccm(1,1,1,l,ns1+ns2-1,nk), &
                     ncpx*ncpy*ncpz,ssvcm(1,1,1,j,ns1+ns2-1,nk),ncpx*ncpy*ncpz,calph,cmat,mstep)
        end do
        call mpi_allreduce(cmat,cmat_a,mstep*mstep,mpi_double_complex,mpi_sum,mpicom_space,mpij)
        !$omp parallel default(shared) private(i,ns2)
        !$omp do
        do i=1,mstep*mstep
          cmat(i,1)=dconjg(cmat_a(i,1))*dx*dy*dz
        end do
        !$omp end parallel
        do ns2= 1,ncol
          call zgemm('n','t',ncpx*ncpy*ncpz,mstep,mstep,-calph,sveccm(1,1,1,j,ns1+ns2-1,nk), &
                     ncpx*ncpy*ncpz,cmat,mstep,calph,sveccm(1,1,1,l,ns1+ns2-1,nk),ncpx*ncpy*ncpz)
          call zgemm('n','t',ncpx*ncpy*ncpz,mstep,mstep,-calph,ssvcm( 1,1,1,j,ns1+ns2-1,nk), &
                     ncpx*ncpy*ncpz,cmat,mstep,calph,ssvcm( 1,1,1,l,ns1+ns2-1,nk),ncpx*ncpy*ncpz)
        end do
        !$omp parallel default(shared) private(ns2)
        do ns2= 1,ncol
          call multc(nprjmx*num_ppcell,mstep,mstep,cspsep(1,1,l,ns1+ns2-1,nk),cspsep(1,1,j,ns1+ns2-1,nk),cmat)
        end do
        !$omp end parallel
      end do
      vecnor=0.0d0
      !$omp parallel default(shared) private(i)
      call orthogonalization_c_05(ncpx*ncpy*ncpz,neigmx,neigmx,ncol,l,sveccm(1,1,1,1,ns1,nk),l,ssvcm(1,1,1,1,ns1,nk),vecnor)
      !$omp barrier
      !$omp single
      call mpi_allreduce(vecnor,tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
      vecnor=tmpall*dx*dy*dz
      if (vecnor<=0.0d0) call stopp('orthogonalization_c_01 (1): wf has negative norm')
      tmp=1.0d0/dsqrt(vecnor)
      !$omp end single
      do ns2= ns1,ns1+ncol-1
        call orthogonalization_c_06(ncpx*ncpy*ncpz,sveccm(1,1,1,l,ns2,nk),tmp)
        call orthogonalization_c_06(ncpx*ncpy*ncpz,ssvcm(1,1,1,l,ns2,nk),tmp)
      end do
      do ns2= ns1,ns1+ncol-1
        !$omp do
        do i=1,nprjmx*num_ppcell
          cspsep(i,1,l,ns2,nk)=cspsep(i,1,l,ns2,nk)*tmp
        end do
        !$omp end do nowait
      end do
      !$omp end parallel
      do l0=l+1,l+mstep-1
        !$omp parallel default(shared) private(i)
        !$omp do
        do i=1,neigmx
          avcmplx(i)=dcmplx(0.0d0,0.0d0)
        end do
        !$omp end parallel
        if (l0-1 > nblas04) then
          do ns2=1,ncol
            call zgemv('c',ncpx*ncpy*ncpz,l0-l,calph,ssvcm(1,1,1,l,ns1+ns2-1,nk),ncpx*ncpy*ncpz &
                ,sveccm(1,1,1,l0,ns1+ns2-1,nk),1,cbeta,avcmplx,1)
          end do
        else
        !$omp parallel default(shared)
          call orthogonalization_c_04(ncpx*ncpy*ncpz,neigmx,neigmx,ncol,l0,sveccm(1,1,1,1,ns1,nk) &
                          ,l,l0-1,ssvcm(1,1,1,1,ns1,nk),avcmplx)
        !$omp end parallel
        end if
        call mpi_allreduce(avcmplx,avcmplxall,(l0-l),mpi_double_complex,mpi_sum,mpicom_space,mpij)
        !$omp parallel default(shared) private(i)
        !$omp do
        do i=1,l0-l
          avcmplx(i)=avcmplxall(i)*dx*dy*dz
        end do
        !$omp end parallel
        do ns2= ns1,ns1+ncol-1
          if (l0-l > nblas03) then
            call zgemv('n',ncpx*ncpy*ncpz,l0-l,-calph,sveccm(1,1,1,l,ns2,nk),ncpx*ncpy*ncpz &
                ,avcmplx,1,calph,sveccm(1,1,1,l0,ns2,nk),1)
            call zgemv('n',ncpx*ncpy*ncpz,l0-l,-calph,ssvcm(1,1,1,l,ns2,nk),ncpx*ncpy*ncpz &
                ,avcmplx,1,calph,ssvcm(1,1,1,l0,ns2,nk),1)
            call zgemv('n',nprjmx*num_ppcell,l0-l,-calph,cspsep(1,1,l,ns2,nk),nprjmx*num_ppcell &
                ,avcmplx,1,calph,cspsep(1,1,l0,ns2,nk),1)
          else
            !$omp parallel default(shared)
            call orthogonalization_c_03(ncpx*ncpy*ncpz,l0-l,l0-l,1,sveccm(1,1,1,l0,ns2,nk),sveccm(1,1,1,l,ns2,nk),avcmplx)
            call orthogonalization_c_03(ncpx*ncpy*ncpz,l0-l,l0-l,1,ssvcm(1,1,1,l0,ns2,nk),ssvcm(1,1,1,l,ns2,nk),avcmplx)
            call multc(nprjmx*num_ppcell,1,l0-l,cspsep(1,1,l0,ns2,nk),cspsep(1,1,l,ns2,nk),avcmplx)
            !$omp end parallel
          end if
        end do
        vecnor=0.0d0
        !$omp parallel default(shared) private(i)
        call orthogonalization_c_05(ncpx*ncpy*ncpz,neigmx,neigmx,ncol,l0,sveccm(1,1,1,1,ns1,nk),l0,ssvcm(1,1,1,1,ns1,nk),vecnor)
        !$omp barrier
        !$omp single
        tmpall=0.0d0
        call mpi_allreduce(vecnor,tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
        vecnor=tmpall*dx*dy*dz
        if (vecnor<=0.0d0) call stopp('orthogonalization_c_01 (2): wf has negative norm')
        tmp=1.0d0/dsqrt(vecnor)
        !$omp end single
        do ns2= ns1,ns1+ncol-1
          call orthogonalization_c_06(ncpx*ncpy*ncpz,sveccm(1,1,1,l0,ns2,nk),tmp)
          call orthogonalization_c_06(ncpx*ncpy*ncpz,ssvcm(1,1,1,l0,ns2,nk),tmp)
        end do
        do ns2= ns1,ns1+ncol-1
          !$omp do
          do i=1,nprjmx*num_ppcell
            cspsep(i,1,l0,ns2,nk)=cspsep(i,1,l0,ns2,nk)*tmp
          end do
        end do
        !$omp end parallel
      end do
      l=l+mstep
    else
      !$omp parallel default(shared) private(i)
      !$omp do
      do i=1,neigmx
        avcmplx(i)=dcmplx(0.0d0,0.0d0)
      end do
      !$omp end parallel
      if (l > nblas04) then
        do ns2=1,ncol
          call zgemv('c',ncpx*ncpy*ncpz,l-1,calph,ssvcm(1,1,1,1,ns1+ns2-1,nk),ncpx*ncpy*ncpz &
              ,sveccm(1,1,1,l,ns1+ns2-1,nk),1,cbeta,avcmplx,1)
        end do
      else
      !$omp parallel default(shared)
        call orthogonalization_c_04(ncpx*ncpy*ncpz,neigmx,neigmx,ncol,l,sveccm(1,1,1,1,ns1,nk) &
                        ,1,l-1,ssvcm(1,1,1,1,ns1,nk),avcmplx)
      !$omp end parallel
      end if
      call mpi_allreduce(avcmplx,avcmplxall,neigmx,mpi_double_complex,mpi_sum,mpicom_space,mpij)
      !$omp parallel default(shared) private(i)
      !$omp do
      do i=1,neigmx
        avcmplx(i)=avcmplxall(i)*dx*dy*dz
      end do
      !$omp end parallel
      do ns2= ns1,ns1+ncol-1
        ! here the noco spin-loop is outside 'orthogonalization_c_03', as input and output arrays ('svec') should not overlap in memory
        if (l-1 > nblas03) then
          call zgemv('n',ncpx*ncpy*ncpz,l-1,-calph,sveccm(1,1,1,1,ns2,nk),ncpx*ncpy*ncpz &
              ,avcmplx,1,calph,sveccm(1,1,1,l,ns2,nk),1)
          call zgemv('n',ncpx*ncpy*ncpz,l-1,-calph,ssvcm(1,1,1,1,ns2,nk),ncpx*ncpy*ncpz &
              ,avcmplx,1,calph,ssvcm(1,1,1,l,ns2,nk),1)
          call zgemv('n',nprjmx*num_ppcell,l-1,-calph,cspsep(1,1,1,ns2,nk),nprjmx*num_ppcell &
              ,avcmplx,1,calph,cspsep(1,1,l,ns2,nk),1)
        else
          !$omp parallel default(shared)
          call orthogonalization_c_03(ncpx*ncpy*ncpz,l-1,l-1,1,sveccm(1,1,1,l,ns2,nk),sveccm(1,1,1,1,ns2,nk),avcmplx(1))
          call orthogonalization_c_03(ncpx*ncpy*ncpz,l-1,l-1,1,ssvcm(1,1,1,l,ns2,nk),ssvcm(1,1,1,1,ns2,nk),avcmplx(1))
          call multc(nprjmx*num_ppcell,1,l-1,cspsep(1,1,l,ns2,nk),cspsep(1,1,1,ns2,nk),avcmplx)
          !$omp end parallel
        end if
      end do
      !$omp parallel default(shared)
      !$omp barrier
      !$omp single
      vecnor=0.0d0
      !$omp end single
      call orthogonalization_c_05(ncpx*ncpy*ncpz,neigmx,neigmx,ncol,l,sveccm(1,1,1,1,ns1,nk) &
                      ,l,ssvcm(1,1,1,1,ns1,nk),vecnor)
      !$omp barrier
      !$omp single
      call mpi_allreduce(vecnor,tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
      vecnor=tmpall*dx*dy*dz
      if (vecnor<=0.0d0) call stopp('orthogonalization_c_01 (3): wf has negative norm')
      tmp=1.0d0/dsqrt(vecnor)
      !$omp end single
      do ns2= ns1,ns1+ncol-1
        call orthogonalization_c_06(ncpx*ncpy*ncpz,sveccm(1,1,1,l,ns2,nk),tmp)
        call orthogonalization_c_06(ncpx*ncpy*ncpz,ssvcm(1,1,1,l,ns2,nk),tmp)
      end do
      do ns2= ns1,ns1+ncol-1
        !$omp do
        do i=1,nprjmx*num_ppcell
          cspsep(i,1,l,ns2,nk)=cspsep(i,1,l,ns2,nk)*tmp
        end do
      end do
      !$omp end parallel
      l=l+1
    end if

  end do !l

  deallocate(cmat,cmat_a)
  deallocate(avcmplx,avcmplxall)
  deallocate(vcccm)
  deallocate(avc)

  return
end subroutine orthogonalization_c_01


subroutine orthogonalization_r_03(nmax,neigmx,lmax,vre,avre,ar)
implicit none
integer, intent(in)    :: nmax,neigmx,lmax
real*8,  intent(in)    :: avre(nmax,neigmx)
real*8,  intent(in)    :: ar(neigmx)
real*8,  intent(inout) :: vre(nmax)
integer i,j

  do j=1,lmax
!$omp do
!ocl norecurrence(vre)
    do i=1,nmax
      vre(i)=vre(i)-ar(j)*avre(i,j)
    end do
!$omp end do nowait
  end do
  return
end subroutine


subroutine orthogonalization_c_03(nmax,neigmx,lmax,ncol,vcm,avcm,ca)
implicit none
integer, intent(in)     :: nmax,neigmx,ncol,lmax
complex*16,intent(in)   :: avcm(nmax,neigmx,ncol)
complex*16,intent(in)   :: ca(neigmx)
complex*16,intent(inout):: vcm(nmax,ncol)
integer ns,i,j

  do ns=1,ncol
    do j=1,lmax
!$omp do
!ocl norecurrence(vcm)
      do i=1,nmax
        vcm(i,ns)=vcm(i,ns)-ca(j)*avcm(i,j,ns)
      end do
!$omp end do nowait
    end do
  end do
  return
end subroutine


subroutine orthogonalization_r_04(nmax,neigmx1,neigmx2,l1,v1re,l2s,l2e,v2re,br)
implicit none
integer, intent(in)    :: nmax,neigmx1,neigmx2,l1,l2s,l2e
real*8,  intent(in)    :: v1re(nmax,neigmx1),v2re(nmax,neigmx2)
real*8,  intent(inout) :: br(neigmx2)
integer i,l
real*8, allocatable::ar(:)

  allocate(ar(neigmx2))
  ar=0.0d0
  do l=l2s,l2e
!$omp do
    do i=1,nmax
      ar(l-l2s+1)=ar(l-l2s+1)+v2re(i,l)*v1re(i,l1)
    end do
!$omp end do nowait
  end do
!$omp critical
  do l=l2s,l2e
    br(l-l2s+1)=br(l-l2s+1)+ar(l-l2s+1)
  end do
!$omp end critical
  deallocate(ar)
  return
end subroutine


subroutine orthogonalization_c_04(nmax,neigmx1,neigmx2,ncol,l1,v1cm,l2s,l2e,v2cm,cb)
implicit none
integer, intent(in)      :: nmax,neigmx1,neigmx2,ncol,l1,l2s,l2e
complex*16,intent(in)    :: v1cm(nmax,neigmx1,ncol),v2cm(nmax,neigmx2,ncol)
complex*16,intent(inout) :: cb(neigmx2)
integer i,ns,l
complex*16,allocatable::ca(:)

  allocate(ca(neigmx2))
  ca=dcmplx(0.0d0,0.0d0)
  do ns=1,ncol
    do l=l2s,l2e
!$omp do
      do i=1,nmax
        ca(l-l2s+1)=ca(l-l2s+1)+dconjg(v2cm(i,l,ns))*v1cm(i,l1,ns)
      end do
!$omp end do nowait
    end do
  end do
!$omp critical
  do l=l2s,l2e
    cb(l-l2s+1)=cb(l-l2s+1)+ca(l-l2s+1)
  end do
!$omp end critical
  deallocate(ca)
  return
end subroutine


subroutine orthogonalization_r_05(nmax,v1re,v2re,b)
implicit none
integer nmax,i
real*8 v1re(nmax),v2re(nmax)
real*8 a,b
  a=0.0d0
!$omp do
  do i=1,nmax
    a=a+v1re(i)*v2re(i)
  end do
!$omo end do nowait
!$omp critical
  b=b+a
!$omp end critical
  return
end subroutine


subroutine orthogonalization_c_05(nmax,neigmx1,neigmx2,ncol,l1,v1cm,l2,v2cm,b)
implicit none
integer,   intent(in)    :: nmax,neigmx1,neigmx2,ncol,l1,l2
complex*16,intent(in)    :: v1cm(nmax,neigmx1,ncol),v2cm(nmax,neigmx2,ncol)
real*8,    intent(inout) :: b
integer i,ns
real*8  a
  a=0.0d0
  do ns=1,ncol
!$omp do
    do i=1,nmax
      a=a+dreal(v1cm(i,l1,ns))*dreal(v2cm(i,l2,ns))+dimag(v1cm(i,l1,ns))*dimag(v2cm(i,l2,ns))
    end do
  end do
!$omo end do nowait
!$omp critical
  b=b+a
!$omp end critical
  return
end subroutine


subroutine orthogonalization_r_06(nmax,v,a)
implicit none
integer nmax,i
real*8 v(nmax)
real*8 a
!cdir nodep
!$omp do
  do i=1,nmax
     v(i)=v(i)*a
  end do
!$omp end do nowait
  return
end subroutine


subroutine orthogonalization_c_06(nmax,v,a)
implicit none
integer nmax,i
complex*16 v(nmax)
real*8 a
!cdir nodep
!$omp do
  do i=1,nmax
     v(i)=v(i)*a
  end do
!$omp end do nowait
  return
end subroutine


subroutine multr(n,m,l,a0,a1,b)
implicit none
integer, intent(in)::n,m,l
real*8, intent(inout)::a0(n,m)
real*8, intent(in)::a1(n,m)
real*8, intent(in)::b(m,m)
integer i,j,k
  do k=1,l
!$omp do
    do j=1,m
      do i=1,n
        a0(i,j)=a0(i,j)-b(j,k)*a1(i,k)
      end do
    end do
!$omp end do nowait
  end do
  return
end subroutine


subroutine multc(n,m,l,ca0,ca1,cb)
implicit none
integer, intent(in)::n,m,l
complex*16, intent(inout)::ca0(n,m)
complex*16, intent(in)::ca1(n,l)
complex*16, intent(in)::cb(m,l)
integer i,j,k
  do k=1,l
!$omp do
    do j=1,m
      do i=1,n
        ca0(i,j)=ca0(i,j)-cb(j,k)*ca1(i,k)
      end do
    end do
!$omp end do nowait
  end do
  return
end subroutine

end module
