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
! **********  force8f.F90 04/27/2023-01  **********

module mod_force
implicit none
contains

! ========================================================================================

subroutine force_eig_r( &
 nperi,numk,nums,neigmx,nf,ncpx,ncpy,ncpz,num_list,natom,num_spe,num_ppcell,nprjmx, & ! <
 key_natpri_in,key_natpri_inps,                                                     & ! <
 lstvec2,latom,lstx,lsty,lstz,nprj,natpri,naps,indspe,natinf,                       & ! <
 xmax,ymax,zmax,                                                                    & ! <
 sss,vnlocp,dij,svecre,sval,fnele,                                                  & ! <
 fatx,faty,fatz)                                                                      ! >
use mod_mpi
use mod_overlap_finitedifference, only:overlap_finitedifference_init,overlap_finitedifference_final &
                                      ,overlap_finitedifference_r,overlap_fdcheck_r
use mod_nonlocaloperation, only:nonlocaloperation_r_01
implicit none
integer, intent(in) ::nperi,numk,nums,neigmx
integer, intent(in) ::nf,ncpx,ncpy,ncpz,num_list
integer, intent(in) ::natom,num_spe,num_ppcell,nprjmx
integer, intent(in) ::key_natpri_in,key_natpri_inps
integer, intent(in) ::lstvec2(num_list,num_ppcell),latom(natom)
integer, intent(in) ::lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell)
integer, intent(in) ::nprj(num_spe)
integer, intent(in) ::natpri(natom),naps(natom),indspe(natom),natinf(natom)
real*8,  intent(in) ::xmax,ymax,zmax
real*8,  intent(in) ::sss(nprjmx,nprjmx,num_spe)
real*8,  intent(in) ::vnlocp(num_list,nprjmx,num_ppcell)
real*8,  intent(in) ::dij(nprjmx,nprjmx,nums,natom)
real*8,  intent(in) ::svecre(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8,  intent(in) ::sval(neigmx,nums,numk),fnele(neigmx,nums,numk)
real*8,  intent(out)::fatx(natom),faty(natom),fatz(natom)
integer ::na,nk,ns,l,lmn,j,i,ix,iy,iz,iaps
real*8  ::dx,dy,dz, fatall
real*8,allocatable::vre2(:,:,:)
real*8,allocatable::avre(:,:,:)
real*8,allocatable::vre(:,:,:)
real*8,allocatable::rpsep(:,:),rdpsep(:,:)
real*8,allocatable::rpsepeig(:,:),rdpsepeig(:,:),rpsepall(:,:)
real*8,allocatable::dijeig(:,:,:,:)
  call overlap_finitedifference_init(ncpx,ncpy,ncpz,nf,1,0)

  dx=2.0d0*xmax/(ncpx*nprocx)
  dy=2.0d0*ymax/(ncpy*nprocy)
  dz=2.0d0*zmax/(ncpz*nprocz)
  allocate(vre2(ncpx,ncpy,ncpz))
  allocate(avre(ncpx,ncpy,ncpz))
  allocate(vre(1-nf:ncpx+nf,1-nf:ncpy+nf,1-nf:ncpz+nf))
  allocate(rpsep(nprjmx,num_ppcell),rdpsep(nprjmx,num_ppcell))
  allocate(rpsepeig(nprjmx*neigmx,num_ppcell),rdpsepeig(nprjmx*neigmx,num_ppcell),rpsepall(nprjmx*neigmx,num_ppcell))
  allocate(dijeig(nprjmx,nprjmx,nums,natom))

! ----------  zero clear  ----------
  fatx=0.0d0
  faty=0.0d0
  fatz=0.0d0
! ----------------------------------

! ----------  nonlocal part of p.p. - valence electron interaction  ----------
  do nk=1,numk
    do ns=1,nums
      do l=1,neigmx
!$omp parallel default(shared) private(i)
!$omp do
        do i=1,ncpx*ncpy*ncpz
          vre2(i,1,1)=svecre(i,1,1,l,ns,nk)
        end do
!$omp do
        do i=1,nprjmx*num_ppcell
          rpsep(i,1)=0.0d0
        end do
        call nonlocaloperation_r_01(natom,nprjmx,num_spe,num_ppcell,num_list,     & ! <
                                    ncpx,ncpy,ncpz,                               & ! <
                                    key_natpri_in,key_natpri_inps,                & ! <
                                    nprj,indspe,natinf,natpri,naps,lstvec2,       & ! <
                                    vnlocp,vre2,                                  & ! <
                                    rpsep)                                          ! X
!$omp end parallel
        do i=1,natom
          na=latom(i)
          if ((natpri(na) .eq. key_natpri_in) .or. (natpri(na) .eq. key_natpri_inps)) then
            iaps=naps(na)
            rpsepeig(nprj(indspe(na))*(l-1)+1:nprj(indspe(na))*l,iaps)=rpsep(1:nprj(indspe(na)),iaps)
          end if
        end do
      end do ! l
      do i=1,natom
        na=latom(i)
        if ((natpri(na) .eq. key_natpri_in) .or. (natpri(na) .eq. key_natpri_inps)) then
          iaps=naps(na)
          call mpi_allreduce(rpsepeig(1,iaps),rpsepall(1,iaps),nprj(indspe(na))*neigmx &
                            ,mpi_double_precision,mpi_sum,mpicom_atom(na),mpij)
          rpsepeig(1:nprj(indspe(na))*neigmx,iaps)=rpsepall(1:nprj(indspe(na))*neigmx,iaps)
        end if
      end do

      do lmn=1,3
        do l=1,neigmx
!$omp parallel default(shared) private(i)
!$omp do
          do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)
            vre(-nf+i,-nf+1,-nf+1)=0.0d0
          end do
!$omp do
          do iz=1,ncpz
          do iy=1,ncpy
          do ix=1,ncpx
            vre(ix,iy,iz)=svecre(ix,iy,iz,l,ns,nk)
          end do
          end do
          end do
!$omp end parallel
          call overlap_finitedifference_r(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,vre)
          call overlap_fdcheck_r(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,vre)
!$omp parallel default(shared)
          if (lmn==1) call force_r_02(ncpx,ncpy,ncpz,nf,dx,vre,avre)
          if (lmn==2) call force_r_03(ncpx,ncpy,ncpz,nf,dy,vre,avre)
          if (lmn==3) call force_r_04(ncpx,ncpy,ncpz,nf,dz,vre,avre)
!$omp do
          do i=1,nprjmx*num_ppcell
            rdpsep(i,1)=0.0d0
          end do
          call nonlocaloperation_r_01(natom,nprjmx,num_spe,num_ppcell,num_list,     & ! <
                                      ncpx,ncpy,ncpz,                               & ! <
                                      key_natpri_in,key_natpri_inps,                & ! <
                                      nprj,indspe,natinf,natpri,naps,lstvec2,       & ! <
                                      vnlocp,avre,                                  & ! <
                                      rdpsep)                                         ! X
!$omp end parallel
          do i=1,natom
            na=latom(i)
            if ((natpri(na) .eq. key_natpri_in) .or. (natpri(na) .eq. key_natpri_inps)) then
              iaps=naps(na)
              rdpsepeig(nprj(indspe(na))*(l-1)+1:nprj(indspe(na))*l,iaps)=rdpsep(1:nprj(indspe(na)),iaps)
            end if
          end do
        end do ! l
        do i=1,natom
          na=latom(i)
          if ((natpri(na) .eq. key_natpri_in) .or. (natpri(na) .eq. key_natpri_inps)) then
            iaps=naps(na)
            call mpi_allreduce(rdpsepeig(1,iaps),rpsepall(1,iaps),nprj(indspe(na))*neigmx &
                              ,mpi_double_precision,mpi_sum,mpicom_atom(na),mpij)
            rdpsepeig(1:nprj(indspe(na))*neigmx,iaps)=rpsepall(1:nprj(indspe(na))*neigmx,iaps)
          end if
        end do

        do l=1,neigmx
          do na=1,natom
            if (natpri(na) .eq. key_natpri_in) then
              iaps=naps(na)
              rpsep(1:nprj(indspe(na)),iaps)=rpsepeig(nprj(indspe(na))*(l-1)+1:nprj(indspe(na))*l,iaps)
              rdpsep(1:nprj(indspe(na)),iaps)=rdpsepeig(nprj(indspe(na))*(l-1)+1:nprj(indspe(na))*l,iaps)
              dijeig(:,:,ns,na)=dij(:,:,ns,na)-sval(l,ns,nk)*sss(:,:,indspe(na))
              do j=1,nprj(indspe(na))
                do i=1,nprj(indspe(na))
                  if (lmn .eq. 1) fatx(na)=fatx(na)-2.0d0*fnele(l,ns,nk)*(rpsep(i,iaps)*rdpsep(j,iaps))*dijeig(i,j,ns,na)
                  if (lmn .eq. 2) faty(na)=faty(na)-2.0d0*fnele(l,ns,nk)*(rpsep(i,iaps)*rdpsep(j,iaps))*dijeig(i,j,ns,na)
                  if (lmn .eq. 3) fatz(na)=fatz(na)-2.0d0*fnele(l,ns,nk)*(rpsep(i,iaps)*rdpsep(j,iaps))*dijeig(i,j,ns,na)
                end do
              end do
            end if
          end do
        end do ! l

      end do ! lmn
    end do ! ns
  end do ! nk
! ----------------------------------------------------------------------------

  call overlap_finitedifference_final

  do na=1,natom
    if (natpri(na) .eq. key_natpri_in) then
      iaps=naps(na)
      fatall= fatx(na)
      call mpi_reduce(fatall,fatx(na),1,mpi_double_precision,mpi_sum,0,mpicom_kpt,mpij)
      fatall= faty(na)
      call mpi_reduce(fatall,faty(na),1,mpi_double_precision,mpi_sum,0,mpicom_kpt,mpij)
      fatall= fatz(na)
      call mpi_reduce(fatall,fatz(na),1,mpi_double_precision,mpi_sum,0,mpicom_kpt,mpij)
    end if
  end do

  deallocate(vre2)
  deallocate(avre)
  deallocate(vre)
  deallocate(rpsep,rdpsep)
  deallocate(rpsepeig,rdpsepeig,rpsepall)
  deallocate(dijeig)

end subroutine force_eig_r


subroutine force_eig_c( &
 nperi,numk,nums,ncol,neigmx,nf,ncpx,ncpy,ncpz,num_list,natom,num_spe,num_ppcell,nprjmx, & ! <
 key_natpri_in,key_natpri_inps,                                                          & ! <
 lstvec2,latom,lstx,lsty,lstz,nprj,natpri,naps,indspe,natinf,                            & ! <
 xmax,ymax,zmax,sss,vnlocp,dij,                                                          & ! <
 natx,naty,natz,skpxx,skpyy,skpzz,                                                       & ! <
 sveccm,sval,fnele,                                                                      & ! <
 fatx,faty,fatz)                                                                           ! >
use mod_mpi
use mod_overlap_finitedifference, only:overlap_finitedifference_init,overlap_finitedifference_final &
                                      ,overlap_finitedifference_c,overlap_fdcheck_c
use mod_nonlocaloperation, only:nonlocaloperation_c_01
implicit none
integer,   intent(in) ::nperi,numk,nums,ncol,neigmx
integer,   intent(in) ::nf,ncpx,ncpy,ncpz,num_list
integer,   intent(in) ::natom,num_spe,num_ppcell,nprjmx
integer,   intent(in) ::key_natpri_in,key_natpri_inps
integer,   intent(in) ::lstvec2(num_list,num_ppcell),latom(natom)
integer,   intent(in) ::lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell)
integer,   intent(in) ::nprj(num_spe)
integer,   intent(in) ::natpri(natom),naps(natom),indspe(natom),natinf(natom)
real*8,    intent(in) ::xmax,ymax,zmax
real*8,    intent(in) ::sss(nprjmx,nprjmx,num_spe)
real*8,    intent(in) ::vnlocp(num_list,nprjmx,num_ppcell)
real*8,    intent(in) ::dij(nprjmx,nprjmx,nums,natom)
integer,   intent(in) ::natx(natom),naty(natom),natz(natom)
real*8,    intent(in) ::skpxx(numk),skpyy(numk),skpzz(numk)
complex*16,intent(in) ::sveccm(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8,    intent(in) ::sval(neigmx,nums,numk),fnele(neigmx,nums,numk)
real*8,    intent(out)::fatx(natom),faty(natom),fatz(natom)
integer               ::na,nk,ns,l,lmn,j,i,ix,iy,iz,iaps
real*8                ::dx,dy,dz, fatall
real*8                ::skpx,skpy,skpz
complex*16,allocatable::vcm2(:,:,:)
complex*16,allocatable::avcm(:,:,:)
complex*16,allocatable::vcm(:,:,:)
complex*16,allocatable::vcccm(:,:)
complex*16,allocatable::cpsep(:,:),cdpsep(:,:)
complex*16,allocatable::cpsepeig(:,:),cdpsepeig(:,:),cpsepall(:,:)
real*8,    allocatable::dijeig(:,:,:,:)
complex*16 cuniti
  cuniti=dcmplx(0.0d0,1.0d0)

  call overlap_finitedifference_init(ncpx,ncpy,ncpz,nf,ncol,1)

  dx=2.0d0*xmax/(ncpx*nprocx)
  dy=2.0d0*ymax/(ncpy*nprocy)
  dz=2.0d0*zmax/(ncpz*nprocz)
  allocate(vcm2(ncpx,ncpy,ncpz))
  allocate(avcm(ncpx,ncpy,ncpz))
  allocate(vcm(1-nf:ncpx+nf,1-nf:ncpy+nf,1-nf:ncpz+nf))
  allocate(vcccm(num_list,ncol))
  allocate(cpsep(nprjmx,num_ppcell),cdpsep(nprjmx,num_ppcell))
  allocate(cpsepeig(nprjmx*neigmx,num_ppcell),cdpsepeig(nprjmx*neigmx,num_ppcell),cpsepall(nprjmx*neigmx,num_ppcell))
  allocate(dijeig(nprjmx,nprjmx,nums*ncol,natom))

! ----------  zero clear  ----------
  fatx=0.0d0
  faty=0.0d0
  fatz=0.0d0
! ----------------------------------

! ----------  nonlocal part of p.p. - valence electron interaction  ----------
  do nk=1,numk
    do ns=1,nums
      skpx=skpxx(nk)
      skpy=skpyy(nk)
      skpz=skpzz(nk)
      do l=1,neigmx
!$omp parallel default(shared) private(i)
!$omp do
        do i=1,ncpx*ncpy*ncpz
          vcm2(i,1,1)=sveccm(i,1,1,l,ns,nk)
        end do
!$omp do
        do i=1,nprjmx*num_ppcell
          cpsep(i,1)=dcmplx(0.0d0,0.0d0)
        end do
        call nonlocaloperation_c_01(natom,nprjmx,num_spe,num_ppcell,num_list,1,ncol,1, & ! <
                                    ncpx,ncpy,ncpz,                                    & ! <
                                    key_natpri_in,key_natpri_inps,                     & ! <
                                    dx,dy,dz,skpx,skpy,skpz,                           & ! <
                                    nprj,indspe,natinf,natpri,naps,lstvec2,            & ! <
                                    lstx,lsty,lstz,natx,naty,natz,                     & ! <
                                    vnlocp,vcm2,                                       & ! <
                                    cpsep,                                             & ! X
                                    vcccm)                                               ! W
!$omp end parallel
        do i=1,natom
          na=latom(i)
          if ((natpri(na) .eq. key_natpri_in) .or. (natpri(na) .eq. key_natpri_inps)) then
            iaps=naps(na)
            cpsepeig(nprj(indspe(na))*(l-1)+1:nprj(indspe(na))*l,iaps)=cpsep(1:nprj(indspe(na)),iaps)
          end if
        end do
      end do
      do i=1,natom
        na=latom(i)
        if ((natpri(na) .eq. key_natpri_in) .or. (natpri(na) .eq. key_natpri_inps)) then
          iaps=naps(na)
          call mpi_allreduce(cpsepeig(1,iaps),cpsepall(1,iaps),nprj(indspe(na))*neigmx &
                            ,mpi_double_complex,mpi_sum,mpicom_atom(na),mpij)
          cpsepeig(1:nprj(indspe(na))*neigmx,iaps)=cpsepall(1:nprj(indspe(na))*neigmx,iaps)
        end if
      end do

      do lmn=1,3
        do l=1,neigmx
!$omp parallel default(shared) private(i,ix,iy,iz)
!$omp do
          do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)
            vcm(-nf+i,-nf+1,-nf+1)=dcmplx(0.0d0,0.0d0)
          end do
!$omp do
          do iz=1,ncpz
          do iy=1,ncpy
          do ix=1,ncpx
            vcm(ix,iy,iz)=sveccm(ix,iy,iz,l,ns,nk)
          end do
          end do
          end do
!$omp end parallel
          call overlap_finitedifference_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,ncol,vcm)
          call overlap_fdcheck_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,ncol,vcm)
!$omp parallel default(shared)
          if (lmn==1) call force_c_02(ncpx,ncpy,ncpz,nf,dx,vcm,avcm)
          if (lmn==2) call force_c_03(ncpx,ncpy,ncpz,nf,dy,vcm,avcm)
          if (lmn==3) call force_c_04(ncpx,ncpy,ncpz,nf,dz,vcm,avcm)
!$omp do
          do i=1,nprjmx*num_ppcell
            cdpsep(i,1)=dcmplx(0.0d0,0.0d0)
          end do
          call nonlocaloperation_c_01(natom,nprjmx,num_spe,num_ppcell,num_list,1,ncol,1, & ! <
                                      ncpx,ncpy,ncpz,                                    & ! <
                                      key_natpri_in,key_natpri_inps,                     & ! <
                                      dx,dy,dz,skpx,skpy,skpz,                           & ! <
                                      nprj,indspe,natinf,natpri,naps,lstvec2,            & ! <
                                      lstx,lsty,lstz,natx,naty,natz,                     & ! <
                                      vnlocp,avcm,                                       & ! <
                                      cdpsep,                                            & ! X
                                      vcccm)                                               ! W
!$omp end parallel
          do i=1,natom
            na=latom(i)
            if ((natpri(na) .eq. key_natpri_in) .or. (natpri(na) .eq. key_natpri_inps)) then
              iaps=naps(na)
              cdpsepeig(nprj(indspe(na))*(l-1)+1:nprj(indspe(na))*l,iaps)=cdpsep(1:nprj(indspe(na)),iaps)
            end if
          end do
        end do ! l
        do i=1,natom
          na=latom(i)
          if ((natpri(na) .eq. key_natpri_in) .or. (natpri(na) .eq. key_natpri_inps)) then
            iaps=naps(na)
            call mpi_allreduce(cdpsepeig(1,iaps),cpsepall(1,iaps),nprj(indspe(na))*neigmx &
                              ,mpi_double_complex,mpi_sum,mpicom_atom(na),mpij)
            cdpsepeig(1:nprj(indspe(na))*neigmx,iaps)=cpsepall(1:nprj(indspe(na))*neigmx,iaps)
          end if
        end do

        do l=1,neigmx
          do na=1,natom
            if (natpri(na) .eq. key_natpri_in) then
              iaps=naps(na)
              cpsep(1:nprj(indspe(na)),iaps)=cpsepeig(nprj(indspe(na))*(l-1)+1:nprj(indspe(na))*l,iaps)
              cdpsep(1:nprj(indspe(na)),iaps)=cdpsepeig(nprj(indspe(na))*(l-1)+1:nprj(indspe(na))*l,iaps)
              dijeig(:,:,ns,na)=dij(:,:,ns,na)-sval(l,ns,nk)*sss(:,:,indspe(na))
              do j=1,nprj(indspe(na))
                do i=1,nprj(indspe(na))
                  if (lmn .eq. 1) fatx(na)=fatx(na)-2.0d0*fnele(l,ns,nk) &
                    *dreal(cpsep(i,iaps)*dconjg(cuniti*skpx*cpsep(i,iaps)+cdpsep(j,iaps)))*dijeig(i,j,ns,na)
                  if (lmn .eq. 2) faty(na)=faty(na)-2.0d0*fnele(l,ns,nk) &
                    *dreal(cpsep(i,iaps)*dconjg(cuniti*skpy*cpsep(i,iaps)+cdpsep(j,iaps)))*dijeig(i,j,ns,na)
                  if (lmn .eq. 3) fatz(na)=fatz(na)-2.0d0*fnele(l,ns,nk) &
                    *dreal(cpsep(i,iaps)*dconjg(cuniti*skpz*cpsep(i,iaps)+cdpsep(j,iaps)))*dijeig(i,j,ns,na)
                end do
              end do
            end if
          end do
        end do ! l

      end do ! lmn
    end do ! ns
  end do ! nk
! ----------------------------------------------------------------------------

  call overlap_finitedifference_final

  do na=1,natom
    if (natpri(na) .eq. key_natpri_in) then
      iaps=naps(na)
      fatall= fatx(na)
      call mpi_reduce(fatall,fatx(na),1,mpi_double_precision,mpi_sum,0,mpicom_kpt,mpij)
      fatall= faty(na)
      call mpi_reduce(fatall,faty(na),1,mpi_double_precision,mpi_sum,0,mpicom_kpt,mpij)
      fatall= fatz(na)
      call mpi_reduce(fatall,fatz(na),1,mpi_double_precision,mpi_sum,0,mpicom_kpt,mpij)
    end if
  end do

  deallocate(vcm2)
  deallocate(avcm)
  deallocate(vcccm)
  deallocate(vcm)
  deallocate(cpsep,cdpsep)
  deallocate(cpsepeig,cdpsepeig,cpsepall)
  deallocate(dijeig)

end subroutine force_eig_c


subroutine force_r_02(ncpx,ncpy,ncpz,nf,dx,vre,avre)
implicit none
integer, intent(in)::ncpx,ncpy,ncpz
integer, intent(in)::nf
real*8, intent(in)::dx
real*8, intent(in)::vre(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf)
real*8, intent(out)::avre(ncpx,ncpy,ncpz)

integer ix,iy,iz
real*8 a1,a2,a3,a4,a5,a6,a7,a8

  if (nf .eq. 1) then
    a1=-0.5d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a1*vre(ix+1,iy,iz) &
                      +a1*vre(ix-1,iy,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 2) then
    a2=1.0d0/12.0d0/dx
    a1=-8.0d0/12.0d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a2*vre(ix+2,iy,iz) &
                      -a1*vre(ix+1,iy,iz) &
                      +a1*vre(ix-1,iy,iz) &
                      +a2*vre(ix-2,iy,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 3) then
    a3=-1.0d0/60.0d0/dx
    a2=9.0d0/60.0d0/dx
    a1=-45.0d0/60.0d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a3*vre(ix+3,iy,iz) &
                      -a2*vre(ix+2,iy,iz) &
                      -a1*vre(ix+1,iy,iz) &
                      +a1*vre(ix-1,iy,iz) &
                      +a2*vre(ix-2,iy,iz) &
                      +a3*vre(ix-3,iy,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 4) then
    a4=3.0d0/840.0d0/dx
    a3=-32.0d0/840.0d0/dx
    a2=168.0d0/840.0d0/dx
    a1=-672.0d0/840.0d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a4*vre(ix+4,iy,iz) &
                      -a3*vre(ix+3,iy,iz) &
                      -a2*vre(ix+2,iy,iz) &
                      -a1*vre(ix+1,iy,iz) &
                      +a1*vre(ix-1,iy,iz) &
                      +a2*vre(ix-2,iy,iz) &
                      +a3*vre(ix-3,iy,iz) &
                      +a4*vre(ix-4,iy,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 5) then
    a5=-7.936507936507937d-4/dx
    a4= 9.920634920634921d-3/dx
    a3=-5.952380952380952d-2/dx
    a2= 0.2380952380952381d0/dx
    a1=-0.8333333333333333d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a5*vre(ix+5,iy,iz) &
                      -a4*vre(ix+4,iy,iz) &
                      -a3*vre(ix+3,iy,iz) &
                      -a2*vre(ix+2,iy,iz) &
                      -a1*vre(ix+1,iy,iz) &
                      +a1*vre(ix-1,iy,iz) &
                      +a2*vre(ix-2,iy,iz) &
                      +a3*vre(ix-3,iy,iz) &
                      +a4*vre(ix-4,iy,iz) &
                      +a5*vre(ix-5,iy,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 6) then
    a6= 1.803751803751804d-4/dx
    a5=-2.597402597402597d-3/dx
    a4= 1.785714285714286d-2/dx
    a3=-7.936507936507937d-2/dx
    a2= 0.2678571428571429d0/dx
    a1=-0.8571428571428571d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a6*vre(ix+6,iy,iz) &
                      -a5*vre(ix+5,iy,iz) &
                      -a4*vre(ix+4,iy,iz) &
                      -a3*vre(ix+3,iy,iz) &
                      -a2*vre(ix+2,iy,iz) &
                      -a1*vre(ix+1,iy,iz) &
                      +a1*vre(ix-1,iy,iz) &
                      +a2*vre(ix-2,iy,iz) &
                      +a3*vre(ix-3,iy,iz) &
                      +a4*vre(ix-4,iy,iz) &
                      +a5*vre(ix-5,iy,iz) &
                      +a6*vre(ix-6,iy,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 7) then
    a7=-4.162504162504163d-5/dx
    a6= 6.798756798756799d-4/dx
    a5=-5.303030303030303d-3/dx
    a4= 2.651515151515152d-2/dx
    a3=-9.722222222222222d-2/dx
    a2= 0.2916666666666667d0/dx
    a1=-0.875d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a7*vre(ix+7,iy,iz) &
                      -a6*vre(ix+6,iy,iz) &
                      -a5*vre(ix+5,iy,iz) &
                      -a4*vre(ix+4,iy,iz) &
                      -a3*vre(ix+3,iy,iz) &
                      -a2*vre(ix+2,iy,iz) &
                      -a1*vre(ix+1,iy,iz) &
                      +a1*vre(ix-1,iy,iz) &
                      +a2*vre(ix-2,iy,iz) &
                      +a3*vre(ix-3,iy,iz) &
                      +a4*vre(ix-4,iy,iz) &
                      +a5*vre(ix-5,iy,iz) &
                      +a6*vre(ix-6,iy,iz) &
                      +a7*vre(ix-7,iy,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 8) then
    a8= 9.712509712509713d-6/dx
    a7=-1.776001776001776d-4/dx
    a6= 1.554001554001554d-3/dx
    a5=-8.702408702408702d-3/dx
    a4= 3.535353535353535d-2/dx
    a3=-0.1131313131313131d0/dx
    a2= 0.3111111111111111d0/dx
    a1=-0.8888888888888889d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a8*vre(ix+8,iy,iz) &
                      -a7*vre(ix+7,iy,iz) &
                      -a6*vre(ix+6,iy,iz) &
                      -a5*vre(ix+5,iy,iz) &
                      -a4*vre(ix+4,iy,iz) &
                      -a3*vre(ix+3,iy,iz) &
                      -a2*vre(ix+2,iy,iz) &
                      -a1*vre(ix+1,iy,iz) &
                      +a1*vre(ix-1,iy,iz) &
                      +a2*vre(ix-2,iy,iz) &
                      +a3*vre(ix-3,iy,iz) &
                      +a4*vre(ix-4,iy,iz) &
                      +a5*vre(ix-5,iy,iz) &
                      +a6*vre(ix-6,iy,iz) &
                      +a7*vre(ix-7,iy,iz) &
                      +a8*vre(ix-8,iy,iz)
    end do
    end do
    end do
  end if
  return
end subroutine force_r_02


subroutine force_r_03(ncpx,ncpy,ncpz,nf,dy,vre,avre)
implicit none
integer, intent(in)::ncpx,ncpy,ncpz
integer, intent(in)::nf
real*8, intent(in)::dy
real*8, intent(in)::vre(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf)
real*8, intent(out)::avre(ncpx,ncpy,ncpz)
integer ix,iy,iz
real*8 a1,a2,a3,a4,a5,a6,a7,a8

  if (nf .eq. 1) then
    a1=-0.5d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a1*vre(ix,iy+1,iz) &
                      +a1*vre(ix,iy-1,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 2) then
    a2=1.0d0/12.0d0/dy
    a1=-8.0d0/12.0d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a2*vre(ix,iy+2,iz) &
                      -a1*vre(ix,iy+1,iz) &
                      +a1*vre(ix,iy-1,iz) &
                      +a2*vre(ix,iy-2,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 3) then
    a3=-1.0d0/60.0d0/dy
    a2=9.0d0/60.0d0/dy
    a1=-45.0d0/60.0d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a3*vre(ix,iy+3,iz) &
                      -a2*vre(ix,iy+2,iz) &
                      -a1*vre(ix,iy+1,iz) &
                      +a1*vre(ix,iy-1,iz) &
                      +a2*vre(ix,iy-2,iz) &
                      +a3*vre(ix,iy-3,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 4) then
    a4=3.0d0/840.0d0/dy
    a3=-32.0d0/840.0d0/dy
    a2=168.0d0/840.0d0/dy
    a1=-672.0d0/840.0d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a4*vre(ix,iy+4,iz) &
                      -a3*vre(ix,iy+3,iz) &
                      -a2*vre(ix,iy+2,iz) &
                      -a1*vre(ix,iy+1,iz) &
                      +a1*vre(ix,iy-1,iz) &
                      +a2*vre(ix,iy-2,iz) &
                      +a3*vre(ix,iy-3,iz) &
                      +a4*vre(ix,iy-4,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 5) then
    a5=-7.936507936507937d-4/dy
    a4= 9.920634920634921d-3/dy
    a3=-5.952380952380952d-2/dy
    a2= 0.2380952380952381d0/dy
    a1=-0.8333333333333333d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a5*vre(ix,iy+5,iz) &
                      -a4*vre(ix,iy+4,iz) &
                      -a3*vre(ix,iy+3,iz) &
                      -a2*vre(ix,iy+2,iz) &
                      -a1*vre(ix,iy+1,iz) &
                      +a1*vre(ix,iy-1,iz) &
                      +a2*vre(ix,iy-2,iz) &
                      +a3*vre(ix,iy-3,iz) &
                      +a4*vre(ix,iy-4,iz) &
                      +a5*vre(ix,iy-5,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 6) then
    a6= 1.803751803751804d-4/dy
    a5=-2.597402597402597d-3/dy
    a4= 1.785714285714286d-2/dy
    a3=-7.936507936507937d-2/dy
    a2= 0.2678571428571429d0/dy
    a1=-0.8571428571428571d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a6*vre(ix,iy+6,iz) &
                      -a5*vre(ix,iy+5,iz) &
                      -a4*vre(ix,iy+4,iz) &
                      -a3*vre(ix,iy+3,iz) &
                      -a2*vre(ix,iy+2,iz) &
                      -a1*vre(ix,iy+1,iz) &
                      +a1*vre(ix,iy-1,iz) &
                      +a2*vre(ix,iy-2,iz) &
                      +a3*vre(ix,iy-3,iz) &
                      +a4*vre(ix,iy-4,iz) &
                      +a5*vre(ix,iy-5,iz) &
                      +a6*vre(ix,iy-6,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 7) then
    a7=-4.162504162504163d-5/dy
    a6= 6.798756798756799d-4/dy
    a5=-5.303030303030303d-3/dy
    a4= 2.651515151515152d-2/dy
    a3=-9.722222222222222d-2/dy
    a2= 0.2916666666666667d0/dy
    a1=-0.875d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a7*vre(ix,iy+7,iz) &
                      -a6*vre(ix,iy+6,iz) &
                      -a5*vre(ix,iy+5,iz) &
                      -a4*vre(ix,iy+4,iz) &
                      -a3*vre(ix,iy+3,iz) &
                      -a2*vre(ix,iy+2,iz) &
                      -a1*vre(ix,iy+1,iz) &
                      +a1*vre(ix,iy-1,iz) &
                      +a2*vre(ix,iy-2,iz) &
                      +a3*vre(ix,iy-3,iz) &
                      +a4*vre(ix,iy-4,iz) &
                      +a5*vre(ix,iy-5,iz) &
                      +a6*vre(ix,iy-6,iz) &
                      +a7*vre(ix,iy-7,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 8) then
    a8= 9.712509712509713d-6/dy
    a7=-1.776001776001776d-4/dy
    a6= 1.554001554001554d-3/dy
    a5=-8.702408702408702d-3/dy
    a4= 3.535353535353535d-2/dy
    a3=-0.1131313131313131d0/dy
    a2= 0.3111111111111111d0/dy
    a1=-0.8888888888888889d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a8*vre(ix,iy+8,iz) &
                      -a7*vre(ix,iy+7,iz) &
                      -a6*vre(ix,iy+6,iz) &
                      -a5*vre(ix,iy+5,iz) &
                      -a4*vre(ix,iy+4,iz) &
                      -a3*vre(ix,iy+3,iz) &
                      -a2*vre(ix,iy+2,iz) &
                      -a1*vre(ix,iy+1,iz) &
                      +a1*vre(ix,iy-1,iz) &
                      +a2*vre(ix,iy-2,iz) &
                      +a3*vre(ix,iy-3,iz) &
                      +a4*vre(ix,iy-4,iz) &
                      +a5*vre(ix,iy-5,iz) &
                      +a6*vre(ix,iy-6,iz) &
                      +a7*vre(ix,iy-7,iz) &
                      +a8*vre(ix,iy-8,iz)
    end do
    end do
    end do
  end if
  return
end subroutine force_r_03


subroutine force_r_04(ncpx,ncpy,ncpz,nf,dz,vre,avre)
implicit none
integer, intent(in)::ncpx,ncpy,ncpz
integer, intent(in)::nf
real*8, intent(in)::dz
real*8, intent(in)::vre(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf)
real*8, intent(out)::avre(ncpx,ncpy,ncpz)
integer ix,iy,iz
real*8 a1,a2,a3,a4,a5,a6,a7,a8

  if (nf .eq. 1) then
    a1=-0.5d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a1*vre(ix,iy,iz+1) &
                      +a1*vre(ix,iy,iz-1)
    end do
    end do
    end do
  end if

  if (nf .eq. 2) then
    a2=1.0d0/12.0d0/dz
    a1=-8.0d0/12.0d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a2*vre(ix,iy,iz+2) &
                      -a1*vre(ix,iy,iz+1) &
                      +a1*vre(ix,iy,iz-1) &
                      +a2*vre(ix,iy,iz-2)
    end do
    end do
    end do
  end if

  if (nf .eq. 3) then
    a3=-1.0d0/60.0d0/dz
    a2=9.0d0/60.0d0/dz
    a1=-45.0d0/60.0d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a3*vre(ix,iy,iz+3) &
                      -a2*vre(ix,iy,iz+2) &
                      -a1*vre(ix,iy,iz+1) &
                      +a1*vre(ix,iy,iz-1) &
                      +a2*vre(ix,iy,iz-2) &
                      +a3*vre(ix,iy,iz-3)
    end do
    end do
    end do
  end if

  if (nf .eq. 4) then
    a4=3.0d0/840.0d0/dz
    a3=-32.0d0/840.0d0/dz
    a2=168.0d0/840.0d0/dz
    a1=-672.0d0/840.0d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a4*vre(ix,iy,iz+4) &
                      -a3*vre(ix,iy,iz+3) &
                      -a2*vre(ix,iy,iz+2) &
                      -a1*vre(ix,iy,iz+1) &
                      +a1*vre(ix,iy,iz-1) &
                      +a2*vre(ix,iy,iz-2) &
                      +a3*vre(ix,iy,iz-3) &
                      +a4*vre(ix,iy,iz-4)
    end do
    end do
    end do
  end if

  if (nf .eq. 5) then
    a5=-7.936507936507937d-4/dz
    a4= 9.920634920634921d-3/dz
    a3=-5.952380952380952d-2/dz
    a2= 0.2380952380952381d0/dz
    a1=-0.8333333333333333d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a5*vre(ix,iy,iz+5) &
                      -a4*vre(ix,iy,iz+4) &
                      -a3*vre(ix,iy,iz+3) &
                      -a2*vre(ix,iy,iz+2) &
                      -a1*vre(ix,iy,iz+1) &
                      +a1*vre(ix,iy,iz-1) &
                      +a2*vre(ix,iy,iz-2) &
                      +a3*vre(ix,iy,iz-3) &
                      +a4*vre(ix,iy,iz-4) &
                      +a5*vre(ix,iy,iz-5)
    end do
    end do
    end do
  end if

  if (nf .eq. 6) then
    a6= 1.803751803751804d-4/dz
    a5=-2.597402597402597d-3/dz
    a4= 1.785714285714286d-2/dz
    a3=-7.936507936507937d-2/dz
    a2= 0.2678571428571429d0/dz
    a1=-0.8571428571428571d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a6*vre(ix,iy,iz+6) &
                      -a5*vre(ix,iy,iz+5) &
                      -a4*vre(ix,iy,iz+4) &
                      -a3*vre(ix,iy,iz+3) &
                      -a2*vre(ix,iy,iz+2) &
                      -a1*vre(ix,iy,iz+1) &
                      +a1*vre(ix,iy,iz-1) &
                      +a2*vre(ix,iy,iz-2) &
                      +a3*vre(ix,iy,iz-3) &
                      +a4*vre(ix,iy,iz-4) &
                      +a5*vre(ix,iy,iz-5) &
                      +a6*vre(ix,iy,iz-6)
    end do
    end do
    end do
  end if

  if (nf .eq. 7) then
    a7=-4.162504162504163d-5/dz
    a6= 6.798756798756799d-4/dz
    a5=-5.303030303030303d-3/dz
    a4= 2.651515151515152d-2/dz
    a3=-9.722222222222222d-2/dz
    a2= 0.2916666666666667d0/dz
    a1=-0.875d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a7*vre(ix,iy,iz+7) &
                      -a6*vre(ix,iy,iz+6) &
                      -a5*vre(ix,iy,iz+5) &
                      -a4*vre(ix,iy,iz+4) &
                      -a3*vre(ix,iy,iz+3) &
                      -a2*vre(ix,iy,iz+2) &
                      -a1*vre(ix,iy,iz+1) &
                      +a1*vre(ix,iy,iz-1) &
                      +a2*vre(ix,iy,iz-2) &
                      +a3*vre(ix,iy,iz-3) &
                      +a4*vre(ix,iy,iz-4) &
                      +a5*vre(ix,iy,iz-5) &
                      +a6*vre(ix,iy,iz-6) &
                      +a7*vre(ix,iy,iz-7)
    end do
    end do
    end do
  end if

  if (nf .eq. 8) then
    a8= 9.712509712509713d-6/dz
    a7=-1.776001776001776d-4/dz
    a6= 1.554001554001554d-3/dz
    a5=-8.702408702408702d-3/dz
    a4= 3.535353535353535d-2/dz
    a3=-0.1131313131313131d0/dz
    a2= 0.3111111111111111d0/dz
    a1=-0.8888888888888889d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avre(ix,iy,iz)=-a8*vre(ix,iy,iz+8) &
                      -a7*vre(ix,iy,iz+7) &
                      -a6*vre(ix,iy,iz+6) &
                      -a5*vre(ix,iy,iz+5) &
                      -a4*vre(ix,iy,iz+4) &
                      -a3*vre(ix,iy,iz+3) &
                      -a2*vre(ix,iy,iz+2) &
                      -a1*vre(ix,iy,iz+1) &
                      +a1*vre(ix,iy,iz-1) &
                      +a2*vre(ix,iy,iz-2) &
                      +a3*vre(ix,iy,iz-3) &
                      +a4*vre(ix,iy,iz-4) &
                      +a5*vre(ix,iy,iz-5) &
                      +a6*vre(ix,iy,iz-6) &
                      +a7*vre(ix,iy,iz-7) &
                      +a8*vre(ix,iy,iz-8)
    end do
    end do
    end do
  end if
  return
end subroutine force_r_04


subroutine force_c_02(ncpx,ncpy,ncpz,nf,dx,vcm,avcm)
implicit none
integer,   intent(in)::ncpx,ncpy,ncpz
integer,   intent(in)::nf
real*8,    intent(in)::dx
complex*16,intent(in)::vcm(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf)
complex*16,intent(out)::avcm(ncpx,ncpy,ncpz)
integer ix,iy,iz
real*8 a1,a2,a3,a4,a5,a6,a7,a8

  if (nf .eq. 1) then
    a1=-0.5d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a1*vcm(ix+1,iy,iz) &
                      +a1*vcm(ix-1,iy,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 2) then
    a2=1.0d0/12.0d0/dx
    a1=-8.0d0/12.0d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a2*vcm(ix+2,iy,iz) &
                      -a1*vcm(ix+1,iy,iz) &
                      +a1*vcm(ix-1,iy,iz) &
                      +a2*vcm(ix-2,iy,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 3) then
    a3=-1.0d0/60.0d0/dx
    a2=9.0d0/60.0d0/dx
    a1=-45.0d0/60.0d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a3*vcm(ix+3,iy,iz) &
                      -a2*vcm(ix+2,iy,iz) &
                      -a1*vcm(ix+1,iy,iz) &
                      +a1*vcm(ix-1,iy,iz) &
                      +a2*vcm(ix-2,iy,iz) &
                      +a3*vcm(ix-3,iy,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 4) then
    a4=3.0d0/840.0d0/dx
    a3=-32.0d0/840.0d0/dx
    a2=168.0d0/840.0d0/dx
    a1=-672.0d0/840.0d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a4*vcm(ix+4,iy,iz) &
                      -a3*vcm(ix+3,iy,iz) &
                      -a2*vcm(ix+2,iy,iz) &
                      -a1*vcm(ix+1,iy,iz) &
                      +a1*vcm(ix-1,iy,iz) &
                      +a2*vcm(ix-2,iy,iz) &
                      +a3*vcm(ix-3,iy,iz) &
                      +a4*vcm(ix-4,iy,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 5) then
    a5=-7.936507936507937d-4/dx
    a4= 9.920634920634921d-3/dx
    a3=-5.952380952380952d-2/dx
    a2= 0.2380952380952381d0/dx
    a1=-0.8333333333333333d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a5*vcm(ix+5,iy,iz) &
                      -a4*vcm(ix+4,iy,iz) &
                      -a3*vcm(ix+3,iy,iz) &
                      -a2*vcm(ix+2,iy,iz) &
                      -a1*vcm(ix+1,iy,iz) &
                      +a1*vcm(ix-1,iy,iz) &
                      +a2*vcm(ix-2,iy,iz) &
                      +a3*vcm(ix-3,iy,iz) &
                      +a4*vcm(ix-4,iy,iz) &
                      +a5*vcm(ix-5,iy,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 6) then
    a6= 1.803751803751804d-4/dx
    a5=-2.597402597402597d-3/dx
    a4= 1.785714285714286d-2/dx
    a3=-7.936507936507937d-2/dx
    a2= 0.2678571428571429d0/dx
    a1=-0.8571428571428571d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a6*vcm(ix+6,iy,iz) &
                      -a5*vcm(ix+5,iy,iz) &
                      -a4*vcm(ix+4,iy,iz) &
                      -a3*vcm(ix+3,iy,iz) &
                      -a2*vcm(ix+2,iy,iz) &
                      -a1*vcm(ix+1,iy,iz) &
                      +a1*vcm(ix-1,iy,iz) &
                      +a2*vcm(ix-2,iy,iz) &
                      +a3*vcm(ix-3,iy,iz) &
                      +a4*vcm(ix-4,iy,iz) &
                      +a5*vcm(ix-5,iy,iz) &
                      +a6*vcm(ix-6,iy,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 7) then
    a7=-4.162504162504163d-5/dx
    a6= 6.798756798756799d-4/dx
    a5=-5.303030303030303d-3/dx
    a4= 2.651515151515152d-2/dx
    a3=-9.722222222222222d-2/dx
    a2= 0.2916666666666667d0/dx
    a1=-0.875d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a7*vcm(ix+7,iy,iz) &
                      -a6*vcm(ix+6,iy,iz) &
                      -a5*vcm(ix+5,iy,iz) &
                      -a4*vcm(ix+4,iy,iz) &
                      -a3*vcm(ix+3,iy,iz) &
                      -a2*vcm(ix+2,iy,iz) &
                      -a1*vcm(ix+1,iy,iz) &
                      +a1*vcm(ix-1,iy,iz) &
                      +a2*vcm(ix-2,iy,iz) &
                      +a3*vcm(ix-3,iy,iz) &
                      +a4*vcm(ix-4,iy,iz) &
                      +a5*vcm(ix-5,iy,iz) &
                      +a6*vcm(ix-6,iy,iz) &
                      +a7*vcm(ix-7,iy,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 8) then
    a8= 9.712509712509713d-6/dx
    a7=-1.776001776001776d-4/dx
    a6= 1.554001554001554d-3/dx
    a5=-8.702408702408702d-3/dx
    a4= 3.535353535353535d-2/dx
    a3=-0.1131313131313131d0/dx
    a2= 0.3111111111111111d0/dx
    a1=-0.8888888888888889d0/dx
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a8*vcm(ix+8,iy,iz) &
                      -a7*vcm(ix+7,iy,iz) &
                      -a6*vcm(ix+6,iy,iz) &
                      -a5*vcm(ix+5,iy,iz) &
                      -a4*vcm(ix+4,iy,iz) &
                      -a3*vcm(ix+3,iy,iz) &
                      -a2*vcm(ix+2,iy,iz) &
                      -a1*vcm(ix+1,iy,iz) &
                      +a1*vcm(ix-1,iy,iz) &
                      +a2*vcm(ix-2,iy,iz) &
                      +a3*vcm(ix-3,iy,iz) &
                      +a4*vcm(ix-4,iy,iz) &
                      +a5*vcm(ix-5,iy,iz) &
                      +a6*vcm(ix-6,iy,iz) &
                      +a7*vcm(ix-7,iy,iz) &
                      +a8*vcm(ix-8,iy,iz)
    end do
    end do
    end do
  end if
  return
end subroutine force_c_02


subroutine force_c_03(ncpx,ncpy,ncpz,nf,dy,vcm,avcm)
implicit none
integer,   intent(in)::ncpx,ncpy,ncpz
integer,   intent(in)::nf
real*8,    intent(in)::dy
complex*16,intent(in)::vcm(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf)
complex*16,intent(out)::avcm(ncpx,ncpy,ncpz)
integer ix,iy,iz
real*8 a1,a2,a3,a4,a5,a6,a7,a8

  if (nf .eq. 1) then
    a1=-0.5d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a1*vcm(ix,iy+1,iz) &
                      +a1*vcm(ix,iy-1,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 2) then
    a2=1.0d0/12.0d0/dy
    a1=-8.0d0/12.0d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a2*vcm(ix,iy+2,iz) &
                      -a1*vcm(ix,iy+1,iz) &
                      +a1*vcm(ix,iy-1,iz) &
                      +a2*vcm(ix,iy-2,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 3) then
    a3=-1.0d0/60.0d0/dy
    a2=9.0d0/60.0d0/dy
    a1=-45.0d0/60.0d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a3*vcm(ix,iy+3,iz) &
                      -a2*vcm(ix,iy+2,iz) &
                      -a1*vcm(ix,iy+1,iz) &
                      +a1*vcm(ix,iy-1,iz) &
                      +a2*vcm(ix,iy-2,iz) &
                      +a3*vcm(ix,iy-3,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 4) then
    a4=3.0d0/840.0d0/dy
    a3=-32.0d0/840.0d0/dy
    a2=168.0d0/840.0d0/dy
    a1=-672.0d0/840.0d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a4*vcm(ix,iy+4,iz) &
                      -a3*vcm(ix,iy+3,iz) &
                      -a2*vcm(ix,iy+2,iz) &
                      -a1*vcm(ix,iy+1,iz) &
                      +a1*vcm(ix,iy-1,iz) &
                      +a2*vcm(ix,iy-2,iz) &
                      +a3*vcm(ix,iy-3,iz) &
                      +a4*vcm(ix,iy-4,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 5) then
    a5=-7.936507936507937d-4/dy
    a4= 9.920634920634921d-3/dy
    a3=-5.952380952380952d-2/dy
    a2= 0.2380952380952381d0/dy
    a1=-0.8333333333333333d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a5*vcm(ix,iy+5,iz) &
                      -a4*vcm(ix,iy+4,iz) &
                      -a3*vcm(ix,iy+3,iz) &
                      -a2*vcm(ix,iy+2,iz) &
                      -a1*vcm(ix,iy+1,iz) &
                      +a1*vcm(ix,iy-1,iz) &
                      +a2*vcm(ix,iy-2,iz) &
                      +a3*vcm(ix,iy-3,iz) &
                      +a4*vcm(ix,iy-4,iz) &
                      +a5*vcm(ix,iy-5,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 6) then
    a6= 1.803751803751804d-4/dy
    a5=-2.597402597402597d-3/dy
    a4= 1.785714285714286d-2/dy
    a3=-7.936507936507937d-2/dy
    a2= 0.2678571428571429d0/dy
    a1=-0.8571428571428571d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a6*vcm(ix,iy+6,iz) &
                      -a5*vcm(ix,iy+5,iz) &
                      -a4*vcm(ix,iy+4,iz) &
                      -a3*vcm(ix,iy+3,iz) &
                      -a2*vcm(ix,iy+2,iz) &
                      -a1*vcm(ix,iy+1,iz) &
                      +a1*vcm(ix,iy-1,iz) &
                      +a2*vcm(ix,iy-2,iz) &
                      +a3*vcm(ix,iy-3,iz) &
                      +a4*vcm(ix,iy-4,iz) &
                      +a5*vcm(ix,iy-5,iz) &
                      +a6*vcm(ix,iy-6,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 7) then
    a7=-4.162504162504163d-5/dy
    a6= 6.798756798756799d-4/dy
    a5=-5.303030303030303d-3/dy
    a4= 2.651515151515152d-2/dy
    a3=-9.722222222222222d-2/dy
    a2= 0.2916666666666667d0/dy
    a1=-0.875d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a7*vcm(ix,iy+7,iz) &
                      -a6*vcm(ix,iy+6,iz) &
                      -a5*vcm(ix,iy+5,iz) &
                      -a4*vcm(ix,iy+4,iz) &
                      -a3*vcm(ix,iy+3,iz) &
                      -a2*vcm(ix,iy+2,iz) &
                      -a1*vcm(ix,iy+1,iz) &
                      +a1*vcm(ix,iy-1,iz) &
                      +a2*vcm(ix,iy-2,iz) &
                      +a3*vcm(ix,iy-3,iz) &
                      +a4*vcm(ix,iy-4,iz) &
                      +a5*vcm(ix,iy-5,iz) &
                      +a6*vcm(ix,iy-6,iz) &
                      +a7*vcm(ix,iy-7,iz)
    end do
    end do
    end do
  end if

  if (nf .eq. 8) then
    a8= 9.712509712509713d-6/dy
    a7=-1.776001776001776d-4/dy
    a6= 1.554001554001554d-3/dy
    a5=-8.702408702408702d-3/dy
    a4= 3.535353535353535d-2/dy
    a3=-0.1131313131313131d0/dy
    a2= 0.3111111111111111d0/dy
    a1=-0.8888888888888889d0/dy
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a8*vcm(ix,iy+8,iz) &
                      -a7*vcm(ix,iy+7,iz) &
                      -a6*vcm(ix,iy+6,iz) &
                      -a5*vcm(ix,iy+5,iz) &
                      -a4*vcm(ix,iy+4,iz) &
                      -a3*vcm(ix,iy+3,iz) &
                      -a2*vcm(ix,iy+2,iz) &
                      -a1*vcm(ix,iy+1,iz) &
                      +a1*vcm(ix,iy-1,iz) &
                      +a2*vcm(ix,iy-2,iz) &
                      +a3*vcm(ix,iy-3,iz) &
                      +a4*vcm(ix,iy-4,iz) &
                      +a5*vcm(ix,iy-5,iz) &
                      +a6*vcm(ix,iy-6,iz) &
                      +a7*vcm(ix,iy-7,iz) &
                      +a8*vcm(ix,iy-8,iz)
    end do
    end do
    end do
  end if
  return
end subroutine force_c_03


subroutine force_c_04(ncpx,ncpy,ncpz,nf,dz,vcm,avcm)
implicit none
integer,   intent(in)::ncpx,ncpy,ncpz
integer,   intent(in)::nf
real*8,    intent(in)::dz
complex*16,intent(in)::vcm(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf)
complex*16,intent(out)::avcm(ncpx,ncpy,ncpz)
integer    ix,iy,iz
real*8 a1,a2,a3,a4,a5,a6,a7,a8

  if (nf .eq. 1) then
    a1=-0.5d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a1*vcm(ix,iy,iz+1) &
                      +a1*vcm(ix,iy,iz-1)
    end do
    end do
    end do
  end if

  if (nf .eq. 2) then
    a2=1.0d0/12.0d0/dz
    a1=-8.0d0/12.0d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a2*vcm(ix,iy,iz+2) &
                      -a1*vcm(ix,iy,iz+1) &
                      +a1*vcm(ix,iy,iz-1) &
                      +a2*vcm(ix,iy,iz-2)
    end do
    end do
    end do
  end if

  if (nf .eq. 3) then
    a3=-1.0d0/60.0d0/dz
    a2=9.0d0/60.0d0/dz
    a1=-45.0d0/60.0d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a3*vcm(ix,iy,iz+3) &
                      -a2*vcm(ix,iy,iz+2) &
                      -a1*vcm(ix,iy,iz+1) &
                      +a1*vcm(ix,iy,iz-1) &
                      +a2*vcm(ix,iy,iz-2) &
                      +a3*vcm(ix,iy,iz-3)
    end do
    end do
    end do
  end if

  if (nf .eq. 4) then
    a4=3.0d0/840.0d0/dz
    a3=-32.0d0/840.0d0/dz
    a2=168.0d0/840.0d0/dz
    a1=-672.0d0/840.0d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a4*vcm(ix,iy,iz+4) &
                      -a3*vcm(ix,iy,iz+3) &
                      -a2*vcm(ix,iy,iz+2) &
                      -a1*vcm(ix,iy,iz+1) &
                      +a1*vcm(ix,iy,iz-1) &
                      +a2*vcm(ix,iy,iz-2) &
                      +a3*vcm(ix,iy,iz-3) &
                      +a4*vcm(ix,iy,iz-4)
    end do
    end do
    end do
  end if

  if (nf .eq. 5) then
    a5=-7.936507936507937d-4/dz
    a4= 9.920634920634921d-3/dz
    a3=-5.952380952380952d-2/dz
    a2= 0.2380952380952381d0/dz
    a1=-0.8333333333333333d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a5*vcm(ix,iy,iz+5) &
                      -a4*vcm(ix,iy,iz+4) &
                      -a3*vcm(ix,iy,iz+3) &
                      -a2*vcm(ix,iy,iz+2) &
                      -a1*vcm(ix,iy,iz+1) &
                      +a1*vcm(ix,iy,iz-1) &
                      +a2*vcm(ix,iy,iz-2) &
                      +a3*vcm(ix,iy,iz-3) &
                      +a4*vcm(ix,iy,iz-4) &
                      +a5*vcm(ix,iy,iz-5)
    end do
    end do
    end do
  end if

  if (nf .eq. 6) then
    a6= 1.803751803751804d-4/dz
    a5=-2.597402597402597d-3/dz
    a4= 1.785714285714286d-2/dz
    a3=-7.936507936507937d-2/dz
    a2= 0.2678571428571429d0/dz
    a1=-0.8571428571428571d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a6*vcm(ix,iy,iz+6) &
                      -a5*vcm(ix,iy,iz+5) &
                      -a4*vcm(ix,iy,iz+4) &
                      -a3*vcm(ix,iy,iz+3) &
                      -a2*vcm(ix,iy,iz+2) &
                      -a1*vcm(ix,iy,iz+1) &
                      +a1*vcm(ix,iy,iz-1) &
                      +a2*vcm(ix,iy,iz-2) &
                      +a3*vcm(ix,iy,iz-3) &
                      +a4*vcm(ix,iy,iz-4) &
                      +a5*vcm(ix,iy,iz-5) &
                      +a6*vcm(ix,iy,iz-6)
    end do
    end do
    end do
  end if

  if (nf .eq. 7) then
    a7=-4.162504162504163d-5/dz
    a6= 6.798756798756799d-4/dz
    a5=-5.303030303030303d-3/dz
    a4= 2.651515151515152d-2/dz
    a3=-9.722222222222222d-2/dz
    a2= 0.2916666666666667d0/dz
    a1=-0.875d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a7*vcm(ix,iy,iz+7) &
                      -a6*vcm(ix,iy,iz+6) &
                      -a5*vcm(ix,iy,iz+5) &
                      -a4*vcm(ix,iy,iz+4) &
                      -a3*vcm(ix,iy,iz+3) &
                      -a2*vcm(ix,iy,iz+2) &
                      -a1*vcm(ix,iy,iz+1) &
                      +a1*vcm(ix,iy,iz-1) &
                      +a2*vcm(ix,iy,iz-2) &
                      +a3*vcm(ix,iy,iz-3) &
                      +a4*vcm(ix,iy,iz-4) &
                      +a5*vcm(ix,iy,iz-5) &
                      +a6*vcm(ix,iy,iz-6) &
                      +a7*vcm(ix,iy,iz-7)
    end do
    end do
    end do
  end if

  if (nf .eq. 8) then
    a8= 9.712509712509713d-6/dz
    a7=-1.776001776001776d-4/dz
    a6= 1.554001554001554d-3/dz
    a5=-8.702408702408702d-3/dz
    a4= 3.535353535353535d-2/dz
    a3=-0.1131313131313131d0/dz
    a2= 0.3111111111111111d0/dz
    a1=-0.8888888888888889d0/dz
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
       avcm(ix,iy,iz)=-a8*vcm(ix,iy,iz+8) &
                      -a7*vcm(ix,iy,iz+7) &
                      -a6*vcm(ix,iy,iz+6) &
                      -a5*vcm(ix,iy,iz+5) &
                      -a4*vcm(ix,iy,iz+4) &
                      -a3*vcm(ix,iy,iz+3) &
                      -a2*vcm(ix,iy,iz+2) &
                      -a1*vcm(ix,iy,iz+1) &
                      +a1*vcm(ix,iy,iz-1) &
                      +a2*vcm(ix,iy,iz-2) &
                      +a3*vcm(ix,iy,iz-3) &
                      +a4*vcm(ix,iy,iz-4) &
                      +a5*vcm(ix,iy,iz-5) &
                      +a6*vcm(ix,iy,iz-6) &
                      +a7*vcm(ix,iy,iz-7) &
                      +a8*vcm(ix,iy,iz-8)
    end do
    end do
    end do
  end if
  return
end subroutine force_c_04

! ========================================================================================

subroutine force(natom,num_spe,num_atcell,num_ppcell_d,num_list_d,nqmx,nradmx,npoint,                 & ! <
                 nmesh,nspv,nperi,nint1dmax,nzmax,jelcalc,                                            & ! <
                 ncpx,ncpy,ncpz,ncpx_d,ncpy_d,ncpz_d,                                                 & ! <
                 new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz,                                     & ! <
                 key_jel_calc,key_natpri_in,key_natpri_inps,                                          & ! <
                 indspe,natpri,natpri_inf,nradct,natprid,napsd,natinfd,natinfd_vloc,lstvecd2,nqctpcc, & ! <
                 lstdx,lstdy,lstdz,                                                                   & ! <
                 veta,psctoff,psftrad,radial,dradial,point,wt,xmax,ymax,zmax,biasx,biasy,biasz,       & ! <
                 cp,coef,                                                                             & ! <
                 atx,aty,atz,                                                                         & ! <
                 dvlocdx_scw,dvlocdy_scw,dvlocdz_scw,                                                 & ! <
                 dvlocdx_hdp,dvlocdy_hdp,dvlocdz_hdp,                                                 & ! <
                 rhotrur,rhosmtr,rhoaugr,drhoaug3ddx,drhoaug3ddy,drhoaug3ddz,                         & ! <
                 rho_coarse,rho_aug_dense,vloc_dense,vh_dense,vx_dense,                               & ! <
                 chrjel,strjel,endjel,                                                                & ! <
                 fatx,faty,fatz,                                                                      & ! X
                 fjel)                                                                                  ! >
use mod_mpi
implicit none
integer,intent(in)   ::natom,num_spe,num_atcell,num_ppcell_d,num_list_d,nqmx,nradmx,npoint
integer,intent(in)   ::nmesh,nspv,nperi,nint1dmax,nzmax,jelcalc
integer,intent(in)   ::ncpx,ncpy,ncpz,ncpx_d,ncpy_d,ncpz_d
integer,intent(in)   ::new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz
integer,intent(in)   ::key_jel_calc,key_natpri_in,key_natpri_inps
integer,intent(in)   ::indspe(natom),natpri(natom),natpri_inf(natom),nradct(num_spe),natprid(natom)
integer,intent(in)   ::napsd(natom),natinfd(natom),natinfd_vloc(natom),nqctpcc(num_spe)
integer,intent(in)   ::lstdx(num_list_d,num_ppcell_d),lstdy(num_list_d,num_ppcell_d),lstdz(num_list_d,num_ppcell_d)
integer,intent(in)   ::lstvecd2(num_list_d,num_ppcell_d)
real*8, intent(in)   ::veta,psctoff,psftrad
real*8, intent(in)   ::xmax,ymax,zmax
real*8, intent(in)   ::biasx,biasy,biasz
real*8, intent(in)   ::radial(nradmx,num_spe),dradial(nradmx,num_spe),point(npoint,3),wt(npoint)
real*8, intent(in)   ::cp(8,num_spe),coef(0:nqmx,0:7,num_spe)
real*8, intent(in)   ::atx(natom),aty(natom),atz(natom)
real*8, intent(in)   ::dvlocdx_scw(ncpx,ncpy,ncpz,natom),dvlocdy_scw(ncpx,ncpy,ncpz,natom),dvlocdz_scw(ncpx,ncpy,ncpz,natom)
real*8, intent(in)   ::dvlocdx_hdp(num_list_d,num_ppcell_d)
real*8, intent(in)   ::dvlocdy_hdp(num_list_d,num_ppcell_d)
real*8, intent(in)   ::dvlocdz_hdp(num_list_d,num_ppcell_d)
real*8, intent(in)   ::rhotrur(nradmx,npoint,nspv,num_atcell),rhosmtr(nradmx,npoint,nspv,num_atcell)
real*8, intent(in)   ::rhoaugr(nradmx,npoint,num_atcell)
real*8, intent(in)   ::drhoaug3ddx(num_list_d,num_ppcell_d)
real*8, intent(in)   ::drhoaug3ddy(num_list_d,num_ppcell_d)
real*8, intent(in)   ::drhoaug3ddz(num_list_d,num_ppcell_d)
real*8, intent(in)   ::rho_coarse(ncpx,ncpy,ncpz)
real*8, intent(in)   ::rho_aug_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
real*8, intent(in)   ::vloc_dense(ncpx_d,ncpy_d,ncpz_d)
real*8, intent(in)   ::vh_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
real*8, intent(in)   ::vx_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh,nspv)
real*8, intent(in)   ::chrjel,strjel,endjel
real*8, intent(inout)::fatx(natom),faty(natom),fatz(natom)
real*8, intent(out)  ::fjel
real*8 dx,dy,dz,ddx,ddy,ddz,fxall,fyall,fzall
integer na,kx,ky,kz,ix,iy,iz,l0,l1,ir,il,ipri,jz,ns,i,nas,nae,natpprcs
real*8 pi,fourpi,omega,omegain,surf,surfin,tmpx,tmpy,tmpz,x,y,z,r,x0,x1,y0,y1,z0,z1,vkx,vky,vkz,vlo0,vlo1,vlo0r,vk2
real*8 cp1,cp2,cp3,cp4,cp5,dr,zz,fjel0,fatxx,fatyy,fatzz,ta,tb,vep0,vep1,vsine,vcosi,t,dt
real*8,allocatable::fatall(:)
real*8 derf

  dx=2.0d0*xmax/dble(ncpx*nprocx)
  dy=2.0d0*ymax/dble(ncpy*nprocy)
  dz=2.0d0*zmax/dble(ncpz*nprocz)
  ddx=dx/dble(nmesh)
  ddy=dy/dble(nmesh)
  ddz=dz/dble(nmesh)
  nas=1
  do i= 0,myrank_glbl-1
    natpprcs=natom/nprocs
    if ( i >= nprocs-mod(natom,nprocs) ) natpprcs=natom/nprocs+1
    nas=nas+natpprcs
  end do
  natpprcs=natom/nprocs
  if ( myrank_glbl >= nprocs-mod(natom,nprocs) ) natpprcs=natom/nprocs+1
  nae=nas+natpprcs-1

  allocate(fatall(natom))

  pi=dacos(-1.0d0)
  fourpi=4.0d0*pi
  omega=8.0d0*xmax*ymax*zmax
  omegain=1.0d0/omega
  surf=4.0d0*xmax*ymax
  surfin=1.0d0/surf

! ----------------------------------------------------------------------------

  do na=1,natom
    tmpx=0.0d0
    tmpy=0.0d0
    tmpz=0.0d0
!$omp parallel default(shared)
    call force_01(natom,nmesh,na,ncpx,ncpy,ncpz,num_ppcell_d,num_list_d, & ! <
                  key_natpri_inps,                                       & ! <
                  dx,dy,dz,                                              & ! <
                  tmpx,tmpy,tmpz,                                        & ! X
                  natprid,napsd,natinfd,natinfd_vloc,lstvecd2,           & ! <
                  lstdx,lstdy,lstdz,                                     & ! <
                  rho_coarse,rho_aug_dense,vh_dense,                     & ! <
                  drhoaug3ddx,drhoaug3ddy,drhoaug3ddz,                   & ! <
                  dvlocdx_scw,dvlocdy_scw,dvlocdz_scw,                   & ! <
                  dvlocdx_hdp,dvlocdy_hdp,dvlocdz_hdp,                   & ! <
                  vloc_dense)                                              ! <
!$omp end parallel
    fatx(na)=fatx(na)+tmpx
    faty(na)=faty(na)+tmpy
    fatz(na)=fatz(na)+tmpz
  end do

! ----------  pcc part  ----------
!$omp parallel default(shared)
  call force_02(natom,num_spe,nspv,nradmx,nqmx,ncpx_d,ncpy_d,ncpz_d, & ! <
                xmax,ymax,zmax,ddx,ddy,ddz,psctoff,psftrad, &          ! <
                nradct,indspe,nqctpcc, &                               ! <
                fatx,faty,fatz, &                                      ! X
                atx,aty,atz, &                                         ! <
                radial,coef,vx_dense)                                  ! <

!  call force_02(natom,nspv,ncpx_d,ncpy_d,ncpz_d,xmax,ymax,zmax,ddx,ddy,ddz,psctoff,psftrad &
!                     ,fatx,faty,fatz,atx,aty,atz,vx_dense)
!$omp end parallel
! --------------------------------

! ----------  jellium part I  ----------
  fjel= 0.0d0
  if (jelcalc==key_jel_calc) then
    do na=1,natom
      cp1=cp(1,indspe(na))
      cp2=cp(2,indspe(na))
      cp3=cp(3,indspe(na))
      cp4=cp(4,indspe(na))
      cp5=cp(5,indspe(na))
      do kz=-nzmax,nzmax
        if (kz**2 .ne. 0) then
          vkz=pi/zmax*kz
          vk2=vkz*vkz
          fjel0=0.0d0
          do iz=1,ncpz_d
          do iy=1,ncpy_d
          do ix=1,ncpx_d
            x=(ix+myrx*ncpx_d)*ddx-0.5d0*ddx-xmax-atx(na)
            y=(iy+myry*ncpy_d)*ddy-0.5d0*ddy-ymax-aty(na)
            z=(iz+myrz*ncpz_d)*ddz-0.5d0*ddz-zmax-atz(na)
            zz=(iz+myrz*ncpz_d)*ddz-0.5d0*ddz-zmax
            if (x .gt.  xmax) x=x-2.0d0*xmax
            if (x .le. -xmax) x=x+2.0d0*xmax
            if (y .gt.  ymax) y=y-2.0d0*ymax
            if (y .le. -ymax) y=y+2.0d0*ymax
            if (z .gt.  zmax) z=z-2.0d0*zmax
            if (z .le. -zmax) z=z+2.0d0*zmax
            r=dsqrt(x*x+y*y+z*z)
            fjel0   =fjel0-4.0d0*pi/omega*chrjel/(strjel-endjel)/vk2*(dcos(vkz*(strjel-zz))-dcos(vkz*(endjel-zz))) &
                   *cp1*(cp4*(cp2/pi)**1.5d0*dexp(-cp2*r*r)+cp5*(cp3/pi)**1.5d0*dexp(-cp3*r*r))*ddx*ddy*ddz
          end do
          end do
          end do
          fjel=fjel+fjel0
          fatz(na)=fatz(na)-fjel0
        end if
      end do
    end do
  end if
! --------------------------------------

! ----------  core - core interaction  ----------
!$omp parallel default(shared) private(l1,x0,y0,z0,x1,y1,z1,x,y,z,r,fatxx,fatyy,fatzz,kx,ky,kz,vkx,vky,vkz,vk2,vlo0)
  do l0=nas,nae
    fatxx=0.0d0
    fatyy=0.0d0
    fatzz=0.0d0
    x0=atx(l0)
    y0=aty(l0)
    z0=atz(l0)
!$omp do
    do l1=1,natom
      x1=atx(l1)
      y1=aty(l1)
      z1=atz(l1)
      select case (nperi)
      case (0)
        x=x0-x1
        y=y0-y1
        z=z0-z1
        r=dsqrt(x*x+y*y+z*z)
        if (r .gt. 1.0d-20) then
          fatxx=fatxx+cp(1,indspe(l0))*cp(1,indspe(l1))*x/r**3
          fatyy=fatyy+cp(1,indspe(l0))*cp(1,indspe(l1))*y/r**3
          fatzz=fatzz+cp(1,indspe(l0))*cp(1,indspe(l1))*z/r**3
        end if
      case (1)
        do kx=-new_pwx,new_pwx
        if (kx**2 .ne. 0) then
          vkx=pi/xmax*kx
          ta=((y0-y1)*(y0-y1)+(z0-z1)*(z0-z1))
          tb=vkx*vkx
          dt=veta/nint1dmax
          vlo0=0.0d0
          vlo1=0.0d0
          do i=1,nint1dmax
            t=i*dt
            vlo0=vlo0+dexp(-ta*t*t-0.25d0*tb/(t*t))/t*dt
            vlo1=vlo1+dexp(-ta*t*t-0.25d0*tb/(t*t))*t*dt
          end do
          vcosi=-dcos(vkx*(x0-x1))/xmax*0.5d0
          vsine=-dsin(vkx*(x0-x1))/xmax*0.5d0
          fatxx=fatxx-2.0d0*vkx*cp(1,indspe(l0))*cp(1,indspe(l1))*vsine*vlo0
          fatyy=fatyy-4.0d0*cp(1,indspe(l0))*cp(1,indspe(l1))*vcosi*vlo1*(y0-y1)
          fatzz=fatzz-4.0d0*cp(1,indspe(l0))*cp(1,indspe(l1))*vcosi*vlo1*(z0-z1)
        end if
        end do
        ta=((y0-y1)*(y0-y1)+(z0-z1)*(z0-z1))
        if (ta .gt. 1.0d-16) then
          vlo0r=-2.0d0*(1.0d0-dexp(-veta*veta*ta))/(ta*2.0d0*xmax)
        else
          vlo0r=-veta*veta/xmax
        end if
        fatyy=fatyy-cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0r*(y0-y1)
        fatzz=fatzz-cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0r*(z0-z1)
        do kx=-new_pwx,new_pwx
          x=kx*2.0d0*xmax+x0-x1
          y=y0-y1
          z=z0-z1
          r=dsqrt(x*x+y*y+z*z)
          if ((r .gt. 1.0d-20) .and.(veta**2*r**2 .lt. 100.0d0))then
            vlo0=-2.0d0*veta/(dexp(veta**2*r**2)*dsqrt(pi)*r)-(1.0d0-derf(veta*r))/r**2
            fatxx=fatxx-cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0*x/r
            fatyy=fatyy-cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0*y/r
            fatzz=fatzz-cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0*z/r
          end if
        end do
      case (2)
        do ky=-new_pwy,new_pwy
        do kx=-new_pwx,new_pwx
          if (kx**2+ky**2 .ne. 0) then
            vkx=pi/xmax*kx
            vky=pi/ymax*ky
            vk2=vkx**2+vky**2
            ta=(z0-z1)
            tb=dsqrt(vk2)
            vep0=-pi/tb/surf*((1.0d0-derf((tb-2.0d0*veta*veta*ta)/(2.0d0*veta)))/dexp(tb*ta) &
                             +(1.0d0-derf((tb+2.0d0*veta*veta*ta)/(2.0d0*veta)))/dexp(-tb*ta))
            vep1=-pi/surf*((1.0d0-derf((tb-2.0d0*veta*veta*ta)/(2.0d0*veta)))*dexp(-tb*ta) &
                          -(1.0d0-derf((tb+2.0d0*veta*veta*ta)/(2.0d0*veta)))*dexp(tb*ta))
            vlo0=dsin(vkx*(x0-x1)+vky*(y0-y1))*vep0
            vlo1=dcos(vkx*(x0-x1)+vky*(y0-y1))*vep1
            fatxx=fatxx-vkx*cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0
            fatyy=fatyy-vky*cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0
            fatzz=fatzz-cp(1,indspe(l0))*cp(1,indspe(l1))*vlo1
          end if
        end do
        end do
        z=dabs(atz(l0)-atz(l1))
        if (z .gt. 1.0d-20) then
          vlo1=-2.0d0*dsqrt(pi)/surf*cp(1,indspe(l0))*cp(1,indspe(l1))*(atz(l0)-atz(l1))/z
          fatzz=fatzz-vlo1*dsqrt(pi)*derf(veta*z)
        end if
        do ky=-new_rsy,new_rsy
        do kx=-new_rsx,new_rsx
          x=kx*2.0d0*xmax+x0-atx(l1)
          y=ky*2.0d0*ymax+y0-aty(l1)
          z=z0-atz(l1)
          r=dsqrt(x*x+y*y+z*z)
          if ((r .gt. 1.0d-20) .and.(veta**2*r**2 .lt. 100.0d0))then
            vlo0=-2.0d0*veta/(dexp(veta**2*r**2)*dsqrt(pi)*r)-(1.0d0-derf(veta*r))/r**2
            fatxx=fatxx-cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0*x/r
            fatyy=fatyy-cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0*y/r
            fatzz=fatzz-cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0*z/r
          end if
        end do
        end do
      case (3)
        do kz=-new_pwz,new_pwz
        do ky=-new_pwy,new_pwy
        do kx=-new_pwx,new_pwx
          if (kx**2+ky**2+kz**2 .ne. 0) then
            vkx=pi/xmax*kx
            vky=pi/ymax*ky
            vkz=pi/zmax*kz
            vk2=vkx**2+vky**2+vkz**2
            if (vk2/(4.0d0*veta**2) .lt. 100.0d0) then
              vlo0=-fourpi*omegain*dsin(vkx*(x0-x1)+vky*(y0-y1)+vkz*(z0-z1))/vk2*dexp(-vk2/(4.0d0*veta**2))
            else
              vlo0=0.0d0
            end if
            fatxx=fatxx-cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0*vkx
            fatyy=fatyy-cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0*vky
            fatzz=fatzz-cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0*vkz
          end if
        end do
        end do
        end do
        do kz=-new_rsz,new_rsz
        do ky=-new_rsy,new_rsy
        do kx=-new_rsx,new_rsx
          x=kx*2.0d0*xmax+x0-atx(l1)
          y=ky*2.0d0*ymax+y0-aty(l1)
          z=kz*2.0d0*zmax+z0-atz(l1)
          r=dsqrt(x*x+y*y+z*z)
          if ((r .gt. 1.0d-20) .and.(veta**2*r**2 .lt. 100.0d0))then
            vlo0=-2.0d0*veta/(dexp(veta**2*r**2)*dsqrt(pi)*r)-(1.0d0-derf(veta*r))/r**2
            fatxx=fatxx-cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0*x/r
            fatyy=fatyy-cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0*y/r
            fatzz=fatzz-cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0*z/r
          end if
        end do
        end do
        end do
      end select
    end do
!$omp critical
  fatx(l0)=fatx(l0)+fatxx
  faty(l0)=faty(l0)+fatyy
  fatz(l0)=fatz(l0)+fatzz
!$omp end critical
!$omp barrier
  end do
!$omp end parallel
! -----------------------------------------------

  call mpi_allreduce(fatx,fatall,natom,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  fatx=fatall
  call mpi_allreduce(faty,fatall,natom,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  faty=fatall
  call mpi_allreduce(fatz,fatall,natom,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  fatz=fatall

! ----------  jellium part II  ----------
  if (jelcalc==key_jel_calc) then
    do na=1,natom
      if (natpri(na) .eq. key_natpri_in) then
        ipri=natpri_inf(na)
        do kz=-nzmax,nzmax
          if (kz**2 .ne. 0) then
            vkz=pi/zmax*kz
            vk2=vkz*vkz
            do ns=1,min(2,nspv)
              do il=1,npoint
                do ir=2,nradct(indspe(na))
                  z=point(il,3)*radial(ir,indspe(na))+atz(na)
                  r=radial(ir,indspe(na))
                  dr=dradial(ir,indspe(na))
                  fjel=fjel+4.0d0*pi/omega*chrjel/(strjel-endjel)/vk2*(dcos(vkz*(strjel-z))-dcos(vkz*(endjel-z))) &
                       *(rhotrur(ir,il,ns,ipri)-rhosmtr(ir,il,ns,ipri)-rhoaugr(ir,il,ipri))*r*r*dr*wt(il)
                       ! ?? nspv-dependent prefactor of rhoaugr needed?
                end do
              end do
            end do
          end if
        end do
      end if
    end do
    do kz=-nzmax,nzmax
      if (kz**2 .ne. 0) then
        vkz=pi/zmax*kz
        vk2=vkz*vkz
        fjel0=0.0d0
        do iz=1,ncpz_d
        do iy=1,ncpy_d
        do ix=1,ncpx_d
          jz=myrz*ncpz_d+iz
          z=jz*ddz-0.5d0*ddz-zmax
          fjel0=fjel0+4.0d0*pi/omega*chrjel/(strjel-endjel)/vk2*(dcos(vkz*(strjel-z))-dcos(vkz*(endjel-z))) &
                       *rho_aug_dense(ix,iy,iz)*ddx*ddy*ddz
        end do
        end do
        end do
        fjel=fjel+fjel0
      end if
    end do
    call mpi_allreduce(fjel,fjel0,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
    fjel=fjel0
  endif
! ---------------------------------------

  do na=1,natom
     fatx(na)=fatx(na)+biasx*cp(1,indspe(na))
     faty(na)=faty(na)+biasy*cp(1,indspe(na))
     fatz(na)=fatz(na)+biasz*cp(1,indspe(na))
  end do
  fxall=0.0d0
  fyall=0.0d0
  fzall=fjel
  do na=1,natom
     fxall=fxall+fatx(na)
     fyall=fyall+faty(na)
     fzall=fzall+fatz(na)
  end do
  do na=1,natom
     fatx(na)=fatx(na)-fxall/dble(natom)
     faty(na)=faty(na)-fyall/dble(natom)
     fatz(na)=fatz(na)-fzall/dble(natom)
  end do

  deallocate(fatall)

end subroutine force


subroutine force_01(natom,nmesh,na,ncpx,ncpy,ncpz,num_ppcell_d,num_list_d, & ! <
                    key_natpri_inps,                                       & ! <
                    dx,dy,dz,                                              & ! <
                    tmpx,tmpy,tmpz,                                        & ! X
                    natprid,napsd,natinfd,natinfd_vloc,lstvecd2,           & ! <
                    lstdx,lstdy,lstdz,                                     & ! <
                    rho_coarse,rho_aug_dense,vh_dense,                     & ! <
                    drhoaug3ddx,drhoaug3ddy,drhoaug3ddz,                   & ! <
                    dvlocdx_scw,dvlocdy_scw,dvlocdz_scw,                   & ! <
                    dvlocdx_hdp,dvlocdy_hdp,dvlocdz_hdp,                   & ! <
                    vloc_dense)                                              ! <
use mod_mpi
implicit none
integer,intent(in)::natom,nmesh,na,ncpx,ncpy,ncpz,num_ppcell_d,num_list_d
integer,intent(in)::key_natpri_inps
real*8,intent(in)::dx,dy,dz
real*8,intent(inout)::tmpx,tmpy,tmpz
integer,intent(in)::natprid(natom),napsd(natom),natinfd(natom),natinfd_vloc(natom)
integer,intent(in)::lstdx(num_list_d,num_ppcell_d),lstdy(num_list_d,num_ppcell_d),lstdz(num_list_d,num_ppcell_d)
integer,intent(in)::lstvecd2(num_list_d,num_ppcell_d)
real*8,intent(in)::rho_coarse(ncpx,ncpy,ncpz)
real*8,intent(in)::rho_aug_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
real*8,intent(in)::vh_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
real*8,intent(in)::vloc_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
real*8,intent(in)::dvlocdx_scw(ncpx,ncpy,ncpz,natom),dvlocdy_scw(ncpx,ncpy,ncpz,natom),dvlocdz_scw(ncpx,ncpy,ncpz,natom)
real*8,intent(in)::dvlocdx_hdp(num_list_d,num_ppcell_d)
real*8,intent(in)::dvlocdy_hdp(num_list_d,num_ppcell_d)
real*8,intent(in)::dvlocdz_hdp(num_list_d,num_ppcell_d)
real*8 drhoaug3ddx(num_list_d,num_ppcell_d)
real*8 drhoaug3ddy(num_list_d,num_ppcell_d)
real*8 drhoaug3ddz(num_list_d,num_ppcell_d)
integer ix,iy,iz,ixyz,i
real*8 tmpx0,tmpy0,tmpz0,dxyz,ddxyz
real*8 ddx,ddy,ddz

  ddx=dx/dble(nmesh)
  ddy=dy/dble(nmesh)
  ddz=dz/dble(nmesh)
  tmpx0=0.0d0
  tmpy0=0.0d0
  tmpz0=0.0d0
! ----------  core - electron attraction  ----------
  dxyz=dx*dy*dz
!$omp do
  do iz=1,ncpz
  do iy=1,ncpy
  do ix=1,ncpx
    tmpx0=tmpx0+dvlocdx_scw(ix,iy,iz,na)*rho_coarse(ix,iy,iz)*dxyz
    tmpy0=tmpy0+dvlocdy_scw(ix,iy,iz,na)*rho_coarse(ix,iy,iz)*dxyz
    tmpz0=tmpz0+dvlocdz_scw(ix,iy,iz,na)*rho_coarse(ix,iy,iz)*dxyz
  end do
  end do
  end do
  if (natprid(na) .eq. key_natpri_inps) then
    ddxyz=ddx*ddy*ddz
!$omp do
    do ixyz=1,natinfd_vloc(na)
      i=lstvecd2(ixyz,napsd(na))
      tmpx0=tmpx0+dvlocdx_hdp(ixyz,napsd(na))*rho_aug_dense(i,1,1)*ddxyz
      tmpy0=tmpy0+dvlocdy_hdp(ixyz,napsd(na))*rho_aug_dense(i,1,1)*ddxyz
      tmpz0=tmpz0+dvlocdz_hdp(ixyz,napsd(na))*rho_aug_dense(i,1,1)*ddxyz
    end do
  end if
! --------------------------------------------------

! ----------  compensation charge - core+electron interaction  ----------
  if (natprid(na) .eq. key_natpri_inps) then
    ddxyz=ddx*ddy*ddz
!$omp do
    do ixyz=1,natinfd(na)
      ix=lstdx(ixyz,napsd(na))
      iy=lstdy(ixyz,napsd(na))
      iz=lstdz(ixyz,napsd(na))
      i=lstvecd2(ixyz,napsd(na))
      tmpx0=tmpx0+(vloc_dense(i,1,1)+vh_dense(i,1,1))*drhoaug3ddx(ixyz,napsd(na))*ddxyz
      tmpy0=tmpy0+(vloc_dense(i,1,1)+vh_dense(i,1,1))*drhoaug3ddy(ixyz,napsd(na))*ddxyz
      tmpz0=tmpz0+(vloc_dense(i,1,1)+vh_dense(i,1,1))*drhoaug3ddz(ixyz,napsd(na))*ddxyz
    end do
  end if
! -----------------------------------------------------------------------
!$omp critical
  tmpx=tmpx+tmpx0
  tmpy=tmpy+tmpy0
  tmpz=tmpz+tmpz0
!$omp end critical
!$omp barrier

end subroutine force_01


subroutine force_02(natom,num_spe,nspv,nradmx,nqmx,ncpx_d,ncpy_d,ncpz_d, & ! <
                    xmax,ymax,zmax,ddx,ddy,ddz,psctoff,psftrad, &          ! <
                    nradct,indspe,nqctpcc, &                               ! <
                    fatx,faty,fatz, &                                      ! X
                    atx,aty,atz, &                                         ! <
                    radial,coef,vx_dense)                                  ! <
use mod_mpi
implicit none
integer, intent(in)::natom,num_spe,nspv,nradmx,nqmx
integer, intent(in)::ncpx_d,ncpy_d,ncpz_d
real*8, intent(in)::xmax,ymax,zmax
real*8, intent(in)::ddx,ddy,ddz
real*8, intent(in)::psctoff,psftrad
integer, intent(in)::nradct(num_spe),indspe(natom),nqctpcc(num_spe)
real*8, intent(inout)::fatx(natom),faty(natom),fatz(natom)
real*8, intent(in)::atx(natom),aty(natom),atz(natom)
real*8, intent(in)::radial(nradmx,num_spe),coef(0:nqmx,0:7,num_spe)
real*8, intent(in)::vx_dense(ncpx_d,ncpy_d,ncpz_d,nspv)
integer na,ns,ix,iy,iz,jx,jy,jz,kx,ky,kz,iq
real*8 x,y,z,r,dbesselj0,qqq,qqqr,qqqrin,dqq,rin,vcc,fourpi,twopicbin,tmp1,rcutss,vlo,ffatx,ffaty,ffatz
real*8 pi

  pi=dacos(-1.0d0)
  twopicbin=1.0d0/(2.0d0*pi)**3
  fourpi=4.0d0*pi

  do na=1,natom
    rcutss=radial(nradct(indspe(na)),indspe(na))*psftrad
    do ns=1,nspv
      do kz=-1,1
      do ky=-1,1
      do kx=-1,1
        ffatx=0.0d0
        ffaty=0.0d0
        ffatz=0.0d0
!$omp do
        do iz=1,ncpz_d
        do iy=1,ncpy_d
        do ix=1,ncpx_d
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=jx*ddx-xmax-0.5d0*ddx-atx(na)+kx*xmax*2.0d0
          y=jy*ddy-ymax-0.5d0*ddy-aty(na)+ky*ymax*2.0d0
          z=jz*ddz-zmax-0.5d0*ddz-atz(na)+kz*zmax*2.0d0
          r=dsqrt(x*x+y*y+z*z)
          tmp1=0.0d0
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
            dqq=pi/rcutss
            rin=1.0d0/r
            do iq=1,nqctpcc(indspe(na))
              qqq=dqq*iq
              qqqr=qqq*r
              qqqrin=rin/qqq
              dbesselj0=-dsin(qqqr)*qqqrin*rin+dcos(qqqr)*rin
              tmp1=tmp1+coef(iq,7,indspe(na))*dbesselj0*qqq*qqq
            end do
            vlo=tmp1*fourpi*dqq*twopicbin/dble(nspv)
            vcc=vx_dense(ix,iy,iz,ns)
            ffatx=ffatx+vlo*vcc*ddx*ddy*ddz*(x/r)
            ffaty=ffaty+vlo*vcc*ddx*ddy*ddz*(y/r)
            ffatz=ffatz+vlo*vcc*ddx*ddy*ddz*(z/r)
          end if
        end do
        end do
        end do
!$omp critical
        fatx(na)=fatx(na)+ffatx
        faty(na)=faty(na)+ffaty
        fatz(na)=fatz(na)+ffatz
!$omp end critical
!$omp barrier
      end do
      end do
      end do
    end do
  end do

end subroutine

! ========================================================================================

end module
