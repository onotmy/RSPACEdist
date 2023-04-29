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
! **********  scf_chargemixing8f.F90 06/20/2018-01  **********

module mod_scf_chargemixing
implicit none
real*8,allocatable::fvec_brydn(:,:,:,:,:),fvectrur_brydn(:,:,:,:,:),fvecsmtr_brydn(:,:,:,:,:)
real*8,allocatable::uvec_brydn(:,:,:,:,:),uvectrur_brydn(:,:,:,:,:),uvecsmtr_brydn(:,:,:,:,:)
real*8,allocatable::gvec_brydn(:,:,:,:)  ,gvectrur_brydn(:,:,:,:)  ,gvecsmtr_brydn(:,:,:,:) 
real*8,allocatable::rhosmt_brydn(:,:,:,:),rhotrur_brydn(:,:,:,:)   ,rhosmtr_brydn(:,:,:,:) 
real*8,allocatable::amat1(:,:),amat3(:)

contains

subroutine scf_chargemixing_initialize( &
 ncpx,ncpy,ncpz,nradmx,nspv,nbrydn,npoint,num_atcell) ! <
implicit none
integer,intent(in)::ncpx,ncpy,ncpz,nradmx,nspv,nbrydn,npoint,num_atcell
  if (nbrydn>1) then
    allocate(fvec_brydn(ncpx,ncpy,ncpz,nspv,nbrydn))
    allocate(fvectrur_brydn(nradmx,npoint,nspv,num_atcell,nbrydn))
    allocate(fvecsmtr_brydn(nradmx,npoint,nspv,num_atcell,nbrydn))
    allocate(uvec_brydn(ncpx,ncpy,ncpz,nspv,2:nbrydn))
    allocate(uvectrur_brydn(nradmx,npoint,nspv,num_atcell,2:nbrydn))
    allocate(uvecsmtr_brydn(nradmx,npoint,nspv,num_atcell,2:nbrydn))
    allocate(gvec_brydn(ncpx,ncpy,ncpz,nspv))
    allocate(gvectrur_brydn(nradmx,npoint,nspv,num_atcell))
    allocate(gvecsmtr_brydn(nradmx,npoint,nspv,num_atcell))
    allocate(rhosmt_brydn(ncpx,ncpy,ncpz,nspv))
    allocate(rhotrur_brydn(nradmx,npoint,nspv,num_atcell))
    allocate(rhosmtr_brydn(nradmx,npoint,nspv,num_atcell))
    allocate(amat1(2:nbrydn,2:nbrydn),amat3(2:nbrydn))
  else
    allocate(fvec_brydn(1,1,1,1,1))
    allocate(fvectrur_brydn(1,1,1,1,1))
    allocate(fvecsmtr_brydn(1,1,1,1,1))
    allocate(uvec_brydn(1,1,1,1,1))
    allocate(uvectrur_brydn(1,1,1,1,1))
    allocate(uvecsmtr_brydn(1,1,1,1,1))
    allocate(gvec_brydn(1,1,1,1))
    allocate(gvectrur_brydn(1,1,1,1))
    allocate(gvecsmtr_brydn(1,1,1,1))
    allocate(rhosmt_brydn(1,1,1,1))
    allocate(rhotrur_brydn(1,1,1,1))
    allocate(rhosmtr_brydn(1,1,1,1))
    allocate(amat1(1,1),amat3(1))
  endif 
end subroutine scf_chargemixing_initialize


subroutine scf_chargemixing_finalize
implicit none
  deallocate(fvec_brydn  ,fvectrur_brydn,fvecsmtr_brydn)
  deallocate(uvec_brydn  ,uvectrur_brydn,uvecsmtr_brydn)
  deallocate(gvec_brydn  ,gvectrur_brydn,gvecsmtr_brydn)
  deallocate(rhosmt_brydn,rhotrur_brydn ,rhosmtr_brydn)
  deallocate(amat1,amat3)
end subroutine scf_chargemixing_finalize


subroutine scf_chargemixing( &
 nomix,                                                      & ! <
 key_natpri_in,key_pp_paw,                                   & ! <
 nbrydn,ndisp,                                               & ! <
 natom,num_spe,num_atcell,nradmx,nspv,npoint,ncpx,ncpy,ncpz, & ! <
 ntyppp,nradct,indspe,natpri,natpri_inf,                     & ! <
 eta,etamag,dx,dy,dz,radial,dradial,wt,                      & ! <
 rhosmt_o,rhotrur_o,rhosmtr_o,                               & ! <
 lbrydn,                                                     & ! X
 rhosmt_i,rhotrur_i,rhosmtr_i)                                 ! X
use mod_mpi, only:myrank_glbl
use mod_stopp
implicit none
logical,parameter:: l_polinp= .false. ! write mag. mom. of input charge ? 
logical,intent(in)   ::nomix
integer,intent(in)   ::key_natpri_in,key_pp_paw
integer,intent(in)   ::nbrydn,ndisp
integer,intent(in)   ::natom,num_spe,num_atcell,nradmx,nspv,npoint
integer,intent(in)   ::ncpx,ncpy,ncpz
integer,intent(in)   ::ntyppp(num_spe),nradct(num_spe)
integer,intent(in)   ::indspe(natom),natpri(natom),natpri_inf(natom)
real*8, intent(in)   ::eta,etamag(2)
real*8, intent(in)   ::dx,dy,dz
real*8, intent(in)   ::radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8, intent(in)   ::wt(npoint)
real*8, intent(in)   ::rhosmt_o(ncpx,ncpy,ncpz,nspv)
real*8, intent(in)   ::rhotrur_o(nradmx,npoint,nspv,num_atcell),rhosmtr_o(nradmx,npoint,nspv,num_atcell)
integer,intent(inout)::lbrydn 
real*8, intent(inout)::rhosmt_i(ncpx,ncpy,ncpz,nspv)
real*8, intent(inout)::rhotrur_i(nradmx,npoint,nspv,num_atcell),rhosmtr_i(nradmx,npoint,nspv,num_atcell)
integer:: ns,ix,na,ipri,il,ir  
real*8 :: tmp,tmpall
real*8,allocatable::tmpmat1(:,:)

  lbrydn= lbrydn +1

  if (lbrydn==1) then

! **********  complete update of charge density for nloop=1  **********
    !$omp parallel default(shared) private(ix,il,ir,ipri)
    do ns=1,nspv
      !$omp do
      do ix=1,ncpx*ncpy*ncpz
        rhosmt_i(ix,1,1,ns)=rhosmt_o(ix,1,1,ns)
      end do
      do na=1,natom
        if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
          ipri=natpri_inf(na)
          !$omp do
          do il=1,npoint
            do ir=1,nradct(indspe(na))-1
              rhotrur_i(ir,il,ns,ipri)=rhotrur_o(ir,il,ns,ipri)
              rhosmtr_i(ir,il,ns,ipri)=rhosmtr_o(ir,il,ns,ipri)
            end do
          end do
        end if
      end do
    end do
    !$omp end parallel
! *********************************************************************

  else

    if ((.not. nomix).and.(nbrydn>0)) then
      if (nbrydn==1) then
        if (myrank_glbl .eq. 0) write(ndisp,*,err=9999) 'mixing ==> straight'
        !$omp parallel default(shared)
        call scf_chargemixing_straight( &
         key_natpri_in,key_pp_paw,                                & ! <
         nspv,                                                    & ! <
         natom,num_spe,num_atcell,nradmx,npoint,ncpx,ncpy,ncpz,   & ! <
         ntyppp,nradct,indspe,natpri,natpri_inf,                  & ! <
         eta,etamag,                                              & ! <
         rhosmt_o,rhotrur_o,rhosmtr_o,                            & ! <
         rhosmt_i,rhotrur_i,rhosmtr_i)                              ! X
        !$omp end parallel
      else
        if (myrank_glbl .eq. 0) write(ndisp,*,err=9999) 'mixing ==> broyden'
        allocate(tmpmat1(2:nbrydn,2:nbrydn))
        if (lbrydn==2) then
          amat1=0.0d0
          amat3=0.0d0
        endif
        !$omp parallel default(shared)
        call scf_chargemixing_broyden( &
         key_natpri_in,key_pp_paw,                                            & ! <
         nspv,                                                                & ! <
         natom,num_spe,num_atcell,nradmx,npoint,nbrydn,lbrydn,ncpx,ncpy,ncpz, & ! <
         ntyppp,nradct,indspe,natpri,natpri_inf,                              & ! <
         eta,etamag,dx,dy,dz,radial,dradial,wt,                               & ! <
         rhosmt_o,rhotrur_o,rhosmtr_o,                                        & ! <
         rhosmt_i,rhotrur_i,rhosmtr_i,                                        & ! X
         tmp,tmpall,tmpmat1)                                                    ! 0
        !$omp end parallel
        deallocate(tmpmat1)
      end if
    endif ! ((.not. nomix).and.(nbrydn>0))

  endif ! (lbrydn==1) else

  ! spin polarization of input charge
  if (l_polinp) call scf_chargemixing_spinpol( &
   '',ndisp,key_pp_paw,key_natpri_in,                        & ! <
   natom,num_spe,num_atcell,nspv,nradmx,npoint,              & ! <
   natpri,natpri_inf,indspe,ntyppp,nradct,wt,radial,dradial, & ! <
   rhotrur_i)                                                  ! <

  return
 9999 call stopp('scf_chargemixing: error writing ndisp')
end subroutine scf_chargemixing


subroutine scf_chargemixing_straight( &
 key_natpri_in,key_pp_paw,                              & ! <
 nspv,                                                  & ! <
 natom,num_spe,num_atcell,nradmx,npoint,ncpx,ncpy,ncpz, & ! <
 ntyppp,nradct,indspe,natpri,natpri_inf,                & ! <
 eta,etamag,                                            & ! <
 rhosmt_o,rhotrur_o,rhosmtr_o,                          & ! <
 rhosmt_i,rhotrur_i,rhosmtr_i)                            ! X
use mod_stopp 
implicit none
integer,intent(in)   ::key_natpri_in,key_pp_paw
integer,intent(in)   ::nspv
integer,intent(in)   ::natom,num_spe,num_atcell,nradmx,npoint,ncpx,ncpy,ncpz
integer,intent(in)   ::ntyppp(num_spe),nradct(num_spe)
integer,intent(in)   ::indspe(natom),natpri(natom),natpri_inf(natom)
real*8, intent(in)   ::eta,etamag(2)
real*8, intent(in)   ::rhosmt_o(ncpx,ncpy,ncpz,nspv)
real*8, intent(in)   ::rhotrur_o(nradmx,npoint,nspv,num_atcell),rhosmtr_o(nradmx,npoint,nspv,num_atcell)
real*8, intent(inout)::rhosmt_i(ncpx,ncpy,ncpz,nspv)
real*8, intent(inout)::rhotrur_i(nradmx,npoint,nspv,num_atcell),rhosmtr_i(nradmx,npoint,nspv,num_atcell)
integer:: ipri,na

  call scf_chargemixing_straight_mix( &
   ncpx*ncpy*ncpz,ncpx*ncpy*ncpz,1,nspv,eta,etamag, & ! <
   rhosmt_o,                                        & ! <
   rhosmt_i)                                          ! X
  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
      ipri = natpri_inf(na)
      call scf_chargemixing_straight_mix( &
       nradmx,nradct(indspe(na))-1,npoint,nspv,eta,etamag, & ! <
       rhotrur_o(1,1,1,ipri),                              & ! <
       rhotrur_i(1,1,1,ipri))                                ! X
      call scf_chargemixing_straight_mix( &
       nradmx,nradct(indspe(na))-1,npoint,nspv,eta,etamag, & ! <
       rhosmtr_o(1,1,1,ipri),                              & ! <
       rhosmtr_i(1,1,1,ipri))                                ! X
    endif
  enddo

end subroutine scf_chargemixing_straight


subroutine scf_chargemixing_broyden( &
 key_natpri_in,key_pp_paw,                                            & ! <
 nspv,                                                                & ! <
 natom,num_spe,num_atcell,nradmx,npoint,nbrydn,lbrydn,ncpx,ncpy,ncpz, & ! <
 ntyppp,nradct,indspe,natpri,natpri_inf,                              & ! <
 eta,etamag,dx,dy,dz,radial,dradial,wt,                               & ! <
 rhosmt_o,rhotrur_o,rhosmtr_o,                                        & ! <
 rhosmt,rhotrur,rhosmtr,                                              & ! X
 tmp,tmpall,tmpmat1)                                                    ! 0
use mod_mpi
use mod_stopp 
implicit none
integer,intent(in)   ::key_natpri_in,key_pp_paw
integer,intent(in)   ::nspv
integer,intent(in)   ::natom,num_spe,num_atcell,nradmx,npoint,nbrydn,lbrydn 
integer,intent(in)   ::ncpx,ncpy,ncpz
integer,intent(in)   ::ntyppp(num_spe),nradct(num_spe)
integer,intent(in)   ::indspe(natom),natpri(natom),natpri_inf(natom)
real*8, intent(in)   ::eta,etamag(2)
real*8, intent(in)   ::dx,dy,dz
real*8, intent(in)   ::radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8, intent(in)   ::wt(npoint)
real*8, intent(in)   ::rhosmt_o(ncpx,ncpy,ncpz,nspv)
real*8, intent(in)   ::rhotrur_o(nradmx,npoint,nspv,num_atcell),rhosmtr_o(nradmx,npoint,nspv,num_atcell)
real*8, intent(inout)::rhosmt(ncpx,ncpy,ncpz,nspv)
real*8, intent(inout)::rhotrur(nradmx,npoint,nspv,num_atcell),rhosmtr(nradmx,npoint,nspv,num_atcell)
! work variables
real*8, intent(inout)::tmp,tmpall
real*8, intent(inout)::tmpmat1(2:nbrydn,2:nbrydn)
real*8 :: tmp0,r,dr
integer:: ns,i,j,ix,ir,il,ipri,na,mbrydn
real*8,allocatable:: pol_o(:)
real*8,allocatable:: tmpmat2(:)
real*8,allocatable:: bvecsmt(:,:,:,:),bvectrur(:,:,:,:),bvecsmtr(:,:,:,:) 

! **********  Broyden mixing according to PRB34, 8391 (1986)  **********

  mbrydn= min( lbrydn-1, nbrydn ) 

  allocate(pol_o(3)) 
  allocate(tmpmat2(2:nbrydn)) 
  allocate(bvecsmt(ncpx,ncpy,ncpz,nspv),bvectrur(nradmx,npoint,nspv,num_atcell),bvecsmtr(nradmx,npoint,nspv,num_atcell))

! ==========  U_m <- G_1 F_{m-1}  ==========
  if (mbrydn>1) then
    do ns=1,nspv
      !$omp do
      do ix=1,ncpx*ncpy*ncpz
        uvec_brydn(ix,1,1,ns,mbrydn)= -gvec_brydn(ix,1,1,ns)
      end do
    end do
    do na=1,natom
      if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
        ipri=natpri_inf(na)
        do ns=1,nspv
          !$omp do
          do il=1,npoint
            do ir=1,nradct(indspe(na))-1
              uvectrur_brydn(ir,il,ns,ipri,mbrydn)= -gvectrur_brydn(ir,il,ns,ipri)
              uvecsmtr_brydn(ir,il,ns,ipri,mbrydn)= -gvecsmtr_brydn(ir,il,ns,ipri)
            end do
          end do
        end do
      end if
    end do
  end if
! ==========================================

! ==========  compute F in Eq.(3) of PRB34,8391(1986), -G_1 F, and m/|m|  ========
  call scf_chargemixing_broyden_vecs( &
   ncpx*ncpy*ncpz,ncpx*ncpy*ncpz,1,nspv,eta,etamag, & ! <
   rhosmt_o,                                        & ! <
   rhosmt,                                          & ! X
   gvec_brydn,fvec_brydn(1,1,1,1,mbrydn))             ! >
  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
      ipri = natpri_inf(na)
      call scf_chargemixing_broyden_vecs( &
       nradmx,nradct(indspe(na))-1,npoint,nspv,eta,etamag,                               & ! <
       rhotrur_o(1,1,1,ipri),                                                            & ! <
       rhotrur(1,1,1,ipri),                                                              & ! X
       gvectrur_brydn(1,1,1,ipri),fvectrur_brydn(1,1,1,ipri,mbrydn))                       ! >
      call scf_chargemixing_broyden_vecs( &
       nradmx,nradct(indspe(na))-1,npoint,nspv,eta,etamag,                               & ! <
       rhosmtr_o(1,1,1,ipri),                                                            & ! <
       rhosmtr(1,1,1,ipri),                                                              & ! X
       gvecsmtr_brydn(1,1,1,ipri),fvecsmtr_brydn(1,1,1,ipri,mbrydn))                       ! >
    end if
  end do
! ================================================================================

  if (mbrydn>1) then 

! ==========  compute V^{T(j)}(F^{(i)}-F^{(i-1)}), V^{T(j)}(F^{(m)}) in Eqs. (7), (6) of PRB34,8391 (1986)  ==========
    !$omp single
    tmpmat1=amat1
    !$omp end single
    do i=2,mbrydn
      !$omp single
      tmp=0.0d0
      !$omp end single
      tmp0=0.0d0
      do ns=1,nspv
        !$omp do
        do ix=1,ncpx*ncpy*ncpz
          tmp0=tmp0+(fvec_brydn(ix,1,1,ns,mbrydn)-fvec_brydn(ix,1,1,ns,mbrydn-1)) &
                   *(fvec_brydn(ix,1,1,ns,     i)-fvec_brydn(ix,1,1,ns,     i-1))*dx*dy*dz
        end do
        do na=1,natom
          if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
            ipri=natpri_inf(na)
            !$omp do
            do il=1,npoint
              do ir=1,nradct(indspe(na))-1
                r=radial(ir,indspe(na))
                dr=dradial(ir,indspe(na))
                tmp0=tmp0+(fvectrur_brydn(ir,il,ns,ipri,mbrydn)-fvectrur_brydn(ir,il,ns,ipri,mbrydn-1)) &
                         *(fvectrur_brydn(ir,il,ns,ipri,     i)-fvectrur_brydn(ir,il,ns,ipri,     i-1))*wt(il)*r*r*dr
                tmp0=tmp0+(fvecsmtr_brydn(ir,il,ns,ipri,mbrydn)-fvecsmtr_brydn(ir,il,ns,ipri,mbrydn-1)) &
                         *(fvecsmtr_brydn(ir,il,ns,ipri,     i)-fvecsmtr_brydn(ir,il,ns,ipri,     i-1))*wt(il)*r*r*dr
              end do
            end do
          end if
        end do
      end do
      !$omp critical
      tmp=tmp+tmp0
      !$omp end critical
      !$omp barrier
      !$omp single
      call mpi_allreduce(tmp,tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
      tmpmat1(mbrydn,i)=tmpall
      tmpmat1(i,mbrydn)=tmpall
      amat1(mbrydn,i)=tmpall
      amat1(i,mbrydn)=tmpall
      !$omp end single
    end do

    do i=2,mbrydn
      !$omp single
      tmp=0.0d0
      !$omp end single
      tmp0=0.0d0
      do ns=1,nspv
        !$omp do
        do ix=1,ncpx*ncpy*ncpz
          tmp0=tmp0+(fvec_brydn(ix,1,1,ns,i)-fvec_brydn(ix,1,1,ns,i-1))*fvec_brydn(ix,1,1,ns,mbrydn)*dx*dy*dz
        end do
        do na=1,natom
          if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
            ipri=natpri_inf(na)
            !$omp do
            do il=1,npoint
              do ir=1,nradct(indspe(na))-1
                r=radial(ir,indspe(na))
                dr=dradial(ir,indspe(na))
                tmp0=tmp0+(fvectrur_brydn(ir,il,ns,ipri,i)-fvectrur_brydn(ir,il,ns,ipri,i-1)) &
                          *fvectrur_brydn(ir,il,ns,ipri,mbrydn)*wt(il)*r*r*dr
                tmp0=tmp0+(fvecsmtr_brydn(ir,il,ns,ipri,i)-fvecsmtr_brydn(ir,il,ns,ipri,i-1)) &
                          *fvecsmtr_brydn(ir,il,ns,ipri,mbrydn)*wt(il)*r*r*dr
              end do
            end do
          end if
        end do
      end do
      !$omp critical
      tmp=tmp+tmp0
      !$omp end critical
      !$omp barrier
      !$omp single
      call mpi_allreduce(tmp,tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
      !$omp end single
      tmpmat2(i)=tmpall
    end do

    !$omp single
    tmp=0.0d0
    !$omp end single
    tmp0=0.0d0
    do ns=1,nspv
      !$omp do
      do ix=1,ncpx*ncpy*ncpz
        tmp0=tmp0+(fvec_brydn(ix,1,1,ns,mbrydn)-fvec_brydn(ix,1,1,ns,mbrydn-1)) &
                 *(fvec_brydn(ix,1,1,ns,mbrydn)-fvec_brydn(ix,1,1,ns,mbrydn-1))*dx*dy*dz
      end do
      do na=1,natom
        if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
          ipri=natpri_inf(na)
          !$omp do
          do il=1,npoint
            do ir=1,nradct(indspe(na))-1
              r=radial(ir,indspe(na))
              dr=dradial(ir,indspe(na))
              tmp0=tmp0+(fvectrur_brydn(ir,il,ns,ipri,mbrydn)-fvectrur_brydn(ir,il,ns,ipri,mbrydn-1)) &
                       *(fvectrur_brydn(ir,il,ns,ipri,mbrydn)-fvectrur_brydn(ir,il,ns,ipri,mbrydn-1))*wt(il)*r*r*dr
              tmp0=tmp0+(fvecsmtr_brydn(ir,il,ns,ipri,mbrydn)-fvecsmtr_brydn(ir,il,ns,ipri,mbrydn-1)) &
                       *(fvecsmtr_brydn(ir,il,ns,ipri,mbrydn)-fvecsmtr_brydn(ir,il,ns,ipri,mbrydn-1))*wt(il)*r*r*dr
            end do
          end do
        end if
      end do
    end do
    !$omp critical
    tmp=tmp+tmp0
    !$omp end critical
    !$omp barrier
    !$omp single
    call mpi_allreduce(tmp,tmpall,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
    amat3(mbrydn)=tmpall
    !$omp end single
    do j=2,mbrydn
      !$omp do
      do i=2,mbrydn
        tmpmat1(j,i)=tmpmat1(j,i)/amat3(j)
      end do
      tmpmat2(j)=tmpmat2(j)/amat3(j)
    end do
! ====================================================================================================================

! ==========  compute U in Eq. (7) of PRB34, 8391 (1986)  ==========
    do ns=1,nspv
      !$omp do
      do ix=1,ncpx*ncpy*ncpz
        uvec_brydn(ix,1,1,ns,mbrydn)= &
         uvec_brydn(ix,1,1,ns,mbrydn) +gvec_brydn(ix,1,1,ns) &
         +rhosmt(ix,1,1,ns) -rhosmt_brydn(ix,1,1,ns)
      end do
    end do
    do j=2,mbrydn-1
      do ns=1,nspv
        !$omp do
        do ix=1,ncpx*ncpy*ncpz
          uvec_brydn(ix,1,1,ns,mbrydn)=uvec_brydn(ix,1,1,ns,mbrydn)-tmpmat1(j,mbrydn)*uvec_brydn(ix,1,1,ns,j)
        end do
      end do
    end do
    do na=1,natom
      if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
        ipri=natpri_inf(na)
        if (mbrydn .ge. 2) then
          do ns=1,nspv
            !$omp do
            do il=1,npoint
              do ir=1,nradct(indspe(na))-1
                uvectrur_brydn(ir,il,ns,ipri,mbrydn)= &
                 uvectrur_brydn(ir,il,ns,ipri,mbrydn) +gvectrur_brydn(ir,il,ns,ipri) &
                 +rhotrur(ir,il,ns,ipri) -rhotrur_brydn(ir,il,ns,ipri)
                uvecsmtr_brydn(ir,il,ns,ipri,mbrydn)= &
                 uvecsmtr_brydn(ir,il,ns,ipri,mbrydn) +gvecsmtr_brydn(ir,il,ns,ipri) &
                 +rhosmtr(ir,il,ns,ipri) -rhosmtr_brydn(ir,il,ns,ipri)
              end do
            end do
          end do
        end if
        do j=2,mbrydn-1
          do ns=1,nspv
            !$omp do
            do il=1,npoint
              do ir=1,nradct(indspe(na))-1
                uvectrur_brydn(ir,il,ns,ipri,mbrydn)= &
                 uvectrur_brydn(ir,il,ns,ipri,mbrydn)-tmpmat1(j,mbrydn)*uvectrur_brydn(ir,il,ns,ipri,j)
                uvecsmtr_brydn(ir,il,ns,ipri,mbrydn)= &
                 uvecsmtr_brydn(ir,il,ns,ipri,mbrydn)-tmpmat1(j,mbrydn)*uvecsmtr_brydn(ir,il,ns,ipri,j)
              end do
            end do
          end do
        end do
      end if
    end do
    !$omp barrier
! ==================================================================

  end if ! (mbrydn>1) 

! ==========  save old rho and compute Broyden correction  ==========
  do ns=1,nspv
    !$omp do
    do ix=1,ncpx*ncpy*ncpz
      rhosmt_brydn(ix,1,1,ns)= rhosmt(ix,1,1,ns)
      bvecsmt(ix,1,1,ns)= gvec_brydn(ix,1,1,ns)
    end do
  end do
  do ns=1,nspv
    do j=2,mbrydn
      !$omp do
      do ix=1,ncpx*ncpy*ncpz
        bvecsmt(ix,1,1,ns)= bvecsmt(ix,1,1,ns)-tmpmat2(j)*uvec_brydn(ix,1,1,ns,j)
      end do
    end do
  end do
  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
      ipri=natpri_inf(na)
      do ns=1,nspv
        !$omp do
        do il=1,npoint
          do ir=1,nradct(indspe(na))-1
            rhotrur_brydn(ir,il,ns,ipri)= rhotrur(ir,il,ns,ipri)
            rhosmtr_brydn(ir,il,ns,ipri)= rhosmtr(ir,il,ns,ipri)
            bvectrur(ir,il,ns,ipri)= gvectrur_brydn(ir,il,ns,ipri)
            bvecsmtr(ir,il,ns,ipri)= gvecsmtr_brydn(ir,il,ns,ipri)
          end do
        end do
        do j=2,mbrydn
          !$omp do
          do il=1,npoint
            do ir=1,nradct(indspe(na))-1
              bvectrur(ir,il,ns,ipri)= bvectrur(ir,il,ns,ipri)-tmpmat2(j)*uvectrur_brydn(ir,il,ns,ipri,j)
              bvecsmtr(ir,il,ns,ipri)= bvecsmtr(ir,il,ns,ipri)-tmpmat2(j)*uvecsmtr_brydn(ir,il,ns,ipri,j)
            end do
          end do
        end do
      end do ! ns
    end if
  end do
! ===================================================================

! ==========  compute new charge density  ==========
  do ns= 1,nspv
    !$omp do
    do ix=1,ncpx*ncpy*ncpz
      rhosmt(ix,1,1,ns)= rhosmt(ix,1,1,ns)+bvecsmt(ix,1,1,ns)  
    end do
  end do
  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
      ipri=natpri_inf(na)
      do ns=1,nspv
        !$omp do
        do il=1,npoint
          do ir=1,nradct(indspe(na))-1
            rhotrur(ir,il,ns,ipri)= rhotrur(ir,il,ns,ipri)+bvectrur(ir,il,ns,ipri)  
            rhosmtr(ir,il,ns,ipri)= rhosmtr(ir,il,ns,ipri)+bvecsmtr(ir,il,ns,ipri)  
          end do
        end do
      end do
    end if
  end do
! ==================================================

! ==========  update history of charge density  ==========
  if (mbrydn .eq. nbrydn) then
    tmpmat2(2)=1.0d0
    do i=3,mbrydn
      tmpmat2(i)=0.0d0
      do j=2,i-1
        tmpmat2(i)=tmpmat2(i)-tmpmat1(j,i)*tmpmat2(j)
      end do
      do ns=1,nspv
        !$omp do
        do ix=1,ncpx*ncpy*ncpz
          uvec_brydn(ix,1,1,ns,i)=uvec_brydn(ix,1,1,ns,i)-tmpmat2(i)*uvec_brydn(ix,1,1,ns,2)
        end do
      end do
    end do
    do i=2,mbrydn
      do ns=1,nspv
        !$omp do
        do ix=1,ncpx*ncpy*ncpz
          fvec_brydn(ix,1,1,ns,i-1)=fvec_brydn(ix,1,1,ns,i)
        end do
        if (i>2) then 
          !$omp do
          do ix=1,ncpx*ncpy*ncpz
            uvec_brydn(ix,1,1,ns,i-1)=uvec_brydn(ix,1,1,ns,i) 
          end do
        end if
      end do
    end do
    !$omp single
    do j=2,mbrydn-1
      do i=2,mbrydn-1
        amat1(i,j)=amat1(i+1,j+1)
      end do
      amat3(j)=amat3(j+1)
    end do
    !$omp end single
    do na=1,natom
      if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
        ipri=natpri_inf(na)
        tmpmat2(2)=1.0d0
        do i=3,mbrydn
          tmpmat2(i)=0.0d0
          do j=2,i-1
            tmpmat2(i)=tmpmat2(i)-tmpmat1(j,i)*tmpmat2(j)
          end do
          do ns=1,nspv
            !$omp do
            do il=1,npoint
              do ir=1,nradct(indspe(na))-1
                uvectrur_brydn(ir,il,ns,ipri,i)=uvectrur_brydn(ir,il,ns,ipri,i)-tmpmat2(i)*uvectrur_brydn(ir,il,ns,ipri,2)
                uvecsmtr_brydn(ir,il,ns,ipri,i)=uvecsmtr_brydn(ir,il,ns,ipri,i)-tmpmat2(i)*uvecsmtr_brydn(ir,il,ns,ipri,2)
              end do
            end do
          end do
        end do
        do i=2,mbrydn
          do ns=1,nspv
            !$omp do
            do il=1,npoint
              do ir=1,nradct(indspe(na))-1
                fvectrur_brydn(ir,il,ns,ipri,i-1)=fvectrur_brydn(ir,il,ns,ipri,i)
                fvecsmtr_brydn(ir,il,ns,ipri,i-1)=fvecsmtr_brydn(ir,il,ns,ipri,i)
              end do
            end do
            if (i>2) then
              !$omp do
              do il=1,npoint
                do ir=1,nradct(indspe(na))-1
                  uvectrur_brydn(ir,il,ns,ipri,i-1)=uvectrur_brydn(ir,il,ns,ipri,i)
                  uvecsmtr_brydn(ir,il,ns,ipri,i-1)=uvecsmtr_brydn(ir,il,ns,ipri,i)
                end do
              end do
            end if
          end do
        end do
        !$omp barrier
      end if
    end do
  end if
! ========================================================

! ==========  convert rho to ususal format  ==========
  call scf_chargemixing_broyden_rho( &
   ncpx*ncpy*ncpz,ncpx*ncpy*ncpz,1,nspv, & ! <
   rhosmt)                                 ! X
  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
      ipri = natpri_inf(na)
      call scf_chargemixing_broyden_rho( &
       nradmx,nradct(indspe(na))-1,npoint,nspv, & ! <
       rhotrur(1,1,1,ipri))                       ! X
      call scf_chargemixing_broyden_rho( &
       nradmx,nradct(indspe(na))-1,npoint,nspv, & ! <
       rhosmtr(1,1,1,ipri))                       ! X
    end if
  end do
! ====================================================

  deallocate(pol_o)
  deallocate(tmpmat2)
  deallocate(bvecsmt,bvectrur,bvecsmtr)

end subroutine scf_chargemixing_broyden


subroutine scf_chargemixing_straight_mix( &
 ndim1,nmax1,nmax2,nspv,eta,etamag, & ! <
 rho_o,                             & ! <
 rho_i)                               ! X
implicit none
integer,intent(in)   :: ndim1,nmax1,nmax2,nspv
real*8, intent(in)   :: eta,etamag(2)
real*8, intent(in)   :: rho_o(ndim1,nmax2,nspv)
real*8, intent(inout):: rho_i(ndim1,nmax2,nspv)
integer:: i1,i2 
real*8 :: chrg_i,chrg_o,chrg
real*8,allocatable:: pol_i(:),pol_o(:),pol(:)

  allocate(pol_i(3),pol_o(3),pol(3))

  select case (nspv) 
  case(1) 
    !$omp do
    do i2=1,nmax2
    do i1=1,nmax1
      rho_i(i1,i2,1)= (1.0d0-eta)*rho_i(i1,i2,1) +eta*rho_o(i1,i2,1)
    enddo
    enddo 
  case(2) 
    !$omp do
    do i2=1,nmax2
    do i1=1,nmax1
      chrg_i  = rho_i(i1,i2,1)+rho_i(i1,i2,nspv)
      pol_i(1)= rho_i(i1,i2,1)-rho_i(i1,i2,nspv)
      chrg_o  = rho_o(i1,i2,1)+rho_o(i1,i2,nspv)
      pol_o(1)= rho_o(i1,i2,1)-rho_o(i1,i2,nspv)
      chrg  = (1.0d0-eta      )*chrg_i   +eta      *chrg_o
      pol(1)= (1.0d0-etamag(1))*pol_i(1) +etamag(1)*pol_o(1)
      rho_i(i1,i2,1   )= 0.5d0*(chrg+pol(1))
      rho_i(i1,i2,nspv)= 0.5d0*(chrg-pol(1))
    enddo
    enddo
  case(4) 
    !$omp do
    do i2=1,nmax2
    do i1=1,nmax1
      chrg_i  = rho_i(i1,i2,1)+rho_i(i1,i2,nspv-2)
      pol_i(1)= rho_i(i1,i2,1)-rho_i(i1,i2,nspv-2)
      pol_i(2)= 2.0d0*rho_i(i1,i2,nspv-1)
      pol_i(3)= 2.0d0*rho_i(i1,i2,nspv-0)
      chrg_o  = rho_o(i1,i2,1)+rho_o(i1,i2,nspv-2)
      pol_o(1)= rho_o(i1,i2,1)-rho_o(i1,i2,nspv-2)
      pol_o(2)= 2.0d0*rho_o(i1,i2,nspv-1)
      pol_o(3)= 2.0d0*rho_o(i1,i2,nspv-0)
      chrg  = (1.0d0-eta      )*chrg_i   +eta      *chrg_o
      pol(1)= (1.0d0-etamag(1))*pol_i(1) +etamag(1)*pol_o(1)
      pol(2)= (1.0d0-etamag(1))*pol_i(2) +etamag(1)*pol_o(2)
      pol(3)= (1.0d0-etamag(1))*pol_i(3) +etamag(1)*pol_o(3)
      rho_i(i1,i2,nspv-3)= 0.5d0*(chrg+pol(1))
      rho_i(i1,i2,nspv-2)= 0.5d0*(chrg-pol(1))
      rho_i(i1,i2,nspv-1)= 0.5d0*pol(2)
      rho_i(i1,i2,nspv-0)= 0.5d0*pol(3)
    enddo
    enddo
  end select 

  deallocate(pol_i,pol_o,pol)

end subroutine scf_chargemixing_straight_mix


subroutine scf_chargemixing_broyden_vecs( &
 ndim1,nmax1,nmax2,nspv,eta,etamag, & ! <
 rho_o,                             & ! <
 rho_i,                             & ! X
 gvec_brydn,fvec_brydn)           ! >
implicit none
integer,intent(in)   :: ndim1,nmax1,nmax2,nspv
real*8, intent(in)   :: eta,etamag(2)
real*8, intent(in)   :: rho_o(ndim1,nmax2,nspv) 
real*8, intent(inout):: rho_i(ndim1,nmax2,nspv)
real*8, intent(out)  :: gvec_brydn(ndim1,nmax2,nspv),fvec_brydn(ndim1,nmax2,nspv)
integer:: i1,i2 
real*8 :: chrg_o,chrg_i
real*8,allocatable:: pol_i(:),pol_o(:),pol_n(:)

  allocate(pol_i(3),pol_o(3),pol_n(3))

  select case (nspv) 
  case(1)
    do i2=1,nmax2
    !$omp do
    do i1=1,nmax1
      fvec_brydn(i1,i2,1)= (rho_o(i1,i2,1)-rho_i(i1,i2,1)) 
      gvec_brydn(i1,i2,1)= (rho_o(i1,i2,1)-rho_i(i1,i2,1))*eta 
    end do
    end do
  case(2)
    do i2=1,nmax2
    !$omp do
    do i1=1,nmax1
      rho_i(i1,i2,1   )= rho_i(i1,i2,1)+rho_i(i1,i2,nspv)
      rho_i(i1,i2,nspv)= rho_i(i1,i2,1)-2.0d0*rho_i(i1,i2,nspv)
      chrg_o  = rho_o(i1,i2,1)+rho_o(i1,i2,nspv)
      pol_o(1)= rho_o(i1,i2,1)-rho_o(i1,i2,nspv)
      fvec_brydn(i1,i2,1   )= (chrg_o  -rho_i(i1,i2,1   ))
      fvec_brydn(i1,i2,nspv)= (pol_o(1)-rho_i(i1,i2,nspv))*etamag(2)
      gvec_brydn(i1,i2,1   )= (chrg_o  -rho_i(i1,i2,1   ))*eta
      gvec_brydn(i1,i2,nspv)= (pol_o(1)-rho_i(i1,i2,nspv))*etamag(1)
    end do
    end do
  case(4)
    do i2=1,nmax2
    !$omp do
    do i1=1,nmax1
      chrg_i  = rho_i(i1,i2,1)+rho_i(i1,i2,nspv-2)
      pol_i(1)= rho_i(i1,i2,1)-rho_i(i1,i2,nspv-2)
      pol_i(2)= 2.0d0*rho_i(i1,i2,nspv-1)
      pol_i(3)= 2.0d0*rho_i(i1,i2,nspv-0)
      chrg_o  = rho_o(i1,i2,1)+rho_o(i1,i2,nspv-2)
      pol_o(1)= rho_o(i1,i2,1)-rho_o(i1,i2,nspv-2)
      pol_o(2)= 2.0d0*rho_o(i1,i2,nspv-1)
      pol_o(3)= 2.0d0*rho_o(i1,i2,nspv-0)
      rho_i(i1,i2,nspv-3)= chrg_i
      rho_i(i1,i2,nspv-2)= pol_i(1) 
      rho_i(i1,i2,nspv-1)= pol_i(2) 
      rho_i(i1,i2,nspv-0)= pol_i(3) 
      fvec_brydn(i1,i2,nspv-3)= (chrg_o  -chrg_i  )
      fvec_brydn(i1,i2,nspv-2)= (pol_o(1)-pol_i(1))*etamag(2)
      fvec_brydn(i1,i2,nspv-1)= (pol_o(2)-pol_i(2))*etamag(2)
      fvec_brydn(i1,i2,nspv-0)= (pol_o(3)-pol_i(3))*etamag(2)
      gvec_brydn(i1,i2,nspv-3)= (chrg_o  -chrg_i  )*eta
      gvec_brydn(i1,i2,nspv-2)= (pol_o(1)-pol_i(1))*etamag(1)
      gvec_brydn(i1,i2,nspv-1)= (pol_o(2)-pol_i(2))*etamag(1)
      gvec_brydn(i1,i2,nspv-0)= (pol_o(3)-pol_i(3))*etamag(1)
    end do
    end do
  end select 

  deallocate(pol_i,pol_o,pol_n)

end subroutine scf_chargemixing_broyden_vecs


subroutine scf_chargemixing_broyden_rho( &
 ndim1,nmax1,nmax2,nspv, & ! <
 rho)                      ! X
implicit none
integer,intent(in)   :: ndim1,nmax1,nmax2,nspv
real*8, intent(inout):: rho(ndim1,nmax2,nspv)
integer:: i1,i2

  select case (nspv)
  case(1)
  case(2) 
    !$omp do
    do i2=1,nmax2
    do i1=1,nmax1
      rho(i1,i2,1   )= 0.5d0*(rho(i1,i2,1)+rho(i1,i2,nspv))
      rho(i1,i2,nspv)= rho(i1,i2,1)-rho(i1,i2,nspv)
    enddo
    enddo
  case(4) 
    !$omp do
    do i2=1,nmax2
    do i1=1,nmax1
      rho(i1,i2,nspv-3)= 0.5d0*(rho(i1,i2,1)+rho(i1,i2,nspv-2))
      rho(i1,i2,nspv-2)= rho(i1,i2,1)-rho(i1,i2,nspv-2)
      rho(i1,i2,nspv-1)= 0.5d0*rho(i1,i2,nspv-1)
      rho(i1,i2,nspv-0)= 0.5d0*rho(i1,i2,nspv-0) 
    enddo
    enddo
  end select

end subroutine scf_chargemixing_broyden_rho


subroutine scf_chargemixing_spinpol( & 
 chinp,ndisp,key_pp_paw,key_natpri_in,                     & ! <
 natom,num_spe,num_atcell,nspv,nradmx,npoint,              & ! <
 natpri,natpri_inf,indspe,ntyppp,nradct,wt,radial,dradial, & ! <
 rhotrur)                                                    ! <
use mod_mpi 
implicit none
character(len=*),intent(in):: chinp
integer,intent(in):: ndisp,key_pp_paw,key_natpri_in
integer,intent(in):: natom,num_spe,num_atcell,nspv,nradmx,npoint
integer,intent(in):: natpri(natom),natpri_inf(natom),indspe(natom)
integer,intent(in):: ntyppp(num_spe),nradct(num_spe)
real*8, intent(in):: wt(npoint)
real*8, intent(in):: radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8, intent(in):: rhotrur(nradmx,npoint,nspv,num_atcell)
integer  :: na,ipri,ispe,il,ir,ns,ns2
character:: chnspv*1
real*8 :: rfac,poltru
real*8,allocatable:: spinpol(:,:),spinpol1(:,:)

  if (nspv>1) then

  allocate(spinpol(max(1,nspv-1),natom),spinpol1(max(1,nspv-1),natom))

  spinpol1(:,:)= 0.0d0 
  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na)==key_natpri_in)) then
      ipri= natpri_inf(na)
      ispe= indspe(na)

      do ns=1,nspv
        poltru= 0.0d0
        do il=1,npoint
          do ir=2,nradct(ispe)-1
            rfac= wt(il)*dradial(ir,ispe)*radial(ir,ispe)**2
            poltru= poltru+rfac*rhotrur(ir,il,ns,ipri)
          end do
        end do
        ns2= max(1,ns-1)
        if (ns==2) then
          poltru= -poltru
        end if
        if (ns>2) then
          poltru= 2.0d0*poltru
        end if
        spinpol1(ns2,na)= spinpol1(ns2,na)+poltru
      end do

    end if ! ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na)==key_natpri_in))
  end do ! na

  call mpi_allreduce(spinpol1,spinpol,max(1,nspv-1)*natom,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  if (myrank_glbl==0) then
    write(chnspv,fmt='(i1)') nspv-1
    do na=1,natom
      write(ndisp,fmt='(3x,"atom #",i4,":  SpinPol",a,"=",'//chnspv//'(1x,e10.4))') &
       na, chinp, (spinpol(ns,na),ns=1,nspv-1)
    enddo
  endif
  
  deallocate(spinpol,spinpol1)

  endif ! (nspv>1)

end subroutine scf_chargemixing_spinpol


end module

