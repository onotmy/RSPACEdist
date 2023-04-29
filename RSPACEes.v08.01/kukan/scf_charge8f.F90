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
! **********  scf_charge8f.F90 06/20/2018-01  **********

module mod_scf_charge
implicit none
contains

! -------------------------------------------------------------------------------------------------

subroutine scf_charge( &
 key_natpri_in,key_pp_paw,key_sym_any,key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp,                            & ! <
 nd,nrc,nsym,numk,neigmx,nums,ncol,nspv,ncpx,ncpy,ncpz,num_spe,natom,num_atcell,num_ppcell,nprjmx,nprmx,lmx,lrhomx, & ! <
 nradmx,npoint,natpri,naps,                                                                                       & ! <
 natpri_inf,indspe,ntyppp,nprj,nlind,noind,nradct,dx,dy,dz,yylm,wt,awf,pwf,radial,dradial,                        & ! <
 fnele,svecre,sveccm,rspsep,cspsep,                                                                               & ! <
 atocc,rhosmt,rhotrur,rhosmtr,spinpol)                                                                              ! >
use mod_mpi
implicit none
integer,   intent(in) :: key_natpri_in,key_pp_paw,key_sym_any,key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp
integer,   intent(in) :: nd,nrc,nsym
integer,   intent(in) :: numk,neigmx,nums,ncol,nspv,ncpx,ncpz,ncpy,num_spe,natom,num_atcell,num_ppcell
integer,   intent(in) :: nprjmx,nprmx,lmx,lrhomx,nradmx,npoint
integer,   intent(in) :: natpri(natom),naps(natom),natpri_inf(natom),indspe(natom)
integer,   intent(in) :: ntyppp(num_spe),nprj(num_spe),nlind(nprjmx,num_spe),noind(nprjmx,num_spe),nradct(num_spe)
real*8,    intent(in) :: dx,dy,dz,yylm(npoint,lrhomx),wt(npoint)
real*8,    intent(in) :: awf(nradmx,nprmx*lmx,num_spe),pwf(nradmx,nprmx*lmx,num_spe)
real*8,    intent(in) :: radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8,    intent(in) :: fnele(neigmx,nums+1-ncol,numk)
real*8,    intent(in) :: svecre(ncpx*(1-nrc)+nrc,ncpy*(1-nrc)+nrc,ncpz*(1-nrc)+nrc,neigmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk*(1-nrc)+nrc)
complex*16,intent(in) :: sveccm(ncpx*nrc-nrc+1,ncpy*nrc-nrc+1,ncpz*nrc-nrc+1,neigmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1)
real*8,    intent(in) :: rspsep(nprjmx*(1-nrc)+nrc,num_ppcell*(1-nrc)+nrc,neigmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk*(1-nrc)+nrc)
complex*16,intent(in) :: cspsep(nprjmx*nrc-nrc+1,num_ppcell*nrc-nrc+1,neigmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1)
real*8,    intent(out):: atocc((nprjmx**2+nprjmx)/2,nspv,num_atcell)
real*8,    intent(out):: rhosmt(ncpx,ncpy,ncpz,nspv)
real*8,    intent(out):: rhotrur(nradmx*nd-nd+1,npoint*nd-nd+1,nspv*nd-nd+1,num_atcell*nd-nd+1)
real*8,    intent(out):: rhosmtr(nradmx*nd-nd+1,npoint*nd-nd+1,nspv*nd-nd+1,num_atcell*nd-nd+1)
real*8,    intent(out):: spinpol(max(1,nspv-1)*nd-nd+1,0:natom*min(1,nspv-1)*nd)
real*8, allocatable::rhoall(:,:,:,:)
real*8, allocatable::rhotrum(:,:,:,:),rhosmtm(:,:,:,:)

! ==========  charge distribution around nucleus (rhotrur: true charge, rhosmtr: smooth charge)  ==

  call scf_charge_atomicocc( &
   key_natpri_in,key_pp_paw,                                                  & ! <
   nrc,natom,num_spe,num_atcell,num_ppcell,nprjmx,nums,ncol,nspv,numk,neigmx, & ! <
   natpri,naps,natpri_inf,indspe,ntyppp,nprj,                                 & ! <
   fnele,rspsep,cspsep,                                                       & ! <
   atocc)                                                                       ! >
  if (myr_kpt==0) then
    allocate( rhoall((nprjmx**2+nprjmx)/2,nspv,num_atcell,1) )
  else
    allocate( rhoall(1,1,1,1) )
  end if
  call mpi_reduce(atocc(1,1,1),rhoall(1,1,1,1),((nprjmx**2+nprjmx)/2)*nspv*num_atcell,mpi_double_precision, &
   mpi_sum,0,mpicom_kpt,mpij)
  if (myr_kpt==0) atocc(:,:,:)= rhoall(:,:,:,1)
  deallocate( rhoall )

!$omp parallel shared(rhotrur,rhosmtr)
  if (myr_kpt==0) call scf_charge_onecenter( &
   key_natpri_in,key_pp_paw,                                                    & ! <
   natom,num_spe,num_atcell,nspv,nradmx,npoint,lmx,lrhomx,nprjmx,nprmx,         & ! <
   natpri,natpri_inf,indspe,ntyppp,nprj,nlind,noind,nradct,yylm,awf,pwf,radial, & ! <
   atocc,                                                                       & ! <
   rhotrur,rhosmtr)                                                               ! >
!$omp end parallel

! ==========  smooth charge density  ==

  if (nrc==0) then
!$omp parallel shared(rhosmt)
    call scf_charge_smt_r( &
     neigmx,nums,nspv,numk,ncpx,ncpy,ncpz, & ! <
     fnele,svecre,                       & ! <
     rhosmt)                               ! >
!$omp end parallel
  else
!$omp parallel shared(rhosmt)
    call scf_charge_smt_c( &
     neigmx,nums,ncol,nspv,numk,ncpx,ncpy,ncpz, & ! <
     fnele,sveccm,                            & ! <
     rhosmt)                                    ! >
!$omp end parallel
  end if
  if (myr_kpt==0) then
    allocate( rhoall(ncpx,ncpy,ncpz,nspv) )
  else
    allocate( rhoall(1,1,1,1) )
  end if
  call mpi_reduce(rhosmt(1,1,1,1),rhoall(1,1,1,1),ncpx*ncpy*ncpz*nspv,mpi_double_precision, &
   mpi_sum,0,mpicom_kpt,mpij)
  if (myr_kpt==0) rhosmt(:,:,:,:)= rhoall(:,:,:,:)
  deallocate( rhoall )

! symmetric operations
  if ((myr_kpt==0) .and. (nsym/=key_sym_any)) then
    call scf_charge_smt_symmetric( &
     key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp, & ! <
     nspv,nsym,ncpx,ncpy,ncpz,                        & ! <
     rhosmt)                                            ! X
  end if
  if ((myr_kpt==0) .and. (nsym/=key_sym_any)) then
  allocate(rhotrum(nradmx,lrhomx,nspv,num_atcell),rhosmtm(nradmx,lrhomx,nspv,num_atcell))
!$omp parallel default(shared)
  call scf_charge_smtr_symmetric(natom,num_spe,num_atcell,nspv,lrhomx,npoint,nsym,nradmx, & ! <
        key_natpri_in,key_pp_paw,                                                         & ! <
        key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp,                                  & ! <
        ntyppp,nradct,indspe,natpri,natpri_inf,rhosmtm,rhotrum,yylm,wt,                   & ! <
        rhosmtr,rhotrur)                                                                    ! X
!$omp end parallel
  deallocate(rhotrum,rhosmtm)
  end if

! ==========  spin polarization  ==

  spinpol(:,:)= 0.0d0

  if ((myr_kpt==0).and.(nspv>1)) then

    call scf_charge_spinpol_onecenter( &
     key_natpri_in,key_pp_paw,                                                                             & ! <
     natom,num_spe,num_atcell,nspv,nradmx,npoint,natpri,natpri_inf,indspe,ntyppp,nradct,wt,radial,dradial, & ! <
     rhotrur,rhosmtr,                                                                                      & ! <
     spinpol)                                                                                                ! X

    call scf_charge_spinpol_smt( &
     nspv,ncpx,ncpy,ncpz,dx,dy,dz, & ! <
     rhosmt,                       & ! <
     spinpol(1,0))                   ! X

    if (myrank_glbl==0) then
      allocate( rhoall(nspv-1,0:natom,1,1) )
    else
      allocate( rhoall(1,1,1,1) )
    end if
    call mpi_reduce(spinpol(1,0),rhoall(1,0,1,1),(nspv-1)*(natom+1),mpi_double_precision, &
     mpi_sum,0,mpicom_space,mpij)
    if (myrank_glbl==0) spinpol(:,:)= rhoall(:,:,1,1)
    deallocate( rhoall )

  end if ! ((myr_kpt==0).and.(nspv>1))

end subroutine scf_charge

! -- occupation of atomic wave functions ----------------------------------------------------------

subroutine scf_charge_atomicocc( &
 key_natpri_in,key_pp_paw,                                                  & ! <
 nrc,natom,num_spe,num_atcell,num_ppcell,nprjmx,nums,ncol,nspv,numk,neigmx, & ! <
 natpri,naps,natpri_inf,indspe,ntyppp,nprj,                                 & ! <
 fnele,rspsep,cspsep,                                                       & ! <
 atocc)                                                                       ! >
implicit none
integer,    intent(in)  :: key_natpri_in,key_pp_paw
integer,    intent(in)  :: nrc,natom,num_spe,num_atcell,num_ppcell,nprjmx,nums,ncol,nspv,numk,neigmx
integer,    intent(in)  :: natpri(natom),naps(natom),natpri_inf(natom),indspe(natom)
integer,    intent(in)  :: ntyppp(num_spe),nprj(num_spe)
real*8,     intent(in)  :: fnele(neigmx,nums+1-ncol,numk)
real*8,     intent(in)  :: rspsep(nprjmx*(1-nrc)+nrc,num_ppcell*(1-nrc)+nrc,neigmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk*(1-nrc)+nrc)
complex*16, intent(in)  :: cspsep(nprjmx*nrc-nrc+1,num_ppcell*nrc-nrc+1,neigmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1)
real*8,     intent(out) :: atocc((nprjmx**2+nprjmx)/2,nspv,num_atcell)
integer    :: na,iaps,ipri,ispe,nk,l,j,i,ij,ns,nsat,nsst, nsmax
real*8     :: tmpr
complex*16 :: tmpc

  atocc(:,:,:)= 0.0d0

  nsmax= max(nums,nspv)

  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na)==key_natpri_in)) then
      iaps= naps(na)
      ipri= natpri_inf(na)
      ispe= indspe(na)

      do nk=1,numk
        do l= 1,neigmx

          ij= 0
          do j=1,nprj(ispe)
            do i=1,j
              ij= ij + 1
              do ns= 1,nsmax
                nsst= min(ns,3-ncol) 
                nsat= min(ns,nspv)
                if (ns<3) then
                  if (nrc==0) then
                    tmpr= rspsep(i,iaps,l,ns,nk)*rspsep(j,iaps,l,ns,nk)
                  else
                    tmpr= dreal(dconjg(cspsep(i,iaps,l,ns,nk))*cspsep(j,iaps,l,ns,nk))
                  end if
                else
                  if (ns==3) then
                    tmpc= dconjg(cspsep(i,iaps,l,1,nk))*cspsep(j,iaps,l,ns-1,nk)  &
                        + dconjg(cspsep(j,iaps,l,1,nk))*cspsep(i,iaps,l,ns-1,nk)
                    tmpr= dreal(tmpc)/2.0d0
                  else
                    tmpr= dimag(tmpc)/2.0d0
                  end if
                end if
                if (i/=j) tmpr= 2.0d0*tmpr
                atocc(ij,nsat,ipri)= atocc(ij,nsat,ipri) + fnele(l,nsst,nk)*tmpr
              end do ! ns

            end do ! i
          end do ! j

        end do ! l
      end do ! nk

    end if
  end do ! na

end subroutine scf_charge_atomicocc

! -- onecenter charge -----------------------------------------------------------------------------

subroutine scf_charge_onecenter( &
 key_natpri_in,key_pp_paw,                                                            & ! <
 natom,num_spe,num_atcell,nspv,nradmx,npoint,lmx,lrhomx,nprjmx,nprmx,                 & ! <
 natpri,natpri_inf,indspe,ntyppp,nprj,nlind,noind,nradct,yylm,awf,pwf,radial,         & ! <
 atocc,                                                                               & ! <
 rhotrur,rhosmtr)                                                                       ! >
implicit none
integer,intent(in) :: key_natpri_in,key_pp_paw
integer,intent(in) :: natom,num_spe,num_atcell,nspv,nradmx,npoint,lmx,lrhomx,nprjmx,nprmx
integer,intent(in) :: natpri(natom),natpri_inf(natom),indspe(natom)
integer,intent(in) :: ntyppp(num_spe),nprj(num_spe),nlind(nprjmx,num_spe),noind(nprjmx,num_spe),nradct(num_spe)
real*8, intent(in) :: yylm(npoint,lrhomx)
real*8, intent(in) :: awf(nradmx,nprmx*lmx,num_spe),pwf(nradmx,nprmx*lmx,num_spe)
real*8, intent(in) :: radial(nradmx,num_spe)
real*8, intent(in) :: atocc((nprjmx**2+nprjmx)/2,nspv,num_atcell)
real*8, intent(out):: rhotrur(nradmx,npoint,nspv,num_atcell)
real*8, intent(out):: rhosmtr(nradmx,npoint,nspv,num_atcell)
integer :: na,ipri,ispe,ns,j,j1,j2,i,i1,i2,ij,il,ir
real*8  :: tmp1,tmp2,tmp3

!$omp do
  do i=1,nradmx*npoint*nspv*num_atcell
    rhotrur(i,1,1,1)= 0.0d0
    rhosmtr(i,1,1,1)= 0.0d0
  end do

  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na)==key_natpri_in)) then
      ipri= natpri_inf(na)
      ispe= indspe(na)

      ij= 0
      do j=1,nprj(ispe)
        j1=nlind(j,ispe)
        j2=noind(j,ispe)
        do i=1,j
          i1=nlind(i,ispe)
          i2=noind(i,ispe)
          ij= ij + 1
          do ns=1,nspv
!$omp do
            do il=1,npoint
              tmp3= yylm(il,i1)*yylm(il,j1)
              do ir=2,nradct(ispe)-1
                tmp1= awf(ir,i2,ispe)*awf(ir,j2,ispe)
                tmp2= pwf(ir,i2,ispe)*pwf(ir,j2,ispe)
                rhotrur(ir,il,ns,ipri)= rhotrur(ir,il,ns,ipri) + tmp1*tmp3*atocc(ij,ns,ipri)
                rhosmtr(ir,il,ns,ipri)= rhosmtr(ir,il,ns,ipri) + tmp2*tmp3*atocc(ij,ns,ipri)
              end do
            end do
          end do

        end do ! i
      end do ! j

      do ns=1,nspv
!$omp do
        do il=1,npoint
          do ir=2,nradct(ispe)
            tmp1=1.0d0/(radial(ir,ispe)*radial(ir,ispe))
            rhotrur(ir,il,ns,ipri)= rhotrur(ir,il,ns,ipri)*tmp1
            rhosmtr(ir,il,ns,ipri)= rhosmtr(ir,il,ns,ipri)*tmp1
          end do
        end do
      end do

    end if ! ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na)==key_natpri_in))
  end do ! na

end subroutine scf_charge_onecenter

! -- smooth charge on real-space grid -------------------------------------------------------------

subroutine scf_charge_smt_r( &
 neigmx,nums,nspv,numk,ncpx,ncpy,ncpz, & ! <
 fnele,svecre,                       & ! <
 rhosmt)                               ! >
implicit none
integer,intent(in)  :: neigmx,nums,nspv,numk,ncpx,ncpy,ncpz
real*8, intent(in)  :: fnele(neigmx,nums,numk)
real*8, intent(in)  :: svecre(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8, intent(out) :: rhosmt(ncpx,ncpy,ncpz,nspv)
real*8 wfre,tmp
integer nk,ns,nsr,l,ix

!$omp do
  do ix=1,ncpx*ncpy*ncpz*nspv
    rhosmt(ix,1,1,1)=0.0d0
  end do

  do ns=1,nums
    nsr= min(ns,nspv)
    do nk=1,numk
      do l=1,neigmx
        tmp=fnele(l,ns,nk)
        !$omp do
        do ix=1,ncpx*ncpy*ncpz
          wfre=svecre(ix,1,1,l,ns,nk)
          rhosmt(ix,1,1,nsr)=rhosmt(ix,1,1,nsr)+tmp*wfre*wfre
        end do
      end do
    end do
  end do

end subroutine scf_charge_smt_r


subroutine scf_charge_smt_c( &
 neigmx,nums,ncol,nspv,numk,ncpx,ncpy,ncpz, & ! <
 fnele,sveccm,                            & ! <
 rhosmt)                                    ! >
implicit none
integer,   intent(in)  :: neigmx,nums,ncol,nspv,numk,ncpx,ncpy,ncpz
real*8,    intent(in)  :: fnele(neigmx,nums+1-ncol,numk)
complex*16,intent(in)  :: sveccm(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8,    intent(out) :: rhosmt(ncpx,ncpy,ncpz,nspv)
real*8 wfre,wfim,tmp
integer nk,ns,nsl,nsr,l,ix

!$omp do
  do ix=1,ncpx*ncpy*ncpz*nspv
    rhosmt(ix,1,1,1)=0.0d0
  end do

  do nk=1,numk
    do ns=1,nums
      if (ncol==1) then
        nsl= ns
      else
        nsl= 1
      end if
      if (nspv>1) then
        nsr= ns
      else
        nsr= 1
      end if
      do l=1,neigmx
        tmp=fnele(l,nsl,nk)
!$omp do
        do ix=1,ncpx*ncpy*ncpz
          wfre=dreal(sveccm(ix,1,1,l,ns,nk))
          wfim=dimag(sveccm(ix,1,1,l,ns,nk))
          rhosmt(ix,1,1,nsr)=rhosmt(ix,1,1,nsr)+tmp*(wfre*wfre+wfim*wfim)
        end do
      end do
    end do

    if (nspv==4) then
      do l=1,neigmx
        tmp=fnele(l,1,nk)
!$omp do
        do ix=1,ncpx*ncpy*ncpz
          wfre= dreal(sveccm(ix,1,1,l,1,nk))*dreal(sveccm(ix,1,1,l,nums,nk))  &
               +dimag(sveccm(ix,1,1,l,1,nk))*dimag(sveccm(ix,1,1,l,nums,nk))
          wfim=-dimag(sveccm(ix,1,1,l,1,nk))*dreal(sveccm(ix,1,1,l,nums,nk))  &
               +dreal(sveccm(ix,1,1,l,1,nk))*dimag(sveccm(ix,1,1,l,nums,nk))
          rhosmt(ix,1,1,nspv-1)=rhosmt(ix,1,1,nspv-1)+tmp*wfre
          rhosmt(ix,1,1,nspv  )=rhosmt(ix,1,1,nspv  )+tmp*wfim
        end do
      end do
    end if

  end do ! nk

end subroutine scf_charge_smt_c


subroutine scf_charge_smt_symmetric( &
 key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp, & ! <
 nspv,nsym,ncpx,ncpy,ncpz,                        & ! <
 rhosmt)                                            ! X
implicit none
integer,intent(in)   :: key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp
integer,intent(in)   :: nspv,nsym,ncpx,ncpy,ncpz
real*8, intent(inout):: rhosmt(ncpx,ncpy,ncpz,nspv)
integer ns,ix,iy,iz
real*8,allocatable :: rhotmp(:,:,:)

  allocate(rhotmp(ncpx,ncpy,ncpz))

!$omp parallel default(shared) private(ns,ix,iy,iz)
  if ((nsym .eq. key_sym_bcc) .or. (nsym .eq. key_sym_fcc)) then
    do ns=1,nspv
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        rhotmp(ix,iy,iz)=0.125d0*(rhosmt(ix,iy,iz,ns)+rhosmt(ncpx-ix+1,iy,iz,ns) &
           +rhosmt(ix,ncpy-iy+1,iz,ns)+rhosmt(ix,iy,ncpz-iz+1,ns) &
           +rhosmt(ncpx-ix+1,ncpy-iy+1,iz,ns)+rhosmt(ncpx-ix+1,iy,ncpz-iz+1,ns) &
           +rhosmt(ix,ncpy-iy+1,ncpz-iz+1,ns)+rhosmt(ncpx-ix+1,ncpy-iy+1,ncpz-iz+1,ns))
      end do
      end do
      end do
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        rhosmt(ix,iy,iz,ns)=1.0d0/6.0d0*(rhotmp(ix,iy,iz)+rhotmp(iy,iz,ix)+rhotmp(iz,ix,iy) &
                                        +rhotmp(iz,iy,ix)+rhotmp(iy,ix,iz)+rhotmp(ix,iz,iy))
      end do
      end do
      end do
    end do
  end if
  if (nsym .eq. key_sym_dia) then
    do ns=1,nspv
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        rhotmp(ix,iy,iz)=0.25d0*(rhosmt(ix,iy,iz,ns)+rhosmt(ncpx-ix+1,ncpy-iy+1,iz,ns) &
           +rhosmt(ix,ncpy-iy+1,ncpz-iz+1,ns)+rhosmt(ncpx-ix+1,iy,ncpz-iz+1,ns))
      end do
      end do
      end do
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        rhosmt(ix,iy,iz,ns)=1.0d0/3.0d0*(rhotmp(ix,iy,iz)+rhotmp(iy,iz,ix)+rhotmp(iz,ix,iy))
      end do
      end do
      end do
    end do
  end if
  if (nsym .eq. key_sym_hcp) then
    do ns=1,nspv
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        rhotmp(ix,iy,iz)=0.25d0*(rhosmt(ix,iy,iz,ns)+rhosmt(ncpx-ix+1,iy,iz,ns) &
           +rhosmt(ix,iy,ncpz-iz+1,ns)+rhosmt(ncpx-ix+1,iy,ncpz-iz+1,ns))
      end do
      end do
      end do
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        rhosmt(ix,iy,iz,ns)=rhotmp(ix,iy,iz)
      end do
      end do
      end do
    end do
  end if
!$omp end parallel

  deallocate( rhotmp )

end subroutine scf_charge_smt_symmetric

! -- spin polarization ----------------------------------------------------------------------------

subroutine scf_charge_spinpol_onecenter( &
 key_natpri_in,key_pp_paw,                                 & ! <
 natom,num_spe,num_atcell,nspv,nradmx,npoint,              & ! <
 natpri,natpri_inf,indspe,ntyppp,nradct,wt,radial,dradial, & ! <
 rhotrur,rhosmtr,                                          & ! <
 spinpol)                                                    ! X
implicit none
integer,intent(in)   :: key_natpri_in,key_pp_paw
integer,intent(in)   :: natom,num_spe,num_atcell,nspv,nradmx,npoint
integer,intent(in)   :: natpri(natom),natpri_inf(natom),indspe(natom)
integer,intent(in)   :: ntyppp(num_spe),nradct(num_spe)
real*8, intent(in)   :: wt(npoint)
real*8, intent(in)   :: radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8, intent(in)   :: rhotrur(nradmx,npoint,nspv,num_atcell)
real*8, intent(in)   :: rhosmtr(nradmx,npoint,nspv,num_atcell)
real*8, intent(inout):: spinpol(max(1,nspv-1),0:natom)
integer:: na,ipri,ispe,il,ir,ns,ns2
real*8 :: rfac,poltru,polsmt

  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na)==key_natpri_in)) then
      ipri= natpri_inf(na)
      ispe= indspe(na)

      do ns=1,nspv
        poltru= 0.0d0
        polsmt= 0.0d0
        do il=1,npoint
          do ir=2,nradct(ispe)
            rfac= wt(il)*dradial(ir,ispe)*radial(ir,ispe)**2
            poltru= poltru+rfac*rhotrur(ir,il,ns,ipri)
            polsmt= polsmt+rfac*rhosmtr(ir,il,ns,ipri)
          end do
        end do
        ns2= max(1,ns-1)
        if (ns==2) then
          poltru= -poltru
          polsmt= -polsmt
        end if
        if (ns>2) then
          poltru= 2.0d0*poltru
          polsmt= 2.0d0*polsmt
        end if
        spinpol(ns2,na)= spinpol(ns2,na)+poltru
        spinpol(ns2,0 )= spinpol(ns2,0 )-polsmt
      end do

      do ns=1,nspv-1
        spinpol(ns,0)= spinpol(ns,0) +spinpol(ns,na)
      end do

    end if ! ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na)==key_natpri_in))
  end do ! na

end subroutine scf_charge_spinpol_onecenter


subroutine scf_charge_spinpol_smt( &
 nspv,ncpx,ncpy,ncpz,dx,dy,dz, & ! <
 rhosmt,                       & ! <
 spinpol)                        ! X
implicit none
integer,intent(in)   :: nspv,ncpx,ncpy,ncpz
real*8, intent(in)   :: dx,dy,dz
real*8, intent(in)   :: rhosmt(ncpx,ncpy,ncpz,nspv)
real*8, intent(inout):: spinpol(max(1,nspv-1))
integer :: ix,ns,ns2
real*8  :: spintmp

  if (nspv>1) then
    do ns= 1,nspv
      spintmp= 0.0d0
      do ix= 1,ncpx*ncpy*ncpz
        spintmp= spintmp + rhosmt(ix,1,1,ns)
      end do
      if (ns==2) spintmp= -spintmp
      if (ns> 2) spintmp= 2.0d0*spintmp
      ns2= max(1,ns-1)
      spinpol(ns2)= spinpol(ns2) + spintmp*dx*dy*dz
    end do
  end if

end subroutine scf_charge_spinpol_smt

! -------------------------------------------------------------------------------------------------


subroutine scf_charge_smtr_symmetric(natom,num_spe,num_atcell,nspv,lrhomx,npoint,nsym,nradmx, & ! <
            key_natpri_in,key_pp_paw,                                                         & ! <
            key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp,                                  & ! <
            ntyppp,nradct,indspe,natpri,natpri_inf,rhosmtm,rhotrum,yylm,wt,                   & ! <
            rhosmtr,rhotrur)                                                                    ! X
! this routine executes symmetric operation of the augmented charge.
use mod_mpi
implicit none
integer, intent(in)::natom,num_spe,num_atcell,nspv,lrhomx,npoint,nsym,nradmx
integer, intent(in)::key_natpri_in &
     ,key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp,key_pp_paw
integer, intent(in)::ntyppp(num_spe),nradct(num_spe)
integer, intent(in)::indspe(natom),natpri(natom),natpri_inf(natom)
real*8, intent(in)::yylm(npoint,lrhomx),wt(npoint)
real*8, intent(out)::rhotrum(nradmx,lrhomx,nspv,num_atcell)
real*8, intent(out)::rhosmtm(nradmx,lrhomx,nspv,num_atcell)
real*8, intent(inout)::rhotrur(nradmx,npoint,nspv,num_atcell)
real*8, intent(inout)::rhosmtr(nradmx,npoint,nspv,num_atcell)
integer na,ipri,ns,la,il,ir

!$omp single
 rhosmtm=0.0d0
 rhotrum=0.0d0
!$omp end single

  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
      ipri=natpri_inf(na)
      do ns=1,nspv
!cdir nodep
!$omp do
        do la=1,lrhomx
          do il=1,npoint
            do ir=2,nradct(indspe(na))-1
              rhosmtm(ir,la,ns,ipri)=rhosmtm(ir,la,ns,ipri)+rhosmtr(ir,il,ns,ipri)*yylm(il,la)*wt(il)
              rhotrum(ir,la,ns,ipri)=rhotrum(ir,la,ns,ipri)+rhotrur(ir,il,ns,ipri)*yylm(il,la)*wt(il)
            end do
          end do
        end do
      end do
    end if
  end do

!$omp single
! ==========  symmetric operations  ==========
  if ((nsym .eq. key_sym_bcc) .or. (nsym .eq. key_sym_fcc) .or. (nsym .eq. key_sym_dia)) then
    rhosmtm(:, 2,:,:)=0.0d0
    rhosmtm(:, 3,:,:)=0.0d0
    rhosmtm(:, 4,:,:)=0.0d0
    rhosmtm(:, 5,:,:)=0.0d0
    rhosmtm(:, 6,:,:)=0.0d0
    rhosmtm(:, 7,:,:)=0.0d0
    rhosmtm(:, 8,:,:)=0.0d0
    rhosmtm(:, 9,:,:)=0.0d0
    rhosmtm(:,10,:,:)=0.0d0
    rhosmtm(:,11,:,:)=0.0d0
    rhosmtm(:,12,:,:)=0.0d0
    rhosmtm(:,13,:,:)=0.0d0
    rhosmtm(:,14,:,:)=0.0d0
    rhosmtm(:,16,:,:)=0.0d0
    rhosmtm(:,18,:,:)=0.0d0
    rhosmtm(:,19,:,:)=0.0d0
    rhosmtm(:,20,:,:)=0.0d0
    rhosmtm(:,22,:,:)=0.0d0
    rhosmtm(:,23,:,:)=0.0d0
    rhosmtm(:,24,:,:)=0.0d0
    rhosmtm(:,25,:,:)=0.0d0
    rhotrum(:, 2,:,:)=0.0d0
    rhotrum(:, 3,:,:)=0.0d0
    rhotrum(:, 4,:,:)=0.0d0
    rhotrum(:, 5,:,:)=0.0d0
    rhotrum(:, 6,:,:)=0.0d0
    rhotrum(:, 7,:,:)=0.0d0
    rhotrum(:, 8,:,:)=0.0d0
    rhotrum(:, 9,:,:)=0.0d0
    rhotrum(:,10,:,:)=0.0d0
    rhotrum(:,11,:,:)=0.0d0
    rhotrum(:,12,:,:)=0.0d0
    rhotrum(:,13,:,:)=0.0d0
    rhotrum(:,14,:,:)=0.0d0
    rhotrum(:,16,:,:)=0.0d0
    rhotrum(:,18,:,:)=0.0d0
    rhotrum(:,19,:,:)=0.0d0
    rhotrum(:,20,:,:)=0.0d0
    rhotrum(:,22,:,:)=0.0d0
    rhotrum(:,23,:,:)=0.0d0
    rhotrum(:,24,:,:)=0.0d0
    rhotrum(:,25,:,:)=0.0d0
  end if
  if ((nsym .eq. key_sym_bcc) .or. (nsym .eq. key_sym_fcc)) then
    rhosmtm(:,15,:,:)=0.0d0
    rhotrum(:,15,:,:)=0.0d0
  end if
  if (nsym .eq. key_sym_hcp) then
    rhosmtm(:, 2,:,:)=0.0d0
    rhosmtm(:, 4,:,:)=0.0d0
    rhosmtm(:, 6,:,:)=0.0d0
    rhosmtm(:, 8,:,:)=0.0d0
    rhosmtm(:, 9,:,:)=0.0d0
    rhosmtm(:,10,:,:)=0.0d0
    rhosmtm(:,11,:,:)=0.0d0
    rhosmtm(:,12,:,:)=0.0d0
    rhosmtm(:,13,:,:)=0.0d0
    rhosmtm(:,15,:,:)=0.0d0
    rhosmtm(:,18,:,:)=0.0d0
    rhosmtm(:,20,:,:)=0.0d0
    rhosmtm(:,22,:,:)=0.0d0
    rhosmtm(:,23,:,:)=0.0d0
    rhosmtm(:,24,:,:)=0.0d0
    rhosmtm(:,25,:,:)=0.0d0
    rhotrum(:, 2,:,:)=0.0d0
    rhotrum(:, 4,:,:)=0.0d0
    rhotrum(:, 6,:,:)=0.0d0
    rhotrum(:, 8,:,:)=0.0d0
    rhotrum(:, 9,:,:)=0.0d0
    rhotrum(:,10,:,:)=0.0d0
    rhotrum(:,11,:,:)=0.0d0
    rhotrum(:,12,:,:)=0.0d0
    rhotrum(:,13,:,:)=0.0d0
    rhotrum(:,15,:,:)=0.0d0
    rhotrum(:,18,:,:)=0.0d0
    rhotrum(:,20,:,:)=0.0d0
    rhotrum(:,22,:,:)=0.0d0
    rhotrum(:,23,:,:)=0.0d0
    rhotrum(:,24,:,:)=0.0d0
    rhotrum(:,25,:,:)=0.0d0
  end if
!     ============================================
!$omp end single

!$omp single
  rhosmtr=0.0d0
  rhotrur=0.0d0
!$omp end single

  do na=1,natom
    if (natpri(na) .eq. key_natpri_in) then
      ipri=natpri_inf(na)
      do ns=1,nspv
        do la=1,lrhomx
!cdir nodep
!$omp do
          do il=1,npoint
            do ir=2,nradct(indspe(na))-1
              rhosmtr(ir,il,ns,ipri)=rhosmtr(ir,il,ns,ipri)+rhosmtm(ir,la,ns,ipri)*yylm(il,la)
              rhotrur(ir,il,ns,ipri)=rhotrur(ir,il,ns,ipri)+rhotrum(ir,la,ns,ipri)*yylm(il,la)
            end do
          end do
        end do
      end do
    end if
  end do

  return
  end subroutine


end module

