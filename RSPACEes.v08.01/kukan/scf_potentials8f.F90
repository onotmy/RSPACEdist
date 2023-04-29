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
! **********  scf_potential8f.f90 01/06/2020-01  **********

module mod_scf_potentials
contains

subroutine scf_potentials( &
 key_natpri_in,key_natpri_inps,key_pp_paw,                                             & ! <
 nmesh,ncpx,ncpy,ncpz,natom,num_atcell,num_spe,nradmx,npoint,lrhomx,lmx,               & ! <
 num_list_d,num_ppcell_d,nspv,                                                         & ! <
 npopshow,nevhist,ndisp,nperi,nf,nfh,ncgres,ncgmin,ncgmax,                             & ! <
 indspe,natpri,natprid,natpri_inf,napsd,natinfd,ndatx,ndaty,ndatz,                     & ! <
 ntyppp,nradct,lpmx,nwexp,lstvecd2,lstdx,lstdy,lstdz,cexco,                            & ! <
 epsvh,xmax,ymax,zmax,dx,dy,dz,tnumele,biasx,biasy,biasz,                              & ! <
 yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2,d2ylm_dtheta_dphi,               & ! <
 wt,point,radial,dradial,cp,rfac,atx,aty,atz,pwei,                                     & ! <
 rhocore,rhopccr,rhotrur,rhosmtr,rhopcc_dense,rhosmt,vloc_coarse,                      & ! <
 vboundx,vboundy,vboundz,                                                              & ! <
 atmpole,rhotrucorer,rhosmt_pcc_dense,rhosmt_pccr,rhoaugr,rhoaug3d,rho_aug_dense,      & ! >
 veff,vh_coarse,vhtrur,vhsmtr,vhaugr,vx,vx_dense,vxctru,vxcsmt,ex_dense,exctru,excsmt, & ! >
 vh_dense)                                                                               ! X
use mod_scf_augcharge,       only:scf_augcharge
use mod_scf_radialmoment,    only:scf_radialmoment
use mod_scf_rhoaugdense,     only:scf_rhoaugdense
use mod_scf_fuzzycellmoment, only:scf_fuzzycellmoment
use mod_scf_vhbound,         only:scf_vhbound
use mod_scf_hartree,         only:scf_hartree
use mod_scf_onecentervh,     only:scf_onecentervh
use mod_scf_rhoxcalc,        only:scf_rhoxcalc_smtpcc, scf_rhoxcalc_onecenter
use mod_scf_vxcalc,          only:scf_vxcalc
use mod_tools,               only:tools_correctpotential
use mod_trans,               only:trans_c2d_smtcharge, trans_d2c_vh, trans_d2c_vxc
use mod_mpi, only:myrx,myry,myrz
implicit none
integer,  intent(in)   :: key_natpri_in,key_natpri_inps,key_pp_paw
integer,  intent(in)   :: nmesh,ncpx,ncpy,ncpz
integer,  intent(in)   :: natom,num_atcell,num_spe,nradmx,npoint,lrhomx,lmx
integer,  intent(in)   :: num_list_d,num_ppcell_d
integer,  intent(in)   :: nspv
integer,  intent(in)   :: npopshow,nevhist,ndisp
integer,  intent(in)   :: nperi,nf,nfh,ncgres,ncgmin,ncgmax
integer,  intent(in)   :: indspe(natom),natpri(natom),natprid(natom),natpri_inf(natom),napsd(natom),natinfd(natom)
integer,  intent(in)   :: ndatx(natom),ndaty(natom),ndatz(natom)
integer,  intent(in)   :: ntyppp(num_spe),nradct(num_spe),lpmx(num_spe),nwexp(num_spe)
integer,  intent(in)   :: lstvecd2(num_list_d,num_ppcell_d)
integer,  intent(in)   :: lstdx(num_list_d,num_ppcell_d),lstdy(num_list_d,num_ppcell_d),lstdz(num_list_d,num_ppcell_d)
character,intent(in)   :: cexco*7
real*8,   intent(in)   :: epsvh
real*8,   intent(in)   :: xmax,ymax,zmax,dx,dy,dz
real*8,   intent(in)   :: tnumele
real*8,   intent(in)   :: biasx,biasy,biasz
real*8,   intent(in)   :: yylm(npoint,lrhomx),dylm_dtheta(npoint,lrhomx)
real*8,   intent(in)   :: d2ylm_dtheta2(npoint,lrhomx),dylm_dphi(npoint,lrhomx)
real*8,   intent(in)   :: d2ylm_dphi2(npoint,lrhomx),d2ylm_dtheta_dphi(npoint,lrhomx)
real*8,   intent(in)   :: wt(npoint),point(npoint,3)
real*8,   intent(in)   :: radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8,   intent(in)   :: cp(8,num_spe),rfac(0:2*(lmx-1),num_spe)
real*8,   intent(in)   :: atx(natom),aty(natom),atz(natom)
real*8,   intent(in)   :: pwei(ncpx,ncpy,ncpz,natom)
real*8,   intent(in)   :: rhocore(nradmx,num_spe)
real*8,   intent(in)   :: rhopccr(nradmx,npoint,num_atcell)
real*8,   intent(in)   :: rhotrur(nradmx,npoint,nspv,num_atcell)
real*8,   intent(in)   :: rhosmtr(nradmx,npoint,nspv,num_atcell)
real*8,   intent(in)   :: rhopcc_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
real*8,   intent(in)   :: rhosmt(ncpx,ncpy,ncpz,nspv)
real*8,   intent(in)   :: vloc_coarse(ncpx,ncpy,ncpz)
real*8,   intent(in)   :: vboundx(-(nfh-1):nfh,ncpy*nmesh,ncpz*nmesh,9*(3-nperi)/3+1,(natom-1)*(3-nperi)/3+1)
real*8,   intent(in)   :: vboundy(ncpx*nmesh,-(nfh-1):nfh,ncpz*nmesh,9*(3-nperi)/2+1,(natom-1)*(3-nperi)/2+1)
real*8,   intent(in)   :: vboundz(ncpx*nmesh,ncpy*nmesh,-(nfh-1):nfh,9*(1-nperi/3)+1,(natom-1)*(1-nperi/3)+1)
real*8,   intent(out)  :: atmpole(10,natom)
real*8,   intent(out)  :: rhotrucorer(nradmx,npoint,nspv,num_atcell)
real*8,   intent(out)  :: rhosmt_pcc_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh,nspv)
real*8,   intent(out)  :: rhosmt_pccr(nradmx,npoint,nspv,num_atcell)
real*8,   intent(out)  :: rhoaugr(nradmx,npoint,num_atcell)
real*8,   intent(out)  :: rhoaug3d(num_list_d,num_ppcell_d)
real*8,   intent(out)  :: rho_aug_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
real*8,   intent(out)  :: veff(ncpx,ncpy,ncpz,nspv)
real*8,   intent(out)  :: vh_coarse(ncpx,ncpy,ncpz)
real*8,   intent(out)  :: vhtrur(nradmx,npoint,num_atcell)
real*8,   intent(out)  :: vhsmtr(nradmx,npoint,num_atcell)
real*8,   intent(out)  :: vhaugr(nradmx,npoint,num_atcell)
real*8,   intent(out)  :: vx(ncpx,ncpy,ncpz,nspv)
real*8,   intent(out)  :: vx_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh,nspv)
real*8,   intent(out)  :: vxctru(nradmx,npoint,nspv,num_atcell)
real*8,   intent(out)  :: vxcsmt(nradmx,npoint,nspv,num_atcell)
real*8,   intent(out)  :: ex_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
real*8,   intent(out)  :: exctru(nradmx,npoint,num_atcell)
real*8,   intent(out)  :: excsmt(nradmx,npoint,num_atcell)
real*8,   intent(inout):: vh_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
integer ::ncpx_d,ncpy_d,ncpz_d
integer ::ix,iy,iz,ns
real*8  ::ddx,ddy,ddz
real*8  ::vdia
real*8,allocatable:: dummy1(:,:),dummy2(:,:),dummy3(:,:)
real*8,allocatable:: rho_dense(:,:,:,:),rhotrum(:,:,:),rhosmtm(:,:,:),rhoaugm(:,:,:)
real*8,allocatable:: boundx(:,:,:),boundy(:,:,:),boundz(:,:,:)


! ************************************************
  ncpx_d=nmesh*ncpx
  ncpy_d=nmesh*ncpy
  ncpz_d=nmesh*ncpz
  ddx=dx/dble(nmesh)
  ddy=dy/dble(nmesh)
  ddz=dz/dble(nmesh)

  allocate( &
   dummy1(1,1),dummy2(1,1),dummy3(1,1), &
   rho_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh,nspv), &
   rhotrum(nradmx,lrhomx,num_atcell), rhosmtm(nradmx,lrhomx,num_atcell), rhoaugm(nradmx,lrhomx,num_atcell) )
  if (nperi/=3) then
    allocate( &
     boundx(-(nfh-1):nfh,ncpy_d,ncpz_d), boundy(ncpx_d,-(nfh-1):nfh,ncpz_d), boundz(ncpx_d,ncpy_d,-(nfh-1):nfh) )
  else
    allocate( boundx(1,1,1), boundy(1,1,1), boundz(1,1,1) )
  endif

! **********  compute augmented charge density **********
  call scf_augcharge( &
   0,key_natpri_in,key_natpri_inps,key_pp_paw,                       & ! <
   nspv,nradmx,npoint,lrhomx,lmx,natom,num_atcell,num_spe,           & ! <
   num_list_d,num_ppcell_d,                                          & ! <
   indspe,natpri,natprid,natpri_inf,napsd,natinfd,ndatx,ndaty,ndatz, & ! <
   ntyppp,nradct,lpmx,nwexp,lstdx,lstdy,lstdz,                       & ! <
   ddx,ddy,ddz,yylm,wt,radial,dradial,rfac,                          & ! <
   atx,aty,atz,rhotrur,rhosmtr,                                      & ! <
   rhoaugr,rhoaug3d, dummy1,dummy2,dummy3)                             ! >
! **********************************************************

! **********  compute charge density on dense grid  **********
  call trans_c2d_smtcharge( &
   nmesh,nspv,nperi,ndisp,nf,ncpx,ncpy,ncpz, & ! <
   rhosmt,                                   & ! <
   rho_dense)                                  ! >
! ************************************************************

! **********  compute Hartree potential  **********
! ==========  compute radial part of charge distribution moment  ==========
!$omp parallel default(shared)
  call scf_radialmoment( &
   key_natpri_in,key_pp_paw,                                 & ! <
   natom,num_spe,num_atcell,nradmx,npoint,lrhomx,nspv,nperi, & ! <
   indspe,natpri,natpri_inf,ntyppp,nradct,lpmx,              & ! <
   xmax,ymax,zmax,tnumele,                                   & ! <
   yylm,wt,                                                  & ! <
   rhotrur,rhosmtr,rhoaugr,                                  & ! <
   rhotrum,rhosmtm,rhoaugm)                                    ! >
!$omp end parallel
! =========================================================================

! ==========  compute charge density for Poisson equation  ==========
!$omp parallel default(shared)
  call scf_rhoaugdense( &
   key_natpri_inps,key_pp_paw,                                      & ! <
   natom,num_spe,nspv,ncpx_d,ncpy_d,ncpz_d,num_list_d,num_ppcell_d, & ! <
   indspe,natprid,napsd,natinfd,ntyppp,lstvecd2,                    & ! <
   rho_dense,rhoaug3d,                                              & ! <
   rho_aug_dense)                                                     ! >
!$omp end parallel
! ===================================================================

! ==========  compute boundary condition of Poisson eq.  ==========
  if ((nperi .lt. 3) .or. (npopshow .eq. 1)) then
!$omp parallel default(shared)
    call scf_fuzzycellmoment( &
     key_natpri_inps,key_pp_paw,nperi,nspv,natom,num_spe,num_ppcell_d,        & ! <
     num_list_d,nmesh,ncpx,ncpy,ncpz,                                         & ! <
     ntyppp,indspe,natprid,natinfd,napsd,ndatx,ndaty,ndatz,lstdx,lstdy,lstdz, & ! <
     dx,dy,dz,xmax,ymax,zmax,                                                 & ! <
     atx,aty,atz,rhosmt,rhoaug3d,pwei,                                        & ! <
     atmpole)                                                                   ! >
!$omp end parallel
  end if
  if (nperi .lt. 3) then
!$omp parallel default(shared)
    call scf_vhbound( &
     natom,num_spe,nperi,indspe,nfh,ncpx_d,ncpy_d,ncpz_d, & ! <
     atmpole,vboundx,vboundy,vboundz,                     & ! <
     boundx,boundy,boundz)                                  ! >
!$omp end parallel
  end if
! =================================================================

  if (nevhist>=0) then
! ==========  solve Poisson equation  ==========
    call scf_hartree( &
     nperi/3,ncpx_d,ncpy_d,ncpz_d,nfh,ndisp,nperi,ncgres,ncgmin,ncgmax,xmax,ymax,zmax,epsvh, & ! <
     rho_aug_dense,boundx,boundy,boundz,                                                     & ! <
     vh_dense)                                                                                 ! X
! ==============================================

! ==========  compensate the divergence of the integration of Coulomb potential  ==========
    if (nperi .eq. 3) then
      call tools_correctpotential(ncpx_d,ncpy_d,ncpz_d,xmax,ymax,zmax,ddx,ddy,ddz,vh_dense)
    end if
! =========================================================================================
  endif

! ==========  compute Hartree potential on coarse grid  ==========
    call trans_d2c_vh(natom,num_spe,nmesh,nperi,ndisp,nf, & ! <
                      ncpx,ncpy,ncpz,                     & ! <
                      xmax,ymax,zmax,                     & ! <
                      indspe,                             & ! <
                      atmpole,cp,atx,aty,atz,             & ! <
                      vh_coarse,                          & ! >
                      vh_dense)                             ! <
! ****************************************************************

! **********  compute one-center Hartree potential  **************
  call scf_onecentervh( &
   key_natpri_in,key_pp_paw,                      & ! <
   natom,num_spe,num_atcell,lrhomx,npoint,nradmx, & ! <
   indspe,natpri,natpri_inf,ntyppp,lpmx,nradct,   & ! <
   yylm,radial,dradial,                           & ! <
   rhotrum,rhosmtm,rhoaugm,                       & ! <
   vhtrur,vhsmtr,vhaugr)                            ! >
! ****************************************************************

! **********  compute Ex. Cor. potential  **********
! ==========  compute charge density for Ex. Cor. potential  ====================
  call scf_rhoxcalc_smtpcc( &
   nspv,ncpx_d,ncpy_d,ncpz_d, & ! <
   rhopcc_dense,rho_dense,    & ! <
   rhosmt_pcc_dense)            ! >

! ==========  compute charge density for one-center Ex. Cor. potential  =========
  call scf_rhoxcalc_onecenter( &
   key_natpri_in,key_pp_paw,                    & ! <
   nspv,natom,num_spe,num_atcell,npoint,nradmx, & ! <
   indspe,natpri,natpri_inf,ntyppp,nradct,      & ! <
   rhocore,rhopccr,rhotrur,rhosmtr,             & ! <
   rhotrucorer,rhosmt_pccr)                       ! >
! ===============================================================================

! ==========  compute Ex. Cor. potential  =======================================
  call scf_vxcalc(  &
   natom,nperi,nfh,num_spe,nspv,npoint,nradmx,num_atcell,                  & ! <
   ncpx_d,ncpy_d,ncpz_d,                                                   & ! <
   key_natpri_in,                                                          & ! <
   ddx,ddy,ddz,                                                            & ! <
   cexco,                                                                  & ! <
   indspe,natpri,natpri_inf,nradct,lrhomx,                                 & ! <
   yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2,d2ylm_dtheta_dphi, & ! <
   point,wt,radial,dradial,                                                & ! <
   rhosmt_pcc_dense,rhotrucorer,rhosmt_pccr,                               & ! <
   vx_dense,vxctru,vxcsmt,ex_dense,exctru,excsmt)                            ! >
! ==================================================

! ==========  compute Ex. Cor. potential on coarse grid  ==========
  call trans_d2c_vxc( &
   nmesh,nspv,nperi,ndisp,nf,ncpx,ncpy,ncpz, & ! <
   vx,                                       & ! <
   vx_dense)                                   ! >
! =================================================================

! **********  sum up contributions of effective potential **********
!$omp parallel default(shared) private(ix,iy,iz,vdia)
  if (nspv<3) then
    do ns=1,nspv
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
      vdia= vloc_coarse(ix,iy,iz) +vh_coarse(ix,iy,iz) &
       + biasx * ( ( dble(myrx*ncpx+ix) -0.5d0 ) * dx - xmax ) &
       + biasy * ( ( dble(myry*ncpy+iy) -0.5d0 ) * dy - ymax ) &
       + biasz * ( ( dble(myrz*ncpz+iz) -0.5d0 ) * dz - zmax )
      veff(ix,iy,iz,ns)= vx(ix,iy,iz,ns)+vdia
    enddo
    enddo
    enddo
    enddo
  else
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
      vdia= vloc_coarse(ix,iy,iz) +vh_coarse(ix,iy,iz) &
       + biasx * ( ( dble(myrx*ncpx+ix) -0.5d0 ) * dx - xmax ) &
       + biasy * ( ( dble(myry*ncpy+iy) -0.5d0 ) * dy - ymax ) &
       + biasz * ( ( dble(myrz*ncpz+iz) -0.5d0 ) * dz - zmax )
      do ns= 1,nspv
        veff(ix,iy,iz,ns)= vx(ix,iy,iz,ns)
        if (ns<3) veff(ix,iy,iz,ns)= veff(ix,iy,iz,ns) +vdia
      enddo
    enddo
    enddo
    enddo
  end if
!$omp end parallel
! ******************************************************************

  deallocate( dummy1,dummy2,dummy3,rho_dense,rhotrum,rhosmtm,rhoaugm, boundx,boundy,boundz )

end subroutine scf_potentials

end module mod_scf_potentials
