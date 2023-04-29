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
! **********  pseudocalc8f.F90 04/27/2023-01  **********

module mod_pseudocalc

implicit none
contains


subroutine pseudocalc(natom,num_spe,nperi,nfdg,npmesh,nmesh,nf,nfh,ndisp,lrhomx,npoint,lsphel,nfiltyp,nqmx,nprjmx,nradmx,nprmx,lmx, & ! <
                      jelcalc,nint1dmax,                                                                                     & ! <
                      ncpx,ncpy,ncpz,ncpx_d,ncpy_d,ncpz_d,npxmax,npymax,npzmax,nxmax,nymax,nzmax,                            & ! <
                      new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz,                                                       & ! <
                      num_atcell,num_ppcell,num_ppcell_d,num_list,num_list_d,                                                & ! <
                      key_pp_paw,key_natpri_in,key_natpri_inps,key_natpri_out,key_jel_calc,                                  & ! <
                      eps,veta,psctoff,psftrad,filpp,rctpcc,chrjel,strjel,endjel,                                            & ! <
                      xmax,ymax,zmax,biasx,biasy,biasz,                                                                      & ! <
                      indspe,nlind,noind,nqct,nqctpcc,ntyppp,nprj,lpmx,nradct,                                               & ! <
                      coef,cp,radial,dradial,potc,awf,pwf,                                                                   & ! <
                      yylm,point,wt,                                                                                         & ! <
                      atx,aty,atz,                                                                                           & ! <
                      natpri,natprid,natpri_inf,natinf,natinfd,natinfd_vloc,naps,napsd,lstvec2,lstvecd2,                     & ! >
                      natx,naty,natz,ndatx,ndaty,ndatz,lstx,lsty,lstz,lstdx,lstdy,lstdz,                                     & ! >
                      qijl,vloc_coarse,vloc_dense,vcorer_all,vcorer,vnlocp,rhopcc_dense,rhopccr,                             & ! >
                      vloc_scw,vloc_hdp,dvlocdx_scw,dvlocdy_scw,dvlocdz_scw,dvlocdx_hdp,dvlocdy_hdp,dvlocdz_hdp,             & ! >
                      vboundx,vboundy,vboundz,vjell)                                                                           ! >
use mod_stopp
use mod_overlap_interpolation, only:overlap_interpolation
use mod_interpolation, only:interpolation
implicit none
integer,intent(in)::natom,num_spe,nperi,nfdg,npmesh,nmesh,nf,nfh,ndisp,lrhomx,npoint,lsphel,nfiltyp,nqmx,nprjmx,nradmx,nprmx,lmx
integer,intent(in)::jelcalc,nint1dmax
integer,intent(in)::ncpx,ncpy,ncpz,ncpx_d,ncpy_d,ncpz_d,npxmax,npymax,npzmax,nxmax,nymax,nzmax
integer,intent(in)::new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz
integer,intent(in)::num_atcell,num_ppcell,num_ppcell_d,num_list,num_list_d
integer,intent(in)::key_pp_paw,key_natpri_in,key_natpri_inps,key_natpri_out,key_jel_calc
real*8, intent(in)::eps,veta,psctoff,psftrad,filpp,rctpcc,chrjel,strjel,endjel
real*8, intent(in)::xmax,ymax,zmax,biasx,biasy,biasz
integer,intent(in)::indspe(natom),nlind(nprjmx,num_spe),noind(nprjmx,num_spe),nqct(num_spe),nqctpcc(num_spe)
integer,intent(in)::ntyppp(num_spe),nprj(num_spe),lpmx(num_spe),nradct(num_spe)
real*8, intent(in)::coef(0:nqmx,0:7,num_spe),cp(8,num_spe)
real*8, intent(in)::radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8, intent(in)::potc(nradmx,num_spe),awf(nradmx,nprmx*lmx,num_spe),pwf(nradmx,nprmx*lmx,num_spe)
real*8, intent(in)::atx(natom),aty(natom),atz(natom)
real*8, intent(in)::yylm(npoint,lrhomx),point(npoint),wt(npoint)
integer,intent(out)::natpri(natom),natprid(natom),natpri_inf(natom),natinf(natom),natinfd(natom),natinfd_vloc(natom)
integer,intent(out)::naps(natom),napsd(natom)
integer,intent(out)::lstvec2(num_list,num_ppcell),lstvecd2(num_list_d,num_ppcell_d)
integer,intent(out)::natx(natom),naty(natom),natz(natom),ndatx(natom),ndaty(natom),ndatz(natom)
integer,intent(out)::lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell)
integer,intent(out)::lstdx(num_list_d,num_ppcell_d),lstdy(num_list_d,num_ppcell_d),lstdz(num_list_d,num_ppcell_d)
real*8, intent(out)::qijl(nprjmx,nprjmx,lrhomx,natom),vloc_coarse(ncpx,ncpy,ncpz),vloc_dense(ncpx_d,ncpy_d,ncpz_d)
real*8, intent(out)::vcorer_all(nradmx,npoint,num_atcell),vcorer(nradmx,npoint,num_atcell)
real*8, intent(out)::vnlocp(num_list,nprjmx,num_ppcell),rhopcc_dense(ncpx_d,ncpy_d,ncpz_d),rhopccr(nradmx,npoint,num_atcell)
real*8, intent(out)::vloc_scw(ncpx,ncpy,ncpz,natom)
real*8, intent(out)::dvlocdx_scw(ncpx,ncpy,ncpz,natom)
real*8, intent(out)::dvlocdy_scw(ncpx,ncpy,ncpz,natom)
real*8, intent(out)::dvlocdz_scw(ncpx,ncpy,ncpz,natom)
real*8, intent(out)::vloc_hdp(num_list_d,num_ppcell_d)
real*8, intent(out)::dvlocdx_hdp(num_list_d,num_ppcell_d)
real*8, intent(out)::dvlocdy_hdp(num_list_d,num_ppcell_d)
real*8, intent(out)::dvlocdz_hdp(num_list_d,num_ppcell_d)
real*8, intent(out)::vjell(ncpz)
real*8, intent(out)::vboundx(-(nfh-1):nfh,ncpy*nmesh,ncpz*nmesh,9*(3-nperi)/3+1,(natom-1)*(3-nperi)/3+1)
real*8, intent(out)::vboundy(ncpx*nmesh,-(nfh-1):nfh,ncpz*nmesh,9*(3-nperi)/2+1,(natom-1)*(3-nperi)/2+1)
real*8, intent(out)::vboundz(ncpx*nmesh,ncpy*nmesh,-(nfh-1):nfh,9*(1-nperi/3)+1,(natom-1)*(1-nperi/3)+1)
integer nxspd,nyspd,nzspd
integer ix,iy,iz,na,la,lla,i,i1,i2,j,j1,j2,nf1
real*8 tmp,sum,sum_tru
real*8 dx,dy,dz
real*8 ddx,ddy,ddz
real*8,allocatable::wxyz(:)
real*8,allocatable::vtmp_pp(:,:,:)
real*8,allocatable::vtmp_pp1(:,:,:)
real*8,allocatable::vtmp_pp2(:,:,:)
real*8,allocatable::vnlspd(:,:,:,:)
real*8,allocatable::vtmp_cc(:,:,:)
real*8,allocatable::vtmp_dd(:,:,:)
real*8,allocatable::rrr(:,:,:),dderfr(:,:,:),ddexpr(:,:,:)

  dx=xmax/nxmax
  dy=ymax/nymax
  dz=zmax/nzmax
  ddx=dx/nmesh
  ddy=dy/nmesh
  ddz=dz/nmesh
  nf1=nf-1
  allocate(wxyz(-nfdg*npmesh:nfdg*npmesh-1))
  allocate(vtmp_cc(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf))
  allocate(vtmp_dd(-(nf*nmesh-1):ncpx_d+nf*nmesh,-(nf*nmesh-1):ncpy_d+nf*nmesh,-(nf*nmesh-1):ncpz_d+nf*nmesh))
  allocate(rrr(ncpx,ncpy,ncpz),dderfr(ncpx,ncpy,ncpz),ddexpr(ncpx,ncpy,ncpz))

! **********  compute double-grid weight  **********
! check for nfdg
  if ((nfdg .lt. 1) .and. (nfdg .gt. 8)) &
   call stopp('error in pseudo. nfdg should be between 1 and 8.')
!$omp parallel default(shared)
  call pseudocalc_01a(nfdg,npmesh,wxyz)
!$omp end parallel
! **************************************************

! **********  compute atom-information vector and list vectors  **********
!$omp parallel default(shared)
  call pseudocalc_01b(nperi,natom,num_spe,nradmx,nmesh,nfdg,ndisp,                                       & ! <
                      num_atcell,num_ppcell,num_ppcell_d,num_list,num_list_d,                            & ! <
                      nxmax,nymax,nzmax,ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                             & ! <
                      key_natpri_in,key_natpri_inps,key_natpri_out,                                      & ! <
                      psctoff,                                                                           & ! <
                      xmax,ymax,zmax,dx,dy,dz,                                                           & ! <
                      indspe,nradct,                                                                     & ! <
                      radial,                                                                            & ! <
                      atx,aty,atz,                                                                       & ! <
                      natx,naty,natz,ndatx,ndaty,ndatz,lstx,lsty,lstz,lstdx,lstdy,lstdz,                 & ! >
                      natinf,natinfd,natinfd_vloc,lstvec2,lstvecd2,natpri,natprid,natpri_inf,naps,napsd)   ! >
!$omp end parallel
!     ************************************************************************

!     **********  compute local potential on coarse grid (vloc_coarse) and dense grid (vloc_dense)  **********
!$omp parallel default(shared) private(i,ix,iy,iz)
!$omp do
  do i=1,ncpx*ncpy*ncpz
    vloc_coarse(i,1,1)=0.0d0
  end do
!$omp do
  do i=1,ncpx_d*ncpy_d*ncpz_d
    vloc_dense(i,1,1)=0.0d0
  end do

! ----------  compute hard local part on dense grid (vloc_hdp) and soft local part on coarse grid (vloc_s)  ----------
  call pseudocalc_02(natom,nradmx,num_spe,nperi,nmesh,nfiltyp,num_list_d,num_ppcell_d,nqmx, & ! <
                     jelcalc,nint1dmax,nzmax,                                               & ! <
                     ncpx,ncpy,ncpz,new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz,        & ! <
                     key_natpri_inps,key_jel_calc,                                          & ! <
                     xmax,ymax,zmax,dx,dy,dz,                                               & ! <
                     eps,psftrad,filpp,veta,chrjel,strjel,endjel,                           & ! <
                     indspe,nradct,natprid,napsd,natinfd_vloc,nqct,                         & ! <
                     lstdx,lstdy,lstdz,ndatx,ndaty,ndatz,                                   & ! <
                     cp,radial,coef,                                                        & ! <
                     atx,aty,atz,                                                           & ! <
                     vloc_scw,dvlocdx_scw,dvlocdy_scw,dvlocdz_scw,                          & ! >
                     vloc_hdp,dvlocdx_hdp,dvlocdy_hdp,dvlocdz_hdp,vjell,                    & ! >
                     rrr,dderfr,ddexpr)                                                       ! W

! ----------  compute the sum of soft local part on coarse grid  ----------
  call pseudocalc_03(natom,num_spe,nradmx,nf,nperi,        & ! <
                     ncpx,ncpy,ncpz,                       & ! <
                     xmax,ymax,zmax,dx,dy,dz,              & ! <
                     indspe,nradct,                        & ! <
                     vloc_scw,vjell,cp,radial,             & ! <
                     atx,aty,atz,                          & ! <
                     vtmp_cc)                                ! >
!$omp end parallel

! ----------  interpolate on dense grid and substitution vtmp_cc=>vtmp_dd=>vloc_dense  ----------
  call overlap_interpolation(ndisp,nperi,ncpx,ncpy,ncpz,nf,nf1,vtmp_cc,1)
  call interpolation(ncpx,ncpy,ncpz,nmesh,vtmp_dd,vtmp_cc,nperi,nf,1,ndisp)
!$omp parallel default(shared) private(i,ix,iy,iz)
!$omp do
  do iz=1,ncpz_d
  do iy=1,ncpy_d
  do ix=1,ncpx_d
    vloc_dense(ix,iy,iz)=vtmp_dd(ix,iy,iz)
  end do
  end do
  end do

! ----------  compute the sum of soft+hard local part on dense grid (vloc_dense)  ----------
  call pseudocalc_04(natom,num_list_d,num_ppcell_d,ncpx_d,ncpy_d,ncpz_d, & ! <
                     key_natpri_inps,                                    & ! <
                     natprid,napsd,natinfd_vloc,lstvecd2,                & ! <
                     vloc_hdp,                                           & ! <
                     vloc_dense)                                           ! X

! ----------  substitution vloc_dense=>vtmp_dd  ----------
!$omp do
  do i=1,(ncpx_d+2*nf*nmesh)*(ncpy_d+2*nf*nmesh)*(ncpz_d+2*nf*nmesh)
    vtmp_dd(-nf*nmesh+i,-nf*nmesh+1,-nf*nmesh+1)=0.0d0
  end do
!$omp do
  do iz=1,ncpz_d
  do iy=1,ncpy_d
  do ix=1,ncpx_d
    vtmp_dd(ix,iy,iz)=vloc_dense(ix,iy,iz)
  end do
  end do
  end do

! ----------  compute boundary condition for the isolated system  ----------
  if (nperi<3) then
  call pseudocalc_05(nperi,natom,num_spe,nf,nmesh,new_pwx,new_pwy,new_rsx,new_rsy,     & ! <
                     ncpx_d,ncpy_d,ncpz_d,nint1dmax,                                   & ! <
                     ddx,ddy,ddz,xmax,ymax,zmax,veta,                                  & ! <
                     indspe,                                                           & ! <
                     atx,aty,atz,                                                      & ! <
                     cp,                                                               & ! <
                     vtmp_dd)                                                            ! X
  end if
!$omp end parallel

! ----------  backwardly interpolate on coarse grid and substitution vtmp_dd=>vtmp_cc=>vloc_coarse  ----------
  call interpolation(ncpx,ncpy,ncpz,nmesh,vtmp_dd,vtmp_cc,nperi,nf,2,ndisp)
  call overlap_interpolation(ndisp,nperi,ncpx,ncpy,ncpz,nf,nf1,vtmp_cc,2)
!$omp parallel default(shared) private(i,ix,iy,iz)
!$omp do
  do iz=1,ncpz
  do iy=1,ncpy
  do ix=1,ncpx
    vloc_coarse(ix,iy,iz)=vtmp_cc(ix,iy,iz)
  end do
  end do
  end do
! ************************************************************************************

! **********  compute pcc charge (rhopcc_dense,rhopccr)  **********
  call pseudocalc_06(natom,nradmx,num_spe,npoint,num_atcell,nfiltyp,nqmx, & ! <
                     ncpx_d,ncpy_d,ncpz_d,                                & ! <
                     key_natpri_in,                                       & ! <
                     eps,psftrad,psctoff,filpp,rctpcc,                    & ! <
                     xmax,ymax,zmax,dx,dy,dz,ddx,ddy,ddz,                 & ! <
                     nradct,indspe,nqctpcc,natpri,natpri_inf,             & ! <
                     point,radial,coef,                                   & ! <
                     atx,aty,atz,                                         & ! <
                     rhopcc_dense,rhopccr)                                  ! >
!$omp end parallel
! *************************************************************

! **********  compute non-local part of pseudopotential (vnlocp)  **********
  nxspd=(npmesh*(npxmax+nfdg)-(npmesh+1)/2)
  nyspd=(npmesh*(npymax+nfdg)-(npmesh+1)/2)
  nzspd=(npmesh*(npzmax+nfdg)-(npmesh+1)/2)
  allocate(vtmp_pp(-npxmax+1:npxmax,-npymax+1:npymax,-npzmax+1:npzmax))
  allocate(vtmp_pp1(-nxspd+1:nxspd,-nyspd+1:nyspd,-npzmax+1:npzmax))
  allocate(vtmp_pp2(-nxspd+1:nxspd,-npymax+1:npymax,-npzmax+1:npzmax))
  allocate(  vnlspd(-nxspd+1:nxspd,-nyspd+1:nyspd,-nzspd+1:nzspd,nprjmx))
!$omp parallel default(shared)
  call pseudocalc_07(natom,num_spe,npmesh,nprjmx,nfdg,nfiltyp,num_list,num_ppcell, & ! <
                     nradmx,nqmx,                                                  & ! <
                     nxmax,nymax,nzmax,npxmax,npymax,npzmax,nxspd,nyspd,nzspd,     & ! <
                     key_natpri_in,key_natpri_inps,                                & ! <
                     eps,psftrad,psctoff,filpp,                                    & ! <
                     dx,dy,dz,                                                     & ! <
                     indspe,nprj,nradct,natpri,naps,natinf,nqct,                   & ! <
                     natx,naty,natz,lstx,lsty,lstz,                                & ! <
                     wxyz,radial,coef,                                             & ! <
                     atx,aty,atz,                                                  & ! <
                     vnlocp,                                                       & ! >
                     vtmp_pp,vtmp_pp1,vtmp_pp2,vnlspd)                               ! X
!$omp end parallel
  deallocate(vtmp_pp,vtmp_pp1,vtmp_pp2,vnlspd)
! **************************************************************************

! **********  compute local potential on radial grid (vloc^1(r))  **********
! check for lrhomx
  if (lrhomx>25) call stopp('error in computation of vloc^1(r): lrhomx must be <= 25!')
!$omp parallel default(shared)
  select case (nperi)
  case (0)
    call pseudocalc_08a(natom,num_spe,nradmx,npoint,num_atcell, & ! <
                        key_natpri_in,key_pp_paw,               & ! <
                        indspe,natpri,natpri_inf,ntyppp,nradct, & ! <
                        cp,radial,point,                        & ! <
                        atx,aty,atz,                            & ! <
                        vcorer)                                   ! >
!$omp barrier
  case (1:3)
! We can omit the computation for the contribution of the boundary value of the Coulomb potential,
! because it is cancelled by substracting \int (v_core) Q_{ij}^L dr from \tilda{D}_{ij}
! Definitions of Q_{ij}^L and \tilda{D}_{ij} are presented in PRB59 1758 (1999).
  call pseudocalc_08b(nperi,natom,num_spe,nradmx,npoint,lsphel,num_atcell,nint1dmax, & ! <
                      new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz,        & ! <
                      key_natpri_in,key_pp_paw,                               & ! <
                      xmax,ymax,zmax,veta,                                    & ! <
                      indspe,natpri,natpri_inf,ntyppp,nradct,                 & ! <
                      cp,radial,point,                                        & ! <
                      atx,aty,atz,                                            & ! <
                      vcorer)                                                   ! >
!$omp barrier
  end select
! ----------  compute the potential from the external electric field  ----------
  call pseudocalc_11(natom,num_spe,nradmx,npoint,num_atcell,nzmax,jelcalc, & ! <
                     key_natpri_in,key_pp_paw,key_jel_calc,                & ! <
                     indspe,ntyppp,nradct,natpri,natpri_inf,               & ! <
                     xmax,ymax,zmax,                                       & ! <
                     biasx,biasy,biasz,                                    & ! <
                     chrjel,strjel,endjel,                                 & ! <
                     point,radial,potc,                                    & ! <
                     vcorer,vcorer_all,                                    & ! X
                     atx,aty,atz)                                            ! <
!$omp end parallel
! **************************************************************************

! **********  compute moments of charges, \int Qij |r|^l Y(\hat{r}) dr in Eq. (26) of PRB59 1758 (1999)  **********
  qijl=0.0d0
  do na=1,natom
    if (ntyppp(indspe(na)) .eq. key_pp_paw) then
      do la=1,(2*lpmx(indspe(na))+1)**2
        if (la .eq. 1) lla=0
        if ((la .ge. 2) .and. (la .le. 4)) lla=1
        if ((la .ge. 5) .and. (la .le. 9)) lla=2
        if ((la .ge. 10) .and. (la .le. 16)) lla=3
        if ((la .ge. 17) .and. (la .le. 25)) lla=4
        do j=1,nprj(indspe(na))
        j1=nlind(j,indspe(na))
        j2=noind(j,indspe(na))
        do i=1,j
          i1=nlind(i,indspe(na))
          i2=noind(i,indspe(na))
          tmp=0.0d0
          sum=0.0d0
!$omp parallel default(shared)
          call pseudocalc_12(natom,nradmx,num_spe,nprmx,lmx,npoint,lrhomx, & ! <
                             na,la,lla,i1,i2,j1,j2,                        & ! <
                             tmp,sum,                                      & ! X
                             indspe,nradct,                                & ! <
                             radial,dradial,awf,pwf,yylm,wt)                 ! <
!$omp end parallel
          qijl(i,j,la,na)=sum*tmp
          qijl(j,i,la,na)=sum*tmp
        end do
        end do
      end do
    end if
  end do
! *****************************************************************************************************************

! **********  compute boundary of hartree potential  **********
!$omp parallel default(shared)
  call pseudocalc_13( &
   natom,num_spe,nperi,nfh,nint1dmax,new_pwx,new_pwy,new_rsx,new_rsy,ncpx_d,ncpy_d,ncpz_d, & ! <
   veta,ddx,ddy,ddz,xmax,ymax,zmax,                                                        & ! <
   indspe,atx,aty,atz,                                                                     & ! <
   vboundx,vboundy,vboundz)                                                                  ! >
!$omp end parallel
! *************************************************************

  deallocate(wxyz)
  deallocate(vtmp_cc,vtmp_dd)
  deallocate(rrr,dderfr,ddexpr)

  return
end subroutine pseudocalc


!this subroutine computes double-grid weight
subroutine pseudocalc_01a(nfdg,npmesh,wxyz)
implicit none
integer, intent(in)::nfdg,npmesh
real*8, intent(out)::wxyz(-nfdg*npmesh:nfdg*npmesh-1)
real*8 tmp,tmp1,r,rmf,rme,rmd,rmb,rmc,rma,rm9,rm8,rm7,rm6,rm5,rm4,rm3,rm2,rm1 &
      ,rp1,rp2,rp3,rp4,rp5,rp6,rp7
integer j,js,je,ntmp1

! **********  compute double-grid weight  **********
  js=-nfdg*npmesh
  je=nfdg*npmesh-1
  tmp1=0.5d0
  ntmp1=1
  if (mod(npmesh,2) .eq. 1) js=-nfdg*npmesh+1
  if (mod(npmesh,2) .eq. 1) tmp1=0.0d0
  if (mod(npmesh,2) .eq. 1) ntmp1=0
  if (nfdg .eq. 1) then
  tmp=1.0d0/npmesh
!$omp do
  do j=js,je
  r=dabs(j*tmp+tmp1*tmp)
  rm1=(r-1.0d0)
  wxyz(j)=rm1/(-1.0d0)
  end do
  end if
  if (nfdg .eq. 2) then
  tmp=1.0d0/npmesh
!$omp do
  do j=js,je
  r=dabs(j*tmp+tmp1*tmp)
  rm3=(r-3.0d0)
  rm2=(r-2.0d0)
  rm1=(r-1.0d0)
  rp1=(r+1.0d0)
  if (abs(2*j+ntmp1) .gt. 2*npmesh) then
    wxyz(j)=rm3*rm2*rm1/(-6.0d0)
  else
    wxyz(j)=rm2*rm1*rp1/(2.0d0)
  end if
  end do
  end if
  if (nfdg .eq. 3) then
  tmp=1.0d0/npmesh
!$omp do
  do j=js,je
  r=dabs(j*tmp+tmp1*tmp)
  rm5=(r-5.0d0)
  rm4=(r-4.0d0)
  rm3=(r-3.0d0)
  rm2=(r-2.0d0)
  rm1=(r-1.0d0)
  rp1=(r+1.0d0)
  rp2=(r+2.0d0)
  if (abs(2*j+ntmp1) .gt. 2*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*2*npmesh) then
    wxyz(j)=rm5*rm4*rm3*rm2*rm1/(-120.0d0)
  else
    wxyz(j)=rm4*rm3*rm2*rm1*rp1/(24.0d0)
  end if
  else
    wxyz(j)=rm3*rm2*rm1*rp1*rp2/(-12.0d0)
  end if
  end do
  end if
  if (nfdg .eq. 4) then
  tmp=1.0d0/npmesh
!$omp do
  do j=js,je
  r=dabs(j*tmp+tmp1*tmp)
  rm7=(r-7.0d0)
  rm6=(r-6.0d0)
  rm5=(r-5.0d0)
  rm4=(r-4.0d0)
  rm3=(r-3.0d0)
  rm2=(r-2.0d0)
  rm1=(r-1.0d0)
  rp1=(r+1.0d0)
  rp2=(r+2.0d0)
  rp3=(r+3.0d0)
  if (abs(2*j+ntmp1) .gt. 2*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*2*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*3*npmesh) then
    wxyz(j)=rm7*rm6*rm5*rm4*rm3*rm2*rm1/(-5040.0d0)
  else
    wxyz(j)=rm6*rm5*rm4*rm3*rm2*rm1*rp1/(720.0d0)
  end if
  else
    wxyz(j)=rm5*rm4*rm3*rm2*rm1*rp1*rp2/(-240.0d0)
  end if
  else
    wxyz(j)=rm4*rm3*rm2*rm1*rp1*rp2*rp3/(144.0d0)
  end if
  end do
  end if
  if (nfdg .eq. 5) then
  tmp=1.0d0/npmesh
!$omp do
  do j=js,je
  r=dabs(j*tmp+tmp1*tmp)
  rm9=(r-9.0d0)
  rm8=(r-8.0d0)
  rm7=(r-7.0d0)
  rm6=(r-6.0d0)
  rm5=(r-5.0d0)
  rm4=(r-4.0d0)
  rm3=(r-3.0d0)
  rm2=(r-2.0d0)
  rm1=(r-1.0d0)
  rp1=(r+1.0d0)
  rp2=(r+2.0d0)
  rp3=(r+3.0d0)
  rp4=(r+4.0d0)
  if (abs(2*j+ntmp1) .gt. 2*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*2*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*3*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*4*npmesh) then
    wxyz(j)=rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1/(-362880.0d0)
  else
    wxyz(j)=rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1/(40320.0d0)
  end if
  else
    wxyz(j)=rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2/(-10080.0d0)
  end if
  else
    wxyz(j)=rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2*rp3/(4320.0d0)
  end if
  else
    wxyz(j)=rm5*rm4*rm3*rm2*rm1*rp1*rp2*rp3*rp4/(-2880.0d0)
  end if
  end do
  end if
  if (nfdg .eq. 6) then
  tmp=1.0d0/npmesh
!$omp do
  do j=js,je
  r=dabs(j*tmp+tmp1*tmp)
  rmb=(r-11.0d0)
  rma=(r-10.0d0)
  rm9=(r-9.0d0)
  rm8=(r-8.0d0)
  rm7=(r-7.0d0)
  rm6=(r-6.0d0)
  rm5=(r-5.0d0)
  rm4=(r-4.0d0)
  rm3=(r-3.0d0)
  rm2=(r-2.0d0)
  rm1=(r-1.0d0)
  rp1=(r+1.0d0)
  rp2=(r+2.0d0)
  rp3=(r+3.0d0)
  rp4=(r+4.0d0)
  rp5=(r+5.0d0)
  if (abs(2*j+ntmp1) .gt. 2*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*2*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*3*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*4*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*5*npmesh) then
    wxyz(j)=rmb*rma*rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1/(-39916800.0d0)
  else
    wxyz(j)=rma*rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1/(3628800.0d0)
  end if
  else
    wxyz(j)=rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2/(-725760.0d0)
  end if
  else
    wxyz(j)=rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2*rp3/(241920.0d0)
  end if
  else
    wxyz(j)=rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2*rp3*rp4/(-120960.0d0)
  end if
  else
    wxyz(j)=rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2*rp3*rp4*rp5/(86400.0d0)
  end if
  end do
  end if
  if (nfdg .eq. 7) then
  tmp=1.0d0/npmesh
!$omp do
  do j=js,je
  r=dabs(j*tmp+tmp1*tmp)
  rmd=(r-13.0d0)
  rmc=(r-12.0d0)
  rmb=(r-11.0d0)
  rma=(r-10.0d0)
  rm9=(r-9.0d0)
  rm8=(r-8.0d0)
  rm7=(r-7.0d0)
  rm6=(r-6.0d0)
  rm5=(r-5.0d0)
  rm4=(r-4.0d0)
  rm3=(r-3.0d0)
  rm2=(r-2.0d0)
  rm1=(r-1.0d0)
  rp1=(r+1.0d0)
  rp2=(r+2.0d0)
  rp3=(r+3.0d0)
  rp4=(r+4.0d0)
  rp5=(r+5.0d0)
  rp6=(r+6.0d0)
  if (abs(2*j+ntmp1) .gt. 2*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*2*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*3*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*4*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*5*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*6*npmesh) then
    wxyz(j)=rmd*rmc*rmb*rma*rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1/(-6227020800.0d0)
  else
    wxyz(j)=rmc*rmb*rma*rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1/(479001600.0d0)
  end if
  else
    wxyz(j)=rmb*rma*rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2/(-79833600.0d0)
  end if
  else
    wxyz(j)=rma*rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2*rp3/(21772800.0d0)
  end if
  else
    wxyz(j)=rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2*rp3*rp4/(-8709120.0d0)
  end if
  else
    wxyz(j)=rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2*rp3*rp4*rp5/(4838400.0d0)
  end if
  else
    wxyz(j)=rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2*rp3*rp4*rp5*rp6/(-3628800.0d0)
  end if
  end do
  end if
  if (nfdg .eq. 8) then
  tmp=1.0d0/npmesh
!$omp do
  do j=js,je
  r=dabs(j*tmp+tmp1*tmp)
  rmf=(r-15.0d0)
  rme=(r-14.0d0)
  rmd=(r-13.0d0)
  rmc=(r-12.0d0)
  rmb=(r-11.0d0)
  rma=(r-10.0d0)
  rm9=(r-9.0d0)
  rm8=(r-8.0d0)
  rm7=(r-7.0d0)
  rm6=(r-6.0d0)
  rm5=(r-5.0d0)
  rm4=(r-4.0d0)
  rm3=(r-3.0d0)
  rm2=(r-2.0d0)
  rm1=(r-1.0d0)
  rp1=(r+1.0d0)
  rp2=(r+2.0d0)
  rp3=(r+3.0d0)
  rp4=(r+4.0d0)
  rp5=(r+5.0d0)
  rp6=(r+6.0d0)
  rp7=(r+7.0d0)
  if (abs(2*j+ntmp1) .gt. 2*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*2*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*3*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*4*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*5*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*6*npmesh) then
  if (abs(2*j+ntmp1) .gt. 2*7*npmesh) then
    wxyz(j)=rmf*rme*rmd*rmc*rmb*rma*rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1/(-1307674368000.0d0)
  else
    wxyz(j)=rme*rmd*rmc*rmb*rma*rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1/(87178291200.0d0)
  end if
  else
    wxyz(j)=rmd*rmc*rmb*rma*rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2/(-12454041600.0d0)
  end if
  else
    wxyz(j)=rmc*rmb*rma*rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2*rp3/(2874009600.0d0)
  end if
  else
    wxyz(j)=rmb*rma*rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2*rp3*rp4/(-958003200.0d0)
  end if
  else
    wxyz(j)=rma*rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2*rp3*rp4*rp5/(435456000.0d0)
  end if
  else
    wxyz(j)=rm9*rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2*rp3*rp4*rp5*rp6/(-261273600.0d0)
  end if
  else
    wxyz(j)=rm8*rm7*rm6*rm5*rm4*rm3*rm2*rm1*rp1*rp2*rp3*rp4*rp5*rp6*rp7/(203212800.0d0)
  end if
  end do
  end if
! **************************************************
  return
end subroutine pseudocalc_01a


!this subroutine computes  atom-information vector and list vectors
subroutine pseudocalc_01b(nperi,natom,num_spe,nradmx,nmesh,nfdg,ndisp,num_atcell,num_ppcell,num_ppcell_d,num_list,      & ! <
                          num_list_d,nxmax,nymax,nzmax,ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                             & ! <
                          key_natpri_in,key_natpri_inps,key_natpri_out,                                                 & ! <
                          psctoff,                                                                                      & ! <
                          xmax,ymax,zmax,dx,dy,dz,                                                                      & ! <
                          indspe,nradct,                                                                                & ! <
                          radial,                                                                                       & ! <
                          atx,aty,atz,                                                                                  & ! <
                          natx,naty,natz,ndatx,ndaty,ndatz,lstx,lsty,lstz,lstdx,lstdy,lstdz,                            & ! >
                          natinf,natinfd,natinfd_vloc,lstvec2,lstvecd2,natpri,natprid,natpri_inf,naps,napsd)              ! >
use mod_stopp
use mod_mpi
implicit none
integer, intent(in)::nperi,natom,num_spe,nradmx,nmesh,nfdg,ndisp
integer, intent(in)::num_atcell,num_ppcell,num_ppcell_d,num_list,num_list_d
integer, intent(in)::nxmax,nymax,nzmax,ncpx,ncpy,ncpz,npxmax,npymax,npzmax
integer, intent(in)::key_natpri_in,key_natpri_inps,key_natpri_out
real*8, intent(in)::xmax,ymax,zmax,dx,dy,dz
real*8, intent(in)::psctoff
integer, intent(in)::indspe(natom),nradct(num_spe)
integer, intent(out)::natinf(natom),natinfd(natom),natinfd_vloc(natom)
integer, intent(out)::lstvec2(num_list,num_ppcell),lstvecd2(num_list_d,num_ppcell_d)
integer, intent(out)::natpri(natom),natprid(natom),natpri_inf(natom),naps(natom),napsd(natom)
integer, intent(out)::natx(natom),naty(natom),natz(natom)
integer, intent(out)::ndatx(natom),ndaty(natom),ndatz(natom)
integer, intent(out)::lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell)
integer, intent(out)::lstdx(num_list_d,num_ppcell_d),lstdy(num_list_d,num_ppcell_d),lstdz(num_list_d,num_ppcell_d)
real*8, intent(in)::radial(nradmx,num_spe)
real*8, intent(in)::atx(natom),aty(natom),atz(natom)
real*8 r,atmtmpx,atmtmpy,atmtmpz,x,y,z
integer i,j,k,l,na,l0,l2,l3,ix,iy,iz,jx,jy,jz,iaps,iapsd,nx,ny,nz
integer ncpx_d,ncpy_d,ncpz_d,npxmax_d,npymax_d,npzmax_d
real*8 ddx,ddy,ddz
logical :: lspec
  ncpx_d=ncpx*nmesh
  ncpy_d=ncpy*nmesh
  ncpz_d=ncpz*nmesh
  npxmax_d=npxmax*nmesh
  npymax_d=npymax*nmesh
  npzmax_d=npzmax*nmesh
  ddx=dx/nmesh
  ddy=dy/nmesh
  ddz=dz/nmesh

! **********  compute atom-information vector  **********
!$omp do
  do na=1,natom
    ! make sure that atom positions were shifted correctly
    if  (( atx(na)-0.5d0*dx < -(xmax+0.1d0*ddx) ).or.( atx(na)-0.5d0*dx > +(xmax+0.1d0*ddx) ) &
     .or.( aty(na)-0.5d0*dy < -(ymax+0.1d0*ddy) ).or.( aty(na)-0.5d0*dy > +(ymax+0.1d0*ddy) ) &
     .or.( atz(na)-0.5d0*dx < -(zmax+0.1d0*ddz) ).or.( atz(na)-0.5d0*dz > +(zmax+0.1d0*ddz) )) then
!      call stopp('pseudocalc_01b: atom position out of bounds!')
      if (myrank_glbl==0) write(ndisp,'(a,i6,a)') 'pseudocalc_01b: atom position',na,'out of bounds!'
    endif
    ! determine nat,ndnat
    natx(na)=int((atx(na)+xmax+0.5d0*dx)/dx)-nxmax
    naty(na)=int((aty(na)+ymax+0.5d0*dy)/dy)-nymax
    natz(na)=int((atz(na)+zmax+0.5d0*dz)/dz)-nzmax
    ndatx(na)=int((atx(na)+xmax+0.5d0*ddx)/ddx)-nxmax*nmesh
    ndaty(na)=int((aty(na)+ymax+0.5d0*ddy)/ddy)-nymax*nmesh
    ndatz(na)=int((atz(na)+zmax+0.5d0*ddz)/ddz)-nzmax*nmesh
    ! take possible rounding errors into account:
!    natx(na)= min(natx(na),nxmax)
!    naty(na)= min(naty(na),nymax)
!    natz(na)= min(natz(na),nzmax)
!    natx(na)= max(natx(na),1-nxmax)
!    naty(na)= max(naty(na),1-nymax)
!    natz(na)= max(natz(na),1-nzmax)
  end do
!$omp do
  do na=1,natom
    i= (natx(na)+nxmax-1)/ncpx
    j= (naty(na)+nymax-1)/ncpy
    k= (natz(na)+nzmax-1)/ncpz
    if ((i==myrx).and.(j==myry).and.(k==myrz)) then
      natpri(na)=key_natpri_in
    else
      natpri(na)=key_natpri_out ! might be changed below
    end if
  end do

!$omp single
  l=0
  natpri_inf(:)=0
  do na=1,natom
    if (natpri(na) .eq. key_natpri_in) then
      l=l+1
      natpri_inf(na)=l
      if (natpri_inf(na) .gt. num_atcell) call stopp('pseudocalc_01b: num_atcell is too small.')
    end if
  end do
  call mpi_reduce(l,na,1,mpi_integer,mpi_max,0,mpicom_space,mpij)
  if (myrank_glbl==0) write(ndisp,*) 'num_atcell   must be >=',na
!$omp end single
! *******************************************************

! **********  compute atom-position informations  **********
  lspec=.false.
  do iz=1,2
  do iy=1,2
  do ix=1,2
  if ((nprocx>=4) .and. (nprocy>=4) .and. (nprocz>=4) .and. &
      (2*npxmax<=ix*ncpx) .and. (2*npymax<=iy*ncpy) .and. (2*npzmax<=iz*ncpz) .and. &
      (mod(nprocx,2*ix)==0) .and. (mod(nprocy,2*iy)==0) .and. (mod(nprocz,2*iz)==0) .and. (.not. lspec)) then
!$omp single
if (myrank_glbl==0) write(ndisp,'(a,3i1,a)') 'special routine for mpi_communicator ',2*ix,2*iy,2*iz,' is used. See pseudocalc_01b.'
!$omp end single
    nx=ix*2
    ny=iy*2
    nz=iz*2
    lspec=.true.
  end if
  end do
  end do
  end do
  if (.not. lspec) then
!$omp single
    if (myrank_glbl==0) write(ndisp,'(a,3i1,a)') 'mpi_communicator is assigned to each atom. See pseudocalc_01b.'
!$omp end single
  end if
  if (lspec) then
!$omp do
    do na=1,natom
      if (natpri(na)==key_natpri_out) then
        i= (natx(na)+nxmax-1-ncpx/(4/nx))
        j= (naty(na)+nymax-1-ncpy/(4/ny))
        k= (natz(na)+nzmax-1-ncpz/(4/nz))
        if (nperi >= 1) then
          do while (i>2*nxmax)
            i=i-2*nxmax
          end do
          do while (i<1)
            i=i+2*nxmax
          end do
        if (nperi >= 2) then
          do while (j>2*nymax)
            j=j-2*nymax
          end do
          do while (j<1)
            j=j+2*nymax
          end do
        if (nperi == 3) then
          do while (k>2*nzmax)
            k=k-2*nzmax
          end do
          do while (k<1)
            k=k+2*nzmax
          end do
        end if
        end if
        end if
        i=(i/((nx/2)*ncpx))*(nx/2)
        j=(j/((ny/2)*ncpy))*(ny/2)
        k=(k/((nz/2)*ncpz))*(nz/2)
        do iz=0,nz-1
        do iy=0,ny-1
        do ix=0,nx-1
          jx=i+ix
          jy=j+iy
          jz=k+iz
          if (jx .ge. nprocx) jx=jx-nprocx
          if (jy .ge. nprocy) jy=jy-nprocy
          if (jz .ge. nprocz) jz=jz-nprocz
          if ((jx==myrx).and.(jy==myry).and.(jz==myrz)) natpri(na)=key_natpri_inps
        end do
        end do
        end do
      end if
    end do
  else
!$omp do
    do na=1,natom
      atmtmpx=atx(na)-(natx(na)*dx-0.5d0*dx)
      atmtmpy=aty(na)-(naty(na)*dy-0.5d0*dy)
      atmtmpz=atz(na)-(natz(na)*dz-0.5d0*dz)
      l0=0
      do iz=-npzmax+1,npzmax
      do iy=-npymax+1,npymax
      do ix=-npxmax+1,npxmax
        if (natpri(na) .eq. key_natpri_out) then
        l=0
        do jz=-nfdg,nfdg
        do jy=-nfdg,nfdg
        do jx=-nfdg,nfdg
          x=(ix+jx)*dx-atmtmpx
          y=(iy+jy)*dy-atmtmpy
          z=(iz+jz)*dz-atmtmpz
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) l=1
        end do
        end do
        end do
        if (l .eq. 1) then
          i=ix+natx(na)+nxmax-myrx*ncpx
          j=iy+naty(na)+nymax-myry*ncpy
          k=iz+natz(na)+nzmax-myrz*ncpz
          if (nperi >= 1) then
            do while (i>2*nxmax)
              i=i-2*nxmax
            end do
            do while (i<1)
              i=i+2*nxmax
            end do
          if (nperi >= 2) then
            do while (j>2*nymax)
              j=j-2*nymax
            end do
            do while (j<1)
              j=j+2*nymax
            end do
          if (nperi == 3) then
            do while (k>2*nzmax)
              k=k-2*nzmax
            end do
            do while (k<1)
              k=k+2*nzmax
            end do
          end if
          end if
          end if
          if (((i-1)*(ncpx-i).ge.0) .and. ((j-1)*(ncpy-j).ge.0) .and. ((k-1)*(ncpz-k).ge.0)) then
            l0=1
          end if
        end if
        end if
      end do
      end do
      end do
      if ((natpri(na) .eq. key_natpri_out) .and. (l0 .eq. 1)) natpri(na)=key_natpri_inps
    end do
  end if

!$omp single
  iaps=0
  naps=0
  do na=1,natom
    if ((natpri(na) .eq. key_natpri_in) .or. (natpri(na) .eq. key_natpri_inps)) then
      iaps=iaps+1
      naps(na)=iaps
      if (naps(na) .gt. num_ppcell) then
        write(ndisp,*) 'error on rank=',myr_space,'naps=',naps(na)
        call stopp('error! num_ppcell is too small.')
      end if
    end if
  end do
  call mpi_allreduce(iaps,na,1,mpi_integer,mpi_max,mpicom_space,mpij)
  if (myrank_glbl==0) write(ndisp,*) 'num_ppcell=',na
!$omp end single
! **********************************************************

! **********  compute list vector for pseudopotential  **********
!$omp do
  do na=1,natom
    if ((natpri(na) .eq. key_natpri_in) .or. (natpri(na) .eq. key_natpri_inps)) then
      iaps=naps(na)
      atmtmpx=atx(na)-(natx(na)*dx-0.5d0*dx)
      atmtmpy=aty(na)-(naty(na)*dy-0.5d0*dy)
      atmtmpz=atz(na)-(natz(na)*dz-0.5d0*dz)
      lstvec2(:,iaps)=0
      lstx(:,iaps)=0
      lsty(:,iaps)=0
      lstz(:,iaps)=0
      l0=0
      do iz=-npzmax+1,npzmax
      do iy=-npymax+1,npymax
      do ix=-npxmax+1,npxmax
        l=0
        do jz=-nfdg,nfdg
        do jy=-nfdg,nfdg
        do jx=-nfdg,nfdg
          x=(ix+jx)*dx-atmtmpx
          y=(iy+jy)*dy-atmtmpy
          z=(iz+jz)*dz-atmtmpz
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) l=1
        end do
        end do
        end do
        if (l .eq. 1) then
          i=ix+natx(na)+nxmax-myrx*ncpx
          j=iy+naty(na)+nymax-myry*ncpy
          k=iz+natz(na)+nzmax-myrz*ncpz
          if (nperi >= 1) then
            do while (i>2*nxmax)
              i=i-2*nxmax
            end do
            do while (i<1)
              i=i+2*nxmax
            end do
          if (nperi >= 2) then
            do while (j>2*nymax)
              j=j-2*nymax
            end do
            do while (j<1)
              j=j+2*nymax
            end do
          if (nperi == 3) then
            do while (k>2*nzmax)
              k=k-2*nzmax
            end do
            do while (k<1)
              k=k+2*nzmax
            end do
          end if
          end if
          end if
          if (((i-1)*(ncpx-i).ge.0) .and. ((j-1)*(ncpy-j).ge.0) .and. ((k-1)*(ncpz-k).ge.0)) then
            l0=l0+1
            l2=(k-1   )* ncpy      * ncpx      +(j-1   )* ncpx      +i
            l3=(iz+npzmax-1)*4*npymax*npxmax+(iy+npymax-1)*2*npxmax+ix+npxmax
            lstvec2(l0,iaps)=l2
            lstx(l0,iaps)=ix
            lsty(l0,iaps)=iy
            lstz(l0,iaps)=iz
          end if
        end if
      end do
      end do
      end do
      natinf(na)=l0

      if (natinf(na) > num_list) then
        write(ndisp,*) 'error on rank=',myr_space,'na=',na,'num_list=',num_list
        call stopp('error! num_list is too small.')
      end if

    else
      natinf(na)=0
    end if
  end do
! ***************************************************************

! **********  check whether the non-local parts of pseudopotential is involved in the subdomain  **********
! key_natpri_inps : ps is involved.
! key_natpri_out  : ps is NOT involved.
!$omp do
  do na=1,natom
    natprid(na)=key_natpri_out
    l0=0
    do iz=-npzmax_d+1,npzmax_d
    do iy=-npymax_d+1,npymax_d
    do ix=-npxmax_d+1,npxmax_d
      i=ix+ndatx(na)+nxmax*nmesh-myrx*ncpx_d
      j=iy+ndaty(na)+nymax*nmesh-myry*ncpy_d
      k=iz+ndatz(na)+nzmax*nmesh-myrz*ncpz_d
      x=ix*ddx-(atx(na)-(ndatx(na)*ddx-0.5d0*ddx))
      y=iy*ddy-(aty(na)-(ndaty(na)*ddy-0.5d0*ddy))
      z=iz*ddz-(atz(na)-(ndatz(na)*ddz-0.5d0*ddz))
      r=dsqrt(x*x+y*y+z*z)
      if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
        if (nperi >= 1) then
          do while (i>2*nxmax*nmesh)
            i=i-2*nxmax*nmesh
          end do
          do while (i<1)
            i=i+2*nxmax*nmesh
          end do
        if (nperi >= 2) then
          do while (j>2*nymax*nmesh)
            j=j-2*nymax*nmesh
          end do
          do while (j<1)
            j=j+2*nymax*nmesh
          end do
        if (nperi == 3) then
          do while (k>2*nzmax*nmesh)
            k=k-2*nzmax*nmesh
          end do
          do while (k<1)
            k=k+2*nzmax*nmesh
          end do
        end if
        end if
        end if
        if (((i-1)*(ncpx_d-i).ge.0) .and. ((j-1)*(ncpy_d-j).ge.0) .and. ((k-1)*(ncpz_d-k).ge.0)) then
          l0=1
        end if
      end if
    end do
    end do
    end do
    if (l0 .eq. 1) natprid(na)=key_natpri_inps
  end do

!$omp single
  iapsd=0
  napsd=0
  do na=1,natom
    if (natprid(na) .eq. key_natpri_inps) then
      iapsd=iapsd+1
      napsd(na)=iapsd
      if (napsd(na) .gt. num_ppcell_d) then
        write(ndisp,*) 'error on rank=',myr_space,'napsd=',napsd(na)
        call stopp('error! num_ppcell_d is too small.')
      end if
    end if
  end do
  call mpi_allreduce(iapsd,na,1,mpi_integer,mpi_max,mpicom_space,mpij)
  if (myrank_glbl==0) write(ndisp,*) 'num_ppcell_d=',na
!$omp end single
! *********************************************************************************************************

! **********  compute list vector for hard local pot. and comp. charge  **********
!$omp do
  do na=1,natom
    if (natprid(na) .eq. key_natpri_inps) then
      iapsd=napsd(na)
      lstvecd2(:,iapsd)=0
      lstdx(:,iapsd)=0
      lstdy(:,iapsd)=0
      lstdz(:,iapsd)=0
      l0=0
      do iz=-npzmax_d+1,npzmax_d
      do iy=-npymax_d+1,npymax_d
      do ix=-npxmax_d+1,npxmax_d
        i=ix+ndatx(na)+nxmax*nmesh-myrx*ncpx_d
        j=iy+ndaty(na)+nymax*nmesh-myry*ncpy_d
        k=iz+ndatz(na)+nzmax*nmesh-myrz*ncpz_d
        x=ix*ddx-(atx(na)-(ndatx(na)*ddx-0.5d0*ddx))
        y=iy*ddy-(aty(na)-(ndaty(na)*ddy-0.5d0*ddy))
        z=iz*ddz-(atz(na)-(ndatz(na)*ddz-0.5d0*ddz))
        r=dsqrt(x*x+y*y+z*z)
        if (r .lt. radial(nradct(indspe(na)),indspe(na))) then
        if (nperi >= 1) then
          do while (i>2*nxmax*nmesh)
            i=i-2*nxmax*nmesh
          end do
          do while (i<1)
            i=i+2*nxmax*nmesh
          end do
        if (nperi >= 2) then
          do while (j>2*nymax*nmesh)
            j=j-2*nymax*nmesh
          end do
          do while (j<1)
            j=j+2*nymax*nmesh
          end do
        if (nperi == 3) then
          do while (k>2*nzmax*nmesh)
            k=k-2*nzmax*nmesh
          end do
          do while (k<1)
            k=k+2*nzmax*nmesh
          end do
        end if
        end if
        end if
          if (((i-1)*(ncpx_d-i).ge.0) .and. ((j-1)*(ncpy_d-j).ge.0) .and. ((k-1)*(ncpz_d-k).ge.0)) then
            l0=l0+1
            l2=(k-1   )* ncpy_d      * ncpx_d      +(j-1   )* ncpx_d      +i
            l3=(iz+npzmax_d-1)*4*npymax_d*npxmax_d+(iy+npymax_d-1)*2*npxmax_d+ix+npxmax_d
            lstvecd2(l0,iapsd)=l2
            lstdx(l0,iapsd)=ix
            lstdy(l0,iapsd)=iy
            lstdz(l0,iapsd)=iz
          end if
        end if
      end do
      end do
      end do
      natinfd(na)=l0
      do iz=-npzmax_d+1,npzmax_d
      do iy=-npymax_d+1,npymax_d
      do ix=-npxmax_d+1,npxmax_d
        i=ix+ndatx(na)+nxmax*nmesh-myrx*ncpx_d
        j=iy+ndaty(na)+nymax*nmesh-myry*ncpy_d
        k=iz+ndatz(na)+nzmax*nmesh-myrz*ncpz_d
        x=ix*ddx-(atx(na)-(ndatx(na)*ddx-0.5d0*ddx))
        y=iy*ddy-(aty(na)-(ndaty(na)*ddy-0.5d0*ddy))
        z=iz*ddz-(atz(na)-(ndatz(na)*ddz-0.5d0*ddz))
        r=dsqrt(x*x+y*y+z*z)
        if ((r .ge. radial(nradct(indspe(na)),indspe(na))) &
          .and. (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff)) then
          if (nperi >= 1) then
            do while (i>2*nxmax*nmesh)
              i=i-2*nxmax*nmesh
            end do
            do while (i<1)
              i=i+2*nxmax*nmesh
            end do
          if (nperi >= 2) then
            do while (j>2*nymax*nmesh)
              j=j-2*nymax*nmesh
            end do
            do while (j<1)
              j=j+2*nymax*nmesh
            end do
          if (nperi == 3) then
            do while (k>2*nzmax*nmesh)
              k=k-2*nzmax*nmesh
            end do
            do while (k<1)
              k=k+2*nzmax*nmesh
            end do
          end if
          end if
          end if
          if (((i-1)*(ncpx_d-i).ge.0) .and. ((j-1)*(ncpy_d-j).ge.0) .and. ((k-1)*(ncpz_d-k).ge.0)) then
            l0=l0+1
            l2=(k-1   )* ncpy_d      * ncpx_d      +(j-1   )* ncpx_d      +i
            l3=(iz+npzmax_d-1)*4*npymax_d*npxmax_d+(iy+npymax_d-1)*2*npxmax_d+ix+npxmax_d
            lstvecd2(l0,iapsd)=l2
            lstdx(l0,iapsd)=ix
            lstdy(l0,iapsd)=iy
            lstdz(l0,iapsd)=iz
          end if
        end if
      end do
      end do
      end do
      natinfd_vloc(na)=l0

      if (natinfd_vloc(na) .gt. num_list_d) then
        write(ndisp,*) 'error on rank=',myr_space,'na=',na,'num_list_d=',num_list_d
        call stopp('error! num_list_d is too small.')
      end if

    else
      natinfd(na)=0
      natinfd_vloc(na)=0
    end if
  end do
! ********************************************************************************

  return
end subroutine pseudocalc_01b


!this subroutine computes hard local part on dense grid (vloc_hdp) and soft local part on coarse grid (vloc_s)
subroutine pseudocalc_02(natom,nradmx,num_spe,nperi,nmesh,nfiltyp,num_list_d,num_ppcell_d,nqmx, & ! <
                         jelcalc,nint1dmax,nzmax,                                               & ! <
                         ncpx,ncpy,ncpz,new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz,        & ! <
                         key_natpri_inps,key_jel_calc,                                          & ! <
                         xmax,ymax,zmax,dx,dy,dz,                                               & ! <
                         eps,psftrad,filpp,veta,chrjel,strjel,endjel,                           & ! <
                         indspe,nradct,natprid,napsd,natinfd_vloc,nqct,                         & ! <
                         lstdx,lstdy,lstdz,ndatx,ndaty,ndatz,                                   & ! <
                         cp,radial,coef,                                                        & ! <
                         atx,aty,atz,                                                           & ! <
                         vloc_scw,dvlocdx_scw,dvlocdy_scw,dvlocdz_scw,                          & ! >
                         vloc_hdp,dvlocdx_hdp,dvlocdy_hdp,dvlocdz_hdp,vjell,                    & ! >
                         rrr,dderfr,ddexpr)                                                       ! W
use mod_mpi,        only: myrx,myry,myrz
use mod_mathfunctions, only: fermidis,expint1
implicit none
integer, intent(in)::natom,nradmx,num_spe,nperi,nmesh,nfiltyp,num_list_d,num_ppcell_d,nqmx
integer, intent(in)::jelcalc,nint1dmax,nzmax
integer, intent(in)::ncpx,ncpy,ncpz,new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz
integer, intent(in)::key_natpri_inps,key_jel_calc
integer, intent(in)::indspe(natom),nradct(num_spe),natprid(natom),napsd(natom),natinfd_vloc(natom),nqct(num_spe)
integer, intent(in)::lstdx(num_list_d,num_ppcell_d),lstdy(num_list_d,num_ppcell_d),lstdz(num_list_d,num_ppcell_d)
integer, intent(in)::ndatx(natom),ndaty(natom),ndatz(natom)
real*8, intent(in)::xmax,ymax,zmax,dx,dy,dz
real*8, intent(in)::eps,psftrad,filpp,veta,chrjel,strjel,endjel
real*8, intent(in)::cp(8,num_spe),radial(nradmx,num_spe),coef(0:nqmx,0:7,num_spe)
real*8, intent(in)::atx(natom),aty(natom),atz(natom)
real*8, intent(out)::vloc_scw(ncpx,ncpy,ncpz,natom)
real*8, intent(out)::dvlocdx_scw(ncpx,ncpy,ncpz,natom)
real*8, intent(out)::dvlocdy_scw(ncpx,ncpy,ncpz,natom)
real*8, intent(out)::dvlocdz_scw(ncpx,ncpy,ncpz,natom)
real*8, intent(out)::vloc_hdp(num_list_d,num_ppcell_d)
real*8, intent(out)::dvlocdx_hdp(num_list_d,num_ppcell_d)
real*8, intent(out)::dvlocdy_hdp(num_list_d,num_ppcell_d)
real*8, intent(out)::dvlocdz_hdp(num_list_d,num_ppcell_d)
real*8, intent(out)::vjell(ncpz)
real*8, intent(inout)::rrr(ncpx,ncpy,ncpz),dderfr(ncpx,ncpy,ncpz),ddexpr(ncpx,ncpy,ncpz)
real*8, allocatable::anum0(:,:,:),anum1(:,:,:)

real*8 pi,twopicbin,omega,omegain,surf,surfin,fourpi,twosqrtpiin,pi05,cp1,rcutss,dqq,tmp, &
       x,y,z,r,rin,tmppo,tmpdx,tmpdy,tmpdz,qqq,qqq2,qqqr,qqqrin,sinqqqr,cosqqqr,besselj0,dbesselj0,dbesselj0qq, &
       vkx,vky,vkz,vk2,tmpdsin,r2in,derfdrr,vlo,vlo0,vlo1,ta,tb,deftp1,deftp2,vep0,vep1,vsine,vcosi,awei,ddx,ddy,ddz,weight,dt,t
real*8 derf
integer na,ix,iy,iz,iapsd,i,ixyz,iq,jx,jy,jz,kx,ky,kz,kpx,kpy,kpz,ierr

  if (nperi==1) allocate(anum0(-new_pwx:new_pwx,ncpy,ncpz),anum1(-new_pwx:new_pwx,ncpy,ncpz))
  pi=dacos(-1.0d0)
  twopicbin=1.0d0/(2.0d0*pi)**3
  omega=xmax*ymax*zmax*8.0d0
  omegain=1.0d0/omega
  surf=xmax*ymax*4.0d0
  surfin=1.0d0/surf
  fourpi=4.0d0*pi
  pi05=dsqrt(pi)
  twosqrtpiin=2.0d0/dsqrt(pi)
  ddx=dx/nmesh
  ddy=dy/nmesh
  ddz=dz/nmesh

!  ----------  zero clear  ----------
!$omp do
  do i=1,ncpx*ncpy*ncpz*natom
       vloc_scw(i,1,1,1)=0.0d0
    dvlocdx_scw(i,1,1,1)=0.0d0
    dvlocdy_scw(i,1,1,1)=0.0d0
    dvlocdz_scw(i,1,1,1)=0.0d0
  end do
!$omp do
  do i=1,num_list_d*num_ppcell_d
       vloc_hdp(i,1)=0.0d0
    dvlocdx_hdp(i,1)=0.0d0
    dvlocdy_hdp(i,1)=0.0d0
    dvlocdz_hdp(i,1)=0.0d0
  end do
!      ----------------------------------

  do na=1,natom
    cp1=cp(1,indspe(na))
    rcutss=radial(nradct(indspe(na)),indspe(na))*psftrad

!  ----------  compute hard local part (vloc_h,dvlocd*_h)  ----------
    if (natprid(na) .eq. key_natpri_inps) then
      iapsd=napsd(na)
      dqq=pi/rcutss
      tmp=fourpi*dqq*twopicbin
!$omp do
      do ixyz=1,natinfd_vloc(na)
        ix=lstdx(ixyz,iapsd)
        iy=lstdy(ixyz,iapsd)
        iz=lstdz(ixyz,iapsd)
        x=ix*ddx-(atx(na)-(ndatx(na)*ddx-0.5d0*ddx))
        y=iy*ddy-(aty(na)-(ndaty(na)*ddy-0.5d0*ddy))
        z=iz*ddz-(atz(na)-(ndatz(na)*ddz-0.5d0*ddz))
        r=dsqrt(x*x+y*y+z*z)
        tmppo=0.0d0
        tmpdx=0.0d0
        tmpdy=0.0d0
        tmpdz=0.0d0
        weight=1.0d0
        if (r .gt. eps) then
          if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,radial(nradct(indspe(na)),indspe(na)), weight)
          rin=1.0d0/r
          do iq=1,nqct(indspe(na))
            qqq=dqq*iq
            qqq2=qqq*qqq
            qqqr=qqq*r
            qqqrin=rin/qqq
            sinqqqr=dsin(qqqr)
            cosqqqr=dcos(qqqr)
            besselj0=sinqqqr*qqqrin
            dbesselj0=-(sinqqqr-qqqr*cosqqqr)*qqqrin*qqqrin*qqq
            dbesselj0qq=dbesselj0*qqq2*rin
            tmppo=tmppo+coef(iq,0,indspe(na))*besselj0*qqq2
            tmpdx=tmpdx+coef(iq,0,indspe(na))*dbesselj0qq*x
            tmpdy=tmpdy+coef(iq,0,indspe(na))*dbesselj0qq*y
            tmpdz=tmpdz+coef(iq,0,indspe(na))*dbesselj0qq*z
          end do
        else
          do iq=1,nqct(indspe(na))
            qqq=dqq*iq
            tmppo=tmppo+coef(iq,0,indspe(na))*qqq*qqq
          end do
        end if
           vloc_hdp(ixyz,iapsd)=tmppo*tmp*weight
        dvlocdx_hdp(ixyz,iapsd)=tmpdx*tmp*weight
        dvlocdy_hdp(ixyz,iapsd)=tmpdy*tmp*weight
        dvlocdz_hdp(ixyz,iapsd)=tmpdz*tmp*weight
      end do
    end if
!  ------------------------------------------------------------------------

!  ----------  compute soft local part (vloc_s,dvlocd*_s)  ----------
    select case (nperi)
    case (0)
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        jx=myrx*ncpx+ix
        jy=myry*ncpy+iy
        jz=myrz*ncpz+iz
        x=(jx*dx-xmax-0.5d0*dx)-atx(na)
        y=(jy*dy-ymax-0.5d0*dy)-aty(na)
        z=(jz*dz-zmax-0.5d0*dz)-atz(na)
        r=dsqrt(x*x+y*y+z*z)
        if (r .lt. radial(nradct(indspe(na)),indspe(na))) then
          tmp=cp(7,indspe(na))*dcos(cp(6,indspe(na))*r)+cp(8,indspe(na))
          derfdrr=-cp(7,indspe(na))*cp(6,indspe(na))*dsin(cp(6,indspe(na))*r)/r
        else
          tmp=-cp1/r
          derfdrr=-tmp/(r*r)
        end if
        vloc_scw(ix,iy,iz,na)=tmp
        dvlocdx_scw(ix,iy,iz,na)=derfdrr*x
        dvlocdy_scw(ix,iy,iz,na)=derfdrr*y
        dvlocdz_scw(ix,iy,iz,na)=derfdrr*z
      end do
      end do
      end do
    case (1)
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do kx=-new_pwx,new_pwx
        if (kx*kx .ne. 0) then
          vkx=pi/xmax*kx
          vk2=vkx**2
          jy=myry*ncpy+iy
          jz=myrz*ncpz+iz
          y=aty(na)-(jy*dy-ymax-0.5d0*dy)
          z=atz(na)-(jz*dz-zmax-0.5d0*dz)
          ta=(y*y+z*z)
          tb=vkx*vkx
          anum0(kx,iy,iz)=0.0d0
          anum1(kx,iy,iz)=0.0d0
          dt=veta/nint1dmax
          do i=1,nint1dmax
            t=i*dt
            anum0(kx,iy,iz)=anum0(kx,iy,iz)+dexp(-ta*t*t-0.25d0*tb/(t*t))/t*dt
            anum1(kx,iy,iz)=anum1(kx,iy,iz)+dexp(-ta*t*t-0.25d0*tb/(t*t))*t*dt
          end do
        end if
      end do
      end do
      end do
      do kx=-new_pwx,new_pwx
        if (kx*kx .ne. 0) then
        vkx=pi/xmax*kx
        vk2=vkx**2
!$omp do
          do iz=1,ncpz
          do iy=1,ncpy
          do ix=1,ncpx
            jx=myrx*ncpx+ix
            jy=myry*ncpy+iy
            jz=myrz*ncpz+iz
            x=(jx*dx-xmax-0.5d0*dx)-atx(na)
            y=(jy*dy-ymax-0.5d0*dy)-aty(na)
            z=(jz*dz-zmax-0.5d0*dz)-atz(na)
            vcosi=dcos(vkx*x)/(xmax*2.0d0)
            vsine=dsin(vkx*x)/(xmax*2.0d0)
            vlo=-2.0d0*cp1*vcosi*anum0(kx,iy,iz)
            vloc_scw(ix,iy,iz,na)=vloc_scw(ix,iy,iz,na)+vlo
            dvlocdx_scw(ix,iy,iz,na)=dvlocdx_scw(ix,iy,iz,na)+2.0d0*cp1*vkx*vsine*anum0(kx,iy,iz)
            dvlocdy_scw(ix,iy,iz,na)=dvlocdy_scw(ix,iy,iz,na)+4.0d0*cp1*vcosi*anum1(kx,iy,iz)*y
            dvlocdz_scw(ix,iy,iz,na)=dvlocdz_scw(ix,iy,iz,na)+4.0d0*cp1*vcosi*anum1(kx,iy,iz)*z
          end do
          end do
          end do
        end if
      end do
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        jx=myrx*ncpx+ix
        jy=myry*ncpy+iy
        jz=myrz*ncpz+iz
        x=(jx*dx-xmax-0.5d0*dx)-atx(na)
        y=(jy*dy-ymax-0.5d0*dy)-aty(na)
        z=(jz*dz-zmax-0.5d0*dz)-atz(na)
        ta=(y*y+z*z)
        vlo=-cp1/xmax*(-0.5d0*expint1(ta*veta*veta,ierr)-0.5d0*dlog(ta))
        vlo0=2.0d0*(1.0d0-dexp(-veta*veta*ta))/(ta*2.0d0*xmax)
        vloc_scw(ix,iy,iz,na)=vloc_scw(ix,iy,iz,na)+vlo
        dvlocdy_scw(ix,iy,iz,na)=dvlocdy_scw(ix,iy,iz,na)+cp1*vlo0*y
        dvlocdz_scw(ix,iy,iz,na)=dvlocdz_scw(ix,iy,iz,na)+cp1*vlo0*z
      end do
      end do
      end do
    case (2)
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        jz=myrz*ncpz+iz
        z=(jz*dz-zmax-0.5d0*dz)-atz(na)
        tmp=dabs(z)
        vlo0=2.0d0*pi05*surfin*cp1*(dexp(-tmp*tmp*veta*veta)/veta+pi05*tmp*derf(veta*tmp))
        vlo1=-2.0d0*pi*surfin*cp1*derf(veta*tmp)*dsign(1.0d0,z)
        vloc_scw(ix,iy,iz,na)=vloc_scw(ix,iy,iz,na)+vlo0
        dvlocdz_scw(ix,iy,iz,na)=dvlocdz_scw(ix,iy,iz,na)-vlo1
      end do
      end do
      end do
      do ky=-new_pwy,new_pwy
      do kx=-new_pwx,new_pwx
        if (kx*kx+ky*ky .ne. 0) then
          vkx=pi/xmax*kx
          vky=pi/ymax*ky
          vk2=vkx**2+vky**2
!$omp do
          do iz=1,ncpz
          do iy=1,ncpy
          do ix=1,ncpx
            jx=myrx*ncpx+ix
            jy=myry*ncpy+iy
            jz=myrz*ncpz+iz
            x=(jx*dx-xmax-0.5d0*dx)-atx(na)
            y=(jy*dy-ymax-0.5d0*dy)-aty(na)
            z=(jz*dz-zmax-0.5d0*dz)-atz(na)
            ta=z
            tb=dsqrt(vk2)
            deftp1=(1.0d0-derf((tb-2.0d0*veta*veta*z)/(2.0d0*veta)))*dexp(-tb*z)
            deftp2=(1.0d0-derf((tb+2.0d0*veta*veta*z)/(2.0d0*veta)))*dexp( tb*z)
            vep0=pi05/tb        *0.500d0*(deftp1+deftp2)
            vep1=pi05/ta        *0.250d0*(deftp1-deftp2)
            vsine=dsin(vkx*x+vky*y)*2.0d0*surfin*pi05
            vcosi=dcos(vkx*x+vky*y)*2.0d0*surfin*pi05
            vloc_scw(ix,iy,iz,na)=vloc_scw(ix,iy,iz,na)-cp1*vcosi*vep0
            dvlocdx_scw(ix,iy,iz,na)=dvlocdx_scw(ix,iy,iz,na)+cp1*vsine*vep0*vkx
            dvlocdy_scw(ix,iy,iz,na)=dvlocdy_scw(ix,iy,iz,na)+cp1*vsine*vep0*vky
            dvlocdz_scw(ix,iy,iz,na)=dvlocdz_scw(ix,iy,iz,na)+cp1*2.0d0*vcosi*vep1*z
          end do
          end do
          end do
        end if
      end do
      end do
    case (3)
      vlo=pi*omegain*cp1/(veta*veta)
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        vloc_scw(ix,iy,iz,na)=vloc_scw(ix,iy,iz,na)+vlo
      end do
      end do
      end do
!      do kz=-new_pwz,new_pwz
      do kz=0,new_pwz
      awei=fourpi*omegain
      if (kz .ne. 0) awei=2.0d0*fourpi*omegain
      do ky=-new_pwy,new_pwy
      do kx=-new_pwx,new_pwx
        if (kx**2+ky**2+kz**2 .ne. 0) then
          vkx=pi/xmax*kx
          vky=pi/ymax*ky
          vkz=pi/zmax*kz
          vk2=vkx*vkx+vky*vky+vkz*vkz
!          tmp=-fourpi*omegain*cp1/vk2*dexp(-vk2/(4.0d0*veta*veta))
          tmp=-awei*cp1/vk2*dexp(-vk2/(4.0d0*veta*veta))
!OCL SIMD,NORECURRENCE(vloc_scw,dvlocdx_scw,dvlocdy_scw,dvlocdz_scw)
!$omp do
          do iz=1,ncpz
          do iy=1,ncpy
          do ix=1,ncpx
            jx=myrx*ncpx+ix
            jy=myry*ncpy+iy
            jz=myrz*ncpz+iz
            x=jx*dx-xmax-0.5d0*dx-atx(na)
            y=jy*dy-ymax-0.5d0*dy-aty(na)
            z=jz*dz-zmax-0.5d0*dz-atz(na)
            tmpdsin=tmp*dsin(vkx*x+vky*y+vkz*z)
            vloc_scw(ix,iy,iz,na)=vloc_scw(ix,iy,iz,na)+tmp*dcos(vkx*x+vky*y+vkz*z)
            dvlocdx_scw(ix,iy,iz,na)=dvlocdx_scw(ix,iy,iz,na)-vkx*tmpdsin
            dvlocdy_scw(ix,iy,iz,na)=dvlocdy_scw(ix,iy,iz,na)-vky*tmpdsin
            dvlocdz_scw(ix,iy,iz,na)=dvlocdz_scw(ix,iy,iz,na)-vkz*tmpdsin
          end do
          end do
          end do
        end if
      end do
      end do
      end do
    end select
    if (nperi > 0) then
      do kpz=-new_rsz,new_rsz
      do kpy=-new_rsy,new_rsy
      do kpx=-new_rsx,new_rsx
!$omp do
        do ixyz=1,ncpx*ncpy*ncpz
          iz=(ixyz-1)/(ncpx*ncpy)+1
          iy=(ixyz-(iz-1)*ncpx*ncpy-1)/ncpx+1
          ix=ixyz-(iz-1)*ncpx*ncpy-(iy-1)*ncpx
          jx=myrx*ncpx+ix
          jy=myry*ncpy+iy
          jz=myrz*ncpz+iz
          x=jx*dx-xmax-0.5d0*dx-atx(na)+kpx*xmax*2.0d0
          y=jy*dy-ymax-0.5d0*dy-aty(na)+kpy*ymax*2.0d0
          z=jz*dz-zmax-0.5d0*dz-atz(na)+kpz*zmax*2.0d0
          rrr(ix,iy,iz)=dsqrt(x*x+y*y+z*z)
        end do
!$omp do
!ocl mfunc(2)
        do ixyz=1,ncpx*ncpy*ncpz
          ddexpr(ixyz,1,1)=dexp(-(veta*rrr(ixyz,1,1))*(veta*rrr(ixyz,1,1)))
          dderfr(ixyz,1,1)=derf(veta*rrr(ixyz,1,1))
        end do
        if ((kpx**2 .le. 1) .and. (kpy**2 .le. 1) .and. (kpz**2 .le. 1)) then
!ocl simd,norecurrence(vloc_scw,dvlocdx_scw,dvlocdy_scw,dvlocdz_scw)
!$omp do
          do ixyz=1,ncpx*ncpy*ncpz
            iz=(ixyz-1)/(ncpx*ncpy)+1
            iy=(ixyz-(iz-1)*ncpx*ncpy-1)/ncpx+1
            ix=ixyz-(iz-1)*ncpx*ncpy-(iy-1)*ncpx
            jx=myrx*ncpx+ix
            jy=myry*ncpy+iy
            jz=myrz*ncpz+iz
            x=jx*dx-xmax-0.5d0*dx-atx(na)+kpx*xmax*2.0d0
            y=jy*dy-ymax-0.5d0*dy-aty(na)+kpy*ymax*2.0d0
            z=jz*dz-zmax-0.5d0*dz-atz(na)+kpz*zmax*2.0d0
            r=rrr(ixyz,1,1)
            rin=1.0d0/r
            r2in=rin*rin
            if (r .lt. radial(nradct(indspe(na)),indspe(na))) then
              tmp=cp1*rin*dderfr(ixyz,1,1) &
                      +(cp(7,indspe(na))*dcos(cp(6,indspe(na))*r)+cp(8,indspe(na)))
              derfdrr=-cp(7,indspe(na))*cp(6,indspe(na))*dsin(cp(6,indspe(na))*r)*rin &
                          +(twosqrtpiin*veta*ddexpr(ixyz,1,1)*rin-dderfr(ixyz,1,1)*r2in)*cp1*rin
            else
              tmp=cp1*rin*(dderfr(ixyz,1,1)-1.0d0)
              derfdrr=cp1*rin*(r2in+twosqrtpiin*veta*ddexpr(ixyz,1,1)*rin-dderfr(ixyz,1,1)*r2in)
            end if
            vloc_scw(ix,iy,iz,na)=vloc_scw(ix,iy,iz,na)+tmp
            dvlocdx_scw(ix,iy,iz,na)=dvlocdx_scw(ix,iy,iz,na)+derfdrr*x
            dvlocdy_scw(ix,iy,iz,na)=dvlocdy_scw(ix,iy,iz,na)+derfdrr*y
            dvlocdz_scw(ix,iy,iz,na)=dvlocdz_scw(ix,iy,iz,na)+derfdrr*z
          end do
        else
!ocl simd,norecurrence(vloc_scw,dvlocdx_scw,dvlocdy_scw,dvlocdz_scw)
!$omp do
          do ixyz=1,ncpx*ncpy*ncpz
            iz=(ixyz-1)/(ncpx*ncpy)+1
            iy=(ixyz-(iz-1)*ncpx*ncpy-1)/ncpx+1
            ix=ixyz-(iz-1)*ncpx*ncpy-(iy-1)*ncpx
            jx=myrx*ncpx+ix
            jy=myry*ncpy+iy
            jz=myrz*ncpz+iz
            x=jx*dx-xmax-0.5d0*dx-atx(na)+kpx*xmax*2.0d0
            y=jy*dy-ymax-0.5d0*dy-aty(na)+kpy*ymax*2.0d0
            z=jz*dz-zmax-0.5d0*dz-atz(na)+kpz*zmax*2.0d0
            r=rrr(ixyz,1,1)
            rin=1.0d0/r
            r2in=rin*rin
            tmp=cp1*rin*(dderfr(ixyz,1,1)-1.0d0)
            derfdrr=cp1*rin*(r2in+twosqrtpiin*veta*ddexpr(ixyz,1,1)*rin-dderfr(ixyz,1,1)*r2in)
            vloc_scw(ix,iy,iz,na)=vloc_scw(ix,iy,iz,na)+tmp
            dvlocdx_scw(ix,iy,iz,na)=dvlocdx_scw(ix,iy,iz,na)+derfdrr*x
            dvlocdy_scw(ix,iy,iz,na)=dvlocdy_scw(ix,iy,iz,na)+derfdrr*y
            dvlocdz_scw(ix,iy,iz,na)=dvlocdz_scw(ix,iy,iz,na)+derfdrr*z
          end do
        end if
      end do
      end do
      end do
    end if
!  ------------------------------------------------------------------------
  end do

! ----------  compute jellium potential  ----------
!$omp do
  do i=1,ncpz
    vjell(i)=0.0d0
  end do
  if (jelcalc==key_jel_calc) then
  do kz=-nzmax,nzmax
     if (kz**2 .ne. 0) then
     vkz=pi/zmax*kz
     vk2=vkz*vkz
!$omp do
     do iz=1,ncpz
        jz=myrz*ncpz+iz
        z=jz*dz-zmax-0.5d0*dz
        vjell(iz)=vjell(iz)-4.0d0*pi/omega*chrjel/(strjel-endjel)/vk2*(dsin(vkz*(strjel-z))-dsin(vkz*(endjel-z)))/vkz
     end do
     end if
  end do
  end if
! -------------------------------------------------
!$omp barrier
  if (nperi==1) deallocate(anum0,anum1)
  return
end subroutine pseudocalc_02


!this subroutine computes the sum of soft local part on coarse grid.
subroutine pseudocalc_03(natom,num_spe,nradmx,nf,nperi,        & ! <
                         ncpx,ncpy,ncpz,                       & ! <
                         xmax,ymax,zmax,dx,dy,dz,              & ! <
                         indspe,nradct,                        & ! <
                         vloc_scw,vjell,cp,radial,             & ! <
                         atx,aty,atz,                          & ! <
                         vtmp_cc)                                ! >
use mod_mpi,              only: myrx,myry,myrz
implicit none
integer, intent(in)::natom,num_spe,nradmx,nf,nperi
integer, intent(in)::ncpx,ncpy,ncpz
real*8, intent(in)::xmax,ymax,zmax
real*8, intent(in)::dx,dy,dz
integer, intent(in)::indspe(natom),nradct(num_spe)
real*8, intent(in)::atx(natom),aty(natom),atz(natom)
real*8, intent(in)::vloc_scw(ncpx,ncpy,ncpz,natom),vjell(ncpz),cp(8,num_spe),radial(nradmx,num_spe)
real*8, intent(out)::vtmp_cc(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf)
real*8 x,y,z,r,pi,cp1
integer ix,iy,iz,na,nf1
  pi=dacos(-1.0d0)
  nf1=nf-1

! ***NOTE*** Since vloc_scw does not have overlap region (OR), vtmp_cc at OR has
!            to be computed here in the case of isolated bc.
!$omp do
  do ix=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)
    vtmp_cc(-nf+ix,-nf+1,-nf+1)=0.0d0
  end do

  select case (nperi)
    case (0)
    do na=1,natom
      cp1=cp(1,indspe(na))
!$omp do
      do iz=-nf1,ncpz+nf
      do iy=-nf1,ncpy+nf
      do ix=-nf1,ncpx+nf
        x=atx(na)-((myrx*ncpx+ix)*dx-xmax-0.5d0*dx)
        y=aty(na)-((myry*ncpy+iy)*dy-ymax-0.5d0*dy)
        z=atz(na)-((myrz*ncpz+iz)*dz-zmax-0.5d0*dz)
        r=dsqrt(x*x+y*y+z*z)
        if (r .lt. radial(nradct(indspe(na)),indspe(na))) then
          vtmp_cc(ix,iy,iz)=vtmp_cc(ix,iy,iz)+(cp(7,indspe(na))*dcos(cp(6,indspe(na))*r)+cp(8,indspe(na)))
        else
          vtmp_cc(ix,iy,iz)=vtmp_cc(ix,iy,iz)-cp1/r
        end if
      end do
      end do
      end do
!$omp end do nowait
    end do
!$omp barrier
  case (1:3)
    do na=1,natom
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        vtmp_cc(ix,iy,iz)=vtmp_cc(ix,iy,iz)+vloc_scw(ix,iy,iz,na)
      end do
      end do
      end do
!$omp end do nowait
    end do
!$omp barrier
  end select
  if (nperi==3) then
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
      vtmp_cc(ix,iy,iz)=vtmp_cc(ix,iy,iz)+vjell(iz)
    end do
    end do
    end do
  end if

  return
end subroutine pseudocalc_03


!this subroutine computes the sum of local parts on dense grid.
subroutine pseudocalc_04(natom,num_list_d,num_ppcell_d,ncpx_d,ncpy_d,ncpz_d, & ! <
                         key_natpri_inps,                                    & ! <
                         natprid,napsd,natinfd_vloc,lstvecd2,                & ! <
                         vloc_hdp,                                           & ! <
                         vloc_dense)                                           ! X
implicit none
integer, intent(in)::natom,num_list_d,num_ppcell_d,ncpx_d,ncpy_d,ncpz_d
integer, intent(in)::key_natpri_inps
integer, intent(in)::natprid(natom),napsd(natom),natinfd_vloc(natom),lstvecd2(num_list_d,num_ppcell_d)
real*8, intent(in)::vloc_hdp(num_list_d,num_ppcell_d)
real*8, intent(inout)::vloc_dense(ncpx_d,ncpy_d,ncpz_d)
integer i,iapsd,na

  do na=1,natom
    if (natprid(na) .eq. key_natpri_inps) then
      iapsd=napsd(na)
!$omp do
      do i=1,natinfd_vloc(na)
        vloc_dense(lstvecd2(i,iapsd),1,1)=vloc_dense(lstvecd2(i,iapsd),1,1)+vloc_hdp(i,iapsd)
      end do
    end if
  end do

  return
end subroutine pseudocalc_04


!this subroutine computes the boundary condition for the isolated system.
subroutine pseudocalc_05(nperi,natom,num_spe,nf,nmesh,new_pwx,new_pwy,new_rsx,new_rsy, & ! <
                         ncpx_d,ncpy_d,ncpz_d,nint1dmax,                               & ! <
                         ddx,ddy,ddz,xmax,ymax,zmax,veta,                              & ! <
                         indspe,                                                       & ! <
                         atx,aty,atz,                                                  & ! <
                         cp,                                                           & ! <
                         vtmp_dd)                                                        ! X
use mod_mpi,              only: myrx,myry,myrz,nprocx,nprocy,nprocz
use mod_mathfunctions,only: expint1
implicit none
integer,intent(in)::nperi,natom,num_spe,nf,nmesh,new_pwx,new_pwy,new_rsx,new_rsy
integer,intent(in)::ncpx_d,ncpy_d,ncpz_d
integer,intent(in)::nint1dmax
integer,intent(in)::indspe(natom)
real*8, intent(in)::ddx,ddy,ddz
real*8, intent(in)::xmax,ymax,zmax
real*8, intent(in)::veta
real*8, intent(in)::atx(natom),aty(natom),atz(natom)
real*8, intent(in)::cp(8,num_spe)
real*8, intent(inout)::vtmp_dd(-(nf*nmesh-1):ncpx_d+nf*nmesh,-(nf*nmesh-1):ncpy_d+nf*nmesh,-(nf*nmesh-1):ncpz_d+nf*nmesh)
real*8, allocatable::anumy(:,:,:),anumz(:,:,:)
real*8 x,y,z,r,rin,cp1,pi,ta,tb,vkx,vky,vk2,tmp,vlo,vlo0,pi05,surf,surfin,deftp1,deftp2,vep0,vep1,vcosi,t,dt
real*8 derf
integer na,ix,iy,iz,jx,jy,jz,kx,ky,kpx,kpy,nf2,nf3,i,ierr

  nf2=nf*nmesh
  nf3=nf2-1
  if (nperi==1) allocate(anumy(-new_pwx:new_pwx,-nf3:nf2,-nf3:ncpz_d+nf2),anumz(-new_pwx:new_pwx,-nf3:ncpy_d+nf2,-nf3:nf2))
  pi=dacos(-1.0d0)
  pi05=dsqrt(pi)
  surf=xmax*ymax*4.0d0
  surfin=1.0d0/surf

  select case (nperi)
  case (0)
    do na=1,natom
      cp1=cp(1,indspe(na))
      if (myrx .eq. 0) then
!$omp do
        do iz=-nf3,ncpz_d+nf2
        do iy=-nf3,ncpy_d+nf2
        do ix=-nf3,0
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=jx*ddx-xmax-0.5d0*ddx-atx(na)
          y=jy*ddy-ymax-0.5d0*ddy-aty(na)
          z=jz*ddz-zmax-0.5d0*ddz-atz(na)
          r=dsqrt(x*x+y*y+z*z)
          rin=1.0d0/r
          vtmp_dd(ix,iy,iz)=vtmp_dd(ix,iy,iz)-cp1*rin
        end do
        end do
        end do
      end if
      if (myrx .eq. nprocx-1) then
!$omp do
        do iz=-nf3,ncpz_d+nf2
        do iy=-nf3,ncpy_d+nf2
        do ix=1,nf2
          jx=myrx*ncpx_d+ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=jx*ddx-xmax-0.5d0*ddx-atx(na)
          y=jy*ddy-ymax-0.5d0*ddy-aty(na)
          z=jz*ddz-zmax-0.5d0*ddz-atz(na)
          r=dsqrt(x*x+y*y+z*z)
          rin=1.0d0/r
          vtmp_dd(ix+ncpx_d,iy,iz)=vtmp_dd(ix+ncpx_d,iy,iz)-cp1*rin
        end do
        end do
        end do
      end if
      if (myry .eq. 0) then
!$omp do
        do iz=-nf3,ncpz_d+nf2
        do iy=-nf3,0
        do ix=-nf3,ncpx_d+nf2
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=jx*ddx-xmax-0.5d0*ddx-atx(na)
          y=jy*ddy-ymax-0.5d0*ddy-aty(na)
          z=jz*ddz-zmax-0.5d0*ddz-atz(na)
          r=dsqrt(x*x+y*y+z*z)
          rin=1.0d0/r
          vtmp_dd(ix,iy,iz)=vtmp_dd(ix,iy,iz)-cp1*rin
        end do
        end do
        end do
      end if
      if (myry .eq. nprocy-1) then
!$omp do
        do iz=-nf3,ncpz_d+nf2
        do iy=1,nf2
        do ix=-nf3,ncpx_d+nf2
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=jx*ddx-xmax-0.5d0*ddx-atx(na)
          y=jy*ddy-ymax-0.5d0*ddy-aty(na)
          z=jz*ddz-zmax-0.5d0*ddz-atz(na)
          r=dsqrt(x*x+y*y+z*z)
          rin=1.0d0/r
          vtmp_dd(ix,iy+ncpy_d,iz)=vtmp_dd(ix,iy+ncpy_d,iz)-cp1*rin
        end do
        end do
        end do
      end if
      if (myrz .eq. 0) then
!$omp do
        do iz=-nf3,0
        do iy=-nf3,ncpy_d+nf2
        do ix=-nf3,ncpx_d+nf2
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=jx*ddx-xmax-0.5d0*ddx-atx(na)
          y=jy*ddy-ymax-0.5d0*ddy-aty(na)
          z=jz*ddz-zmax-0.5d0*ddz-atz(na)
          r=dsqrt(x*x+y*y+z*z)
          rin=1.0d0/r
          vtmp_dd(ix,iy,iz)=vtmp_dd(ix,iy,iz)-cp1*rin
        end do
        end do
        end do
      end if
      if (myrz .eq. nprocz-1) then
!$omp do
        do iz=1,nf2
        do iy=-nf3,ncpy_d+nf2
        do ix=-nf3,ncpx_d+nf2
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+ncpz_d+iz
          x=jx*ddx-xmax-0.5d0*ddx-atx(na)
          y=jy*ddy-ymax-0.5d0*ddy-aty(na)
          z=jz*ddz-zmax-0.5d0*ddz-atz(na)
          r=dsqrt(x*x+y*y+z*z)
          rin=1.0d0/r
          vtmp_dd(ix,iy,iz+ncpz_d)=vtmp_dd(ix,iy,iz+ncpz_d)-cp1*rin
        end do
        end do
        end do
      end if
    end do
  case (1)
    do na=1,natom
      cp1=cp(1,indspe(na))
      if (myry .eq. 0) then
!$omp do
        do iz=-nf3,ncpz_d+nf2
        do iy=-nf3,0
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
            vkx=pi/xmax*kx
            vk2=vkx**2
            jy=myry*ncpy_d+iy
            jz=myrz*ncpz_d+iz
            y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
            z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
            ta=(y*y+z*z)
            tb=vkx*vkx
            anumy(kx,iy,iz)=0.0d0
            dt=veta/nint1dmax
            do i=1,nint1dmax
              t=i*dt
              anumy(kx,iy,iz)=anumy(kx,iy,iz)+dexp(-ta*t*t-0.25d0*tb/(t*t))/t*dt
            end do
          end if
        end do
        end do
        end do
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
          vkx=pi/xmax*kx
          vk2=vkx**2
!$omp do
            do iz=-nf3,ncpz_d+nf2
            do iy=-nf3,0
            do ix=-nf3,ncpx_d+nf2
              jx=myrx*ncpx_d+ix
              jy=myry*ncpy_d+iy
              jz=myrz*ncpz_d+iz
              x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
              y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
              z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
              vcosi=dcos(vkx*x)/(xmax*2.0d0)
              vlo=-2.0d0*cp1*vcosi*anumy(kx,iy,iz)
              vtmp_dd(ix,iy,iz)=vtmp_dd(ix,iy,iz)+vlo
            end do
            end do
            end do
          end if
        end do
!$omp do
        do iz=-nf3,ncpz_d+nf2
        do iy=-nf3,0
        do ix=-nf3,ncpx_d+nf2
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
          y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
          z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
          ta=(y*y+z*z)
          vlo=-cp1/xmax*(-0.5d0*expint1(ta*veta*veta,ierr)-0.5d0*dlog(ta))
          vlo0=2.0d0*(1.0d0-dexp(-veta*veta*ta))/(ta*2.0d0*xmax)
          vtmp_dd(ix,iy,iz)=vtmp_dd(ix,iy,iz)+vlo
        end do
        end do
        end do
        do kpx=-new_rsx,new_rsx
!ocl simd,norecurrence(vtmp_dd)
!$omp do
          do iz=-nf3,ncpz_d+nf2
          do iy=-nf3,0
          do ix=-nf3,ncpx_d+nf2
            jx=myrx*ncpx_d+ix
            jy=myry*ncpy_d+iy
            jz=myrz*ncpz_d+iz
            x=jx*ddx-xmax-0.5d0*ddx-atx(na)+kpx*xmax*2.0d0
            y=jy*ddy-ymax-0.5d0*ddy-aty(na)
            z=jz*ddz-zmax-0.5d0*ddz-atz(na)
            r=dsqrt(x*x+y*y+z*z)
            rin=1.0d0/r
            tmp=cp1*rin*(derf(veta*r)-1.0d0)
            vtmp_dd(ix,iy,iz)=vtmp_dd(ix,iy,iz)+tmp
          end do
          end do
          end do
        end do
      end if
      if (myry .eq. nprocy-1) then
!$omp do
        do iz=-nf3,ncpz_d+nf2
        do iy=1,nf2
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
            vkx=pi/xmax*kx
            vk2=vkx**2
            jy=myry*ncpy_d+ncpy_d+iy
            jz=myrz*ncpz_d+iz
            y=aty(na)-(jy*ddy-ymax-0.5d0*ddy)
            z=atz(na)-(jz*ddz-zmax-0.5d0*ddz)
            ta=(y*y+z*z)
            tb=vkx*vkx
            anumy(kx,iy,iz)=0.0d0
            dt=veta/nint1dmax
            do i=1,nint1dmax
              t=i*dt
              anumy(kx,iy,iz)=anumy(kx,iy,iz)+dexp(-ta*t*t-0.25d0*tb/(t*t))/t*dt
            end do
          end if
        end do
        end do
        end do
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
          vkx=pi/xmax*kx
          vk2=vkx**2
!$omp do
            do iz=-nf3,ncpz_d+nf2
            do iy=1,nf2
            do ix=-nf3,ncpx_d+nf2
              jx=myrx*ncpx_d+ix
              jy=myry*ncpy_d+ncpy_d+iy
              jz=myrz*ncpz_d+iz
              x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
              y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
              z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
              vcosi=dcos(vkx*x)/(xmax*2.0d0)
              vlo=-2.0d0*cp1*vcosi*anumy(kx,iy,iz)
              vtmp_dd(ix,iy+ncpy_d,iz)=vtmp_dd(ix,iy+ncpy_d,iz)+vlo
            end do
            end do
            end do
          end if
        end do
!$omp do
        do iz=-nf3,ncpz_d+nf2
        do iy=1,nf2
        do ix=-nf3,ncpx_d+nf2
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
          y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
          z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
          ta=(y*y+z*z)
          vlo=-cp1/xmax*(-0.5d0*expint1(ta*veta*veta,ierr)-0.5d0*dlog(ta))
          vlo0=2.0d0*(1.0d0-dexp(-veta*veta*ta))/(ta*2.0d0*xmax)
          vtmp_dd(ix,iy+ncpy_d,iz)=vtmp_dd(ix,iy+ncpy_d,iz)+vlo
        end do
        end do
        end do
        do kpx=-new_rsx,new_rsx
!ocl simd,norecurrence(vtmp_dd)
!$omp do
          do iz=-nf3,ncpz_d+nf2
          do iy=1,nf2
          do ix=-nf3,ncpx_d+nf2
            jx=myrx*ncpx_d+ix
            jy=myry*ncpy_d+ncpy_d+iy
            jz=myrz*ncpz_d+iz
            x=jx*ddx-xmax-0.5d0*ddx-atx(na)+kpx*xmax*2.0d0
            y=jy*ddy-ymax-0.5d0*ddy-aty(na)
            z=jz*ddz-zmax-0.5d0*ddz-atz(na)
            r=dsqrt(x*x+y*y+z*z)
            rin=1.0d0/r
            tmp=cp1*rin*(derf(veta*r)-1.0d0)
            vtmp_dd(ix,iy+ncpy_d,iz)=vtmp_dd(ix,iy+ncpy_d,iz)+tmp
          end do
          end do
          end do
        end do
      end if
      if (myrz .eq. 0) then
!$omp do
        do iz=-nf3,0
        do iy=-nf3,ncpy_d+nf2
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
            vkx=pi/xmax*kx
            vk2=vkx**2
            jy=myry*ncpy_d+iy
            jz=myrz*ncpz_d+iz
            y=aty(na)-(jy*ddy-ymax-0.5d0*ddy)
            z=atz(na)-(jz*ddz-zmax-0.5d0*ddz)
            ta=(y*y+z*z)
            tb=vkx*vkx
            anumz(kx,iy,iz)=0.0d0
            dt=veta/nint1dmax
            do i=1,nint1dmax
              t=i*dt
              anumz(kx,iy,iz)=anumz(kx,iy,iz)+dexp(-ta*t*t-0.25d0*tb/(t*t))/t*dt
            end do
          end if
        end do
        end do
        end do
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
          vkx=pi/xmax*kx
          vk2=vkx**2
!$omp do
            do iz=-nf3,0
            do iy=-nf3,ncpy_d+nf2
            do ix=-nf3,ncpx_d+nf2
              jx=myrx*ncpx_d+ix
              jy=myry*ncpy_d+iy
              jz=myrz*ncpz_d+iz
              x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
              y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
              z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
              vcosi=dcos(vkx*x)/(xmax*2.0d0)
              vlo=-2.0d0*cp1*vcosi*anumz(kx,iy,iz)
              vtmp_dd(ix,iy,iz)=vtmp_dd(ix,iy,iz)+vlo
            end do
            end do
            end do
          end if
        end do
!$omp do
        do iz=-nf3,0
        do iy=-nf3,ncpy_d+nf2
        do ix=-nf3,ncpx_d+nf2
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
          y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
          z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
          ta=(y*y+z*z)
          vlo=-cp1/xmax*(-0.5d0*expint1(ta*veta*veta,ierr)-0.5d0*dlog(ta))
          vlo0=2.0d0*(1.0d0-dexp(-veta*veta*ta))/(ta*2.0d0*xmax)
          vtmp_dd(ix,iy,iz)=vtmp_dd(ix,iy,iz)+vlo
        end do
        end do
        end do
        do kpx=-new_rsx,new_rsx
!ocl simd,norecurrence(vtmp_dd)
!$omp do
          do iz=-nf3,0
          do iy=-nf3,ncpy_d+nf2
          do ix=-nf3,ncpx_d+nf2
            jx=myrx*ncpx_d+ix
            jy=myry*ncpy_d+iy
            jz=myrz*ncpz_d+iz
            x=jx*ddx-xmax-0.5d0*ddx-atx(na)+kpx*xmax*2.0d0
            y=jy*ddy-ymax-0.5d0*ddy-aty(na)
            z=jz*ddz-zmax-0.5d0*ddz-atz(na)
            r=dsqrt(x*x+y*y+z*z)
            rin=1.0d0/r
            tmp=cp1*rin*(derf(veta*r)-1.0d0)
            vtmp_dd(ix,iy,iz)=vtmp_dd(ix,iy,iz)+tmp
          end do
          end do
          end do
        end do
      end if
      if (myrz .eq. nprocz-1) then
!$omp do
        do iz=1,nf2
        do iy=-nf3,ncpy_d+nf2
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
            vkx=pi/xmax*kx
            vk2=vkx**2
            jy=myry*ncpy_d+iy
            jz=myrz*ncpz_d+ncpz_d+iz
            y=aty(na)-(jy*ddy-ymax-0.5d0*ddy)
            z=atz(na)-(jz*ddz-zmax-0.5d0*ddz)
            ta=(y*y+z*z)
            tb=vkx*vkx
            anumz(kx,iy,iz)=0.0d0
            dt=veta/nint1dmax
            do i=1,nint1dmax
              t=i*dt
              anumz(kx,iy,iz)=anumz(kx,iy,iz)+dexp(-ta*t*t-0.25d0*tb/(t*t))/t*dt
            end do
          end if
        end do
        end do
        end do
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
          vkx=pi/xmax*kx
          vk2=vkx**2
!$omp do
            do iz=1,nf2
            do iy=-nf3,ncpy_d+nf2
            do ix=-nf3,ncpx_d+nf2
              jx=myrx*ncpx_d+ix
              jy=myry*ncpy_d+iy
              jz=myrz*ncpz_d+ncpz_d+iz
              x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
              y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
              z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
              vcosi=dcos(vkx*x)/(xmax*2.0d0)
              vlo=-2.0d0*cp1*vcosi*anumz(kx,iy,iz)
              vtmp_dd(ix,iy,iz+ncpz_d)=vtmp_dd(ix,iy,iz+ncpz_d)+vlo
            end do
            end do
            end do
          end if
        end do
!$omp do
        do iz=1,nf2
        do iy=-nf3,ncpy_d+nf2
        do ix=-nf3,ncpx_d+nf2
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+ncpz_d+iz
          x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
          y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
          z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
          ta=(y*y+z*z)
          vlo=-cp1/xmax*(-0.5d0*expint1(ta*veta*veta,ierr)-0.5d0*dlog(ta))
          vlo0=2.0d0*(1.0d0-dexp(-veta*veta*ta))/(ta*2.0d0*xmax)
          vtmp_dd(ix,iy,iz+ncpz_d)=vtmp_dd(ix,iy,iz+ncpz_d)+vlo
        end do
        end do
        end do
        do kpx=-new_rsx,new_rsx
!ocl simd,norecurrence(vtmp_dd)
!$omp do
          do iz=1,nf2
          do iy=-nf3,ncpy_d+nf2
          do ix=-nf3,ncpx_d+nf2
            jx=myrx*ncpx_d+ix
            jy=myry*ncpy_d+iy
            jz=myrz*ncpz_d+ncpz_d+iz
            x=jx*ddx-xmax-0.5d0*ddx-atx(na)+kpx*xmax*2.0d0
            y=jy*ddy-ymax-0.5d0*ddy-aty(na)
            z=jz*ddz-zmax-0.5d0*ddz-atz(na)
            r=dsqrt(x*x+y*y+z*z)
            rin=1.0d0/r
            tmp=cp1*rin*(derf(veta*r)-1.0d0)
            vtmp_dd(ix,iy,iz+ncpz_d)=vtmp_dd(ix,iy,iz+ncpz_d)+tmp
          end do
          end do
          end do
        end do
      end if
    end do
  case (2)
    do na=1,natom
      cp1=cp(1,indspe(na))
      if (myrz .eq. 0) then
!$omp do
        do iz=-nf3,0
        do iy=-nf3,ncpy_d+nf2
        do ix=-nf3,ncpx_d+nf2
          jz=myrz*ncpz_d+iz
          z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
          tmp=dabs(z)
          vlo0=2.0d0*pi05*surfin*cp1*(dexp(-tmp*tmp*veta*veta)/veta+pi05*tmp*derf(tmp*veta))
          vtmp_dd(ix,iy,iz)=vtmp_dd(ix,iy,iz)+vlo0
        end do
        end do
        end do
        do ky=-new_pwy,new_pwy
        do kx=-new_pwx,new_pwx
          if (kx*kx+ky*ky .ne. 0) then
            vkx=pi/xmax*kx
            vky=pi/ymax*ky
            vk2=vkx**2+vky**2
!$omp do
            do iz=-nf3,0
            do iy=-nf3,ncpy_d+nf2
            do ix=-nf3,ncpx_d+nf2
              jx=myrx*ncpx_d+ix
              jy=myry*ncpy_d+iy
              jz=myrz*ncpz_d+iz
              x=jx*ddx-xmax-0.5d0*ddx-atx(na)
              y=jy*ddy-ymax-0.5d0*ddy-aty(na)
              z=jz*ddz-zmax-0.5d0*ddz-atz(na)
              ta=z
              tb=dsqrt(vk2)
              deftp1=(1.0d0-derf((tb-2.0d0*veta*veta*z)/(2.0d0*veta)))*dexp(-tb*z)
              deftp2=(1.0d0-derf((tb+2.0d0*veta*veta*z)/(2.0d0*veta)))*dexp( tb*z)
              vep0=pi05/tb        *0.500d0*(deftp1+deftp2)
              vep1=pi05/ta        *0.250d0*(deftp1-deftp2)
              vcosi=dcos(vkx*x+vky*y)*2.0d0*surfin*pi05
              vtmp_dd(ix,iy,iz)=vtmp_dd(ix,iy,iz)-cp1*vcosi*vep0
            end do
            end do
            end do
          end if
        end do
        end do
        do kpy=-new_rsy,new_rsy
        do kpx=-new_rsx,new_rsx
!ocl simd,norecurrence(vtmp_dd)
!$omp do
          do iz=-nf3,0
          do iy=-nf3,ncpy_d+nf2
          do ix=-nf3,ncpx_d+nf2
            jx=myrx*ncpx_d+ix
            jy=myry*ncpy_d+iy
            jz=myrz*ncpz_d+iz
            x=jx*ddx-xmax-0.5d0*ddx-atx(na)+kpx*xmax*2.0d0
            y=jy*ddy-ymax-0.5d0*ddy-aty(na)+kpy*ymax*2.0d0
            z=jz*ddz-zmax-0.5d0*ddz-atz(na)
            r=dsqrt(x*x+y*y+z*z)
            rin=1.0d0/r
            tmp=cp1*rin*(derf(veta*r)-1.0d0)
            vtmp_dd(ix,iy,iz)=vtmp_dd(ix,iy,iz)+tmp
          end do
          end do
          end do
        end do
        end do
!$omp barrier
      end if
      if (myrz .eq. nprocz-1) then
!$omp do
        do iz=1,nf2
        do iy=-nf3,ncpy_d+nf2
        do ix=-nf3,ncpx_d+nf2
          jz=myrz*ncpz_d+ncpz_d+iz
          z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
          tmp=dabs(z)
          vlo0=2.0d0*pi05*surfin*cp1*(dexp(-tmp*tmp*veta*veta)/veta+pi05*tmp*derf(tmp*veta))
          vtmp_dd(ix,iy,iz+ncpz_d)=vtmp_dd(ix,iy,iz+ncpz_d)+vlo0
        end do
        end do
        end do
        do ky=-new_pwy,new_pwy
        do kx=-new_pwx,new_pwx
          if (kx*kx+ky*ky .ne. 0) then
            vkx=pi/xmax*kx
            vky=pi/ymax*ky
            vk2=vkx**2+vky**2
!$omp do
            do iz=1,nf2
            do iy=-nf3,ncpy_d+nf2
            do ix=-nf3,ncpx_d+nf2
              jx=myrx*ncpx_d+ix
              jy=myry*ncpy_d+iy
              jz=myrz*ncpz_d+ncpz_d+iz
              x=jx*ddx-xmax-0.5d0*ddx-atx(na)
              y=jy*ddy-ymax-0.5d0*ddy-aty(na)
              z=jz*ddz-zmax-0.5d0*ddz-atz(na)
              ta=z
              tb=dsqrt(vk2)
              deftp1=(1.0d0-derf((tb-2.0d0*veta*veta*z)/(2.0d0*veta)))*dexp(-tb*z)
              deftp2=(1.0d0-derf((tb+2.0d0*veta*veta*z)/(2.0d0*veta)))*dexp( tb*z)
              vep0=pi05/tb        *0.500d0*(deftp1+deftp2)
              vep1=pi05/ta        *0.250d0*(deftp1-deftp2)
              vcosi=dcos(vkx*x+vky*y)*2.0d0*surfin*pi05
              vtmp_dd(ix,iy,iz+ncpz_d)=vtmp_dd(ix,iy,iz+ncpz_d)-cp1*vcosi*vep0
            end do
            end do
            end do
          end if
        end do
        end do
        do kpy=-new_rsy,new_rsy
        do kpx=-new_rsx,new_rsx
!ocl simd,norecurrence(vtmp_dd)
!$omp do
          do iz=1,nf2
          do iy=-nf3,ncpy_d+nf2
          do ix=-nf3,ncpx_d+nf2
            jx=myrx*ncpx_d+ix
            jy=myry*ncpy_d+iy
            jz=myrz*ncpz_d+ncpz_d+iz
            x=jx*ddx-xmax-0.5d0*ddx-atx(na)+kpx*xmax*2.0d0
            y=jy*ddy-ymax-0.5d0*ddy-aty(na)+kpy*ymax*2.0d0
            z=jz*ddz-zmax-0.5d0*ddz-atz(na)
            r=dsqrt(x*x+y*y+z*z)
            rin=1.0d0/r
            tmp=cp1*rin*(derf(veta*r)-1.0d0)
            vtmp_dd(ix,iy,iz+ncpz_d)=vtmp_dd(ix,iy,iz+ncpz_d)+tmp
          end do
          end do
          end do
        end do
        end do
!$omp barrier
      end if
    end do
  end select

  if (nperi==1) deallocate(anumy,anumz)
  return
end subroutine pseudocalc_05


!this subroutine computes pcc charge.
subroutine pseudocalc_06(natom,nradmx,num_spe,npoint,num_atcell,nfiltyp,nqmx, & ! <
                         ncpx_d,ncpy_d,ncpz_d,                                & ! <
                         key_natpri_in,                                       & ! <
                         eps,psftrad,psctoff,filpp,rctpcc,                    & ! <
                         xmax,ymax,zmax,dx,dy,dz,ddx,ddy,ddz,                 & ! <
                         nradct,indspe,nqctpcc,natpri,natpri_inf,             & ! <
                         point,radial,coef,                                   & ! <
                         atx,aty,atz,                                         & ! <
                         rhopcc_dense,rhopccr)                                  ! >
use mod_mathfunctions, only: fermidis
use mod_mpi,           only: myrx,myry,myrz
implicit none
integer, intent(in)::natom,nradmx,num_spe,npoint,num_atcell,nfiltyp,nqmx
integer, intent(in)::key_natpri_in
integer, intent(in)::ncpx_d,ncpy_d,ncpz_d
real*8, intent(in)::eps,psftrad,psctoff,filpp,rctpcc
real*8, intent(in)::xmax,ymax,zmax
real*8, intent(in)::dx,dy,dz,ddx,ddy,ddz
integer, intent(in)::nradct(num_spe),indspe(natom),nqctpcc(num_spe),natpri(natom),natpri_inf(natom)
real*8, intent(in)::point(npoint,3),radial(nradmx,num_spe),coef(0:nqmx,0:7,num_spe)
real*8, intent(in)::atx(natom),aty(natom),atz(natom)
real*8, intent(out)::rhopcc_dense(ncpx_d,ncpy_d,ncpz_d),rhopccr(nradmx,npoint,num_atcell)
real*8 pi,twopicbin,fourpi,rcutss,dqq,tmp,x,y,z,r,rin,qqq,qqqr,qqqrin,tmpx,tmpy,tmpz,tmppcc,besselj0,weight,psctoff2,rcut
integer ix,iy,iz,ipri,il,na,kxmx,kymx,kzmx,kx,ky,kz,jx,jy,jz,ir,iq,na0,na1

  pi=dacos(-1.0d0)
  twopicbin=1.0d0/(2.0d0*pi)**3
  fourpi=4.0d0*pi

  psctoff2= psctoff
  if (nfiltyp==1) psctoff2= psctoff2*rctpcc

!  ----------  zero clear  ----------
!$omp do
  do ix=1,ncpx_d*ncpy_d*ncpz_d
    rhopcc_dense(ix,1,1)=0.0d0
  end do
!$omp do
  do ir=1,nradmx*npoint*num_atcell
    rhopccr(ir,1,1)=0.0d0
  end do
!  ----------------------------------

  do na=1,natom
    tmppcc=0.0d0
    do iq=1,nqctpcc(indspe(na))
      tmppcc=tmppcc+coef(iq,7,indspe(na))*coef(iq,7,indspe(na))
    end do
    if (tmppcc .gt. 1.0d-15) then
      rcutss=radial(nradct(indspe(na)),indspe(na))*psftrad
      dqq=pi/rcutss
      tmp=fourpi*dqq*twopicbin
      rcut= radial(nradct(indspe(na)),indspe(na))*psctoff2
      kxmx= idint(rcut/(2.0d0*xmax)+1.0d0+(dx+ddx)/4.0d0)
      kymx= idint(rcut/(2.0d0*ymax)+1.0d0+(dy+ddy)/4.0d0)
      kzmx= idint(rcut/(2.0d0*zmax)+1.0d0+(dz+ddz)/4.0d0)
      do kz=-kzmx,kzmx
      do ky=-kymx,kymx
      do kx=-kxmx,kxmx
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
          if (r .lt. rcut) then
            weight=1.0d0
            if (r .gt. eps) then
              if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,rctpcc*radial(nradct(indspe(na)),indspe(na)), weight)
              rin=1.0d0/r
              tmppcc=0.0d0
              do iq=1,nqctpcc(indspe(na))
                qqq=dqq*iq
                qqqrin=rin/qqq
                besselj0=dsin(qqq*r)*qqqrin
                tmppcc=tmppcc+coef(iq,7,indspe(na))*besselj0*qqq*qqq
              end do
            else
              tmppcc=0.0d0
              do iq=1,nqctpcc(indspe(na))
                qqq=dqq*iq
                tmppcc=tmppcc+coef(iq,7,indspe(na))*qqq*qqq
              end do
            end if
            rhopcc_dense(ix,iy,iz)=rhopcc_dense(ix,iy,iz)+tmppcc*tmp*weight
          end if
        end do
        end do
        end do
!$omp end do nowait
      end do
      end do
      end do
    end if
  end do
!$omp barrier

  do na0=1,natom
    tmppcc=0.0d0
    do iq=1,nqctpcc(indspe(na0))
      tmppcc=tmppcc+coef(iq,7,indspe(na0))*coef(iq,7,indspe(na0))
    end do
    if (tmppcc .gt. 1.0d-15) then
      if (natpri(na0) .eq. key_natpri_in) then
        ipri=natpri_inf(na0)
        do na1=1,natom
          rcutss=radial(nradct(indspe(na1)),indspe(na1))*psftrad
          dqq=pi/rcutss
          tmp=fourpi*dqq*twopicbin
          rcut= 2.0d0*radial(nradct(indspe(na1)),indspe(na1))*psctoff2
          kxmx= idint(rcut/(2.0d0*xmax)+1.0d0)
          kymx= idint(rcut/(2.0d0*ymax)+1.0d0)
          kzmx= idint(rcut/(2.0d0*zmax)+1.0d0)
          do kz=-kzmx,kzmx
          do ky=-kymx,kymx
          do kx=-kxmx,kxmx
            tmpx=atx(na0)-atx(na1)+kx*xmax*2.0d0
            tmpy=aty(na0)-aty(na1)+ky*ymax*2.0d0
            tmpz=atz(na0)-atz(na1)+kz*zmax*2.0d0
            r=dsqrt(tmpx*tmpx+tmpy*tmpy+tmpz*tmpz)
            if (r .lt. rcut) then
!$omp do
            do il=1,npoint
            do ir=2,nradct(indspe(na0))
              x=point(il,1)*radial(ir,indspe(na0))+tmpx
              y=point(il,2)*radial(ir,indspe(na0))+tmpy
              z=point(il,3)*radial(ir,indspe(na0))+tmpz
              r=dsqrt(x*x+y*y+z*z)
              if (r .lt. radial(nradct(indspe(na1)),indspe(na1))*psctoff2) then
                weight=1.0d0
                if (r .gt. eps) then
                  if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,rctpcc*radial(nradct(indspe(na1)),indspe(na1)), weight)
                  rin=1.0d0/r
                  tmppcc=0.0d0
                  do iq=1,nqctpcc(indspe(na1))
                    qqq=dqq*iq
                    qqqr=qqq*r
                    qqqrin=rin/qqq
                    besselj0=dsin(qqqr)*qqqrin
                    tmppcc=tmppcc+coef(iq,7,indspe(na1))*besselj0*qqq*qqq
                  end do
                else
                  tmppcc=0.0d0
                  do iq=1,nqctpcc(indspe(na1))
                    qqq=dqq*iq
                    tmppcc=tmppcc+coef(iq,7,indspe(na1))*qqq*qqq
                  end do
                end if
                rhopccr(ir,il,ipri)=rhopccr(ir,il,ipri)+tmppcc*tmp*weight
              end if
            end do
            end do
!$omp end do nowait
            end if
          end do
          end do
          end do
        end do
      end if
    end if
  end do
!$omp barrier

  return
end subroutine pseudocalc_06


!this subroutine computes non-local parts of pseudopotential.
subroutine pseudocalc_07(natom,num_spe,npmesh,nprjmx,nfdg,nfiltyp,num_list,num_ppcell, & ! <
                         nradmx,nqmx,                                                  & ! <
                         nxmax,nymax,nzmax,npxmax,npymax,npzmax,nxspd,nyspd,nzspd,     & ! <
                         key_natpri_in,key_natpri_inps,                                & ! <
                         eps,psftrad,psctoff,filpp,                                    & ! <
                         dx,dy,dz,                                                     & ! <
                         indspe,nprj,nradct,natpri,naps,natinf,nqct,                   & ! <
                         natx,naty,natz,lstx,lsty,lstz,                                & ! <
                         wxyz,radial,coef,                                             & ! <
                         atx,aty,atz,                                                  & ! <
                         vnlocp,                                                       & ! >
                         vtmp_pp,vtmp_pp1,vtmp_pp2,vnlspd)                               ! X
use mod_mathfunctions, only: fermidis
implicit none
integer, intent(in)::natom,num_spe,npmesh,nprjmx,nfdg,nfiltyp,num_list,num_ppcell,nradmx,nqmx
integer, intent(in)::nxmax,nymax,nzmax,npxmax,npymax,npzmax
integer, intent(in)::nxspd,nyspd,nzspd
integer, intent(in)::key_natpri_in,key_natpri_inps
real*8, intent(in)::eps,psftrad,psctoff,filpp
real*8, intent(in)::dx,dy,dz
integer, intent(in)::indspe(natom),nprj(num_spe),nradct(num_spe),natpri(natom),naps(natom),natinf(natom),nqct(num_spe)
integer, intent(in)::natx(natom),naty(natom),natz(natom)
integer, intent(in)::lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell)
real*8, intent(in)::wxyz(-nfdg*npmesh:nfdg*npmesh-1),radial(nradmx,num_spe),coef(0:nqmx,0:7,num_spe)
real*8, intent(in)::atx(natom),aty(natom),atz(natom)
real*8, intent(inout)::vtmp_pp(-npxmax+1:npxmax,-npymax+1:npymax,-npzmax+1:npzmax)
real*8, intent(inout)::vtmp_pp1(-nxspd+1:nxspd,-nyspd+1:nyspd,-npzmax+1:npzmax)
real*8, intent(inout)::vtmp_pp2(-nxspd+1:nxspd,-npymax+1:npymax,-npzmax+1:npzmax)
real*8, intent(inout)::vnlspd(-nxspd+1:nxspd,-nyspd+1:nyspd,-nzspd+1:nzspd,nprjmx)
real*8, intent(out)::vnlocp(num_list,nprjmx,num_ppcell)
integer npxmax_d,npymax_d,npzmax_d,npatx,npaty,npatz
real*8 dpx,dpy,dpz
real*8 pi,sqfourpi,sq3fourpi,sq1516pi,sq1504pi,sq0516pi
real*8 twopicbin,fourpi,rcutss,tmp,x,y,z,r,dqq,rin &
      ,tmpspd1,tmpspd2,tmpspd3,tmpspd4,tmpspd5,tmpspd6,qqq,qqq2,qqqr,qqqr2,qqqrin,qqqr2in,qqqr3in &
      ,sinqqqr,cosqqqr,besselj0,besselj1,besselj2,dpxyz,weight
real*8 ylm01,ylm02,ylm03,ylm04,ylm05,ylm06,ylm07,ylm08,ylm09
integer iaps,na,ix,iy,iz,kx,ky,kz,iq,ixyz,kkx,kky,kkz,jx,jy,jz,la,js,je


  pi=dacos(-1.0d0)
  sqfourpi=dsqrt(1.0d0/4.0d0/pi)
  sq3fourpi=dsqrt(3.0d0/4.0d0/pi)
  sq1516pi=dsqrt(15.0d0/16.0d0/pi)
  sq1504pi=dsqrt(15.0d0/4.0d0/pi)
  sq0516pi=dsqrt( 5.0d0/16.0d0/pi)
  dpx=dx/npmesh
  dpy=dy/npmesh
  dpz=dz/npmesh
  dpxyz=dpx*dpy*dpz
  npxmax_d=npxmax*npmesh
  npymax_d=npymax*npmesh
  npzmax_d=npzmax*npmesh
  twopicbin=1.0d0/(2.0d0*pi)**3
  fourpi=4.0d0*pi

!$omp do
  do ix=1,num_list*nprjmx*num_ppcell
    vnlocp(ix,1,1)=0.0d0
  enddo
  do na=1,natom
!$omp do
    do ix=1,8*nxspd*nyspd*nzspd*nprjmx
      vnlspd(-nxspd+ix,-nyspd+1,-nzspd+1,1)=0.0d0
    enddo
    if ((natpri(na) .eq. key_natpri_in) .or. (natpri(na) .eq. key_natpri_inps)) then
      npatx=int((atx(na)+nxmax*dx+0.5d0*dpx)/dpx)-nxmax*npmesh
      npaty=int((aty(na)+nymax*dy+0.5d0*dpy)/dpy)-nymax*npmesh
      npatz=int((atz(na)+nzmax*dz+0.5d0*dpz)/dpz)-nzmax*npmesh

      iaps=naps(na)
! ==========  assign nonlocal pseudopotential  ==========
      rcutss=radial(nradct(indspe(na)),indspe(na))*psftrad
! ==========  for nprj= 1  ==========
      if (nprj(indspe(na)) .eq. 1) then
        dqq=pi/rcutss
        tmp=fourpi*dqq*twopicbin
!$omp do
        do iz=-npzmax_d+1,npzmax_d
        do iy=-npymax_d+1,npymax_d
        do ix=-npxmax_d+1,npxmax_d
          kx=ix+(npatx-npmesh*natx(na))
          ky=iy+(npaty-npmesh*naty(na))
          kz=iz+(npatz-npmesh*natz(na))
          x=ix*dpx-(atx(na)-(npatx*dpx-0.5d0*dpx))
          y=iy*dpy-(aty(na)-(npaty*dpy-0.5d0*dpy))
          z=iz*dpz-(atz(na)-(npatz*dpz-0.5d0*dpz))
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
            weight=1.0d0
            if (r .gt. eps) then
              if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,radial(nradct(indspe(na)),indspe(na)), weight)
              rin=1.0d0/r
              ylm01=1.0d0
              tmpspd1=0.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                qqqr=qqq*r
                qqqrin=rin/qqq
                sinqqqr=dsin(qqqr)
                besselj0=sinqqqr*qqqrin
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*besselj0*qqq*qqq
              end do
            else
              ylm01=1.0d0
              tmpspd1=0.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*qqq*qqq
              end do
            end if
            vnlspd(kx,ky,kz, 1)=tmpspd1*tmp*dpxyz*sqfourpi*ylm01*weight
          end if
        end do
        end do
        end do
      end if
! ===================================
! ==========  for nprj= 2  ==========
      if (nprj(indspe(na)) .eq. 2) then
        dqq=pi/rcutss
        tmp=fourpi*dqq*twopicbin
!$omp do
        do iz=-npzmax_d+1,npzmax_d
        do iy=-npymax_d+1,npymax_d
        do ix=-npxmax_d+1,npxmax_d
          kx=ix+(npatx-npmesh*natx(na))
          ky=iy+(npaty-npmesh*naty(na))
          kz=iz+(npatz-npmesh*natz(na))
          x=ix*dpx-(atx(na)-(npatx*dpx-0.5d0*dpx))
          y=iy*dpy-(aty(na)-(npaty*dpy-0.5d0*dpy))
          z=iz*dpz-(atz(na)-(npatz*dpz-0.5d0*dpz))
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
            tmpspd1=0.0d0
            tmpspd2=0.0d0
            weight=1.0d0
            if (r .gt. eps) then
              if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,radial(nradct(indspe(na)),indspe(na)), weight)
              rin=1.0d0/r
              ylm01=1.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                qqq2=qqq*qqq
                qqqr=qqq*r
                qqqrin=rin/qqq
                qqqr2in=qqqrin*qqqrin
                sinqqqr=dsin(qqqr)
                cosqqqr=dcos(qqqr)
                besselj0=sinqqqr*qqqrin
                besselj1=(sinqqqr-qqqr*cosqqqr)*qqqr2in
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*besselj0*qqq2
                tmpspd2=tmpspd2+coef(iq,2,indspe(na))*besselj0*qqq2
              end do
            else
              ylm01=1.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*qqq*qqq
                tmpspd2=tmpspd2+coef(iq,2,indspe(na))*qqq*qqq
              end do
            end if
            vnlspd(kx,ky,kz, 1)=tmpspd1*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 2)=tmpspd2*tmp*dpxyz*sqfourpi*ylm01*weight
          end if
        end do
        end do
        end do
      end if
!  ===================================
! ==========  for nprj= 4  ==========
      if (nprj(indspe(na)) .eq. 4) then
        dqq=pi/rcutss
        tmp=fourpi*dqq*twopicbin
!$omp do
        do iz=-npzmax_d+1,npzmax_d
        do iy=-npymax_d+1,npymax_d
        do ix=-npxmax_d+1,npxmax_d
          kx=ix+(npatx-npmesh*natx(na))
          ky=iy+(npaty-npmesh*naty(na))
          kz=iz+(npatz-npmesh*natz(na))
          x=ix*dpx-(atx(na)-(npatx*dpx-0.5d0*dpx))
          y=iy*dpy-(aty(na)-(npaty*dpy-0.5d0*dpy))
          z=iz*dpz-(atz(na)-(npatz*dpz-0.5d0*dpz))
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
            tmpspd1=0.0d0
            tmpspd3=0.0d0
            weight=1.0d0
            if (r .gt. eps) then
              if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,radial(nradct(indspe(na)),indspe(na)), weight)
              rin=1.0d0/r
              ylm01=1.0d0
              ylm02=rin*x
              ylm03=rin*y
              ylm04=rin*z
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                qqq2=qqq*qqq
                qqqr=qqq*r
                qqqrin=rin/qqq
                qqqr2in=qqqrin*qqqrin
                sinqqqr=dsin(qqqr)
                cosqqqr=dcos(qqqr)
                besselj0=sinqqqr*qqqrin
                besselj1=(sinqqqr-qqqr*cosqqqr)*qqqr2in
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*besselj0*qqq2
                tmpspd3=tmpspd3+coef(iq,3,indspe(na))*besselj1*qqq2
              end do
            else
              ylm01=1.0d0
              ylm02=0.0d0
              ylm03=0.0d0
              ylm04=0.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*qqq*qqq
              end do
            end if
            vnlspd(kx,ky,kz, 1)=tmpspd1*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 2)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 3)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 4)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm04*weight
          end if
        end do
        end do
        end do
      end if
! ===================================
! ==========  for nprj= 5  ==========
      if (nprj(indspe(na)) .eq. 5) then
        dqq=pi/rcutss
        tmp=fourpi*dqq*twopicbin
!$omp do
        do iz=-npzmax_d+1,npzmax_d
        do iy=-npymax_d+1,npymax_d
        do ix=-npxmax_d+1,npxmax_d
          kx=ix+(npatx-npmesh*natx(na))
          ky=iy+(npaty-npmesh*naty(na))
          kz=iz+(npatz-npmesh*natz(na))
          x=ix*dpx-(atx(na)-(npatx*dpx-0.5d0*dpx))
          y=iy*dpy-(aty(na)-(npaty*dpy-0.5d0*dpy))
          z=iz*dpz-(atz(na)-(npatz*dpz-0.5d0*dpz))
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
            tmpspd1=0.0d0
            tmpspd2=0.0d0
            tmpspd3=0.0d0
            weight=1.0d0
            if (r .gt. eps) then
              if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,radial(nradct(indspe(na)),indspe(na)), weight)
              rin=1.0d0/r
              ylm01=1.0d0
              ylm02=rin*x
              ylm03=rin*y
              ylm04=rin*z
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                qqq2=qqq*qqq
                qqqr=qqq*r
                qqqrin=rin/qqq
                qqqr2in=qqqrin*qqqrin
                sinqqqr=dsin(qqqr)
                cosqqqr=dcos(qqqr)
                besselj0=sinqqqr*qqqrin
                besselj1=(sinqqqr-qqqr*cosqqqr)*qqqr2in
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*besselj0*qqq2
                tmpspd2=tmpspd2+coef(iq,2,indspe(na))*besselj0*qqq2
                tmpspd3=tmpspd3+coef(iq,3,indspe(na))*besselj1*qqq2
              end do
            else
              ylm01=1.0d0
              ylm02=0.0d0
              ylm03=0.0d0
              ylm04=0.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*qqq*qqq
                tmpspd2=tmpspd2+coef(iq,2,indspe(na))*qqq*qqq
              end do
            end if
            vnlspd(kx,ky,kz, 1)=tmpspd1*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 2)=tmpspd2*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 3)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 4)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 5)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm04*weight
          end if
        end do
        end do
        end do
      end if
! ===================================
! ==========  for nprj= 7  ==========
      if (nprj(indspe(na)) .eq. 7) then
        dqq=pi/rcutss
        tmp=fourpi*dqq*twopicbin
!$omp do
        do iz=-npzmax_d+1,npzmax_d
        do iy=-npymax_d+1,npymax_d
        do ix=-npxmax_d+1,npxmax_d
          kx=ix+(npatx-npmesh*natx(na))
          ky=iy+(npaty-npmesh*naty(na))
          kz=iz+(npatz-npmesh*natz(na))
          x=ix*dpx-(atx(na)-(npatx*dpx-0.5d0*dpx))
          y=iy*dpy-(aty(na)-(npaty*dpy-0.5d0*dpy))
          z=iz*dpz-(atz(na)-(npatz*dpz-0.5d0*dpz))
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
            tmpspd1=0.0d0
            tmpspd3=0.0d0
            tmpspd4=0.0d0
            weight=1.0d0
            if (r .gt. eps) then
              if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,radial(nradct(indspe(na)),indspe(na)), weight)
              rin=1.0d0/r
              ylm01=1.0d0
              ylm02=rin*x
              ylm03=rin*y
              ylm04=rin*z
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                qqq2=qqq*qqq
                qqqr=qqq*r
                qqqrin=rin/qqq
                qqqr2in=qqqrin*qqqrin
                sinqqqr=dsin(qqqr)
                cosqqqr=dcos(qqqr)
                besselj0=sinqqqr*qqqrin
                besselj1=(sinqqqr-qqqr*cosqqqr)*qqqr2in
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*besselj0*qqq2
                tmpspd3=tmpspd3+coef(iq,3,indspe(na))*besselj1*qqq2
                tmpspd4=tmpspd4+coef(iq,4,indspe(na))*besselj1*qqq2
              end do
            else
              ylm01=1.0d0
              ylm02=0.0d0
              ylm03=0.0d0
              ylm04=0.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*qqq*qqq
              end do
            end if
            vnlspd(kx,ky,kz, 1)=tmpspd1*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 2)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 3)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 4)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm04*weight
            vnlspd(kx,ky,kz, 5)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 6)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 7)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm04*weight
          end if
        end do
        end do
        end do
      end if
! ===================================
! ==========  for nprj= 8  ==========
      if (nprj(indspe(na)) .eq. 8) then
        dqq=pi/rcutss
        tmp=fourpi*dqq*twopicbin
!$omp do
        do iz=-npzmax_d+1,npzmax_d
        do iy=-npymax_d+1,npymax_d
        do ix=-npxmax_d+1,npxmax_d
          kx=ix+(npatx-npmesh*natx(na))
          ky=iy+(npaty-npmesh*naty(na))
          kz=iz+(npatz-npmesh*natz(na))
          x=ix*dpx-(atx(na)-(npatx*dpx-0.5d0*dpx))
          y=iy*dpy-(aty(na)-(npaty*dpy-0.5d0*dpy))
          z=iz*dpz-(atz(na)-(npatz*dpz-0.5d0*dpz))
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
            tmpspd1=0.0d0
            tmpspd2=0.0d0
            tmpspd3=0.0d0
            tmpspd4=0.0d0
            weight=1.0d0
            if (r .gt. eps) then
              if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,radial(nradct(indspe(na)),indspe(na)), weight)
              rin=1.0d0/r
              ylm01=1.0d0
              ylm02=rin*x
              ylm03=rin*y
              ylm04=rin*z
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                qqq2=qqq*qqq
                qqqr=qqq*r
                qqqrin=rin/qqq
                qqqr2in=qqqrin*qqqrin
                sinqqqr=dsin(qqqr)
                cosqqqr=dcos(qqqr)
                besselj0=sinqqqr*qqqrin
                besselj1=(sinqqqr-qqqr*cosqqqr)*qqqr2in
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*besselj0*qqq2
                tmpspd2=tmpspd2+coef(iq,2,indspe(na))*besselj0*qqq2
                tmpspd3=tmpspd3+coef(iq,3,indspe(na))*besselj1*qqq2
                tmpspd4=tmpspd4+coef(iq,4,indspe(na))*besselj1*qqq2
              end do
            else
              ylm01=1.0d0
              ylm02=0.0d0
              ylm03=0.0d0
              ylm04=0.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*qqq*qqq
                tmpspd2=tmpspd2+coef(iq,2,indspe(na))*qqq*qqq
              end do
            end if
            vnlspd(kx,ky,kz, 1)=tmpspd1*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 2)=tmpspd2*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 3)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 4)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 5)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm04*weight
            vnlspd(kx,ky,kz, 6)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 7)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 8)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm04*weight
          end if
        end do
        end do
        end do
      end if
! ===================================
! ==========  for nprj= 9  ==========
      if (nprj(indspe(na)) .eq. 9) then
        dqq=pi/rcutss
        tmp=fourpi*dqq*twopicbin
!$omp do
        do iz=-npzmax_d+1,npzmax_d
        do iy=-npymax_d+1,npymax_d
        do ix=-npxmax_d+1,npxmax_d
          kx=ix+(npatx-npmesh*natx(na))
          ky=iy+(npaty-npmesh*naty(na))
          kz=iz+(npatz-npmesh*natz(na))
          x=ix*dpx-(atx(na)-(npatx*dpx-0.5d0*dpx))
          y=iy*dpy-(aty(na)-(npaty*dpy-0.5d0*dpy))
          z=iz*dpz-(atz(na)-(npatz*dpz-0.5d0*dpz))
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
            tmpspd1=0.0d0
            tmpspd3=0.0d0
            tmpspd5=0.0d0
            weight=1.0d0
            if (r .gt. eps) then
              if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,radial(nradct(indspe(na)),indspe(na)), weight)
              rin=1.0d0/r
              ylm01=1.0d0
              ylm02=rin*x
              ylm03=rin*y
              ylm04=rin*z
              ylm05= rin*rin*(x*x-y*y)
              ylm06=-rin*rin*z*x
              ylm07= rin*rin*(3.0d0*z*z-r*r)
              ylm08=-rin*rin*y*z
              ylm09= rin*rin*x*y
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                qqq2=qqq*qqq
                qqqr=qqq*r
                qqqr2=qqqr*qqqr
                qqqrin=rin/qqq
                qqqr2in=qqqrin*qqqrin
                qqqr3in=qqqr2in*qqqrin
                sinqqqr=dsin(qqqr)
                cosqqqr=dcos(qqqr)
                besselj0=sinqqqr*qqqrin
                besselj1=(sinqqqr-qqqr*cosqqqr)*qqqr2in
                besselj2=((3.0d0-qqqr2)*sinqqqr-3.0d0*qqqr*cosqqqr)*qqqr3in
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*besselj0*qqq2
                tmpspd3=tmpspd3+coef(iq,3,indspe(na))*besselj1*qqq2
                tmpspd5=tmpspd5+coef(iq,5,indspe(na))*besselj2*qqq2
              end do
            else
              ylm01=1.0d0
              ylm02=0.0d0
              ylm03=0.0d0
              ylm04=0.0d0
              ylm05=0.0d0
              ylm06=0.0d0
              ylm07=0.0d0
              ylm08=0.0d0
              ylm09=0.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*qqq*qqq
              end do
            end if
            vnlspd(kx,ky,kz, 1)=tmpspd1*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 2)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 3)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 4)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm04*weight
            vnlspd(kx,ky,kz, 5)=tmpspd5*tmp*dpxyz*sq1516pi*ylm05*weight
            vnlspd(kx,ky,kz, 6)=tmpspd5*tmp*dpxyz*sq1504pi*ylm06*weight
            vnlspd(kx,ky,kz, 7)=tmpspd5*tmp*dpxyz*sq0516pi*ylm07*weight
            vnlspd(kx,ky,kz, 8)=tmpspd5*tmp*dpxyz*sq1504pi*ylm08*weight
            vnlspd(kx,ky,kz, 9)=tmpspd5*tmp*dpxyz*sq1504pi*ylm09*weight
          end if
        end do
        end do
        end do
      end if
! ===================================
! ==========  for nprj=10  ==========
      if (nprj(indspe(na)) .eq. 10) then
        dqq=pi/rcutss
        tmp=fourpi*dqq*twopicbin
!$omp do
        do iz=-npzmax_d+1,npzmax_d
        do iy=-npymax_d+1,npymax_d
        do ix=-npxmax_d+1,npxmax_d
          kx=ix+(npatx-npmesh*natx(na))
          ky=iy+(npaty-npmesh*naty(na))
          kz=iz+(npatz-npmesh*natz(na))
          x=ix*dpx-(atx(na)-(npatx*dpx-0.5d0*dpx))
          y=iy*dpy-(aty(na)-(npaty*dpy-0.5d0*dpy))
          z=iz*dpz-(atz(na)-(npatz*dpz-0.5d0*dpz))
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
            tmpspd1=0.0d0
            tmpspd2=0.0d0
            tmpspd3=0.0d0
            tmpspd5=0.0d0
            weight=1.0d0
            if (r .gt. eps) then
              if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,radial(nradct(indspe(na)),indspe(na)), weight)
              rin=1.0d0/r
              ylm01=1.0d0
              ylm02=rin*x
              ylm03=rin*y
              ylm04=rin*z
              ylm05= rin*rin*(x*x-y*y)
              ylm06=-rin*rin*z*x
              ylm07= rin*rin*(3.0d0*z*z-r*r)
              ylm08=-rin*rin*y*z
              ylm09= rin*rin*x*y
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                qqq2=qqq*qqq
                qqqr=qqq*r
                qqqr2=qqqr*qqqr
                qqqrin=rin/qqq
                qqqr2in=qqqrin*qqqrin
                qqqr3in=qqqr2in*qqqrin
                sinqqqr=dsin(qqqr)
                cosqqqr=dcos(qqqr)
                besselj0=sinqqqr*qqqrin
                besselj1=(sinqqqr-qqqr*cosqqqr)*qqqr2in
                besselj2=((3.0d0-qqqr2)*sinqqqr-3.0d0*qqqr*cosqqqr)*qqqr3in
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*besselj0*qqq2
                tmpspd2=tmpspd2+coef(iq,2,indspe(na))*besselj0*qqq2
                tmpspd3=tmpspd3+coef(iq,3,indspe(na))*besselj1*qqq2
                tmpspd5=tmpspd5+coef(iq,5,indspe(na))*besselj2*qqq2
              end do
            else
              ylm01=1.0d0
              ylm02=0.0d0
              ylm03=0.0d0
              ylm04=0.0d0
              ylm05=0.0d0
              ylm06=0.0d0
              ylm07=0.0d0
              ylm08=0.0d0
              ylm09=0.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*qqq*qqq
                tmpspd2=tmpspd2+coef(iq,2,indspe(na))*qqq*qqq
              end do
            end if
            vnlspd(kx,ky,kz, 1)=tmpspd1*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 2)=tmpspd2*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 3)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 4)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 5)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm04*weight
            vnlspd(kx,ky,kz, 6)=tmpspd5*tmp*dpxyz*sq1516pi*ylm05*weight
            vnlspd(kx,ky,kz, 7)=tmpspd5*tmp*dpxyz*sq1504pi*ylm06*weight
            vnlspd(kx,ky,kz, 8)=tmpspd5*tmp*dpxyz*sq0516pi*ylm07*weight
            vnlspd(kx,ky,kz, 9)=tmpspd5*tmp*dpxyz*sq1504pi*ylm08*weight
            vnlspd(kx,ky,kz,10)=tmpspd5*tmp*dpxyz*sq1504pi*ylm09*weight
          end if
        end do
        end do
        end do
      end if
! ===================================
! ==========  for nprj=12  ==========
      if (nprj(indspe(na)) .eq. 12) then
        dqq=pi/rcutss
        tmp=fourpi*dqq*twopicbin
!$omp do
        do iz=-npzmax_d+1,npzmax_d
        do iy=-npymax_d+1,npymax_d
        do ix=-npxmax_d+1,npxmax_d
          kx=ix+(npatx-npmesh*natx(na))
          ky=iy+(npaty-npmesh*naty(na))
          kz=iz+(npatz-npmesh*natz(na))
          x=ix*dpx-(atx(na)-(npatx*dpx-0.5d0*dpx))
          y=iy*dpy-(aty(na)-(npaty*dpy-0.5d0*dpy))
          z=iz*dpz-(atz(na)-(npatz*dpz-0.5d0*dpz))
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
            tmpspd1=0.0d0
            tmpspd3=0.0d0
            tmpspd4=0.0d0
            tmpspd5=0.0d0
            weight=1.0d0
            if (r .gt. eps) then
              if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,radial(nradct(indspe(na)),indspe(na)), weight)
              rin=1.0d0/r
              ylm01=1.0d0
              ylm02=rin*x
              ylm03=rin*y
              ylm04=rin*z
              ylm05= rin*rin*(x*x-y*y)
              ylm06=-rin*rin*z*x
              ylm07= rin*rin*(3.0d0*z*z-r*r)
              ylm08=-rin*rin*y*z
              ylm09= rin*rin*x*y
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                qqq2=qqq*qqq
                qqqr=qqq*r
                qqqr2=qqqr*qqqr
                qqqrin=rin/qqq
                qqqr2in=qqqrin*qqqrin
                qqqr3in=qqqr2in*qqqrin
                sinqqqr=dsin(qqqr)
                cosqqqr=dcos(qqqr)
                besselj0=sinqqqr*qqqrin
                besselj1=(sinqqqr-qqqr*cosqqqr)*qqqr2in
                besselj2=((3.0d0-qqqr2)*sinqqqr-3.0d0*qqqr*cosqqqr)*qqqr3in
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*besselj0*qqq2
                tmpspd3=tmpspd3+coef(iq,3,indspe(na))*besselj1*qqq2
                tmpspd4=tmpspd4+coef(iq,4,indspe(na))*besselj1*qqq2
                tmpspd5=tmpspd5+coef(iq,5,indspe(na))*besselj2*qqq2
              end do
            else
              ylm01=1.0d0
              ylm02=0.0d0
              ylm03=0.0d0
              ylm04=0.0d0
              ylm05=0.0d0
              ylm06=0.0d0
              ylm07=0.0d0
              ylm08=0.0d0
              ylm09=0.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*qqq*qqq
              end do
            end if
            vnlspd(kx,ky,kz, 1)=tmpspd1*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 2)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 3)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 4)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm04*weight
            vnlspd(kx,ky,kz, 5)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 6)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 7)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm04*weight
            vnlspd(kx,ky,kz, 8)=tmpspd5*tmp*dpxyz*sq1516pi*ylm05*weight
            vnlspd(kx,ky,kz, 9)=tmpspd5*tmp*dpxyz*sq1504pi*ylm06*weight
            vnlspd(kx,ky,kz,10)=tmpspd5*tmp*dpxyz*sq0516pi*ylm07*weight
            vnlspd(kx,ky,kz,11)=tmpspd5*tmp*dpxyz*sq1504pi*ylm08*weight
            vnlspd(kx,ky,kz,12)=tmpspd5*tmp*dpxyz*sq1504pi*ylm09*weight
          end if
        end do
        end do
        end do
      end if
! ===================================
! ==========  for nprj=13  ==========
      if (nprj(indspe(na)) .eq. 13) then
        dqq=pi/rcutss
        tmp=fourpi*dqq*twopicbin
!$omp do
        do iz=-npzmax_d+1,npzmax_d
        do iy=-npymax_d+1,npymax_d
        do ix=-npxmax_d+1,npxmax_d
          kx=ix+(npatx-npmesh*natx(na))
          ky=iy+(npaty-npmesh*naty(na))
          kz=iz+(npatz-npmesh*natz(na))
          x=ix*dpx-(atx(na)-(npatx*dpx-0.5d0*dpx))
          y=iy*dpy-(aty(na)-(npaty*dpy-0.5d0*dpy))
          z=iz*dpz-(atz(na)-(npatz*dpz-0.5d0*dpz))
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
            tmpspd1=0.0d0
            tmpspd2=0.0d0
            tmpspd3=0.0d0
            tmpspd4=0.0d0
            tmpspd5=0.0d0
            weight=1.0d0
            if (r .gt. eps) then
              if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,radial(nradct(indspe(na)),indspe(na)), weight)
              rin=1.0d0/r
              ylm01=1.0d0
              ylm02=rin*x
              ylm03=rin*y
              ylm04=rin*z
              ylm05= rin*rin*(x*x-y*y)
              ylm06=-rin*rin*z*x
              ylm07= rin*rin*(3.0d0*z*z-r*r)
              ylm08=-rin*rin*y*z
              ylm09= rin*rin*x*y
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                qqq2=qqq*qqq
                qqqr=qqq*r
                qqqr2=qqqr*qqqr
                qqqrin=rin/qqq
                qqqr2in=qqqrin*qqqrin
                qqqr3in=qqqr2in*qqqrin
                sinqqqr=dsin(qqqr)
                cosqqqr=dcos(qqqr)
                besselj0=sinqqqr*qqqrin
                besselj1=(sinqqqr-qqqr*cosqqqr)*qqqr2in
                besselj2=((3.0d0-qqqr2)*sinqqqr-3.0d0*qqqr*cosqqqr)*qqqr3in
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*besselj0*qqq2
                tmpspd2=tmpspd2+coef(iq,2,indspe(na))*besselj0*qqq2
                tmpspd3=tmpspd3+coef(iq,3,indspe(na))*besselj1*qqq2
                tmpspd4=tmpspd4+coef(iq,4,indspe(na))*besselj1*qqq2
                tmpspd5=tmpspd5+coef(iq,5,indspe(na))*besselj2*qqq2
              end do
            else
              ylm01=1.0d0
              ylm02=0.0d0
              ylm03=0.0d0
              ylm04=0.0d0
              ylm05=0.0d0
              ylm06=0.0d0
              ylm07=0.0d0
              ylm08=0.0d0
              ylm09=0.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*qqq*qqq
                tmpspd2=tmpspd2+coef(iq,2,indspe(na))*qqq*qqq
              end do
            end if
            vnlspd(kx,ky,kz, 1)=tmpspd1*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 2)=tmpspd2*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 3)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 4)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 5)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm04*weight
            vnlspd(kx,ky,kz, 6)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 7)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 8)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm04*weight
            vnlspd(kx,ky,kz, 9)=tmpspd5*tmp*dpxyz*sq1516pi*ylm05*weight
            vnlspd(kx,ky,kz,10)=tmpspd5*tmp*dpxyz*sq1504pi*ylm06*weight
            vnlspd(kx,ky,kz,11)=tmpspd5*tmp*dpxyz*sq0516pi*ylm07*weight
            vnlspd(kx,ky,kz,12)=tmpspd5*tmp*dpxyz*sq1504pi*ylm08*weight
            vnlspd(kx,ky,kz,13)=tmpspd5*tmp*dpxyz*sq1504pi*ylm09*weight
          end if
        end do
        end do
        end do
      end if
! ===================================
! ==========  for nprj=14  ==========
      if (nprj(indspe(na)) .eq. 14) then
        dqq=pi/rcutss
        tmp=fourpi*dqq*twopicbin
!$omp do
        do iz=-npzmax_d+1,npzmax_d
        do iy=-npymax_d+1,npymax_d
        do ix=-npxmax_d+1,npxmax_d
          kx=ix+(npatx-npmesh*natx(na))
          ky=iy+(npaty-npmesh*naty(na))
          kz=iz+(npatz-npmesh*natz(na))
          x=ix*dpx-(atx(na)-(npatx*dpx-0.5d0*dpx))
          y=iy*dpy-(aty(na)-(npaty*dpy-0.5d0*dpy))
          z=iz*dpz-(atz(na)-(npatz*dpz-0.5d0*dpz))
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
            tmpspd1=0.0d0
            tmpspd3=0.0d0
            tmpspd5=0.0d0
            tmpspd6=0.0d0
            weight=1.0d0
            if (r .gt. eps) then
              if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,radial(nradct(indspe(na)),indspe(na)), weight)
              rin=1.0d0/r
              ylm01=1.0d0
              ylm02=rin*x
              ylm03=rin*y
              ylm04=rin*z
              ylm05= rin*rin*(x*x-y*y)
              ylm06=-rin*rin*z*x
              ylm07= rin*rin*(3.0d0*z*z-r*r)
              ylm08=-rin*rin*y*z
              ylm09= rin*rin*x*y
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                qqq2=qqq*qqq
                qqqr=qqq*r
                qqqr2=qqqr*qqqr
                qqqrin=rin/qqq
                qqqr2in=qqqrin*qqqrin
                qqqr3in=qqqr2in*qqqrin
                sinqqqr=dsin(qqqr)
                cosqqqr=dcos(qqqr)
                besselj0=sinqqqr*qqqrin
                besselj1=(sinqqqr-qqqr*cosqqqr)*qqqr2in
                besselj2=((3.0d0-qqqr2)*sinqqqr-3.0d0*qqqr*cosqqqr)*qqqr3in
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*besselj0*qqq2
                tmpspd3=tmpspd3+coef(iq,3,indspe(na))*besselj1*qqq2
                tmpspd5=tmpspd5+coef(iq,5,indspe(na))*besselj2*qqq2
                tmpspd6=tmpspd6+coef(iq,6,indspe(na))*besselj2*qqq2
              end do
            else
              ylm01=1.0d0
              ylm02=0.0d0
              ylm03=0.0d0
              ylm04=0.0d0
              ylm05=0.0d0
              ylm06=0.0d0
              ylm07=0.0d0
              ylm08=0.0d0
              ylm09=0.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*qqq*qqq
              end do
            end if
            vnlspd(kx,ky,kz, 1)=tmpspd1*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 2)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 3)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 4)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm04*weight
            vnlspd(kx,ky,kz, 5)=tmpspd5*tmp*dpxyz*sq1516pi*ylm05*weight
            vnlspd(kx,ky,kz, 6)=tmpspd5*tmp*dpxyz*sq1504pi*ylm06*weight
            vnlspd(kx,ky,kz, 7)=tmpspd5*tmp*dpxyz*sq0516pi*ylm07*weight
            vnlspd(kx,ky,kz, 8)=tmpspd5*tmp*dpxyz*sq1504pi*ylm08*weight
            vnlspd(kx,ky,kz, 9)=tmpspd5*tmp*dpxyz*sq1504pi*ylm09*weight
            vnlspd(kx,ky,kz,10)=tmpspd6*tmp*dpxyz*sq1516pi*ylm05*weight
            vnlspd(kx,ky,kz,11)=tmpspd6*tmp*dpxyz*sq1504pi*ylm06*weight
            vnlspd(kx,ky,kz,12)=tmpspd6*tmp*dpxyz*sq0516pi*ylm07*weight
            vnlspd(kx,ky,kz,13)=tmpspd6*tmp*dpxyz*sq1504pi*ylm08*weight
            vnlspd(kx,ky,kz,14)=tmpspd6*tmp*dpxyz*sq1504pi*ylm09*weight
          end if
        end do
        end do
        end do
      end if
! ===================================
! ==========  for nprj=15  ==========
      if (nprj(indspe(na)) .eq. 15) then
        dqq=pi/rcutss
        tmp=fourpi*dqq*twopicbin
!$omp do
        do iz=-npzmax_d+1,npzmax_d
        do iy=-npymax_d+1,npymax_d
        do ix=-npxmax_d+1,npxmax_d
          kx=ix+(npatx-npmesh*natx(na))
          ky=iy+(npaty-npmesh*naty(na))
          kz=iz+(npatz-npmesh*natz(na))
          x=ix*dpx-(atx(na)-(npatx*dpx-0.5d0*dpx))
          y=iy*dpy-(aty(na)-(npaty*dpy-0.5d0*dpy))
          z=iz*dpz-(atz(na)-(npatz*dpz-0.5d0*dpz))
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
            tmpspd1=0.0d0
            tmpspd2=0.0d0
            tmpspd3=0.0d0
            tmpspd5=0.0d0
            tmpspd6=0.0d0
            weight=1.0d0
            if (r .gt. eps) then
              if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,radial(nradct(indspe(na)),indspe(na)), weight)
              rin=1.0d0/r
              ylm01=1.0d0
              ylm02=rin*x
              ylm03=rin*y
              ylm04=rin*z
              ylm05= rin*rin*(x*x-y*y)
              ylm06=-rin*rin*z*x
              ylm07= rin*rin*(3.0d0*z*z-r*r)
              ylm08=-rin*rin*y*z
              ylm09= rin*rin*x*y
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                qqq2=qqq*qqq
                qqqr=qqq*r
                qqqr2=qqqr*qqqr
                qqqrin=rin/qqq
                qqqr2in=qqqrin*qqqrin
                qqqr3in=qqqr2in*qqqrin
                sinqqqr=dsin(qqqr)
                cosqqqr=dcos(qqqr)
                besselj0=sinqqqr*qqqrin
                besselj1=(sinqqqr-qqqr*cosqqqr)*qqqr2in
                besselj2=((3.0d0-qqqr2)*sinqqqr-3.0d0*qqqr*cosqqqr)*qqqr3in
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*besselj0*qqq2
                tmpspd2=tmpspd2+coef(iq,2,indspe(na))*besselj0*qqq2
                tmpspd3=tmpspd3+coef(iq,3,indspe(na))*besselj1*qqq2
                tmpspd5=tmpspd5+coef(iq,5,indspe(na))*besselj2*qqq2
                tmpspd6=tmpspd6+coef(iq,6,indspe(na))*besselj2*qqq2
              end do
            else
              ylm01=1.0d0
              ylm02=0.0d0
              ylm03=0.0d0
              ylm04=0.0d0
              ylm05=0.0d0
              ylm06=0.0d0
              ylm07=0.0d0
              ylm08=0.0d0
              ylm09=0.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*qqq*qqq
                tmpspd2=tmpspd2+coef(iq,2,indspe(na))*qqq*qqq
              end do
            end if
            vnlspd(kx,ky,kz, 1)=tmpspd1*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 2)=tmpspd2*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 3)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 4)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 5)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm04*weight
            vnlspd(kx,ky,kz, 6)=tmpspd5*tmp*dpxyz*sq1516pi*ylm05*weight
            vnlspd(kx,ky,kz, 7)=tmpspd5*tmp*dpxyz*sq1504pi*ylm06*weight
            vnlspd(kx,ky,kz, 8)=tmpspd5*tmp*dpxyz*sq0516pi*ylm07*weight
            vnlspd(kx,ky,kz, 9)=tmpspd5*tmp*dpxyz*sq1504pi*ylm08*weight
            vnlspd(kx,ky,kz,10)=tmpspd5*tmp*dpxyz*sq1504pi*ylm09*weight
            vnlspd(kx,ky,kz,11)=tmpspd6*tmp*dpxyz*sq1516pi*ylm05*weight
            vnlspd(kx,ky,kz,12)=tmpspd6*tmp*dpxyz*sq1504pi*ylm06*weight
            vnlspd(kx,ky,kz,13)=tmpspd6*tmp*dpxyz*sq0516pi*ylm07*weight
            vnlspd(kx,ky,kz,14)=tmpspd6*tmp*dpxyz*sq1504pi*ylm08*weight
            vnlspd(kx,ky,kz,15)=tmpspd6*tmp*dpxyz*sq1504pi*ylm09*weight
          end if
        end do
        end do
        end do
      end if
! ===================================
! ==========  for nprj=17  ==========
      if (nprj(indspe(na)) .eq. 17) then
        dqq=pi/rcutss
        tmp=fourpi*dqq*twopicbin
!$omp do
        do iz=-npzmax_d+1,npzmax_d
        do iy=-npymax_d+1,npymax_d
        do ix=-npxmax_d+1,npxmax_d
          kx=ix+(npatx-npmesh*natx(na))
          ky=iy+(npaty-npmesh*naty(na))
          kz=iz+(npatz-npmesh*natz(na))
          x=ix*dpx-(atx(na)-(npatx*dpx-0.5d0*dpx))
          y=iy*dpy-(aty(na)-(npaty*dpy-0.5d0*dpy))
          z=iz*dpz-(atz(na)-(npatz*dpz-0.5d0*dpz))
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
            tmpspd1=0.0d0
            tmpspd2=0.0d0
            tmpspd3=0.0d0
            tmpspd4=0.0d0
            tmpspd5=0.0d0
            tmpspd6=0.0d0
            weight=1.0d0
            if (r .gt. eps) then
              if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,radial(nradct(indspe(na)),indspe(na)), weight)
              rin=1.0d0/r
              ylm01=1.0d0
              ylm02=rin*x
              ylm03=rin*y
              ylm04=rin*z
              ylm05= rin*rin*(x*x-y*y)
              ylm06=-rin*rin*z*x
              ylm07= rin*rin*(3.0d0*z*z-r*r)
              ylm08=-rin*rin*y*z
              ylm09= rin*rin*x*y
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                qqq2=qqq*qqq
                qqqr=qqq*r
                qqqr2=qqqr*qqqr
                qqqrin=rin/qqq
                qqqr2in=qqqrin*qqqrin
                qqqr3in=qqqr2in*qqqrin
                sinqqqr=dsin(qqqr)
                cosqqqr=dcos(qqqr)
                besselj0=sinqqqr*qqqrin
                besselj1=(sinqqqr-qqqr*cosqqqr)*qqqr2in
                besselj2=((3.0d0-qqqr2)*sinqqqr-3.0d0*qqqr*cosqqqr)*qqqr3in
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*besselj0*qqq2
                tmpspd3=tmpspd3+coef(iq,3,indspe(na))*besselj1*qqq2
                tmpspd4=tmpspd4+coef(iq,4,indspe(na))*besselj1*qqq2
                tmpspd5=tmpspd5+coef(iq,5,indspe(na))*besselj2*qqq2
                tmpspd6=tmpspd6+coef(iq,6,indspe(na))*besselj2*qqq2
              end do
            else
              ylm01=1.0d0
              ylm02=0.0d0
              ylm03=0.0d0
              ylm04=0.0d0
              ylm05=0.0d0
              ylm06=0.0d0
              ylm07=0.0d0
              ylm08=0.0d0
              ylm09=0.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*qqq*qqq
              end do
            end if
            vnlspd(kx,ky,kz, 1)=tmpspd1*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 2)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 3)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 4)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm04*weight
            vnlspd(kx,ky,kz, 5)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 6)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 7)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm04*weight
            vnlspd(kx,ky,kz, 8)=tmpspd5*tmp*dpxyz*sq1516pi*ylm05*weight
            vnlspd(kx,ky,kz, 9)=tmpspd5*tmp*dpxyz*sq1504pi*ylm06*weight
            vnlspd(kx,ky,kz,10)=tmpspd5*tmp*dpxyz*sq0516pi*ylm07*weight
            vnlspd(kx,ky,kz,11)=tmpspd5*tmp*dpxyz*sq1504pi*ylm08*weight
            vnlspd(kx,ky,kz,12)=tmpspd5*tmp*dpxyz*sq1504pi*ylm09*weight
            vnlspd(kx,ky,kz,13)=tmpspd6*tmp*dpxyz*sq1516pi*ylm05*weight
            vnlspd(kx,ky,kz,14)=tmpspd6*tmp*dpxyz*sq1504pi*ylm06*weight
            vnlspd(kx,ky,kz,15)=tmpspd6*tmp*dpxyz*sq0516pi*ylm07*weight
            vnlspd(kx,ky,kz,16)=tmpspd6*tmp*dpxyz*sq1504pi*ylm08*weight
            vnlspd(kx,ky,kz,17)=tmpspd6*tmp*dpxyz*sq1504pi*ylm09*weight
          end if
        end do
        end do
        end do
      end if
! ===================================
! ==========  for nprj=18  ==========
      if (nprj(indspe(na)) .eq. 18) then
        dqq=pi/rcutss
        tmp=fourpi*dqq*twopicbin
!$omp do
        do iz=-npzmax_d+1,npzmax_d
        do iy=-npymax_d+1,npymax_d
        do ix=-npxmax_d+1,npxmax_d
          kx=ix+(npatx-npmesh*natx(na))
          ky=iy+(npaty-npmesh*naty(na))
          kz=iz+(npatz-npmesh*natz(na))
          x=ix*dpx-(atx(na)-(npatx*dpx-0.5d0*dpx))
          y=iy*dpy-(aty(na)-(npaty*dpy-0.5d0*dpy))
          z=iz*dpz-(atz(na)-(npatz*dpz-0.5d0*dpz))
          r=dsqrt(x*x+y*y+z*z)
          if (r .lt. radial(nradct(indspe(na)),indspe(na))*psctoff) then
            tmpspd1=0.0d0
            tmpspd2=0.0d0
            tmpspd3=0.0d0
            tmpspd4=0.0d0
            tmpspd5=0.0d0
            tmpspd6=0.0d0
            weight=1.0d0
            if (r .gt. eps) then
              if (nfiltyp .eq. 1) call fermidis(1.0d0,r,filpp,radial(nradct(indspe(na)),indspe(na)), weight)
              rin=1.0d0/r
              ylm01=1.0d0
              ylm02=rin*x
              ylm03=rin*y
              ylm04=rin*z
              ylm05= rin*rin*(x*x-y*y)
              ylm06=-rin*rin*z*x
              ylm07= rin*rin*(3.0d0*z*z-r*r)
              ylm08=-rin*rin*y*z
              ylm09= rin*rin*x*y
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                qqq2=qqq*qqq
                qqqr=qqq*r
                qqqr2=qqqr*qqqr
                qqqrin=rin/qqq
                qqqr2in=qqqrin*qqqrin
                qqqr3in=qqqr2in*qqqrin
                sinqqqr=dsin(qqqr)
                cosqqqr=dcos(qqqr)
                besselj0=sinqqqr*qqqrin
                besselj1=(sinqqqr-qqqr*cosqqqr)*qqqr2in
                besselj2=((3.0d0-qqqr2)*sinqqqr-3.0d0*qqqr*cosqqqr)*qqqr3in
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*besselj0*qqq2
                tmpspd2=tmpspd2+coef(iq,2,indspe(na))*besselj0*qqq2
                tmpspd3=tmpspd3+coef(iq,3,indspe(na))*besselj1*qqq2
                tmpspd4=tmpspd4+coef(iq,4,indspe(na))*besselj1*qqq2
                tmpspd5=tmpspd5+coef(iq,5,indspe(na))*besselj2*qqq2
                tmpspd6=tmpspd6+coef(iq,6,indspe(na))*besselj2*qqq2
              end do
            else
              ylm01=1.0d0
              ylm02=0.0d0
              ylm03=0.0d0
              ylm04=0.0d0
              ylm05=0.0d0
              ylm06=0.0d0
              ylm07=0.0d0
              ylm08=0.0d0
              ylm09=0.0d0
              do iq=1,nqct(indspe(na))
                qqq=dqq*iq
                tmpspd1=tmpspd1+coef(iq,1,indspe(na))*qqq*qqq
                tmpspd2=tmpspd2+coef(iq,2,indspe(na))*qqq*qqq
              end do
            end if
            vnlspd(kx,ky,kz, 1)=tmpspd1*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 2)=tmpspd2*tmp*dpxyz*sqfourpi*ylm01*weight
            vnlspd(kx,ky,kz, 3)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 4)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 5)=tmpspd3*tmp*dpxyz*sq3fourpi*ylm04*weight
            vnlspd(kx,ky,kz, 6)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm02*weight
            vnlspd(kx,ky,kz, 7)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm03*weight
            vnlspd(kx,ky,kz, 8)=tmpspd4*tmp*dpxyz*sq3fourpi*ylm04*weight
            vnlspd(kx,ky,kz, 9)=tmpspd5*tmp*dpxyz*sq1516pi*ylm05*weight
            vnlspd(kx,ky,kz,10)=tmpspd5*tmp*dpxyz*sq1504pi*ylm06*weight
            vnlspd(kx,ky,kz,11)=tmpspd5*tmp*dpxyz*sq0516pi*ylm07*weight
            vnlspd(kx,ky,kz,12)=tmpspd5*tmp*dpxyz*sq1504pi*ylm08*weight
            vnlspd(kx,ky,kz,13)=tmpspd5*tmp*dpxyz*sq1504pi*ylm09*weight
            vnlspd(kx,ky,kz,14)=tmpspd6*tmp*dpxyz*sq1516pi*ylm05*weight
            vnlspd(kx,ky,kz,15)=tmpspd6*tmp*dpxyz*sq1504pi*ylm06*weight
            vnlspd(kx,ky,kz,16)=tmpspd6*tmp*dpxyz*sq0516pi*ylm07*weight
            vnlspd(kx,ky,kz,17)=tmpspd6*tmp*dpxyz*sq1504pi*ylm08*weight
            vnlspd(kx,ky,kz,18)=tmpspd6*tmp*dpxyz*sq1504pi*ylm09*weight
          end if
        end do
        end do
        end do
      end if
! ===================================

      js=-nfdg*npmesh
      je=nfdg*npmesh-1
      if (mod(npmesh,2) .eq. 1) js=-nfdg*npmesh+1
      do la=1,nprj(indspe(na))
!$omp do
        do ix=1,8*npxmax*npymax*npzmax
          vtmp_pp(-npxmax+ix,-npymax+1,-npzmax+1)=0.0d0
        end do
!$omp do
        do ix=1,8*nxspd*nyspd*npzmax
          vtmp_pp1(-nxspd+ix,-nyspd+1,-npzmax+1)=0.0d0
        end do
!$omp do
        do ix=1,8*nxspd*npymax*npzmax
          vtmp_pp2(-nxspd+ix,-npymax+1,-npzmax+1)=0.0d0
        end do
!$omp do
        do iz=-npzmax+1,npzmax
        do iy=-nyspd+1,nyspd
        do ix=-nxspd+1,nxspd
          kkz=iz*npmesh-(npmesh-1)/2
          do jz=js,je
            kz=kkz+jz
            vtmp_pp1(ix,iy,iz)=vtmp_pp1(ix,iy,iz)+vnlspd(ix,iy,kz,la)*wxyz(jz)
          end do
        end do
        end do
        end do
!$omp do
        do iz=-npzmax+1,npzmax
        do iy=-npymax+1,npymax
        do ix=-nxspd+1,nxspd
          kky=iy*npmesh-(npmesh-1)/2
          do jy=js,je
            ky=kky+jy
            vtmp_pp2(ix,iy,iz)=vtmp_pp2(ix,iy,iz)+vtmp_pp1(ix,ky,iz)*wxyz(jy)
          end do
        end do
        end do
        end do
!$omp do
        do iz=-npzmax+1,npzmax
        do iy=-npymax+1,npymax
        do ix=-npxmax+1,npxmax
          kkx=ix*npmesh-(npmesh-1)/2
          do jx=js,je
            kx=kkx+jx
            vtmp_pp(ix,iy,iz)=vtmp_pp(ix,iy,iz)+vtmp_pp2(kx,iy,iz)*wxyz(jx)
          end do
        end do
        end do
        end do
!$omp do
        do ixyz=1,natinf(na)
          ix=lstx(ixyz,iaps)
          iy=lsty(ixyz,iaps)
          iz=lstz(ixyz,iaps)
          vnlocp(ixyz,la,iaps)=vtmp_pp(ix,iy,iz)
        end do
      end do
! =======================================================
    end if
  end do

  return
end subroutine pseudocalc_07


!this subroutine computes local parts of pseudopotential on radial grid.
subroutine pseudocalc_08a(natom,num_spe,nradmx,npoint,num_atcell, & ! <
                          key_natpri_in,key_pp_paw,               & ! <
                          indspe,natpri,natpri_inf,ntyppp,nradct, & ! <
                          cp,radial,point,                        & ! <
                          atx,aty,atz,                            & ! <
                          vcorer)                                   ! >
implicit none
integer, intent(in)::natom,num_spe,nradmx,npoint,num_atcell
integer, intent(in)::key_natpri_in,key_pp_paw
integer, intent(in)::indspe(natom),natpri(natom),natpri_inf(natom),ntyppp(num_spe),nradct(num_spe)
real*8, intent(in)::cp(8,num_spe),radial(nradmx,num_spe)
real*8, intent(in)::point(npoint,3)
real*8, intent(in)::atx(natom),aty(natom),atz(natom)
real*8, intent(out)::vcorer(nradmx,npoint,num_atcell)
real*8 x,y,z,r,tmp,cp1
integer na0,na1,ipri,il,ir

!$omp do
  do ir=1,nradmx*npoint*num_atcell
    vcorer(ir,1,1)=0.0d0
  end do

  do na0=1,natom
    if ((ntyppp(indspe(na0)) .eq. key_pp_paw) .and. (natpri(na0) .eq. key_natpri_in)) then
      ipri=natpri_inf(na0)
      do na1=1,natom
! vcorer does not include Coulomb potential of the owner of augmented sphere.
        if (na1 .ne. na0) then
          cp1=cp(1,indspe(na1))
!$omp do
          do il=1,npoint
            do ir=2,nradct(indspe(na0))
              x=point(il,1)*radial(ir,indspe(na0))+atx(na0)
              y=point(il,2)*radial(ir,indspe(na0))+aty(na0)
              z=point(il,3)*radial(ir,indspe(na0))+atz(na0)
              r=dsqrt((x-atx(na1))*(x-atx(na1))+(y-aty(na1))*(y-aty(na1))+(z-atz(na1))*(z-atz(na1)))
              tmp=-cp1/r
              vcorer(ir,il,ipri)=vcorer(ir,il,ipri)+tmp
            end do
          end do
!$omp end do nowait
        end if
      end do
    end if
  end do
!$omp barrier

  return
end subroutine pseudocalc_08a


!this subroutine computes local parts of pseudopotential on the radial grid of the augmented sphere.
subroutine pseudocalc_08b(nperi,natom,num_spe,nradmx,npoint,lsphel,num_atcell,nint1dmax, & ! <
                          new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz,    & ! <
                          key_natpri_in,key_pp_paw,                           & ! <
                          xmax,ymax,zmax,veta,                                & ! <
                          indspe,natpri,natpri_inf,ntyppp,nradct,             & ! <
                          cp,radial,point,                                    & ! <
                          atx,aty,atz,                                        & ! <
                          vcorer)                                               ! >
use mod_mpi
implicit none
integer, intent(in)::nperi,natom,num_spe,nradmx,npoint,lsphel,num_atcell,nint1dmax
integer, intent(in)::new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz
integer, intent(in)::key_natpri_in,key_pp_paw
integer, intent(in)::indspe(natom),natpri(natom),natpri_inf(natom),ntyppp(num_spe),nradct(num_spe)
real*8, intent(in)::xmax,ymax,zmax,veta
real*8, intent(in)::cp(8,num_spe),radial(nradmx,num_spe)
real*8, intent(in)::point(npoint,3)
real*8, intent(in)::atx(natom),aty(natom),atz(natom)
real*8, intent(out)::vcorer(nradmx,npoint,num_atcell)
real*8, allocatable::anum0(:,:,:)
real*8 derf
real*8 pi,pi05,omega,omegain,surf,surfin,avcharge,x,y,z,x0,y0,z0,z1,r,rcut,cp1, &
       tmp,theta,rin,r2in,t,dt,dr
integer, parameter ::nt=200,nr=200

integer na,na0,na1,ipri,il,ir,kx,ky,kz,kpx,kpy,kpz,i,j,mcirc

  if (nperi==1) allocate(anum0(-new_pwx:new_pwx,nradmx,npoint))
  pi=dacos(-1.0d0)
  pi05=dsqrt(pi)

!$omp do
  do ir=1,nradmx*npoint*num_atcell
    vcorer(ir,1,1)=0.0d0
  end do

  select case (nperi)
  case(1)
    do na0=1,natom
      if ((ntyppp(indspe(na0)) .eq. key_pp_paw) .and. (natpri(na0) .eq. key_natpri_in)) then
        ipri=natpri_inf(na0)
        do na1=1,natom
          cp1=cp(1,indspe(na1))
          y0=aty(na1)-aty(na0)
          z0=atz(na1)-atz(na0)
          if (y0*y0+z0*z0 .lt. radial(nradct(indspe(na0)),indspe(na0))*radial(nradct(indspe(na0)),indspe(na0))) then
            dr=dsqrt(radial(nradct(indspe(na0)),indspe(na0))**2-y0*y0-z0*z0)/nr
!$omp do
            do il=1,npoint
              do ir=2,nradct(indspe(na0))
                x=point(il,1)*radial(ir,indspe(na0))
                y=point(il,2)*radial(ir,indspe(na0))
                z=point(il,3)*radial(ir,indspe(na0))
                do i=-nr,nr
                  x0=i*dr-0.5d0*dr
                  r=dsqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0))
                  vcorer(ir,il,ipri)=vcorer(ir,il,ipri)+cp1/(2.0d0*xmax)*dr/r
                end do
              end do
            end do
          end if
        end do !na1
      end if
    end do !na0
  case(2)
    surf=xmax*ymax*4.0d0
    surfin=1.0d0/surf
    mcirc=(lsphel+1)/2
    do na0=1,natom
      if ((ntyppp(indspe(na0)) .eq. key_pp_paw) .and. (natpri(na0) .eq. key_natpri_in)) then
        ipri=natpri_inf(na0)
        do na1=1,natom
          cp1=cp(1,indspe(na1))
          z0=atz(na1)-atz(na0)
          if (z0*z0 .lt. radial(nradct(indspe(na0)),indspe(na0))*radial(nradct(indspe(na0)),indspe(na0))) then
            dt=pi/nt
            dr=dsqrt(radial(nradct(indspe(na0)),indspe(na0))*radial(nradct(indspe(na0)),indspe(na0))-z0*z0)/nr
!$omp do
            do ir=2,nradct(indspe(na0))
              do il=1,npoint/(2*mcirc)
                z=point(il,3)*radial(ir,indspe(na0))
                tmp=0.0d0
                x=point(il,1)*radial(ir,indspe(na0))
                y=point(il,2)*radial(ir,indspe(na0))
                z=point(il,3)*radial(ir,indspe(na0))
                theta=dimag(log(dcmplx(x,y)))
                do j=1,nt
                  do i=1,nr
                    x0=(i*dr)*dcos(theta+j*dt-0.5d0*dt)
                    y0=(i*dr)*dsin(theta+j*dt-0.5d0*dt)
                    r=dsqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0))
                    tmp=tmp+2.0d0*cp1*surfin*dr*dt*(i*dr)/r
                  end do
                end do
                do i=0,2*mcirc-1
                  vcorer(ir,il+i*npoint/(2*mcirc),ipri)=vcorer(ir,il+i*npoint/(2*mcirc),ipri)+tmp
                end do
              end do
            end do
          end if
        end do !na1
      end if
    end do !na0
  case(3)
    omega=xmax*ymax*zmax*8.0d0
    omegain=1.0d0/omega

    avcharge=0.0d0
    do na=1,natom
      cp1=cp(1,indspe(na))
      avcharge=avcharge+cp1*omegain
    end do

    do na=1,natom
      if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
        ipri=natpri_inf(na)
        cp1=cp(1,indspe(na))
!$omp do
        do il=1,npoint
        do ir=2,nradct(indspe(na))
          r=radial(ir,indspe(na))
          rcut=radial(nradct(indspe(na)),indspe(na))
! vcorer does not include Coulomb potential of the owner of augmented sphere.
          vcorer(ir,il,ipri)=avcharge/6.0d0*(4.0d0*pi)*(rcut**2-r**2)+cp1/rcut
        end do
        end do
!$omp end do nowait
      end if
    end do
!$omp barrier
  end select
  if (nperi==1) deallocate(anum0)
  return
end subroutine pseudocalc_08b


!this subroutine computes the potential from the external electric field.
subroutine pseudocalc_11(natom,num_spe,nradmx,npoint,num_atcell,nzmax,jelcalc, & ! <
                         key_natpri_in,key_pp_paw,key_jel_calc,                & ! <
                         indspe,ntyppp,nradct,natpri,natpri_inf,               & ! <
                         xmax,ymax,zmax,                                       & ! <
                         biasx,biasy,biasz,                                    & ! <
                         chrjel,strjel,endjel,                                 & ! <
                         point,radial,potc,                                    & ! <
                         vcorer,vcorer_all,                                    & ! X
                         atx,aty,atz)                                            ! <
implicit none
integer,intent(in)::natom,num_spe,nradmx,npoint,num_atcell,nzmax,jelcalc
integer,intent(in)::key_natpri_in,key_pp_paw,key_jel_calc
integer,intent(in)::indspe(natom),ntyppp(num_spe),nradct(num_spe),natpri(natom),natpri_inf(natom)
real*8, intent(in)::xmax,ymax,zmax
real*8, intent(in)::biasx,biasy,biasz
real*8, intent(in)::chrjel,strjel,endjel
real*8, intent(in)::point(npoint,3),radial(nradmx,num_spe),potc(nradmx,num_spe)
real*8, intent(inout)::vcorer(nradmx,npoint,num_atcell)
real*8, intent(out)::vcorer_all(nradmx,npoint,num_atcell)
real*8, intent(in)::atx(natom),aty(natom),atz(natom)
real*8 x,y,z,omega,vkz,vk2,pi
integer na,ipri,il,ir,kz
  pi=dacos(-1.0d0)

  omega=8.0d0*xmax*ymax*zmax
  if (jelcalc==key_jel_calc) then
    do na=1,natom
      if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
        ipri=natpri_inf(na)
        do kz=-nzmax,nzmax
          if (kz**2 .ne. 0) then
            vkz=pi/zmax*kz
            vk2=vkz*vkz
!$omp do
            do il=1,npoint
              do ir=2,nradct(indspe(na))
                z=point(il,3)*radial(ir,indspe(na))+atz(na)
                vcorer(ir,il,ipri)=vcorer(ir,il,ipri) &
                 -4.0d0*pi/omega*chrjel/(strjel-endjel)/vk2*(dsin(vkz*(strjel-z))-dsin(vkz*(endjel-z)))/vkz
              end do
            end do
!$omp end do nowait
          end if
        end do
      end if
    end do
!$omp barrier
  end if

  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
      ipri=natpri_inf(na)
!$omp do
      do il=1,npoint
        do ir=2,nradct(indspe(na))
          x=point(il,1)*radial(ir,indspe(na))+atx(na)
          y=point(il,2)*radial(ir,indspe(na))+aty(na)
          z=point(il,3)*radial(ir,indspe(na))+atz(na)
! vcorer does not include Coulomb potential of the owner of augmented sphere.
! vcorer_all includes Coulomb potential of all atoms.
          vcorer(ir,il,ipri)=vcorer(ir,il,ipri)+biasx*x+biasy*y+biasz*z
          vcorer_all(ir,il,ipri)=vcorer(ir,il,ipri)+potc(ir,indspe(na))
        end do
      end do
!$omp end do nowait
    end if
  end do
!$omp barrier

  return
end subroutine pseudocalc_11


!this subroutine computes moments of charges, \int Qij |r|^l Y(\hat{r}) dr in Eq. (26) of PRB59 1758 (1999).
subroutine pseudocalc_12(natom,nradmx,num_spe,nprmx,lmx,npoint,lrhomx, & ! <
                         na,la,lla,i1,i2,j1,j2,                        & ! <
                         tmp,sum,                                      & ! X
                         indspe,nradct,                                & ! <
                         radial,dradial,awf,pwf,yylm,wt)                 ! <
implicit none
integer,intent(in)   ::natom,nradmx,num_spe,nprmx,lmx
integer,intent(in)   ::na,la,lla,i1,i2,j1,j2
integer,intent(in)   ::npoint,lrhomx
integer,intent(in)   ::indspe(natom),nradct(num_spe)
real*8, intent(in)   ::radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8, intent(in)   ::awf(nradmx,nprmx*lmx,num_spe),pwf(nradmx,nprmx*lmx,num_spe)
real*8, intent(in)   ::yylm(npoint,lrhomx),wt(npoint)
real*8, intent(inout)::tmp,sum
integer il,ir
real*8  tmp0,sum0,r,dr

  tmp0=0.0d0
!$omp do
  do il=1,npoint
    tmp0=tmp0+yylm(il,i1)*yylm(il,j1)*yylm(il,la)*wt(il)
  end do
!$omp end do nowait
  sum0=0.0d0
!$omp do
  do ir=2,nradct(indspe(na))
    r=radial(ir,indspe(na))
    dr=dradial(ir,indspe(na))
    sum0=sum0+(awf(ir,i2,indspe(na))*awf(ir,j2,indspe(na)) &
              -pwf(ir,i2,indspe(na))*pwf(ir,j2,indspe(na)))*dr*r**lla
  end do
!$omp end do nowait
!$omp critical
  tmp=tmp+tmp0
  sum=sum+sum0
!$omp end critical

  return
end subroutine pseudocalc_12


subroutine pseudocalc_13( &
 natom,num_spe,nperi,nfh,nint1dmax,new_pwx,new_pwy,new_rsx,new_rsy,ncpx_d,ncpy_d,ncpz_d, & ! <
 veta,ddx,ddy,ddz,xmax,ymax,zmax,                                                        & ! <
 indspe,atx,aty,atz,                                                                     & ! <
 vboundx,vboundy,vboundz)                                                                  ! >
use mod_mpi,       only: myrx,myry,myrz,nprocx,nprocy,nprocz
use mod_mathfunctions,only: expint1
implicit none
integer, intent(in)::natom,num_spe,nperi,nfh,nint1dmax,new_pwx,new_pwy,new_rsx,new_rsy
integer, intent(in)::ncpx_d,ncpy_d,ncpz_d
integer, intent(in)::indspe(natom)
real*8, intent(in)::veta
real*8, intent(in)::ddx,ddy,ddz
real*8, intent(in)::xmax,ymax,zmax
real*8, intent(in)::atx(natom),aty(natom),atz(natom)
real*8, intent(out)::vboundx(-(nfh-1):nfh,ncpy_d,ncpz_d,9*(3-nperi)/3+1,(natom-1)*(3-nperi)/3+1)
real*8, intent(out)::vboundy(ncpx_d,-(nfh-1):nfh,ncpz_d,9*(3-nperi)/2+1,(natom-1)*(3-nperi)/2+1)
real*8, intent(out)::vboundz(ncpx_d,ncpy_d,-(nfh-1):nfh,9*(1-nperi/3)+1,(natom-1)*(1-nperi/3)+1)
real*8, allocatable::anumy(:,:,:,:),anumz(:,:,:,:)
integer na,ix,iy,iz,jx,jy,jz,kx,ky,i,ierr
real*8 x,y,z,r,rin,r2in,r3in,r5in,derfdrr,pi,surfin,pisurf,vkx,vky,vk2,ta,tb,deftp1,deftp2,pi05,t,dt
real*8 vep0,vep1,vep2,vlo,vlo0,vlo1,vlo2,vlo0r,vlo1r,vlo2r,vlo2rrr,vlo12r,vlorr,vsine,vcosi,zab,derftmp,dexptmp,deftmp0,dextmp0
real*8 derf

  if (nperi==1) allocate(anumy(-new_pwx:new_pwx,-(nfh-1):nfh,ncpz_d,3),anumz(-new_pwx:new_pwx,ncpy_d,-(nfh-1):nfh,3))
  pi=dacos(-1.0d0)
  pi05=dsqrt(pi)
  surfin=1.0d0/(4.0d0*xmax*ymax)

!$omp single
  vboundx(:,:,:,:,:)=0.0d0
  vboundy(:,:,:,:,:)=0.0d0
  vboundz(:,:,:,:,:)=0.0d0
!$omp end single
!$omp barrier

  select case (nperi)
  case (0)
    do na=1,natom
    if (myrx==0) then
!$omp do
      do iz=1,ncpz_d
      do iy=1,ncpy_d
      do ix=-(nfh-1),0
        jx=myrx*ncpx_d+ix
        jy=myry*ncpy_d+iy
        jz=myrz*ncpz_d+iz
        x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
        y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
        z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
        r=dsqrt(x*x+y*y+z*z)
        rin=1.0d0/r
        r2in=rin*rin
        r3in=r2in*rin
        r5in=r3in*r2in
        vboundx(ix,iy,iz, 1,na)=rin
        derfdrr=1.0d0*r3in
        vboundx(ix,iy,iz, 2,na)=derfdrr*x
        vboundx(ix,iy,iz, 3,na)=derfdrr*y
        vboundx(ix,iy,iz, 4,na)=derfdrr*z
        vboundx(ix,iy,iz, 5,na)=(3.0d0*x*x*r5in-r3in)
        vboundx(ix,iy,iz, 6,na)=(3.0d0*y*y*r5in-r3in)
        vboundx(ix,iy,iz, 7,na)=(3.0d0*z*z*r5in-r3in)
        vboundx(ix,iy,iz, 8,na)=3.0d0*x*y*r5in
        vboundx(ix,iy,iz, 9,na)=3.0d0*y*z*r5in
        vboundx(ix,iy,iz,10,na)=3.0d0*z*x*r5in
      end do
      end do
      end do
    end if
    if (myrx==nprocx-1) then
!$omp do
      do iz=1,ncpz_d
      do iy=1,ncpy_d
      do ix=1,nfh
        jx=myrx*ncpx_d+ncpx_d+ix
        jy=myry*ncpy_d+iy
        jz=myrz*ncpz_d+iz
        x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
        y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
        z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
        r=dsqrt(x*x+y*y+z*z)
        rin=1.0d0/r
        r2in=rin*rin
        r3in=r2in*rin
        r5in=r3in*r2in
        vboundx(ix,iy,iz, 1,na)=rin
        derfdrr=1.0d0*r3in
        vboundx(ix,iy,iz, 2,na)=derfdrr*x
        vboundx(ix,iy,iz, 3,na)=derfdrr*y
        vboundx(ix,iy,iz, 4,na)=derfdrr*z
        vboundx(ix,iy,iz, 5,na)=(3.0d0*x*x*r5in-r3in)
        vboundx(ix,iy,iz, 6,na)=(3.0d0*y*y*r5in-r3in)
        vboundx(ix,iy,iz, 7,na)=(3.0d0*z*z*r5in-r3in)
        vboundx(ix,iy,iz, 8,na)=3.0d0*x*y*r5in
        vboundx(ix,iy,iz, 9,na)=3.0d0*y*z*r5in
        vboundx(ix,iy,iz,10,na)=3.0d0*z*x*r5in
      end do
      end do
      end do
    end if
    if (myry==0) then
!$omp do
      do iz=1,ncpz_d
      do iy=-(nfh-1),0
      do ix=1,ncpx_d
        jx=myrx*ncpx_d+ix
        jy=myry*ncpy_d+iy
        jz=myrz*ncpz_d+iz
        x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
        y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
        z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
        r=dsqrt(x*x+y*y+z*z)
        rin=1.0d0/r
        r2in=rin*rin
        r3in=r2in*rin
        r5in=r3in*r2in
        vboundy(ix,iy,iz, 1,na)=rin
        derfdrr=1.0d0*r3in
        vboundy(ix,iy,iz, 2,na)=derfdrr*x
        vboundy(ix,iy,iz, 3,na)=derfdrr*y
        vboundy(ix,iy,iz, 4,na)=derfdrr*z
        vboundy(ix,iy,iz, 5,na)=(3.0d0*x*x*r5in-r3in)
        vboundy(ix,iy,iz, 6,na)=(3.0d0*y*y*r5in-r3in)
        vboundy(ix,iy,iz, 7,na)=(3.0d0*z*z*r5in-r3in)
        vboundy(ix,iy,iz, 8,na)=3.0d0*x*y*r5in
        vboundy(ix,iy,iz, 9,na)=3.0d0*y*z*r5in
        vboundy(ix,iy,iz,10,na)=3.0d0*z*x*r5in
      end do
      end do
      end do
    end if
    if (myry==nprocy-1) then
!$omp do
      do iz=1,ncpz_d
      do iy=1,nfh
      do ix=1,ncpx_d
        jx=myrx*ncpx_d+ix
        jy=myry*ncpy_d+ncpy_d+iy
        jz=myrz*ncpz_d+iz
        x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
        y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
        z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
        r=dsqrt(x*x+y*y+z*z)
        rin=1.0d0/r
        r2in=rin*rin
        r3in=r2in*rin
        r5in=r3in*r2in
        vboundy(ix,iy,iz, 1,na)=rin
        derfdrr=1.0d0*r3in
        vboundy(ix,iy,iz, 2,na)=derfdrr*x
        vboundy(ix,iy,iz, 3,na)=derfdrr*y
        vboundy(ix,iy,iz, 4,na)=derfdrr*z
        vboundy(ix,iy,iz, 5,na)=(3.0d0*x*x*r5in-r3in)
        vboundy(ix,iy,iz, 6,na)=(3.0d0*y*y*r5in-r3in)
        vboundy(ix,iy,iz, 7,na)=(3.0d0*z*z*r5in-r3in)
        vboundy(ix,iy,iz, 8,na)=3.0d0*x*y*r5in
        vboundy(ix,iy,iz, 9,na)=3.0d0*y*z*r5in
        vboundy(ix,iy,iz,10,na)=3.0d0*z*x*r5in
      end do
      end do
      end do
    end if
    if (myrz==0) then
!$omp do
      do iz=-(nfh-1),0
      do iy=1,ncpy_d
      do ix=1,ncpx_d
        jx=myrx*ncpx_d+ix
        jy=myry*ncpy_d+iy
        jz=myrz*ncpz_d+iz
        x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
        y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
        z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
        r=dsqrt(x*x+y*y+z*z)
        rin=1.0d0/r
        r2in=rin*rin
        r3in=r2in*rin
        r5in=r3in*r2in
        vboundz(ix,iy,iz, 1,na)=rin
        derfdrr=1.0d0*r3in
        vboundz(ix,iy,iz, 2,na)=derfdrr*x
        vboundz(ix,iy,iz, 3,na)=derfdrr*y
        vboundz(ix,iy,iz, 4,na)=derfdrr*z
        vboundz(ix,iy,iz, 5,na)=(3.0d0*x*x*r5in-r3in)
        vboundz(ix,iy,iz, 6,na)=(3.0d0*y*y*r5in-r3in)
        vboundz(ix,iy,iz, 7,na)=(3.0d0*z*z*r5in-r3in)
        vboundz(ix,iy,iz, 8,na)=3.0d0*x*y*r5in
        vboundz(ix,iy,iz, 9,na)=3.0d0*y*z*r5in
        vboundz(ix,iy,iz,10,na)=3.0d0*z*x*r5in
      end do
      end do
      end do
    end if
    if (myrz==nprocz-1) then
!$omp do
      do iz=1,nfh
      do iy=1,ncpy_d
      do ix=1,ncpx_d
        jx=myrx*ncpx_d+ix
        jy=myry*ncpy_d+iy
        jz=myrz*ncpz_d+ncpz_d+iz
        x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
        y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
        z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
        r=dsqrt(x*x+y*y+z*z)
        rin=1.0d0/r
        r2in=rin*rin
        r3in=r2in*rin
        r5in=r3in*r2in
        vboundz(ix,iy,iz, 1,na)=rin
        derfdrr=1.0d0*r3in
        vboundz(ix,iy,iz, 2,na)=derfdrr*x
        vboundz(ix,iy,iz, 3,na)=derfdrr*y
        vboundz(ix,iy,iz, 4,na)=derfdrr*z
        vboundz(ix,iy,iz, 5,na)=(3.0d0*x*x*r5in-r3in)
        vboundz(ix,iy,iz, 6,na)=(3.0d0*y*y*r5in-r3in)
        vboundz(ix,iy,iz, 7,na)=(3.0d0*z*z*r5in-r3in)
        vboundz(ix,iy,iz, 8,na)=3.0d0*x*y*r5in
        vboundz(ix,iy,iz, 9,na)=3.0d0*y*z*r5in
        vboundz(ix,iy,iz,10,na)=3.0d0*z*x*r5in
      end do
      end do
      end do
    end if
    end do
  case (1)
    do na=1,natom
      if (myry==0) then
!$omp do
        do iz=1,ncpz_d
        do iy=-(nfh-1),0
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
            vkx=pi/xmax*kx
            vk2=vkx**2
            jy=myry*ncpy_d+iy
            jz=myrz*ncpz_d+iz
            y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
            z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
            ta=(y*y+z*z)
            tb=vkx*vkx
            anumy(kx,iy,iz,1)=0.0d0
            anumy(kx,iy,iz,2)=0.0d0
            anumy(kx,iy,iz,3)=0.0d0
            dt=veta/nint1dmax
            do i=1,nint1dmax
              t=i*dt
              anumy(kx,iy,iz,1)=anumy(kx,iy,iz,1)+dexp(-ta*t*t-0.25d0*tb/(t*t))/t*dt
              anumy(kx,iy,iz,2)=anumy(kx,iy,iz,2)+dexp(-ta*t*t-0.25d0*tb/(t*t))*t*dt
              anumy(kx,iy,iz,3)=anumy(kx,iy,iz,3)+dexp(-ta*t*t-0.25d0*tb/(t*t))*t*t*t*dt
            end do
          end if
        end do
        end do
        end do
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
            vkx=pi/xmax*kx
            vk2=vkx**2
!$omp do
            do iz=1,ncpz_d
            do iy=-(nfh-1),0
            do ix=1,ncpx_d
              jx=myrx*ncpx_d+ix
              jy=myry*ncpy_d+iy
              jz=myrz*ncpz_d+iz
              x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
              y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
              z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
              vcosi=dcos(vkx*x)/(xmax*2.0d0)
              vsine=dsin(vkx*x)/(xmax*2.0d0)
              vboundy(ix,iy,iz, 1,na)=vboundy(ix,iy,iz, 1,na)+2.0d0*vcosi*anumy(kx,iy,iz,1)
              vboundy(ix,iy,iz, 2,na)=vboundy(ix,iy,iz, 2,na)+2.0d0*vkx*vsine*anumy(kx,iy,iz,1)
              vboundy(ix,iy,iz, 3,na)=vboundy(ix,iy,iz, 3,na)+4.0d0*vcosi*anumy(kx,iy,iz,2)*y
              vboundy(ix,iy,iz, 4,na)=vboundy(ix,iy,iz, 4,na)+4.0d0*vcosi*anumy(kx,iy,iz,2)*z
              vboundy(ix,iy,iz, 5,na)=vboundy(ix,iy,iz, 5,na)-2.0d0*vkx*vkx*vcosi*anumy(kx,iy,iz,1)
              vboundy(ix,iy,iz, 6,na)=vboundy(ix,iy,iz, 6,na)-4.0d0*vcosi*(anumy(kx,iy,iz,2)-2.0d0*y*y*anumy(kx,iy,iz,3))
              vboundy(ix,iy,iz, 7,na)=vboundy(ix,iy,iz, 7,na)-4.0d0*vcosi*(anumy(kx,iy,iz,2)-2.0d0*z*z*anumy(kx,iy,iz,3))
              vboundy(ix,iy,iz, 8,na)=vboundy(ix,iy,iz, 8,na)+4.0d0*vkx*vsine*anumy(kx,iy,iz,2)*y
              vboundy(ix,iy,iz, 9,na)=vboundy(ix,iy,iz, 9,na)+8.0d0*vcosi*anumy(kx,iy,iz,3)*y*z
              vboundy(ix,iy,iz,10,na)=vboundy(ix,iy,iz,10,na)+4.0d0*vkx*vsine*anumy(kx,iy,iz,2)*z
            end do
            end do
            end do
          end if
        end do
!$omp do
        do iz=1,ncpz_d
        do iy=-(nfh-1),0
        do ix=1,ncpx_d
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
          y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
          z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
          ta=(y*y+z*z)
          vlo=(-0.5d0*expint1(ta*veta*veta,ierr)-0.5d0*dlog(ta))/xmax
          vlo0r=2.0d0*(1.0d0-dexp(-veta*veta*ta))/(ta*2.0d0*xmax)
          vlorr=-2.0d0*(1.0d0-dexp(-veta*veta*ta))/(ta*ta*2.0d0*xmax) &
                +2.0d0*veta*veta*dexp(-veta*veta*ta)/(ta*2.0d0*xmax)
          vboundy(ix,iy,iz, 1,na)=vboundy(ix,iy,iz, 1,na)+vlo
          vboundy(ix,iy,iz, 3,na)=vboundy(ix,iy,iz, 3,na)+vlo0r*y
          vboundy(ix,iy,iz, 4,na)=vboundy(ix,iy,iz, 4,na)+vlo0r*z
          vboundy(ix,iy,iz, 6,na)=vboundy(ix,iy,iz, 6,na)-(vlo0r+2.0d0*y*y*vlorr)
          vboundy(ix,iy,iz, 7,na)=vboundy(ix,iy,iz, 7,na)-(vlo0r+2.0d0*z*z*vlorr)
          vboundy(ix,iy,iz, 9,na)=vboundy(ix,iy,iz, 9,na)-2.0d0*y*z*vlorr
        end do
        end do
        end do
        do kx=-new_rsx,new_rsx
!$omp do
          do iz=1,ncpz_d
          do iy=-(nfh-1),0
          do ix=1,ncpx_d
            jx=myrx*ncpx_d+ix
            jy=myry*ncpy_d+iy
            jz=myrz*ncpz_d+iz
            x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)+kx*xmax*2.0d0
            y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
            z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
            r=dsqrt(x*x+y*y+z*z)
            rin=1.0d0/r
            r2in=rin*rin
            r3in=r2in*rin
            derftmp=(1.0d0-derf(veta*r))
            dexptmp=dexp(-veta*veta*r*r)/pi05
            vlo0=-2.0d0*veta*dexptmp*rin-derftmp*r2in
            vlo1=4.0d0*veta*dexptmp*(veta*veta+r2in)+2.0d0*derftmp*r3in
            vlo2=-vlo0
            vlo0r=vlo0*rin
            vlo1r=vlo1*r2in
            vlo2r=vlo2*r3in
            vlo2rrr=vlo2*rin
            vlo12r=(vlo1r+vlo2r)
            vboundy(ix,iy,iz, 1,na)=vboundy(ix,iy,iz, 1,na)+rin*(-derf(veta*r)+1.0d0)
            vboundy(ix,iy,iz, 2,na)=vboundy(ix,iy,iz, 2,na)-vlo0r*x
            vboundy(ix,iy,iz, 3,na)=vboundy(ix,iy,iz, 3,na)-vlo0r*y
            vboundy(ix,iy,iz, 4,na)=vboundy(ix,iy,iz, 4,na)-vlo0r*z
            vboundy(ix,iy,iz, 5,na)=vboundy(ix,iy,iz, 5,na)+(x*x*vlo12r-vlo2rrr)
            vboundy(ix,iy,iz, 6,na)=vboundy(ix,iy,iz, 6,na)+(y*y*vlo12r-vlo2rrr)
            vboundy(ix,iy,iz, 7,na)=vboundy(ix,iy,iz, 7,na)+(z*z*vlo12r-vlo2rrr)
            vboundy(ix,iy,iz, 8,na)=vboundy(ix,iy,iz, 8,na)+x*y*vlo12r
            vboundy(ix,iy,iz, 9,na)=vboundy(ix,iy,iz, 9,na)+y*z*vlo12r
            vboundy(ix,iy,iz,10,na)=vboundy(ix,iy,iz,10,na)+z*x*vlo12r
          end do
          end do
          end do
        end do
      end if
      if (myry==nprocy-1) then
!$omp do
        do iz=1,ncpz_d
        do iy=1,nfh
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
            vkx=pi/xmax*kx
            vk2=vkx**2
            jy=myry*ncpy_d+ncpy_d+iy
            jz=myrz*ncpz_d+iz
            y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
            z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
            ta=(y*y+z*z)
            tb=vkx*vkx
            anumy(kx,iy,iz,1)=0.0d0
            anumy(kx,iy,iz,2)=0.0d0
            anumy(kx,iy,iz,3)=0.0d0
            dt=veta/nint1dmax
            do i=1,nint1dmax
              t=i*dt
              anumy(kx,iy,iz,1)=anumy(kx,iy,iz,1)+dexp(-ta*t*t-0.25d0*tb/(t*t))/t*dt
              anumy(kx,iy,iz,2)=anumy(kx,iy,iz,2)+dexp(-ta*t*t-0.25d0*tb/(t*t))*t*dt
              anumy(kx,iy,iz,3)=anumy(kx,iy,iz,3)+dexp(-ta*t*t-0.25d0*tb/(t*t))*t*t*t*dt
            end do
          end if
        end do
        end do
        end do
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
            vkx=pi/xmax*kx
            vk2=vkx**2
!$omp do
            do iz=1,ncpz_d
            do iy=1,nfh
            do ix=1,ncpx_d
              jx=myrx*ncpx_d+ix
              jy=myry*ncpy_d+ncpy_d+iy
              jz=myrz*ncpz_d+iz
              x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
              y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
              z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
              vcosi=dcos(vkx*x)/(xmax*2.0d0)
              vsine=dsin(vkx*x)/(xmax*2.0d0)
              vboundy(ix,iy,iz, 1,na)=vboundy(ix,iy,iz, 1,na)+2.0d0*vcosi*anumy(kx,iy,iz,1)
              vboundy(ix,iy,iz, 2,na)=vboundy(ix,iy,iz, 2,na)+2.0d0*vkx*vsine*anumy(kx,iy,iz,1)
              vboundy(ix,iy,iz, 3,na)=vboundy(ix,iy,iz, 3,na)+4.0d0*vcosi*anumy(kx,iy,iz,2)*y
              vboundy(ix,iy,iz, 4,na)=vboundy(ix,iy,iz, 4,na)+4.0d0*vcosi*anumy(kx,iy,iz,2)*z
              vboundy(ix,iy,iz, 5,na)=vboundy(ix,iy,iz, 5,na)-2.0d0*vkx*vkx*vcosi*anumy(kx,iy,iz,1)
              vboundy(ix,iy,iz, 6,na)=vboundy(ix,iy,iz, 6,na)-4.0d0*vcosi*(anumy(kx,iy,iz,2)-2.0d0*y*y*anumy(kx,iy,iz,3))
              vboundy(ix,iy,iz, 7,na)=vboundy(ix,iy,iz, 7,na)-4.0d0*vcosi*(anumy(kx,iy,iz,2)-2.0d0*z*z*anumy(kx,iy,iz,3))
              vboundy(ix,iy,iz, 8,na)=vboundy(ix,iy,iz, 8,na)+4.0d0*vkx*vsine*anumy(kx,iy,iz,2)*y
              vboundy(ix,iy,iz, 9,na)=vboundy(ix,iy,iz, 9,na)+8.0d0*vcosi*anumy(kx,iy,iz,3)*y*z
              vboundy(ix,iy,iz,10,na)=vboundy(ix,iy,iz,10,na)+4.0d0*vkx*vsine*anumy(kx,iy,iz,2)*z
            end do
            end do
            end do
          end if
        end do
!$omp do
        do iz=1,ncpz_d
        do iy=1,nfh
        do ix=1,ncpx_d
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
          y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
          z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
          ta=(y*y+z*z)
          vlo=(-0.5d0*expint1(ta*veta*veta,ierr)-0.5d0*dlog(ta))/xmax
          vlo0r=2.0d0*(1.0d0-dexp(-veta*veta*ta))/(ta*2.0d0*xmax)
          vlorr=-2.0d0*(1.0d0-dexp(-veta*veta*ta))/(ta*ta*2.0d0*xmax) &
                +2.0d0*veta*veta*dexp(-veta*veta*ta)/(ta*2.0d0*xmax)
          vboundy(ix,iy,iz, 1,na)=vboundy(ix,iy,iz, 1,na)+vlo
          vboundy(ix,iy,iz, 3,na)=vboundy(ix,iy,iz, 3,na)+vlo0r*y
          vboundy(ix,iy,iz, 4,na)=vboundy(ix,iy,iz, 4,na)+vlo0r*z
          vboundy(ix,iy,iz, 6,na)=vboundy(ix,iy,iz, 6,na)-(vlo0r+2.0d0*y*y*vlorr)
          vboundy(ix,iy,iz, 7,na)=vboundy(ix,iy,iz, 7,na)-(vlo0r+2.0d0*z*z*vlorr)
          vboundy(ix,iy,iz, 9,na)=vboundy(ix,iy,iz, 9,na)-2.0d0*y*z*vlorr
        end do
        end do
        end do
        do kx=-new_rsx,new_rsx
!$omp do
          do iz=1,ncpz_d
          do iy=1,nfh
          do ix=1,ncpx_d
            jx=myrx*ncpx_d+ix
            jy=myry*ncpy_d+ncpy_d+iy
            jz=myrz*ncpz_d+iz
            x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)+kx*xmax*2.0d0
            y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
            z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
            r=dsqrt(x*x+y*y+z*z)
            rin=1.0d0/r
            r2in=rin*rin
            r3in=r2in*rin
            derftmp=(1.0d0-derf(veta*r))
            dexptmp=dexp(-veta*veta*r*r)/pi05
            vlo0=-2.0d0*veta*dexptmp*rin-derftmp*r2in
            vlo1=4.0d0*veta*dexptmp*(veta*veta+r2in)+2.0d0*derftmp*r3in
            vlo2=-vlo0
            vlo0r=vlo0*rin
            vlo1r=vlo1*r2in
            vlo2r=vlo2*r3in
            vlo2rrr=vlo2*rin
            vlo12r=(vlo1r+vlo2r)
            vboundy(ix,iy,iz, 1,na)=vboundy(ix,iy,iz, 1,na)+rin*(-derf(veta*r)+1.0d0)
            vboundy(ix,iy,iz, 2,na)=vboundy(ix,iy,iz, 2,na)-vlo0r*x
            vboundy(ix,iy,iz, 3,na)=vboundy(ix,iy,iz, 3,na)-vlo0r*y
            vboundy(ix,iy,iz, 4,na)=vboundy(ix,iy,iz, 4,na)-vlo0r*z
            vboundy(ix,iy,iz, 5,na)=vboundy(ix,iy,iz, 5,na)+(x*x*vlo12r-vlo2rrr)
            vboundy(ix,iy,iz, 6,na)=vboundy(ix,iy,iz, 6,na)+(y*y*vlo12r-vlo2rrr)
            vboundy(ix,iy,iz, 7,na)=vboundy(ix,iy,iz, 7,na)+(z*z*vlo12r-vlo2rrr)
            vboundy(ix,iy,iz, 8,na)=vboundy(ix,iy,iz, 8,na)+x*y*vlo12r
            vboundy(ix,iy,iz, 9,na)=vboundy(ix,iy,iz, 9,na)+y*z*vlo12r
            vboundy(ix,iy,iz,10,na)=vboundy(ix,iy,iz,10,na)+z*x*vlo12r
          end do
          end do
          end do
        end do
      end if
      if (myrz==0) then
!$omp do
        do iz=-(nfh-1),0
        do iy=1,ncpy_d
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
            vkx=pi/xmax*kx
            vk2=vkx**2
            jy=myry*ncpy_d+iy
            jz=myrz*ncpz_d+iz
            y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
            z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
            ta=(y*y+z*z)
            tb=vkx*vkx
            anumz(kx,iy,iz,1)=0.0d0
            anumz(kx,iy,iz,2)=0.0d0
            anumz(kx,iy,iz,3)=0.0d0
            dt=veta/nint1dmax
            do i=1,nint1dmax
              t=i*dt
              anumz(kx,iy,iz,1)=anumz(kx,iy,iz,1)+dexp(-ta*t*t-0.25d0*tb/(t*t))/t*dt
              anumz(kx,iy,iz,2)=anumz(kx,iy,iz,2)+dexp(-ta*t*t-0.25d0*tb/(t*t))*t*dt
              anumz(kx,iy,iz,3)=anumz(kx,iy,iz,3)+dexp(-ta*t*t-0.25d0*tb/(t*t))*t*t*t*dt
            end do
          end if
        end do
        end do
        end do
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
            vkx=pi/xmax*kx
            vk2=vkx**2
!$omp do
            do iz=-(nfh-1),0
            do iy=1,ncpy_d
            do ix=1,ncpx_d
              jx=myrx*ncpx_d+ix
              jy=myry*ncpy_d+iy
              jz=myrz*ncpz_d+iz
              x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
              y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
              z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
              vcosi=dcos(vkx*x)/(xmax*2.0d0)
              vsine=dsin(vkx*x)/(xmax*2.0d0)
              vboundz(ix,iy,iz, 1,na)=vboundz(ix,iy,iz, 1,na)+2.0d0*vcosi*anumz(kx,iy,iz,1)
              vboundz(ix,iy,iz, 2,na)=vboundz(ix,iy,iz, 2,na)+2.0d0*vkx*vsine*anumz(kx,iy,iz,1)
              vboundz(ix,iy,iz, 3,na)=vboundz(ix,iy,iz, 3,na)+4.0d0*vcosi*anumz(kx,iy,iz,2)*y
              vboundz(ix,iy,iz, 4,na)=vboundz(ix,iy,iz, 4,na)+4.0d0*vcosi*anumz(kx,iy,iz,2)*z
              vboundz(ix,iy,iz, 5,na)=vboundz(ix,iy,iz, 5,na)-2.0d0*vkx*vkx*vcosi*anumz(kx,iy,iz,1)
              vboundz(ix,iy,iz, 6,na)=vboundz(ix,iy,iz, 6,na)-4.0d0*vcosi*(anumz(kx,iy,iz,2)-2.0d0*y*y*anumz(kx,iy,iz,3))
              vboundz(ix,iy,iz, 7,na)=vboundz(ix,iy,iz, 7,na)-4.0d0*vcosi*(anumz(kx,iy,iz,2)-2.0d0*z*z*anumz(kx,iy,iz,3))
              vboundz(ix,iy,iz, 8,na)=vboundz(ix,iy,iz, 8,na)+4.0d0*vkx*vsine*anumz(kx,iy,iz,2)*y
              vboundz(ix,iy,iz, 9,na)=vboundz(ix,iy,iz, 9,na)+8.0d0*vcosi*anumz(kx,iy,iz,3)*y*z
              vboundz(ix,iy,iz,10,na)=vboundz(ix,iy,iz,10,na)+4.0d0*vkx*vsine*anumz(kx,iy,iz,2)*z
            end do
            end do
            end do
          end if
        end do
!$omp do
        do iz=-(nfh-1),0
        do iy=1,ncpy_d
        do ix=1,ncpx_d
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
          y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
          z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
          ta=(y*y+z*z)
          vlo=(-0.5d0*expint1(ta*veta*veta,ierr)-0.5d0*dlog(ta))/xmax
          vlo0r=2.0d0*(1.0d0-dexp(-veta*veta*ta))/(ta*2.0d0*xmax)
          vlorr=-2.0d0*(1.0d0-dexp(-veta*veta*ta))/(ta*ta*2.0d0*xmax) &
                +2.0d0*veta*veta*dexp(-veta*veta*ta)/(ta*2.0d0*xmax)
          vboundz(ix,iy,iz, 1,na)=vboundz(ix,iy,iz, 1,na)+vlo
          vboundz(ix,iy,iz, 3,na)=vboundz(ix,iy,iz, 3,na)+vlo0r*y
          vboundz(ix,iy,iz, 4,na)=vboundz(ix,iy,iz, 4,na)+vlo0r*z
          vboundz(ix,iy,iz, 6,na)=vboundz(ix,iy,iz, 6,na)-(vlo0r+2.0d0*y*y*vlorr)
          vboundz(ix,iy,iz, 7,na)=vboundz(ix,iy,iz, 7,na)-(vlo0r+2.0d0*z*z*vlorr)
          vboundz(ix,iy,iz, 9,na)=vboundz(ix,iy,iz, 9,na)-2.0d0*y*z*vlorr
        end do
        end do
        end do
        do kx=-new_rsx,new_rsx
!$omp do
          do iz=-(nfh-1),0
          do iy=1,ncpy_d
          do ix=1,ncpx_d
            jx=myrx*ncpx_d+ix
            jy=myry*ncpy_d+iy
            jz=myrz*ncpz_d+iz
            x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)+kx*xmax*2.0d0
            y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
            z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
            r=dsqrt(x*x+y*y+z*z)
            rin=1.0d0/r
            r2in=rin*rin
            r3in=r2in*rin
            derftmp=(1.0d0-derf(veta*r))
            dexptmp=dexp(-veta*veta*r*r)/pi05
            vlo0=-2.0d0*veta*dexptmp*rin-derftmp*r2in
            vlo1=4.0d0*veta*dexptmp*(veta*veta+r2in)+2.0d0*derftmp*r3in
            vlo2=-vlo0
            vlo0r=vlo0*rin
            vlo1r=vlo1*r2in
            vlo2r=vlo2*r3in
            vlo2rrr=vlo2*rin
            vlo12r=(vlo1r+vlo2r)
            vboundz(ix,iy,iz, 1,na)=vboundz(ix,iy,iz, 1,na)+rin*(-derf(veta*r)+1.0d0)
            vboundz(ix,iy,iz, 2,na)=vboundz(ix,iy,iz, 2,na)-vlo0r*x
            vboundz(ix,iy,iz, 3,na)=vboundz(ix,iy,iz, 3,na)-vlo0r*y
            vboundz(ix,iy,iz, 4,na)=vboundz(ix,iy,iz, 4,na)-vlo0r*z
            vboundz(ix,iy,iz, 5,na)=vboundz(ix,iy,iz, 5,na)+(x*x*vlo12r-vlo2rrr)
            vboundz(ix,iy,iz, 6,na)=vboundz(ix,iy,iz, 6,na)+(y*y*vlo12r-vlo2rrr)
            vboundz(ix,iy,iz, 7,na)=vboundz(ix,iy,iz, 7,na)+(z*z*vlo12r-vlo2rrr)
            vboundz(ix,iy,iz, 8,na)=vboundz(ix,iy,iz, 8,na)+x*y*vlo12r
            vboundz(ix,iy,iz, 9,na)=vboundz(ix,iy,iz, 9,na)+y*z*vlo12r
            vboundz(ix,iy,iz,10,na)=vboundz(ix,iy,iz,10,na)+z*x*vlo12r
          end do
          end do
          end do
        end do
      end if
      if (myrz==nprocz-1) then
!$omp do
        do iz=1,nfh
        do iy=1,ncpy_d
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
            vkx=pi/xmax*kx
            vk2=vkx**2
            jy=myry*ncpy_d+iy
            jz=myrz*ncpz_d+ncpz_d+iz
            y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
            z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
            ta=(y*y+z*z)
            tb=vkx*vkx
            anumz(kx,iy,iz,1)=0.0d0
            anumz(kx,iy,iz,2)=0.0d0
            anumz(kx,iy,iz,3)=0.0d0
            dt=veta/nint1dmax
            do i=1,nint1dmax
              t=i*dt
              anumz(kx,iy,iz,1)=anumz(kx,iy,iz,1)+dexp(-ta*t*t-0.25d0*tb/(t*t))/t*dt
              anumz(kx,iy,iz,2)=anumz(kx,iy,iz,2)+dexp(-ta*t*t-0.25d0*tb/(t*t))*t*dt
              anumz(kx,iy,iz,3)=anumz(kx,iy,iz,3)+dexp(-ta*t*t-0.25d0*tb/(t*t))*t*t*t*dt
            end do
          end if
        end do
        end do
        end do
        do kx=-new_pwx,new_pwx
          if (kx*kx .ne. 0) then
            vkx=pi/xmax*kx
            vk2=vkx**2
!$omp do
            do iz=1,nfh
            do iy=1,ncpy_d
            do ix=1,ncpx_d
              jx=myrx*ncpx_d+ix
              jy=myry*ncpy_d+iy
              jz=myrz*ncpz_d+ncpz_d+iz
              x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
              y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
              z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
              vcosi=dcos(vkx*x)/(xmax*2.0d0)
              vsine=dsin(vkx*x)/(xmax*2.0d0)
              vboundz(ix,iy,iz, 1,na)=vboundz(ix,iy,iz, 1,na)+2.0d0*vcosi*anumz(kx,iy,iz,1)
              vboundz(ix,iy,iz, 2,na)=vboundz(ix,iy,iz, 2,na)+2.0d0*vkx*vsine*anumz(kx,iy,iz,1)
              vboundz(ix,iy,iz, 3,na)=vboundz(ix,iy,iz, 3,na)+4.0d0*vcosi*anumz(kx,iy,iz,2)*y
              vboundz(ix,iy,iz, 4,na)=vboundz(ix,iy,iz, 4,na)+4.0d0*vcosi*anumz(kx,iy,iz,2)*z
              vboundz(ix,iy,iz, 5,na)=vboundz(ix,iy,iz, 5,na)-2.0d0*vkx*vkx*vcosi*anumz(kx,iy,iz,1)
              vboundz(ix,iy,iz, 6,na)=vboundz(ix,iy,iz, 6,na)-4.0d0*vcosi*(anumz(kx,iy,iz,2)-2.0d0*y*y*anumz(kx,iy,iz,3))
              vboundz(ix,iy,iz, 7,na)=vboundz(ix,iy,iz, 7,na)-4.0d0*vcosi*(anumz(kx,iy,iz,2)-2.0d0*z*z*anumz(kx,iy,iz,3))
              vboundz(ix,iy,iz, 8,na)=vboundz(ix,iy,iz, 8,na)+4.0d0*vkx*vsine*anumz(kx,iy,iz,2)*y
              vboundz(ix,iy,iz, 9,na)=vboundz(ix,iy,iz, 9,na)+8.0d0*vcosi*anumz(kx,iy,iz,3)*y*z
              vboundz(ix,iy,iz,10,na)=vboundz(ix,iy,iz,10,na)+4.0d0*vkx*vsine*anumz(kx,iy,iz,2)*z
            end do
            end do
            end do
          end if
        end do
!$omp do
        do iz=1,nfh
        do iy=1,ncpy_d
        do ix=1,ncpx_d
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+ncpz_d+iz
          x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
          y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
          z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
          ta=(y*y+z*z)
          vlo=(-0.5d0*expint1(ta*veta*veta,ierr)-0.5d0*dlog(ta))/xmax
          vlo0r=2.0d0*(1.0d0-dexp(-veta*veta*ta))/(ta*2.0d0*xmax)
          vlorr=-2.0d0*(1.0d0-dexp(-veta*veta*ta))/(ta*ta*2.0d0*xmax) &
                +2.0d0*veta*veta*dexp(-veta*veta*ta)/(ta*2.0d0*xmax)
          vboundz(ix,iy,iz, 1,na)=vboundz(ix,iy,iz, 1,na)+vlo
          vboundz(ix,iy,iz, 3,na)=vboundz(ix,iy,iz, 3,na)+vlo0r*y
          vboundz(ix,iy,iz, 4,na)=vboundz(ix,iy,iz, 4,na)+vlo0r*z
          vboundz(ix,iy,iz, 6,na)=vboundz(ix,iy,iz, 6,na)-(vlo0r+2.0d0*y*y*vlorr)
          vboundz(ix,iy,iz, 7,na)=vboundz(ix,iy,iz, 7,na)-(vlo0r+2.0d0*z*z*vlorr)
          vboundz(ix,iy,iz, 9,na)=vboundz(ix,iy,iz, 9,na)-2.0d0*y*z*vlorr
        end do
        end do
        end do
        do kx=-new_rsx,new_rsx
!$omp do
          do iz=1,nfh
          do iy=1,ncpy_d
          do ix=1,ncpx_d
            jx=myrx*ncpx_d+ix
            jy=myry*ncpy_d+iy
            jz=myrz*ncpz_d+ncpz_d+iz
            x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)+kx*xmax*2.0d0
            y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
            z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
            r=dsqrt(x*x+y*y+z*z)
            rin=1.0d0/r
            r2in=rin*rin
            r3in=r2in*rin
            derftmp=(1.0d0-derf(veta*r))
            dexptmp=dexp(-veta*veta*r*r)/pi05
            vlo0=-2.0d0*veta*dexptmp*rin-derftmp*r2in
            vlo1=4.0d0*veta*dexptmp*(veta*veta+r2in)+2.0d0*derftmp*r3in
            vlo2=-vlo0
            vlo0r=vlo0*rin
            vlo1r=vlo1*r2in
            vlo2r=vlo2*r3in
            vlo2rrr=vlo2*rin
            vlo12r=(vlo1r+vlo2r)
            vboundz(ix,iy,iz, 1,na)=vboundz(ix,iy,iz, 1,na)+rin*(-derf(veta*r)+1.0d0)
            vboundz(ix,iy,iz, 2,na)=vboundz(ix,iy,iz, 2,na)-vlo0r*x
            vboundz(ix,iy,iz, 3,na)=vboundz(ix,iy,iz, 3,na)-vlo0r*y
            vboundz(ix,iy,iz, 4,na)=vboundz(ix,iy,iz, 4,na)-vlo0r*z
            vboundz(ix,iy,iz, 5,na)=vboundz(ix,iy,iz, 5,na)+(x*x*vlo12r-vlo2rrr)
            vboundz(ix,iy,iz, 6,na)=vboundz(ix,iy,iz, 6,na)+(y*y*vlo12r-vlo2rrr)
            vboundz(ix,iy,iz, 7,na)=vboundz(ix,iy,iz, 7,na)+(z*z*vlo12r-vlo2rrr)
            vboundz(ix,iy,iz, 8,na)=vboundz(ix,iy,iz, 8,na)+x*y*vlo12r
            vboundz(ix,iy,iz, 9,na)=vboundz(ix,iy,iz, 9,na)+y*z*vlo12r
            vboundz(ix,iy,iz,10,na)=vboundz(ix,iy,iz,10,na)+z*x*vlo12r
          end do
          end do
          end do
        end do
      end if
    end do
  case (2)
    do na=1,natom
    if (myrz==0) then
      pisurf=2.0d0/(4.0d0*xmax*ymax)*dsqrt(pi)
      do ky=-new_pwy,new_pwy
      do kx=-new_pwx,new_pwx
      if (kx*kx+ky*ky .ne. 0) then
         vkx=pi/xmax*kx
         vky=pi/ymax*ky
         vk2=vkx**2+vky**2
!$omp do
        do iz=-(nfh-1),0
        do iy=1,ncpy_d
        do ix=1,ncpx_d
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
          y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
          z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
          ta=z
          tb=dsqrt(vk2)
          deftp1=(1.0d0-derf((tb-2.0d0*veta*veta*z)/(2.0d0*veta)))*dexp(-tb*z)
          deftp2=(1.0d0-derf((tb+2.0d0*veta*veta*z)/(2.0d0*veta)))*dexp( tb*z)
          vep0=pi05/tb        *0.500d0*(deftp1+deftp2)
          vep1=pi05/ta        *0.250d0*(deftp1-deftp2)
          vep2=pi05/(ta*ta*ta)*0.125d0*(deftp1*(1.0d0+ta*tb)-deftp2*(1.0d0-ta*tb)) &
               -0.5d0*veta/(ta*ta)*dexp(-(tb*tb+4.0d0*ta*ta*veta*veta*veta*veta)/(4.0d0*veta*veta))
          vsine=dsin(vkx*x+vky*y)*pisurf
          vcosi=dcos(vkx*x+vky*y)*pisurf
          vboundz(ix,iy,iz, 1,na)=vboundz(ix,iy,iz, 1,na)+vcosi*vep0
          vboundz(ix,iy,iz, 2,na)=vboundz(ix,iy,iz, 2,na)+vsine*vep0*vkx
          vboundz(ix,iy,iz, 3,na)=vboundz(ix,iy,iz, 3,na)+vsine*vep0*vky
          vboundz(ix,iy,iz, 4,na)=vboundz(ix,iy,iz, 4,na)+2.0d0*vcosi*vep1*z
          vboundz(ix,iy,iz, 5,na)=vboundz(ix,iy,iz, 5,na)-vcosi*vep0*vkx*vkx
          vboundz(ix,iy,iz, 6,na)=vboundz(ix,iy,iz, 6,na)-vcosi*vep0*vky*vky
          vboundz(ix,iy,iz, 7,na)=vboundz(ix,iy,iz, 7,na)-vcosi*(2.0d0*vep1-4.0d0*z*z*vep2)
          vboundz(ix,iy,iz, 8,na)=vboundz(ix,iy,iz, 8,na)-vcosi*vep0*vkx*vky
          vboundz(ix,iy,iz, 9,na)=vboundz(ix,iy,iz, 9,na)+2.0d0*vsine*vep1*vky*z
          vboundz(ix,iy,iz,10,na)=vboundz(ix,iy,iz,10,na)+2.0d0*vsine*vep1*z*vkx
        end do
        end do
        end do
      end if
      end do
      end do
!$omp do
      do iz=-(nfh-1),0
      do iy=1,ncpy_d
      do ix=1,ncpx_d
        jz=myrz*ncpz_d+iz
        z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
        zab=dabs(z)
        deftmp0=derf(veta*zab)
        dextmp0=dexp(-zab*zab*veta*veta)
        vlo0=2.0d0*pi05*surfin*(dextmp0/veta+pi05*zab*deftmp0)
        vlo1=-2.0d0*pi*surfin*deftmp0*dsign(1.0d0,z)
        vlo2=4.0d0*pi05*veta*surfin*dextmp0
        vboundz(ix,iy,iz, 1,na)=vboundz(ix,iy,iz, 1,na)-vlo0
        vboundz(ix,iy,iz, 4,na)=vboundz(ix,iy,iz, 4,na)-vlo1
        vboundz(ix,iy,iz, 7,na)=vboundz(ix,iy,iz, 7,na)-vlo2
      end do
      end do
      end do
      do ky=-new_rsy,new_rsy
      do kx=-new_rsx,new_rsx
!$omp do
        do iz=-(nfh-1),0
        do iy=1,ncpy_d
        do ix=1,ncpx_d
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+iz
          x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)+kx*xmax*2.0d0
          y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)+ky*ymax*2.0d0
          z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
          r=dsqrt(x*x+y*y+z*z)
          rin=1.0d0/r
          r2in=rin*rin
          r3in=r2in*rin
          derftmp=(1.0d0-derf(veta*r))
          dexptmp=dexp(-veta*veta*r*r)/pi05
          vlo0=-2.0d0*veta*dexptmp*rin-derftmp*r2in
          vlo1=4.0d0*veta*dexptmp*(veta*veta+r2in)+2.0d0*derftmp*r3in
          vlo2=-vlo0
          vlo0r=vlo0*rin
          vlo1r=vlo1*r2in
          vlo2r=vlo2*r3in
          vlo2rrr=vlo2*rin
          vlo12r=(vlo1r+vlo2r)
          vboundz(ix,iy,iz, 1,na)=vboundz(ix,iy,iz, 1,na)+rin*(-derf(veta*r)+1.0d0)
          vboundz(ix,iy,iz, 2,na)=vboundz(ix,iy,iz, 2,na)-vlo0r*x
          vboundz(ix,iy,iz, 3,na)=vboundz(ix,iy,iz, 3,na)-vlo0r*y
          vboundz(ix,iy,iz, 4,na)=vboundz(ix,iy,iz, 4,na)-vlo0r*z
          vboundz(ix,iy,iz, 5,na)=vboundz(ix,iy,iz, 5,na)+(x*x*vlo12r-vlo2rrr)
          vboundz(ix,iy,iz, 6,na)=vboundz(ix,iy,iz, 6,na)+(y*y*vlo12r-vlo2rrr)
          vboundz(ix,iy,iz, 7,na)=vboundz(ix,iy,iz, 7,na)+(z*z*vlo12r-vlo2rrr)
          vboundz(ix,iy,iz, 8,na)=vboundz(ix,iy,iz, 8,na)+x*y*vlo12r
          vboundz(ix,iy,iz, 9,na)=vboundz(ix,iy,iz, 9,na)+y*z*vlo12r
          vboundz(ix,iy,iz,10,na)=vboundz(ix,iy,iz,10,na)+z*x*vlo12r
        end do
        end do
        end do
      end do
      end do
    end if
    if (myrz==nprocz-1) then
      pisurf=2.0d0/(4.0d0*xmax*ymax)*dsqrt(pi)
      do ky=-new_pwy,new_pwy
      do kx=-new_pwx,new_pwx
      if (kx*kx+ky*ky .ne. 0) then
         vkx=pi/xmax*kx
         vky=pi/ymax*ky
         vk2=vkx**2+vky**2
!$omp do
        do iz=1,nfh
        do iy=1,ncpy_d
        do ix=1,ncpx_d
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+ncpz_d+iz
          x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)
          y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)
          z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
          ta=z
          tb=dsqrt(vk2)
          deftp1=(1.0d0-derf((tb-2.0d0*veta*veta*z)/(2.0d0*veta)))*dexp(-tb*z)
          deftp2=(1.0d0-derf((tb+2.0d0*veta*veta*z)/(2.0d0*veta)))*dexp( tb*z)
          vep0=pi05/tb        *0.500d0*(deftp1+deftp2)
          vep1=pi05/ta        *0.250d0*(deftp1-deftp2)
          vep2=pi05/(ta*ta*ta)*0.125d0*(deftp1*(1.0d0+ta*tb)-deftp2*(1.0d0-ta*tb)) &
               -0.5d0*veta/(ta*ta)*dexp(-(tb*tb+4.0d0*ta*ta*veta*veta*veta*veta)/(4.0d0*veta*veta))
          vsine=dsin(vkx*x+vky*y)*pisurf
          vcosi=dcos(vkx*x+vky*y)*pisurf
          vboundz(ix,iy,iz, 1,na)=vboundz(ix,iy,iz, 1,na)+vcosi*vep0
          vboundz(ix,iy,iz, 2,na)=vboundz(ix,iy,iz, 2,na)+vsine*vep0*vkx
          vboundz(ix,iy,iz, 3,na)=vboundz(ix,iy,iz, 3,na)+vsine*vep0*vky
          vboundz(ix,iy,iz, 4,na)=vboundz(ix,iy,iz, 4,na)+2.0d0*vcosi*vep1*z
          vboundz(ix,iy,iz, 5,na)=vboundz(ix,iy,iz, 5,na)-vcosi*vep0*vkx*vkx
          vboundz(ix,iy,iz, 6,na)=vboundz(ix,iy,iz, 6,na)-vcosi*vep0*vky*vky
          vboundz(ix,iy,iz, 7,na)=vboundz(ix,iy,iz, 7,na)-vcosi*(2.0d0*vep1-4.0d0*z*z*vep2)
          vboundz(ix,iy,iz, 8,na)=vboundz(ix,iy,iz, 8,na)-vcosi*vep0*vkx*vky
          vboundz(ix,iy,iz, 9,na)=vboundz(ix,iy,iz, 9,na)+2.0d0*vsine*vep1*vky*z
          vboundz(ix,iy,iz,10,na)=vboundz(ix,iy,iz,10,na)+2.0d0*vsine*vep1*z*vkx
        end do
        end do
        end do
      end if
      end do
      end do
!$omp do
      do iz=1,nfh
      do iy=1,ncpy_d
      do ix=1,ncpx_d
        jz=myrz*ncpz_d+ncpz_d+iz
        z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
        zab=dabs(z)
        deftmp0=derf(veta*zab)
        dextmp0=dexp(-zab*zab*veta*veta)
        vlo0=2.0d0*pi05*surfin*(dextmp0/veta+pi05*zab*deftmp0)
        vlo1=-2.0d0*pi*surfin*deftmp0*dsign(1.0d0,z)
        vlo2=4.0d0*pi05*veta*surfin*dextmp0
        vboundz(ix,iy,iz, 1,na)=vboundz(ix,iy,iz, 1,na)-vlo0
        vboundz(ix,iy,iz, 4,na)=vboundz(ix,iy,iz, 4,na)-vlo1
        vboundz(ix,iy,iz, 7,na)=vboundz(ix,iy,iz, 7,na)-vlo2
      end do
      end do
      end do
      do ky=-new_rsy,new_rsy
      do kx=-new_rsx,new_rsx
!$omp do
        do iz=1,nfh
        do iy=1,ncpy_d
        do ix=1,ncpx_d
          jx=myrx*ncpx_d+ix
          jy=myry*ncpy_d+iy
          jz=myrz*ncpz_d+ncpz_d+iz
          x=(jx*ddx-xmax-0.5d0*ddx)-atx(na)+kx*xmax*2.0d0
          y=(jy*ddy-ymax-0.5d0*ddy)-aty(na)+ky*ymax*2.0d0
          z=(jz*ddz-zmax-0.5d0*ddz)-atz(na)
          r=dsqrt(x*x+y*y+z*z)
          rin=1.0d0/r
          r2in=rin*rin
          r3in=r2in*rin
          derftmp=(1.0d0-derf(veta*r))
          dexptmp=dexp(-veta*veta*r*r)/pi05
          vlo0=-2.0d0*veta*dexptmp*rin-derftmp*r2in
          vlo1=4.0d0*veta*dexptmp*(veta*veta+r2in)+2.0d0*derftmp*r3in
          vlo2=-vlo0
          vlo0r=vlo0*rin
          vlo1r=vlo1*r2in
          vlo2r=vlo2*r3in
          vlo2rrr=vlo2*rin
          vlo12r=(vlo1r+vlo2r)
          vboundz(ix,iy,iz, 1,na)=vboundz(ix,iy,iz, 1,na)+rin*(-derf(veta*r)+1.0d0)
          vboundz(ix,iy,iz, 2,na)=vboundz(ix,iy,iz, 2,na)-vlo0r*x
          vboundz(ix,iy,iz, 3,na)=vboundz(ix,iy,iz, 3,na)-vlo0r*y
          vboundz(ix,iy,iz, 4,na)=vboundz(ix,iy,iz, 4,na)-vlo0r*z
          vboundz(ix,iy,iz, 5,na)=vboundz(ix,iy,iz, 5,na)+(x*x*vlo12r-vlo2rrr)
          vboundz(ix,iy,iz, 6,na)=vboundz(ix,iy,iz, 6,na)+(y*y*vlo12r-vlo2rrr)
          vboundz(ix,iy,iz, 7,na)=vboundz(ix,iy,iz, 7,na)+(z*z*vlo12r-vlo2rrr)
          vboundz(ix,iy,iz, 8,na)=vboundz(ix,iy,iz, 8,na)+x*y*vlo12r
          vboundz(ix,iy,iz, 9,na)=vboundz(ix,iy,iz, 9,na)+y*z*vlo12r
          vboundz(ix,iy,iz,10,na)=vboundz(ix,iy,iz,10,na)+z*x*vlo12r
        end do
        end do
        end do
      end do
      end do
    end if
    end do
  end select
  if (nperi==1) deallocate(anumy,anumz)

return
end subroutine pseudocalc_13


subroutine pseudocalc_broadcast( &
 nprjmx,natom,num_list,num_ppcell,                                           & ! <
 natpri,natpri_inf,naps,natinf,natx,naty,natz,lstx,lsty,lstz,lstvec2,vnlocp)   ! X
use mod_mpi
implicit none
integer, intent(in)    :: nprjmx,natom,num_list,num_ppcell
integer, intent(inout) :: natpri(natom)
integer, intent(inout) :: natpri_inf(natom)
integer, intent(inout) :: naps(natom)
integer, intent(inout) :: natinf(natom)
integer, intent(inout) :: natx(natom),naty(natom),natz(natom)
integer, intent(inout) :: lstx(num_list,num_ppcell)
integer, intent(inout) :: lsty(num_list,num_ppcell)
integer, intent(inout) :: lstz(num_list,num_ppcell)
integer, intent(inout) :: lstvec2(num_list,num_ppcell)
real*8,  intent(inout) :: vnlocp(num_list,nprjmx,num_ppcell)

  call mpi_bcast(natpri,natom,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(natpri_inf,natom,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(naps,natom,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(natinf,natom,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(natx,natom,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(naty,natom,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(natz,natom,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(lstx,num_list*num_ppcell,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(lsty,num_list*num_ppcell,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(lstz,num_list*num_ppcell,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(lstvec2,num_list*num_ppcell,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(vnlocp,num_list*nprjmx*num_ppcell,mpi_double_precision,0,mpicom_kpt,mpij)

end subroutine pseudocalc_broadcast


end module
