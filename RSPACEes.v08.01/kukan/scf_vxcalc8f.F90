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
! **********  scf_vxcalc8f.F90 01/12/2016-01  **********

module mod_scf_vxcalc
implicit none
contains

subroutine scf_vxcalc(  &
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
use mod_vxpot, only:vxpot_lda_vwn,vxpot_lda_vwns,vxpot_lda_pz,vxpot_lda_pzs &
                   ,vxpot_gga_pw91,vxpot_gga_pw91s,vxpot_gga_pbe,vxpot_gga_pbes
use mod_vxpot_ggaxyz
use mod_vxpot_ggartp
use mod_stopp
integer,  intent(in) :: natom,nperi,nfh,num_spe,nspv,npoint,nradmx,lrhomx,num_atcell
integer,  intent(in) :: ncpx_d,ncpy_d,ncpz_d
integer,  intent(in) :: key_natpri_in
real*8,   intent(in) :: ddx,ddy,ddz
character,intent(in) :: cexco*7
integer,  intent(in) :: indspe(natom),natpri(natom),natpri_inf(natom),nradct(num_spe)
real*8,   intent(in) :: yylm(npoint,lrhomx),dylm_dtheta(npoint,lrhomx)
real*8,   intent(in) :: d2ylm_dtheta2(npoint,lrhomx),dylm_dphi(npoint,lrhomx)
real*8,   intent(in) :: d2ylm_dphi2(npoint,lrhomx),d2ylm_dtheta_dphi(npoint,lrhomx)
real*8,   intent(in) :: wt(npoint),point(npoint,3)
real*8,   intent(in) :: radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8,   intent(in) :: rhosmt_pcc_dense(ncpx_d,ncpy_d,ncpz_d,nspv)
real*8,   intent(in) :: rhotrucorer(nradmx,npoint,nspv,num_atcell)
real*8,   intent(in) :: rhosmt_pccr(nradmx,npoint,nspv,num_atcell)
real*8,   intent(out):: vx_dense(ncpx_d,ncpy_d,ncpz_d,nspv)
real*8,   intent(out):: vxctru(nradmx,npoint,nspv,num_atcell)
real*8,   intent(out):: vxcsmt(nradmx,npoint,nspv,num_atcell)
real*8,   intent(out):: ex_dense(ncpx_d,ncpy_d,ncpz_d)
real*8,   intent(out):: exctru(nradmx,npoint,num_atcell)
real*8,   intent(out):: excsmt(nradmx,npoint,num_atcell)
integer :: nspin,nrad,na,ipri,ispe
real*8, allocatable :: abgr(:,:,:,:),ggabg(:,:,:,:),ggr(:,:,:,:)
real*8, allocatable ::agr(:,:,:),gggr(:,:,:),g2r(:,:,:)
real*8, allocatable::rhodia(:,:,:),vxcdia(:,:,:),mag(:,:,:)
real*8, allocatable::rhototal(:,:,:)


  ! *********************************************************************
  select case (trim(cexco))

  case ('vwn') ! l(s)da

    nspin= min(nspv,2)

    if (nspv<4) then
      select case (nspin)
      case (1) ! no spin
        !$omp parallel default(shared)
        call vxpot_lda_vwn(ncpx_d*ncpy_d*ncpz_d,rhosmt_pcc_dense, vx_dense, ex_dense)
        !$omp end parallel
      case (2) ! spin free
        !$omp parallel default(shared)
        call vxpot_lda_vwns(ncpx_d*ncpy_d*ncpz_d,rhosmt_pcc_dense, vx_dense, ex_dense)
        !$omp end parallel
      end select
    else
      allocate(rhodia(ncpx_d*ncpy_d,ncpz_d,2),vxcdia(ncpx_d*ncpy_d,ncpz_d,2),mag(ncpx_d*ncpy_d,ncpz_d,3))
      !$omp parallel default(shared)
      call scf_vxcalc_rhomagloc(ncpx_d*ncpy_d,ncpx_d*ncpy_d,ncpz_d,rhosmt_pcc_dense, rhodia,mag)
      call vxpot_lda_vwns(ncpx_d*ncpy_d*ncpz_d,rhodia,vxcdia,ex_dense)
      call scf_vxcalc_rhomagglo(ncpx_d*ncpy_d,ncpx_d*ncpy_d,ncpz_d,vxcdia,mag, vx_dense)
      !$omp end parallel
      deallocate(rhodia,vxcdia,mag)
    end if

    do na= 1,natom
      if (natpri(na)==key_natpri_in) then
        ipri=natpri_inf(na)
        nrad= nradct(indspe(na))

        if (nspv<4) then
          select case (nspin)
          case (1) ! no spin
            !$omp parallel default(shared)
            call vxpot_lda_vwn(nradmx*npoint,rhosmt_pccr(1,1,1,ipri),vxcsmt(1,1,1,ipri),excsmt(1,1,ipri))
            call vxpot_lda_vwn(nradmx*npoint,rhotrucorer(1,1,1,ipri),vxctru(1,1,1,ipri),exctru(1,1,ipri))
            !$omp end parallel
          case (2) ! spin free
            !$omp parallel default(shared)
            call vxpot_lda_vwns(nradmx*npoint,rhosmt_pccr(1,1,1,ipri),vxcsmt(1,1,1,ipri),excsmt(1,1,ipri))
            call vxpot_lda_vwns(nradmx*npoint,rhotrucorer(1,1,1,ipri),vxctru(1,1,1,ipri),exctru(1,1,ipri))
            !$omp end parallel
          end select
        else
          allocate(rhodia(nradmx,npoint,2),vxcdia(nradmx,npoint,2),mag(nrad,npoint,3))
          !$omp parallel default(shared)
          call scf_vxcalc_rhomagloc(nradmx,nrad,npoint,rhosmt_pccr(1,1,1,ipri), rhodia,mag)
          call vxpot_lda_vwns(nradmx*npoint,rhodia,vxcdia,excsmt(1,1,ipri))
          call scf_vxcalc_rhomagglo(nradmx,nrad,npoint,vxcdia,mag, vxcsmt(1,1,1,ipri))
          call scf_vxcalc_rhomagloc(nradmx,nrad,npoint,rhotrucorer(1,1,1,ipri), rhodia,mag)
          call vxpot_lda_vwns(nradmx*npoint,rhodia,vxcdia,exctru(1,1,ipri))
          call scf_vxcalc_rhomagglo(nradmx,nrad,npoint,vxcdia,mag, vxctru(1,1,1,ipri))
          !$omp end parallel
          deallocate(rhodia,vxcdia,mag)
        end if
      end if
    end do

  case ('pz') ! l(s)da

    nspin= min(nspv,2)

    if (nspv<4) then
      select case (nspin)
      case (1) ! no spin
        !$omp parallel default(shared)
        call vxpot_lda_pz(ncpx_d*ncpy_d*ncpz_d,rhosmt_pcc_dense, vx_dense, ex_dense)
        !$omp end parallel
      case (2) ! spin free
        !$omp parallel default(shared)
        call vxpot_lda_pzs(ncpx_d*ncpy_d*ncpz_d,rhosmt_pcc_dense, vx_dense, ex_dense)
        !$omp end parallel
      end select
    else
      allocate(rhodia(ncpx_d*ncpy_d,ncpz_d,2),vxcdia(ncpx_d*ncpy_d,ncpz_d,2),mag(ncpx_d*ncpy_d,ncpz_d,3))
      !$omp parallel default(shared)
      call scf_vxcalc_rhomagloc(ncpx_d*ncpy_d,ncpx_d*ncpy_d,ncpz_d,rhosmt_pcc_dense, rhodia,mag)
      call vxpot_lda_pzs(ncpx_d*ncpy_d*ncpz_d,rhodia,vxcdia,ex_dense)
      call scf_vxcalc_rhomagglo(ncpx_d*ncpy_d,ncpx_d*ncpy_d,ncpz_d,vxcdia,mag, vx_dense)
      !$omp end parallel
      deallocate(rhodia,vxcdia,mag)
    end if

    do na= 1,natom
      if (natpri(na)==key_natpri_in) then
        ipri=natpri_inf(na)
        nrad= nradct(indspe(na))

        if (nspv<4) then
          select case (nspin)
          case (1) ! no spin
            !$omp parallel default(shared)
            call vxpot_lda_pz(nradmx*npoint,rhosmt_pccr(1,1,1,ipri),vxcsmt(1,1,1,ipri),excsmt(1,1,ipri))
            call vxpot_lda_pz(nradmx*npoint,rhotrucorer(1,1,1,ipri),vxctru(1,1,1,ipri),exctru(1,1,ipri))
            !$omp end parallel
          case (2) ! spin free
            !$omp parallel default(shared)
            call vxpot_lda_pzs(nradmx*npoint,rhosmt_pccr(1,1,1,ipri),vxcsmt(1,1,1,ipri),excsmt(1,1,ipri))
            call vxpot_lda_pzs(nradmx*npoint,rhotrucorer(1,1,1,ipri),vxctru(1,1,1,ipri),exctru(1,1,ipri))
            !$omp end parallel
          end select
        else
          allocate(rhodia(nradmx,npoint,2),vxcdia(nradmx,npoint,2),mag(nrad,npoint,3))
          !$omp parallel default(shared)
          call scf_vxcalc_rhomagloc(nradmx,nrad,npoint,rhosmt_pccr(1,1,1,ipri), rhodia,mag)
          call vxpot_lda_pzs(nradmx*npoint,rhodia,vxcdia,excsmt(1,1,ipri))
          call scf_vxcalc_rhomagglo(nradmx,nrad,npoint,vxcdia,mag, vxcsmt(1,1,1,ipri))
          call scf_vxcalc_rhomagloc(nradmx,nrad,npoint,rhotrucorer(1,1,1,ipri), rhodia,mag)
          call vxpot_lda_pzs(nradmx*npoint,rhodia,vxcdia,exctru(1,1,ipri))
          call scf_vxcalc_rhomagglo(nradmx,nrad,npoint,vxcdia,mag, vxctru(1,1,1,ipri))
          !$omp end parallel
          deallocate(rhodia,vxcdia,mag)
        end if
      end if
    end do

  case ('pw91') ! GGA(PW91)

    nspin= min(nspv,2)

    if (nspv<4) then
      select case (nspin)
      case (1) ! no spin
        allocate(abgr(ncpx_d,ncpy_d,ncpz_d,nspin),ggabg(ncpx_d,ncpy_d,ncpz_d,nspin),ggr(ncpx_d,ncpy_d,ncpz_d,nspin))
        call vxpot_ggaxyz(nperi,nfh,ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,rhosmt_pcc_dense,abgr,ggabg,ggr)
        !$omp parallel default(shared)
        call vxpot_gga_pw91(ncpx_d*ncpy_d*ncpz_d,rhosmt_pcc_dense,abgr,ggabg,ggr,vx_dense,ex_dense)
        !$omp end parallel
      case (2) ! spin free
        allocate(abgr(ncpx_d,ncpy_d,ncpz_d,nspin+1),ggabg(ncpx_d,ncpy_d,ncpz_d,nspin+1),ggr(ncpx_d,ncpy_d,ncpz_d,nspin+1))
        allocate(rhototal(ncpx_d,ncpy_d,ncpz_d))
        rhototal(:,:,:)=rhosmt_pcc_dense(:,:,:,1)+rhosmt_pcc_dense(:,:,:,2)
        call vxpot_ggaxyz(nperi,nfh,ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,rhosmt_pcc_dense(1,1,1,1) &
                         ,abgr(1,1,1,1),ggabg(1,1,1,1),ggr(1,1,1,1))
        call vxpot_ggaxyz(nperi,nfh,ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,rhosmt_pcc_dense(1,1,1,2) &
                         ,abgr(1,1,1,2),ggabg(1,1,1,2),ggr(1,1,1,2))
        call vxpot_ggaxyz(nperi,nfh,ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,rhototal &
                         ,abgr(1,1,1,3),ggabg(1,1,1,3),ggr(1,1,1,3))
        !$omp parallel default(shared)
        call vxpot_gga_pw91s(ncpx_d*ncpy_d*ncpz_d,rhosmt_pcc_dense,abgr,ggabg,ggr,vx_dense,ex_dense)
        !$omp end parallel
      end select
    else
      allocate(abgr(ncpx_d,ncpy_d,ncpz_d,nspin+1),ggabg(ncpx_d,ncpy_d,ncpz_d,nspin+1),ggr(ncpx_d,ncpy_d,ncpz_d,nspin+1))
      allocate(rhototal(ncpx_d*ncpy_d,ncpz_d,1))
      allocate(rhodia(ncpx_d*ncpy_d,ncpz_d,2),vxcdia(ncpx_d*ncpy_d,ncpz_d,2),mag(ncpx_d*ncpy_d,ncpz_d,3))
      rhototal(:,:,1)=rhodia(:,:,1)+rhodia(:,:,2)
      !$omp parallel default(shared)
      call scf_vxcalc_rhomagloc(ncpx_d*ncpy_d,ncpx_d*ncpy_d,ncpz_d,rhosmt_pcc_dense,rhodia,mag)
      !$omp end parallel
      call vxpot_ggaxyz(nperi,nfh,ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,rhodia(1,1,1),abgr(1,1,1,1),ggabg(1,1,1,1),ggr(1,1,1,1))
      call vxpot_ggaxyz(nperi,nfh,ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,rhodia(1,1,2),abgr(1,1,1,2),ggabg(1,1,1,2),ggr(1,1,1,2))
      call vxpot_ggaxyz(nperi,nfh,ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,rhototal(1,1,1),abgr(1,1,1,3),ggabg(1,1,1,3),ggr(1,1,1,3))
      !$omp parallel default(shared)
      call vxpot_gga_pw91s(ncpx_d*ncpy_d*ncpz_d,rhodia,abgr,ggabg,ggr,vxcdia,ex_dense)
      call scf_vxcalc_rhomagglo(ncpx_d*ncpy_d,ncpx_d*ncpy_d,ncpz_d,vxcdia,mag,vx_dense)
      !$omp end parallel
      deallocate(rhodia,vxcdia,mag)
    end if
    deallocate(abgr,ggabg,ggr)
    if (nspin==2) deallocate (rhototal)

    select case (nspin)
    case(1)
      allocate(agr(nradmx,npoint,nspin),gggr(nradmx,npoint,nspin),g2r(nradmx,npoint,nspin))
    case(2)
      allocate(agr(nradmx,npoint,nspin+1),gggr(nradmx,npoint,nspin+1),g2r(nradmx,npoint,nspin+1))
      allocate(rhototal(nradmx,npoint,1))
    case(4)
      allocate(agr(nradmx,npoint,nspin+1),gggr(nradmx,npoint,nspin+1),g2r(nradmx,npoint,nspin+1))
      allocate(rhototal(nradmx,npoint,1))
      allocate(rhodia(nradmx,npoint,2),vxcdia(nradmx,npoint,2),mag(nrad,npoint,3))
    end select
    do na= 1,natom
      if (natpri(na)==key_natpri_in) then
        ipri=natpri_inf(na)
        ispe=indspe(na)
        nrad= nradct(indspe(na))

        if (nspv<4) then
          select case (nspin)
          case (1) ! no spin
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhosmt_pccr(1,1,1,ipri)                           & ! <
                              ,agr,gggr,g2r)                                                                                 ! >
            !$omp parallel default(shared)
            call vxpot_gga_pw91(nradmx*npoint,rhosmt_pccr(1,1,1,ipri),agr,gggr,g2r,vxcsmt(1,1,1,ipri),excsmt(1,1,ipri))
            !$omp end parallel
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhotrucorer(1,1,1,ipri)                           & ! <
                              ,agr,gggr,g2r)                                                                                 ! >
            !$omp parallel default(shared)
            call vxpot_gga_pw91(nradmx*npoint,rhotrucorer(1,1,1,ipri),agr,gggr,g2r,vxctru(1,1,1,ipri),exctru(1,1,ipri))
            !$omp end parallel
          case (2) ! spin free
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhosmt_pccr(1,1,1,ipri)                           & ! <
                              ,agr(1,1,1),gggr(1,1,1),g2r(1,1,1))                                                            ! >
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhosmt_pccr(1,1,2,ipri)                           & ! <
                              ,agr(1,1,2),gggr(1,1,2),g2r(1,1,2))                                                            ! >
            rhototal(:,:,1)=rhosmt_pccr(:,:,1,ipri)+rhosmt_pccr(:,:,2,ipri)
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhototal(1,1,1)                                   & ! <
                              ,agr(1,1,3),gggr(1,1,3),g2r(1,1,3))                                                            ! >
            !$omp parallel default(shared)
            call vxpot_gga_pw91s(nradmx*npoint,rhosmt_pccr(1,1,1,ipri),agr,gggr,g2r,vxcsmt(1,1,1,ipri),excsmt(1,1,ipri))
            !$omp end parallel
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhotrucorer(1,1,1,ipri)                           & ! <
                              ,agr(1,1,1),gggr(1,1,1),g2r(1,1,1))                                                            ! >
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhotrucorer(1,1,2,ipri)                           & ! <
                              ,agr(1,1,2),gggr(1,1,2),g2r(1,1,2))                                                            ! >
            rhototal(:,:,1)=rhotrucorer(:,:,1,ipri)+rhotrucorer(:,:,2,ipri)
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhototal(1,1,1)                                   & ! <
                              ,agr(1,1,3),gggr(1,1,3),g2r(1,1,3))                                                            ! >
            !$omp parallel default(shared)
            call vxpot_gga_pw91s(nradmx*npoint,rhotrucorer(1,1,1,ipri),agr,gggr,g2r,vxctru(1,1,1,ipri),exctru(1,1,ipri))
            !$omp end parallel
          end select
        else
          !$omp parallel default(shared)
          call scf_vxcalc_rhomagloc(nradmx,nrad,npoint,rhosmt_pccr(1,1,1,ipri),rhodia,mag)
          !$omp end parallel
          call vxpot_ggartp(num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                           ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhodia(1,1,1)                                     & ! <
                           ,agr(1,1,1),gggr(1,1,1),g2r(1,1,1))                                                            ! >
          call vxpot_ggartp(num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                           ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhodia(1,1,2)                                     & ! <
                           ,agr(1,1,2),gggr(1,1,2),g2r(1,1,2))                                                            ! >
          rhototal(:,:,1)=rhodia(:,:,1)+rhodia(:,:,2)
          call vxpot_ggartp(num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                           ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhototal(1,1,1)                                   & ! <
                           ,agr(1,1,3),gggr(1,1,3),g2r(1,1,3))                                                            ! >
          !$omp parallel default(shared)
          call vxpot_gga_pw91s(nradmx*npoint,rhodia,agr,gggr,g2r,vxcdia,excsmt(1,1,ipri))
          call scf_vxcalc_rhomagglo(nradmx,nrad,npoint,vxcdia,mag, vxcsmt(1,1,1,ipri))
          call scf_vxcalc_rhomagloc(nradmx,nrad,npoint,rhotrucorer(1,1,1,ipri), rhodia,mag)
          !$omp end parallel
          call vxpot_ggartp(num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                           ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhodia(1,1,1)                                     & ! <
                           ,agr(1,1,1),gggr(1,1,1),g2r(1,1,1))                                                            ! >
          call vxpot_ggartp(num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                           ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhodia(1,1,2)                                     & ! <
                           ,agr(1,1,2),gggr(1,1,2),g2r(1,1,2))                                                            ! >
          rhototal(:,:,1)=rhodia(:,:,1)+rhodia(:,:,2)
          call vxpot_ggartp(num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                           ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhototal(1,1,1)                                   & ! <
                           ,agr(1,1,3),gggr(1,1,3),g2r(1,1,3))                                                            ! >
          !$omp parallel default(shared)
          call vxpot_gga_pw91s(nradmx*npoint,rhodia,agr,gggr,g2r,vxcdia,exctru(1,1,ipri))
          call scf_vxcalc_rhomagglo(nradmx,nrad,npoint,vxcdia,mag, vxctru(1,1,1,ipri))
          !$omp end parallel
        end if
      end if
    end do
    deallocate(agr,gggr,g2r)
    if (nspin==2) deallocate (rhototal)
    if (nspin==4) deallocate(rhodia,vxcdia,mag)

  case ('pbe') ! GGA(PBE)

    nspin= min(nspv,2)

    if (nspv<4) then
      select case (nspin)
      case (1) ! no spin
        allocate(abgr(ncpx_d,ncpy_d,ncpz_d,nspin),ggabg(ncpx_d,ncpy_d,ncpz_d,nspin),ggr(ncpx_d,ncpy_d,ncpz_d,nspin))
        call vxpot_ggaxyz(nperi,nfh,ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,rhosmt_pcc_dense,abgr,ggabg,ggr)
        !$omp parallel default(shared)
        call vxpot_gga_pbe(ncpx_d*ncpy_d*ncpz_d,rhosmt_pcc_dense,abgr,ggabg,ggr,vx_dense,ex_dense)
        !$omp end parallel
      case (2) ! spin free
        allocate(abgr(ncpx_d,ncpy_d,ncpz_d,nspin+1),ggabg(ncpx_d,ncpy_d,ncpz_d,nspin+1),ggr(ncpx_d,ncpy_d,ncpz_d,nspin+1))
        allocate(rhototal(ncpx_d,ncpy_d,ncpz_d))
        rhototal(:,:,:)=rhosmt_pcc_dense(:,:,:,1)+rhosmt_pcc_dense(:,:,:,2)
        call vxpot_ggaxyz(nperi,nfh,ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,rhosmt_pcc_dense(1,1,1,1) &
                         ,abgr(1,1,1,1),ggabg(1,1,1,1),ggr(1,1,1,1))
        call vxpot_ggaxyz(nperi,nfh,ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,rhosmt_pcc_dense(1,1,1,2) &
                         ,abgr(1,1,1,2),ggabg(1,1,1,2),ggr(1,1,1,2))
        call vxpot_ggaxyz(nperi,nfh,ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,rhototal &
                         ,abgr(1,1,1,3),ggabg(1,1,1,3),ggr(1,1,1,3))
        !$omp parallel default(shared)
        call vxpot_gga_pbes(ncpx_d*ncpy_d*ncpz_d,rhosmt_pcc_dense,abgr,ggabg,ggr,vx_dense,ex_dense)
        !$omp end parallel
      end select
    else
      allocate(abgr(ncpx_d,ncpy_d,ncpz_d,nspin+1),ggabg(ncpx_d,ncpy_d,ncpz_d,nspin+1),ggr(ncpx_d,ncpy_d,ncpz_d,nspin+1))
      allocate(rhototal(ncpx_d*ncpy_d,ncpz_d,1))
      allocate(rhodia(ncpx_d*ncpy_d,ncpz_d,2),vxcdia(ncpx_d*ncpy_d,ncpz_d,2),mag(ncpx_d*ncpy_d,ncpz_d,3))
      rhototal(:,:,1)=rhodia(:,:,1)+rhodia(:,:,2)
      !$omp parallel default(shared)
      call scf_vxcalc_rhomagloc(ncpx_d*ncpy_d,ncpx_d*ncpy_d,ncpz_d,rhosmt_pcc_dense,rhodia,mag)
      !$omp end parallel
      call vxpot_ggaxyz(nperi,nfh,ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,rhodia(1,1,1),abgr(1,1,1,1),ggabg(1,1,1,1),ggr(1,1,1,1))
      call vxpot_ggaxyz(nperi,nfh,ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,rhodia(1,1,2),abgr(1,1,1,2),ggabg(1,1,1,2),ggr(1,1,1,2))
      call vxpot_ggaxyz(nperi,nfh,ncpx_d,ncpy_d,ncpz_d,ddx,ddy,ddz,rhototal(1,1,1),abgr(1,1,1,3),ggabg(1,1,1,3),ggr(1,1,1,3))
      !$omp parallel default(shared)
      call vxpot_gga_pbes(ncpx_d*ncpy_d*ncpz_d,rhodia,abgr,ggabg,ggr,vxcdia,ex_dense)
      call scf_vxcalc_rhomagglo(ncpx_d*ncpy_d,ncpx_d*ncpy_d,ncpz_d,vxcdia,mag,vx_dense)
      !$omp end parallel
      deallocate(rhodia,vxcdia,mag)
    end if
    deallocate(abgr,ggabg,ggr)
    if (nspin==2) deallocate (rhototal)

    select case (nspin)
    case(1)
      allocate(agr(nradmx,npoint,nspin),gggr(nradmx,npoint,nspin),g2r(nradmx,npoint,nspin))
    case(2)
      allocate(agr(nradmx,npoint,nspin+1),gggr(nradmx,npoint,nspin+1),g2r(nradmx,npoint,nspin+1))
      allocate(rhototal(nradmx,npoint,1))
    case(4)
      allocate(agr(nradmx,npoint,nspin+1),gggr(nradmx,npoint,nspin+1),g2r(nradmx,npoint,nspin+1))
      allocate(rhototal(nradmx,npoint,1))
      allocate(rhodia(nradmx,npoint,2),vxcdia(nradmx,npoint,2),mag(nrad,npoint,3))
    end select
    do na= 1,natom
      if (natpri(na)==key_natpri_in) then
        ipri=natpri_inf(na)
        ispe=indspe(na)
        nrad= nradct(indspe(na))

        if (nspv<4) then
          select case (nspin)
          case (1) ! no spin
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhosmt_pccr(1,1,1,ipri)                           & ! <
                              ,agr,gggr,g2r)                                                                                 ! >
            !$omp parallel default(shared)
            call vxpot_gga_pbe(nradmx*npoint,rhosmt_pccr(1,1,1,ipri),agr,gggr,g2r,vxcsmt(1,1,1,ipri),excsmt(1,1,ipri))
            !$omp end parallel
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhotrucorer(1,1,1,ipri)                           & ! <
                              ,agr,gggr,g2r)                                                                                 ! >
            !$omp parallel default(shared)
            call vxpot_gga_pbe(nradmx*npoint,rhotrucorer(1,1,1,ipri),agr,gggr,g2r,vxctru(1,1,1,ipri),exctru(1,1,ipri))
            !$omp end parallel
          case (2) ! spin free
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhosmt_pccr(1,1,1,ipri)                           & ! <
                              ,agr(1,1,1),gggr(1,1,1),g2r(1,1,1))                                                            ! >
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhosmt_pccr(1,1,2,ipri)                           & ! <
                              ,agr(1,1,2),gggr(1,1,2),g2r(1,1,2))                                                            ! >
            rhototal(:,:,1)=rhosmt_pccr(:,:,1,ipri)+rhosmt_pccr(:,:,2,ipri)
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhototal(1,1,1)                                   & ! <
                              ,agr(1,1,3),gggr(1,1,3),g2r(1,1,3))                                                            ! >
            !$omp parallel default(shared)
            call vxpot_gga_pbes(nradmx*npoint,rhosmt_pccr(1,1,1,ipri),agr,gggr,g2r,vxcsmt(1,1,1,ipri),excsmt(1,1,ipri))
            !$omp end parallel
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhotrucorer(1,1,1,ipri)                           & ! <
                              ,agr(1,1,1),gggr(1,1,1),g2r(1,1,1))                                                            ! >
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhotrucorer(1,1,2,ipri)                           & ! <
                              ,agr(1,1,2),gggr(1,1,2),g2r(1,1,2))                                                            ! >
            rhototal(:,:,1)=rhotrucorer(:,:,1,ipri)+rhotrucorer(:,:,2,ipri)
            call vxpot_ggartp (num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                              ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhototal(1,1,1)                                   & ! <
                              ,agr(1,1,3),gggr(1,1,3),g2r(1,1,3))                                                            ! >
            !$omp parallel default(shared)
            call vxpot_gga_pbes(nradmx*npoint,rhotrucorer(1,1,1,ipri),agr,gggr,g2r,vxctru(1,1,1,ipri),exctru(1,1,ipri))
            !$omp end parallel
          end select
        else
          !$omp parallel default(shared)
          call scf_vxcalc_rhomagloc(nradmx,nrad,npoint,rhosmt_pccr(1,1,1,ipri),rhodia,mag)
          !$omp end parallel
          call vxpot_ggartp(num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                           ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhodia(1,1,1)                                     & ! <
                           ,agr(1,1,1),gggr(1,1,1),g2r(1,1,1))                                                            ! >
          call vxpot_ggartp(num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                           ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhodia(1,1,2)                                     & ! <
                           ,agr(1,1,2),gggr(1,1,2),g2r(1,1,2))                                                            ! >
          rhototal(:,:,1)=rhodia(:,:,1)+rhodia(:,:,2)
          call vxpot_ggartp(num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                           ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhototal(1,1,1)                                   & ! <
                           ,agr(1,1,3),gggr(1,1,3),g2r(1,1,3))                                                            ! >
          !$omp parallel default(shared)
          call vxpot_gga_pbes(nradmx*npoint,rhodia,agr,gggr,g2r,vxcdia,excsmt(1,1,ipri))
          call scf_vxcalc_rhomagglo(nradmx,nrad,npoint,vxcdia,mag, vxcsmt(1,1,1,ipri))
          call scf_vxcalc_rhomagloc(nradmx,nrad,npoint,rhotrucorer(1,1,1,ipri), rhodia,mag)
          !$omp end parallel
          call vxpot_ggartp(num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                           ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhodia(1,1,1)                                     & ! <
                           ,agr(1,1,1),gggr(1,1,1),g2r(1,1,1))                                                            ! >
          call vxpot_ggartp(num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                           ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhodia(1,1,2)                                     & ! <
                           ,agr(1,1,2),gggr(1,1,2),g2r(1,1,2))                                                            ! >
          rhototal(:,:,1)=rhodia(:,:,1)+rhodia(:,:,2)
          call vxpot_ggartp(num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                           ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhototal(1,1,1)                                   & ! <
                           ,agr(1,1,3),gggr(1,1,3),g2r(1,1,3))                                                            ! >
          !$omp parallel default(shared)
          call vxpot_gga_pbes(nradmx*npoint,rhodia,agr,gggr,g2r,vxcdia,exctru(1,1,ipri))
          call scf_vxcalc_rhomagglo(nradmx,nrad,npoint,vxcdia,mag, vxctru(1,1,1,ipri))
          !$omp end parallel
        end if
      end if
    end do
    deallocate(agr,gggr,g2r)
    if (nspin==2) deallocate (rhototal)
    if (nspin==4) deallocate(rhodia,vxcdia,mag)

  case default
    call stopp ('scf_vxcalc:  unknown xc potential')
  end select
  ! *********************************************************************

end subroutine scf_vxcalc


subroutine scf_vxcalc_rhomagloc (n1mx,n1,n2,rho,rhodia,mag)
implicit none
integer,intent(in) :: n1mx,n1,n2
real*8, intent(in) :: rho(n1mx,n2,4)
real*8, intent(out):: rhodia(n1mx,n2,2),mag(n1,n2,3)
real*8 :: eps_m ; parameter(eps_m=1.0d-10)
integer i1,i2

!$omp do
  do i2=1,n2
  do i1=1,n1
    rhodia(i1,i2,1)= rho(i1,i2,1)+rho(i1,i2,2)
    mag(i1,i2,3)= rho(i1,i2,1)-rho(i1,i2,2)
    mag(i1,i2,1)= 2.0d0*rho(i1,i2,3)
    mag(i1,i2,2)= 2.0d0*rho(i1,i2,4)
    rhodia(i1,i2,2)= mag(i1,i2,1)**2+mag(i1,i2,2)**2+mag(i1,i2,3)**2
    rhodia(i1,i2,2)= dsqrt(rhodia(i1,i2,2))
    if (rhodia(i1,i2,2)<eps_m) then
      rhodia(i1,i2,2)= 0.0d0
      mag(i1,i2,1)= 0.0d0
      mag(i1,i2,2)= 0.0d0
      mag(i1,i2,3)= 0.0d0
    else
      mag(i1,i2,1)= mag(i1,i2,1)/rhodia(i1,i2,2)
      mag(i1,i2,2)= mag(i1,i2,2)/rhodia(i1,i2,2)
      mag(i1,i2,3)= mag(i1,i2,3)/rhodia(i1,i2,2)
    end if
    rhodia(i1,i2,1)= (rhodia(i1,i2,1)+rhodia(i1,i2,2))/2.0d0
    rhodia(i1,i2,2)= rhodia(i1,i2,1)-rhodia(i1,i2,2)
  end do
  end do
end subroutine scf_vxcalc_rhomagloc


subroutine scf_vxcalc_rhomagglo (n1mx,n1,n2,vxcdia,mag,potxc)
implicit none
integer,intent(in) :: n1mx,n1,n2
real*8, intent(in) :: vxcdia(n1mx,n2,2),mag(n1,n2,3)
real*8, intent(out):: potxc(n1mx,n2,4)
real*8 :: pot1, pot2
integer i1,i2

!$omp do
  do i2=1,n2
  do i1=1,n1
    pot1= (vxcdia(i1,i2,1)+vxcdia(i1,i2,2))/2.0d0
    pot2= (vxcdia(i1,i2,1)-vxcdia(i1,i2,2))/2.0d0
    potxc(i1,i2,1)= pot1+pot2*mag(i1,i2,3)
    potxc(i1,i2,2)= pot1-pot2*mag(i1,i2,3)
    potxc(i1,i2,3)= pot2*mag(i1,i2,1)
    potxc(i1,i2,4)= pot2*mag(i1,i2,2)
  end do
  end do

end subroutine scf_vxcalc_rhomagglo


end module
