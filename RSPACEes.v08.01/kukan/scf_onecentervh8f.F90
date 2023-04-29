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
! **********  scf_onecentervh8e.f90 11/19/2013-01  **********

module mod_scf_onecentervh
implicit none
contains

subroutine scf_onecentervh( &
 key_natpri_in,key_pp_paw,                      & ! <
 natom,num_spe,num_atcell,lrhomx,npoint,nradmx, & ! <
 indspe,natpri,natpri_inf,ntyppp,lpmx,nradct,   & ! <
 yylm,radial,dradial,                           & ! <
 rhotrum,rhosmtm,rhoaugm,                       & ! <
 vhtrur,vhsmtr,vhaugr)                            ! >
use mod_stopp
use mod_radialhartree
implicit none
integer,intent(in) ::key_natpri_in,key_pp_paw
integer,intent(in) ::natom,num_spe,num_atcell,lrhomx,npoint,nradmx
integer,intent(in) ::indspe(natom),natpri(natom),natpri_inf(natom),ntyppp(num_spe),lpmx(num_spe),nradct(num_spe)
real*8, intent(in) ::yylm(npoint,lrhomx),radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8, intent(in) ::rhotrum(nradmx,lrhomx,num_atcell)
real*8, intent(in) ::rhosmtm(nradmx,lrhomx,num_atcell)
real*8, intent(in) ::rhoaugm(nradmx,lrhomx,num_atcell)
real*8, intent(out)::vhtrur(nradmx,npoint,num_atcell)
real*8, intent(out)::vhsmtr(nradmx,npoint,num_atcell)
real*8, intent(out)::vhaugr(nradmx,npoint,num_atcell)
real*8,allocatable::vtmpm(:,:,:)
integer na,ipri,ispe
  if (lrhomx > 25) call stopp ('scf_onecentervh: lrhomx must be <= 25!')
! We can omit the computation for the contribution of the boundary value of the Coulomb potential,
! because it is cancelled by substracting \int (v_h) Q_{ij}^L dr from \tilda{D}_{ij}
! Definitions of Q_{ij}^L and \tilda{D}_{ij} are presented in PRB59 1758 (1999).
  allocate(vtmpm(nradmx,lrhomx,num_atcell))
! zero clear is carried out in radialhartree_01
!  vhsmtr=0.0d0
!  vhtrur=0.0d0
!  vhaugr=0.0d0
  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
      ipri=natpri_inf(na)
      ispe=indspe(na)
!$omp parallel default(shared)
    call radialhartree_01(nradmx,lrhomx,nradct(ispe),lpmx(ispe),rhosmtm(1,1,ipri),vtmpm(1,1,ipri),radial(1,ispe),dradial(1,ispe))
    call radialhartree_02(nradmx,lrhomx,nradct(ispe),npoint,lpmx(ispe),vhsmtr(1,1,ipri),vtmpm(1,1,ipri),yylm)
    call radialhartree_01(nradmx,lrhomx,nradct(ispe),lpmx(ispe),rhotrum(1,1,ipri),vtmpm(1,1,ipri),radial(1,ispe),dradial(1,ispe))
    call radialhartree_02(nradmx,lrhomx,nradct(ispe),npoint,lpmx(ispe),vhtrur(1,1,ipri),vtmpm(1,1,ipri),yylm)
    call radialhartree_01(nradmx,lrhomx,nradct(ispe),lpmx(ispe),rhoaugm(1,1,ipri),vtmpm(1,1,ipri),radial(1,ispe),dradial(1,ispe))
    call radialhartree_02(nradmx,lrhomx,nradct(ispe),npoint,lpmx(ispe),vhaugr(1,1,ipri),vtmpm(1,1,ipri),yylm)
!$omp end parallel
    end if
  end do
  deallocate(vtmpm)

end subroutine

end module
