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
! **********  scf_rhoxcalc8e.f90 10/04/2016-01  **********

module mod_scf_rhoxcalc
implicit none
contains


subroutine scf_rhoxcalc_smtpcc(  &
 nspv,ncpx_d,ncpy_d,ncpz_d, & ! <
 rhopcc_dense,rho_dense,    & ! <
 rhosmt_pcc_dense)            ! >
integer,intent(in) ::nspv,ncpx_d,ncpy_d,ncpz_d
real*8, intent(in) ::rhopcc_dense(ncpx_d,ncpy_d,ncpz_d),rho_dense(ncpx_d,ncpy_d,ncpz_d,nspv)
real*8, intent(out)::rhosmt_pcc_dense(ncpx_d,ncpy_d,ncpz_d,nspv)
integer ns,ix
real*8  rspv

!$omp parallel default(shared) private(ns,ix,rspv)
  rspv= dble(min(2,nspv))
  do ns=1,nspv
    if (ns<3) then
!$omp do
      do ix=1,ncpx_d*ncpy_d*ncpz_d
        rhosmt_pcc_dense(ix,1,1,ns)=rho_dense(ix,1,1,ns)+rhopcc_dense(ix,1,1)/rspv
      end do
    else
!$omp do
      do ix=1,ncpx_d*ncpy_d*ncpz_d
        rhosmt_pcc_dense(ix,1,1,ns)=rho_dense(ix,1,1,ns)
      end do
    end if
  end do
!$omp end parallel
  return
end subroutine scf_rhoxcalc_smtpcc


subroutine scf_rhoxcalc_onecenter( &
 key_natpri_in,key_pp_paw,                    & ! <
 nspv,natom,num_spe,num_atcell,npoint,nradmx, & ! <
 indspe,natpri,natpri_inf,ntyppp,nradct,      & ! <
 rhocore,rhopccr,rhotrur,rhosmtr,             & ! <
 rhotrucorer,rhosmt_pccr)                       ! >
use mod_mpi
integer, intent(in)  :: key_natpri_in,key_pp_paw
integer, intent(in)  :: nspv,natom,num_spe,num_atcell,npoint,nradmx
integer, intent(in)  :: indspe(natom),natpri(natom),natpri_inf(natom),ntyppp(num_spe),nradct(num_spe)
real*8,  intent(in)  :: rhocore(nradmx,num_spe),rhopccr(nradmx,npoint,num_atcell)
real*8,  intent(in)  :: rhotrur(nradmx,npoint,nspv,num_atcell),rhosmtr(nradmx,npoint,nspv,num_atcell)
real*8,  intent(out) :: rhotrucorer(nradmx,npoint,nspv,num_atcell),rhosmt_pccr(nradmx,npoint,nspv,num_atcell)
integer:: na,ipri,ns,il,ir
real*8 :: rspv

!$omp parallel default(shared) private(na,ipri,ns,il,ir,rspv)

rspv= dble(min(2,nspv))

do na=1,natom
  if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
    ipri=natpri_inf(na)
    do ns=1,nspv
      if (ns<3) then
!$omp do
      do il=1,npoint
        do ir=2,nradct(indspe(na))-1
          rhotrucorer(ir,il,ns,ipri)=rhotrur(ir,il,ns,ipri)+(rhocore(ir,indspe(na))+rhopccr(ir,il,ipri))/rspv
          rhosmt_pccr(ir,il,ns,ipri)=rhosmtr(ir,il,ns,ipri)+rhopccr(ir,il,ipri)/rspv
        end do
        do ir=nradct(indspe(na)),nradmx
          rhotrucorer(ir,il,ns,ipri)=0.0d0
          rhosmt_pccr(ir,il,ns,ipri)=0.0d0
        end do
      end do
      else
!$omp do
      do il=1,npoint
        do ir=2,nradct(indspe(na))-1
          rhotrucorer(ir,il,ns,ipri)=rhotrur(ir,il,ns,ipri)
          rhosmt_pccr(ir,il,ns,ipri)=rhosmtr(ir,il,ns,ipri)
        end do
        do ir=nradct(indspe(na)),nradmx
          rhotrucorer(ir,il,ns,ipri)=0.0d0
          rhosmt_pccr(ir,il,ns,ipri)=0.0d0
        end do
      end do
      end if
    end do
  end if
end do

!$omp end parallel

end subroutine scf_rhoxcalc_onecenter


end module
