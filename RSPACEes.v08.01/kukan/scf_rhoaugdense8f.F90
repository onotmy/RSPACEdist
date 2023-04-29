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
! **********  scf_rhoaugdense8e.f90 11/19/2013-01  **********

module mod_scf_rhoaugdense
implicit none

contains

subroutine scf_rhoaugdense( &
 key_natpri_inps,key_pp_paw,                                      & ! <
 natom,num_spe,nspv,ncpx_d,ncpy_d,ncpz_d,num_list_d,num_ppcell_d, & ! <
 indspe,natprid,napsd,natinfd,ntyppp,lstvecd2,                    & ! <
 rho_dense,rhoaug3d,                                              & ! <
 rho_aug_dense)                                                     ! >
implicit none
integer,intent(in) ::key_natpri_inps,key_pp_paw
integer,intent(in) ::natom,num_spe,nspv
integer,intent(in) ::ncpx_d,ncpy_d,ncpz_d,num_list_d,num_ppcell_d
integer,intent(in) ::indspe(natom),natprid(natom),napsd(natom),natinfd(natom)
integer,intent(in) ::ntyppp(num_spe)
integer,intent(in) ::lstvecd2(num_list_d,num_ppcell_d)
real*8, intent(in) ::rho_dense(ncpx_d,ncpy_d,ncpz_d,nspv)
real*8, intent(in) ::rhoaug3d(num_list_d,num_ppcell_d)
real*8, intent(out)::rho_aug_dense(ncpx_d,ncpy_d,ncpz_d)
integer na,iapsd,i

!$omp do
  do i=1,ncpx_d*ncpy_d*ncpz_d
    rho_aug_dense(i,1,1)= rho_dense(i,1,1,1)
  end do
  if (nspv>1) then
!$omp do
    do i=1,ncpx_d*ncpy_d*ncpz_d
      rho_aug_dense(i,1,1)= rho_aug_dense(i,1,1) +rho_dense(i,1,1,2)
    end do
  end if

  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natprid(na) .eq. key_natpri_inps)) then
      iapsd=napsd(na)
!$omp do
      do i=1,natinfd(na)
        rho_aug_dense(lstvecd2(i,iapsd),1,1)=rho_aug_dense(lstvecd2(i,iapsd),1,1)+rhoaug3d(i,iapsd)
      end do
!$omp barrier
    end if
  end do

  end subroutine


end module
