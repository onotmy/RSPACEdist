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
!     **********  scf_diffcharge8f.F90 09/21/2012-01  **********

module mod_scf_diffcharge
implicit none

contains

subroutine scf_diffcharge(ncpx,ncpy,ncpz,nspv,rhosmt_o,rhosmt_i,diff_charge_all)
use mod_mpi
implicit none
integer,intent(in) ::ncpx,ncpy,ncpz,nspv
real*8, intent(in) ::rhosmt_o(ncpx,ncpy,ncpz,nspv)
real*8, intent(in) ::rhosmt_i(ncpx,ncpy,ncpz,nspv)
real*8, intent(out)::diff_charge_all
integer ix,ns
real*8 diff_charge

  diff_charge_all= 0.0d0
!$omp parallel default(shared) private(ns,ix,diff_charge)
  diff_charge=0.0d0
  do ns=nspv,1,-1
!$omp do
    do ix=1,ncpx*ncpy*ncpz
      diff_charge=diff_charge+(rhosmt_i(ix,1,1,ns)-rhosmt_o(ix,1,1,ns))**2
    end do
    if (ns==3) diff_charge= diff_charge*2.0d0
  end do
!$omp critical
  diff_charge_all= diff_charge_all +diff_charge
!$omp end critical
!$omp end parallel

  diff_charge= diff_charge_all
  call mpi_reduce(diff_charge,diff_charge_all,1,mpi_double_precision,mpi_sum,0,mpicom_space,mpij)

  end subroutine

end module
