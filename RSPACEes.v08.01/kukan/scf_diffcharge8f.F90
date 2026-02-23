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

! In the magnetic case, the distance is defined as { (n_a-n_b)^2 + (m_a-m_b)^2 }/2 .
! Thereby, m can have (x,y,z)-components. Note, that
! ( Up_a - Up_b )^2 + ( Down_a - Down_b )^2
!  = { [ ( Up_a + Down_a ) - ( Up_b + Down_b ) ]^2 + [ ( Up_a - Down_a ) - ( Up_b - Down_b ) ]^2 }/2
!  = { [       n_a         -       n_b         ]^2 + [       mz_a        -       mz_b        ]^2 }/2 .

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
