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
! **********  module8f.F90 08/21/2018-01  **********

module mod_mpi
implicit none
include 'mpif.h'
integer mpistat(mpi_status_size)
integer mpicom_space,mpicom_kpt,mpicom_rhpt,mpicom_grn,mpicom_ene
integer,allocatable::mpicom_atom(:)
integer mpij
integer ndim_ke
integer,allocatable::numkproc(:),ntilecomm(:)
integer myrank_glbl,myr_space,myrx,myry,myrz,myr_kpt,myr_rhpt,myr_grn,myr_ene
integer nprocworld,nprocw,nprocs,nprocx,nprocy,nprocz,nprock,nprocsk,nprocg,nprocrhpt,nprocene
integer nprocx_inv,nprocy_inv,nprocz_inv,nprocg_inv
real*8 comcount,comtime

integer ncontext_ss,ncontext_frac

end module mod_mpi


subroutine comput_spacerank
use mod_mpi
implicit none
  myrz=myr_space/(nprocx*nprocy)
  myry=(myr_space-myrz*nprocx*nprocy)/nprocx
  myrx=myr_space-myrz*nprocx*nprocy-myry*nprocx
return
end subroutine comput_spacerank


module mod_stopp
implicit none
integer:: ndisp_stopp

contains

subroutine stopfile (looptime)
! useful for mpi debugging
implicit none
integer, intent(in) :: looptime
logical :: stopnow
integer :: jfile
  stopnow= .false.
  do while (.not. stopnow)
    inquire(file='stopfile',exist=stopnow)
    if (stopnow) then
      jfile= 10
      do while (stopnow)
        jfile= jfile+1
        inquire(unit=jfile,opened=stopnow)
      enddo
      open(jfile,file='stopfile',form='formatted')
      read(jfile,fmt=*) stopnow
      close(jfile)
    endif
    if (.not. stopnow) call sleep(looptime)
  enddo
  open(jfile,file='stopfile',form='formatted',status='replace')
  write(jfile,fmt='(l1)') .false.
  close(jfile)
  call stopp('stopfile')
end subroutine stopfile

subroutine stopp_initialize (ndisp)
! this subroutine is necessaray as stopp should be very low in the makefile hirarchy
implicit none
integer,intent(in) :: ndisp
  ndisp_stopp= ndisp
end subroutine stopp_initialize

subroutine stopp0 (cherror)
use mod_mpi
  character(len=*), intent(in) :: cherror
  if (myrank_glbl==0) then
    call stopp (cherror)
  endif
  call mpi_barrier(mpi_comm_world,mpij)
end subroutine stopp0

subroutine stopp (cherror)
use mod_mpi
implicit none
  character(len=*), intent(in) :: cherror
  if (myrank_glbl==0) then
    write(ndisp_stopp,fmt='(1x,a)') cherror
    if (ndisp_stopp/=6) close(ndisp_stopp)
  else
    write(6,*) cherror
  endif
  call mpi_abort(mpi_comm_world,mpij)
  stop
end subroutine stopp

end module mod_stopp


module mod_disptime
use mod_mpi
contains

subroutine starttime(ts)
implicit none
real*8, intent(out)::ts
  ts=mpi_wtime()
end subroutine

subroutine endtime(ndisp,ts,chara)
implicit none
integer, intent(in) ::ndisp
real*8, intent(in)::ts
  character(len=*), intent(in) :: chara
  real*8 te
  te=mpi_wtime()
  if (myrank_glbl .eq. 0) write(ndisp,*) chara,te-ts,'(sec)'
end subroutine

end module mod_disptime


module mod_proccount
contains

subroutine proccount(ndisp)
use mod_mpi
implicit none
integer, intent(in) :: ndisp
integer :: nmpi,n

  if (myrank_glbl==0) then
    write(ndisp,*)
    call mpi_comm_size(mpi_comm_world,nmpi,mpij)
    n= 0 
    !$omp parallel default(shared) 
    !$omp critical 
    n= n+1 
    !$omp end critical
    !$omp end parallel 
    write(ndisp,fmt='(" # MPI processes =",i6,2x,",",2x,"# OMP threads =",i4)') nmpi,n
    write(ndisp,*)
  endif

end subroutine proccount

end module mod_proccount

