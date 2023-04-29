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
! **********  potout8f.F90 12/05/2022-01 **********

subroutine potout(chdir,ndisp,nspv,nums,ncpx,ncpy,ncpz,myrank_glbl,nprocx,nprocy,nprocz,mpicom_space,mpi_status_size, &
                  mpi_double_precision,mpistat,ferm,veff)
implicit none
character, intent(in) :: chdir*200
integer,   intent(in) :: ndisp,nspv,nums,ncpx,ncpy,ncpz,myrank_glbl,nprocx,nprocy,nprocz
integer,   intent(in) :: mpicom_space,mpi_status_size,mpi_double_precision
integer               :: mpistat(mpi_status_size)
real*8,    intent(in) :: veff(ncpx,ncpy,ncpz,nspv)
real*8,   allocatable :: vrecv(:,:,:),veff_all(:,:,:)
integer ix,iy,iz,jx,jy,jz,j,ns,myrx,myry,myrz,mx,my,mz,ierr
real*8 ferm,tmp
character :: fname*200

  allocate (vrecv(ncpx,ncpy,ncpz))
  if (myrank_glbl == 0) allocate(veff_all(ncpx*nprocx,ncpy*nprocy,ncpz*nprocz))

  myrz=myrank_glbl/(nprocx*nprocy)
  myry=(myrank_glbl-myrz*nprocx*nprocy)/nprocx
  myrx=myrank_glbl-myrz*nprocx*nprocy-myry*nprocx

  if (myrank_glbl == 0) then
    fname='Potential.txt'
    if (len_trim(chdir) > 0) fname=trim(chdir)//'/'//fname
    open(90,file=fname,form='formatted')
  end if

  do ns=1,nums
    if (myrank_glbl /= 0) call mpi_send(veff(1,1,1,ns),ncpx*ncpy*ncpz,mpi_double_precision,0,1,mpicom_space,ierr)

    if (myrank_glbl == 0) then
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        veff_all(ix,iy,iz)=veff(ix,iy,iz,ns)
      end do
      end do
      end do

      do j=1,nprocx*nprocy*nprocz-1
        mz=j/(nprocx*nprocy)
        my=(j-mz*nprocx*nprocy)/nprocx
        mx=j-mz*nprocx*nprocy-my*nprocx
        call mpi_recv(vrecv,ncpx*ncpy*ncpz,mpi_double_precision,j,1,mpicom_space,mpistat,ierr)
        do iz=1,ncpz
        do iy=1,ncpy
        do ix=1,ncpx
          jx=mx*ncpx+ix
          jy=my*ncpy+iy
          jz=mz*ncpz+iz
          veff_all(jx,jy,jz)=vrecv(ix,iy,iz)
        end do
        end do
        end do
      end do
      do ix=1,ncpx*nprocx*ncpy*nprocy*ncpz*nprocz
        write(90,'(d25.16)') veff_all(ix,1,1)
      end do

      tmp=0.0d0
      do iy=1,ncpy*nprocy
      do ix=1,ncpx*nprocx
        tmp=tmp+veff_all(ix,iy,1)/(ncpy*nprocy*ncpx*nprocx)
      end do
      end do
      write(ndisp,*) 'E_F from the bottom of V_eff at the boundary for spin (left)',ns
      write(ndisp,*) ferm-tmp,'(hartree)'
      tmp=0.0d0
      do iy=1,ncpy*nprocy
      do ix=1,ncpx*nprocx
        tmp=tmp+veff_all(ix,iy,ncpz*nprocz)/(ncpy*nprocy*ncpx*nprocx)
      end do
      end do
      write(ndisp,*) 'E_F from the bottom of V_eff at the boundary for spin (right)',ns
      write(ndisp,*) ferm-tmp,'(hartree)'
      end if ! (myrank_glbl == 0)
  end do !ns

  if (myrank_glbl == 0) close(90)
  if (myrank_glbl == 0) deallocate(veff_all)

  deallocate (vrecv)

  return
end subroutine
