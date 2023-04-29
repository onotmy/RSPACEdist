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
! **********  reorderstates8f.F90 06/20/2018-01  **********

module mod_reorderstates
implicit none

contains

subroutine reorderstates(nrc,num_ppcell,nprjmx,nei000gmx,nums,ncol,numk,ncpx,ncpy,ncpz,ns1,nk,    & ! <
                        sval,residual_states,rspsep,cspsep,svecre,sveccm,ssvre,ssvcm,hsvre,hsvcm)   ! X
use mod_mpi, only: myrank_glbl
implicit none
integer,    intent(in)   ::nrc,num_ppcell,nprjmx,nei000gmx,nums,ncol,numk,ncpx,ncpy,ncpz,ns1,nk
real*8,     intent(inout)::sval(nei000gmx,nums+1-ncol,numk),residual_states(nei000gmx,nums+1-ncol,numk)
real*8,     intent(inout)::rspsep(nprjmx*(1-nrc)+nrc,num_ppcell*(1-nrc)+nrc,nei000gmx*(1-nrc)+nrc, &
                           nums*(1-nrc)+nrc,numk*(1-nrc)+nrc)
complex*16, intent(inout)::cspsep(nprjmx*nrc-nrc+1,num_ppcell*nrc-nrc+1,nei000gmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1)
real*8,     intent(inout)::svecre(ncpx*(1-nrc)+nrc,ncpy*(1-nrc)+nrc,ncpz*(1-nrc)+nrc &
                                  ,nei000gmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk*(1-nrc)+nrc)
complex*16, intent(inout)::sveccm(ncpx*nrc-nrc+1,ncpy*nrc-nrc+1,ncpz*nrc-nrc+1,nei000gmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1)
real*8,     intent(inout)::ssvre(ncpx*(1-nrc)+nrc,ncpy*(1-nrc)+nrc,ncpz*(1-nrc)+nrc &
                                  ,nei000gmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk*(1-nrc)+nrc)
complex*16, intent(inout)::ssvcm(ncpx*nrc-nrc+1,ncpy*nrc-nrc+1,ncpz*nrc-nrc+1,nei000gmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1)
real*8,     intent(inout)::hsvre(ncpx*(1-nrc)+nrc,ncpy*(1-nrc)+nrc,ncpz*(1-nrc)+nrc &
                                  ,nei000gmx*(1-nrc)+nrc)
complex*16, intent(inout)::hsvcm(ncpx*nrc-nrc+1,ncpy*nrc-nrc+1,ncpz*nrc-nrc+1,nei000gmx*nrc-nrc+1,ncol*nrc-nrc+1)
integer i,j, ns,ns2
real*8,allocatable::rpsep(:,:,:,:),sval_tmp(:),residual_tmp(:)
complex*16,allocatable::cpsep(:,:,:,:)
integer,allocatable::lank(:)

  allocate(lank(nei000gmx),sval_tmp(nei000gmx),residual_tmp(nei000gmx))
  if (nrc==0) then
    allocate(rpsep(nprjmx,num_ppcell,nei000gmx,ncol))
  else
    allocate(cpsep(nprjmx,num_ppcell,nei000gmx,ncol))
  endif

  lank(1)=1
  do j=2,nei000gmx
    lank(j)=1
    do i=1,j-1
      if (sval(i,ns1,nk) .gt. sval(j,ns1,nk)) then
        lank(i)=lank(i)+1
      else
        lank(j)=lank(j)+1
      end if
    end do
  end do

  !$omp parallel default(shared) private(j,ns,ns2)

  if (nrc==0) then
    !$omp do
    do j=1,nei000gmx
      sval_tmp(lank(j))=sval(j,ns1,nk)
      residual_tmp(lank(j))=residual_states(j,ns1,nk)
      rpsep(:,:,lank(j),1)=rspsep(:,:,j,ns1,nk)
      hsvre(:,:,:,lank(j))=svecre(:,:,:,j,ns1,nk)
    end do
    !$omp do
    do j=1,nei000gmx
      sval(j,ns1,nk)=sval_tmp(j)
      residual_states(j,ns1,nk)=residual_tmp(j)
      rspsep(:,:,j,ns1,nk)=rpsep(:,:,j,1)
      svecre(:,:,:,j,ns1,nk)=hsvre(:,:,:,j)
    enddo  
    !$omp do
    do j=1,nei000gmx
      hsvre(:,:,:,lank(j))=ssvre(:,:,:,j,ns1,nk)
    end do
    !$omp do
    do j=1,nei000gmx
      ssvre(:,:,:,j,ns1,nk)=hsvre(:,:,:,j)
    end do
  end if

  if (nrc==1) then
    !$omp do
    do j=1,nei000gmx
      sval_tmp(lank(j))=sval(j,ns1,nk)
      residual_tmp(lank(j))=residual_states(j,ns1,nk)
    end do
    !$omp do
    do j=1,nei000gmx
      sval(j,ns1,nk)=sval_tmp(j)
      residual_states(j,ns1,nk)=residual_tmp(j)
    end do
    do ns= 1,ncol
      ns2= ns1+ns-1
      !$omp do
      do j=1,nei000gmx
        cpsep(:,:,lank(j),ns)=cspsep(:,:,j,ns2,nk)
        hsvcm(:,:,:,lank(j),ns)=sveccm(:,:,:,j,ns2,nk)
      end do
      !$omp do
      do j=1,nei000gmx
        cspsep(:,:,j,ns2,nk)=cpsep(:,:,j,ns)
        sveccm(:,:,:,j,ns2,nk)=hsvcm(:,:,:,j,ns)
      enddo 
      !$omp do
      do j=1,nei000gmx
        hsvcm(:,:,:,lank(j),ns)=ssvcm(:,:,:,j,ns2,nk)
      end do
      !$omp do
      do j=1,nei000gmx
        ssvcm(:,:,:,j,ns2,nk)=hsvcm(:,:,:,j,ns)
      end do 
    end do
  end if

  !$omp end parallel

  if (nrc==0) then
    deallocate(rpsep)
  else
    deallocate(cpsep)
  endif
  deallocate(lank,sval_tmp,residual_tmp)

  return
end subroutine


end module
