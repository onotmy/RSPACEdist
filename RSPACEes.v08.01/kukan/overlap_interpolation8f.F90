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
! **********  overlap_interpolation8f.F90 04/18/2023-01  **********

module mod_overlap_interpolation
implicit none

contains

subroutine overlap_interpolation(ndisp,nperi,ncpa,ncpb,ncpc,nf,nf1,v_cc,nflag)
use mod_mpi,myra=>myrx,myrb=>myry,myrc=>myrz,nproca=>nprocx,nprocb=>nprocy,nprocc=>nprocz
implicit none
integer ndisp,nperi,ncpa,ncpb,ncpc,nf,nf1
integer ia,ib,ic
integer iov,nflag
integer idests0,idests1,idestr0,idestr1
real*8 stime,etime,weightcom
integer req1,req2
real*8 v_cc(-nf1:ncpa+nf,-nf1:ncpb+nf,-nf1:ncpc+nf)
real*8,allocatable::&
       vsenda2(:,:,:),vsenda3(:,:,:) &
      ,vrecva2(:,:,:),vrecva3(:,:,:) &
      ,vsendb2(:,:,:),vsendb3(:,:,:) &
      ,vrecvb2(:,:,:),vrecvb3(:,:,:) &
      ,vsendc2(:,:,:),vsendc3(:,:,:) &
      ,vrecvc2(:,:,:),vrecvc3(:,:,:)

  allocate(vsenda2(nf,-nf1:ncpb+nf,-nf1:ncpc+nf),vsenda3(nf,-nf1:ncpb+nf,-nf1:ncpc+nf) &
          ,vrecva2(nf,-nf1:ncpb+nf,-nf1:ncpc+nf),vrecva3(nf,-nf1:ncpb+nf,-nf1:ncpc+nf) &
          ,vsendb2(-nf1:ncpa+nf,nf,-nf1:ncpc+nf),vsendb3(-nf1:ncpa+nf,nf,-nf1:ncpc+nf) &
          ,vrecvb2(-nf1:ncpa+nf,nf,-nf1:ncpc+nf),vrecvb3(-nf1:ncpa+nf,nf,-nf1:ncpc+nf) &
          ,vsendc2(-nf1:ncpa+nf,-nf1:ncpb+nf,nf),vsendc3(-nf1:ncpa+nf,-nf1:ncpb+nf,nf) &
          ,vrecvc2(-nf1:ncpa+nf,-nf1:ncpb+nf,nf),vrecvc3(-nf1:ncpa+nf,-nf1:ncpb+nf,nf))

  if ((nflag .ne. 1) .and. (nflag .ne. 2)) then
    if (myrank_glbl .eq. 0) write(ndisp,*) 'error in overlap_interpolation. nflag should be 1 or 2!'
      call mpi_barrier(mpi_comm_world,mpij)
      call mpi_finalize(mpij)
    stop
  end if

  if (nflag .eq. 1) then
  if (nperi .eq. 0) then
! **********  communicate in x direction  **********
  if (nproca .ne. 1) then
  do ic=-nf1,ncpc+nf
  do ib=-nf1,ncpb+nf
  do iov=1,nf
     vsenda2(iov,ib,ic)=v_cc(ncpa-nf+iov,ib,ic)
     vsenda3(iov,ib,ic)=v_cc(        iov,ib,ic)
  end do
  end do
  end do
  end if

  stime=mpi_wtime()
  if (nproca .ne. 1) then
  weightcom=2.0d0*nf
  if (nproca .ge. 4) weightcom=4.0d0*nf
  comcount=comcount+2.0d0*nf*(ncpb+2*nf)*(ncpc+2*nf)*weightcom
  vrecva2=0.0d0
  vrecva3=0.0d0
  idests0=myr_space+1
  idests1=myr_space-1
  idestr0=myr_space-1
  idestr1=myr_space+1
  if (myra .eq. nproca-1) idests0=myr_space-(nproca-1)
  if (myra .eq.        0) idests1=myr_space+(nproca-1)
  if (myra .eq.        0) idestr0=myr_space+(nproca-1)
  if (myra .eq. nproca-1) idestr1=myr_space-(nproca-1)

  if (mod(myra,2) .eq. 0) then
    call mpi_isend(vsenda2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests0,0,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr1,1,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if (mod(myra,2) .eq. 1) then
    call mpi_isend(vsenda3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests1,1,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr0,0,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if ((mod(myra,2) .eq. 0) .and. (myra .ne. 0)) then
    call mpi_isend(vsenda3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests1,2,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr0,3,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if ((mod(myra,2) .eq. 1) .and. (myra .ne. nproca-1)) then
    call mpi_isend(vsenda2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests0,3,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr1,2,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  end if
  etime=mpi_wtime()
  comtime=comtime+(etime-stime)

  if (nproca .ne. 1) then
  if ((myra .ne. 0) .and. (myra .ne. nproca-1)) then
  do ic=-nf1,ncpc+nf
  do ib=-nf1,ncpb+nf
  do iov=1,nf
     v_cc( -nf+iov,ib,ic)=vrecva2(iov,ib,ic)
     v_cc(ncpa+iov,ib,ic)=vrecva3(iov,ib,ic)
  end do
  end do
  end do
  end if
  if (myra .eq. 0) then
  do ic=-nf1,ncpc+nf
  do ib=-nf1,ncpb+nf
  do iov=1,nf
     v_cc(ncpa+iov,ib,ic)=vrecva3(iov,ib,ic)
  end do
  end do
  end do
  end if
  if (myra .eq. nproca-1) then
  do ic=-nf1,ncpc+nf
  do ib=-nf1,ncpb+nf
  do iov=1,nf
     v_cc( -nf+iov,ib,ic)=vrecva2(iov,ib,ic)
  end do
  end do
  end do
  end if
  end if
! **************************************************
  else
! **********  communicate in x direction  **********
  do ic=-nf1,ncpc+nf
  do ib=-nf1,ncpb+nf
  do iov=1,nf
     vsenda2(iov,ib,ic)=v_cc(ncpa-nf+iov,ib,ic)
     vrecva2(iov,ib,ic)=v_cc(ncpa-nf+iov,ib,ic)
     vsenda3(iov,ib,ic)=v_cc(        iov,ib,ic)
     vrecva3(iov,ib,ic)=v_cc(        iov,ib,ic)
  end do
  end do
  end do

  stime=mpi_wtime()
  if (nproca .ne. 1) then
  weightcom=4.0d0*nf
  comcount=comcount+2.0d0*nf*(ncpb+2*nf)*(ncpc+2*nf)*weightcom
  vrecva2=0.0d0
  vrecva3=0.0d0
  idests0=myr_space+1
  idests1=myr_space-1
  idestr0=myr_space-1
  idestr1=myr_space+1
  if (myra .eq. nproca-1) idests0=myr_space-(nproca-1)
  if (myra .eq.        0) idests1=myr_space+(nproca-1)
  if (myra .eq.        0) idestr0=myr_space+(nproca-1)
  if (myra .eq. nproca-1) idestr1=myr_space-(nproca-1)
  if (mod(myra,2) .eq. 0) then
    call mpi_isend(vsenda2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests0,0,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr1,1,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
    call mpi_isend(vsenda3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests1,2,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr0,3,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if (mod(myra,2) .eq. 1) then
    call mpi_isend(vsenda3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests1,1,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr0,0,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
    call mpi_isend(vsenda2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests0,3,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr1,2,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  end if
  etime=mpi_wtime()
  comtime=comtime+(etime-stime)

  do ic=-nf1,ncpc+nf
  do ib=-nf1,ncpb+nf
  do iov=1,nf
     v_cc( -nf+iov,ib,ic)=vrecva2(iov,ib,ic)
     v_cc(ncpa+iov,ib,ic)=vrecva3(iov,ib,ic)
  end do
  end do
  end do
! **************************************************
  end if

  if (nperi .le. 1) then
! **********  communicate in y direction  **********
  if (nprocb .ne. 1) then
  do ic=-nf1,ncpc+nf
  do iov=1,nf
  do ia=-nf1,ncpa+nf
     vsendb2(ia,iov,ic)=v_cc(ia,ncpb-nf+iov,ic)
     vsendb3(ia,iov,ic)=v_cc(ia,        iov,ic)
  end do
  end do
  end do
  end if

  stime=mpi_wtime()
  if (nprocb .ne. 1) then
  weightcom=2.0d0*nf
  if (nprocb .ge. 4) weightcom=4.0d0*nf
  comcount=comcount+2.0d0*(ncpa+2*nf)*nf*(ncpc+2*nf)*weightcom
  vrecvb2=0.0d0
  vrecvb3=0.0d0
  idests0=myr_space+nproca
  idests1=myr_space-nproca
  idestr0=myr_space-nproca
  idestr1=myr_space+nproca
  if (myrb .eq. nprocb-1) idests0=myr_space-(nprocb-1)*nproca
  if (myrb .eq.        0) idests1=myr_space+(nprocb-1)*nproca
  if (myrb .eq.        0) idestr0=myr_space+(nprocb-1)*nproca
  if (myrb .eq. nprocb-1) idestr1=myr_space-(nprocb-1)*nproca
  if (mod(myrb,2) .eq. 0) then
    call mpi_isend(vsendb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests0,4,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr1,5,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if (mod(myrb,2) .eq. 1) then
    call mpi_isend(vsendb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests1,5,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr0,4,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if ((mod(myrb,2) .eq. 0) .and. (myrb .ne. 0)) then
    call mpi_isend(vsendb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests1,6,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr0,7,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if ((mod(myrb,2) .eq. 1) .and. (myrb .ne. nprocb-1)) then
    call mpi_isend(vsendb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests0,7,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr1,6,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  end if
  etime=mpi_wtime()
  comtime=comtime+(etime-stime)

  if (nprocb .ne. 1) then
  if ((myrb .ne. 0) .and. (myrb .ne. nprocb-1)) then
  do ic=-nf1,ncpc+nf
  do iov=1,nf
  do ia=-nf1,ncpa+nf
     v_cc(ia, -nf+iov,ic)=vrecvb2(ia,iov,ic)
     v_cc(ia,ncpb+iov,ic)=vrecvb3(ia,iov,ic)
  end do
  end do
  end do
  end if
  if (myrb .eq. 0) then
  do ic=-nf1,ncpc+nf
  do iov=1,nf
  do ia=-nf1,ncpa+nf
     v_cc(ia,ncpb+iov,ic)=vrecvb3(ia,iov,ic)
  end do
  end do
  end do
  end if
  if (myrb .eq. nprocb-1) then
  do ic=-nf1,ncpc+nf
  do iov=1,nf
  do ia=-nf1,ncpa+nf
     v_cc(ia, -nf+iov,ic)=vrecvb2(ia,iov,ic)
  end do
  end do
  end do
  end if
  end if
! **************************************************
  else
! **********  communicate in y direction  **********
  do ic=-nf1,ncpc+nf
  do iov=1,nf
  do ia=-nf1,ncpa+nf
     vsendb2(ia,iov,ic)=v_cc(ia,ncpb-nf+iov,ic)
     vrecvb2(ia,iov,ic)=v_cc(ia,ncpb-nf+iov,ic)
     vsendb3(ia,iov,ic)=v_cc(ia,        iov,ic)
     vrecvb3(ia,iov,ic)=v_cc(ia,        iov,ic)
  end do
  end do
  end do

  stime=mpi_wtime()
  if (nprocb .ne. 1) then
  weightcom=4.0d0*nf
  comcount=comcount+2.0d0*(ncpa+2*nf)*nf*(ncpc+2*nf)*weightcom
  vrecvb2=0.0d0
  vrecvb3=0.0d0
  idests0=myr_space+nproca
  idests1=myr_space-nproca
  idestr0=myr_space-nproca
  idestr1=myr_space+nproca
  if (myrb .eq. nprocb-1) idests0=myr_space-(nprocb-1)*nproca
  if (myrb .eq.        0) idests1=myr_space+(nprocb-1)*nproca
  if (myrb .eq.        0) idestr0=myr_space+(nprocb-1)*nproca
  if (myrb .eq. nprocb-1) idestr1=myr_space-(nprocb-1)*nproca
  if (mod(myrb,2) .eq. 0) then
    call mpi_isend(vsendb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests0,4,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr1,5,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
    call mpi_isend(vsendb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests1,6,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr0,7,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if (mod(myrb,2) .eq. 1) then
    call mpi_isend(vsendb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests1,5,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr0,4,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
    call mpi_isend(vsendb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests0,7,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr1,6,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  end if
  etime=mpi_wtime()
  comtime=comtime+(etime-stime)

  do ic=-nf1,ncpc+nf
  do iov=1,nf
  do ia=-nf1,ncpa+nf
     v_cc(ia, -nf+iov,ic)=vrecvb2(ia,iov,ic)
     v_cc(ia,ncpb+iov,ic)=vrecvb3(ia,iov,ic)
  end do
  end do
  end do
! **************************************************
  end if

  if (nperi .le. 2) then
! **********  communicate in z direction  **********
  if (nprocc .ne. 1) then
  do iov=1,nf
  do ib=-nf1,ncpb+nf
  do ia=-nf1,ncpa+nf
     vsendc2(ia,ib,iov)=v_cc(ia,ib,ncpc-nf+iov)
     vsendc3(ia,ib,iov)=v_cc(ia,ib,        iov)
  end do
  end do
  end do
  end if

  stime=mpi_wtime()
  if (nprocc .ne. 1) then
  weightcom=2.0d0*nf
  if (nprocc .ge. 4) weightcom=4.0d0*nf
  comcount=comcount+2.0d0*(ncpa+2*nf)*(ncpb+2*nf)*nf*weightcom
  vrecvc2=0.0d0
  vrecvc3=0.0d0
  idests0=myr_space+nproca*nprocb
  idests1=myr_space-nproca*nprocb
  idestr0=myr_space-nproca*nprocb
  idestr1=myr_space+nproca*nprocb
  if (myrc .eq. nprocc-1) idests0=myr_space-(nprocc-1)*nproca*nprocb
  if (myrc .eq.        0) idests1=myr_space+(nprocc-1)*nproca*nprocb
  if (myrc .eq.        0) idestr0=myr_space+(nprocc-1)*nproca*nprocb
  if (myrc .eq. nprocc-1) idestr1=myr_space-(nprocc-1)*nproca*nprocb
  if (mod(myrc,2) .eq. 0) then
    call mpi_isend(vsendc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests0,8,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr1,9,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if (mod(myrc,2) .eq. 1) then
    call mpi_isend(vsendc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests1,9,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr0,8,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if ((mod(myrc,2) .eq. 0) .and. (myrc .ne. 0)) then
    call mpi_isend(vsendc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests1,10,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr0,11,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if ((mod(myrc,2) .eq. 1) .and. (myrc .ne. nprocc-1)) then
    call mpi_isend(vsendc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests0,11,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr1,10,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  end if
  etime=mpi_wtime()
  comtime=comtime+(etime-stime)

  if (nprocc .ne. 1) then
  if ((myrc .ne. 0) .and. (myrc .ne. nprocc-1)) then
  do iov=1,nf
  do ib=-nf1,ncpb+nf
  do ia=-nf1,ncpa+nf
     v_cc(ia,ib, -nf+iov)=vrecvc2(ia,ib,iov)
     v_cc(ia,ib,ncpc+iov)=vrecvc3(ia,ib,iov)
  end do
  end do
  end do
  end if
  if (myrc .eq. 0) then
  do iov=1,nf
  do ib=-nf1,ncpb+nf
  do ia=-nf1,ncpa+nf
     v_cc(ia,ib,ncpc+iov)=vrecvc3(ia,ib,iov)
  end do
  end do
  end do
  end if
  if (myrc .eq. nprocc-1) then
  do iov=1,nf
  do ib=-nf1,ncpb+nf
  do ia=-nf1,ncpa+nf
     v_cc(ia,ib, -nf+iov)=vrecvc2(ia,ib,iov)
  end do
  end do
  end do
  end if
  end if
! **************************************************
  else
! **********  communicate in z direction  **********
  do iov=1,nf
  do ib=-nf1,ncpb+nf
  do ia=-nf1,ncpa+nf
     vsendc2(ia,ib,iov)=v_cc(ia,ib,ncpc-nf+iov)
     vrecvc2(ia,ib,iov)=v_cc(ia,ib,ncpc-nf+iov)
     vsendc3(ia,ib,iov)=v_cc(ia,ib,        iov)
     vrecvc3(ia,ib,iov)=v_cc(ia,ib,        iov)
  end do
  end do
  end do

  stime=mpi_wtime()
  if (nprocc .ne. 1) then
  weightcom=4.0d0*nf
  comcount=comcount+2.0d0*(ncpa+2*nf)*(ncpb+2*nf)*nf*weightcom
  vrecvc2=0.0d0


  vrecvc3=0.0d0
  idests0=myr_space+nproca*nprocb
  idests1=myr_space-nproca*nprocb
  idestr0=myr_space-nproca*nprocb
  idestr1=myr_space+nproca*nprocb
  if (myrc .eq. nprocc-1) idests0=myr_space-(nprocc-1)*nproca*nprocb
  if (myrc .eq.        0) idests1=myr_space+(nprocc-1)*nproca*nprocb
  if (myrc .eq.        0) idestr0=myr_space+(nprocc-1)*nproca*nprocb
  if (myrc .eq. nprocc-1) idestr1=myr_space-(nprocc-1)*nproca*nprocb
  if (mod(myrc,2) .eq. 0) then
    call mpi_isend(vsendc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests0,8,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr1,9,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
    call mpi_isend(vsendc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests1,10,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr0,11,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if (mod(myrc,2) .eq. 1) then
    call mpi_isend(vsendc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests1,9,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr0,8,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
    call mpi_isend(vsendc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests0,11,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr1,10,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  end if
  etime=mpi_wtime()
  comtime=comtime+(etime-stime)

  do iov=1,nf
  do ib=-nf1,ncpb+nf
  do ia=-nf1,ncpa+nf
     v_cc(ia,ib, -nf+iov)=vrecvc2(ia,ib,iov)
     v_cc(ia,ib,ncpc+iov)=vrecvc3(ia,ib,iov)
  end do
  end do
  end do
! **************************************************
  end if
  end if

  if (nflag .eq. 2) then
  if (nperi .eq. 0) then
! **********  communicate in x direction  **********
  do ic=-nf+1,ncpc+nf
  do ib=-nf+1,ncpb+nf
  do iov=1,nf
    vsenda2(iov,ib,ic)=v_cc( -nf+iov,ib,ic)
    vsenda3(iov,ib,ic)=v_cc(ncpa+iov,ib,ic)
  end do
  end do
  end do

  stime=mpi_wtime()
  if (nproca .ne. 1) then
  weightcom=2.0d0*nf
  if (nproca .ge. 4) weightcom=4.0d0*nf
  comcount=comcount+2.0d0*nf*(ncpb+2*nf)*(ncpc+2*nf)*weightcom
  vrecva2=0.0d0
  vrecva3=0.0d0
  idests0=myr_space-1
  idests1=myr_space+1
  idestr0=myr_space+1
  idestr1=myr_space-1
  if (myra .eq.        0) idests0=myr_space+(nproca-1)
  if (myra .eq. nproca-1) idests1=myr_space-(nproca-1)
  if (myra .eq. nproca-1) idestr0=myr_space-(nproca-1)
  if (myra .eq.        0) idestr1=myr_space+(nproca-1)

  if ((mod(myra,2) .eq. 0) .and. (myra .ne. 0)) then
    call mpi_isend(vsenda2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests0,0,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr1,1,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if ((mod(myra,2) .eq. 1) .and. (myra .ne. nproca-1)) then
    call mpi_isend(vsenda3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests1,1,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr0,0,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if (mod(myra,2) .eq. 0) then
    call mpi_isend(vsenda3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests1,2,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr0,3,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if (mod(myra,2) .eq. 1) then
    call mpi_isend(vsenda2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests0,3,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr1,2,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  end if
  etime=mpi_wtime()
  comtime=comtime+(etime-stime)

  if (nproca .ne. 1) then
    do ic=-nf+1,ncpc+nf
    do ib=-nf+1,ncpb+nf
    do iov=1,nf
      v_cc(        iov,ib,ic)=v_cc(        iov,ib,ic)+vrecva3(iov,ib,ic)
      v_cc(ncpa-nf+iov,ib,ic)=v_cc(ncpa-nf+iov,ib,ic)+vrecva2(iov,ib,ic)
    end do
    end do
    end do
  end if
! **************************************************
  else
! **********  communicate in x direction  **********
  do ic=-nf+1,ncpc+nf
  do ib=-nf+1,ncpb+nf
  do iov=1,nf
     vsenda2(iov,ib,ic)=v_cc( -nf+iov,ib,ic)
     vsenda3(iov,ib,ic)=v_cc(ncpa+iov,ib,ic)
  end do
  end do
  end do

  stime=mpi_wtime()
  if (nproca .ne. 1) then
  weightcom=2.0d0*nf
  if (nproca .ge. 4) weightcom=4.0d0*nf
  comcount=comcount+2.0d0*nf*(ncpb+2*nf)*(ncpc+2*nf)*weightcom
  vrecva2=0.0d0
  vrecva3=0.0d0
  idests0=myr_space-1
  idests1=myr_space+1
  idestr0=myr_space+1
  idestr1=myr_space-1
  if (myra .eq.        0) idests0=myr_space+(nproca-1)
  if (myra .eq. nproca-1) idests1=myr_space-(nproca-1)
  if (myra .eq. nproca-1) idestr0=myr_space-(nproca-1)
  if (myra .eq.        0) idestr1=myr_space+(nproca-1)

  if (mod(myra,2) .eq. 0) then
    call mpi_isend(vsenda2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests0,0,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr1,1,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
    call mpi_isend(vsenda3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests1,2,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr0,3,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if (mod(myra,2) .eq. 1) then
    call mpi_isend(vsenda3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests1,1,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr0,0,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
    call mpi_isend(vsenda2,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idests0,3,mpicom_space,req1,mpij)
    call mpi_irecv(vrecva3,nf*(ncpb+2*nf)*(ncpc+2*nf),mpi_double_precision,idestr1,2,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  end if
  etime=mpi_wtime()
  comtime=comtime+(etime-stime)

  if (nproca .eq. 1) then
    vrecva3=vsenda3
    vrecva2=vsenda2
  end if
  do ic=-nf+1,ncpc+nf
  do ib=-nf+1,ncpb+nf
  do iov=1,nf
    v_cc(        iov,ib,ic)=v_cc(        iov,ib,ic)+vrecva3(iov,ib,ic)
    v_cc(ncpa-nf+iov,ib,ic)=v_cc(ncpa-nf+iov,ib,ic)+vrecva2(iov,ib,ic)
  end do
  end do
  end do
! **************************************************
  end if

  if (nperi .le. 1) then
! **********  communicate in y direction  **********
  do ic=-nf+1,ncpc+nf
  do iov=1,nf
  do ia=-nf+1,ncpa+nf
    vsendb2(ia,iov,ic)=v_cc(ia, -nf+iov,ic)
    vsendb3(ia,iov,ic)=v_cc(ia,ncpb+iov,ic)
  end do
  end do
  end do

  stime=mpi_wtime()
  if (nprocb .ne. 1) then
  weightcom=2.0d0*nf
  if (nprocb .ge. 4) weightcom=4.0d0*nf
  comcount=comcount+2.0d0*(ncpa+2*nf)*nf*(ncpc+2*nf)*weightcom
  vrecvb2=0.0d0
  vrecvb3=0.0d0
  idests0=myr_space-nproca
  idests1=myr_space+nproca
  idestr0=myr_space+nproca
  idestr1=myr_space-nproca
  if (myrb .eq.        0) idests0=myr_space+(nprocb-1)*nproca
  if (myrb .eq. nprocb-1) idests1=myr_space-(nprocb-1)*nproca
  if (myrb .eq. nprocb-1) idestr0=myr_space-(nprocb-1)*nproca
  if (myrb .eq.        0) idestr1=myr_space+(nprocb-1)*nproca
  if ((mod(myrb,2) .eq. 0) .and. (myrb .ne. 0)) then
    call mpi_isend(vsendb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests0,4,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr1,5,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if ((mod(myrb,2) .eq. 1) .and. (myrb .ne. nprocb-1)) then
    call mpi_isend(vsendb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests1,5,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr0,4,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if (mod(myrb,2) .eq. 0) then
    call mpi_isend(vsendb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests1,6,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr0,7,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if (mod(myrb,2) .eq. 1) then
    call mpi_isend(vsendb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests0,7,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr1,6,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  end if
  etime=mpi_wtime()
  comtime=comtime+(etime-stime)

  if (nprocb .ne. 1) then
    do ic=-nf+1,ncpc+nf
    do iov=1,nf
    do ia=-nf+1,ncpa+nf
      v_cc(ia,        iov,ic)=v_cc(ia,        iov,ic)+vrecvb3(ia,iov,ic)
      v_cc(ia,ncpb-nf+iov,ic)=v_cc(ia,ncpb-nf+iov,ic)+vrecvb2(ia,iov,ic)
    end do
    end do
    end do
  end if
! **************************************************
  else
! **********  communicate in y direction  **********
  do ic=-nf+1,ncpc+nf
  do iov=1,nf
  do ia=-nf+1,ncpa+nf
     vsendb2(ia,iov,ic)=v_cc(ia, -nf+iov,ic)
     vsendb3(ia,iov,ic)=v_cc(ia,ncpb+iov,ic)
  end do
  end do
  end do

  stime=mpi_wtime()
  if (nprocb .ne. 1) then
  weightcom=2.0d0*nf
  if (nprocb .ge. 4) weightcom=4.0d0*nf
  comcount=comcount+2.0d0*(ncpa+2*nf)*nf*(ncpc+2*nf)*weightcom
  vrecvb2=0.0d0
  vrecvb3=0.0d0
  idests0=myr_space-nproca
  idests1=myr_space+nproca
  idestr0=myr_space+nproca
  idestr1=myr_space-nproca
  if (myrb .eq.        0) idests0=myr_space+(nprocb-1)*nproca
  if (myrb .eq. nprocb-1) idests1=myr_space-(nprocb-1)*nproca
  if (myrb .eq. nprocb-1) idestr0=myr_space-(nprocb-1)*nproca
  if (myrb .eq.        0) idestr1=myr_space+(nprocb-1)*nproca
  if (mod(myrb,2) .eq. 0) then
    call mpi_isend(vsendb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests0,4,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr1,5,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
    call mpi_isend(vsendb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests1,6,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr0,7,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if (mod(myrb,2) .eq. 1) then
    call mpi_isend(vsendb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests1,5,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr0,4,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
    call mpi_isend(vsendb2,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idests0,7,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvb3,(ncpa+2*nf)*nf*(ncpc+2*nf),mpi_double_precision,idestr1,6,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  end if
  etime=mpi_wtime()
  comtime=comtime+(etime-stime)

  if (nprocb .eq. 1) then
    vrecvb3=vsendb3
    vrecvb2=vsendb2
  end if
  do ic=-nf+1,ncpc+nf
  do iov=1,nf
  do ia=-nf+1,ncpa+nf
     v_cc(ia,        iov,ic)=v_cc(ia,        iov,ic)+vrecvb3(ia,iov,ic)
     v_cc(ia,ncpb-nf+iov,ic)=v_cc(ia,ncpb-nf+iov,ic)+vrecvb2(ia,iov,ic)
  end do
  end do
  end do
! **************************************************
  end if

  if (nperi .le. 2) then
! **********  communicate in z direction  **********
  do iov=1,nf
  do ib=-nf+1,ncpb+nf
  do ia=-nf+1,ncpa+nf
    vsendc2(ia,ib,iov)=v_cc(ia,ib, -nf+iov)
    vsendc3(ia,ib,iov)=v_cc(ia,ib,ncpc+iov)
  end do
  end do
  end do

  stime=mpi_wtime()
  if (nprocc .ne. 1) then
  weightcom=2.0d0*nf
  if (nprocc .ge. 4) weightcom=4.0d0*nf
  comcount=comcount+2.0d0*(ncpa+2*nf)*(ncpb+2*nf)*nf*weightcom
  vrecvc2=0.0d0
  vrecvc3=0.0d0
  idests0=myr_space-nproca*nprocb
  idests1=myr_space+nproca*nprocb
  idestr0=myr_space+nproca*nprocb
  idestr1=myr_space-nproca*nprocb
  if (myrc .eq.        0) idests0=myr_space+(nprocc-1)*nproca*nprocb
  if (myrc .eq. nprocc-1) idests1=myr_space-(nprocc-1)*nproca*nprocb
  if (myrc .eq. nprocc-1) idestr0=myr_space-(nprocc-1)*nproca*nprocb
  if (myrc .eq.        0) idestr1=myr_space+(nprocc-1)*nproca*nprocb
  if ((mod(myrc,2) .eq. 0) .and. (myrc .ne. 0)) then
    call mpi_isend(vsendc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests0,8,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr1,9,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if ((mod(myrc,2) .eq. 1) .and. (myrc .ne. nprocc-1)) then
    call mpi_isend(vsendc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests1,9,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr0,8,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if (mod(myrc,2) .eq. 0) then
    call mpi_isend(vsendc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests1,10,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr0,11,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if (mod(myrc,2) .eq. 1) then
    call mpi_isend(vsendc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests0,11,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr1,10,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  end if
  etime=mpi_wtime()
  comtime=comtime+(etime-stime)

  if (nprocc .ne. 1) then
    do iov=1,nf
    do ib=-nf+1,ncpb+nf
    do ia=-nf+1,ncpa+nf
      v_cc(ia,ib,        iov)=v_cc(ia,ib,        iov)+vrecvc3(ia,ib,iov)
      v_cc(ia,ib,ncpc-nf+iov)=v_cc(ia,ib,ncpc-nf+iov)+vrecvc2(ia,ib,iov)
    end do
    end do
    end do
  end if
! **************************************************
  else
! **********  communicate in z direction  **********
  do iov=1,nf
  do ib=-nf+1,ncpb+nf
  do ia=-nf+1,ncpa+nf
    vsendc2(ia,ib,iov)=v_cc(ia,ib, -nf+iov)
    vsendc3(ia,ib,iov)=v_cc(ia,ib,ncpc+iov)
  end do
  end do
  end do

  stime=mpi_wtime()
  if (nprocc .ne. 1) then
  weightcom=2.0d0*nf
  if (nprocc .ge. 4) weightcom=4.0d0*nf
  comcount=comcount+2.0d0*(ncpa+2*nf)*(ncpb+2*nf)*nf*weightcom
  vrecvc2=0.0d0
  vrecvc3=0.0d0
  idests0=myr_space-nproca*nprocb
  idests1=myr_space+nproca*nprocb
  idestr0=myr_space+nproca*nprocb
  idestr1=myr_space-nproca*nprocb
  if (myrc .eq.        0) idests0=myr_space+(nprocc-1)*nproca*nprocb
  if (myrc .eq. nprocc-1) idests1=myr_space-(nprocc-1)*nproca*nprocb
  if (myrc .eq. nprocc-1) idestr0=myr_space-(nprocc-1)*nproca*nprocb
  if (myrc .eq.        0) idestr1=myr_space+(nprocc-1)*nproca*nprocb
  if (mod(myrc,2) .eq. 0) then
    call mpi_isend(vsendc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests0,8,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr1,9,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
    call mpi_isend(vsendc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests1,10,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr0,11,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  if (mod(myrc,2) .eq. 1) then
    call mpi_isend(vsendc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests1,9,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr0,8,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
    call mpi_isend(vsendc2,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idests0,11,mpicom_space,req1,mpij)
    call mpi_irecv(vrecvc3,(ncpa+2*nf)*(ncpb+2*nf)*nf,mpi_double_precision,idestr1,10,mpicom_space,req2,mpij)
    call mpi_wait(req1,mpistat,mpij)
    call mpi_wait(req2,mpistat,mpij)
  end if
  end if
  etime=mpi_wtime()
  comtime=comtime+(etime-stime)

  if (nprocc .eq. 1) then
    vrecvc3=vsendc3
    vrecvc2=vsendc2
  end if
  do iov=1,nf
  do ib=-nf+1,ncpb+nf
  do ia=-nf+1,ncpa+nf
    v_cc(ia,ib,        iov)=v_cc(ia,ib,        iov)+vrecvc3(ia,ib,iov)
    v_cc(ia,ib,ncpc-nf+iov)=v_cc(ia,ib,ncpc-nf+iov)+vrecvc2(ia,ib,iov)
  end do
  end do
  end do
! **************************************************
  end if
  end if

  deallocate(vsenda2,vsenda3,vrecva2,vrecva3,vsendb2,vsendb3,vrecvb2,vrecvb3,vsendc2,vsendc3,vrecvc2,vrecvc3)

return
end subroutine


end module
