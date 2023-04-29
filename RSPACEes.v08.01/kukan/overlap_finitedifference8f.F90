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
!     **********  overlap_finitedifference8f.F90 12/08/2016-01  **********

      module mod_overlap_finitedifference
      implicit none
      real*8,allocatable::&
             vsenda0(:),vsenda1(:) &
            ,vrecva0(:),vrecva1(:) &
            ,vsendb0(:),vsendb1(:) &
            ,vrecvb0(:),vrecvb1(:) &
            ,vsendc0(:),vsendc1(:) &
            ,vrecvc0(:),vrecvc1(:)
      integer ireqa1,ireqa2,ireqa3,ireqa4,ireqa5,ireqa6,ireqa7,ireqa8
      integer ireqb1,ireqb2,ireqb3,ireqb4,ireqb5,ireqb6,ireqb7,ireqb8
      integer ireqc1,ireqc2,ireqc3,ireqc4,ireqc5,ireqc6,ireqc7,ireqc8
      real*8 stime,etime

      contains

      subroutine overlap_finitedifference_init(ncpa,ncpb,ncpc,nf,ncol,nrc)
      implicit none
      integer ncpa,ncpb,ncpc,nf,ncol,nrc
      allocate(vsenda0(nf*ncpb*ncpc*ncol*(nrc+1)),vsenda1(nf*ncpb*ncpc*ncol*(nrc+1)) &
              ,vrecva0(nf*ncpb*ncpc*ncol*(nrc+1)),vrecva1(nf*ncpb*ncpc*ncol*(nrc+1)) &
              ,vsendb0(ncpa*nf*ncpc*ncol*(nrc+1)),vsendb1(ncpa*nf*ncpc*ncol*(nrc+1)) &
              ,vrecvb0(ncpa*nf*ncpc*ncol*(nrc+1)),vrecvb1(ncpa*nf*ncpc*ncol*(nrc+1)) &
              ,vsendc0(ncpa*ncpb*nf*ncol*(nrc+1)),vsendc1(ncpa*ncpb*nf*ncol*(nrc+1)) &
              ,vrecvc0(ncpa*ncpb*nf*ncol*(nrc+1)),vrecvc1(ncpa*ncpb*nf*ncol*(nrc+1)))
      end subroutine


      subroutine overlap_finitedifference_final
      deallocate(vsenda0,vsenda1,vrecva0,vrecva1,vsendb0,vsendb1,vrecvb0,vrecvb1,vsendc0,vsendc1,vrecvc0,vrecvc1)
      end subroutine


      subroutine overlap_finitedifference_r(nperi,ncpa,ncpb,ncpc,nf,nf1,nsnd,v_cc)
      use mod_mpi,myra=>myrx,myrb=>myry,myrc=>myrz,nproca=>nprocx,nprocb=>nprocy,nprocc=>nprocz
      implicit none
      integer nperi,ncpa,ncpb,ncpc,nf,nf1,nsnd
      real*8 v_cc(-nf1:ncpa+nf,-nf1:ncpb+nf,-nf1:ncpc+nf)
      integer ia,ib,ic,i
      integer iov
      integer idests0,idests1,idestr0,idestr1

      real*8 weightcom

      stime=mpi_wtime()

      if (nperi .eq. 0) then
!     **********  communicate in a direction  **********
      if (nproca .ne. 1) then
      i=0
      do ic=1,ncpc
      do ib=1,ncpb
      do iov=1,nsnd
         i=i+1
         vsenda0(i)=v_cc(ncpa-nsnd+iov,ib,ic)
         vsenda1(i)=v_cc(          iov,ib,ic)
      end do
      end do
      end do
      end if

      if (nproca .ne. 1) then
      weightcom=2.0d0*nsnd
      if (nproca .ge. 4) weightcom=4.0d0*nsnd
      comcount=comcount+2.0d0*nsnd*ncpb*ncpc*weightcom
      vrecva0=0.0d0
      vrecva1=0.0d0
      idests0=myr_space+1
      idests1=myr_space-1
      idestr0=myr_space-1
      idestr1=myr_space+1
      if (myra .eq. nproca-1) idests0=idests0-nproca
      if (myra .eq.        0) idests1=idests1+nproca
      if (myra .eq.        0) idestr0=idestr0+nproca
      if (myra .eq. nproca-1) idestr1=idestr1-nproca

      if (mod(myra,2) .eq. 0) then
      call mpi_isend(vsenda0,nsnd*ncpb*ncpc,mpi_double_precision,idests0,0,mpicom_space,ireqa1,mpij)
      call mpi_irecv(vrecva1,nsnd*ncpb*ncpc,mpi_double_precision,idestr1,1,mpicom_space,ireqa2,mpij)
      end if
      if (mod(myra,2) .eq. 1) then
      call mpi_isend(vsenda1,nsnd*ncpb*ncpc,mpi_double_precision,idests1,1,mpicom_space,ireqa3,mpij)
      call mpi_irecv(vrecva0,nsnd*ncpb*ncpc,mpi_double_precision,idestr0,0,mpicom_space,ireqa4,mpij)
      end if
      if ((mod(myra,2) .eq. 0) .and. (myra .ne. 0)) then
      call mpi_isend(vsenda1,nsnd*ncpb*ncpc,mpi_double_precision,idests1,2,mpicom_space,ireqa5,mpij)
      call mpi_irecv(vrecva0,nsnd*ncpb*ncpc,mpi_double_precision,idestr0,3,mpicom_space,ireqa6,mpij)
      end if
      if ((mod(myra,2) .eq. 1) .and. (myra .ne. nproca-1)) then
      call mpi_isend(vsenda0,nsnd*ncpb*ncpc,mpi_double_precision,idests0,3,mpicom_space,ireqa7,mpij)
      call mpi_irecv(vrecva1,nsnd*ncpb*ncpc,mpi_double_precision,idestr1,2,mpicom_space,ireqa8,mpij)
      end if
      end if
!     **************************************************
      else
!     **********  communicate in a direction  **********
      i=0
      do ic=1,ncpc
      do ib=1,ncpb
      do iov=1,nsnd
         i=i+1
         vsenda0(i)=v_cc(ncpa-nsnd+iov,ib,ic)
         vrecva0(i)=v_cc(ncpa-nsnd+iov,ib,ic)
         vsenda1(i)=v_cc(          iov,ib,ic)
         vrecva1(i)=v_cc(          iov,ib,ic)
      end do
      end do
      end do

      if (nproca .ne. 1) then
      weightcom=4.0d0*nsnd
      comcount=comcount+2.0d0*nsnd*ncpc*ncpb*weightcom
      vrecva0=0.0d0
      vrecva1=0.0d0
      idests0=myr_space+1
      idests1=myr_space-1
      idestr0=myr_space-1
      idestr1=myr_space+1
      if (myra .eq. nproca-1) idests0=idests0-nproca
      if (myra .eq.        0) idests1=idests1+nproca
      if (myra .eq.        0) idestr0=idestr0+nproca
      if (myra .eq. nproca-1) idestr1=idestr1-nproca
      if (mod(myra,2) .eq. 0) then
      call mpi_isend(vsenda0,nsnd*ncpb*ncpc,mpi_double_precision,idests0,0,mpicom_space,ireqa1,mpij)
      call mpi_irecv(vrecva1,nsnd*ncpb*ncpc,mpi_double_precision,idestr1,1,mpicom_space,ireqa2,mpij)
      call mpi_isend(vsenda1,nsnd*ncpb*ncpc,mpi_double_precision,idests1,2,mpicom_space,ireqa3,mpij)
      call mpi_irecv(vrecva0,nsnd*ncpb*ncpc,mpi_double_precision,idestr0,3,mpicom_space,ireqa4,mpij)
      end if
      if (mod(myra,2) .eq. 1) then
      call mpi_isend(vsenda1,nsnd*ncpb*ncpc,mpi_double_precision,idests1,1,mpicom_space,ireqa5,mpij)
      call mpi_irecv(vrecva0,nsnd*ncpb*ncpc,mpi_double_precision,idestr0,0,mpicom_space,ireqa6,mpij)
      call mpi_isend(vsenda0,nsnd*ncpb*ncpc,mpi_double_precision,idests0,3,mpicom_space,ireqa7,mpij)
      call mpi_irecv(vrecva1,nsnd*ncpb*ncpc,mpi_double_precision,idestr1,2,mpicom_space,ireqa8,mpij)
      end if
      end if
!     **************************************************
      end if

      if (nperi .le. 1) then
!     **********  communicate in b direction  **********
      if (nprocb .ne. 1) then
      i=0
      do ic=1,ncpc
      do iov=1,nsnd
      do ia=1,ncpa
         i=i+1
         vsendb0(i)=v_cc(ia,ncpb-nsnd+iov,ic)
         vsendb1(i)=v_cc(ia,          iov,ic)
      end do
      end do
      end do
      end if

      if (nprocb .ne. 1) then
      weightcom=2.0d0*nsnd
      if (nprocb .ge. 4) weightcom=4.0d0*nsnd
      comcount=comcount+2.0d0*nsnd*ncpa*ncpb*weightcom
      vrecvb0=0.0d0
      vrecvb1=0.0d0
      idests0=myr_space+nproca
      idests1=myr_space-nproca
      idestr0=myr_space-nproca
      idestr1=myr_space+nproca
      if (myrb .eq. nprocb-1) idests0=idests0-nprocb*nproca
      if (myrb .eq.        0) idests1=idests1+nprocb*nproca
      if (myrb .eq.        0) idestr0=idestr0+nprocb*nproca
      if (myrb .eq. nprocb-1) idestr1=idestr1-nprocb*nproca
      if (mod(myrb,2) .eq. 0) then
      call mpi_isend(vsendb0,ncpa*nsnd*ncpc,mpi_double_precision,idests0,4,mpicom_space,ireqb1,mpij)
      call mpi_irecv(vrecvb1,ncpa*nsnd*ncpc,mpi_double_precision,idestr1,5,mpicom_space,ireqb2,mpij)
      end if
      if (mod(myrb,2) .eq. 1) then
      call mpi_isend(vsendb1,ncpa*nsnd*ncpc,mpi_double_precision,idests1,5,mpicom_space,ireqb3,mpij)
      call mpi_irecv(vrecvb0,ncpa*nsnd*ncpc,mpi_double_precision,idestr0,4,mpicom_space,ireqb4,mpij)
      end if
      if ((mod(myrb,2) .eq. 0) .and. (myrb .ne. 0)) then
      call mpi_isend(vsendb1,ncpa*nsnd*ncpc,mpi_double_precision,idests1,6,mpicom_space,ireqb5,mpij)
      call mpi_irecv(vrecvb0,ncpa*nsnd*ncpc,mpi_double_precision,idestr0,7,mpicom_space,ireqb6,mpij)
      end if
      if ((mod(myrb,2) .eq. 1) .and. (myrb .ne. nprocb-1)) then
      call mpi_isend(vsendb0,ncpa*nsnd*ncpc,mpi_double_precision,idests0,7,mpicom_space,ireqb7,mpij)
      call mpi_irecv(vrecvb1,ncpa*nsnd*ncpc,mpi_double_precision,idestr1,6,mpicom_space,ireqb8,mpij)
      end if
      end if
!     **************************************************
      else
!     **********  communicate in b direction  **********
      i=0
      do ic=1,ncpc
      do iov=1,nsnd
      do ia=1,ncpa
         i=i+1
         vsendb0(i)=v_cc(ia,ncpb-nsnd+iov,ic)
         vrecvb0(i)=v_cc(ia,ncpb-nsnd+iov,ic)
         vsendb1(i)=v_cc(ia,          iov,ic)
         vrecvb1(i)=v_cc(ia,          iov,ic)
      end do
      end do
      end do

      if (nprocb .ne. 1) then
      weightcom=4.0d0*nsnd
      comcount=comcount+2.0d0*nsnd*ncpc*ncpa*weightcom
      vrecvb0=0.0d0
      vrecvb1=0.0d0
      idests0=myr_space+nproca
      idests1=myr_space-nproca
      idestr0=myr_space-nproca
      idestr1=myr_space+nproca
      if (myrb .eq. nprocb-1) idests0=idests0-nprocb*nproca
      if (myrb .eq.        0) idests1=idests1+nprocb*nproca
      if (myrb .eq.        0) idestr0=idestr0+nprocb*nproca
      if (myrb .eq. nprocb-1) idestr1=idestr1-nprocb*nproca
      if (mod(myrb,2) .eq. 0) then
      call mpi_isend(vsendb0,ncpa*nsnd*ncpc,mpi_double_precision,idests0,4,mpicom_space,ireqb1,mpij)
      call mpi_irecv(vrecvb1,ncpa*nsnd*ncpc,mpi_double_precision,idestr1,5,mpicom_space,ireqb2,mpij)
      call mpi_isend(vsendb1,ncpa*nsnd*ncpc,mpi_double_precision,idests1,6,mpicom_space,ireqb3,mpij)
      call mpi_irecv(vrecvb0,ncpa*nsnd*ncpc,mpi_double_precision,idestr0,7,mpicom_space,ireqb4,mpij)
      end if
      if (mod(myrb,2) .eq. 1) then
      call mpi_isend(vsendb1,ncpa*nsnd*ncpc,mpi_double_precision,idests1,5,mpicom_space,ireqb5,mpij)
      call mpi_irecv(vrecvb0,ncpa*nsnd*ncpc,mpi_double_precision,idestr0,4,mpicom_space,ireqb6,mpij)
      call mpi_isend(vsendb0,ncpa*nsnd*ncpc,mpi_double_precision,idests0,7,mpicom_space,ireqb7,mpij)
      call mpi_irecv(vrecvb1,ncpa*nsnd*ncpc,mpi_double_precision,idestr1,6,mpicom_space,ireqb8,mpij)
      end if
      end if
!     **************************************************
      end if

      if (nperi .le. 2) then
!     **********  communicate in c direction  **********
      if (nprocc .ne. 1) then
      i=0
      do iov=1,nsnd
      do ib=1,ncpb
      do ia=1,ncpa
         i=i+1
         vsendc0(i)=v_cc(ia,ib,ncpc-nsnd+iov)
         vsendc1(i)=v_cc(ia,ib,          iov)
      end do
      end do
      end do
      end if

      if (nprocc .ne. 1) then
      weightcom=2.0d0*nsnd
      if (nprocc .ge. 4) weightcom=4.0d0*nsnd
      comcount=comcount+2.0d0*nsnd*ncpb*ncpa*weightcom
      vrecvc0=0.0d0
      vrecvc1=0.0d0
      idests0=myr_space+nproca*nprocb
      idests1=myr_space-nproca*nprocb
      idestr0=myr_space-nproca*nprocb
      idestr1=myr_space+nproca*nprocb
      if (myrc .eq. nprocc-1) idests0=idests0-nprocc*nproca*nprocb
      if (myrc .eq.        0) idests1=idests1+nprocc*nproca*nprocb
      if (myrc .eq.        0) idestr0=idestr0+nprocc*nproca*nprocb
      if (myrc .eq. nprocc-1) idestr1=idestr1-nprocc*nproca*nprocb
      if (mod(myrc,2) .eq. 0) then
      call mpi_isend(vsendc0,ncpa*ncpb*nsnd,mpi_double_precision,idests0,8,mpicom_space,ireqc1,mpij)
      call mpi_irecv(vrecvc1,ncpa*ncpb*nsnd,mpi_double_precision,idestr1,9,mpicom_space,ireqc2,mpij)
      end if
      if (mod(myrc,2) .eq. 1) then
      call mpi_isend(vsendc1,ncpa*ncpb*nsnd,mpi_double_precision,idests1,9,mpicom_space,ireqc3,mpij)
      call mpi_irecv(vrecvc0,ncpa*ncpb*nsnd,mpi_double_precision,idestr0,8,mpicom_space,ireqc4,mpij)
      end if
      if ((mod(myrc,2) .eq. 0) .and. (myrc .ne. 0)) then
      call mpi_isend(vsendc1,ncpa*ncpb*nsnd,mpi_double_precision,idests1,10,mpicom_space,ireqc5,mpij)
      call mpi_irecv(vrecvc0,ncpa*ncpb*nsnd,mpi_double_precision,idestr0,11,mpicom_space,ireqc6,mpij)
      end if
      if ((mod(myrc,2) .eq. 1) .and. (myrc .ne. nprocc-1)) then
      call mpi_isend(vsendc0,ncpa*ncpb*nsnd,mpi_double_precision,idests0,11,mpicom_space,ireqc7,mpij)
      call mpi_irecv(vrecvc1,ncpa*ncpb*nsnd,mpi_double_precision,idestr1,10,mpicom_space,ireqc8,mpij)
      end if
      end if
!     **************************************************
      else
!     **********  communicate in c direction  **********
      i=0
      do iov=1,nsnd
      do ib=1,ncpb
      do ia=1,ncpa
         i=i+1
         vsendc0(i)=v_cc(ia,ib,ncpc-nsnd+iov)
         vrecvc0(i)=v_cc(ia,ib,ncpc-nsnd+iov)
         vsendc1(i)=v_cc(ia,ib,          iov)
         vrecvc1(i)=v_cc(ia,ib,          iov)
      end do
      end do
      end do

      if (nprocc .ne. 1) then
      weightcom=4.0d0*nsnd
      comcount=comcount+2.0d0*nsnd*ncpb*ncpa*weightcom
      vrecvc0=0.0d0
      vrecvc1=0.0d0
      idests0=myr_space+nproca*nprocb
      idests1=myr_space-nproca*nprocb
      idestr0=myr_space-nproca*nprocb
      idestr1=myr_space+nproca*nprocb
      if (myrc .eq. nprocc-1) idests0=idests0-nprocc*nproca*nprocb
      if (myrc .eq.        0) idests1=idests1+nprocc*nproca*nprocb
      if (myrc .eq.        0) idestr0=idestr0+nprocc*nproca*nprocb
      if (myrc .eq. nprocc-1) idestr1=idestr1-nprocc*nproca*nprocb
      if (mod(myrc,2) .eq. 0) then
      call mpi_isend(vsendc0,ncpa*ncpb*nsnd,mpi_double_precision,idests0,8,mpicom_space,ireqc1,mpij)
      call mpi_irecv(vrecvc1,ncpa*ncpb*nsnd,mpi_double_precision,idestr1,9,mpicom_space,ireqc2,mpij)
      call mpi_isend(vsendc1,ncpa*ncpb*nsnd,mpi_double_precision,idests1,10,mpicom_space,ireqc3,mpij)
      call mpi_irecv(vrecvc0,ncpa*ncpb*nsnd,mpi_double_precision,idestr0,11,mpicom_space,ireqc4,mpij)
      end if
      if (mod(myrc,2) .eq. 1) then
      call mpi_isend(vsendc1,ncpa*ncpb*nsnd,mpi_double_precision,idests1,9,mpicom_space,ireqc5,mpij)
      call mpi_irecv(vrecvc0,ncpa*ncpb*nsnd,mpi_double_precision,idestr0,8,mpicom_space,ireqc6,mpij)
      call mpi_isend(vsendc0,ncpa*ncpb*nsnd,mpi_double_precision,idests0,11,mpicom_space,ireqc7,mpij)
      call mpi_irecv(vrecvc1,ncpa*ncpb*nsnd,mpi_double_precision,idestr1,10,mpicom_space,ireqc8,mpij)
      end if
      end if
!     **************************************************
      end if

      return
      end subroutine


      subroutine overlap_finitedifference_c(nperi,ncpa,ncpb,ncpc,nf,nf1,nsnd,ncol,v_cc)
      use mod_mpi,myra=>myrx,myrb=>myry,myrc=>myrz,nproca=>nprocx,nprocb=>nprocy,nprocc=>nprocz
      implicit none
      integer nperi,ncpa,ncpb,ncpc,nf,nf1,nsnd,ncol
      complex*16 v_cc(-nf1:ncpa+nf,-nf1:ncpb+nf,-nf1:ncpc+nf,ncol)
      integer ia,ib,ic,is,i
      integer iov
      integer idests0,idests1,idests2,idests3,idestr0,idestr1,idestr2,idestr3

      real*8 weightcom

      stime=mpi_wtime()

      if (nperi .eq. 0) then
!     **********  communicate in a direction  **********
      if (nproca .ne. 1) then
      i=0
      do is=1,ncol
      do ic=1,ncpc
      do ib=1,ncpb
      do iov=1,nsnd
         i=i+1
         vsenda0(i                    )=dreal(v_cc(ncpa-nsnd+iov,ib,ic,is))
         vsenda1(i                    )=dreal(v_cc(          iov,ib,ic,is))
         vsenda0(i+nsnd*ncpb*ncpc*ncol)=dimag(v_cc(ncpa-nsnd+iov,ib,ic,is))
         vsenda1(i+nsnd*ncpb*ncpc*ncol)=dimag(v_cc(          iov,ib,ic,is))
      end do
      end do
      end do
      end do
      end if

      if (nproca .ne. 1) then
      weightcom=2.0d0*nsnd
      if (nproca .ge. 4) weightcom=4.0d0*nsnd
      comcount=comcount+2.0d0*nsnd*ncpb*ncpc*weightcom*2*ncol
      vrecva0=0.0d0
      vrecva1=0.0d0
      idests0=myr_space+1
      idests1=myr_space-1
      idestr0=myr_space-1
      idestr1=myr_space+1
      if (myra .eq. nproca-1) idests0=idests0-nproca
      if (myra .eq.        0) idests1=idests1+nproca
      if (myra .eq.        0) idestr0=idestr0+nproca
      if (myra .eq. nproca-1) idestr1=idestr1-nproca

      if (mod(myra,2) .eq. 0) then
      call mpi_isend(vsenda0,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idests0,0,mpicom_space,ireqa1,mpij)
      call mpi_irecv(vrecva1,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idestr1,1,mpicom_space,ireqa2,mpij)
      end if
      if (mod(myra,2) .eq. 1) then
      call mpi_isend(vsenda1,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idests1,1,mpicom_space,ireqa3,mpij)
      call mpi_irecv(vrecva0,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idestr0,0,mpicom_space,ireqa4,mpij)
      end if
      if ((mod(myra,2) .eq. 0) .and. (myra .ne. 0)) then
      call mpi_isend(vsenda1,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idests1,2,mpicom_space,ireqa5,mpij)
      call mpi_irecv(vrecva0,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idestr0,3,mpicom_space,ireqa6,mpij)
      end if
      if ((mod(myra,2) .eq. 1) .and. (myra .ne. nproca-1)) then
      call mpi_isend(vsenda0,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idests0,3,mpicom_space,ireqa7,mpij)
      call mpi_irecv(vrecva1,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idestr1,2,mpicom_space,ireqa8,mpij)
      end if
      end if
!     **************************************************
      else
!     **********  communicate in a direction  **********
      i=0
      do is=1,ncol
      do ic=1,ncpc
      do ib=1,ncpb
      do iov=1,nsnd
         i=i+1
         vsenda0(i                    )=dreal(v_cc(ncpa-nsnd+iov,ib,ic,is))
         vsenda1(i                    )=dreal(v_cc(          iov,ib,ic,is))
         vsenda0(i+nsnd*ncpb*ncpc*ncol)=dimag(v_cc(ncpa-nsnd+iov,ib,ic,is))
         vsenda1(i+nsnd*ncpb*ncpc*ncol)=dimag(v_cc(          iov,ib,ic,is))
         vrecva0(i                    )=dreal(v_cc(ncpa-nsnd+iov,ib,ic,is))
         vrecva0(i+nsnd*ncpb*ncpc*ncol)=dimag(v_cc(ncpa-nsnd+iov,ib,ic,is))
         vrecva1(i                    )=dreal(v_cc(          iov,ib,ic,is))
         vrecva1(i+nsnd*ncpb*ncpc*ncol)=dimag(v_cc(          iov,ib,ic,is))
      end do
      end do
      end do
      end do

      if (nproca .ne. 1) then
      weightcom=4.0d0*nsnd
      comcount=comcount+2.0d0*nsnd*ncpc*ncpb*weightcom*2*ncol
      vrecva0=0.0d0
      vrecva1=0.0d0
      idests0=myr_space+1
      idests1=myr_space-1
      idestr0=myr_space-1
      idestr1=myr_space+1
      if (myra .eq. nproca-1) idests0=idests0-nproca
      if (myra .eq.        0) idests1=idests1+nproca
      if (myra .eq.        0) idestr0=idestr0+nproca
      if (myra .eq. nproca-1) idestr1=idestr1-nproca
      if (mod(myra,2) .eq. 0) then
      call mpi_isend(vsenda0,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idests0,0,mpicom_space,ireqa1,mpij)
      call mpi_irecv(vrecva1,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idestr1,1,mpicom_space,ireqa2,mpij)
      call mpi_isend(vsenda1,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idests1,2,mpicom_space,ireqa3,mpij)
      call mpi_irecv(vrecva0,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idestr0,3,mpicom_space,ireqa4,mpij)
      end if
      if (mod(myra,2) .eq. 1) then
      call mpi_isend(vsenda1,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idests1,1,mpicom_space,ireqa5,mpij)
      call mpi_irecv(vrecva0,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idestr0,0,mpicom_space,ireqa6,mpij)
      call mpi_isend(vsenda0,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idests0,3,mpicom_space,ireqa7,mpij)
      call mpi_irecv(vrecva1,nsnd*ncpb*ncpc*2*ncol,mpi_double_precision,idestr1,2,mpicom_space,ireqa8,mpij)
      end if
      end if
!     **************************************************
      end if

      if (nperi .le. 1) then
!     **********  communicate in b direction  **********
      if (nprocb .ne. 1) then
      i=0
      do is=1,ncol
      do ic=1,ncpc
      do iov=1,nsnd
      do ia=1,ncpa
         i=i+1
         vsendb0(i                    )=dreal(v_cc(ia,ncpb-nsnd+iov,ic,is))
         vsendb1(i                    )=dreal(v_cc(ia,          iov,ic,is))
         vsendb0(i+ncpa*nsnd*ncpc*ncol)=dimag(v_cc(ia,ncpb-nsnd+iov,ic,is))
         vsendb1(i+ncpa*nsnd*ncpc*ncol)=dimag(v_cc(ia,          iov,ic,is))
      end do
      end do
      end do
      end do
      end if

      if (nprocb .ne. 1) then
      weightcom=2.0d0*nsnd
      if (nprocb .ge. 4) weightcom=4.0d0*nsnd
      comcount=comcount+2.0d0*nsnd*ncpa*ncpb*weightcom*2*ncol
      vrecvb0=0.0d0
      vrecvb1=0.0d0
      idests0=myr_space+nproca
      idests1=myr_space-nproca
      idestr0=myr_space-nproca
      idestr1=myr_space+nproca
      if (myrb .eq. nprocb-1) idests0=idests0-nprocb*nproca
      if (myrb .eq.        0) idests1=idests1+nprocb*nproca
      if (myrb .eq.        0) idestr0=idestr0+nprocb*nproca
      if (myrb .eq. nprocb-1) idestr1=idestr1-nprocb*nproca
      if (mod(myrb,2) .eq. 0) then
      call mpi_isend(vsendb0,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idests0,4,mpicom_space,ireqb1,mpij)
      call mpi_irecv(vrecvb1,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idestr1,5,mpicom_space,ireqb2,mpij)
      end if
      if (mod(myrb,2) .eq. 1) then
      call mpi_isend(vsendb1,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idests1,5,mpicom_space,ireqb3,mpij)
      call mpi_irecv(vrecvb0,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idestr0,4,mpicom_space,ireqb4,mpij)
      end if
      if ((mod(myrb,2) .eq. 0) .and. (myrb .ne. 0)) then
      call mpi_isend(vsendb1,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idests1,6,mpicom_space,ireqb5,mpij)
      call mpi_irecv(vrecvb0,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idestr0,7,mpicom_space,ireqb6,mpij)
      end if
      if ((mod(myrb,2) .eq. 1) .and. (myrb .ne. nprocb-1)) then
      call mpi_isend(vsendb0,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idests0,7,mpicom_space,ireqb7,mpij)
      call mpi_irecv(vrecvb1,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idestr1,6,mpicom_space,ireqb8,mpij)
      end if
      end if
!     **************************************************
      else
!     **********  communicate in b direction  **********
      i=0
      do is=1,ncol
      do ic=1,ncpc
      do iov=1,nsnd
      do ia=1,ncpa
         i=i+1
         vsendb0(i                    )=dreal(v_cc(ia,ncpb-nsnd+iov,ic,is))
         vsendb1(i                    )=dreal(v_cc(ia,          iov,ic,is))
         vsendb0(i+ncpa*nsnd*ncpc*ncol)=dimag(v_cc(ia,ncpb-nsnd+iov,ic,is))
         vsendb1(i+ncpa*nsnd*ncpc*ncol)=dimag(v_cc(ia,          iov,ic,is))
         vrecvb0(i                    )=dreal(v_cc(ia,ncpb-nsnd+iov,ic,is))
         vrecvb1(i                    )=dreal(v_cc(ia,          iov,ic,is))
         vrecvb0(i+ncpa*nsnd*ncpc*ncol)=dimag(v_cc(ia,ncpb-nsnd+iov,ic,is))
         vrecvb1(i+ncpa*nsnd*ncpc*ncol)=dimag(v_cc(ia,          iov,ic,is))
      end do
      end do
      end do
      end do

      if (nprocb .ne. 1) then
      weightcom=4.0d0*nsnd
      comcount=comcount+2.0d0*nsnd*ncpc*ncpa*weightcom*2*ncol
      vrecvb0=0.0d0
      vrecvb1=0.0d0
      idests0=myr_space+nproca
      idests1=myr_space-nproca
      idestr0=myr_space-nproca
      idestr1=myr_space+nproca
      if (myrb .eq. nprocb-1) idests0=idests0-nprocb*nproca
      if (myrb .eq.        0) idests1=idests1+nprocb*nproca
      if (myrb .eq.        0) idestr0=idestr0+nprocb*nproca
      if (myrb .eq. nprocb-1) idestr1=idestr1-nprocb*nproca
      if (mod(myrb,2) .eq. 0) then
      call mpi_isend(vsendb0,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idests0,4,mpicom_space,ireqb1,mpij)
      call mpi_irecv(vrecvb1,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idestr1,5,mpicom_space,ireqb2,mpij)
      call mpi_isend(vsendb1,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idests1,6,mpicom_space,ireqb3,mpij)
      call mpi_irecv(vrecvb0,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idestr0,7,mpicom_space,ireqb4,mpij)
      end if
      if (mod(myrb,2) .eq. 1) then
      call mpi_isend(vsendb1,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idests1,5,mpicom_space,ireqb5,mpij)
      call mpi_irecv(vrecvb0,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idestr0,4,mpicom_space,ireqb6,mpij)
      call mpi_isend(vsendb0,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idests0,7,mpicom_space,ireqb7,mpij)
      call mpi_irecv(vrecvb1,ncpa*nsnd*ncpc*2*ncol,mpi_double_precision,idestr1,6,mpicom_space,ireqb8,mpij)
      end if
      end if
!     **************************************************
      end if

      if (nperi .le. 2) then
!     **********  communicate in c direction  **********
      if (nprocc .ne. 1) then
      i=0
      do is=1,ncol
      do iov=1,nsnd
      do ib=1,ncpb
      do ia=1,ncpa
         i=i+1
         vsendc0(i                    )=dreal(v_cc(ia,ib,ncpc-nsnd+iov,is))
         vsendc1(i                    )=dreal(v_cc(ia,ib,          iov,is))
         vsendc0(i+ncpa*ncpb*nsnd*ncol)=dimag(v_cc(ia,ib,ncpc-nsnd+iov,is))
         vsendc1(i+ncpa*ncpb*nsnd*ncol)=dimag(v_cc(ia,ib,          iov,is))
      end do
      end do
      end do
      end do
      end if

      if (nprocc .ne. 1) then
      weightcom=2.0d0*nsnd
      if (nprocc .ge. 4) weightcom=4.0d0*nsnd
      comcount=comcount+2.0d0*nsnd*ncpb*ncpa*weightcom*2*ncol
      vrecvc0=0.0d0
      vrecvc1=0.0d0
      idests0=myr_space+nproca*nprocb
      idests1=myr_space-nproca*nprocb
      idestr0=myr_space-nproca*nprocb
      idestr1=myr_space+nproca*nprocb
      if (myrc .eq. nprocc-1) idests0=idests0-nprocc*nproca*nprocb
      if (myrc .eq.        0) idests1=idests1+nprocc*nproca*nprocb
      if (myrc .eq.        0) idestr0=idestr0+nprocc*nproca*nprocb
      if (myrc .eq. nprocc-1) idestr1=idestr1-nprocc*nproca*nprocb
      if (mod(myrc,2) .eq. 0) then
      call mpi_isend(vsendc0,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idests0,8,mpicom_space,ireqc1,mpij)
      call mpi_irecv(vrecvc1,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idestr1,9,mpicom_space,ireqc2,mpij)
      end if
      if (mod(myrc,2) .eq. 1) then
      call mpi_isend(vsendc1,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idests1,9,mpicom_space,ireqc3,mpij)
      call mpi_irecv(vrecvc0,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idestr0,8,mpicom_space,ireqc4,mpij)
      end if
      if ((mod(myrc,2) .eq. 0) .and. (myrc .ne. 0)) then
      call mpi_isend(vsendc1,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idests1,10,mpicom_space,ireqc5,mpij)
      call mpi_irecv(vrecvc0,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idestr0,11,mpicom_space,ireqc6,mpij)
      end if
      if ((mod(myrc,2) .eq. 1) .and. (myrc .ne. nprocc-1)) then
      call mpi_isend(vsendc0,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idests0,11,mpicom_space,ireqc7,mpij)
      call mpi_irecv(vrecvc1,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idestr1,10,mpicom_space,ireqc8,mpij)
      end if
      end if
!     **************************************************
      else
!     **********  communicate in c direction  **********
      i=0
      do is=1,ncol
      do iov=1,nsnd
      do ib=1,ncpb
      do ia=1,ncpa
         i=i+1
         vsendc0(i                    )=dreal(v_cc(ia,ib,ncpc-nsnd+iov,is))
         vsendc1(i                    )=dreal(v_cc(ia,ib,          iov,is))
         vsendc0(i+ncpa*ncpb*nsnd*ncol)=dimag(v_cc(ia,ib,ncpc-nsnd+iov,is))
         vsendc1(i+ncpa*ncpb*nsnd*ncol)=dimag(v_cc(ia,ib,          iov,is))
         vrecvc0(i                    )=dreal(v_cc(ia,ib,ncpc-nsnd+iov,is))
         vrecvc1(i                    )=dreal(v_cc(ia,ib,          iov,is))
         vrecvc0(i+ncpa*ncpb*nsnd*ncol)=dimag(v_cc(ia,ib,ncpc-nsnd+iov,is))
         vrecvc1(i+ncpa*ncpb*nsnd*ncol)=dimag(v_cc(ia,ib,          iov,is))
      end do
      end do
      end do
      end do

      if (nprocc .ne. 1) then
      weightcom=4.0d0*nsnd
      comcount=comcount+2.0d0*nsnd*ncpb*ncpa*weightcom*2*ncol
      vrecvc0=0.0d0
      vrecvc1=0.0d0
      idests0=myr_space+nproca*nprocb
      idests1=myr_space-nproca*nprocb
      idestr0=myr_space-nproca*nprocb
      idestr1=myr_space+nproca*nprocb
      if (myrc .eq. nprocc-1) idests0=idests0-nprocc*nproca*nprocb
      if (myrc .eq.        0) idests1=idests1+nprocc*nproca*nprocb
      if (myrc .eq.        0) idestr0=idestr0+nprocc*nproca*nprocb
      if (myrc .eq. nprocc-1) idestr1=idestr1-nprocc*nproca*nprocb
      if (mod(myrc,2) .eq. 0) then
      call mpi_isend(vsendc0,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idests0,8,mpicom_space,ireqc1,mpij)
      call mpi_irecv(vrecvc1,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idestr1,9,mpicom_space,ireqc2,mpij)
      call mpi_isend(vsendc1,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idests1,10,mpicom_space,ireqc3,mpij)
      call mpi_irecv(vrecvc0,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idestr0,11,mpicom_space,ireqc4,mpij)
      end if
      if (mod(myrc,2) .eq. 1) then
      call mpi_isend(vsendc1,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idests1,9,mpicom_space,ireqc5,mpij)
      call mpi_irecv(vrecvc0,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idestr0,8,mpicom_space,ireqc6,mpij)
      call mpi_isend(vsendc0,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idests0,11,mpicom_space,ireqc7,mpij)
      call mpi_irecv(vrecvc1,ncpa*ncpb*nsnd*2*ncol,mpi_double_precision,idestr1,10,mpicom_space,ireqc8,mpij)
      end if
      end if
!     **************************************************
      end if

      return
      end subroutine


      subroutine overlap_fdcheck_r(nperi,ncpa,ncpb,ncpc,nf,nf1,nsnd,v_cc)
      use mod_mpi,myra=>myrx,myrb=>myry,myrc=>myrz,nproca=>nprocx,nprocb=>nprocy,nprocc=>nprocz
      implicit none
      integer nperi,ncpa,ncpb,ncpc,nf,nf1,nsnd
      real*8 v_cc(-nf1:ncpa+nf,-nf1:ncpb+nf,-nf1:ncpc+nf)
      integer ia,ib,ic,i
      integer iov

      if (nperi .eq. 0) then
        if (nproca .ne. 1) then
          if (mod(myra,2) .eq. 0) then
          call mpi_wait(ireqa1,mpistat,mpij)
          call mpi_wait(ireqa2,mpistat,mpij)
          end if
          if (mod(myra,2) .eq. 1) then
          call mpi_wait(ireqa3,mpistat,mpij)
          call mpi_wait(ireqa4,mpistat,mpij)
          end if
          if ((mod(myra,2) .eq. 0) .and. (myra .ne. 0)) then
          call mpi_wait(ireqa5,mpistat,mpij)
          call mpi_wait(ireqa6,mpistat,mpij)
          end if
          if ((mod(myra,2) .eq. 1) .and. (myra .ne. nproca-1)) then
          call mpi_wait(ireqa7,mpistat,mpij)
          call mpi_wait(ireqa8,mpistat,mpij)
          end if
          i=0
          do ic=1,ncpc
          do ib=1,ncpb
          do iov=1,nsnd
             i=i+1
             v_cc(-nsnd+iov,ib,ic)=vrecva0(i)
             v_cc( ncpa+iov,ib,ic)=vrecva1(i)
          end do
          end do
          end do
        end if
      else
        if (nproca .ne. 1) then
          if (mod(myra,2) .eq. 0) then
          call mpi_wait(ireqa1,mpistat,mpij)
          call mpi_wait(ireqa2,mpistat,mpij)
          call mpi_wait(ireqa3,mpistat,mpij)
          call mpi_wait(ireqa4,mpistat,mpij)
          end if
          if (mod(myra,2) .eq. 1) then
          call mpi_wait(ireqa5,mpistat,mpij)
          call mpi_wait(ireqa6,mpistat,mpij)
          call mpi_wait(ireqa7,mpistat,mpij)
          call mpi_wait(ireqa8,mpistat,mpij)
          end if
        end if
        i=0
        do ic=1,ncpc
        do ib=1,ncpb
        do iov=1,nsnd
           i=i+1
           v_cc(-nsnd+iov,ib,ic)=vrecva0(i)
           v_cc( ncpa+iov,ib,ic)=vrecva1(i)
        end do
        end do
        end do
      end if

      if (nperi .le. 1) then
        if (nprocb .ne. 1) then
          if (mod(myrb,2) .eq. 0) then
          call mpi_wait(ireqb1,mpistat,mpij)
          call mpi_wait(ireqb2,mpistat,mpij)
          end if
          if (mod(myrb,2) .eq. 1) then
          call mpi_wait(ireqb3,mpistat,mpij)
          call mpi_wait(ireqb4,mpistat,mpij)
          end if
          if ((mod(myrb,2) .eq. 0) .and. (myrb .ne. 0)) then
          call mpi_wait(ireqb5,mpistat,mpij)
          call mpi_wait(ireqb6,mpistat,mpij)
          end if
          if ((mod(myrb,2) .eq. 1) .and. (myrb .ne. nprocb-1)) then
          call mpi_wait(ireqb7,mpistat,mpij)
          call mpi_wait(ireqb8,mpistat,mpij)
          end if
          i=0
          do ic=1,ncpc
          do iov=1,nsnd
          do ia=1,ncpa
             i=i+1
             v_cc(ia,-nsnd+iov,ic)=vrecvb0(i)
             v_cc(ia, ncpb+iov,ic)=vrecvb1(i)
          end do
          end do
          end do
        end if
      else
        if (nprocb .ne. 1) then
          if (mod(myrb,2) .eq. 0) then
          call mpi_wait(ireqb1,mpistat,mpij)
          call mpi_wait(ireqb2,mpistat,mpij)
          call mpi_wait(ireqb3,mpistat,mpij)
          call mpi_wait(ireqb4,mpistat,mpij)
          end if
          if (mod(myrb,2) .eq. 1) then
          call mpi_wait(ireqb5,mpistat,mpij)
          call mpi_wait(ireqb6,mpistat,mpij)
          call mpi_wait(ireqb7,mpistat,mpij)
          call mpi_wait(ireqb8,mpistat,mpij)
          end if
        end if
        i=0
        do ic=1,ncpc
        do iov=1,nsnd
        do ia=1,ncpa
           i=i+1
           v_cc(ia,-nsnd+iov,ic)=vrecvb0(i)
           v_cc(ia, ncpb+iov,ic)=vrecvb1(i)
        end do
        end do
        end do
      end if

      if (nperi .le. 2) then
        if (nprocc .ne. 1) then
          if (mod(myrc,2) .eq. 0) then
          call mpi_wait(ireqc1,mpistat,mpij)
          call mpi_wait(ireqc2,mpistat,mpij)
          end if
          if (mod(myrc,2) .eq. 1) then
          call mpi_wait(ireqc3,mpistat,mpij)
          call mpi_wait(ireqc4,mpistat,mpij)
          end if
          if ((mod(myrc,2) .eq. 0) .and. (myrc .ne. 0)) then
          call mpi_wait(ireqc5,mpistat,mpij)
          call mpi_wait(ireqc6,mpistat,mpij)
          end if
          if ((mod(myrc,2) .eq. 1) .and. (myrc .ne. nprocc-1)) then
          call mpi_wait(ireqc7,mpistat,mpij)
          call mpi_wait(ireqc8,mpistat,mpij)
          end if
          i=0
          do iov=1,nsnd
          do ib=1,ncpb
          do ia=1,ncpa
             i=i+1
             v_cc(ia,ib,-nsnd+iov)=vrecvc0(i)
             v_cc(ia,ib, ncpc+iov)=vrecvc1(i)
          end do
          end do
          end do
        end if
      else
        if (nprocc .ne. 1) then
          if (mod(myrc,2) .eq. 0) then
          call mpi_wait(ireqc1,mpistat,mpij)
          call mpi_wait(ireqc2,mpistat,mpij)
          call mpi_wait(ireqc3,mpistat,mpij)
          call mpi_wait(ireqc4,mpistat,mpij)
          end if
          if (mod(myrc,2) .eq. 1) then
          call mpi_wait(ireqc5,mpistat,mpij)
          call mpi_wait(ireqc6,mpistat,mpij)
          call mpi_wait(ireqc7,mpistat,mpij)
          call mpi_wait(ireqc8,mpistat,mpij)
          end if
        end if
        i=0
        do iov=1,nsnd
        do ib=1,ncpb
        do ia=1,ncpa
           i=i+1
           v_cc(ia,ib,-nsnd+iov)=vrecvc0(i)
           v_cc(ia,ib, ncpc+iov)=vrecvc1(i)
        end do
        end do
        end do
      end if
      etime=mpi_wtime()
      comtime=comtime+(etime-stime)

      return
      end subroutine


      subroutine overlap_fdcheck_c(nperi,ncpa,ncpb,ncpc,nf,nf1,nsnd,ncol,v_cc)
      use mod_mpi,myra=>myrx,myrb=>myry,myrc=>myrz,nproca=>nprocx,nprocb=>nprocy,nprocc=>nprocz
      implicit none
      integer nperi,ncpa,ncpb,ncpc,nf,nf1,nsnd,ncol
      complex*16 v_cc(-nf1:ncpa+nf,-nf1:ncpb+nf,-nf1:ncpc+nf,ncol)
      integer ia,ib,ic,is,i
      integer iov

      if (nperi .eq. 0) then
        if (nproca .ne. 1) then
          if (mod(myra,2) .eq. 0) then
          call mpi_wait(ireqa1,mpistat,mpij)
          call mpi_wait(ireqa2,mpistat,mpij)
          end if
          if (mod(myra,2) .eq. 1) then
          call mpi_wait(ireqa3,mpistat,mpij)
          call mpi_wait(ireqa4,mpistat,mpij)
          end if
          if ((mod(myra,2) .eq. 0) .and. (myra .ne. 0)) then
          call mpi_wait(ireqa5,mpistat,mpij)
          call mpi_wait(ireqa6,mpistat,mpij)
          end if
          if ((mod(myra,2) .eq. 1) .and. (myra .ne. nproca-1)) then
          call mpi_wait(ireqa7,mpistat,mpij)
          call mpi_wait(ireqa8,mpistat,mpij)
          end if
          i=0
          do is=1,ncol
          do ic=1,ncpc
          do ib=1,ncpb
          do iov=1,nsnd
             i=i+1
             v_cc(-nsnd+iov,ib,ic,is)=dcmplx(vrecva0(i),vrecva0(i+nsnd*ncpb*ncpc*ncol))
             v_cc( ncpa+iov,ib,ic,is)=dcmplx(vrecva1(i),vrecva1(i+nsnd*ncpb*ncpc*ncol))
          end do
          end do
          end do
          end do
        end if
      else
        if (nproca .ne. 1) then
          if (mod(myra,2) .eq. 0) then
          call mpi_wait(ireqa1,mpistat,mpij)
          call mpi_wait(ireqa2,mpistat,mpij)
          call mpi_wait(ireqa3,mpistat,mpij)
          call mpi_wait(ireqa4,mpistat,mpij)
          end if
          if (mod(myra,2) .eq. 1) then
          call mpi_wait(ireqa5,mpistat,mpij)
          call mpi_wait(ireqa6,mpistat,mpij)
          call mpi_wait(ireqa7,mpistat,mpij)
          call mpi_wait(ireqa8,mpistat,mpij)
          end if
        end if
        i=0
        do is=1,ncol
        do ic=1,ncpc
        do ib=1,ncpb
        do iov=1,nsnd
           i=i+1
           v_cc(-nsnd+iov,ib,ic,is)=dcmplx(vrecva0(i),vrecva0(i+nsnd*ncpb*ncpc*ncol))
           v_cc( ncpa+iov,ib,ic,is)=dcmplx(vrecva1(i),vrecva1(i+nsnd*ncpb*ncpc*ncol))
        end do
        end do
        end do
        end do
      end if

      if (nperi .le. 1) then
        if (nprocb .ne. 1) then
          if (mod(myrb,2) .eq. 0) then
          call mpi_wait(ireqb1,mpistat,mpij)
          call mpi_wait(ireqb2,mpistat,mpij)
          end if
          if (mod(myrb,2) .eq. 1) then
          call mpi_wait(ireqb3,mpistat,mpij)
          call mpi_wait(ireqb4,mpistat,mpij)
          end if
          if ((mod(myrb,2) .eq. 0) .and. (myrb .ne. 0)) then
          call mpi_wait(ireqb5,mpistat,mpij)
          call mpi_wait(ireqb6,mpistat,mpij)
          end if
          if ((mod(myrb,2) .eq. 1) .and. (myrb .ne. nprocb-1)) then
          call mpi_wait(ireqb7,mpistat,mpij)
          call mpi_wait(ireqb8,mpistat,mpij)
          end if
          i=0
          do is=1,ncol
          do ic=1,ncpc
          do iov=1,nsnd
          do ia=1,ncpa
             i=i+1
             v_cc(ia,-nsnd+iov,ic,is)=dcmplx(vrecvb0(i),vrecvb0(i+ncpa*nsnd*ncpc*ncol))
             v_cc(ia, ncpb+iov,ic,is)=dcmplx(vrecvb1(i),vrecvb1(i+ncpa*nsnd*ncpc*ncol))
          end do
          end do
          end do
          end do
        end if
      else
        if (nprocb .ne. 1) then
          if (mod(myrb,2) .eq. 0) then
          call mpi_wait(ireqb1,mpistat,mpij)
          call mpi_wait(ireqb2,mpistat,mpij)
          call mpi_wait(ireqb3,mpistat,mpij)
          call mpi_wait(ireqb4,mpistat,mpij)
          end if
          if (mod(myrb,2) .eq. 1) then
          call mpi_wait(ireqb5,mpistat,mpij)
          call mpi_wait(ireqb6,mpistat,mpij)
          call mpi_wait(ireqb7,mpistat,mpij)
          call mpi_wait(ireqb8,mpistat,mpij)
          end if
        end if
        i=0
        do is=1,ncol
        do ic=1,ncpc
        do iov=1,nsnd
        do ia=1,ncpa
           i=i+1
           v_cc(ia,-nsnd+iov,ic,is)=dcmplx(vrecvb0(i),vrecvb0(i+ncpa*nsnd*ncpc*ncol))
           v_cc(ia, ncpb+iov,ic,is)=dcmplx(vrecvb1(i),vrecvb1(i+ncpa*nsnd*ncpc*ncol))
        end do
        end do
        end do
        end do
      end if

      if (nperi .le. 2) then
        if (nprocc .ne. 1) then
          if (mod(myrc,2) .eq. 0) then
          call mpi_wait(ireqc1,mpistat,mpij)
          call mpi_wait(ireqc2,mpistat,mpij)
          end if
          if (mod(myrc,2) .eq. 1) then
          call mpi_wait(ireqc3,mpistat,mpij)
          call mpi_wait(ireqc4,mpistat,mpij)
          end if
          if ((mod(myrc,2) .eq. 0) .and. (myrc .ne. 0)) then
          call mpi_wait(ireqc5,mpistat,mpij)
          call mpi_wait(ireqc6,mpistat,mpij)
          end if
          if ((mod(myrc,2) .eq. 1) .and. (myrc .ne. nprocc-1)) then
          call mpi_wait(ireqc7,mpistat,mpij)
          call mpi_wait(ireqc8,mpistat,mpij)
          end if
          i=0
          do is=1,ncol
          do iov=1,nsnd
          do ib=1,ncpb
          do ia=1,ncpa
             i=i+1
             v_cc(ia,ib,-nsnd+iov,is)=dcmplx(vrecvc0(i),vrecvc0(i+ncpa*ncpb*nsnd*ncol))
             v_cc(ia,ib, ncpc+iov,is)=dcmplx(vrecvc1(i),vrecvc1(i+ncpa*ncpb*nsnd*ncol))
          end do
          end do
          end do
          end do
        end if
      else
        if (nprocc .ne. 1) then
          if (mod(myrc,2) .eq. 0) then
          call mpi_wait(ireqc1,mpistat,mpij)
          call mpi_wait(ireqc2,mpistat,mpij)
          call mpi_wait(ireqc3,mpistat,mpij)
          call mpi_wait(ireqc4,mpistat,mpij)
          end if
          if (mod(myrc,2) .eq. 1) then
          call mpi_wait(ireqc5,mpistat,mpij)
          call mpi_wait(ireqc6,mpistat,mpij)
          call mpi_wait(ireqc7,mpistat,mpij)
          call mpi_wait(ireqc8,mpistat,mpij)
          end if
        end if
        i=0
        do is=1,ncol
        do iov=1,nsnd
        do ib=1,ncpb
        do ia=1,ncpa
           i=i+1
           v_cc(ia,ib,-nsnd+iov,is)=dcmplx(vrecvc0(i),vrecvc0(i+ncpa*ncpb*nsnd*ncol))
           v_cc(ia,ib, ncpc+iov,is)=dcmplx(vrecvc1(i),vrecvc1(i+ncpa*ncpb*nsnd*ncol))
        end do
        end do
        end do
        end do
      end if
      etime=mpi_wtime()
      comtime=comtime+(etime-stime)

      return
      end subroutine


      end module
