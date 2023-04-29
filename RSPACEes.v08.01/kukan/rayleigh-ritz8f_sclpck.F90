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
! **********  rayleighritz8f_sclpck.F90 04/26/2023-01  **********

module mod_rayleighritz
use mod_stopp
implicit none
integer ndesca(9)
integer ncontext,nprow,npcol,myrow,mycol,mxllda,mxlocc,mb,lwork,lrwork,liwork
integer,allocatable::lrow(:),lcol(:),jrow(:),jcol(:)

contains


subroutine rayleighritz_blacs_ini(neigmx,nrc)
use mod_mpi
implicit none
integer, intent(in)::neigmx,nrc
complex*16,allocatable::ca(:,:),cz(:,:),cwork(:)
real*8,allocatable::ra(:,:),rz(:,:),eig(:),rwork(:)
integer,allocatable::iwork(:)
integer i,j,k,ierr
integer numroc

  allocate(lrow(neigmx),lcol(neigmx),jrow(neigmx),jcol(neigmx))
  nprow=int(dsqrt(dfloat(nprocs)))
  npcol=nprocs/nprow
  mb=4
  mxllda=(neigmx-1)/nprow+1+mod((neigmx-1)/nprow+1,mb)
  mxlocc=(neigmx-1)/npcol+1+mod((neigmx-1)/npcol+1,mb)
  mxllda=((mxllda-1)/mb+1)*mb
  mxlocc=((mxlocc-1)/mb+1)*mb
!  call sl_init(ncontext,nprow,npcol)
  ncontext=mpicom_space
  call blacs_gridinit(ncontext,'Col',nprow,npcol)
  call blacs_gridinfo(ncontext,nprow,npcol,myrow,mycol)
  if(myrow /= -1) then
    call descinit(ndesca,neigmx,neigmx,mb,mb,0,0,ncontext,mxllda,ierr)
    if (nrc == 0) then
      lwork=max(1+6*neigmx+2*numroc(neigmx,mb,myrow,0,nprow)*numroc(neigmx,mb,mycol,0,npcol), &
         3*neigmx+max(mb*(numroc(neigmx,mb,myrow,0,nprow)+1),3*mb))+2*neigmx
      liwork=7*neigmx+8*npcol+2
      allocate(ra(mxllda,mxlocc),rz(mxllda,mxlocc),eig(neigmx),rwork(lwork),iwork(liwork))
!     lwork is multiplied by 10 and lower triangle is refferd to avoid problem in ohtaka 04/13/2023.
!      call pdsyevd('v','u',neigmx,ra,1,1,ndesca,eig,rz,1,1,ndesca,rwork,-1,iwork,liwork,ierr)
      call pdsyevd('v','l',neigmx,ra,1,1,ndesca,eig,rz,1,1,ndesca,rwork,-1,iwork,liwork,ierr)
      lwork=nint(rwork(1))
      deallocate(ra,rz,eig,rwork,iwork)
    else
      lwork=neigmx+(numroc(max(neigmx,mb,2),mb,0,0,nprow)+numroc(max(neigmx,mb,2),mb,0,0,npcol)+mb)*mb
      lrwork=1+9*neigmx+3*numroc(neigmx,mb,myrow,0,nprow)*numroc(neigmx,mb,mycol,0,npcol)
      liwork=7*neigmx+8*npcol+2
      allocate(ca(mxllda,mxlocc),cz(mxllda,mxlocc),eig(neigmx),cwork(lwork),rwork(lrwork),iwork(liwork))
      ca=dcmplx(0.0d0,0.0d0)
      do k=1,neigmx
        call pzelset(ca,k,k,ndesca,dcmplx(1.0d0,0.0d0))
      end do
      call pzheevd('v','u',neigmx,ca,1,1,ndesca,eig,cz,1,1,ndesca,cwork,-1,rwork,lrwork,iwork,liwork,ierr)
      lwork=nint(dreal(cwork(1)))
      lrwork=nint(rwork(1))
      deallocate(ca,cz,eig,cwork,rwork,iwork)
    end if
  end if
  i=0
  j=0
  do k=1,neigmx
    i=i+1
    if (i>mb) then
      i=1
      j=j+1
    end if
    if (j>=nprow) j=0
    lrow(k)=j
  end do
  i=0
  j=0
  do k=1,neigmx
    i=i+1
    if (i>mb) then
      i=1
      j=j+1
    end if
    if (j>=npcol) j=0
    lcol(k)=j
  end do
  do k=0,nprow-1
    i=0
    do j=1,neigmx
      if (k==lrow(j)) then
        i=i+1
        jrow(j)=i
      end if
    end do
  end do
  do k=0,npcol-1
    i=0
    do j=1,neigmx
      if (k==lcol(j)) then
        i=i+1
        jcol(j)=i
      end if
    end do
  end do
  call mpi_bcast(lrow,neigmx,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(lcol,neigmx,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(jrow,neigmx,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(jcol,neigmx,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(nprow,1,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(npcol,1,mpi_integer,0,mpicom_space,mpij)
return
end subroutine


subroutine rayleighritz_blacs_fin
implicit none
  deallocate(lrow,lcol,jrow,jcol)
  if(myrow /= -1) then
    call blacs_gridexit(ncontext)
  end if
  call blacs_exit(1)
return
end subroutine


subroutine rayleighritz_r( &
 natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,nspv, & ! <
 num_list,ncpx,ncpy,ncpz,                               & ! <
 nperi,nk,ns,nf,ndisp,                                  & ! <
 key_natpri_in,key_natpri_inps,                         & ! <
 dx,dy,dz,                                              & ! <
 nprj,natinf,lstvec2,indspe,natpri,naps,ntyppp,         & ! <
 dij,veff,vnlocp,                                       & ! <
 l_rrb,                                                 & ! <
 rspsep,svecre,ssvre,                                   & ! X
 hsvre,                                                 & ! K
 sval)                                                    ! >
use mod_mpi
use mod_stopp
use mod_disptime
use mod_overlap_finitedifference, only:overlap_finitedifference_init,overlap_finitedifference_final, &
                                       overlap_finitedifference_r,overlap_fdcheck_r
implicit none
integer, intent(in)   ::natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,nspv
integer, intent(in)   ::num_list,ncpx,ncpy,ncpz
integer, intent(in)   ::nperi,nk,ns,nf,ndisp
integer, intent(in)   ::key_natpri_in,key_natpri_inps
real*8,  intent(in)   ::dx,dy,dz
integer, intent(in)   ::nprj(num_spe),natinf(natom),lstvec2(num_list,num_ppcell)
integer, intent(in)   ::indspe(natom),natpri(natom),naps(natom),ntyppp(num_spe)
real*8,  intent(in)   ::dij(nprjmx,nprjmx,nums,natom)
real*8,  intent(in)   ::veff(ncpx,ncpy,ncpz,nspv)
real*8,  intent(in)   ::vnlocp(num_list,nprjmx,num_ppcell)
real*8,  intent(inout)::svecre(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8,  intent(inout)::rspsep(nprjmx,num_ppcell,neigmx,nums,numk)
real*8,  intent(inout)::ssvre(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8,  intent(inout)::hsvre(ncpx,ncpy,ncpz,neigmx)
real*8,  intent(out)  ::sval(neigmx,nums,numk)
logical, intent(in)   ::l_rrb
real*8    ,allocatable::ra(:,:),rz(:,:),rzall(:,:),rwork(:)           ! for isolated system
real*8    ,allocatable::eig(:)
real*8    ,allocatable::rspsep_rr(:,:,:) ! for isolated system
integer   ,allocatable::iwork(:)
real*8    ,allocatable::sssvre(:,:,:,:)
integer   ::ierr,i,j,k,jj,kk
real*8 ::stime0
integer   ::mstep

  mstep=32
  if (neigmx .gt. 1024) mstep=64
  if (neigmx .gt. 2048) mstep=128
  if (neigmx .gt. 4096) mstep=256

  if (myrank_glbl==0) write(ndisp,*) 'diagonalization in the subspace was implemented.'

  call starttime(stime0)
  allocate(ra(mxllda,mxlocc),eig(neigmx))
  allocate(sssvre(ncpx,ncpy,ncpz,neigmx),rspsep_rr(nprjmx,num_ppcell,neigmx))
  allocate(rz(mstep,mstep),rzall(mstep,mstep))
  if (l_rrb) then
    call rayleighritz_hamiltonian_r( &
     natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,nspv, & ! <
     num_list,ncpx,ncpy,ncpz,                               & ! <
     nperi,nk,ns,nf,                                        & ! <
     key_natpri_in,key_natpri_inps,                         & ! <
     dx,dy,dz,                                              & ! <
     nprj,natinf,lstvec2,indspe,natpri,naps,                & ! <
     dij,veff,vnlocp,svecre,rspsep,                         & ! <
     hsvre)                                                   ! >
  end if ! l_rrb

! ==========  setup matrix  ==========
  k=0
  do while (k<neigmx)
    if (neigmx-k>=mstep) then
      j=0
      do while (j<neigmx)
        if (neigmx-j>=mstep) then
          call dgemm('t','n',mstep,mstep,ncpx*ncpy*ncpz,1.0d0,svecre(1,1,1,j+1,ns,nk) &
                    ,ncpx*ncpy*ncpz,hsvre(1,1,1,k+1),ncpx*ncpy*ncpz,0.0d0,rz,mstep)
          call mpi_allreduce(rz,rzall,mstep*mstep,mpi_double_precision,mpi_sum,mpicom_space,mpij)
          do kk=1,mstep,mb
            do jj=1,mstep,mb
              if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) &
                  ra(jrow(j+jj):jrow(j+jj)+mb-1,jcol(k+kk):jcol(k+kk)+mb-1)=rzall(jj:jj+mb-1,kk:kk+mb-1)*dx*dy*dz
            end do
          end do
          j=j+mstep
        else
          call dgemm('t','n',neigmx-j,mstep,ncpx*ncpy*ncpz,1.0d0,svecre(1,1,1,j+1,ns,nk) &
                    ,ncpx*ncpy*ncpz,hsvre(1,1,1,k+1),ncpx*ncpy*ncpz,0.0d0,rz,mstep)
          do kk=1,mstep
            call mpi_allreduce(rz(1,kk),rzall(1,kk),neigmx-j,mpi_double_precision,mpi_sum,mpicom_space,mpij)
          end do
          do kk=1,mstep,mb
            do jj=1,neigmx-j
              if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) &
                ra(jrow(j+jj),jcol(k+kk):jcol(k+kk)+mb-1)=rzall(jj,kk:kk+mb-1)*dx*dy*dz
            end do
          end do
          j=neigmx
        end if
      end do
      k=k+mstep
    else
      j=0
      do while (j<neigmx)
        if (neigmx-j>=mstep) then
          call dgemm('t','n',mstep,neigmx-k,ncpx*ncpy*ncpz,1.0d0,svecre(1,1,1,j+1,ns,nk) &
                    ,ncpx*ncpy*ncpz,hsvre(1,1,1,k+1),ncpx*ncpy*ncpz,0.0d0,rz,mstep)
          call mpi_allreduce(rz(1,1),rzall(1,1),mstep*(neigmx-k),mpi_double_precision,mpi_sum,mpicom_space,mpij)
          do kk=1,neigmx-k
            do jj=1,mstep,mb
              if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) &
                ra(jrow(j+jj):jrow(j+jj)+mb-1,jcol(k+kk))=rzall(jj:jj+mb-1,kk)*dx*dy*dz
            end do
          end do
          j=j+mstep
        else
          call dgemm('t','n',neigmx-j,neigmx-k,ncpx*ncpy*ncpz,1.0d0,svecre(1,1,1,j+1,ns,nk) &
                    ,ncpx*ncpy*ncpz,hsvre(1,1,1,k+1),ncpx*ncpy*ncpz,0.0d0,rz,mstep)
          do kk=1,neigmx-k
            call mpi_allreduce(rz(1,kk),rzall(1,kk),neigmx-j,mpi_double_precision,mpi_sum,mpicom_space,mpij)
          end do
          do kk=1,neigmx-k
            do jj=1,neigmx-j
              if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) ra(jrow(j+jj),jcol(k+kk))=rzall(jj,kk)*dx*dy*dz
            end do
          end do
          j=neigmx
        end if
      end do
      k=neigmx
    end if
  end do
  deallocate(rz,rzall)
  allocate(rz(mxllda,mxlocc))
! ====================================
  call endtime(ndisp,stime0,'[TI_sub] setup matrix :')

  call starttime(stime0)
! ==========  eigenproblem  ==========
  ierr=0
  if (myrow /= -1) then
    allocate(rwork(lwork),iwork(liwork))
!    call pdsyevd('v','u',neigmx,ra,1,1,ndesca,eig,rz,1,1,ndesca,rwork,lwork,iwork,liwork,ierr)
    call pdsyevd('v','l',neigmx,ra,1,1,ndesca,eig,rz,1,1,ndesca,rwork,lwork,iwork,liwork,ierr)
    deallocate(rwork,iwork)
    sval(:,ns,nk)=eig(:)
  end if
  call mpi_bcast(ierr,1,mpi_integer,0,mpicom_space,mpij)
  if (ierr /= 0) call stopp('rayleighritz: error occurs in pdsyevd')
  call mpi_bcast(sval(1,ns,nk),neigmx,mpi_double_precision,0,mpicom_space,mpij)
! ====================================
  call endtime(ndisp,stime0,'[TI_sub] sub-space diag. :')

  call starttime(stime0)
  deallocate(ra)
  allocate(ra(mstep,mstep),rzall(mstep,mstep))
! ==========  update vectors  ==========
!$omp parallel default(shared) private(i)
!$omp do
  do i=1,ncpx*ncpy*ncpz*neigmx
    hsvre(i,1,1,1)=0.0d0
    sssvre(i,1,1,1)=0.0d0
  end do
!$omp do
  do i=1,nprjmx*num_ppcell*neigmx
    rspsep_rr(i,1,1)=0.0d0
  end do
!$omp end parallel
  k=0
  do while (k<neigmx)
    if (neigmx-k>=mstep) then
      j=0
      do while (j<neigmx)
!$omp parallel default(shared) private(i)
!$omp do
        do i=1,mstep*mstep
          rzall(i,1)=0.0d0
        end do
!$omp end parallel
        if (neigmx-j>=mstep) then
          do kk=1,mstep,mb
            do jj=1,mstep,mb
              if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) &
                  rzall(jj:jj+mb-1,kk:kk+mb-1)=rz(jrow(j+jj):jrow(j+jj)+mb-1,jcol(k+kk):jcol(k+kk)+mb-1)
            end do
          end do
          call mpi_allreduce(rzall,ra,mstep*mstep,mpi_double_precision,mpi_sum,mpicom_space,mpij)
          call dgemm('n','n',ncpx*ncpy*ncpz,mstep,mstep,1.0d0,svecre(1,1,1,j+1,ns,nk) &
                    ,ncpx*ncpy*ncpz,ra,mstep,1.0d0,hsvre(1,1,1,k+1),ncpx*ncpy*ncpz)
          call dgemm('n','n',ncpx*ncpy*ncpz,mstep,mstep,1.0d0,ssvre(1,1,1,j+1,ns,nk) &
                    ,ncpx*ncpy*ncpz,ra,mstep,1.0d0,sssvre(1,1,1,k+1),ncpx*ncpy*ncpz)
          call dgemm('n','n',nprjmx*num_ppcell,mstep,mstep,1.0d0,rspsep(1,1,j+1,ns,nk) &
                    ,nprjmx*num_ppcell,ra,mstep,1.0d0,rspsep_rr(1,1,k+1),nprjmx*num_ppcell)
          j=j+mstep
        else
          do kk=1,mstep,mb
            do jj=1,neigmx-j
              if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) &
                rzall(jj,kk:kk+mb-1)=rz(jrow(j+jj),jcol(k+kk):jcol(k+kk)+mb-1)
            end do
          end do
          do kk=1,mstep
            call mpi_allreduce(rzall(1,kk),ra(1,kk),neigmx-j,mpi_double_precision,mpi_sum,mpicom_space,mpij)
          end do
          call dgemm('n','n',ncpx*ncpy*ncpz,mstep,neigmx-j,1.0d0,svecre(1,1,1,j+1,ns,nk) &
                    ,ncpx*ncpy*ncpz,ra,mstep,1.0d0,hsvre(1,1,1,k+1),ncpx*ncpy*ncpz)
          call dgemm('n','n',ncpx*ncpy*ncpz,mstep,neigmx-j,1.0d0,ssvre(1,1,1,j+1,ns,nk) &
                    ,ncpx*ncpy*ncpz,ra,mstep,1.0d0,sssvre(1,1,1,k+1),ncpx*ncpy*ncpz)
          call dgemm('n','n',nprjmx*num_ppcell,mstep,neigmx-j,1.0d0,rspsep(1,1,j+1,ns,nk) &
                    ,nprjmx*num_ppcell,ra,mstep,1.0d0,rspsep_rr(1,1,k+1),nprjmx*num_ppcell)
          j=neigmx
        end if
      end do
      k=k+mstep
    else
      j=0
      do while (j<neigmx)
!$omp parallel default(shared) private(i)
!$omp do
        do i=1,mstep*mstep
          rzall(i,1)=0.0d0
        end do
!$omp end parallel
        if (neigmx-j>=mstep) then
          do kk=1,neigmx-k
            do jj=1,mstep,mb
              if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) &
                rzall(jj:jj+mb-1,kk)=rz(jrow(j+jj):jrow(j+jj)+mb-1,jcol(k+kk))
            end do
          end do
          call mpi_allreduce(rzall(1,1),ra(1,1),mstep*(neigmx-k),mpi_double_precision,mpi_sum,mpicom_space,mpij)
          call dgemm('n','n',ncpx*ncpy*ncpz,neigmx-k,mstep,1.0d0,svecre(1,1,1,j+1,ns,nk) &
                    ,ncpx*ncpy*ncpz,ra,mstep,1.0d0,hsvre(1,1,1,k+1),ncpx*ncpy*ncpz)
          call dgemm('n','n',ncpx*ncpy*ncpz,neigmx-k,mstep,1.0d0,ssvre(1,1,1,j+1,ns,nk) &
                    ,ncpx*ncpy*ncpz,ra,mstep,1.0d0,sssvre(1,1,1,k+1),ncpx*ncpy*ncpz)
          call dgemm('n','n',nprjmx*num_ppcell,neigmx-k,mstep,1.0d0,rspsep(1,1,j+1,ns,nk) &
                    ,nprjmx*num_ppcell,ra,mstep,1.0d0,rspsep_rr(1,1,k+1),nprjmx*num_ppcell)
          j=j+mstep
        else
          do kk=1,neigmx-k
            do jj=1,neigmx-j
              if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) rzall(jj,kk)=rz(jrow(j+jj),jcol(k+kk))
            end do
          end do
           do kk=1,neigmx-k
            call mpi_allreduce(rzall(1,kk),ra(1,kk),neigmx-j,mpi_double_precision,mpi_sum,mpicom_space,mpij)
          end do
          call dgemm('n','n',ncpx*ncpy*ncpz,neigmx-k,neigmx-j,1.0d0,svecre(1,1,1,j+1,ns,nk) &
                    ,ncpx*ncpy*ncpz,ra,mstep,1.0d0,hsvre(1,1,1,k+1),ncpx*ncpy*ncpz)
          call dgemm('n','n',ncpx*ncpy*ncpz,neigmx-k,neigmx-j,1.0d0,ssvre(1,1,1,j+1,ns,nk) &
                    ,ncpx*ncpy*ncpz,ra,mstep,1.0d0,sssvre(1,1,1,k+1),ncpx*ncpy*ncpz)
          call dgemm('n','n',nprjmx*num_ppcell,neigmx-k,neigmx-j,1.0d0,rspsep(1,1,j+1,ns,nk) &
                    ,nprjmx*num_ppcell,ra,mstep,1.0d0,rspsep_rr(1,1,k+1),nprjmx*num_ppcell)
          j=neigmx
        end if
      end do
      k=neigmx
    end if
  end do
!$omp parallel default(shared) private(i)
!$omp do
  do i=1,ncpx*ncpy*ncpz*neigmx
    svecre(i,1,1,1,ns,nk)=hsvre(i,1,1,1)
    ssvre(i,1,1,1,ns,nk)=sssvre(i,1,1,1)
  end do
!$omp do
  do i=1,nprjmx*num_ppcell*neigmx
    rspsep(i,1,1,ns,nk)=rspsep_rr(i,1,1)
  end do
!$omp end parallel
  deallocate(ra,rzall)
  deallocate(sssvre,rspsep_rr)
! ======================================
  call endtime(ndisp,stime0,'[TI_sub] update vectors :')

  deallocate(rz,eig)

end subroutine rayleighritz_r


subroutine rayleighritz_c( &
 nso,natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,ncol,nspv, & ! <
 num_list,ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                   & ! <
 nperi,nk,ns1,nf,ndisp,                                          & ! <
 key_natpri_in,key_natpri_inps,key_soc_calc,                     & ! <
 dx,dy,dz,xmax,ymax,zmax,skpxx,skpyy,skpzz,                      & ! <
 nprj,natinf,lstvec2,indspe,natpri,naps,natsoc,ntyppp,           & ! <
 lstx,lsty,lstz,natx,naty,natz,                                  & ! <
 dij,dijsoc,veff,vnlocp,                                         & ! <
 l_rrb,                                                          & ! <
 cspsep,sveccm,ssvcm,                                            & ! X
 hsvcm,                                                          & ! K
 sval)                                                             ! >
use mod_mpi
use mod_stopp
use mod_disptime
use mod_overlap_finitedifference, only:overlap_finitedifference_init,overlap_finitedifference_final, &
                                       overlap_finitedifference_c,overlap_fdcheck_c
implicit none
integer,    intent(in)   ::nso,natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,ncol,nspv
integer,    intent(in)   ::num_list,ncpx,ncpy,ncpz,npxmax,npymax,npzmax
integer,    intent(in)   ::nperi,nk,ns1,nf,ndisp
integer,    intent(in)   ::key_natpri_in,key_natpri_inps,key_soc_calc
real*8,     intent(in)   ::dx,dy,dz,xmax,ymax,zmax,skpxx(numk),skpyy(numk),skpzz(numk)
integer,    intent(in)   ::nprj(num_spe),natinf(natom),lstvec2(num_list,num_ppcell)
integer,    intent(in)   ::indspe(natom),natpri(natom),naps(natom),natsoc(natom),ntyppp(num_spe)
integer,    intent(in)   ::lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell)
integer,    intent(in)   ::natx(natom),naty(natom),natz(natom)
real*8,     intent(in)   ::dij(nprjmx,nprjmx,nums*ncol,natom)
real*8,     intent(in)   ::dijsoc(nprjmx*nso-nso+1,nprjmx*nso-nso+1,3*nso-nso+1,natom*nso-nso+1)
real*8,     intent(in)   ::veff(ncpx,ncpy,ncpz,nspv)
real*8,     intent(in)   ::vnlocp(num_list,nprjmx,num_ppcell)
complex*16, intent(inout)::sveccm(ncpx,ncpy,ncpz,neigmx,nums,numk)
complex*16, intent(inout)::cspsep(nprjmx,num_ppcell,neigmx,nums,numk)
complex*16, intent(inout)::hsvcm(ncpx,ncpy,ncpz,neigmx,ncol)
complex*16, intent(inout)::ssvcm(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8,     intent(out)  ::sval(neigmx,nums+1-ncol,numk)
logical,intent(in) ::l_rrb
complex*16,allocatable::ca(:,:),cz(:,:),czall(:,:),cwork(:) ! for periodic system
real*8    ,allocatable::eig(:)
complex*16,allocatable::cspsep_rr(:,:,:,:) ! for periodic system
integer   ,allocatable::iwork(:)
complex*16,allocatable::sssvcm(:,:,:,:,:)
real*8,    allocatable::rwork(:)
integer   ::ierr,i,j,k,jj,kk,ns,ns2
real*8    ::stime0
complex*16::ctmp0,ctmp1
integer   ::mstep

  mstep=64
  if (neigmx .gt. 1024) mstep=64
  if (neigmx .gt. 2048) mstep=128
  if (neigmx .gt. 4096) mstep=256

  if (myrank_glbl==0) write(ndisp,*) 'diagonalization in the subspace was implemented.'

  call starttime(stime0)
! ==========  setup matrix  ==========
  allocate(ca(mxllda,mxlocc),eig(neigmx))
  allocate(sssvcm(ncpx,ncpy,ncpz,neigmx,ncol),cspsep_rr(nprjmx,num_ppcell,neigmx,ncol))
  allocate(cz(mstep,mstep),czall(mstep,mstep))
  if (l_rrb) then
    call rayleighritz_hamiltonian_c( &
     nso,natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,ncol,nspv, & ! <
     num_list,ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                   & ! <
     nperi,nk,ns1,nf,                                                & ! <
     key_natpri_in,key_natpri_inps,key_soc_calc,                     & ! <
     dx,dy,dz,xmax,ymax,zmax,skpxx,skpyy,skpzz,                      & ! <
     nprj,natinf,lstvec2,indspe,natpri,naps,natsoc,                  & ! <
     lstx,lsty,lstz,natx,naty,natz,                                  & ! <
     dij,dijsoc,veff,vnlocp,sveccm,cspsep,                           & ! <
     hsvcm)                                                            ! >
  end if ! l_rrb

  k=0
  do while (k<neigmx)
    if (neigmx-k>=mstep) then
      j=0
      do while (j<neigmx)
        if (neigmx-j>=mstep) then
          cz(:,:)=dcmplx(0.0d0,0.0d0)
          do ns= 1,ncol
            ns2= ns1+ns-1
            ctmp0=dcmplx(dx*dy*dz,0.0d0)
            ctmp1=dcmplx(1.0d0,0.0d0)
            call zgemm('c','n',mstep,mstep,ncpx*ncpy*ncpz,ctmp0,sveccm(1,1,1,j+1,ns2,nk) &
                      ,ncpx*ncpy*ncpz,hsvcm(1,1,1,k+1,ns),ncpx*ncpy*ncpz,ctmp1,cz,mstep)
          end do
          call mpi_allreduce(cz,czall,mstep*mstep,mpi_double_complex,mpi_sum,mpicom_space,mpij)
          do kk=1,mstep,mb
            do jj=1,mstep,mb
              if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) &
                  ca(jrow(j+jj):jrow(j+jj)+mb-1,jcol(k+kk):jcol(k+kk)+mb-1)=czall(jj:jj+mb-1,kk:kk+mb-1)
            end do
          end do
          j=j+mstep
        else
          cz(:,:)=dcmplx(0.0d0,0.0d0)
          do ns= 1,ncol
            ns2= ns1+ns-1
            ctmp0=dcmplx(dx*dy*dz,0.0d0)
            ctmp1=dcmplx(1.0d0,0.0d0)
            call zgemm('c','n',neigmx-j,mstep,ncpx*ncpy*ncpz,ctmp0,sveccm(1,1,1,j+1,ns2,nk) &
                      ,ncpx*ncpy*ncpz,hsvcm(1,1,1,k+1,ns),ncpx*ncpy*ncpz,ctmp1,cz,mstep)
          end do
          do kk=1,mstep
            call mpi_allreduce(cz(1,kk),czall(1,kk),neigmx-j,mpi_double_complex,mpi_sum,mpicom_space,mpij)
          end do
          do kk=1,mstep,mb
            do jj=1,neigmx-j
              if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) &
                ca(jrow(j+jj),jcol(k+kk):jcol(k+kk)+mb-1)=czall(jj,kk:kk+mb-1)
            end do
          end do
          j=neigmx
        end if
      end do
      k=k+mstep
    else
      j=0
      do while (j<neigmx)
        if (neigmx-j>=mstep) then
          cz(:,:)=dcmplx(0.0d0,0.0d0)
          do ns= 1,ncol
            ns2= ns1+ns-1
            ctmp0=dcmplx(dx*dy*dz,0.0d0)
            ctmp1=dcmplx(1.0d0,0.0d0)
            call zgemm('c','n',mstep,neigmx-k,ncpx*ncpy*ncpz,ctmp0,sveccm(1,1,1,j+1,ns2,nk) &
                      ,ncpx*ncpy*ncpz,hsvcm(1,1,1,k+1,ns),ncpx*ncpy*ncpz,ctmp1,cz,mstep)
          end do
          call mpi_allreduce(cz(1,1),czall(1,1),mstep*(neigmx-k),mpi_double_complex,mpi_sum,mpicom_space,mpij)
          do kk=1,neigmx-k
            do jj=1,mstep,mb
              if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) &
                ca(jrow(j+jj):jrow(j+jj)+mb-1,jcol(k+kk))=czall(jj:jj+mb-1,kk)
            end do
          end do
          j=j+mstep
        else
          cz(:,:)=dcmplx(0.0d0,0.0d0)
          do ns= 1,ncol
            ns2= ns1+ns-1
            ctmp0=dcmplx(dx*dy*dz,0.0d0)
            ctmp1=dcmplx(1.0d0,0.0d0)
            call zgemm('c','n',neigmx-j,neigmx-k,ncpx*ncpy*ncpz,ctmp0,sveccm(1,1,1,j+1,ns2,nk) &
                      ,ncpx*ncpy*ncpz,hsvcm(1,1,1,k+1,ns),ncpx*ncpy*ncpz,ctmp1,cz,mstep)
          end do
          do kk=1,neigmx-k
            call mpi_allreduce(cz(1,kk),czall(1,kk),neigmx-j,mpi_double_complex,mpi_sum,mpicom_space,mpij)
          end do
          do kk=1,neigmx-k
            do jj=1,neigmx-j
              if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) ca(jrow(j+jj),jcol(k+kk))=czall(jj,kk)
            end do
          end do
          j=neigmx
        end if
      end do
      k=neigmx
    end if
  end do
  deallocate(cz,czall)
  allocate(cz(mxllda,mxlocc))
! ====================================
  call endtime(ndisp,stime0,'[TI_sub] setup matrix :')

  call starttime(stime0)
! ==========  eigenproblem  ==========
  ierr=0
  if(myrow /= -1) then
    allocate(cwork(lwork),rwork(lrwork),iwork(liwork))
!    call pzheevd('v','u',neigmx,ca,1,1,ndesca,eig,cz,1,1,ndesca,cwork,lwork,rwork,lrwork,iwork,liwork,ierr)
    call pzheevd('v','l',neigmx,ca,1,1,ndesca,eig,cz,1,1,ndesca,cwork,lwork,rwork,lrwork,iwork,liwork,ierr)
    deallocate(cwork,rwork,iwork)
    sval(:,ns1,nk)=eig(:)
  end if
  call mpi_bcast(ierr,1,mpi_integer,0,mpicom_space,mpij)
  if (ierr /= 0) call stopp('rayleighritz: error occurs in pzheevd')
  call mpi_bcast(sval(1,ns1,nk),neigmx,mpi_double_precision,0,mpicom_space,mpij)
! ====================================
  call endtime(ndisp,stime0,'[TI_sub] sub-space diag. :')

  call starttime(stime0)
  deallocate(ca)
  allocate(ca(mstep,mstep),czall(mstep,mstep))
! ==========  update vectors  ==========
  do ns= ns1,ns1+ncol-1
!$omp parallel default(shared) private(i)
!$omp do
    do i=1,ncpx*ncpy*ncpz*neigmx
      hsvcm(i,1,1,1,ns-ns1+1)=dcmplx(0.0d0,0.0d0)
      sssvcm(i,1,1,1,ns-ns1+1)=dcmplx(0.0d0,0.0d0)
    end do
!$omp do
    do i=1,nprjmx*num_ppcell*neigmx
      cspsep_rr(i,1,1,ns-ns1+1)=dcmplx(0.0d0,0.0d0)
    end do
!$omp end parallel
    k=0
    do while (k<neigmx)
      if (neigmx-k>=mstep) then
        j=0
        do while (j<neigmx)
!$omp parallel default(shared) private(i)
!$omp do
          do i=1,mstep*mstep
            czall(i,1)=dcmplx(0.0d0,0.0d0)
          end do
!$omp end parallel
          if (neigmx-j>=mstep) then
            do kk=1,mstep,mb
              do jj=1,mstep,mb
                if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) &
                    czall(jj:jj+mb-1,kk:kk+mb-1)=cz(jrow(j+jj):jrow(j+jj)+mb-1,jcol(k+kk):jcol(k+kk)+mb-1)
              end do
            end do
            call mpi_allreduce(czall,ca,mstep*mstep,mpi_double_complex,mpi_sum,mpicom_space,mpij)
            call zgemm('n','n',ncpx*ncpy*ncpz,mstep,mstep,dcmplx(1.0d0,0.0d0),sveccm(1,1,1,j+1,ns,nk) &
                      ,ncpx*ncpy*ncpz,ca,mstep,dcmplx(1.0d0,0.0d0),hsvcm(1,1,1,k+1,ns-ns1+1),ncpx*ncpy*ncpz)
            call zgemm('n','n',ncpx*ncpy*ncpz,mstep,mstep,dcmplx(1.0d0,0.0d0),ssvcm(1,1,1,j+1,ns,nk) &
                      ,ncpx*ncpy*ncpz,ca,mstep,dcmplx(1.0d0,0.0d0),sssvcm(1,1,1,k+1,ns-ns1+1),ncpx*ncpy*ncpz)
            call zgemm('n','n',nprjmx*num_ppcell,mstep,mstep,dcmplx(1.0d0,0.0d0),cspsep(1,1,j+1,ns,nk) &
                      ,nprjmx*num_ppcell,ca,mstep,dcmplx(1.0d0,0.0d0),cspsep_rr(1,1,k+1,ns-ns1+1),nprjmx*num_ppcell)
            j=j+mstep
          else
            do kk=1,mstep,mb
              do jj=1,neigmx-j
                if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) &
                  czall(jj,kk:kk+mb-1)=cz(jrow(j+jj),jcol(k+kk):jcol(k+kk)+mb-1)
              end do
            end do
            do kk=1,mstep
              call mpi_allreduce(czall(1,kk),ca(1,kk),neigmx-j,mpi_double_complex,mpi_sum,mpicom_space,mpij)
            end do
            call zgemm('n','n',ncpx*ncpy*ncpz,mstep,neigmx-j,dcmplx(1.0d0,0.0d0),sveccm(1,1,1,j+1,ns,nk) &
                      ,ncpx*ncpy*ncpz,ca,mstep,dcmplx(1.0d0,0.0d0),hsvcm(1,1,1,k+1,ns-ns1+1),ncpx*ncpy*ncpz)
            call zgemm('n','n',ncpx*ncpy*ncpz,mstep,neigmx-j,dcmplx(1.0d0,0.0d0),ssvcm(1,1,1,j+1,ns,nk) &
                      ,ncpx*ncpy*ncpz,ca,mstep,dcmplx(1.0d0,0.0d0),sssvcm(1,1,1,k+1,ns-ns1+1),ncpx*ncpy*ncpz)
            call zgemm('n','n',nprjmx*num_ppcell,mstep,neigmx-j,dcmplx(1.0d0,0.0d0),cspsep(1,1,j+1,ns,nk) &
                      ,nprjmx*num_ppcell,ca,mstep,dcmplx(1.0d0,0.0d0),cspsep_rr(1,1,k+1,ns-ns1+1),nprjmx*num_ppcell)
            j=neigmx
          end if
        end do
        k=k+mstep
      else
        j=0
        do while (j<neigmx)
!$omp parallel default(shared) private(i)
!$omp do
          do i=1,mstep*mstep
            czall(i,1)=dcmplx(0.0d0,0.0d0)
          end do
!$omp end parallel
          if (neigmx-j>=mstep) then
            do kk=1,neigmx-k
              do jj=1,mstep,mb
                if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) &
                  czall(jj:jj+mb-1,kk)=cz(jrow(j+jj):jrow(j+jj)+mb-1,jcol(k+kk))
              end do
            end do
            call mpi_allreduce(czall(1,1),ca(1,1),mstep*(neigmx-k),mpi_double_complex,mpi_sum,mpicom_space,mpij)
            call zgemm('n','n',ncpx*ncpy*ncpz,neigmx-k,mstep,dcmplx(1.0d0,0.0d0),sveccm(1,1,1,j+1,ns,nk) &
                      ,ncpx*ncpy*ncpz,ca,mstep,dcmplx(1.0d0,0.0d0),hsvcm(1,1,1,k+1,ns-ns1+1),ncpx*ncpy*ncpz)
            call zgemm('n','n',ncpx*ncpy*ncpz,neigmx-k,mstep,dcmplx(1.0d0,0.0d0),ssvcm(1,1,1,j+1,ns,nk) &
                      ,ncpx*ncpy*ncpz,ca,mstep,dcmplx(1.0d0,0.0d0),sssvcm(1,1,1,k+1,ns-ns1+1),ncpx*ncpy*ncpz)
            call zgemm('n','n',nprjmx*num_ppcell,neigmx-k,mstep,dcmplx(1.0d0,0.0d0),cspsep(1,1,j+1,ns,nk) &
                      ,nprjmx*num_ppcell,ca,mstep,dcmplx(1.0d0,0.0d0),cspsep_rr(1,1,k+1,ns-ns1+1),nprjmx*num_ppcell)
            j=j+mstep
          else
            do kk=1,neigmx-k
              do jj=1,neigmx-j
                if (myr_space == lrow(j+jj)+lcol(k+kk)*nprow) czall(jj,kk)=cz(jrow(j+jj),jcol(k+kk))
              end do
            end do
             do kk=1,neigmx-k
              call mpi_allreduce(czall(1,kk),ca(1,kk),neigmx-j,mpi_double_complex,mpi_sum,mpicom_space,mpij)
            end do
            call zgemm('n','n',ncpx*ncpy*ncpz,neigmx-k,neigmx-j,dcmplx(1.0d0,0.0d0),sveccm(1,1,1,j+1,ns,nk) &
                      ,ncpx*ncpy*ncpz,ca,mstep,dcmplx(1.0d0,0.0d0),hsvcm(1,1,1,k+1,ns-ns1+1),ncpx*ncpy*ncpz)
            call zgemm('n','n',ncpx*ncpy*ncpz,neigmx-k,neigmx-j,dcmplx(1.0d0,0.0d0),ssvcm(1,1,1,j+1,ns,nk) &
                      ,ncpx*ncpy*ncpz,ca,mstep,dcmplx(1.0d0,0.0d0),sssvcm(1,1,1,k+1,ns-ns1+1),ncpx*ncpy*ncpz)
            call zgemm('n','n',nprjmx*num_ppcell,neigmx-k,neigmx-j,dcmplx(1.0d0,0.0d0),cspsep(1,1,j+1,ns,nk) &
                      ,nprjmx*num_ppcell,ca,mstep,dcmplx(1.0d0,0.0d0),cspsep_rr(1,1,k+1,ns-ns1+1),nprjmx*num_ppcell)
            j=neigmx
          end if
        end do
        k=neigmx
      end if
    end do
!$omp parallel default(shared) private(i)
!$omp do
    do i=1,ncpx*ncpy*ncpz*neigmx
      sveccm(i,1,1,1,ns,nk)=hsvcm(i,1,1,1,ns-ns1+1)
      ssvcm(i,1,1,1,ns,nk)=sssvcm(i,1,1,1,ns-ns1+1)
    end do
!$omp do
    do i=1,nprjmx*num_ppcell*neigmx
      cspsep(i,1,1,ns,nk)=cspsep_rr(i,1,1,ns-ns1+1)
    end do
!$omp end parallel
  end do !ns
  deallocate(ca,czall)
  deallocate(sssvcm,cspsep_rr)
! ======================================
  call endtime(ndisp,stime0,'[TI_sub] update vectors :')

  deallocate(cz,eig)

end subroutine rayleighritz_c


subroutine rayleighritzsoc_r( &
 natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,nspv, & ! <
 num_list,ncpx,ncpy,ncpz,                               & ! <
 nperi,nk,nf,ndisp,                                     & ! <
 key_natpri_in,key_natpri_inps,key_soc_calc,            & ! <
 dx,dy,dz,                                              & ! <
 nprj,natinf,lstvec2,indspe,natpri,naps,natsoc,ntyppp,  & ! <
 dij,dijsoc,veff,vnlocp,                                & ! <
 rspsep,cspsepsoc,svecre,svecsoc,                       & ! X
 hsvre,                                                 & ! K
 svalsoc)                                                 ! >
use mod_mpi
use mod_stopp
use mod_disptime
use mod_overlap_finitedifference, only:overlap_finitedifference_init,overlap_finitedifference_final, &
                                       overlap_finitedifference_r,overlap_fdcheck_r
implicit none
integer,    intent(in)   ::natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,nspv
integer,    intent(in)   ::num_list,ncpx,ncpy,ncpz
integer,    intent(in)   ::nperi,nk,nf,ndisp
integer,    intent(in)   ::key_natpri_in,key_natpri_inps,key_soc_calc
real*8,     intent(in)   ::dx,dy,dz
integer,    intent(in)   ::nprj(num_spe),natinf(natom),lstvec2(num_list,num_ppcell)
integer,    intent(in)   ::indspe(natom),natpri(natom),naps(natom),natsoc(natom),ntyppp(num_spe)
real*8,     intent(in)   ::dij(nprjmx,nprjmx,nums,natom)
real*8,     intent(in)   ::dijsoc(nprjmx,nprjmx,3,natom)
real*8,     intent(in)   ::veff(ncpx,ncpy,ncpz,nspv)
real*8,     intent(in)   ::vnlocp(num_list,nprjmx,num_ppcell)
real*8,     intent(inout)::svecre(ncpx,ncpy,ncpz,neigmx,nums,numk)
complex*16, intent(out)  ::svecsoc(ncpx,ncpy,ncpz,2*neigmx,2,numk)
real*8,     intent(inout)::rspsep(nprjmx,num_ppcell,neigmx,nums,numk)
complex*16, intent(out)  ::cspsepsoc(nprjmx,num_ppcell,2*neigmx,2,numk)
real*8,     intent(inout)::hsvre(ncpx,ncpy,ncpz,neigmx)
real*8,     intent(out)  ::svalsoc(2*neigmx,1,numk)
integer ::ierr,i,ns,ns0,ns2,l,l2
real*8 ::stime0
real*8,    allocatable::raqr(:,:)
complex*16,allocatable::caqrsoc(:,:,:,:)
complex*16,allocatable::cwork(:)
real*8,    allocatable::rwork(:)
  allocate(raqr(neigmx,neigmx))
  allocate(caqrsoc(neigmx,2,neigmx,2)) 
  allocate(cwork(neigmx*max(6,nprjmx*num_ppcell)),rwork(2*neigmx*3))

  if (myrank_glbl==0) write(ndisp,*) 'diagonalization in the soc subspace was implemented.'

  call starttime(stime0)
! ==========  setup matrix  ==========
  do ns=1,nums
    raqr=0.0d0
    call rayleighritz_hamiltonian_r( &
     natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,nspv, & ! <
     num_list,ncpx,ncpy,ncpz,                               & ! <
     nperi,nk,ns,nf,                                        & ! <
     key_natpri_in,key_natpri_inps,                         & ! <
     dx,dy,dz,                                              & ! <
     nprj,natinf,lstvec2,indspe,natpri,naps,                & ! <
     dij,veff,vnlocp,svecre,rspsep,                         & ! <
     hsvre)                                                   ! >
    call dgemm('t','n',neigmx,neigmx,ncpx*ncpy*ncpz,dx*dy*dz,svecre(1,1,1,1,ns,nk), &
               ncpx*ncpy*ncpz,hsvre,ncpx*ncpy*ncpz,1.0d0,raqr,neigmx)
    do l= 1,neigmx  
      caqrsoc(l:neigmx,ns,l,ns)= dcmplx(raqr(l:neigmx,l),0.0d0)
    enddo  
  end do
!$omp parallel default(shared)
  call rayleighritz_matrixsoc_r( &
   natom,num_spe,num_ppcell,nprjmx,nums,neigmx,numk,nk, & ! <
   key_natpri_in,key_soc_calc,                          & ! <
   indspe,natpri,naps,natsoc,nprj,                      & ! <
   rspsep,dijsoc,                                       & ! <
   caqrsoc)                                               ! X
!$omp end parallel
  do ns=1,nums
    do i= 1,neigmx    
      call mpi_allreduce(caqrsoc(1,1,i,ns),cwork(1),2*neigmx,mpi_double_complex,mpi_sum,mpicom_space,mpij)
      caqrsoc(:,1,i,ns)= cwork(1:neigmx)
      caqrsoc(:,2,i,ns)= cwork(neigmx+1:2*neigmx)
    end do 
  end do 
! ====================================
  call endtime(ndisp,stime0,'[TI_sub] soc setup matrix :')

  call starttime(stime0)
! ==========  eigenproblem  ==========
  call zheev('v','l',2*neigmx,caqrsoc,2*neigmx,svalsoc(1,1,nk),cwork,2*neigmx*3,rwork,ierr)
  if (ierr/=0) call stopp('rayleighritz_diag_c: error in soc matrix diagonalization')
! ====================================
  call endtime(ndisp,stime0,'[TI_sub] soc sub-space diag. :')

  call starttime(stime0)
! ==========  update vectors  ==========
  do ns= 1,2
    do l= 1,neigmx
      do ns2= 1,2
        ns0= min(nums,ns2)
        !$omp parallel default(shared) private(i,l2)
        !$omp do
        do i= 1,ncpx*ncpy*ncpz
          svecsoc(i,1,1,l+neigmx*(ns-1),ns2,nk)= dcmplx(0.0d0,0.0d0)
        enddo
        !$omp do
        do i= 1,nprjmx*num_ppcell
          cspsepsoc(i,1,l+neigmx*(ns-1),ns2,nk)= dcmplx(0.0d0,0.0d0)
        enddo
        !$omp do
        do l2= 1,neigmx
          do i= 1,ncpx*ncpy*ncpz
            svecsoc(i,1,1,l+neigmx*(ns-1),ns2,nk)= &
             svecsoc(i,1,1,l+neigmx*(ns-1),ns2,nk) +svecre(i,1,1,l2,ns0,nk)*caqrsoc(l2,ns2,l,ns)
          enddo
          do i= 1,nprjmx*num_ppcell
            cspsepsoc(i,1,l+neigmx*(ns-1),ns2,nk)= &
             cspsepsoc(i,1,l+neigmx*(ns-1),ns2,nk) +rspsep(i,1,l2,ns0,nk)*caqrsoc(l2,ns2,l,ns)
          enddo
        enddo
        !$omp end parallel
      enddo
    enddo
  enddo
! ======================================
  call endtime(ndisp,stime0,'[TI_sub] soc update vectors :')

  deallocate(raqr,caqrsoc,rwork,cwork)

end subroutine rayleighritzsoc_r


subroutine rayleighritzsoc_c( &
 nso,natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,ncol,nspv, & ! <
 num_list,ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                   & ! <
 nperi,nk,nf,ndisp,                                              & ! <
 key_natpri_in,key_natpri_inps,key_soc_calc,                     & ! <
 dx,dy,dz,xmax,ymax,zmax,skpxx,skpyy,skpzz,                      & ! <
 nprj,natinf,lstvec2,indspe,natpri,naps,natsoc,ntyppp,           & ! <
 lstx,lsty,lstz,natx,naty,natz,                                  & ! <
 dij,dijsoc,veff,vnlocp,                                         & ! <
 cspsep,cspsepsoc,sveccm,svecsoc,                                & ! X
 hsvcm,                                                          & ! K
 svalsoc)                                                          ! >
use mod_mpi
use mod_stopp
use mod_disptime
use mod_overlap_finitedifference, only:overlap_finitedifference_init,overlap_finitedifference_final, &
                                       overlap_finitedifference_c,overlap_fdcheck_c
implicit none
integer,    intent(in)   ::nso,natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,ncol,nspv
integer,    intent(in)   ::num_list,ncpx,ncpy,ncpz,npxmax,npymax,npzmax
integer,    intent(in)   ::nperi,nk,nf,ndisp
integer,    intent(in)   ::key_natpri_in,key_natpri_inps,key_soc_calc
real*8,     intent(in)   ::dx,dy,dz,xmax,ymax,zmax,skpxx(numk),skpyy(numk),skpzz(numk)
integer,    intent(in)   ::nprj(num_spe),natinf(natom),lstvec2(num_list,num_ppcell)
integer,    intent(in)   ::indspe(natom),natpri(natom),naps(natom),natsoc(natom),ntyppp(num_spe)
integer,    intent(in)   ::lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell)
integer,    intent(in)   ::natx(natom),naty(natom),natz(natom)
real*8,     intent(in)   ::dij(nprjmx,nprjmx,nums*ncol,natom)
real*8,     intent(in)   ::dijsoc(nprjmx,nprjmx,3,natom)
real*8,     intent(in)   ::veff(ncpx,ncpy,ncpz,nspv)
real*8,     intent(in)   ::vnlocp(num_list,nprjmx,num_ppcell)
complex*16, intent(inout)::sveccm(ncpx,ncpy,ncpz,neigmx,nums,numk)
complex*16, intent(out)  ::svecsoc(ncpx,ncpy,ncpz,2*neigmx,1,numk)
complex*16, intent(inout)::cspsep(nprjmx,num_ppcell,neigmx,nums,numk)
complex*16, intent(out)  ::cspsepsoc(nprjmx,num_ppcell,2*neigmx,1,numk)
complex*16, intent(inout)::hsvcm(ncpx,ncpy,ncpz,neigmx,ncol)
real*8,     intent(out)  ::svalsoc(2*neigmx,1,numk)
integer   ::ierr,i,ns,ns0,ns2,l
real*8    ::stime0
complex*16::ca,cb
complex*16,allocatable::caqr(:,:)
complex*16,allocatable::caqrsoc(:,:,:,:)
complex*16,allocatable::cwork(:)
real*8,    allocatable::rwork(:)
  allocate(caqrsoc(neigmx,2,neigmx,2)) 
  allocate(caqr(neigmx,neigmx))
  allocate(cwork(neigmx*max(6,nprjmx*num_ppcell)),rwork(2*neigmx*3))

  if (myrank_glbl==0) write(ndisp,*) 'diagonalization in the soc subspace was implemented.'

  call starttime(stime0)
! ==========  setup matrix  ==========
  do ns=1,nums
    caqr=dcmplx(0.0d0,0.0d0)
    call rayleighritz_hamiltonian_c( &
     nso,natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,ncol,nspv, & ! <
     num_list,ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                   & ! <
     nperi,nk,ns,nf,                                                 & ! <
     key_natpri_in,key_natpri_inps,key_soc_calc,                     & ! <
     dx,dy,dz,xmax,ymax,zmax,skpxx,skpyy,skpzz,                      & ! <
     nprj,natinf,lstvec2,indspe,natpri,naps,natsoc,                  & ! <
     lstx,lsty,lstz,natx,naty,natz,                                  & ! <
     dij,dijsoc,veff,vnlocp,sveccm,cspsep,                           & ! <
     hsvcm)                                                            ! >
    ca=dcmplx(dx*dy*dz,0.0d0)
    cb=dcmplx(1.0d0,0.0d0)
    call zgemm('c','n',neigmx,neigmx,ncpx*ncpy*ncpz,ca,sveccm(1,1,1,1,ns,nk),ncpx*ncpy*ncpz, &
               hsvcm,ncpx*ncpy*ncpz,cb,caqr,neigmx)
    do l= 1,neigmx
      caqrsoc(l:neigmx,ns,l,ns)= caqr(l:neigmx,l)
    end do
  end do
!$omp parallel default(shared)
  call rayleighritz_matrixsoc_c( &
   natom,num_spe,num_ppcell,nprjmx,nums,neigmx,numk,nk, & ! <
   key_natpri_in,key_soc_calc,                          & ! <
   indspe,natpri,naps,natsoc,nprj,                      & ! <
   cspsep,dijsoc,                                       & ! <
   caqrsoc)                                               ! X
!$omp end parallel
  do ns=1,nums
    do i= 1,neigmx
      call mpi_allreduce(caqrsoc(1,1,i,ns),cwork(1),2*neigmx,mpi_double_complex,mpi_sum,mpicom_space,mpij)
      caqrsoc(:,1,i,ns)= cwork(1:neigmx)
      caqrsoc(:,2,i,ns)= cwork(neigmx+1:2*neigmx)
    end do 
  end do 
! ====================================
  call endtime(ndisp,stime0,'[TI_sub] soc setup matrix :')

  call starttime(stime0)
! ==========  eigenproblem  ==========
  call zheev('v','l',2*neigmx,caqrsoc,2*neigmx,svalsoc(1,1,nk),cwork,2*neigmx*3,rwork,ierr)
  if (ierr/=0) call stopp('rayleighritz_diag_c: error in soc matrix diagonalization')
! ====================================
  call endtime(ndisp,stime0,'[TI_sub] soc sub-space diag. :')

  call starttime(stime0)
! ==========  update vectors  ==========
  do ns= 1,2 
    do l= 1,neigmx
      do ns2= 1,2
        ns0= min(nums,ns2)
        call zgemm('n','n',ncpx*ncpy*ncpz,1,neigmx, &
         dcmplx(1.0d0,0.0d0),sveccm(1,1,1,1,ns0,nk),ncpx*ncpy*ncpz,caqrsoc(1,ns2,l,ns),neigmx, &
         dcmplx(0.0d0,0.0d0),svecsoc(1,1,1,l+neigmx*(ns-1),ns2,nk),ncpx*ncpy*ncpz)
        call zgemm('n','n',nprjmx*num_ppcell,1,neigmx, &
         dcmplx(1.0d0,0.0d0),cspsep(1,1,1,ns0,nk),nprjmx*num_ppcell,caqrsoc(1,ns2,l,ns),neigmx, &
         dcmplx(0.0d0,0.0d0),cspsepsoc(1,1,l+neigmx*(ns-1),ns2,nk),nprjmx*num_ppcell)
      enddo
    enddo
  enddo
! ======================================
  call endtime(ndisp,stime0,'[TI_sub] soc update vectors :')

  deallocate(caqr,caqrsoc,rwork,cwork)

end subroutine rayleighritzsoc_c


subroutine rayleighritz_hamiltonian_r( &
 natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,nspv, & ! <
 num_list,ncpx,ncpy,ncpz,                               & ! <
 nperi,nk,ns,nf,                                        & ! <
 key_natpri_in,key_natpri_inps,                         & ! <
 dx,dy,dz,                                              & ! <
 nprj,natinf,lstvec2,indspe,natpri,naps,                & ! <
 dij,veff,vnlocp,svecre,rspsep,                         & ! <
 hsvre)                                                   ! >
use mod_overlap_finitedifference, only:overlap_finitedifference_init,overlap_finitedifference_final, &
                                       overlap_finitedifference_r,overlap_fdcheck_r
use mod_nonlocaloperation, only:nonlocaloperation_r_02
use mod_kslaplacian, only:kslaplacian_r,kslaplacian_initialize,kslaplacian_finalize
implicit none
integer, intent(in)  ::natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,nspv
integer, intent(in)  ::num_list,ncpx,ncpy,ncpz
integer, intent(in)  ::nperi,nk,ns,nf
integer, intent(in)  ::key_natpri_in,key_natpri_inps
real*8,  intent(in)  ::dx,dy,dz
integer, intent(in)  ::nprj(num_spe),natinf(natom),lstvec2(num_list,num_ppcell)
integer, intent(in)  ::indspe(natom),natpri(natom),naps(natom)
real*8,  intent(in)  ::rspsep(nprjmx,num_ppcell,neigmx,nums,numk)
real*8,  intent(in)  ::dij(nprjmx,nprjmx,nums,natom)
real*8,  intent(in)  ::veff(ncpx,ncpy,ncpz,nspv)
real*8,  intent(in)  ::vnlocp(num_list,nprjmx,num_ppcell)
real*8,  intent(in)  ::svecre(ncpx,ncpy,ncpz,neigmx,nums,numk)
real*8,  intent(out) ::hsvre(ncpx,ncpy,ncpz,neigmx)
integer ::l,ix,iy,iz,nsv,i
real*8,    allocatable::vre(:,:,:)
real*8,    allocatable::avr(:)

  allocate(vre(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf))
  allocate(avr(num_list))

  call overlap_finitedifference_init(ncpx,ncpy,ncpz,nf,1,0)
  call kslaplacian_initialize(nf)

  nsv= min(ns,nspv)

!$omp parallel default(shared) private(i,ix,iy,iz)
  do l= 1,neigmx
!$omp do
    do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)
      vre(-nf+i,-nf+1,-nf+1)=0.0d0
    end do
!$omp do
    do i=1,ncpx*ncpy*ncpz
      hsvre(i,1,1,l)=0.0d0
    end do
!$omp do
    do iz=1,ncpz
    do iy=1,ncpy
    do ix=1,ncpx
      vre(ix,iy,iz)=svecre(ix,iy,iz,l,ns,nk)
    end do
    end do
    end do
!$omp single
    call overlap_finitedifference_r(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,vre)
!$omp end single
    call nonlocaloperation_r_02(1,natom,num_spe,nums,nprjmx,ns,num_list,num_ppcell, & ! <
                                ncpx,ncpy,ncpz,                                     & ! <
                                key_natpri_in,key_natpri_inps,                      & ! <
                                dx,dy,dz,                                           & ! <
                                nprj,indspe,natinf,lstvec2,natpri,naps,             & ! <
                                vnlocp,rspsep(1,1,l,ns,nk),dij,                     & ! <
                                hsvre(1,1,1,l),                                     & ! X
                                avr)                                                  ! W
!$omp barrier
!$omp single
    call overlap_fdcheck_r(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,vre)
!$omp end single
    call kslaplacian_r(ncpx,ncpy,ncpz,neigmx,nf,dx,dy,dz,vre,veff(1,1,1,nsv),l,hsvre)
!$omp barrier
  end do ! l
!$omp end parallel

  deallocate( vre )
  deallocate(avr)
  call overlap_finitedifference_final
  call kslaplacian_finalize
  return
end subroutine rayleighritz_hamiltonian_r


subroutine rayleighritz_hamiltonian_c( &
 nso,natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,ncol,nspv, & ! <
 num_list,ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                   & ! <
 nperi,nk,ns1,nf,                                                & ! <
 key_natpri_in,key_natpri_inps,key_soc_calc,                     & ! <
 dx,dy,dz,xmax,ymax,zmax,skpxx,skpyy,skpzz,                      & ! <
 nprj,natinf,lstvec2,indspe,natpri,naps,natsoc,                  & ! <
 lstx,lsty,lstz,natx,naty,natz,                                  & ! <
 dij,dijsoc,veff,vnlocp,sveccm,cspsep,                           & ! <
 hsvcm)                                                            ! >
use mod_mpi
use mod_overlap_finitedifference, only:overlap_finitedifference_init,overlap_finitedifference_final, &
                                       overlap_finitedifference_c,overlap_fdcheck_c
use mod_nonlocaloperation, only:nonlocaloperation_c_02
use mod_kslaplacian, only:kslaplacian_c,kslaplacian_initialize,kslaplacian_finalize
implicit none
integer,    intent(in)  ::nso,natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,ncol,nspv
integer,    intent(in)  ::num_list,ncpx,ncpy,ncpz,npxmax,npymax,npzmax
integer,    intent(in)  ::nperi,nk,ns1,nf
integer,    intent(in)  ::key_natpri_in,key_natpri_inps,key_soc_calc
real*8,     intent(in)  ::dx,dy,dz,xmax,ymax,zmax,skpxx(numk),skpyy(numk),skpzz(numk)
integer,    intent(in)  ::nprj(num_spe),natinf(natom),lstvec2(num_list,num_ppcell)
integer,    intent(in)  ::indspe(natom),natpri(natom),naps(natom),natsoc(natom)
integer,    intent(in)  ::lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell)
integer,    intent(in)  ::natx(natom),naty(natom),natz(natom)
real*8,     intent(in)  ::dij(nprjmx,nprjmx,nums*ncol,natom)
real*8,     intent(in)  ::dijsoc(nprjmx*nso-nso+1,nprjmx*nso-nso+1,3*nso-nso+1,natom*nso-nso+1)
real*8,     intent(in)  ::veff(ncpx,ncpy,ncpz,nspv)
real*8,     intent(in)  ::vnlocp(num_list,nprjmx,num_ppcell)
complex*16, intent(in)  ::sveccm(ncpx,ncpy,ncpz,neigmx,nums,numk)
complex*16, intent(in)  ::cspsep(nprjmx,num_ppcell,neigmx,nums,numk)
complex*16, intent(out) ::hsvcm(ncpx,ncpy,ncpz,neigmx,ncol)
integer ::nvef,l,ns2,nsv,ix,iy,iz,i
complex*16,allocatable::vcm(:,:,:,:)
complex*16,allocatable::workc(:,:,:,:)
complex*16,allocatable::cpsep(:,:,:)
complex*16,allocatable::vcccm(:)
complex*16,allocatable::avc(:,:)

  allocate(vcm(-(nf-1):ncpx+nf,-(nf-1):ncpy+nf,-(nf-1):ncpz+nf,ncol))
  allocate(workc(ncpx,ncpy,ncpz,ncol))
  allocate(cpsep(nprjmx,num_ppcell,ncol))
  allocate(vcccm(num_list))
  allocate(avc(num_list,ncol))
  call overlap_finitedifference_init(ncpx,ncpy,ncpz,nf,ncol,1)
  call kslaplacian_initialize(nf)

  nvef=1
  if (ncol==2) nvef=nspv
  nsv= min(ns1,nspv)

!$omp parallel default(shared) private(i,ix,iy,iz)
  do l=1,neigmx
    do ns2= 1,ncol
!$omp do
      do i=1,(ncpx+2*nf)*(ncpy+2*nf)*(ncpz+2*nf)
        vcm(-nf+i,-nf+1,-nf+1,ns2)=dcmplx(0.0d0,0.0d0)
      end do
!$omp do
      do i=1,ncpx*ncpy*ncpz
        hsvcm(i,1,1,l,ns2)=dcmplx(0.0d0,0.0d0)
      end do
!$omp do
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        vcm(ix,iy,iz,ns2)=sveccm(ix,iy,iz,l,ns1+ns2-1,nk)
      end do
      end do
      end do
    end do ! ns2
!$omp single
    call overlap_finitedifference_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,ncol,vcm)
!$omp end single
    do ns2= 1,ncol
!$omp do
      do i=1,nprjmx*num_ppcell
        cpsep(i,1,ns2)= cspsep(i,1,l,ns1+ns2-1,nk)
      end do
    end do
    call nonlocaloperation_c_02(1,nso,natom,num_spe,neigmx,nums,nprjmx,ncol,ns1,l,num_list,num_ppcell, & ! <
                                ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                                   & ! <
                                key_natpri_in,key_natpri_inps,key_soc_calc,                            & ! <
                                dx,dy,dz,skpxx(nk),skpyy(nk),skpzz(nk),                                & ! <
                                nprj,indspe,natinf,lstvec2,natpri,naps,natsoc,                         & ! <
                                lstx,lsty,lstz,natx,naty,natz,                                         & ! <
                                vnlocp,cpsep,dij,dijsoc,                                               & ! <
                                hsvcm,                                                                 & ! X
                                vcccm,avc)                                                               ! W
!$omp barrier
!$omp single
    call overlap_fdcheck_c(nperi,ncpx,ncpy,ncpz,nf,nf-1,nf,ncol,vcm)
!$omp end single
    call kslaplacian_c(ncpx,ncpy,ncpz,neigmx,ncol,nvef,nf,dx,dy,dz,xmax,ymax,zmax,skpxx(nk),skpyy(nk),skpzz(nk) &
                       ,vcm,veff(1,1,1,nsv),l,hsvcm,workc)
!$omp barrier
  end do ! l
!$omp end parallel

  deallocate(vcm,cpsep)
  deallocate(workc)
  deallocate(vcccm,avc)
  call overlap_finitedifference_final
  call kslaplacian_finalize
  return
end subroutine rayleighritz_hamiltonian_c


subroutine rayleighritz_matrixsoc_r( &
 natom,num_spe,num_ppcell,nprjmx,nums,neigmx,numk,nk, & ! <
 key_natpri_in,key_soc_calc,                          & ! <
 indspe,natpri,naps,natsoc,nprj,                      & ! <
 rspsep,dijsoc,                                       & ! <
 caqrsoc)                                                ! X
implicit none
integer,   intent(in)   ::natom,num_spe,num_ppcell,nprjmx,nums,neigmx,numk,nk
integer,   intent(in)   ::key_natpri_in,key_soc_calc
integer,   intent(in)   ::indspe(natom),natpri(natom),naps(natom),natsoc(natom),nprj(num_spe)
real*8,    intent(in)   ::rspsep(nprjmx,num_ppcell,neigmx,nums,numk)
real*8,    intent(in)   ::dijsoc(nprjmx,nprjmx,3,natom)
complex*16,intent(inout)::caqrsoc(neigmx,2,neigmx,2)
integer :: na,iaps,iprj,ns,ns0,l2,l1,i2,i1
real*8  :: sig

  !$omp do
  do l2= 1,neigmx
    caqrsoc(1:neigmx,2,l2,1)= dcmplx(0.0d0,0.0d0)
  enddo 

  do na= 1,natom
  if ((natpri(na)==key_natpri_in) .and. (natsoc(na)==key_soc_calc)) then
    iaps= naps(na)
    iprj= nprj(indspe(na))
    do ns= 2,1,-1
      ns0= min(ns,nums)
      sig= dble(3-2*ns)
      !$omp do
      do l2= 1,neigmx
      do l1= l2,neigmx
        caqrsoc(l1,ns,l2,ns)= caqrsoc(l1,ns0,l2,ns0)
        do i2= 1,iprj
        do i1= 1,iprj
          caqrsoc(l1,ns,l2,ns)= caqrsoc(l1,ns,l2,ns) &
           + dcmplx( 0.0d0, rspsep(i1,iaps,l1,ns0,nk)*rspsep(i2,iaps,l2,ns0,nk) *sig*dijsoc(i1,i2,1,na) )
        enddo
        enddo
      enddo
      enddo
    enddo
    !$omp do
    do l2= 1,neigmx
    do l1= 1,neigmx
      do i2= 1,iprj
      do i1= 1,iprj
        caqrsoc(l1,2,l2,1)= caqrsoc(l1,2,l2,1) &
         + dcmplx( rspsep(i1,iaps,l1,nums,nk)*rspsep(i2,iaps,l2,1,nk), 0.0d0 ) &
          *dcmplx( dijsoc(i1,i2,2,na), dijsoc(i1,i2,3,na) )
      enddo
      enddo
    enddo
    enddo
  endif
  enddo ! na

end subroutine rayleighritz_matrixsoc_r


subroutine rayleighritz_matrixsoc_c( &
 natom,num_spe,num_ppcell,nprjmx,nums,neigmx,numk,nk, & ! <
 key_natpri_in,key_soc_calc,                          & ! <
 indspe,natpri,naps,natsoc,nprj,                      & ! <
 cspsep,dijsoc,                                       & ! <
 caqrsoc)                                               ! X
implicit none
integer,   intent(in)   ::natom,num_spe,num_ppcell,nprjmx,nums,neigmx,numk,nk
integer,   intent(in)   ::key_natpri_in,key_soc_calc
integer,   intent(in)   ::indspe(natom),natpri(natom),naps(natom),natsoc(natom),nprj(num_spe)
complex*16,intent(in)   ::cspsep(nprjmx,num_ppcell,neigmx,nums,numk)
real*8,    intent(in)   ::dijsoc(nprjmx,nprjmx,3,natom)
complex*16,intent(inout)::caqrsoc(neigmx,2,neigmx,2)
integer :: na,iaps,iprj,ns,ns0,l2,l1,i2,i1
real*8  :: sig

  !$omp do
  do l2= 1,neigmx
    caqrsoc(1:neigmx,2,l2,1)= dcmplx(0.0d0,0.0d0)
  enddo 

  do na= 1,natom
  if ((natpri(na)==key_natpri_in) .and. (natsoc(na)==key_soc_calc)) then
    iaps= naps(na)
    iprj= nprj(indspe(na))
    do ns= 2,1,-1
      ns0= min(ns,nums)
      sig= dble(3-2*ns)
      !$omp do 
      do l2= 1,neigmx
      do l1= l2,neigmx
        caqrsoc(l1,ns,l2,ns)= caqrsoc(l1,ns0,l2,ns0)
        do i2= 1,iprj
        do i1= 1,iprj
          caqrsoc(l1,ns,l2,ns)= caqrsoc(l1,ns,l2,ns) & 
           + dconjg(cspsep(i1,iaps,l1,ns0,nk))*cspsep(i2,iaps,l2,ns0,nk) &
            *dcmplx( 0.0d0, sig*dijsoc(i1,i2,1,na) ) 
        enddo
        enddo
      enddo 
      enddo 
    enddo
    !$omp do 
    do l2= 1,neigmx
    do l1= 1,neigmx
      do i2= 1,iprj
      do i1= 1,iprj
        caqrsoc(l1,2,l2,1)= caqrsoc(l1,2,l2,1) &
         + dconjg(cspsep(i1,iaps,l1,nums,nk))*cspsep(i2,iaps,l2,1,nk) &
          *dcmplx( dijsoc(i1,i2,2,na), dijsoc(i1,i2,3,na) )
      enddo
      enddo
    enddo
    enddo
  endif
  enddo ! na

end subroutine rayleighritz_matrixsoc_c

end module
 
