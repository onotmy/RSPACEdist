! **********  spinflip.f90 06/28/2013-01  **********

subroutine tools_spinflip( &
 nrc,nums,nspv,ncpx,ncpy,ncpz,nradmx,npoint,num_atcell,neigmx,numk, & ! <
 svecre,sveccm,rhosmt,rhotrur,rhosmtr)                              ! >
use mod_stopp
implicit none
integer,   intent(in)   ::nrc,nums,nspv,ncpx,ncpy,ncpz,nradmx,npoint,num_atcell,neigmx,numk
real*8,    intent(inout)::svecre(ncpx*(1-nrc)+nrc,ncpy*(1-nrc)+nrc,ncpz*(1-nrc)+nrc, &
                                 neigmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk*(1-nrc)+nrc)
complex*16,intent(inout)::sveccm(ncpx*nrc-nrc+1,ncpy*nrc-nrc+1,ncpz*nrc-nrc+1, &
                                 neigmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1)
real*8,    intent(inout)::rhosmt(ncpx,ncpy,ncpz,nspv)
real*8,    intent(inout)::rhotrur(nradmx,npoint,nspv,num_atcell),rhosmtr(nradmx,npoint,nspv,num_atcell)
integer::nk,na 
real*8,    allocatable::rbuf(:,:,:,:)
complex*16,allocatable::cbuf(:,:,:)

  if (nums/=0) then
    if (nums/=2) call stopp('tools_spinflip: nums /= 2')
    if (nrc==0) then
      allocate(rbuf(ncpx,ncpy,ncpz,neigmx))
      do nk= 1,numk
        rbuf(:,:,:,:)       = svecre(:,:,:,:,1,nk)
        svecre(:,:,:,:,1,nk)= svecre(:,:,:,:,2,nk)
        svecre(:,:,:,:,2,nk)= rbuf(:,:,:,:)
      enddo
      deallocate(rbuf)
    else
      allocate(cbuf(ncpx,ncpy,ncpz,neigmx))
      do nk= 1,numk
        cbuf(:,:,:,:)          = sveccm(:,:,:,:,1   ,nk)
        sveccm(:,:,:,:,1   ,nk)= sveccm(:,:,:,:,nums,nk)
        sveccm(:,:,:,:,nums,nk)= cbuf(:,:,:,:)
      enddo
      deallocate(cbuf)
    endif
  endif

  if (nspv/=0) then
    if (nspv/=2) call stopp('tools_spinflip: nspv /= 2')
    allocate(rbuf(ncpx,ncpy,ncpz,1))
    rbuf(:,:,:,1)     = rhosmt(:,:,:,1   )
    rhosmt(:,:,:,1   )= rhosmt(:,:,:,nums)
    rhosmt(:,:,:,nums)= rbuf(:,:,:,1)
    deallocate(rbuf)
    allocate(rbuf(nradmx,npoint,1,1))
    do na= 1,num_atcell
      rbuf(:,:,1,1)       = rhotrur(:,:,1   ,na)
      rhotrur(:,:,1   ,na)= rhotrur(:,:,nspv,na)
      rhotrur(:,:,nspv,na)= rhobuf(:,:,1,1)
      rbuf(:,:,1,1)       = rhosmtr(:,:,1   ,na)
      rhosmtr(:,:,1   ,na)= rhosmtr(:,:,nspv,na)
      rhosmtr(:,:,nspv,na)= rhobuf(:,:,1,1)
    enddo
    deallocate(rbuf)
  endif 

end subroutine toops_spinflip 

