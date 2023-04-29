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
! **********  read_ppfile8f.F90 12/05/2022-01  **********

module mod_read_ppfile
implicit none
contains


subroutine read_ppfile( &
 chdir,cexco,                                 & ! <
 ndisp,num_spe,nradmx,nprmx,lmx,              & ! <
 key_pp_ncps,key_pp_paw,                      & ! <
 mspe,                                        & ! <
 gmaxqp,                                      & ! >
 ntyppp,nradps,nradct,npr,lpmx,nrprj,         & ! >
 cp,radial,dradial,                           & ! >
 potc,awf,pwf,rhocore,rhopcc,grwein,          & ! >
 sss0,akv0,prj,gmax,aeeig,aepot)                ! >
use mod_mpi
use mod_stopp
use mod_tools, only: tools_bisec,tools_deri_grdch
use mod_vxpot
implicit none
character,intent(in)::chdir*200,cexco*7
integer,intent(in)::ndisp,num_spe,nradmx,nprmx,lmx
integer,intent(in)::key_pp_ncps,key_pp_paw
integer,intent(in)::mspe(num_spe)
real*8, intent(out)::gmaxqp
integer,intent(out)::ntyppp(num_spe)
integer,intent(out)::nradps(num_spe),nradct(num_spe)
integer,intent(out)::npr(0:lmx-1,num_spe),lpmx(num_spe)
integer,intent(out)::nrprj(num_spe)
real*8, intent(out)::cp(8,num_spe)
real*8, intent(out)::radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8, intent(out)::potc(nradmx,num_spe),awf(nradmx,nprmx*lmx,num_spe),pwf(nradmx,nprmx*lmx,num_spe)
real*8, intent(out)::rhocore(nradmx,num_spe),rhopcc(nradmx,num_spe)
real*8, intent(out)::grwein(15)
real*8, intent(out)::sss0(nprmx*lmx,nprmx*lmx,num_spe),akv0(nprmx*lmx,nprmx*lmx,num_spe),prj(nradmx,nprmx*lmx,num_spe)
real*8, intent(out)::gmax(num_spe)
real*8, intent(out)::aeeig(nprmx*lmx,num_spe),aepot(nradmx,num_spe)
real*8, parameter :: r_eps=1.0d-15

integer  :: ispe,i,j,ir,ierr
real*8   :: ar0,ar1,ar2,ar3,ar4,ar5,ar6,uabgr,uggabg,uggr
character:: cpsfiles*13,cfilename*7,chara*1,cpptype*4,cherror*28,fname*200
real*8   :: pi
integer, allocatable:: irprj(:)
real*8, allocatable::ppp(:),vx(:),ex(:),drr(:),ddrr(:),abgr(:),ggabg(:),ggr(:)
allocate(irprj(nprmx*lmx),ppp(nradmx),vx(nradmx),ex(nradmx),drr(nradmx),ddrr(nradmx),abgr(nradmx),ggabg(nradmx),ggr(nradmx))
pi=dacos(-1.0d0)

! ==========  zero clear  ==========
  sss0(:,:,:)=0.0d0
  akv0(:,:,:)=0.0d0
  prj(:,:,:)=0.0d0
  npr(:,:)=0
  radial(:,:)=0.0d0
  dradial(:,:)=0.0d0
  potc(:,:)=0.0d0
  awf(:,:,:)=0.0d0
  pwf(:,:,:)=0.0d0
  rhocore(:,:)=0.0d0
  rhopcc(:,:)=0.0d0
! ==================================

  if (myrank_glbl .eq. 0) then

  gmaxqp=2.0d0 !this should be var gmaxqp as read in readpa
  grwein( 1)= 4.49d0
  grwein( 2)= 5.76d0
  grwein( 3)= 6.99d0
  grwein( 4)= 8.18d0
  grwein( 5)= 9.36d0
  grwein( 6)=10.51d0
  grwein( 7)=11.66d0
  grwein( 8)=12.79d0
  grwein( 9)=13.92d0
  grwein(10)=15.03d0
  grwein(11)=16.14d0
  grwein(12)=17.25d0
  grwein(13)=18.35d0
  grwein(14)=19.45d0
  grwein(15)=20.54d0

  do ispe=1,num_spe
!   ----------  prepeare string for error messages  ----------
    write(cherror,fmt='("error in pp.-file for N=",i2,":",1x)') mspe(ispe)
!   ==========  read pseusopotentials and wave functions  ==========
    cpsfiles(1:6)='pspaw/'
    write(cfilename,fmt='("paw.",i3.3)') mspe(ispe)
    cpsfiles(7:13)=cfilename
    fname=cpsfiles
    if (len_trim(chdir) > 0) fname=trim(chdir)//'/'//fname
    open(10,file=fname)
!   ----------  read header  ----------
    read(10,*,err=9999) chara
    read(10,*,err=9999) chara
    read(10,*,err=9999) chara
    read(10,*,err=9999) chara
    read(10,*,err=9999) chara
    read(10,*,err=9999) chara
    read(10,*,err=9999) chara
    read(10,*,err=9999) chara
    read(10,*,err=9999) chara
    read(10,*,err=9999) chara
!   ----------  mesh type: exponential  ----------
    read(10,*,err=9999) chara
    read(10,*,err=9999) cpptype
    cpptype= adjustl(cpptype)
    do j= 1,len(cpptype)
      if ((cpptype(j:j)>='A').and.(cpptype(j:j)<='Z')) &
       cpptype(j:j)=char(ichar(cpptype(j:j))+ichar('a')-ichar('A'))
    end do
    select case (trim(cpptype))
    case ('ncps')
      ntyppp(ispe)= key_pp_ncps
    case ('paw' )
      ntyppp(ispe)= key_pp_paw
    case default
      call stopp(cherror//'type of pseudopotential is unknown')
    end select
    read(10,*,err=9999) chara
    read(10,*,err=9999) j,i
    if (j .ne. mspe(ispe)) call stopp(cherror//'wrong atomic number')
    cp(1,ispe)=dfloat(i)
    read(10,*,err=9999) chara
    read(10,*,err=9999) i
    lpmx(ispe)= i-1
    read(10,*,err=9999) chara
    read(10,*,err=9999) i
    if ((lpmx(ispe)>lmx-1).or.(i>nprmx).or.(min(lpmx(ispe)+1,i)<0)) &
     call stopp(cherror//'number of projectors out of bounds')
    nrprj(ispe)= i*(lpmx(ispe)+1)
    do j= 1,nrprj(ispe)
      irprj(j)= ( (j-1)/i )*nprmx +mod(j-1,i) +1
    end do
!   ----------  # of dyadics for 0<=l<=lpmx  ----------
    read(10,*,err=9999) chara
    npr(:,ispe)= 0
    do j=0,lpmx(ispe)
      read(10,*,err=9999) i,npr(j,ispe)
      if (( i/=j ).or.( npr(j,ispe)<0 ).or.( npr(j,ispe)*(lpmx(ispe)+1)>nrprj(ispe) )) &
       call stopp(cherror//'unexpected number of dyadics')
    end do
!   ----------  core parameters from bhs table  ----------
    read(10,*,err=9999) chara
    read(10,*,err=9999) cp(2,ispe),cp(3,ispe),ar0,cp(4,ispe),cp(5,ispe)
    read(10,*,err=9999) chara
    read(10,*,err=9999) i,gmax(ispe)
    read(10,*,err=9999) chara
    read(10,*,err=9999) chara
    read(10,*,err=9999) chara
    read(10,*,err=9999) chara
    read(10,*,err=9999) chara
    read(10,*,err=9999) chara
!   ----------  parameters for radial mesh  ----------
    read(10,*,err=9999) chara
    read(10,*,err=9999) nradps(ispe),nradct(ispe)
    if (nradct(ispe)>nradmx) call stopp(cherror//'nradmx is too small')
    if (nradps(ispe)<nradct(ispe)) call stopp(cherror//'nradps < nradct')
!   ----------  (PS-radial mesh) r, drdi, vloc,vhxc  ----------
    read(10,*,err=9999) chara
    do ir=1,nradps(ispe)
      read(10,*,err=9999) ar0,ar1,ar2,ar3,ar4,ar5,ar6  ! (in Ry.)
      radial(ir,ispe)= ar0
      dradial(ir,ispe)=ar1
      potc(ir,ispe)=(ar2-ar4)*0.5d0
      aepot(ir,ispe)=ar5*0.5d0
      ppp(ir)=ar6
    end do
    if (ntyppp(ispe)==key_pp_paw) then
      select case (trim(cexco))
      case ('vwn')
        call vxpot_lda_vwn(nradmx,ppp,vx,ex)
      case ('pz')
        call vxpot_lda_pz(nradmx,ppp,vx,ex)
      case ('pw91')
        call tools_deri_grdch(nradmx,radial(1,ispe),dradial(1,ispe),ppp,drr,ddrr)
        do ir=1,nradps(ispe)
          abgr(ir)=dsqrt(drr(ir)**2)
          ggabg(ir)=drr(ir)**2*ddrr(ir)/abgr(ir)
          ggr(ir)=ddrr(ir)+2.0d0*drr(ir)/radial(ir,ispe)
        end do
        call vxpot_gga_pw91(nradps(ispe),ppp,abgr,ggabg,ggr,vx,ex)
      case ('pbe')
        call tools_deri_grdch(nradmx,radial(1,ispe),dradial(1,ispe),ppp,drr,ddrr)
        do ir=1,nradps(ispe)
          abgr(ir)=dsqrt(drr(ir)**2)
          ggabg(ir)=drr(ir)**2*ddrr(ir)/abgr(ir)
          ggr(ir)=ddrr(ir)+2.0d0*drr(ir)/radial(ir,ispe)
        end do
        call vxpot_gga_pbe(nradps(ispe),ppp,abgr,ggabg,ggr,vx,ex)
      case default
        call stopp ('read_ppfile:  unknown xc potential')
      end select
      do ir=1,nradps(ispe)
        potc(ir,ispe)=potc(ir,ispe)-vx(ir)
      end do
    end if
!   ----------  core-densities n_c, \tilde{n}_c  ----------
    read(10,*,err=9999) chara
    read(10,*,err=9999) i,ar0
    j=0
    do ir=1,nradps(ispe)
      read(10,*,err=9999) ar0,ar1,ar2,ar3  ! <= r,drdi,4.0*pi*r*r*n_{core-densities}(r),4.0*pi*r*r*n_{pcc}(r)
      if (ar0>r_eps) then
        ar0= ar0*ar0*4.0d0*pi
        rhopcc(ir,ispe)=ar3/ar0
        if (ir<=nradct(ispe)) rhocore(ir,ispe)=(ar2-ar3)/ar0
      else
        j= max(j,ir)
      end if
    end do
    if (j> 1) call stopp(cherror//'problem with radial grid of core densities')
    if (j==1) then
      rhopcc( 1,ispe)= rhopcc( 2,ispe)
      rhocore(1,ispe)= rhocore(2,ispe)
    end if
!   ----------  for all dyadics: ps-part.wave,ae-part.wave, projectors  ----------
    read(10,*,err=9999) chara
    do j=1,nrprj(ispe)
      read(10,*,err=9999) chara
      read(10,*,err=9999) chara
      read(10,*,err=9999) chara
      do ir=1,nradct(ispe)
        read(10,*,err=9999) ar1,ar2,ar3  ! <= r*R_{pwf}(r),r*R_{awf}(r),r*prj(r), where R is radial part of wavefunction
        pwf(ir,irprj(j),ispe)=ar1
        awf(ir,irprj(j),ispe)=ar2
        prj(ir,irprj(j),ispe)=ar3
      end do
    end do
    read(10,*,err=9999) chara
!   ---------  charge deficit matrix qvd  ----------
    read(10,*,err=9999) chara
    do j=1,nrprj(ispe)
      read(10,*,err=9999) (sss0(irprj(i),irprj(j),ispe),i=1,nrprj(ispe))
    end do
    read(10,*,err=9999) chara
!   ----------  matrix dion (kinetic  + loc.psp energy)  ----------
    read(10,*,err=9999) chara
    do j=1,nrprj(ispe)
      read(10,*,err=9999) (akv0(irprj(i),irprj(j),ispe),i=1,nrprj(ispe))  ! (in Ry.)
    end do
    akv0(:,:,ispe)=akv0(:,:,ispe)*0.5d0
    read(10,*,err=9999) chara
!   ----------  matrix dhxc (int Q_ij(r) V_hxc(r) dr)  ----------
    read(10,*,err=9999) chara
    do j=1,nrprj(ispe)
      read(10,*,err=9999) chara
    end do
    read(10,*,err=9999) chara
!   ----------  all-electron eigenvalues and potential (only for soc)  ----------
    read(10,*,err=9999) chara
    do j=0,lpmx(ispe)
      read(10,*,err=9999) i, aeeig(j*nprmx+1:j*nprmx+npr(j,ispe),ispe)
    end do
!   -----------------------------------------------------------------------------
    close(10)

    call tools_bisec(cp(1,ispe),radial(nradct(ispe),ispe),1.0d-12,cp(6,ispe),cp(7,ispe),cp(8,ispe),ierr)

!   ----------  matrix dhxc (int Q_ij(r) V_hxc(r) dr)  ----------
!   ! read(10,*,err=9999) chara
!   ! do j=1,nrprj(ispe)
!   !   read(10,*,err=9999) (aa,i=1,nrprj(ispe))
!   ! end do
!   ! read(10,*,err=9999) chara
!   ---------  AE-potential (only needed for log. derivatives))  ----------
!   ! read(10,*,err=9999) chara
!   ! read(10,*,err=9999) chara
!   ! read(10,*,err=9999) chara
!   ! do j=1,nradps(ispe)
!   !   read(10,*,err=9999) aa
!   ! end do
!   ----------  rhocps  ----------
!   ! read(10,*,err=9999) chara
!   ! read(10,*,err=9999) (rhopcc(ir,ispe),ir=1,nradps(ispe))  ! <= 4.0*pi*r*r*n_{pcc}(r)
!   !  close(10)
!
!   ----------  reduce number of partial waves (if applicable)  ----------
!   ! i= 0
!   ! lpmx(ispe)= 0
!   ! do j=0,lmx-1
!   !   i= max(i,npr(j,ispe))
!   !   if (npr(j,ispe)/=0) lpmx(ispe)= j
!   ! end do
!   ! nrprj(ispe)= i*(lpmx(ispe)+1)

!   ================================================================

  end do ! ispe

  end if ! (myrank_glbl .eq. 0)

  call mpi_bcast(ntyppp,num_spe,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(potc,nradmx*num_spe,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(awf,nradmx*nprmx*lmx*num_spe,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(pwf,nradmx*nprmx*lmx*num_spe,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(rhocore,nradmx*num_spe,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(rhopcc,nradmx*num_spe,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(nradps,num_spe,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(nradct,num_spe,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(npr,lmx*num_spe,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(radial,nradmx*num_spe,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(dradial,nradmx*num_spe,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(cp,8*num_spe,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(lpmx,num_spe,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(aeeig,nprmx*lmx*num_spe,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(aepot,nradmx*num_spe,mpi_double_precision,0,mpicom_space,mpij)

  deallocate(irprj,ppp,vx,ex,drr,ddrr,abgr,ggabg,ggr)

return
9999 continue 
  write(ndisp,*) 'Pseudopotential data is wrong!'
  call stopp('error occurs while reading '//cfilename)
end subroutine


end module

