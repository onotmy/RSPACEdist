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
! **********  rw_wfpot8f.F90 12/05/2022-01  **********

module mod_rw_wfpot
contains


subroutine read_wf(  &
 nr1In,nr2In,chdir,nrc,nums,ncol,neigmx,numk,ncpx,ncpy,ncpz,ksym, & ! <
 key_ksym_inv2,                                                 & ! <
 svecre,sveccm,sval)                                              ! >
use mod_mpi
use mod_stopp
implicit none
integer,     intent(in) ::nr1In,nr2In
character(*),intent(in) ::chdir
integer,     intent(in) ::nrc,nums,ncol,neigmx,numk,ncpx,ncpy,ncpz,ksym
integer,     intent(in) ::key_ksym_inv2 
real*8,      intent(out)::svecre(ncpx*(1-nrc)+nrc,ncpy*(1-nrc)+nrc,ncpz*(1-nrc)+nrc, &
                                 neigmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk*(1-nrc)+nrc)
complex*16,  intent(out)::sveccm(ncpx*nrc-nrc+1,ncpy*nrc-nrc+1,ncpz*nrc-nrc+1, &
                                 neigmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1)
real*8,      intent(out)::sval(neigmx,nums+1-ncol,numk)
integer  ::nr1,nr2
character::fname*20
integer  ::version1,iok
integer  ::l,l1,lspin(2),ns,nk
integer  ::nrc1,nums1,ncol1,neigmx1,numk1,ncpx1,ncpy1,ncpz1,ksym1
real*8,    allocatable :: sval1(:,:,:), svecr1(:,:,:)
complex*16,allocatable :: svecc1(:,:,:)
integer,   allocatable :: maplist(:,:)

integer,parameter :: wfc_version= 5

  nr1= nr1In
  nr2= nr2In
  if (nr1<0) nr1= myr_space
  if (nr2<0) nr2= myr_kpt
  call get_fname(chdir,'wfc',nr1,nr2,fname)

  open(20,file=trim(chdir)//trim(fname),iostat=iok,form='unformatted',action='read')
  if (iok/=0) call stopp('cannot open wfc.* file(s)')
  read(20,err=999) version1
  if (version1/=wfc_version) call stopp('wrong wfc version')
  read(20,err=999) nrc1,nums1,ncol1,neigmx1,numk1,ncpx1,ncpy1,ncpz1
  read(20,err=999) ! nrx,nry,nrz,myr_kpt  not needed here
  read(20,err=999) ksym1

  if ( &
   (ncpx1/=ncpx).or.(ncpy1/=ncpy).or.(ncpz1/=ncpz) &
   .or. (nrc1>nrc).or.(nums1>nums).or.(ncol1>ncol) &
   .or. ( (numk1/=numk).and.((numk1*2/=numk).or.(ksym/=key_ksym_inv2)) ) &
   ) call stopp('read_wf: unexpected array dimensions in file '//trim(chdir)//trim(fname))
  if (neigmx1*ncol<neigmx*ncol1) call stopp('read_wf: not enough eigenstates in file '//trim(chdir)//trim(fname))

  if ((nrc1==nrc).and.(nums1==nums).and.(ncol1==ncol).and.(neigmx1==neigmx).and.(numk1==numk)) then
  ! wf arrays in code and file have the same format
    read(20,err=999) sval
    do nk=1,numk
      do ns=1,nums
        do l=1,neigmx
          if (nrc==0) read(20,err=999) svecre(:,:,:,l,ns,nk)
          if (nrc==1) read(20,err=999) sveccm(:,:,:,l,ns,nk)
        enddo
      enddo
    enddo

  else
  ! mapping of wf arrays from file to code

    allocate( sval1(neigmx1,nums1+1-ncol1,numk1) )
    if ((nrc1==0).or.(nrc==0)) allocate( svecr1(ncpx1,ncpy1,ncpz1) )
    if ((nrc1==1).or.(nrc==1)) allocate( svecc1(ncpx1,ncpy1,ncpz1) )

    read(20,err=999) sval1

    if ((ncol1==ncol).or.(nums1<nums)) then
    ! all cases apart from mapping collinear to noncollinear

      do nk=1,numk1
        do ns=1,nums1
          do l1=1,neigmx1
            if (ncol1==ncol) then
              l= l1
            else
            ! mapping nonmagnetic to noncollinear
              l= 2*l1-1
            endif

            if ((l<=neigmx).and.(ns<=nums+1-ncol)) sval(l,ns,nk)= sval1(l1,ns,nk)
            if (nrc1==0) read(20,err=999) svecr1
            if (nrc1==1) read(20,err=999) svecc1
            if (l<=neigmx) then
              if (nrc1<nrc) svecc1(:,:,:)= dcmplx(svecr1(:,:,:),0.0d0)
              if (nrc==0) svecre(:,:,:,l,ns,nk)= svecr1(:,:,:)
              if (nrc==1) sveccm(:,:,:,l,ns,nk)= svecc1(:,:,:)
              if (nums1<nums) then
                if (ncol1==ncol) then
                ! mapping nonmagnetic to collinear
                  sval(l,nums,nk)= sval(l,1,nk)
                  if (nrc==0) svecre(:,:,:,l,nums,nk)= svecre(:,:,:,l,1,nk)
                  if (nrc==1) sveccm(:,:,:,l,nums,nk)= sveccm(:,:,:,l,1,nk)
                else
                ! mapping nonmagnetic to noncollinear
                  sveccm(:,:,:,l,nums,nk)= dcmplx(0.0d0,0.0d0)
                  if (l+1<=neigmx) then
                    sval(l+1,1,nk)= sval(l,1,nk)
                    sveccm(:,:,:,l+1,nums,nk)= sveccm(:,:,:,l,1,nk)
                    sveccm(:,:,:,l+1,1   ,nk)= dcmplx(0.0d0,0.0d0)
                  endif
                endif
              endif
            endif ! (l<=neigmx)

          enddo ! l1
        enddo ! ns
      enddo ! nk

    else
    ! mapping magnetic collinear to noncollinear

      allocate( maplist(nums1,neigmx1) )
      do nk=1,numk1
        if (myr_space==0) then
          maplist(:,:)= 0
          l= 1
          lspin(:)= 1
          do while (l<=neigmx)
            ns= 0
            if (lspin(1)>neigmx1) ns= 2
            if (lspin(2)>neigmx1) ns= 1
            if (ns==0) then
              if ( sval1(lspin(1),1,nk) <= sval1(lspin(2),2,nk) ) then
                ns= 1
              else
                ns= 2
              endif
            endif
            maplist(ns,lspin(ns))= l
            lspin(ns)= lspin(ns) + 1
            l= l+1
          enddo
        endif
        call mpi_bcast(maplist,nums1*neigmx1,mpi_integer,0,mpicom_space,mpij)
        do ns=1,nums1
          do l1=1,neigmx1
            l= maplist(ns,l1)
            if (nrc1==0) then
              read(20,err=999) svecr1
              svecc1(:,:,:)= dcmplx(svecr1(:,:,:),0.0d0)
            else
              read(20,err=999) svecc1
            endif
            if (l>0) then
              sval(l,1,nk)= sval1(l1,ns,nk)
              sveccm(:,:,:,l,ns,nk)= svecc1(:,:,:)
              sveccm(:,:,:,l,3-ns,nk)= dcmplx(0.0d0,0.0d0)
            endif
          enddo
        enddo

      enddo ! nk
      deallocate( maplist )

    endif ! ((ncol1==ncol).or.(nums1<nums)) else

    deallocate( sval1 )
    if ((nrc1==0).or.(nrc==0)) deallocate( svecr1 )
    if ((nrc1==1).or.(nrc==1)) deallocate( svecc1 )

    if (numk1*2==numk) then
    ! mapping of k-points (assuming appropriate k-point set to exploit time inversion symmetry)
      if (ksym/=key_ksym_inv2) call stopp('read_wf: wrong ksym for k-point mapping')
      do nk= 1,numk1
        sval(:,:,2*numk1+1-nk)= sval(:,:,nk)
        if (nrc==0) svecre(:,:,:,:,:,2*numk1+1-nk)= svecre(:,:,:,:,:,nk)
        if (nrc==1) sveccm(:,:,:,:,:,2*numk1+1-nk)= dconjg(sveccm(:,:,:,:,:,nk))
      enddo
    endif

  endif ! mapping of wf arrays from file to code
  close(20)

  return
  999 call stopp('error found in wfc.*')
end subroutine read_wf


subroutine read_rho( &
 nrIn,chdir,nevhist,nspv,ncpx,ncpy,ncpz,nradmx,npoint,num_atcell, & ! <
 lbrydn,na_atcell,natcell_inf,rhosmt,rhotrur,rhosmtr)               ! >
use mod_mpi
use mod_stopp
implicit none
integer,     intent(in) ::nrIn
character(*),intent(in) ::chdir
integer,     intent(in) ::nevhist,nspv,ncpx,ncpy,ncpz,nradmx,npoint,num_atcell
integer,     intent(out)::lbrydn,na_atcell,natcell_inf(num_atcell)
real*8,      intent(out)::rhosmt(ncpx,ncpy,ncpz,nspv)
real*8,      intent(out)::rhotrur(nradmx,npoint,nspv,num_atcell),rhosmtr(nradmx,npoint,nspv,num_atcell)
integer  ::nr
character::fname*20
logical  ::histyn 
integer  ::version1,iok
integer  ::nspv1,ncpx1,ncpy1,ncpz1,nradmx1,npoint1,num_atcell1
integer  ::ns,ns1,ix,iy,iz,ir,il,na
real*8   ::fctr
integer,allocatable:: natcell_inf1(:)
real*8, allocatable:: rhosmt1(:,:,:,:),rhotrur1(:,:,:,:),rhosmtr1(:,:,:,:)

integer,parameter :: rho_version= 3

  histyn= ((nevhist==1).or.(nevhist==3).or.(nevhist==-1).or.(nevhist==4))
  if (histyn) then
    nr= nrIn
    if (nr<0) nr=myr_space
    if (nevhist==4) then 
      call get_fname(chdir,'r_o',nr,-1,fname)
    else 
      call get_fname(chdir,'rho',nr,-1,fname)
    endif 
    if (myrank_glbl==0) inquire(file=trim(chdir)//trim(fname),exist=histyn)
    call mpi_bcast(histyn,1,mpi_logical,0,mpicom_space,mpij)
    if ((nevhist==-1) .and. (.not. histyn)) call stopp('read_rho: nevhist=-1, but no density file rho.*')
    if (histyn) then

      lbrydn= 1
      open(20,file=trim(chdir)//trim(fname),iostat=iok,form='unformatted',action='read')
      if (iok/=0) call stopp('cannot open rho.* file(s)')
      read(20,err=999) version1
      if (version1/=rho_version) call stopp('wrong rho version')
      read(20,err=999) nspv1,ncpx1,ncpy1,ncpz1,nradmx1,npoint1,num_atcell1
      read(20,err=999) ! nrx,nry,nrz  not needed here
      allocate(natcell_inf1(num_atcell1))
      read(20,err=999) na_atcell,natcell_inf1
      if ((nspv1>nspv).and.(nspv1==2)) &
       call stopp('read_rho:  mapping of (mag dens) to (non-mag dens) not implemented')
      if ((nspv1>nspv).and.(nspv1==4)) &
       call stopp('read_rho:  mapping of (noco dens) to (col dens) not implemented')
      if ( (nspv1>nspv).or.(ncpx1/=ncpx).or.(ncpy1/=ncpy).or.(ncpz1/=ncpz).or. &
       (nradmx1/=nradmx).or.(npoint1/=npoint)) &
       call stopp('read_rho: unexpected array dimensions in density file '//trim(chdir)//trim(fname))
      if (na_atcell>num_atcell) &
       call stopp('read_rho:  num_atcell is too small')
      natcell_inf(1:na_atcell)= natcell_inf1(1:na_atcell)
      deallocate(natcell_inf1)
      if (nspv1==nspv) then
        read(20,err=999) rhosmt
        if (num_atcell1<=num_atcell) then
          read(20,err=999) rhotrur(:,:,:,1:num_atcell1)
          read(20,err=999) rhosmtr(:,:,:,1:num_atcell1)
        else
          allocate(rhotrur1(nradmx,npoint,nspv,num_atcell1))
          read(20,err=999) rhotrur1
          rhotrur(:,:,:,1:na_atcell)= rhotrur1(:,:,:,1:na_atcell) 
          read(20,err=999) rhotrur1
          rhosmtr(:,:,:,1:na_atcell)= rhotrur1(:,:,:,1:na_atcell) 
          deallocate(rhotrur1)
        endif
      else
        allocate( &
         rhosmt1(ncpx,ncpy,ncpz,nspv1), &
         rhotrur1(nradmx,npoint,nspv1,num_atcell1), rhosmtr1(nradmx,npoint,nspv1,num_atcell1) )
        read(20,err=999) rhosmt1
        read(20,err=999) rhotrur1
        read(20,err=999) rhosmtr1
        do ns= 1,nspv
          ns1= min(ns,nspv1)
          fctr=1.0d0
          if (nspv1==1) fctr=0.5d0
          if (ns>2)     fctr=0.0d0
          do iz=1,ncpz
          do iy=1,ncpy
          do ix=1,ncpx
            rhosmt(ix,iy,iz,ns)= rhosmt1(ix,iy,iz,ns1)*fctr
          enddo
          enddo
          enddo
          do na=1,na_atcell
            do il=1,npoint
              do ir=1,nradmx
                rhotrur(ir,il,ns,na)= rhotrur1(ir,il,ns1,na)*fctr
                rhosmtr(ir,il,ns,na)= rhosmtr1(ir,il,ns1,na)*fctr
              enddo
            enddo
          enddo
        enddo
        deallocate( rhosmt1,rhotrur1,rhosmtr1 )
      endif
      close(20)

    endif
  endif
  if (.not. histyn) then
    lbrydn= 0 
    ! rhosmt(:,:,:,:) = 0.0d0
    ! rhotrur(:,:,:,:)= 0.0d0
    ! rhosmtr(:,:,:,:)= 0.0d0
  endif

  return
  999 call stopp('error found in rho.*')
end subroutine read_rho


subroutine read_vht(  &
 nrIn,chdir,nevhist,ncpx,ncpy,ncpz,nmesh, & ! <
 vh_dense)                                  ! >
use mod_mpi
use mod_stopp
implicit none
integer,     intent(in) ::nrIn
character(*),intent(in) ::chdir
integer,     intent(in) ::nevhist,ncpx,ncpy,ncpz,nmesh
real*8,      intent(out)::vh_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
integer  ::nr
character::fname*20
logical  ::histyn
integer  ::version1,iok
integer  ::ncpx1,ncpy1,ncpz1,nmesh1

integer,parameter :: vht_version= 3

  histyn= ((nevhist==2).or.(nevhist==3).or.(nevhist==-1))
  if (histyn) then
    nr= nrIn
    if (nr<0) nr=myr_space
    call get_fname(chdir,'vht',nr,-1,fname)
    if (myrank_glbl==0) inquire(file=trim(chdir)//trim(fname),exist=histyn)
    call mpi_bcast(histyn,1,mpi_logical,0,mpicom_space,mpij)
    if (histyn) then

      open(20,file=trim(chdir)//trim(fname),iostat=iok,form='unformatted',action='read')
      if (iok/=0) call stopp('cannot open vht.* file(s)')
      read(20,err=999) version1
      if (version1/=vht_version) call stopp('wrong vht version')
      read(20,err=999) ncpx1,ncpy1,ncpz1,nmesh1
      read(20,err=999) ! nrx,nry,nrz  not needed here
      if ((ncpx1/=ncpx).or.(ncpy1/=ncpy).or.(ncpz1/=ncpz).or.(nmesh1/=nmesh)) &
       call stopp('read_vht: unexpected array dimensions in potential file '//trim(chdir)//trim(fname))
      read(20,err=999) vh_dense
      close(20)

    endif
  endif
  if (.not. histyn) vh_dense(:,:,:)= 0.0d0

  return
  999 call stopp('error found in vht.*')
end subroutine read_vht 


subroutine write_wf(  &
 nr1In,nr2In,chdir,nrc,nums,ncol,neigmx,numk,ncpx,ncpy,ncpz,nrxIn,nryIn,nrzIn,ksym, & ! <
 svecre,sveccm,sval)                                                                  ! <
use mod_mpi,only: myr_space,myr_kpt,myrx,myry,myrz
use mod_stopp
implicit none
integer,     intent(in)::nr1In,nr2In
character(*),intent(in)::chdir
integer,     intent(in)::nrc,nums,ncol,neigmx,numk,ncpx,ncpy,ncpz,nrxIn,nryIn,nrzIn,ksym
real*8,      intent(in)::svecre(ncpx*(1-nrc)+nrc,ncpy*(1-nrc)+nrc,ncpz*(1-nrc)+nrc, &
                                 neigmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk*(1-nrc)+nrc)
complex*16,  intent(in)::sveccm(ncpx*nrc-nrc+1,ncpy*nrc-nrc+1,ncpz*nrc-nrc+1, &
                                 neigmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1)
real*8,      intent(in)::sval(neigmx,nums+1-ncol,numk)
integer  ::nr1,nr2,nrx,nry,nrz, l,ns,nk
character::fname*20

integer,parameter :: wfc_version= 5

  nr1= nr1In
  nr2= nr2In
  if (nr1<0) nr1= myr_space
  if (nr2<0) nr2= myr_kpt
  nrx= nrxIn
  nry= nryIn
  nrz= nrzIn
  if (nrx<0) nrx= myrx
  if (nry<0) nry= myry
  if (nrz<0) nrz= myrz

  call get_fname(chdir,'wfc',nr1,nr2,fname)
  open(20,file=trim(chdir)//trim(fname),form='unformatted',status='replace')
  write(20) wfc_version
  write(20) nrc,nums,ncol,neigmx,numk,ncpx,ncpy,ncpz
  write(20) nrx,nry,nrz,myr_kpt
  write(20) ksym
  write(20) sval
  do nk=1,numk
    do ns=1,nums
      do l=1,neigmx
        if (nrc==0) write(20) svecre(:,:,:,l,ns,nk)
        if (nrc==1) write(20) sveccm(:,:,:,l,ns,nk)
      enddo
    enddo
  enddo
  close(20)

end subroutine write_wf


subroutine write_rhovht(  &
 nrIn,chdir,nrho,nele,npre,nevhist,nspv,ncpx,ncpy,ncpz,nmesh,nradmx,npoint,num_atcell,nrxIn,nryIn,nrzIn, & ! <
 natom,key_natpri_in,natpri,                                                                             & ! <
 rhosmt_i,rhotrur_i,rhosmtr_i,rhosmt_o,rhotrur_o,rhosmtr_o,vh_dense)                                       ! <
use mod_mpi,only: myr_space,myrx,myry,myrz
use mod_stopp
implicit none
integer,     intent(in)::nrIn
character(*),intent(in)::chdir
integer,     intent(in)::nrho,nele,npre,nevhist
integer,     intent(in)::nspv,ncpx,ncpy,ncpz,nmesh,nradmx,npoint,num_atcell,nrxIn,nryIn,nrzIn
integer,     intent(in)::natom,key_natpri_in,natpri(natom)
real*8,      intent(in)::rhosmt_i(ncpx,ncpy,ncpz,nspv)
real*8,      intent(in)::rhotrur_i(nradmx,npoint,nspv,num_atcell),rhosmtr_i(nradmx,npoint,nspv,num_atcell)
real*8,      intent(in)::rhosmt_o(ncpx,ncpy,ncpz,nspv)
real*8,      intent(in)::rhotrur_o(nradmx,npoint,nspv,num_atcell),rhosmtr_o(nradmx,npoint,nspv,num_atcell)
real*8,      intent(in)::vh_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
character::fname*20
integer  ::nr,nrx,nry,nrz,irho 
integer  ::na,na_atcell,num_atcell1
integer,allocatable:: natcell_inf(:)

logical,parameter :: l_samedim= .false.
integer,parameter :: rho_version= 3
integer,parameter :: vht_version= 3

  nr = nrIn
  if (nr <0) nr = myr_space
  nrx= nrxIn
  nry= nryIn
  nrz= nrzIn
  if (nrx<0) nrx= myrx
  if (nry<0) nry= myry
  if (nrz<0) nrz= myrz

  if ((nevhist==1).or.(nevhist==3)) then
    do irho= 1,nrho 
      if (irho==1) then
        call get_fname(chdir,'rho',nr,-1,fname)
      else
        call get_fname(chdir,'r_o',nr,-1,fname)
      endif 
      open(20,file=trim(chdir)//trim(fname),form='unformatted',status='replace')
      if (npre==2) then
        close(20,status='delete')
      else
        allocate( natcell_inf(num_atcell) )
        natcell_inf(:)= 0
        na_atcell= 0
        do na= 1,natom
          if (natpri(na)==key_natpri_in) then
            na_atcell= na_atcell +1
            if (na_atcell>num_atcell) call stopp('write_rhovht: na_atcell > num_atcell')
            natcell_inf(na_atcell)= na
          endif
        enddo
        if (l_samedim) then
          num_atcell1= num_atcell
        else
          num_atcell1= na_atcell
        endif
        write(20) rho_version
        write(20) nspv,ncpx,ncpy,ncpz,nradmx,npoint,num_atcell1
        write(20) nrx,nry,nrz
        write(20) na_atcell,natcell_inf
        if (irho==1) then  
          write(20) rhosmt_i
          write(20) rhotrur_i(:,:,:,1:num_atcell1)
          write(20) rhosmtr_i(:,:,:,1:num_atcell1)
        else
          write(20) rhosmt_o
          write(20) rhotrur_o(:,:,:,1:num_atcell1)
          write(20) rhosmtr_o(:,:,:,1:num_atcell1)
        endif 
        deallocate( natcell_inf )
        close(20)
      endif
    enddo 
  endif

  if ((nevhist==2).or.(nevhist==3)) then
    call get_fname(chdir,'vht',nr,-1,fname)
    open(20,file=trim(chdir)//trim(fname),form='unformatted',status='replace')
    if (npre==2) then
      close(20,status='delete')
    else
      write(20) vht_version
      write(20) ncpx,ncpy,ncpz,nmesh
      write(20) nrx,nry,nrz
      write(20) vh_dense
      close(20)
    endif
  endif

  if (nele==1) then
    call get_fname(chdir,'ele',nr,-1,fname)
    open(20,file=trim(chdir)//trim(fname),form='unformatted',status='replace')
    write(20) rhosmt_o 
    close(20)
  endif

end subroutine write_rhovht


subroutine filenumbers_rhovht( &
 chdir,       & ! <
 nrrho,nrpot)   ! >
use mod_stopp
implicit none
character(*),intent(in) :: chdir
integer,     intent(out):: nrrho,nrpot
character:: fname*20
logical  :: fileyn

  nrrho=0
  fileyn= .true.
  do while (fileyn)
    call get_fname(chdir,'rho',nrrho,-1,fname)
    inquire(file=trim(chdir)//trim(fname),exist=fileyn)
    if (fileyn) nrrho= nrrho+1
  enddo

  nrpot=0
  fileyn= .true.
  do while (fileyn)
    call get_fname(chdir,'vht',nrpot,-1,fname)
    inquire(file=trim(chdir)//trim(fname),exist=fileyn)
    if (fileyn) nrpot= nrpot+1
  enddo

end subroutine filenumbers_rhovht


subroutine get_fname( &
 chdir,chnam,nr1,nr2, & ! <
 fname)                 ! >
use mod_stopp
implicit none
character(*),intent(in) ::chdir,chnam
integer,     intent(in) ::nr1,nr2
character(*),intent(out)::fname
integer  :: flen
character:: chfmt*30

integer,parameter :: digit1=5, digit2=3

  if (nr2>=0) then
    chfmt= '(    a,".",iX.X,".",iX.X)'
    write(chfmt(13:15),fmt='(i1,".",i1)') digit1,digit1
    write(chfmt(22:24),fmt='(i1,".",i1)') digit2,digit2
    flen= len_trim(chnam) +2 +digit1+digit2
  else
    chfmt= '(    a,".",iX.X)'
    write(chfmt(13:15),fmt='(i1,".",i1)') digit1,digit1
    flen= len_trim(chnam) +1 +digit1
  endif
  if (len_trim(chdir)>0) then
    chfmt(2:5)= '"/",'
    flen= flen +1
  endif

  if (len(fname)<flen) call stopp('get_fname: string too short')
  if ((nr1>=10**digit1).or.(nr2>=10**digit2)) call stopp('get_fname: not enough digits')

  if (nr2>=0) then
    write(fname,fmt=trim(chfmt)) trim(chnam),nr1,nr2
  else
    write(fname,fmt=trim(chfmt)) trim(chnam),nr1
  endif

end subroutine get_fname


subroutine get_wfpotdir( &
 ndisp,                     & ! <
 scfctr2,chdirinp,chdirout)   ! >
use mod_mpi
implicit none
integer,     intent(in) ::ndisp
logical,     intent(out)::scfctr2
character(*),intent(out)::chdirinp,chdirout
logical:: fileyn
integer:: ifileend
character(len(chdirinp)+10):: chdir

  if (myrank_glbl==0) then
    scfctr2 = .true.  
    inquire(file='wfcpot.dir',exist=fileyn)
    if (fileyn) then
      open(20,file='wfcpot.dir',form='formatted',action='read')
      ifileend= 0
      do while (ifileend==0)
        read(20,fmt='(a)',iostat=ifileend) chdir
        if (ifileend==0) then
          if (index(chdir,'inp:')==1) then
            chdir(1:4)= '    '
            chdirinp= adjustl(chdir)
          endif
          if (index(chdir,'out:')==1) then
            chdir(1:4)= '    '
            chdirout= adjustl(chdir)
          endif
        endif
      enddo
      close(20)
      if (len_trim(chdirinp)/=0) write(ndisp,fmt='(1x,a,1x,a)') 'input directory: ',trim(chdirinp)
      if (trim(chdirout)=='no output') then
        scfctr2= .false.
        chdirout= ''
        write(ndisp,fmt='(1x,a)') 'no output of wfc.* rho.* vht.*'
      else 
        if (len_trim(chdirout)/=0) write(ndisp,fmt='(1x,a,1x,a)') 'output directory:',trim(chdirout)
      endif
    endif
  endif
  call mpi_bcast(scfctr2,1,mpi_logical,0,mpi_comm_world,mpij)
  call mpi_bcast(chdirinp,len(chdirinp),mpi_character,0,mpi_comm_world,mpij)
  call mpi_bcast(chdirout,len(chdirout),mpi_character,0,mpi_comm_world,mpij)

end subroutine get_wfpotdir


end module

