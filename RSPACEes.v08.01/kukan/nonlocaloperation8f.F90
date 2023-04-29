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
! **********  nonlocaloperation8f.F90 06/28/2013-01  **********

module mod_nonlocaloperation
use mod_mpi, only: nprocx,nprocy,nprocz
implicit none
contains


subroutine nonlocaloperation_r_01(natom,nprjmx,num_spe,num_ppcell,num_list,     & ! <
                                  ncpx,ncpy,ncpz,                               & ! <
                                  key_natpri_in,key_natpri_inps,                & ! <
                                  nprj,indspe,natinf,natpri,naps,lstvec2,       & ! <
                                  vnlocp,vre2,                                  & ! <
                                  rpsep)                                          ! X
implicit none
integer,intent(in)   ::natom,nprjmx,num_spe,num_ppcell,num_list,ncpx,ncpy,ncpz
integer,intent(in)   ::key_natpri_in,key_natpri_inps
integer,intent(in)   ::nprj(num_spe),indspe(natom),natinf(natom),lstvec2(num_list,num_ppcell),natpri(natom),naps(natom)
real*8, intent(in)   ::vnlocp(num_list,nprjmx,num_ppcell),vre2(ncpx,ncpy,ncpz)
real*8, intent(inout)::rpsep(nprjmx,num_ppcell)
real*8,allocatable::tmpvre(:)
integer:: na,i,ii,iaps,iatinf,iprj

  allocate(tmpvre(nprjmx))

  do na=1,natom
  if ((natpri(na) .eq. key_natpri_in) .or. (natpri(na) .eq. key_natpri_inps)) then
    iaps  = naps(na)
    iatinf= natinf(na)
    iprj  = nprj(indspe(na))

    tmpvre(1:iprj)=0.0d0

    do i=1,iprj
!$omp do
      do ii= 1,iatinf
        tmpvre(i)= tmpvre(i) +vnlocp(ii,i,iaps)*vre2(lstvec2(ii,iaps),1,1)
      end do
!$omp end do nowait
    end do

!$omp critical
    do i=1,iprj
      rpsep(i,iaps)= rpsep(i,iaps) +tmpvre(i)
    end do
!$omp end critical

  end if
  end do

  deallocate(tmpvre)
  return
end subroutine nonlocaloperation_r_01


subroutine nonlocaloperation_c_01(natom,nprjmx,num_spe,num_ppcell,num_list,neigmx,ncol,l, & ! <
                                  ncpx,ncpy,ncpz,                                         & ! <
                                  key_natpri_in,key_natpri_inps,                          & ! <
                                  dx,dy,dz,skpx,skpy,skpz,                                & ! <
                                  nprj,indspe,natinf,natpri,naps,lstvec2,                 & ! <
                                  lstx,lsty,lstz,natx,naty,natz,                          & ! <
                                  vnlocp,vcm2,                                            & ! <
                                  cpsep,                                                  & ! X
                                  vcccm)                                                    ! W
implicit none
integer,   intent(in)   ::natom,nprjmx,num_spe,num_ppcell,num_list,neigmx,ncol,l
integer,   intent(in)   ::ncpx,ncpy,ncpz
integer,   intent(in)   ::key_natpri_in,key_natpri_inps
integer,   intent(in)   ::nprj(num_spe),indspe(natom),natinf(natom),lstvec2(num_list,num_ppcell),natpri(natom),naps(natom)
integer,   intent(in)   ::lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell)
integer,   intent(in)   ::natx(natom),naty(natom),natz(natom)
real*8,    intent(in)   ::vnlocp(num_list,nprjmx,num_ppcell)
complex*16,intent(in)   ::vcm2(ncpx*ncpy*ncpz,neigmx,ncol)
real*8,    intent(in)   ::dx,dy,dz
real*8,    intent(in)   ::skpx,skpy,skpz
complex*16,intent(inout)::cpsep(nprjmx,num_ppcell,ncol)
complex*16,intent(inout)::vcccm(num_list,ncol)
complex*16,allocatable::tmpvcm(:,:)
real*8,allocatable::dskl(:),dsnl(:),dcsl(:)
integer:: na,i,j,ii,iaps,iatinf,iprj,ix0,iy0,iz0,ix,iy,iz,ns
real*8 :: xx,yy,zz

  allocate(tmpvcm(nprjmx,ncol))
  allocate(dskl(num_list),dsnl(num_list),dcsl(num_list))
  do na=1,natom
  if ((natpri(na) .eq. key_natpri_in) .or. (natpri(na) .eq. key_natpri_inps)) then
    iaps  = naps(na)
    iatinf= natinf(na)
    iprj  = nprj(indspe(na))
    ix0   = natx(na)
    iy0   = naty(na)
    iz0   = natz(na)

    tmpvcm(1:iprj,:)=dcmplx(0.0d0,0.0d0)

    do ns=1,ncol
!$omp do
      do i= 1,iatinf
! ==========  multiply phase factor to wave function  ==========
        ix= lstx(i,iaps)+ix0
        iy= lsty(i,iaps)+iy0
        iz= lstz(i,iaps)+iz0
        xx= ix*dx-0.5d0*dx
        yy= iy*dy-0.5d0*dy
        zz= iz*dz-0.5d0*dz
        dskl(i)= skpx*xx+skpy*yy+skpz*zz
! ==============================================================
      end do
!$omp do
      do i= 1,iatinf
        dcsl(i)=dcos(dskl(i))
        dsnl(i)=dsin(dskl(i))
      end do
!$omp do
      do i= 1,iatinf
        ii= lstvec2(i,iaps)
        vcccm(i,ns)= dcmplx(dcsl(i),dsnl(i))*vcm2(ii,l,ns)
      end do
    end do

    do ns= 1,ncol
      do j=1,iprj
!$omp do
        do i=1,iatinf
          tmpvcm(j,ns)= tmpvcm(j,ns)+vnlocp(i,j,iaps)*vcccm(i,ns)
        end do
!$omp end do nowait
      end do
    end do

!$omp critical
    do ns= 1,ncol
      do j= 1,iprj
        cpsep(j,iaps,ns)= cpsep(j,iaps,ns) +tmpvcm(j,ns)
      end do
    end do
!$omp end critical

  end if
  end do

  deallocate(tmpvcm)
  deallocate(dskl,dsnl,dcsl)

  return
end subroutine nonlocaloperation_c_01


subroutine nonlocaloperation_r_02(nhs,natom,num_spe,nums,nprjmx,ns,num_list,num_ppcell, & ! <
                                  ncpx,ncpy,ncpz,                                       & ! <
                                  key_natpri_in,key_natpri_inps,                        & ! <
                                  dx,dy,dz,                                             & ! <
                                  nprj,indspe,natinf,lstvec2,natpri,naps,               & ! <
                                  vnlocp,rpsep,dijsss,                                  & ! <
                                  avre,                                                 & ! X
                                  avr)                                                    ! W
use mod_stopp
implicit none
integer,intent(in)   ::nhs
integer,intent(in)   ::natom,num_spe,nums,nprjmx,ns,num_list,num_ppcell
integer,intent(in)   ::ncpx,ncpy,ncpz
integer,intent(in)   ::key_natpri_in,key_natpri_inps
real*8, intent(in)   ::dx,dy,dz
integer,intent(in)   ::nprj(num_spe),indspe(natom),natinf(natom),lstvec2(num_list,num_ppcell)
integer,intent(in)   ::natpri(natom),naps(natom)
real*8, intent(in)   ::vnlocp(num_list,nprjmx,num_ppcell),rpsep(nprjmx,num_ppcell)
real*8, intent(in)   ::dijsss(nprjmx,nprjmx,nums*nhs-nhs+1,natom*nhs+num_spe*(1-nhs))
real*8, intent(inout)::avre(ncpx,ncpy,ncpz)
real*8, intent(inout)::avr(num_list)
real*8, allocatable::tmpvre(:)
integer:: na,j,iaps,iatinf,ispe,iprj,i,jj
real*8 :: dxyzin

  if ((nhs/=1).and.(nhs/=0)) call stopp('error occurs in nonlocaloperation_r_02')
  if ((nhs==0).and.(ns/=1)) call stopp('nonlocaloperation_r_02: nhs=0 and ns/=1')

  allocate(tmpvre(nprjmx))

  dxyzin=1.0d0/(dx*dy*dz)

  do na=1,natom
  if ((natpri(na) .eq. key_natpri_in) .or. (natpri(na) .eq. key_natpri_inps)) then
    iaps  = naps(na)
    iatinf= natinf(na)
    ispe  = indspe(na)
    iprj  = nprj(indspe(na))

    if (nhs==1) then
      do i=1,iprj
        tmpvre(i)= 0.0d0
        do j= 1,iprj
          tmpvre(i)= tmpvre(i) +rpsep(j,iaps)*dijsss(j,i,ns,na)
        end do
        tmpvre(i)= tmpvre(i) *dxyzin
      end do
    else
      do i=1,iprj
        tmpvre(i)= 0.0d0
        do j= 1,iprj
          tmpvre(i)= tmpvre(i) +rpsep(j,iaps)*dijsss(j,i,ns,ispe)
        end do
        tmpvre(i)= tmpvre(i) *dxyzin
      end do
    end if

!$omp do
!ocl norecurrence(avr)
    do i=1,iatinf
      avr(i)=tmpvre(1)*vnlocp(i,1,iaps)
    end do
    do j=2,iprj
!$omp do
!ocl norecurrence(avr)
      do i=1,iatinf
        avr(i)= avr(i)+tmpvre(j)*vnlocp(i,j,iaps)
      end do
!$omp end do nowait
    end do

!$omp do
!ocl norecurrence(avre)
    do i=1,iatinf
      jj=lstvec2(i,iaps)
      avre(jj,1,1)=avre(jj,1,1)+avr(i)
    end do
!$omp end do nowait
!$omp barrier

  end if
  end do

  deallocate(tmpvre)

  return
end subroutine nonlocaloperation_r_02


subroutine nonlocaloperation_c_02(nhs,nso,natom,num_spe,neigmx,nums,nprjmx,ncol,ns1,l,num_list,num_ppcell, & ! <
                                  ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                                     & ! <
                                  key_natpri_in,key_natpri_inps,key_soc_calc,                              & ! <
                                  dx,dy,dz,skpx,skpy,skpz,                                                 & ! <
                                  nprj,indspe,natinf,lstvec2,natpri,naps,natsoc,                           & ! <
                                  lstx,lsty,lstz,natx,naty,natz,                                           & ! <
                                  vnlocp,cpsep,dijsss,dijsoc,                                              & ! <
                                  avcm,                                                                    & ! X
                                  vcccm,avc)                                                                 ! W
use mod_stopp
implicit none
integer,   intent(in)   :: nhs,nso
integer,   intent(in)   :: natom,num_spe,neigmx,nums,nprjmx,ncol,ns1,l,num_list,num_ppcell
integer,   intent(in)   :: ncpx,ncpy,ncpz,npxmax,npymax,npzmax
integer,   intent(in)   :: key_natpri_in,key_natpri_inps,key_soc_calc
real*8,    intent(in)   :: dx,dy,dz,skpx,skpy,skpz
integer,   intent(in)   :: nprj(num_spe),indspe(natom),natinf(natom),lstvec2(num_list,num_ppcell)
integer,   intent(in)   :: natpri(natom),naps(natom),natsoc(natom*nhs-nhs+1)
integer,   intent(in)   :: lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell)
integer,   intent(in)   :: natx(natom),naty(natom),natz(natom)
real*8,    intent(in)   :: vnlocp(num_list,nprjmx,num_ppcell)
complex*16,intent(in)   :: cpsep(nprjmx,num_ppcell,ncol)
real*8,    intent(in)   :: dijsss(nprjmx,nprjmx,nums*ncol*nhs-nhs+1,natom*nhs+num_spe*(1-nhs))
real*8,    intent(in)   :: dijsoc(nprjmx*nso-nso+1,nprjmx*nso-nso+1,3*nso-nso+1,natom*nso-nso+1)
complex*16,intent(inout):: avcm(ncpx*ncpy*ncpz,neigmx,ncol)
complex*16,intent(inout):: avc(num_list,ncol)
complex*16,intent(inout):: vcccm(num_list)
complex*16,allocatable::ctmpv(:,:),ctmpv2(:,:)
real*8,allocatable::dskl(:),dcsl(:),dsnl(:)
integer:: n0,na,i,j,iaps,iatinf,ispe,iprj,ix0,iy0,iz0,ix,iy,iz,ns2
real*8 :: dxyzin,xx,yy,zz
logical :: lg_overlap

  if ((max(nhs,nso)>1).or.(min(nhs,nso)<0)) call stopp('error occurs in nonlocaloperation_c_02')
  if ((nhs==0).and.(ns1/=1)) call stopp('nonlocaloperation_c_02: nhs=0 and ns1/=1')

  allocate(ctmpv(nprjmx,ncol),ctmpv2(nprjmx,ncol))
  allocate(dskl(num_list),dcsl(num_list),dsnl(num_list))

  n0= max(0,ncol-2) ! always zero, only to avoid constant array indices

  lg_overlap=((ncpx*nprocx.lt.2*npxmax).or.(ncpy*nprocy.lt.2*npymax).or.(ncpz*nprocz.lt.2*npzmax))

  dxyzin= 1.0d0/(dx*dy*dz)

  do na=1,natom
  if ((natpri(na) .eq. key_natpri_in) .or. (natpri(na) .eq. key_natpri_inps)) then
    iaps  = naps(na)
    iatinf= natinf(na)
    ispe  = indspe(na)
    iprj  = nprj(indspe(na))
    ix0   = natx(na)
    iy0   = naty(na)
    iz0   = natz(na)

    do ns2= 1,ncol
      do i= 1,iprj
        ctmpv2(i,ns2)=dcmplx(0.0d0,0.0d0)
        ctmpv(i,ns2)=cpsep(i,iaps,ns2)*dxyzin
      end do
    end do

    if (nhs==1) then
      do j= 1,iprj
        do i= 1,iprj
          if (ncol==1) then
            ctmpv2(i,1)= ctmpv2(i,1) +ctmpv(j,1)*dijsss(i,j,ns1,na)
          else
            ctmpv2(i,n0+1)= ctmpv2(i,n0+1) &
             +ctmpv(j,n0+1)*dijsss(i,j,n0+1,na) +ctmpv(j,n0+2)*dcmplx(dijsss(i,j,n0+3,na),-dijsss(i,j,n0+4,na))
            ctmpv2(i,n0+2)= ctmpv2(i,n0+2) &
             +ctmpv(j,n0+2)*dijsss(i,j,n0+2,na) +ctmpv(j,n0+1)*dcmplx(dijsss(i,j,n0+3,na),+dijsss(i,j,n0+4,na))
            if ((nso==1).and.(natsoc(na)==key_soc_calc)) then
              ctmpv2(i,n0+1)= ctmpv2(i,n0+1) &
               +ctmpv(j,n0+1)*dcmplx(0.0d0,dijsoc(i,j,n0+1,na)) -ctmpv(j,n0+2)*dcmplx(dijsoc(i,j,n0+2,na),-dijsoc(i,j,n0+3,na))
              ctmpv2(i,n0+2)= ctmpv2(i,n0+2) &
               -ctmpv(j,n0+2)*dcmplx(0.0d0,dijsoc(i,j,n0+1,na)) +ctmpv(j,n0+1)*dcmplx(dijsoc(i,j,n0+2,na),+dijsoc(i,j,n0+3,na))
            endif 
          end if
        end do
      end do
    else
      do ns2= 1,ncol
        do j= 1,iprj
          do i= 1,iprj
            ctmpv2(i,ns2)= ctmpv2(i,ns2) +ctmpv(j,ns2)*dijsss(i,j,1,ispe)
          end do
        end do
      end do
    end if

    do ns2=1,ncol
!$omp do
      do i=1,iatinf
        avc(i,ns2)=ctmpv2(1,ns2)*vnlocp(i,1,iaps)
      end do
!$omp end do nowait
      do j=2,iprj
!$omp do
!ocl norecurrence(avc)
        do i=1,iatinf
          avc(i,ns2)=avc(i,ns2)+ctmpv2(j,ns2)*vnlocp(i,j,iaps)
        end do
!$omp end do nowait
      end do
    end do

    do ns2= 1,ncol
!$omp do
      do i=1,iatinf
        ix= lstx(i,iaps)+ix0
        iy= lsty(i,iaps)+iy0
        iz= lstz(i,iaps)+iz0
        xx= ix*dx-0.5d0*dx
        yy= iy*dy-0.5d0*dy
        zz= iz*dz-0.5d0*dz
        dskl(i)= skpx*xx+skpy*yy+skpz*zz
      end do
!$omp end do nowait
!$omp do
      do i=1,iatinf
        dcsl(i)=dcos(dskl(i))
        dsnl(i)=dsin(dskl(i))
      end do
!$omp end do nowait
      if (.not. lg_overlap) then
!$omp do
!ocl norecurrence(avcm)
        do i=1,iatinf
          j=lstvec2(i,iaps)
          avcm(j,l,ns2)=avcm(j,l,ns2)+dcmplx(dcsl(i),-dsnl(i))*avc(i,ns2)
        end do
!$omp end do nowait
      else
!$omp do
        do i=1,iatinf
          vcccm(i)=dcmplx(dcsl(i),-dsnl(i))*avc(i,ns2)
        end do
!$omp single
        do i=1,iatinf
          j=lstvec2(i,iaps)
          avcm(j,l,ns2)= avcm(j,l,ns2)+vcccm(i)
        end do
!$omp end single
      end if
    end do
!$omp barrier

  end if
  end do

  deallocate(ctmpv,ctmpv2)
  deallocate(dskl,dcsl,dsnl)

  return
end subroutine nonlocaloperation_c_02


end module
