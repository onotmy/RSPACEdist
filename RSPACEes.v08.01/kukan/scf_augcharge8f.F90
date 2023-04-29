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
! **********  scf_augcharge8e.f90 01/05/2020-01  **********

module mod_scf_augcharge
implicit none

contains

subroutine scf_augcharge( &
 nrd,key_natpri_in,key_natpri_inps,key_pp_paw,                     & ! <
 nspv,nradmx,npoint,lrhomx,lmx,natom,num_atcell,num_spe,           & ! <
 num_list_d,num_ppcell_d,                                          & ! <
 indspe,natpri,natprid,natpri_inf,napsd,natinfd,ndatx,ndaty,ndatz, & ! <
 ntyppp,nradct,lpmx,nwexp,lstdx,lstdy,lstdz,                       & ! <
 ddx,ddy,ddz,yylm,wt,radial,dradial,rfac,                          & ! <
 atx,aty,atz,rhotrur,rhosmtr,                                      & ! <
 rhoaugr,rhoaug3d, drhoaug3ddx,drhoaug3ddy,drhoaug3ddz)              ! >
use mod_mpi
use mod_stopp
implicit none
integer,intent(in) :: nrd
integer,intent(in) :: key_natpri_in,key_natpri_inps,key_pp_paw
integer,intent(in) :: nspv,nradmx,npoint,lrhomx,lmx,natom,num_atcell,num_spe
integer,intent(in) :: num_list_d,num_ppcell_d
integer,intent(in) :: indspe(natom),natpri(natom),natprid(natom),natpri_inf(natom),napsd(natom),natinfd(natom)
integer,intent(in) :: ndatx(natom),ndaty(natom),ndatz(natom)
integer,intent(in) :: ntyppp(num_spe),nradct(num_spe),lpmx(num_spe),nwexp(num_spe)
integer,intent(in) :: lstdx(num_list_d,num_ppcell_d),lstdy(num_list_d,num_ppcell_d),lstdz(num_list_d,num_ppcell_d)
real*8, intent(in) :: ddx,ddy,ddz
real*8, intent(in) :: yylm(npoint,lrhomx),wt(npoint)
real*8, intent(in) :: radial(nradmx,num_spe),dradial(nradmx,num_spe),rfac(0:2*(lmx-1),num_spe)
real*8, intent(in) :: atx(natom),aty(natom),atz(natom)
real*8, intent(in) :: rhotrur(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) :: rhosmtr(nradmx,npoint,nspv,num_atcell)
real*8, intent(out):: rhoaugr(nradmx,npoint,num_atcell)
real*8, intent(out):: rhoaug3d(num_list_d,num_ppcell_d)
real*8, intent(out):: drhoaug3ddx(num_list_d*nrd-nrd+1,num_ppcell_d*nrd-nrd+1)
real*8, intent(out):: drhoaug3ddy(num_list_d*nrd-nrd+1,num_ppcell_d*nrd-nrd+1)
real*8, intent(out):: drhoaug3ddz(num_list_d*nrd-nrd+1,num_ppcell_d*nrd-nrd+1)
real*8,allocatable::rhoaugmom(:,:),rhoaugmoma(:,:)
  allocate(rhoaugmom(lrhomx,natom),rhoaugmoma(lrhomx,natom))

  if (lrhomx > 25) call stopp ('error in augcharge. lrhomx must be <= 25.')
  if ((nrd/=0).and.(nrd/=1)) call stopp ('error in augcharge. nrd should be 0 or 1.')

!$omp parallel default(shared)

  call augcharge_01( &
   key_natpri_in,key_pp_paw,                           & ! <
   nspv,nradmx,npoint,lrhomx,natom,num_atcell,num_spe, & ! <
   indspe,natpri,natpri_inf,ntyppp,nradct,lpmx,        & ! <
   yylm,wt,radial,dradial,                             & ! <
   rhotrur,rhosmtr,                                    & ! <
   rhoaugmom)                                            ! >

!$omp single
  call mpi_allreduce(rhoaugmom,rhoaugmoma,lrhomx*natom,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  rhoaugmom=rhoaugmoma
!$omp end single

  call augcharge_02( &
   key_natpri_in,key_pp_paw,                          & ! <
   nradmx,npoint,lrhomx,lmx,natom,num_atcell,num_spe, & ! <
   indspe,natpri,natpri_inf,ntyppp,nradct,lpmx,nwexp, & ! <
   yylm,radial,rfac,                                  & ! <
   rhoaugmom,                                         & ! <
   rhoaugr)                                             ! >
!$omp barrier

  if (nrd==1) then
    call augcharge_03( &
     key_natpri_inps,key_pp_paw,                                                                & ! <
     nradmx,lrhomx,lmx,natom,num_spe,num_list_d,num_ppcell_d,                                   & ! <
     indspe,natprid,napsd,natinfd,ndatx,ndaty,ndatz,ntyppp,nradct,lpmx,nwexp,lstdx,lstdy,lstdz, & ! <
     ddx,ddy,ddz,radial,rfac,atx,aty,atz,rhoaugmom,                                             & ! <
     rhoaug3d,drhoaug3ddx,drhoaug3ddy,drhoaug3ddz)                                                ! >
  end if

  if (nrd==0) then
    call augcharge_04( &
     key_natpri_inps,key_pp_paw,                                                                & ! <
     nradmx,lrhomx,lmx,natom,num_spe,num_list_d,num_ppcell_d,                                   & ! <
     indspe,natprid,napsd,natinfd,ndatx,ndaty,ndatz,ntyppp,nradct,lpmx,nwexp,lstdx,lstdy,lstdz, & ! <
     ddx,ddy,ddz,radial,rfac,atx,aty,atz,rhoaugmom,                                             & ! <
     rhoaug3d)                                                                                    ! >
  end if

!$omp end parallel

  deallocate(rhoaugmom,rhoaugmoma)
  return
end subroutine


subroutine augcharge_01( &
 key_natpri_in,key_pp_paw,                           & ! <
 nspv,nradmx,npoint,lrhomx,natom,num_atcell,num_spe, & ! <
 indspe,natpri,natpri_inf,ntyppp,nradct,lpmx,        & ! <
 yylm,wt,radial,dradial,                             & ! <
 rhotrur,rhosmtr,                                    & ! <
 rhoaugmom)                                            ! >
use mod_mpi
implicit none
integer,intent(in)   :: key_natpri_in,key_pp_paw
integer,intent(in)   :: nspv,nradmx,npoint,lrhomx,natom,num_atcell,num_spe
integer,intent(in)   :: indspe(natom),natpri(natom),natpri_inf(natom),nradct(num_spe),lpmx(num_spe),ntyppp(num_spe)
real*8, intent(in)   :: yylm(npoint,lrhomx),wt(npoint)
real*8, intent(in)   :: radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8, intent(in)   :: rhotrur(nradmx,npoint,nspv,num_atcell),rhosmtr(nradmx,npoint,nspv,num_atcell)
real*8, intent(out)  :: rhoaugmom(lrhomx,natom)
integer na,ipri,ispe,il,ir,i,i_end
real*8, allocatable::sum(:),tmp(:)
real*8  r,r2,r3,r4,r5,r6,tmp0,tmpr2,tmpr3,tmpr4,tmpr5,tmpr6

  allocate(sum(25),tmp(25))

!$omp do
  do i=1,lrhomx*natom
    rhoaugmom(i,1)=0.0d0
  end do

  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
      ipri=natpri_inf(na)
      ispe=indspe(na)
      if (lpmx(ispe)==0) i_end=1
      if (lpmx(ispe)==1) i_end=9
      if (lpmx(ispe)==2) i_end=25
      do i=1,i_end
        sum(i)=0.0d0
      end do
!$omp do
      do il=1,npoint
        do i=1,i_end
          tmp(i)=yylm(il,i)*wt(il)
        end do
        if (nspv>1) then
          select case (lpmx(ispe))
          case (0)
            do ir=2,nradct(ispe)-1
              tmp0=(rhotrur(ir,il,1,ipri)-rhosmtr(ir,il,1,ipri)+rhotrur(ir,il,2,ipri)-rhosmtr(ir,il,2,ipri))*dradial(ir,ispe)
              r=radial(ir,ispe)
              r2=r*r
              tmpr2=tmp0*r2
              sum(1)=sum(1)+tmpr2*tmp(1)
            end do
          case (1)
            do ir=2,nradct(ispe)-1
              tmp0=(rhotrur(ir,il,1,ipri)-rhosmtr(ir,il,1,ipri)+rhotrur(ir,il,2,ipri)-rhosmtr(ir,il,2,ipri))*dradial(ir,ispe)
              r=radial(ir,ispe)
              r2=r*r
              r3=r2*r
              r4=r3*r
              tmpr2=tmp0*r2
              tmpr3=tmp0*r3
              tmpr4=tmp0*r4
              sum(1)=sum(1)+tmpr2*tmp(1)
              sum(2)=sum(2)+tmpr3*tmp(2)
              sum(3)=sum(3)+tmpr3*tmp(3)
              sum(4)=sum(4)+tmpr3*tmp(4)
              sum(5)=sum(5)+tmpr4*tmp(5)
              sum(6)=sum(6)+tmpr4*tmp(6)
              sum(7)=sum(7)+tmpr4*tmp(7)
              sum(8)=sum(8)+tmpr4*tmp(8)
              sum(9)=sum(9)+tmpr4*tmp(9)
            end do
          case(2)
            do ir=2,nradct(ispe)-1
              tmp0=(rhotrur(ir,il,1,ipri)-rhosmtr(ir,il,1,ipri)+rhotrur(ir,il,2,ipri)-rhosmtr(ir,il,2,ipri))*dradial(ir,ispe)
              r=radial(ir,ispe)
              r2=r*r
              r3=r2*r
              r4=r3*r
              r5=r4*r
              r6=r5*r
              tmpr2=tmp0*r2
              tmpr3=tmp0*r3
              tmpr4=tmp0*r4
              tmpr5=tmp0*r5
              tmpr6=tmp0*r6
              sum( 1)=sum( 1)+tmpr2*tmp( 1)
              sum( 2)=sum( 2)+tmpr3*tmp( 2)
              sum( 3)=sum( 3)+tmpr3*tmp( 3)
              sum( 4)=sum( 4)+tmpr3*tmp( 4)
              sum( 5)=sum( 5)+tmpr4*tmp( 5)
              sum( 6)=sum( 6)+tmpr4*tmp( 6)
              sum( 7)=sum( 7)+tmpr4*tmp( 7)
              sum( 8)=sum( 8)+tmpr4*tmp( 8)
              sum( 9)=sum( 9)+tmpr4*tmp( 9)
              sum(10)=sum(10)+tmpr5*tmp(10)
              sum(11)=sum(11)+tmpr5*tmp(11)
              sum(12)=sum(12)+tmpr5*tmp(12)
              sum(13)=sum(13)+tmpr5*tmp(13)
              sum(14)=sum(14)+tmpr5*tmp(14)
              sum(15)=sum(15)+tmpr5*tmp(15)
              sum(16)=sum(16)+tmpr5*tmp(16)
              sum(17)=sum(17)+tmpr6*tmp(17)
              sum(18)=sum(18)+tmpr6*tmp(18)
              sum(19)=sum(19)+tmpr6*tmp(19)
              sum(20)=sum(20)+tmpr6*tmp(20)
              sum(21)=sum(21)+tmpr6*tmp(21)
              sum(22)=sum(22)+tmpr6*tmp(22)
              sum(23)=sum(23)+tmpr6*tmp(23)
              sum(24)=sum(24)+tmpr6*tmp(24)
              sum(25)=sum(25)+tmpr6*tmp(25)
            end do
          end select
        else
          select case (lpmx(ispe))
          case (0)
            do ir=2,nradct(ispe)-1
              tmp0=(rhotrur(ir,il,1,ipri)-rhosmtr(ir,il,1,ipri))*dradial(ir,ispe)
              r=radial(ir,ispe)
              r2=r*r
              tmpr2=tmp0*r2
              sum(1)=sum(1)+tmpr2*tmp(1)
            end do
          case(1)
            do ir=2,nradct(ispe)-1
              tmp0=(rhotrur(ir,il,1,ipri)-rhosmtr(ir,il,1,ipri))*dradial(ir,ispe)
              r=radial(ir,ispe)
              r2=r*r
              r3=r2*r
              r4=r3*r
              tmpr2=tmp0*r2
              tmpr3=tmp0*r3
              tmpr4=tmp0*r4
              sum(1)=sum(1)+tmpr2*tmp(1)
              sum(2)=sum(2)+tmpr3*tmp(2)
              sum(3)=sum(3)+tmpr3*tmp(3)
              sum(4)=sum(4)+tmpr3*tmp(4)
              sum(5)=sum(5)+tmpr4*tmp(5)
              sum(6)=sum(6)+tmpr4*tmp(6)
              sum(7)=sum(7)+tmpr4*tmp(7)
              sum(8)=sum(8)+tmpr4*tmp(8)
              sum(9)=sum(9)+tmpr4*tmp(9)
            end do
          case(2)
            do ir=2,nradct(ispe)-1
              tmp0=(rhotrur(ir,il,1,ipri)-rhosmtr(ir,il,1,ipri))*dradial(ir,ispe)
              r=radial(ir,ispe)
              r2=r*r
              r3=r2*r
              r4=r3*r
              r5=r4*r
              r6=r5*r
              tmpr2=tmp0*r2
              tmpr3=tmp0*r3
              tmpr4=tmp0*r4
              tmpr5=tmp0*r5
              tmpr6=tmp0*r6
              sum( 1)=sum( 1)+tmpr2*tmp( 1)
              sum( 2)=sum( 2)+tmpr3*tmp( 2)
              sum( 3)=sum( 3)+tmpr3*tmp( 3)
              sum( 4)=sum( 4)+tmpr3*tmp( 4)
              sum( 5)=sum( 5)+tmpr4*tmp( 5)
              sum( 6)=sum( 6)+tmpr4*tmp( 6)
              sum( 7)=sum( 7)+tmpr4*tmp( 7)
              sum( 8)=sum( 8)+tmpr4*tmp( 8)
              sum( 9)=sum( 9)+tmpr4*tmp( 9)
              sum(10)=sum(10)+tmpr5*tmp(10)
              sum(11)=sum(11)+tmpr5*tmp(11)
              sum(12)=sum(12)+tmpr5*tmp(12)
              sum(13)=sum(13)+tmpr5*tmp(13)
              sum(14)=sum(14)+tmpr5*tmp(14)
              sum(15)=sum(15)+tmpr5*tmp(15)
              sum(16)=sum(16)+tmpr5*tmp(16)
              sum(17)=sum(17)+tmpr6*tmp(17)
              sum(18)=sum(18)+tmpr6*tmp(18)
              sum(19)=sum(19)+tmpr6*tmp(19)
              sum(20)=sum(20)+tmpr6*tmp(20)
              sum(21)=sum(21)+tmpr6*tmp(21)
              sum(22)=sum(22)+tmpr6*tmp(22)
              sum(23)=sum(23)+tmpr6*tmp(23)
              sum(24)=sum(24)+tmpr6*tmp(24)
              sum(25)=sum(25)+tmpr6*tmp(25)
            end do
          end select
        end if
      end do
!$omp critical
      rhoaugmom( 1,na)=rhoaugmom( 1,na)+sum( 1)
      if (lpmx(ispe)>0) then
        rhoaugmom( 2,na)=rhoaugmom( 2,na)+sum( 2)
        rhoaugmom( 3,na)=rhoaugmom( 3,na)+sum( 3)
        rhoaugmom( 4,na)=rhoaugmom( 4,na)+sum( 4)
        rhoaugmom( 5,na)=rhoaugmom( 5,na)+sum( 5)
        rhoaugmom( 6,na)=rhoaugmom( 6,na)+sum( 6)
        rhoaugmom( 7,na)=rhoaugmom( 7,na)+sum( 7)
        rhoaugmom( 8,na)=rhoaugmom( 8,na)+sum( 8)
        rhoaugmom( 9,na)=rhoaugmom( 9,na)+sum( 9)
      end if
      if (lpmx(ispe)>1) then
        rhoaugmom(10,na)=rhoaugmom(10,na)+sum(10)
        rhoaugmom(11,na)=rhoaugmom(11,na)+sum(11)
        rhoaugmom(12,na)=rhoaugmom(12,na)+sum(12)
        rhoaugmom(13,na)=rhoaugmom(13,na)+sum(13)
        rhoaugmom(14,na)=rhoaugmom(14,na)+sum(14)
        rhoaugmom(15,na)=rhoaugmom(15,na)+sum(15)
        rhoaugmom(16,na)=rhoaugmom(16,na)+sum(16)
        rhoaugmom(17,na)=rhoaugmom(17,na)+sum(17)
        rhoaugmom(18,na)=rhoaugmom(18,na)+sum(18)
        rhoaugmom(19,na)=rhoaugmom(19,na)+sum(19)
        rhoaugmom(20,na)=rhoaugmom(20,na)+sum(20)
        rhoaugmom(21,na)=rhoaugmom(21,na)+sum(21)
        rhoaugmom(22,na)=rhoaugmom(22,na)+sum(22)
        rhoaugmom(23,na)=rhoaugmom(23,na)+sum(23)
        rhoaugmom(24,na)=rhoaugmom(24,na)+sum(24)
        rhoaugmom(25,na)=rhoaugmom(25,na)+sum(25)
      end if
!$omp end critical
!$omp barrier
    end if
  end do

  deallocate(sum,tmp)
  return
end subroutine


subroutine augcharge_02( &
 key_natpri_in,key_pp_paw,                          & ! <
 nradmx,npoint,lrhomx,lmx,natom,num_atcell,num_spe, & ! <
 indspe,natpri,natpri_inf,ntyppp,nradct,lpmx,nwexp, & ! <
 yylm,radial,rfac,                                  & ! <
 rhoaugmom,                                         & ! <
 rhoaugr)                                             ! >
use mod_mpi
implicit none
integer,intent(in) :: key_natpri_in,key_pp_paw
integer,intent(in) :: nradmx,npoint,lrhomx,lmx,natom,num_atcell,num_spe
integer,intent(in) :: indspe(natom),natpri(natom),natpri_inf(natom)
integer,intent(in) :: ntyppp(num_spe),nradct(num_spe),lpmx(num_spe),nwexp(num_spe)
real*8, intent(in) :: yylm(npoint,lrhomx),radial(nradmx,num_spe),rfac(0:2*(lmx-1),num_spe)
real*8, intent(in) :: rhoaugmom(lrhomx,natom)
real*8, intent(out):: rhoaugr(nradmx,npoint,num_atcell)
integer na,ipri,il,ir,ispe,i,i_end
real*8, allocatable::tmp(:)
real*8  aug0,aug1,aug2,aug3,aug4,augbaser0,augbaser1,augbaser2,augbaser3,augbaser4 &
       ,rcut,rcutin,rcut2in,rcut3in,rcut4in,r,r2,r3,r4

  allocate(tmp(25))

  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
      ispe=indspe(na)
      if (lpmx(ispe)==0) i_end=1
      if (lpmx(ispe)==1) i_end=9
      if (lpmx(ispe)==2) i_end=25
      rcut=radial(nradct(ispe),ispe)
      rcutin=1.0d0/rcut
      rcut2in=rcutin*rcutin
      rcut3in=rcut2in*rcutin
      rcut4in=rcut3in*rcutin
      ipri=natpri_inf(na)
!!$omp do
!      do i=1,nradmx*npoint
!        rhoaugr(i,1,ipri)=0.0d0
!      end do
!$omp do
      do il=1,npoint
        do i=1,i_end
          tmp(i)=rhoaugmom(i,na)*yylm(il,i)
        end do
        select case(lpmx(ispe))
        case(0)
          do ir=2,nradct(ispe)-1
            r=radial(ir,ispe)
            r2=r*r
            aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
            aug3=aug4*(1.0d0-r2*rcut2in)
            aug2=aug3*(1.0d0-r2*rcut2in)
            aug1=aug2*(1.0d0-r2*rcut2in)
            aug0=aug1*(1.0d0-r2*rcut2in)
            augbaser0=rfac(0,ispe)*aug0
            rhoaugr(ir,il,ipri)= tmp(1)*augbaser0
          end do
        case(1)
          do ir=2,nradct(ispe)-1
            r=radial(ir,ispe)
            r2=r*r
            aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
            aug3=aug4*(1.0d0-r2*rcut2in)
            aug2=aug3*(1.0d0-r2*rcut2in)
            aug1=aug2*(1.0d0-r2*rcut2in)
            aug0=aug1*(1.0d0-r2*rcut2in)
            augbaser0=rfac(0,ispe)*aug0
            augbaser1=rfac(1,ispe)*r *rcutin *aug1
            augbaser2=rfac(2,ispe)*r2*rcut2in*aug2
            rhoaugr(ir,il,ipri)=   tmp( 1)*augbaser0 &
                +tmp( 2)*augbaser1+tmp( 3)*augbaser1+tmp( 4)*augbaser1 &
                +tmp( 5)*augbaser2+tmp( 6)*augbaser2+tmp( 7)*augbaser2 &
                +tmp( 8)*augbaser2+tmp( 9)*augbaser2
          end do
        case(2)
          do ir=2,nradct(ispe)-1
            r=radial(ir,ispe)
            r2=r*r
            r3=r2*r
            r4=r3*r
            aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
            aug3=aug4*(1.0d0-r2*rcut2in)
            aug2=aug3*(1.0d0-r2*rcut2in)
            aug1=aug2*(1.0d0-r2*rcut2in)
            aug0=aug1*(1.0d0-r2*rcut2in)
            augbaser0=rfac(0,ispe)*aug0
            augbaser1=rfac(1,ispe)*r *rcutin *aug1
            augbaser2=rfac(2,ispe)*r2*rcut2in*aug2
            augbaser3=rfac(3,ispe)*r3*rcut3in*aug3
            augbaser4=rfac(4,ispe)*r4*rcut4in*aug4
            rhoaugr(ir,il,ipri)=   tmp( 1)*augbaser0 &
                +tmp( 2)*augbaser1+tmp( 3)*augbaser1+tmp( 4)*augbaser1 &
                +tmp( 5)*augbaser2+tmp( 6)*augbaser2+tmp( 7)*augbaser2 &
                +tmp( 8)*augbaser2+tmp( 9)*augbaser2+tmp(10)*augbaser3 &
                +tmp(11)*augbaser3+tmp(12)*augbaser3+tmp(13)*augbaser3 &
                +tmp(14)*augbaser3+tmp(15)*augbaser3+tmp(16)*augbaser3 &
                +tmp(17)*augbaser4+tmp(18)*augbaser4+tmp(19)*augbaser4 &
                +tmp(20)*augbaser4+tmp(21)*augbaser4+tmp(22)*augbaser4 &
                +tmp(23)*augbaser4+tmp(24)*augbaser4+tmp(25)*augbaser4
          end do
        end select
      end do
    end if
  end do

  deallocate(tmp)
  return
end subroutine


subroutine augcharge_03( &
 key_natpri_inps,key_pp_paw,                                                                & ! <
 nradmx,lrhomx,lmx,natom,num_spe,num_list_d,num_ppcell_d,                                   & ! <
 indspe,natprid,napsd,natinfd,ndatx,ndaty,ndatz,ntyppp,nradct,lpmx,nwexp,lstdx,lstdy,lstdz, & ! <
 ddx,ddy,ddz,radial,rfac,atx,aty,atz,rhoaugmom,                                             & ! <
 rhoaug3d,drhoaug3ddx,drhoaug3ddy,drhoaug3ddz)                                                ! >
use mod_mpi
implicit none
integer,intent(in) :: key_natpri_inps,key_pp_paw
integer,intent(in) :: nradmx,lrhomx,lmx,natom,num_spe,num_list_d,num_ppcell_d
integer,intent(in) :: indspe(natom),natprid(natom),napsd(natom),natinfd(natom)
integer,intent(in) :: ndatx(natom),ndaty(natom),ndatz(natom)
integer,intent(in) :: ntyppp(num_spe),nradct(num_spe),lpmx(num_spe),nwexp(num_spe)
integer,intent(in) :: lstdx(num_list_d,num_ppcell_d),lstdy(num_list_d,num_ppcell_d),lstdz(num_list_d,num_ppcell_d)
real*8, intent(in) :: ddx,ddy,ddz
real*8, intent(in) :: radial(nradmx,num_spe),rfac(0:2*(lmx-1),num_spe)
real*8, intent(in) :: atx(natom),aty(natom),atz(natom)
real*8, intent(in) :: rhoaugmom(lrhomx,natom)
real*8, intent(out):: rhoaug3d(num_list_d,num_ppcell_d)
real*8, intent(out):: drhoaug3ddx(num_list_d,num_ppcell_d)
real*8, intent(out):: drhoaug3ddy(num_list_d,num_ppcell_d)
real*8, intent(out):: drhoaug3ddz(num_list_d,num_ppcell_d)
integer ix,iy,iz,na,iapsd,ispe,i
real*8 rcut,rcutin,rcut2in,rcut3in,rcut4in,rcut5in,rcut6in &
      ,r,r2,r3,r4,r5,r6,rin,r2in,r3in,r4in,r5in,r6in
real*8 x,y,z,x2,y2,z2,x3,y3,z3,xy,yz,zx,xyz
real*8 ylm01,ylm02,ylm03,ylm04,ylm05,ylm06,ylm07,ylm08,ylm09,ylm10 &
      ,ylm11,ylm12,ylm13,ylm14,ylm15,ylm16,ylm17,ylm18,ylm19,ylm20 &
      ,ylm21,ylm22,ylm23,ylm24,ylm25
real*8 dylmdx01,dylmdx02,dylmdx03,dylmdx04,dylmdx05,dylmdx06,dylmdx07,dylmdx08,dylmdx09,dylmdx10 &
      ,dylmdx11,dylmdx12,dylmdx13,dylmdx14,dylmdx15,dylmdx16,dylmdx17,dylmdx18,dylmdx19,dylmdx20 &
      ,dylmdx21,dylmdx22,dylmdx23,dylmdx24,dylmdx25
real*8 dylmdy01,dylmdy02,dylmdy03,dylmdy04,dylmdy05,dylmdy06,dylmdy07,dylmdy08,dylmdy09,dylmdy10 &
      ,dylmdy11,dylmdy12,dylmdy13,dylmdy14,dylmdy15,dylmdy16,dylmdy17,dylmdy18,dylmdy19,dylmdy20 &
      ,dylmdy21,dylmdy22,dylmdy23,dylmdy24,dylmdy25
real*8 dylmdz01,dylmdz02,dylmdz03,dylmdz04,dylmdz05,dylmdz06,dylmdz07,dylmdz08,dylmdz09,dylmdz10 &
      ,dylmdz11,dylmdz12,dylmdz13,dylmdz14,dylmdz15,dylmdz16,dylmdz17,dylmdz18,dylmdz19,dylmdz20 &
      ,dylmdz21,dylmdz22,dylmdz23,dylmdz24,dylmdz25
real*8 aug0,aug1,aug2,aug3,aug4,aug5
real*8 augbase3d0,augbase3d1,augbase3d2,augbase3d3,augbase3d4
real*8 daugbase3ddr0,daugbase3ddr1,daugbase3ddr2,daugbase3ddr3,daugbase3ddr4
real*8 tmp
real*8 atmtmpx,atmtmpy,atmtmpz
real*8 pi,sqfourpi,sq3fourpi,sq1516pi,sq0516pi,sq1504pi,sq0716pi,sq2132pi,sq10516pi &
      ,sq3532pi,sq10504pi,sq09256pi,sq4532pi,sq4564pi,sq31532pi,sq315256pi,sq4516pi,sq31516pi

  pi=4.0d0*datan(1.0d0)
  sqfourpi=dsqrt(1.0d0/4.0d0/pi)
  sq3fourpi=dsqrt(3.0d0/4.0d0/pi)
  sq1516pi=dsqrt(15.0d0/16.0d0/pi)
  sq0516pi=dsqrt( 5.0d0/16.0d0/pi)
  sq1504pi=dsqrt(15.0d0/4.0d0/pi)
  sq0716pi=dsqrt(7.0d0/16.0d0/pi)
  sq2132pi=dsqrt(21.0d0/32.0d0/pi)
  sq10516pi=dsqrt(105.0d0/16.0d0/pi)
  sq3532pi=dsqrt(35.0d0/32.0d0/pi)
  sq10504pi=dsqrt(105.0d0/4.0d0/pi)
  sq09256pi=dsqrt(9.0d0/256.0d0/pi)
  sq4532pi=dsqrt(45.0d0/32.0d0/pi)
  sq4564pi=dsqrt(45.0d0/64.0d0/pi)
  sq31532pi=dsqrt(315.0d0/32.0d0/pi)
  sq315256pi=dsqrt(315.0d0/256.0d0/pi)
  sq4516pi=dsqrt(45.0d0/16.0d0/pi)
  sq31516pi=dsqrt(315.0d0/16.0d0/pi)

!$omp do
  do i=1,num_list_d*num_ppcell_d
       rhoaug3d(i,1)=0.0d0
    drhoaug3ddx(i,1)=0.0d0
    drhoaug3ddy(i,1)=0.0d0
    drhoaug3ddz(i,1)=0.0d0
  end do

  do na=1,natom
  if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natprid(na) .eq. key_natpri_inps)) then
    iapsd=napsd(na)
    ispe=indspe(na)
    rcut=radial(nradct(ispe),ispe)
    rcutin=1.0d0/rcut
    rcut2in=rcutin*rcutin
    rcut3in=rcut2in*rcutin
    rcut4in=rcut3in*rcutin
    rcut5in=rcut4in*rcutin
    rcut6in=rcut5in*rcutin

    select case(lpmx(ispe))
    case(0)
!$omp do
      do i=1,natinfd(na)
        ix=lstdx(i,iapsd)
        iy=lstdy(i,iapsd)
        iz=lstdz(i,iapsd)
        atmtmpx=atx(na)-(ndatx(na)*ddx-0.5d0*ddx)
        atmtmpy=aty(na)-(ndaty(na)*ddy-0.5d0*ddy)
        atmtmpz=atz(na)-(ndatz(na)*ddz-0.5d0*ddz)
        x=ix*ddx-atmtmpx
        y=iy*ddy-atmtmpy
        z=iz*ddz-atmtmpz
        x2=x*x
        y2=y*y
        z2=z*z
        r2=x2+y2+z2
        r=dsqrt(r2)
!        if (r .le. radial(nradct(ispe),ispe)) then
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)
          aug5=aug4/(1.0d0-r2*rcut2in)
          ylm01= sqfourpi
          augbase3d0= rfac(0,ispe)*aug0
          daugbase3ddr0= -2.0d0*rfac(0,ispe)*dfloat(nwexp(ispe))*r*rcut2in*aug1
          rhoaug3d(i,iapsd)= rhoaugmom( 1,na)*ylm01*augbase3d0
          tmp              = rhoaugmom( 1,na)*ylm01*daugbase3ddr0
          drhoaug3ddx(i,iapsd)=tmp*(x/r)
          drhoaug3ddy(i,iapsd)=tmp*(y/r)
          drhoaug3ddz(i,iapsd)=tmp*(z/r)
!        end if
      end do
    case(1)
!$omp do
      do i=1,natinfd(na)
        ix=lstdx(i,iapsd)
        iy=lstdy(i,iapsd)
        iz=lstdz(i,iapsd)
        atmtmpx=atx(na)-(ndatx(na)*ddx-0.5d0*ddx)
        atmtmpy=aty(na)-(ndaty(na)*ddy-0.5d0*ddy)
        atmtmpz=atz(na)-(ndatz(na)*ddz-0.5d0*ddz)
        x=ix*ddx-atmtmpx
        y=iy*ddy-atmtmpy
        z=iz*ddz-atmtmpz
        x2=x*x
        y2=y*y
        z2=z*z
        r2=x2+y2+z2
        r=dsqrt(r2)
!        if (r .le. radial(nradct(ispe),ispe)) then
          xy=x*y
          yz=y*z
          zx=z*x
          r3=r2*r
          rin=1.0d0/r
          r2in=rin*rin
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)
          aug5=aug4/(1.0d0-r2*rcut2in)
          ylm01= sqfourpi
          ylm02= sq3fourpi*x*rin
          ylm03= sq3fourpi*y*rin
          ylm04= sq3fourpi*z*rin
          ylm05= sq1516pi*(x2-y2)*r2in
          ylm06=-sq1504pi*zx*r2in
          ylm07= sq0516pi*(3.0d0*z2-r2)*r2in
          ylm08=-sq1504pi*yz*r2in
          ylm09= sq1504pi*xy*r2in
          augbase3d0= rfac(0,ispe)*aug0
          daugbase3ddr0= -2.0d0*rfac(0,ispe)*dfloat(nwexp(ispe))*r*rcut2in*aug1
          augbase3d1=rfac(1,ispe)*r *rcutin *aug1
          augbase3d2=rfac(2,ispe)*r2*rcut2in*aug2
          daugbase3ddr1=rfac(1,ispe)         *rcutin *aug1 &
                 -2.0d0*rfac(1,ispe)*dfloat(nwexp(ispe)-1)  *r2*rcut3in*aug2
          daugbase3ddr2=rfac(2,ispe)*2.0d0*r *rcut2in*aug2 &
                 -2.0d0*rfac(2,ispe)*dfloat(nwexp(ispe)-2)  *r3*rcut4in*aug3
          rhoaug3d(i,iapsd)= rhoaugmom( 1,na)*ylm01*augbase3d0 &
                            +rhoaugmom( 2,na)*ylm02*augbase3d1 &
                            +rhoaugmom( 3,na)*ylm03*augbase3d1 &
                            +rhoaugmom( 4,na)*ylm04*augbase3d1 &
                            +rhoaugmom( 5,na)*ylm05*augbase3d2 &
                            +rhoaugmom( 6,na)*ylm06*augbase3d2 &
                            +rhoaugmom( 7,na)*ylm07*augbase3d2 &
                            +rhoaugmom( 8,na)*ylm08*augbase3d2 &
                            +rhoaugmom( 9,na)*ylm09*augbase3d2
          tmp= rhoaugmom( 1,na)*ylm01*daugbase3ddr0 &
              +rhoaugmom( 2,na)*ylm02*daugbase3ddr1 &
              +rhoaugmom( 3,na)*ylm03*daugbase3ddr1 &
              +rhoaugmom( 4,na)*ylm04*daugbase3ddr1 &
              +rhoaugmom( 5,na)*ylm05*daugbase3ddr2 &
              +rhoaugmom( 6,na)*ylm06*daugbase3ddr2 &
              +rhoaugmom( 7,na)*ylm07*daugbase3ddr2 &
              +rhoaugmom( 8,na)*ylm08*daugbase3ddr2 &
              +rhoaugmom( 9,na)*ylm09*daugbase3ddr2
          drhoaug3ddx(i,iapsd)=tmp*(x/r)
          drhoaug3ddy(i,iapsd)=tmp*(y/r)
          drhoaug3ddz(i,iapsd)=tmp*(z/r)
!        end if
      end do
    case(2)
!$omp do
      do i=1,natinfd(na)
        ix=lstdx(i,iapsd)
        iy=lstdy(i,iapsd)
        iz=lstdz(i,iapsd)
        atmtmpx=atx(na)-(ndatx(na)*ddx-0.5d0*ddx)
        atmtmpy=aty(na)-(ndaty(na)*ddy-0.5d0*ddy)
        atmtmpz=atz(na)-(ndatz(na)*ddz-0.5d0*ddz)
        x=ix*ddx-atmtmpx
        y=iy*ddy-atmtmpy
        z=iz*ddz-atmtmpz
        x2=x*x
        y2=y*y
        z2=z*z
        r2=x2+y2+z2
        r=dsqrt(r2)
!        if (r .le. radial(nradct(ispe),ispe)) then
          xy=x*y
          yz=y*z
          zx=z*x
          xyz=xy*z
          r3=r2*r
          r4=r3*r
          r5=r4*r
          r6=r5*r
          rin=1.0d0/r
          r2in=rin*rin
          r3in=r2in*rin
          r4in=r3in*rin
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)
          aug5=aug4/(1.0d0-r2*rcut2in)
          ylm01= sqfourpi
          ylm02= sq3fourpi*x*rin
          ylm03= sq3fourpi*y*rin
          ylm04= sq3fourpi*z*rin
          ylm05= sq1516pi*(x2-y2)*r2in
          ylm06=-sq1504pi*zx*r2in
          ylm07= sq0516pi*(3.0d0*z2-r2)*r2in
          ylm08=-sq1504pi*yz*r2in
          ylm09= sq1504pi*xy*r2in
          ylm10= sq0716pi*z*(5.0d0*z2-3.0d0*r2)*r3in
          ylm11=-sq2132pi*x*(5.0d0*z2-r2)*r3in
          ylm12= sq10516pi*z*(x2-y2)*r3in
          ylm13=-sq3532pi*x*(x2-3.0d0*y2)*r3in
          ylm14=-sq2132pi*y*(5.0d0*z2-r2)*r3in
          ylm15= sq10504pi*xyz*r3in
          ylm16=-sq3532pi*y*(3.0d0*x2-y2)*r3in
          ylm17= sq09256pi*(35.0d0*z2*z2-30.0d0*z2*r2+3.0d0*r2*r2)*r4in
          ylm18= sq4532pi*x*z*(7.0d0*z2-3.0d0*r2)*r4in
          ylm19= sq4564pi*(x2-y2)*(7.0d0*z2-r2)*r4in
          ylm20= sq31532pi*zx*(x2-3.0d0*y2)*r4in
          ylm21= sq315256pi*(x2*x2-6.0d0*x2*y2+y2*y2)*r4in
          ylm22= sq4532pi*yz*(7.0d0*z2-3.0d0*r2)*r4in
          ylm23= sq4516pi*xy*(7.0d0*z2-r2)*r4in
          ylm24= sq31532pi*yz*(y2-3.0d0*x2)*r4in
          ylm25= sq31516pi*xy*(x2-y2)*r4in
          augbase3d0= rfac(0,ispe)*aug0
          augbase3d1=rfac(1,ispe)*r *rcutin *aug1
          augbase3d2=rfac(2,ispe)*r2*rcut2in*aug2
          augbase3d3=rfac(3,ispe)*r3*rcut3in*aug3
          augbase3d4=rfac(4,ispe)*r4*rcut4in*aug4
          daugbase3ddr0= -2.0d0*rfac(0,ispe)*dfloat(nwexp(ispe))*r*rcut2in*aug1
          daugbase3ddr1=rfac(1,ispe)         *rcutin *aug1 &
                 -2.0d0*rfac(1,ispe)*dfloat(nwexp(ispe)-1)  *r2*rcut3in*aug2
          daugbase3ddr2=rfac(2,ispe)*2.0d0*r *rcut2in*aug2 &
                 -2.0d0*rfac(2,ispe)*dfloat(nwexp(ispe)-2)  *r3*rcut4in*aug3
          daugbase3ddr3=rfac(3,ispe)*3.0d0*r2*rcut3in*aug3 &
                 -2.0d0*rfac(3,ispe)*dfloat(nwexp(ispe)-3)  *r4*rcut5in*aug4
          daugbase3ddr4=rfac(4,ispe)*4.0d0*r3*rcut4in*aug4 &
                 -2.0d0*rfac(4,ispe)*dfloat(nwexp(ispe)-4)  *r5*rcut6in*aug5
          rhoaug3d(i,iapsd)= rhoaugmom( 1,na)*ylm01*augbase3d0 &
                            +rhoaugmom( 2,na)*ylm02*augbase3d1 &
                            +rhoaugmom( 3,na)*ylm03*augbase3d1 &
                            +rhoaugmom( 4,na)*ylm04*augbase3d1 &
                            +rhoaugmom( 5,na)*ylm05*augbase3d2 &
                            +rhoaugmom( 6,na)*ylm06*augbase3d2 &
                            +rhoaugmom( 7,na)*ylm07*augbase3d2 &
                            +rhoaugmom( 8,na)*ylm08*augbase3d2 &
                            +rhoaugmom( 9,na)*ylm09*augbase3d2 &
                            +rhoaugmom(10,na)*ylm10*augbase3d3 &
                            +rhoaugmom(11,na)*ylm11*augbase3d3 &
                            +rhoaugmom(12,na)*ylm12*augbase3d3 &
                            +rhoaugmom(13,na)*ylm13*augbase3d3 &
                            +rhoaugmom(14,na)*ylm14*augbase3d3 &
                            +rhoaugmom(15,na)*ylm15*augbase3d3 &
                            +rhoaugmom(16,na)*ylm16*augbase3d3 &
                            +rhoaugmom(17,na)*ylm17*augbase3d4 &
                            +rhoaugmom(18,na)*ylm18*augbase3d4 &
                            +rhoaugmom(19,na)*ylm19*augbase3d4 &
                            +rhoaugmom(20,na)*ylm20*augbase3d4 &
                            +rhoaugmom(21,na)*ylm21*augbase3d4 &
                            +rhoaugmom(22,na)*ylm22*augbase3d4 &
                            +rhoaugmom(23,na)*ylm23*augbase3d4 &
                            +rhoaugmom(24,na)*ylm24*augbase3d4 &
                            +rhoaugmom(25,na)*ylm25*augbase3d4
          tmp= rhoaugmom( 1,na)*ylm01*daugbase3ddr0 &
              +rhoaugmom( 2,na)*ylm02*daugbase3ddr1 &
              +rhoaugmom( 3,na)*ylm03*daugbase3ddr1 &
              +rhoaugmom( 4,na)*ylm04*daugbase3ddr1 &
              +rhoaugmom( 5,na)*ylm05*daugbase3ddr2 &
              +rhoaugmom( 6,na)*ylm06*daugbase3ddr2 &
              +rhoaugmom( 7,na)*ylm07*daugbase3ddr2 &
              +rhoaugmom( 8,na)*ylm08*daugbase3ddr2 &
              +rhoaugmom( 9,na)*ylm09*daugbase3ddr2 &
              +rhoaugmom(10,na)*ylm10*daugbase3ddr3 &
              +rhoaugmom(11,na)*ylm11*daugbase3ddr3 &
              +rhoaugmom(12,na)*ylm12*daugbase3ddr3 &
              +rhoaugmom(13,na)*ylm13*daugbase3ddr3 &
              +rhoaugmom(14,na)*ylm14*daugbase3ddr3 &
              +rhoaugmom(15,na)*ylm15*daugbase3ddr3 &
              +rhoaugmom(16,na)*ylm16*daugbase3ddr3 &
              +rhoaugmom(17,na)*ylm17*daugbase3ddr4 &
              +rhoaugmom(18,na)*ylm18*daugbase3ddr4 &
              +rhoaugmom(19,na)*ylm19*daugbase3ddr4 &
              +rhoaugmom(20,na)*ylm20*daugbase3ddr4 &
              +rhoaugmom(21,na)*ylm21*daugbase3ddr4 &
              +rhoaugmom(22,na)*ylm22*daugbase3ddr4 &
              +rhoaugmom(23,na)*ylm23*daugbase3ddr4 &
              +rhoaugmom(24,na)*ylm24*daugbase3ddr4 &
              +rhoaugmom(25,na)*ylm25*daugbase3ddr4
          drhoaug3ddx(i,iapsd)=tmp*(x/r)
          drhoaug3ddy(i,iapsd)=tmp*(y/r)
          drhoaug3ddz(i,iapsd)=tmp*(z/r)
!        end if
      end do
    end select

    select case(lpmx(ispe))
    case(0)
!$omp do
      do i=1,natinfd(na)
        ix=lstdx(i,iapsd)
        iy=lstdy(i,iapsd)
        iz=lstdz(i,iapsd)
        atmtmpx=atx(na)-(ndatx(na)*ddx-0.5d0*ddx)
        atmtmpy=aty(na)-(ndaty(na)*ddy-0.5d0*ddy)
        atmtmpz=atz(na)-(ndatz(na)*ddz-0.5d0*ddz)
        x=ix*ddx-atmtmpx
        y=iy*ddy-atmtmpy
        z=iz*ddz-atmtmpz
        x2=x*x
        y2=y*y
        z2=z*z
        r2=x2+y2+z2
        r=dsqrt(r2)
!        if (r .le. radial(nradct(ispe),ispe)) then
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)
          dylmdx01= 0.0d0
          augbase3d0=rfac(0,ispe)*aug0
          drhoaug3ddx(i,iapsd)=drhoaug3ddx(i,iapsd) +rhoaugmom(1,na)*dylmdx01*augbase3d0
!        end if
      end do
    case(1)
!$omp do
      do i=1,natinfd(na)
        ix=lstdx(i,iapsd)
        iy=lstdy(i,iapsd)
        iz=lstdz(i,iapsd)
        atmtmpx=atx(na)-(ndatx(na)*ddx-0.5d0*ddx)
        atmtmpy=aty(na)-(ndaty(na)*ddy-0.5d0*ddy)
        atmtmpz=atz(na)-(ndatz(na)*ddz-0.5d0*ddz)
        x=ix*ddx-atmtmpx
        y=iy*ddy-atmtmpy
        z=iz*ddz-atmtmpz
        x2=x*x
        y2=y*y
        z2=z*z
        r2=x2+y2+z2
        r=dsqrt(r2)
!        if (r .le. radial(nradct(ispe),ispe)) then
          xyz=x*y*z
          xy=x*y
          zx=x*z
          rin=1.0d0/r
          r2in=rin*rin
          r3in=r2in*rin
          r4in=r3in*rin
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)
          dylmdx01= 0.0d0
          dylmdx02= sq3fourpi*(r2-x2)*r3in
          dylmdx03= sq3fourpi*(-xy)*r3in
          dylmdx04= sq3fourpi*(-zx)*r3in
          dylmdx05= sq1516pi*2.0d0*x*(r2-x2+y2)*r4in
          dylmdx06=-sq1504pi*z*(r2-2.0d0*x2)*r4in
          dylmdx07= sq0516pi*(-6.0d0*x*z2)*r4in
          dylmdx08= sq1504pi*2.0d0*xyz*r4in
          dylmdx09= sq1504pi*y*(r2-2.0d0*x2)*r4in
          augbase3d0=rfac(0,ispe)*aug0
          augbase3d1=rfac(1,ispe)*r *rcutin *aug1
          augbase3d2=rfac(2,ispe)*r2*rcut2in*aug2
          drhoaug3ddx(i,iapsd)=drhoaug3ddx(i,iapsd) &
                       +rhoaugmom( 1,na)*dylmdx01*augbase3d0 &
                       +rhoaugmom( 2,na)*dylmdx02*augbase3d1 &
                       +rhoaugmom( 3,na)*dylmdx03*augbase3d1 &
                       +rhoaugmom( 4,na)*dylmdx04*augbase3d1 &
                       +rhoaugmom( 5,na)*dylmdx05*augbase3d2 &
                       +rhoaugmom( 6,na)*dylmdx06*augbase3d2 &
                       +rhoaugmom( 7,na)*dylmdx07*augbase3d2 &
                       +rhoaugmom( 8,na)*dylmdx08*augbase3d2 &
                       +rhoaugmom( 9,na)*dylmdx09*augbase3d2
!        end if
      end do
    case(2)
!$omp do
      do i=1,natinfd(na)
        ix=lstdx(i,iapsd)
        iy=lstdy(i,iapsd)
        iz=lstdz(i,iapsd)
        atmtmpx=atx(na)-(ndatx(na)*ddx-0.5d0*ddx)
        atmtmpy=aty(na)-(ndaty(na)*ddy-0.5d0*ddy)
        atmtmpz=atz(na)-(ndatz(na)*ddz-0.5d0*ddz)
        x=ix*ddx-atmtmpx
        y=iy*ddy-atmtmpy
        z=iz*ddz-atmtmpz
        x2=x*x
        y2=y*y
        z2=z*z
        r2=x2+y2+z2
        r=dsqrt(r2)
!        if (r .le. radial(nradct(ispe),ispe)) then
          xyz=x*y*z
          xy=x*y
          zx=x*z
          x3=x2*x
          y3=y2*y
          z3=z2*z
          yz=y*z
          r3=r2*r
          r4=r3*r
          r5=r4*r
          r6=r5*r
          rin=1.0d0/r
          r2in=rin*rin
          r3in=r2in*rin
          r4in=r3in*rin
          r5in=r4in*rin
          r6in=r5in*rin
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)
          dylmdx01= 0.0d0
          dylmdx02= sq3fourpi*(r2-x2)*r3in
          dylmdx03= sq3fourpi*(-xy)*r3in
          dylmdx04= sq3fourpi*(-zx)*r3in
          dylmdx05= sq1516pi*2.0d0*x*(r2-x2+y2)*r4in
          dylmdx06=-sq1504pi*z*(r2-2.0d0*x2)*r4in
          dylmdx07= sq0516pi*(-6.0d0*x*z2)*r4in
          dylmdx08= sq1504pi*2.0d0*xyz*r4in
          dylmdx09= sq1504pi*y*(r2-2.0d0*x2)*r4in
          dylmdx10= sq0716pi*(3.0d0*zx*r2-15.0d0*x*z3)*r5in
          dylmdx11=-sq2132pi*(5.0d0*z2*r2-r4+x2*r2-15.0d0*x2*z2)*r5in
          dylmdx12= sq10516pi*(2.0d0*zx*r2-3.0d0*x3*z+3.0d0*x*y2*z)*r5in
          dylmdx13=-sq3532pi*(3.0d0*x2*r2-3.0d0*y2*r2-3.0d0*x2*x2+9.0d0*x2*y2)*r5in
          dylmdx14=-sq2132pi*(xy*r2-15.0d0*xy*z2)*r5in
          dylmdx15= sq10504pi*(yz*r2-3.0d0*x2*yz)*r5in
          dylmdx16=-sq3532pi*(6.0d0*xy*r2-9.0d0*x3*y+3.0d0*x*y3)*r5in
          dylmdx17= sq09256pi*(-140.0d0*x*z2*z2+60.0d0*x*z2*r2)*r6in
          dylmdx18= sq4532pi*(7.0d0*z3*r2-3.0d0*z*r4+6.0d0*x2*z*r2-28.0d0*x2*z3)*r6in
          dylmdx19= sq4564pi*(2.0d0*x*(7.0d0*z2-r2)*r4in-2.0d0*r*(x2-y2)*x*r5in-4.0d0*(x2-y2)*(7.0d0*z2-r2)*x*r6in)
          dylmdx20= sq31532pi*(3.0d0*x2*z*r2-3.0d0*y2*z*r2-4.0d0*(x2-3.0d0*y2)*x2*z)*r6in
          dylmdx21= sq315256pi*((4.0d0*x3-12.0d0*x*y2)*r4in-4.0d0*(x2*x2-6.0d0*x2*y2+y2*y2)*x*r6in)
          dylmdx22= sq4532pi*(-6.0d0*xyz*r*r5in-4.0d0*(7.0d0*z2-3.0d0*r2)*xyz*r6in)
          dylmdx23= sq4516pi*((7.0d0*z2-r2)*y*r4in-2.0d0*x2*y*r4in-4.0d0*(7.0d0*z2-r2)*x2*y*r6in)
          dylmdx24= sq31532pi*(-6.0d0*xyz*r4in-4.0d0*(y2-3.0d0*x2)*xyz*r6in)
          dylmdx25= sq31516pi*(2.0d0*x2*y*r4in+(x2-y2)*y*r4in-4.0d0*(x2-y2)*x2*y*r6in)
          augbase3d0=rfac(0,ispe)*aug0
          augbase3d1=rfac(1,ispe)*r *rcutin *aug1
          augbase3d2=rfac(2,ispe)*r2*rcut2in*aug2
          augbase3d3=rfac(3,ispe)*r3*rcut3in*aug3
          augbase3d4=rfac(4,ispe)*r4*rcut4in*aug4
            drhoaug3ddx(i,iapsd)=drhoaug3ddx(i,iapsd) &
                         +rhoaugmom( 1,na)*dylmdx01*augbase3d0 &
                         +rhoaugmom( 2,na)*dylmdx02*augbase3d1 &
                         +rhoaugmom( 3,na)*dylmdx03*augbase3d1 &
                         +rhoaugmom( 4,na)*dylmdx04*augbase3d1 &
                         +rhoaugmom( 5,na)*dylmdx05*augbase3d2 &
                         +rhoaugmom( 6,na)*dylmdx06*augbase3d2 &
                         +rhoaugmom( 7,na)*dylmdx07*augbase3d2 &
                         +rhoaugmom( 8,na)*dylmdx08*augbase3d2 &
                         +rhoaugmom( 9,na)*dylmdx09*augbase3d2 &
                         +rhoaugmom(10,na)*dylmdx10*augbase3d3 &
                         +rhoaugmom(11,na)*dylmdx11*augbase3d3 &
                         +rhoaugmom(12,na)*dylmdx12*augbase3d3 &
                         +rhoaugmom(13,na)*dylmdx13*augbase3d3 &
                         +rhoaugmom(14,na)*dylmdx14*augbase3d3 &
                         +rhoaugmom(15,na)*dylmdx15*augbase3d3 &
                         +rhoaugmom(16,na)*dylmdx16*augbase3d3 &
                         +rhoaugmom(17,na)*dylmdx17*augbase3d4 &
                         +rhoaugmom(18,na)*dylmdx18*augbase3d4 &
                         +rhoaugmom(19,na)*dylmdx19*augbase3d4 &
                         +rhoaugmom(20,na)*dylmdx20*augbase3d4 &
                         +rhoaugmom(21,na)*dylmdx21*augbase3d4 &
                         +rhoaugmom(22,na)*dylmdx22*augbase3d4 &
                         +rhoaugmom(23,na)*dylmdx23*augbase3d4 &
                         +rhoaugmom(24,na)*dylmdx24*augbase3d4 &
                         +rhoaugmom(25,na)*dylmdx25*augbase3d4
!        end if
      end do
    end select

    select case(lpmx(ispe))
    case(0)
!$omp do
      do i=1,natinfd(na)
        ix=lstdx(i,iapsd)
        iy=lstdy(i,iapsd)
        iz=lstdz(i,iapsd)
        atmtmpx=atx(na)-(ndatx(na)*ddx-0.5d0*ddx)
        atmtmpy=aty(na)-(ndaty(na)*ddy-0.5d0*ddy)
        atmtmpz=atz(na)-(ndatz(na)*ddz-0.5d0*ddz)
        x=ix*ddx-atmtmpx
        y=iy*ddy-atmtmpy
        z=iz*ddz-atmtmpz
        x2=x*x
        y2=y*y
        z2=z*z
        r2=x2+y2+z2
        r=dsqrt(r2)
!        if (r .le. radial(nradct(ispe),ispe)) then
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)
          dylmdy01= 0.0d0
          augbase3d0=rfac(0,ispe)*aug0
          drhoaug3ddy(i,iapsd)=drhoaug3ddy(i,iapsd) +rhoaugmom( 1,na)*dylmdy01*augbase3d0
!        end if
      end do
    case(1)
!$omp do
      do i=1,natinfd(na)
        ix=lstdx(i,iapsd)
        iy=lstdy(i,iapsd)
        iz=lstdz(i,iapsd)
        atmtmpx=atx(na)-(ndatx(na)*ddx-0.5d0*ddx)
        atmtmpy=aty(na)-(ndaty(na)*ddy-0.5d0*ddy)
        atmtmpz=atz(na)-(ndatz(na)*ddz-0.5d0*ddz)
        x=ix*ddx-atmtmpx
        y=iy*ddy-atmtmpy
        z=iz*ddz-atmtmpz
        x2=x*x
        y2=y*y
        z2=z*z
        r2=x2+y2+z2
        r=dsqrt(r2)
!        if (r .le. radial(nradct(ispe),ispe)) then
          xy=x*y
          yz=y*z
          r3=r2*r
          r4=r3*r
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)
          dylmdy01= 0.0d0
          dylmdy02= dsqrt(3.0d0/(4.0d0*pi))*(-xy)/r3
          dylmdy03= dsqrt(3.0d0/(4.0d0*pi))*(r2-y2)/r3
          dylmdy04= dsqrt(3.0d0/(4.0d0*pi))*(-yz)/r3
          dylmdy05= dsqrt(15.0d0/16.0d0/pi)*(-2.0d0*y)*(r2+x2-y2)/r4
          dylmdy06=-dsqrt(15.0d0/4.0d0/pi)*(-2.0d0*x*y*z)/r4
          dylmdy07= dsqrt(5.0d0/16.0d0/pi)*(-6.0d0*y*z2)/r4
          dylmdy08= dsqrt(15.0d0/4.0d0/pi)*(-z*(r2-2.0d0*y2)/r4)
          dylmdy09= dsqrt(15.0d0/4.0d0/pi)*x*(r2-2.0d0*y2)/r4
          augbase3d0=rfac(0,ispe)*aug0
          augbase3d1=rfac(1,ispe)*r *rcutin *aug1
          augbase3d2=rfac(2,ispe)*r2*rcut2in*aug2
          drhoaug3ddy(i,iapsd)=drhoaug3ddy(i,iapsd) &
                       +rhoaugmom( 1,na)*dylmdy01*augbase3d0 &
                       +rhoaugmom( 2,na)*dylmdy02*augbase3d1 &
                       +rhoaugmom( 3,na)*dylmdy03*augbase3d1 &
                       +rhoaugmom( 4,na)*dylmdy04*augbase3d1 &
                       +rhoaugmom( 5,na)*dylmdy05*augbase3d2 &
                       +rhoaugmom( 6,na)*dylmdy06*augbase3d2 &
                       +rhoaugmom( 7,na)*dylmdy07*augbase3d2 &
                       +rhoaugmom( 8,na)*dylmdy08*augbase3d2 &
                       +rhoaugmom( 9,na)*dylmdy09*augbase3d2
!        end if
      end do
    case(2)
!$omp do
      do i=1,natinfd(na)
        ix=lstdx(i,iapsd)
        iy=lstdy(i,iapsd)
        iz=lstdz(i,iapsd)
        atmtmpx=atx(na)-(ndatx(na)*ddx-0.5d0*ddx)
        atmtmpy=aty(na)-(ndaty(na)*ddy-0.5d0*ddy)
        atmtmpz=atz(na)-(ndatz(na)*ddz-0.5d0*ddz)
        x=ix*ddx-atmtmpx
        y=iy*ddy-atmtmpy
        z=iz*ddz-atmtmpz
        x2=x*x
        y2=y*y
        z2=z*z
        r2=x2+y2+z2
        r=dsqrt(r2)
!        if (r .le. radial(nradct(ispe),ispe)) then
          xy=x*y
          yz=y*z
          r3=r2*r
          r4=r3*r
          x3=x2*x
          y3=y2*y
          z3=z2*z
          xyz=x*y*z
          zx=x*z
          r3=r2*r
          r4=r3*r
          r5=r4*r
          r6=r5*r
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)
          dylmdy01= 0.0d0
          dylmdy02= dsqrt(3.0d0/(4.0d0*pi))*(-xy)/r3
          dylmdy03= dsqrt(3.0d0/(4.0d0*pi))*(r2-y2)/r3
          dylmdy04= dsqrt(3.0d0/(4.0d0*pi))*(-yz)/r3
          dylmdy05= dsqrt(15.0d0/16.0d0/pi)*(-2.0d0*y)*(r2+x2-y2)/r4
          dylmdy06=-dsqrt(15.0d0/4.0d0/pi)*(-2.0d0*x*y*z)/r4
          dylmdy07= dsqrt(5.0d0/16.0d0/pi)*(-6.0d0*y*z2)/r4
          dylmdy08= dsqrt(15.0d0/4.0d0/pi)*(-z*(r2-2.0d0*y2)/r4)
          dylmdy09= dsqrt(15.0d0/4.0d0/pi)*x*(r2-2.0d0*y2)/r4
          dylmdy10= dsqrt(7.0d0/16.0d0/pi)*(3.0d0*yz*r2-15.0d0*y*z3)/r5
          dylmdy11=-dsqrt(21.0d0/32.0d0/pi)*(xy*r2-15.0d0*xy*z2)/r5
          dylmdy12= dsqrt(105.0d0/16.0d0/pi)*(-2.0d0*yz*r2-3.0d0*x2*yz+3.0d0*y3*z)/r5
          dylmdy13=-dsqrt(35.0d0/32.0d0/pi)*(-6.0d0*xy*r2-3.0d0*x3*y+9.0d0*x*y3)/r5
          dylmdy14=-dsqrt(21.0d0/32.0d0/pi)*(5.0d0*z2*r2-r4+y2*r2-15.0d0*y2*z2)/r5
          dylmdy15= dsqrt(105.0d0/4.0d0/pi)*(zx*r2-3.0d0*zx*y2)/r5
          dylmdy16=-dsqrt(35.0d0/32.0d0/pi)*(3.0d0*x2*r2-3.0d0*y2*r2-9.0d0*x2*y2+3.0d0*y2*y2)/r5
          dylmdy17= dsqrt(9.0d0/256.0d0/pi)*(-140.0d0*y*z2*z2+60.0d0*y*z2*r2)/r6
          dylmdy18= dsqrt(45.0d0/32.0d0/pi)*(6.0d0*xyz*r2-28.0d0*xyz*z2)/r6
          dylmdy19= dsqrt(45.0d0/64.0d0/pi)*(-14.0d0*y*z2*r2+2.0d0*y*r4+2.0d0*x2*y*r2 &
                                                                     -2.0d0*y3*r2-28.0d0*x2*y*z2+28.0d0*y3*z2)/r6
          dylmdy20= dsqrt(315.0d0/32.0d0/pi)*(-6.0d0*xyz*r2-4.0d0*x3*yz+12.0d0*xyz*y2)/r6
          dylmdy21= dsqrt(315.0d0/256.0d0/pi)*(-12.0d0*x2*y*r2+4.0d0*y3*r2-4.0d0*x3*xy+24.0d0*x2*y3-4.0d0*y2*y3)/r6
          dylmdy22= dsqrt(45.0d0/32.0d0/pi)*(6.0d0*y2*z*r2+7.0d0*z3*r2-3.0d0*z*r4-28.0d0*y2*z3)/r6
          dylmdy23= dsqrt(45.0d0/16.0d0/pi)*(2.0d0*x*y2*r2+7.0d0*x*z2*r2-x*r4-28.0d0*x*y2*z2)/r6
          dylmdy24= dsqrt(315.0d0/32.0d0/pi)*(3.0d0*y2*z*r2-3.0d0*x2*z*r2-4.0d0*y2*y2*z+12.0d0*x2*y2*z)/r6
          dylmdy25= dsqrt(315.0d0/16.0d0/pi)*(-3.0d0*x*y2*r2+x3*r2-4.0d0*x3*y2+4.0d0*x*y2*y2)/r6
          augbase3d0=rfac(0,ispe)*aug0
          augbase3d1=rfac(1,ispe)*r *rcutin *aug1
          augbase3d2=rfac(2,ispe)*r2*rcut2in*aug2
          augbase3d3=rfac(3,ispe)*r3*rcut3in*aug3
          augbase3d4=rfac(4,ispe)*r4*rcut4in*aug4
          drhoaug3ddy(i,iapsd)=drhoaug3ddy(i,iapsd) &
                       +rhoaugmom( 1,na)*dylmdy01*augbase3d0 &
                       +rhoaugmom( 2,na)*dylmdy02*augbase3d1 &
                       +rhoaugmom( 3,na)*dylmdy03*augbase3d1 &
                       +rhoaugmom( 4,na)*dylmdy04*augbase3d1 &
                       +rhoaugmom( 5,na)*dylmdy05*augbase3d2 &
                       +rhoaugmom( 6,na)*dylmdy06*augbase3d2 &
                       +rhoaugmom( 7,na)*dylmdy07*augbase3d2 &
                       +rhoaugmom( 8,na)*dylmdy08*augbase3d2 &
                       +rhoaugmom( 9,na)*dylmdy09*augbase3d2 &
                       +rhoaugmom(10,na)*dylmdy10*augbase3d3 &
                       +rhoaugmom(11,na)*dylmdy11*augbase3d3 &
                       +rhoaugmom(12,na)*dylmdy12*augbase3d3 &
                       +rhoaugmom(13,na)*dylmdy13*augbase3d3 &
                       +rhoaugmom(14,na)*dylmdy14*augbase3d3 &
                       +rhoaugmom(15,na)*dylmdy15*augbase3d3 &
                       +rhoaugmom(16,na)*dylmdy16*augbase3d3 &
                       +rhoaugmom(17,na)*dylmdy17*augbase3d4 &
                       +rhoaugmom(18,na)*dylmdy18*augbase3d4 &
                       +rhoaugmom(19,na)*dylmdy19*augbase3d4 &
                       +rhoaugmom(20,na)*dylmdy20*augbase3d4 &
                       +rhoaugmom(21,na)*dylmdy21*augbase3d4 &
                       +rhoaugmom(22,na)*dylmdy22*augbase3d4 &
                       +rhoaugmom(23,na)*dylmdy23*augbase3d4 &
                       +rhoaugmom(24,na)*dylmdy24*augbase3d4 &
                       +rhoaugmom(25,na)*dylmdy25*augbase3d4
!        end if
      end do
    end select

    select case(lpmx(ispe))
    case(0)
!$omp do
      do i=1,natinfd(na)
        ix=lstdx(i,iapsd)
        iy=lstdy(i,iapsd)
        iz=lstdz(i,iapsd)
        atmtmpx=atx(na)-(ndatx(na)*ddx-0.5d0*ddx)
        atmtmpy=aty(na)-(ndaty(na)*ddy-0.5d0*ddy)
        atmtmpz=atz(na)-(ndatz(na)*ddz-0.5d0*ddz)
        x=ix*ddx-atmtmpx
        y=iy*ddy-atmtmpy
        z=iz*ddz-atmtmpz
        x2=x*x
        y2=y*y
        z2=z*z
        r2=x2+y2+z2
        r=dsqrt(r2)
!        if (r .le. radial(nradct(ispe),ispe)) then
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)
          dylmdz01= 0.0d0
          augbase3d0=rfac(0,ispe)*aug0
          drhoaug3ddz(i,iapsd)=drhoaug3ddz(i,iapsd) +rhoaugmom( 1,na)*dylmdz01*augbase3d0
!        end if
      end do
    case(1)
!$omp do
    do i=1,natinfd(na)
      ix=lstdx(i,iapsd)
      iy=lstdy(i,iapsd)
      iz=lstdz(i,iapsd)
      atmtmpx=atx(na)-(ndatx(na)*ddx-0.5d0*ddx)
      atmtmpy=aty(na)-(ndaty(na)*ddy-0.5d0*ddy)
      atmtmpz=atz(na)-(ndatz(na)*ddz-0.5d0*ddz)
      x=ix*ddx-atmtmpx
      y=iy*ddy-atmtmpy
      z=iz*ddz-atmtmpz
      x2=x*x
      y2=y*y
      z2=z*z
      r2=x2+y2+z2
      r=dsqrt(r2)
!      if (r .le. radial(nradct(ispe),ispe)) then
        yz=y*z
        zx=x*z
        xyz=x*y*z
        rin=1.0d0/r
        r2in=rin*rin
        r3in=r2in*rin
        r4in=r3in*rin
        aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
        aug3=aug4*(1.0d0-r2*rcut2in)
        aug2=aug3*(1.0d0-r2*rcut2in)
        aug1=aug2*(1.0d0-r2*rcut2in)
        aug0=aug1*(1.0d0-r2*rcut2in)
        dylmdz01= 0.0d0
        dylmdz02= sq3fourpi*(-zx)*r3in
        dylmdz03= sq3fourpi*(-yz)*r3in
        dylmdz04= sq3fourpi*(r2-z2)*r3in
        dylmdz05= sq1516pi*(-2.0d0*z)*(x2-y2)*r4in
        dylmdz06=-sq1504pi*(r2-2.0d0*z2)*x*r4in
        dylmdz07= sq0516pi*6.0d0*z*(r2-z2)*r4in
        dylmdz08= sq1504pi*(-y*(r2-2.0d0*z2)*r4in)
        dylmdz09= sq1504pi*(-2.0d0*xyz)*r4in
        augbase3d0=rfac(0,ispe)*aug0
        augbase3d1=rfac(1,ispe)*r *rcutin *aug1
        augbase3d2=rfac(2,ispe)*r2*rcut2in*aug2
        drhoaug3ddz(i,iapsd)=drhoaug3ddz(i,iapsd) &
                       +rhoaugmom( 1,na)*dylmdz01*augbase3d0 &
                       +rhoaugmom( 2,na)*dylmdz02*augbase3d1 &
                       +rhoaugmom( 3,na)*dylmdz03*augbase3d1 &
                       +rhoaugmom( 4,na)*dylmdz04*augbase3d1 &
                       +rhoaugmom( 5,na)*dylmdz05*augbase3d2 &
                       +rhoaugmom( 6,na)*dylmdz06*augbase3d2 &
                       +rhoaugmom( 7,na)*dylmdz07*augbase3d2 &
                       +rhoaugmom( 8,na)*dylmdz08*augbase3d2 &
                       +rhoaugmom( 9,na)*dylmdz09*augbase3d2
!      end if
    end do
    case(2)
!$omp do
      do i=1,natinfd(na)
        ix=lstdx(i,iapsd)
        iy=lstdy(i,iapsd)
        iz=lstdz(i,iapsd)
        atmtmpx=atx(na)-(ndatx(na)*ddx-0.5d0*ddx)
        atmtmpy=aty(na)-(ndaty(na)*ddy-0.5d0*ddy)
        atmtmpz=atz(na)-(ndatz(na)*ddz-0.5d0*ddz)
        x=ix*ddx-atmtmpx
        y=iy*ddy-atmtmpy
        z=iz*ddz-atmtmpz
        x2=x*x
        y2=y*y
        z2=z*z
        r2=x2+y2+z2
        r=dsqrt(r2)
!        if (r .le. radial(nradct(ispe),ispe)) then
          yz=y*z
          zx=x*z
          xyz=x*y*z
          x3=x2*x
          y3=y2*y
          z3=z2*z
          xy=x*y
          r3=r2*r
          r4=r3*r
          rin=1.0d0/r
          r2in=rin*rin
          r3in=r2in*rin
          r4in=r3in*rin
          r5in=r4in*rin
          r6in=r5in*rin
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)
          dylmdz01= 0.0d0
          dylmdz02= sq3fourpi*(-zx)*r3in
          dylmdz03= sq3fourpi*(-yz)*r3in
          dylmdz04= sq3fourpi*(r2-z2)*r3in
          dylmdz05= sq1516pi*(-2.0d0*z)*(x2-y2)*r4in
          dylmdz06=-sq1504pi*(r2-2.0d0*z2)*x*r4in
          dylmdz07= sq0516pi*6.0d0*z*(r2-z2)*r4in
          dylmdz08= sq1504pi*(-y*(r2-2.0d0*z2)*r4in)
          dylmdz09= sq1504pi*(-2.0d0*xyz)*r4in
          dylmdz10= sq0716pi*(18.0d0*z2*r2-15.0d0*z2*z2-3.0d0*r4)*r5in
          dylmdz11=-sq2132pi*(11.0d0*zx*r2-15.0d0*x*z3)*r5in
          dylmdz12= sq10516pi*(x2-y2)*(r2-3.0d0*z2)*r5in
          dylmdz13=-sq3532pi*(-3.0d0*(x2-3.0d0*y2)*zx*r5in)
          dylmdz14=-sq2132pi*(11.0d0*yz*r2-15.0d0*y*z3)*r5in
          dylmdz15= sq10504pi*(xy*r2-3.0d0*xy*z2)*r5in
          dylmdz16=-sq3532pi*(-3.0d0*(3.0d0*x2-y2)*yz*r5in)
          dylmdz17= sq09256pi*(200.0d0*z3*r2-60.0d0*z*r4-140.0d0*z2*z3)*r6in
          dylmdz18= sq4532pi*(27.0d0*x*z2*r2-3.0d0*x*r4-28.0d0*x*z2*z2)*r6in
          dylmdz19= sq4564pi*(16.0d0*x2*z*r2-16.0d0*y2*z*r2-28.0d0*x2*z3+28.0d0*y2*z3)*r6in
          dylmdz20= sq31532pi*(x3*r2-3.0d0*x*y2*r2-4.0d0*x3*z2+12.0d0*x*y2*z2)*r6in
          dylmdz21= sq315256pi*(-4.0d0*x2*x2*z+24.0d0*x2*y2*z-4.0d0*y2*y2*z)*r6in
          dylmdz22= sq4532pi*(27.0d0*y*z2*r2-3.0d0*y*r4-28.0d0*y*z2*z2)*r6in
          dylmdz23= sq4516pi*(16.0d0*xyz*r2-28.0d0*xyz*z2)*r6in
          dylmdz24= sq31532pi*(y3*r2-3.0d0*x2*y*r2-4.0d0*y3*z2+12.0d0*x2*y*z2)*r6in
          dylmdz25= sq31516pi*(-4.0d0*x3*yz+4.0d0*zx*y3)*r6in
          augbase3d0=rfac(0,ispe)*aug0
          augbase3d1=rfac(1,ispe)*r *rcutin *aug1
          augbase3d2=rfac(2,ispe)*r2*rcut2in*aug2
          augbase3d3=rfac(3,ispe)*r3*rcut3in*aug3
          augbase3d4=rfac(4,ispe)*r4*rcut4in*aug4
          drhoaug3ddz(i,iapsd)=drhoaug3ddz(i,iapsd) &
                         +rhoaugmom( 1,na)*dylmdz01*augbase3d0 &
                         +rhoaugmom( 2,na)*dylmdz02*augbase3d1 &
                         +rhoaugmom( 3,na)*dylmdz03*augbase3d1 &
                         +rhoaugmom( 4,na)*dylmdz04*augbase3d1 &
                         +rhoaugmom( 5,na)*dylmdz05*augbase3d2 &
                         +rhoaugmom( 6,na)*dylmdz06*augbase3d2 &
                         +rhoaugmom( 7,na)*dylmdz07*augbase3d2 &
                         +rhoaugmom( 8,na)*dylmdz08*augbase3d2 &
                         +rhoaugmom( 9,na)*dylmdz09*augbase3d2 &
                         +rhoaugmom(10,na)*dylmdz10*augbase3d3 &
                         +rhoaugmom(11,na)*dylmdz11*augbase3d3 &
                         +rhoaugmom(12,na)*dylmdz12*augbase3d3 &
                         +rhoaugmom(13,na)*dylmdz13*augbase3d3 &
                         +rhoaugmom(14,na)*dylmdz14*augbase3d3 &
                         +rhoaugmom(15,na)*dylmdz15*augbase3d3 &
                         +rhoaugmom(16,na)*dylmdz16*augbase3d3 &
                         +rhoaugmom(17,na)*dylmdz17*augbase3d4 &
                         +rhoaugmom(18,na)*dylmdz18*augbase3d4 &
                         +rhoaugmom(19,na)*dylmdz19*augbase3d4 &
                         +rhoaugmom(20,na)*dylmdz20*augbase3d4 &
                         +rhoaugmom(21,na)*dylmdz21*augbase3d4 &
                         +rhoaugmom(22,na)*dylmdz22*augbase3d4 &
                         +rhoaugmom(23,na)*dylmdz23*augbase3d4 &
                         +rhoaugmom(24,na)*dylmdz24*augbase3d4 &
                         +rhoaugmom(25,na)*dylmdz25*augbase3d4
!        end if
      end do
    end select
  end if
  end do

  return
end subroutine


subroutine augcharge_04( &
 key_natpri_inps,key_pp_paw,                                                                & ! <
 nradmx,lrhomx,lmx,natom,num_spe,num_list_d,num_ppcell_d,                                   & ! <
 indspe,natprid,napsd,natinfd,ndatx,ndaty,ndatz,ntyppp,nradct,lpmx,nwexp,lstdx,lstdy,lstdz, & ! <
 ddx,ddy,ddz,radial,rfac,atx,aty,atz,rhoaugmom,                                             & ! <
 rhoaug3d)                                                                                    ! >
use mod_mpi
implicit none
integer,intent(in) :: key_natpri_inps,key_pp_paw
integer,intent(in) :: nradmx,lrhomx,lmx,natom,num_spe,num_list_d,num_ppcell_d
integer,intent(in) :: indspe(natom),natprid(natom),napsd(natom),natinfd(natom)
integer,intent(in) :: ndatx(natom),ndaty(natom),ndatz(natom)
integer,intent(in) :: ntyppp(num_spe),nradct(num_spe),lpmx(num_spe),nwexp(num_spe)
integer,intent(in) :: lstdx(num_list_d,num_ppcell_d),lstdy(num_list_d,num_ppcell_d),lstdz(num_list_d,num_ppcell_d)
real*8, intent(in) :: ddx,ddy,ddz
real*8, intent(in) :: radial(nradmx,num_spe),rfac(0:2*(lmx-1),num_spe)
real*8, intent(in) :: atx(natom),aty(natom),atz(natom)
real*8, intent(in) :: rhoaugmom(lrhomx,natom)
real*8, intent(out):: rhoaug3d(num_list_d,num_ppcell_d)
real*8 ylm01,ylm02,ylm03,ylm04,ylm05,ylm06,ylm07,ylm08,ylm09,ylm10 &
      ,ylm11,ylm12,ylm13,ylm14,ylm15,ylm16,ylm17,ylm18,ylm19,ylm20 &
      ,ylm21,ylm22,ylm23,ylm24,ylm25
real*8 aug0,aug1,aug2,aug3,aug4
real*8 augbase3d0,augbase3d1,augbase3d2,augbase3d3,augbase3d4
real*8 rcut,rcutin,rcut2in,rcut3in,rcut4in &
      ,r,r2,r3,r4,rin,r2in,r3in,r4in
real*8 x,y,z,x2,y2,z2,xy,yz,zx,xyz
real*8 atmtmpx,atmtmpy,atmtmpz
real*8 pi,sqfourpi,sq3fourpi,sq1516pi,sq0516pi,sq1504pi,sq0716pi,sq2132pi,sq10516pi &
      ,sq3532pi,sq10504pi,sq09256pi,sq4532pi,sq4564pi,sq31532pi,sq315256pi,sq4516pi,sq31516pi
integer na,iapsd,ispe,i,ix,iy,iz

  pi=4.0d0*datan(1.0d0)
  sqfourpi=dsqrt(1.0d0/4.0d0/pi)
  sq3fourpi=dsqrt(3.0d0/4.0d0/pi)
  sq1516pi=dsqrt(15.0d0/16.0d0/pi)
  sq0516pi=dsqrt( 5.0d0/16.0d0/pi)
  sq1504pi=dsqrt(15.0d0/4.0d0/pi)
  sq0716pi=dsqrt(7.0d0/16.0d0/pi)
  sq2132pi=dsqrt(21.0d0/32.0d0/pi)
  sq10516pi=dsqrt(105.0d0/16.0d0/pi)
  sq3532pi=dsqrt(35.0d0/32.0d0/pi)
  sq10504pi=dsqrt(105.0d0/4.0d0/pi)
  sq09256pi=dsqrt(9.0d0/256.0d0/pi)
  sq4532pi=dsqrt(45.0d0/32.0d0/pi)
  sq4564pi=dsqrt(45.0d0/64.0d0/pi)
  sq31532pi=dsqrt(315.0d0/32.0d0/pi)
  sq315256pi=dsqrt(315.0d0/256.0d0/pi)
  sq4516pi=dsqrt(45.0d0/16.0d0/pi)
  sq31516pi=dsqrt(315.0d0/16.0d0/pi)

!$omp do
  do i=1,num_list_d*num_ppcell_d
    rhoaug3d(i,1)=0.0d0
  end do

  do na=1,natom
  if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natprid(na) .eq. key_natpri_inps)) then
    iapsd=napsd(na)
    ispe=indspe(na)
    rcut=radial(nradct(ispe),ispe)
    rcutin=1.0d0/rcut
    rcut2in=rcutin*rcutin
    rcut3in=rcut2in*rcutin
    rcut4in=rcut3in*rcutin
    select case(lpmx(ispe))
    case(0)
!$omp do
      do i=1,natinfd(na)
        ix=lstdx(i,iapsd)
        iy=lstdy(i,iapsd)
        iz=lstdz(i,iapsd)
        atmtmpx=atx(na)-(ndatx(na)*ddx-0.5d0*ddx)
        atmtmpy=aty(na)-(ndaty(na)*ddy-0.5d0*ddy)
        atmtmpz=atz(na)-(ndatz(na)*ddz-0.5d0*ddz)
        x=ix*ddx-atmtmpx
        y=iy*ddy-atmtmpy
        z=iz*ddz-atmtmpz
        x2=x*x
        y2=y*y
        z2=z*z
        r2=x2+y2+z2
        r=dsqrt(r2)
!        if (r .le. radial(nradct(ispe),ispe)) then
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)

          ylm01= sqfourpi
          augbase3d0= rfac(0,ispe)*aug0
          rhoaug3d(i,iapsd)= rhoaugmom( 1,na)*ylm01*augbase3d0
!        end if
      end do
    case(1)
!$omp do
      do i=1,natinfd(na)
        ix=lstdx(i,iapsd)
        iy=lstdy(i,iapsd)
        iz=lstdz(i,iapsd)
        atmtmpx=atx(na)-(ndatx(na)*ddx-0.5d0*ddx)
        atmtmpy=aty(na)-(ndaty(na)*ddy-0.5d0*ddy)
        atmtmpz=atz(na)-(ndatz(na)*ddz-0.5d0*ddz)
        x=ix*ddx-atmtmpx
        y=iy*ddy-atmtmpy
        z=iz*ddz-atmtmpz
        x2=x*x
        y2=y*y
        z2=z*z
        r2=x2+y2+z2
        r=dsqrt(r2)
!        if (r .le. radial(nradct(ispe),ispe)) then
          xy=x*y
          yz=y*z
          zx=z*x
          rin=1.0d0/r
          r2in=rin*rin
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)
          ylm01= sqfourpi
          ylm02= sq3fourpi*x*rin
          ylm03= sq3fourpi*y*rin
          ylm04= sq3fourpi*z*rin
          ylm05= sq1516pi*(x2-y2)*r2in
          ylm06=-sq1504pi*zx*r2in
          ylm07= sq0516pi*(3.0d0*z2-r2)*r2in
          ylm08=-sq1504pi*yz*r2in
          ylm09= sq1504pi*xy*r2in
          augbase3d0= rfac(0,ispe)*aug0
          augbase3d1=rfac(1,ispe)*r *rcutin *aug1
          augbase3d2=rfac(2,ispe)*r2*rcut2in*aug2
          rhoaug3d(i,iapsd)= rhoaugmom( 1,na)*ylm01*augbase3d0 &
                            +rhoaugmom( 2,na)*ylm02*augbase3d1 &
                            +rhoaugmom( 3,na)*ylm03*augbase3d1 &
                            +rhoaugmom( 4,na)*ylm04*augbase3d1 &
                            +rhoaugmom( 5,na)*ylm05*augbase3d2 &
                            +rhoaugmom( 6,na)*ylm06*augbase3d2 &
                            +rhoaugmom( 7,na)*ylm07*augbase3d2 &
                            +rhoaugmom( 8,na)*ylm08*augbase3d2 &
                            +rhoaugmom( 9,na)*ylm09*augbase3d2
!        end if
      end do
    case(2)
!$omp do
      do i=1,natinfd(na)
        ix=lstdx(i,iapsd)
        iy=lstdy(i,iapsd)
        iz=lstdz(i,iapsd)
        atmtmpx=atx(na)-(ndatx(na)*ddx-0.5d0*ddx)
        atmtmpy=aty(na)-(ndaty(na)*ddy-0.5d0*ddy)
        atmtmpz=atz(na)-(ndatz(na)*ddz-0.5d0*ddz)
        x=ix*ddx-atmtmpx
        y=iy*ddy-atmtmpy
        z=iz*ddz-atmtmpz
        x2=x*x
        y2=y*y
        z2=z*z
        r2=x2+y2+z2
        r=dsqrt(r2)
!        if (r .le. radial(nradct(ispe),ispe)) then
          xy=x*y
          yz=y*z
          zx=z*x
          xyz=xy*z
          r3=r2*r
          r4=r3*r
          rin=1.0d0/r
          r2in=rin*rin
          r3in=r2in*rin
          r4in=r3in*rin
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)
          ylm01= sqfourpi
          ylm02= sq3fourpi*x*rin
          ylm03= sq3fourpi*y*rin
          ylm04= sq3fourpi*z*rin
          ylm05= sq1516pi*(x2-y2)*r2in
          ylm06=-sq1504pi*zx*r2in
          ylm07= sq0516pi*(3.0d0*z2-r2)*r2in
          ylm08=-sq1504pi*yz*r2in
          ylm09= sq1504pi*xy*r2in
          ylm10= sq0716pi*z*(5.0d0*z2-3.0d0*r2)*r3in
          ylm11=-sq2132pi*x*(5.0d0*z2-r2)*r3in
          ylm12= sq10516pi*z*(x2-y2)*r3in
          ylm13=-sq3532pi*x*(x2-3.0d0*y2)*r3in
          ylm14=-sq2132pi*y*(5.0d0*z2-r2)*r3in
          ylm15= sq10504pi*xyz*r3in
          ylm16=-sq3532pi*y*(3.0d0*x2-y2)*r3in
          ylm17= sq09256pi*(35.0d0*z2*z2-30.0d0*z2*r2+3.0d0*r2*r2)*r4in
          ylm18= sq4532pi*x*z*(7.0d0*z2-3.0d0*r2)*r4in
          ylm19= sq4564pi*(x2-y2)*(7.0d0*z2-r2)*r4in
          ylm20= sq31532pi*zx*(x2-3.0d0*y2)*r4in
          ylm21= sq315256pi*(x2*x2-6.0d0*x2*y2+y2*y2)*r4in
          ylm22= sq4532pi*yz*(7.0d0*z2-3.0d0*r2)*r4in
          ylm23= sq4516pi*xy*(7.0d0*z2-r2)*r4in
          ylm24= sq31532pi*yz*(y2-3.0d0*x2)*r4in
          ylm25= sq31516pi*xy*(x2-y2)*r4in
          augbase3d0= rfac(0,ispe)*aug0
          augbase3d1=rfac(1,ispe)*r *rcutin *aug1
          augbase3d2=rfac(2,ispe)*r2*rcut2in*aug2
          augbase3d3=rfac(3,ispe)*r3*rcut3in*aug3
          augbase3d4=rfac(4,ispe)*r4*rcut4in*aug4
          rhoaug3d(i,iapsd)= rhoaugmom( 1,na)*ylm01*augbase3d0 &
                            +rhoaugmom( 2,na)*ylm02*augbase3d1 &
                            +rhoaugmom( 3,na)*ylm03*augbase3d1 &
                            +rhoaugmom( 4,na)*ylm04*augbase3d1 &
                            +rhoaugmom( 5,na)*ylm05*augbase3d2 &
                            +rhoaugmom( 6,na)*ylm06*augbase3d2 &
                            +rhoaugmom( 7,na)*ylm07*augbase3d2 &
                            +rhoaugmom( 8,na)*ylm08*augbase3d2 &
                            +rhoaugmom( 9,na)*ylm09*augbase3d2 &
                            +rhoaugmom(10,na)*ylm10*augbase3d3 &
                            +rhoaugmom(11,na)*ylm11*augbase3d3 &
                            +rhoaugmom(12,na)*ylm12*augbase3d3 &
                            +rhoaugmom(13,na)*ylm13*augbase3d3 &
                            +rhoaugmom(14,na)*ylm14*augbase3d3 &
                            +rhoaugmom(15,na)*ylm15*augbase3d3 &
                            +rhoaugmom(16,na)*ylm16*augbase3d3 &
                            +rhoaugmom(17,na)*ylm17*augbase3d4 &
                            +rhoaugmom(18,na)*ylm18*augbase3d4 &
                            +rhoaugmom(19,na)*ylm19*augbase3d4 &
                            +rhoaugmom(20,na)*ylm20*augbase3d4 &
                            +rhoaugmom(21,na)*ylm21*augbase3d4 &
                            +rhoaugmom(22,na)*ylm22*augbase3d4 &
                            +rhoaugmom(23,na)*ylm23*augbase3d4 &
                            +rhoaugmom(24,na)*ylm24*augbase3d4 &
                            +rhoaugmom(25,na)*ylm25*augbase3d4
!        end if
      end do
    end select
  end if
  end do

  return
end subroutine


end module
