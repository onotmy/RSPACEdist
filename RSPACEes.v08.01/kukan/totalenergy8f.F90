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
! **********  totalenergy8f.F90 04/27/2023-01  **********

module mod_totalenergy
implicit none
contains


subroutine totalenergy_space( &
 key_natpri_in,key_pp_paw,key_jel_calc,                                       & ! <
 key_polcon_atoms,key_polcon_asa,key_polcon2_size,                            & ! <
 idim,natom,num_spe,nperi,nmesh,nspv,npoint,ncpx,ncpy,ncpz,num_atcell,nradmx, & ! <
 jelcalc,nint1dmax,nzmax,new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz,     & ! <
 natpri,natpri_inf,indspe,ntyppp,nradct,npolcon,npolcon2,                     & ! <
 wt,radial,dradial,cp,veta,tnumele,dx,dy,dz,xmax,ymax,zmax,                   & ! <
 biasx,biasy,biasz,chrjel,endjel,strjel,atx,aty,atz,                          & ! <
 rhosmt,rho_aug_dense,rhosmt_pcc_dense,                                       & ! <
 rhotrur,rhosmtr,rhoaugr,rhosmt_pccr,rhotrucorer,spinpol,polconb,             & ! <
 vh_dense,vx,vhtrur,vhsmtr,vhaugr,vxctru,vxcsmt,ex_dense,exctru,excsmt,       & ! <
 eneco,eneha,eneex,eneat,eneof,enejel,enebc)                                    ! >
use mod_mpi
implicit none
integer,intent(in) ::key_natpri_in,key_pp_paw,key_jel_calc
integer,intent(in) ::key_polcon_atoms,key_polcon_asa,key_polcon2_size
integer,intent(in) ::idim,natom,num_spe,nperi,nmesh,nspv,npoint,ncpx,ncpy,ncpz,num_atcell,nradmx
integer,intent(in) ::jelcalc,nint1dmax,nzmax,new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz
integer,intent(in) ::natpri(natom),natpri_inf(natom),indspe(natom),ntyppp(num_spe),nradct(num_spe)
integer,intent(in) ::npolcon,npolcon2(0,idim)
real*8, intent(in) ::wt(npoint),radial(nradmx,num_spe),dradial(nradmx,num_spe),cp(8,num_spe)
real*8, intent(in) ::veta,tnumele
real*8, intent(in) ::dx,dy,dz,xmax,ymax,zmax,biasx,biasy,biasz,chrjel,endjel,strjel
real*8, intent(in) ::atx(natom),aty(natom),atz(natom)
real*8, intent(in) ::rhosmt(ncpx,ncpy,ncpz,nspv)
real*8, intent(in) ::rho_aug_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
real*8, intent(in) ::rhosmt_pcc_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh,nspv)
real*8, intent(in) ::vh_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
real*8, intent(in) ::vx(ncpx,ncpy,ncpz,nspv)
real*8, intent(in) ::ex_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
real*8, intent(in) ::rhotrur(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::rhosmtr(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::rhoaugr(nradmx,npoint,num_atcell)
real*8, intent(in) ::rhosmt_pccr(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::rhotrucorer(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::spinpol(max(1,nspv-1),0:min(1,nspv-1)*natom) 
real*8, intent(in) ::polconb(3,idim)
real*8, intent(in) ::vhtrur(nradmx,npoint,num_atcell)
real*8, intent(in) ::vhsmtr(nradmx,npoint,num_atcell)
real*8, intent(in) ::vhaugr(nradmx,npoint,num_atcell)
real*8, intent(in) ::vxctru(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::vxcsmt(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::exctru(nradmx,npoint,num_atcell)
real*8, intent(in) ::excsmt(nradmx,npoint,num_atcell)
real*8, intent(out)::eneco,eneha,eneex,eneat,eneof,enejel,enebc
integer :: na
real*8  :: enehaa,eneexa
real*8  :: eneall

!  ==========  compensate double counting  ==========
!  -H1 : the 1st term of E_{dc}^1 in Eq. (48) of PRB59 1758 (1999)
!   H2 : the 1st term of \tilde{E}_{dc}^1 in Eq. (48) of PRB59 1758 (1999)
!  -EX1: the 3rd term of E_{dc}^1 in Eq. (48) of PRB59 1758 (1999)
!   EX2: the 2nd term of E_{dc}^1 in Eq. (48) of PRB59 1758 (1999)
!   EX3: the 3rd term of \tilda{E}_{dc}^1 in Eq. (48) of PRB59 1758 (1999)
!  -EX4: the 2nd term of \tilda{E}_{dc}^1 in Eq. (48) of PRB59 1758 (1999)
  !$omp parallel default(shared)
  call totalenergy_space_01( &
   nspv,nmesh,ncpx,ncpy,ncpz,dx,dy,dz,                         & ! <
   rhosmt,rho_aug_dense,rhosmt_pcc_dense,vh_dense,vx,ex_dense, & ! <
   eneha,eneex)                                                  ! >
  !$omp end parallel
  !$omp parallel default(shared)
  call totalenergy_space_02( &
   key_natpri_in,key_pp_paw,                                                                         & ! <
   natom,num_spe,num_atcell,nspv,npoint,nradmx,                                                      & ! <
   indspe,natpri,natpri_inf,ntyppp,nradct,wt,radial,dradial,                                         & ! <
   rhotrur,rhosmtr,rhoaugr,rhosmt_pccr,rhotrucorer,vhtrur,vhsmtr,vhaugr,vxctru,vxcsmt,exctru,excsmt, & ! <
   enehaa,eneexa)                                                                                      ! >
  !$omp end parallel
  eneha=eneha+enehaa
  eneex=eneex+eneexa
! ==================================================

! ==========  ewald energy  ==========
  !$omp parallel default(shared)
  call totalenergy_space_03( &
   natom,num_spe,nperi,nint1dmax,indspe,new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz,xmax,ymax,zmax,veta,cp, & ! <
   atx,aty,atz,                                                                                       & ! <
   eneat)                                                                                               ! >
  !$omp end parallel
! ====================================

! ==========  field correction  ==========
  eneco=0.0d0
  do na=1,natom
    eneco=eneco-(biasx*atx(na)+biasy*aty(na)+biasz*atz(na))*cp(1,indspe(na))
  end do
! ========================================

! ==========  energy offset of pseudopotential  ==========
  !$omp parallel default(shared)
  call totalenergy_space_05( &
   natom,num_spe,nperi,npoint,num_atcell,nradmx, & ! <
   key_pp_paw,key_natpri_in,                     & ! <
   tnumele,xmax,ymax,zmax,                       & ! <
   natpri,natpri_inf,indspe,ntyppp,nradct,       & ! <
   wt,radial,dradial,vhsmtr,vhaugr,vhtrur,       & ! <
   eneof)                                          ! >
  !$omp end parallel
! ========================================================

! ==========  energy between jellium and ion  ==========
  !$omp parallel default(shared)
  call totalenergy_space_04( &
   jelcalc,nperi,natom,num_spe,indspe,ncpx,ncpy,ncpz,nmesh,nzmax,key_jel_calc, & ! <
   xmax,ymax,zmax,dx,dy,dz,cp,chrjel,endjel,strjel,                            & ! <
   atx,aty,atz,                                                                & ! <
   enejel)                                                                       ! >
  !$omp end parallel
! ======================================================

! ==========  energy correction due to parallel constraining B-field  ==========
  if (nspv>1) then 
    call totalenergy_space_06( &
    key_polcon_atoms,key_polcon_asa,key_polcon2_size,key_natpri_in, & ! <
    natom,nspv,npolcon,npolcon2,natpri,                             & ! <
    spinpol,polconb,                                                & ! <
    enebc)                                                            ! >
  else
    enebc= 0.0d0
  endif
! =============================================================================


    eneall= eneha
    call mpi_reduce(eneall,eneha,1,mpi_double_precision,mpi_sum,0,mpicom_space,mpij)
    eneall= eneex
    call mpi_reduce(eneall,eneex,1,mpi_double_precision,mpi_sum,0,mpicom_space,mpij)
    eneall= eneof
    call mpi_reduce(eneall,eneof,1,mpi_double_precision,mpi_sum,0,mpicom_space,mpij)
    if (jelcalc==key_jel_calc) then
      eneall= enejel
      call mpi_reduce(eneall,enejel,1,mpi_double_precision,mpi_sum,0,mpicom_space,mpij)
    end if
    if ((npolcon==key_polcon_atoms).or.(npolcon==key_polcon_asa)) then
      eneall= enebc
      call mpi_reduce(eneall,enebc,1,mpi_double_precision,mpi_sum,0,mpicom_space,mpij)
    endif 

end subroutine totalenergy_space


subroutine totalenergy_eig( &
 nums,ncol,numk,neigmx,nwskk,nwskptot,tf,sval,fnele, & ! <
 eneel,eneth)                                        ! >
use mod_mpi
implicit none
integer, intent(in) :: nums,ncol,numk,neigmx
integer, intent(in) :: nwskk(numk),nwskptot
real*8,  intent(in) :: tf, sval(neigmx,nums+1-ncol,numk), fnele(neigmx,nums+1-ncol,numk)
real*8,  intent(out):: eneel, eneth
integer :: nk,ns,l
real*8  :: tmpfac,tmpele,tmp1,tmp2
real*8  :: eneall

! ==========  band-strucure energy  ==========
  eneel=0.0d0
  do nk=1,numk
  do ns=1,nums+1-ncol
  do l=1,neigmx
     eneel=eneel+fnele(l,ns,nk)*sval(l,ns,nk)
  end do
  end do
  end do
! ============================================

! ==========  helmholtz free energy  ==========
  eneth=0.0d0
  do nk=1,numk
    tmpfac= dble(nwskptot*nums)/dble(2*nwskk(nk))
    do ns=1,nums+1-ncol
    do l=1,neigmx
      tmpele=fnele(l,ns,nk)*tmpfac
      tmp1=0.0d0
      tmp2=0.0d0
      if (tmpele .gt. 1.0d-8) tmp1=dlog(tmpele)
      if (1.0d0-tmpele .gt. 1.0d-8) tmp2=dlog(1.0d0-tmpele)
      eneth=eneth+(tmpele*tmp1+(1.0d0-tmpele)*tmp2)/tmpfac
    end do
    end do
  end do
! =============================================

    eneall= eneel
    call mpi_reduce(eneall,eneel,1,mpi_double_precision,mpi_sum,0,mpicom_kpt,mpij)
    eneall= eneth
    call mpi_reduce(eneall,eneth,1,mpi_double_precision,mpi_sum,0,mpicom_kpt,mpij)

  if (myrank_glbl==0) eneth=eneth*tf

end subroutine totalenergy_eig

! --------------------------------------------------------------------------------------------------

subroutine totalenergy_space_01( &
 nspv,nmesh,ncpx,ncpy,ncpz,dx,dy,dz,                         & ! <
 rhosmt,rho_aug_dense,rhosmt_pcc_dense,vh_dense,vx,ex_dense, & ! <
 enehaa,eneexa)                                                ! >
! -H2 : the 1st term of \tilda{E}_{dc} in Eq. (48) of PRB59 1758 (1999)
! -EX3: the 3rd term of \tilda{E}_{dc} in Eq. (48) of PRB59 1758 (1999)
!  EX4: the 2nd term of \tilda{E}_{dc} in Eq. (48) of PRB59 1758 (1999)
implicit none
integer,intent(in) ::nspv,nmesh,ncpx,ncpy,ncpz
real*8, intent(in) ::dx,dy,dz
real*8, intent(in) ::rhosmt(ncpx,ncpy,ncpz,nspv)
real*8, intent(in) ::rho_aug_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
real*8, intent(in) ::rhosmt_pcc_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh,nspv)
real*8, intent(in) ::vh_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
real*8, intent(in) ::vx(ncpx,ncpy,ncpz,nspv)
real*8, intent(in) ::ex_dense(ncpx*nmesh,ncpy*nmesh,ncpz*nmesh)
real*8, intent(out)::enehaa,eneexa
integer ncpx_d,ncpy_d,ncpz_d
real*8 ddx,ddy,ddz
integer ns,ix
real*8 eneha,eneex,dxyz,ddxyz

  ncpx_d=ncpx*nmesh
  ncpy_d=ncpy*nmesh
  ncpz_d=ncpz*nmesh
  ddx=dx/dble(nmesh)
  ddy=dy/dble(nmesh)
  ddz=dz/dble(nmesh)
  dxyz=dx*dy*dz
  ddxyz=ddx*ddy*ddz

  !$omp single
  enehaa=0.0d0
  eneexa=0.0d0
  !$omp end single
  !$omp barrier

  eneha=0.0d0
  !$omp do
  do ix=1,ncpx_d*ncpy_d*ncpz_d
!     ----------  subtract derivative of EH of smooth charge (-H2)  ----------
    eneha=eneha-0.5d0*vh_dense(ix,1,1)*rho_aug_dense(ix,1,1)*ddxyz
  end do
  !$omp critical
  enehaa=enehaa+eneha
  !$omp end critical
  !$omp barrier

  do ns=1,nspv
    eneex=0.0d0
    !$omp do
    do ix=1,ncpx*ncpy*ncpz
!     ----------  subtract derivative of EXC of smooth charge (-EX3)  ----------
      eneex=eneex-rhosmt(ix,1,1,ns)*vx(ix,1,1,ns)*dxyz
    end do
    if (ns<3) then
      !$omp do
      do ix=1,ncpx_d*ncpy_d*ncpz_d
!     ----------  add EXC of smooth charge (EX4)  ----------
        eneex=eneex+rhosmt_pcc_dense(ix,1,1,ns)*ex_dense(ix,1,1)*ddxyz
      end do
    else
      eneex= eneex*2.0d0
    end if
    !$omp critical
    eneexa=eneexa+eneex
    !$omp end critical
    !$omp barrier
  end do

end subroutine totalenergy_space_01


subroutine totalenergy_space_02( &
 key_natpri_in,key_pp_paw,                                                                         & ! <
 natom,num_spe,num_atcell,nspv,npoint,nradmx,                                                      & ! <
 indspe,natpri,natpri_inf,ntyppp,nradct,wt,radial,dradial,                                         & ! <
 rhotrur,rhosmtr,rhoaugr,rhosmt_pccr,rhotrucorer,vhtrur,vhsmtr,vhaugr,vxctru,vxcsmt,exctru,excsmt, & ! <
 enehaa,eneexa)                                                                                      ! >
implicit none
integer,intent(in) ::key_natpri_in,key_pp_paw
integer,intent(in) ::natom,num_spe,num_atcell,nspv,npoint,nradmx
integer,intent(in) ::indspe(natom),natpri(natom),natpri_inf(natom),ntyppp(num_spe),nradct(num_spe)
real*8, intent(in) ::wt(npoint),radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8, intent(in) ::rhotrur(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::rhosmtr(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::rhoaugr(nradmx,npoint,num_atcell)
real*8, intent(in) ::rhosmt_pccr(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::rhotrucorer(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::vhtrur(nradmx,npoint,num_atcell)
real*8, intent(in) ::vhsmtr(nradmx,npoint,num_atcell)
real*8, intent(in) ::vhaugr(nradmx,npoint,num_atcell)
real*8, intent(in) ::vxctru(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::vxcsmt(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::exctru(nradmx,npoint,num_atcell)
real*8, intent(in) ::excsmt(nradmx,npoint,num_atcell)
real*8, intent(out)::enehaa,eneexa
integer na,ipri,ns,il,ir
real*8 eneha0,eneex0,eneex,r,dr

  !$omp single
  enehaa=0.0d0
  eneexa=0.0d0
  !$omp end single
  !$omp barrier

  eneha0=0.0d0
  eneex0=0.0d0
  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na)==key_natpri_in)) then
      ipri=natpri_inf(na)
      do ns=1,min(2,nspv)
        do il=1,npoint
          !$omp do
          do ir=2,nradct(indspe(na))-1
            r=radial(ir,indspe(na))
            dr=dradial(ir,indspe(na))
!     ----------  subtract derivative of EH of true charge (-H1)  ----------
            eneha0=eneha0-0.5d0*(vhtrur(ir,il,ipri)                  ) &
                  *(rhotrur(ir,il,ns,ipri)                     )*r*r*dr*wt(il)
!     ----------  add derivative of EH of smooth charge (H2)  ----------
            eneha0=eneha0+0.5d0*(vhsmtr(ir,il,ipri)+ vhaugr(ir,il,ipri)) &
                  *(rhosmtr(ir,il,ns,ipri)+rhoaugr(ir,il,ipri)/dble(min(2,nspv)) )*r*r*dr*wt(il)
          end do
        end do
      end do
      do ns=1,nspv
        do il=1,npoint
          !$omp do
          do ir=2,nradct(indspe(na))-1
            eneex= 0.0d0
!     ----------  subtract derivative of EXC of true charge (-EX1)  ----------
            eneex=eneex-vxctru(ir,il,ns,ipri)*rhotrur(ir,il,ns,ipri)
!     ----------  add derivative of EXC of smooth charge (EX3)  ----------
            eneex=eneex+vxcsmt(ir,il,ns,ipri)*rhosmtr(ir,il,ns,ipri)
            if (ns<3) then
!     ----------  add EXC of true charge (EX2)  ----------
              eneex=eneex+exctru(ir,il,ipri)*rhotrucorer(ir,il,ns,ipri)
!     ----------  subtract EXC of smooth charge (-EX4)  ----------
              eneex=eneex-excsmt(ir,il,ipri)*rhosmt_pccr(ir,il,ns,ipri)
            else
              eneex= eneex*2.0d0
            end if
            r=radial(ir,indspe(na))
            dr=dradial(ir,indspe(na))
            eneex0= eneex0 + eneex*r*r*dr*wt(il)
          end do
        end do
      end do
    end if
  end do
  !$omp critical
  enehaa=enehaa+eneha0
  eneexa=eneexa+eneex0
  !$omp end critical
  !$omp barrier

end subroutine totalenergy_space_02


subroutine totalenergy_space_03( &
 natom,num_spe,nperi,nint1dmax,indspe,new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz,xmax,ymax,zmax,veta,cp, & ! <
 atx,aty,atz,                                                                                       & ! <
 eneat)                                                                                               ! >
use mod_mpi
use mod_mathfunctions, only: expint1
implicit none
integer,intent(in) :: natom,num_spe,nperi,nint1dmax,indspe(natom)
integer,intent(in) :: new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz
real*8, intent(in) :: xmax,ymax,zmax
real*8, intent(in) :: veta
real*8, intent(in) :: cp(8,num_spe)
real*8, intent(in) :: atx(natom),aty(natom),atz(natom)
real*8, intent(out):: eneat
real*8  :: pi,fourpi,omega,omegain,surf,surfin,eneat0
real*8  :: x,y,z,x0,y0,z0,x1,y1,z1,r
real*8  :: vlo0,vlo1,vlo2,vlo3,vep,delta
real*8  :: vkx,vky,vkz,vk2,ta,tb,t,dt
integer :: na,na1,l0,l1,nas,nae,i,natpprcs
integer :: kx,ky,kz,ierr

real*8 derf

  pi=dacos(-1.0d0)
  fourpi=4.0d0*pi
  omega=8.0d0*xmax*ymax*zmax
  omegain=1.0d0/omega
  surf=4.0d0*xmax*ymax
  surfin=1.0d0/surf

  !$omp single
  eneat=0.0d0
  !$omp end single
  !$omp barrier

  nas=1
  do i= 0,myrank_glbl-1
    natpprcs=natom/nprocs
    if ( i >= nprocs-mod(natom,nprocs) ) natpprcs=natom/nprocs+1
    nas=nas+natpprcs
  end do
  natpprcs=natom/nprocs
  if ( myrank_glbl >= nprocs-mod(natom,nprocs) ) natpprcs=natom/nprocs+1
  nae=nas+natpprcs-1

  eneat0=0.0d0
  select case(nperi)
  case (0)
    do na=nas,nae
      x0=atx(na)
      y0=aty(na)
      z0=atz(na)
      !$omp do
      do na1=1,natom
        x1=atx(na1)
        y1=aty(na1)
        z1=atz(na1)
        x=x0-x1
        y=y0-y1
        z=z0-z1
        r=dsqrt(x*x+y*y+z*z)
        if (r .gt. 1.0d-20) then
          eneat0=eneat0+0.5d0*cp(1,indspe(na))*cp(1,indspe(na1))/r
        end if
      end do
    end do
  case (1)
    do l0=nas,nae
      x0=atx(l0)
      y0=aty(l0)
      z0=atz(l0)
      !$omp do
      do l1=1,natom
        x1=atx(l1)
        y1=aty(l1)
        z1=atz(l1)
        do kx=-new_pwx,new_pwx
          if (kx**2 .ne. 0) then
          vkx=pi/xmax*kx
          ta=((y0-y1)*(y0-y1)+(z0-z1)*(z0-z1))
          tb=vkx*vkx
          dt=veta/nint1dmax
          vlo0=0.0d0
          do i=1,nint1dmax
            t=i*dt
            vlo0=vlo0+dexp(-ta*t*t-0.25d0*tb/(t*t))/t*dt
          end do
          vlo0=0.5d0/xmax*dcos(vkx*(x0-x1))*vlo0
          eneat0=eneat0+cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0
          end if
        end do
        ta=((y0-y1)*(y0-y1)+(z0-z1)*(z0-z1))
        delta=0.0d0
        if (l0 .eq. l1) delta=1.0d0
        if (ta .lt. 1.0d-16) then
          vlo1=0.5d0/xmax*(0.5d0*0.577215664901532d0+0.5d0*dlog(veta**2))
        else
          vlo1=0.5d0/xmax*(-0.5d0*expint1(ta*veta*veta,ierr)-0.5d0*dlog(ta))
        end if
        vlo2=-delta*veta/dsqrt(pi)
        eneat0=eneat0+cp(1,indspe(l0))*cp(1,indspe(l1))*(vlo1+vlo2)
        do kx=-new_rsx,new_rsx
          x=kx*2.0d0*xmax+x0-x1
          y=y0-y1
          z=z0-z1
          r=dsqrt(x*x+y*y+z*z)
          if (r .gt. 1.0d-20) then
            vlo0=(1.0d0-derf(veta*r))/r
            eneat0=eneat0+0.5d0*cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0
          end if
        end do
      end do
    end do
  case (2)
    do l0=nas,nae
      x0=atx(l0)
      y0=aty(l0)
      z0=atz(l0)
      !$omp do
      do l1=1,natom
        x1=atx(l1)
        y1=aty(l1)
        z1=atz(l1)
        do ky=-new_pwy,new_pwy
        do kx=-new_pwx,new_pwx
          if (kx**2+ky**2 .ne. 0) then
          vkx=pi/xmax*kx
          vky=pi/ymax*ky
          vk2=vkx**2+vky**2
          ta=(z0-z1)
          tb=dsqrt(vk2)
          vep=pi/tb/surf*((1.0d0-derf((tb-2.0d0*veta*veta*ta)/(2.0d0*veta)))/dexp(tb*ta) &
                         +(1.0d0-derf((tb+2.0d0*veta*veta*ta)/(2.0d0*veta)))/dexp(-tb*ta))
          vlo0=dcos(vkx*(x0-x1)+vky*(y0-y1))*vep
          eneat0=eneat0+0.5d0*cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0
          end if
        end do
        end do
        delta=0.0d0
        if (l0 .eq. l1) delta=1.0d0
        z=dabs(atz(l0)-atz(l1))
        vlo1=-2.0d0*dsqrt(pi)/surf*(-dsqrt(pi)*z+dexp(-z**2*veta**2)/veta+dsqrt(pi)*z*derf(z*veta))
        vlo2=-delta*2.0d0*veta/dsqrt(pi)
        vlo3= 2.0d0*dsqrt(pi)/surf*(-dsqrt(pi)*z)
        eneat0=eneat0+0.5d0*cp(1,indspe(l0))*cp(1,indspe(l1))*(vlo1+vlo2+vlo3)
        do ky=-new_rsy,new_rsy
        do kx=-new_rsx,new_rsx
          x=kx*2.0d0*xmax+x0-atx(l1)
          y=ky*2.0d0*ymax+y0-aty(l1)
          z=z0-atz(l1)
          r=dsqrt(x*x+y*y+z*z)
          if (r .gt. 1.0d-20) then
            vlo0=(1.0d0-derf(veta*r))/r
            eneat0=eneat0+0.5d0*cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0
          end if
        end do
        end do
      end do
    end do
  case (3)
    do l0=nas,nae
      x0=atx(l0)
      y0=aty(l0)
      z0=atz(l0)
      !$omp do
      do l1=1,natom
        x1=atx(l1)
        y1=aty(l1)
        z1=atz(l1)
        do kz=-new_pwz,new_pwz
        do ky=-new_pwy,new_pwy
        do kx=-new_pwx,new_pwx
          if (kx**2+ky**2+kz**2 .ne. 0) then
            vkx=pi/xmax*kx
            vky=pi/ymax*ky
            vkz=pi/zmax*kz
            vk2=vkx**2+vky**2+vkz**2
            vlo0=fourpi*omegain*dcos(vkx*(x0-x1)+vky*(y0-y1)+vkz*(z0-z1))/vk2*dexp(-vk2/(4.0d0*veta**2))
            eneat0=eneat0+0.5d0*cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0
          end if
        end do
        end do
        end do
        delta=0.0d0
        if (l0 .eq. l1) delta=1.0d0
        vlo1=-pi/omega/veta**2
        vlo2=-delta*2.0d0*veta/dsqrt(pi)
        eneat0=eneat0+0.5d0*cp(1,indspe(l0))*cp(1,indspe(l1))*(vlo1+vlo2)
        do kz=-new_rsz,new_rsz
        do ky=-new_rsy,new_rsy
        do kx=-new_rsx,new_rsx
          x=kx*2.0d0*xmax+x0-atx(l1)
          y=ky*2.0d0*ymax+y0-aty(l1)
          z=kz*2.0d0*zmax+z0-atz(l1)
          r=dsqrt(x*x+y*y+z*z)
          if (r .gt. 1.0d-20) then
            vlo0=(1.0d0-derf(veta*r))/r
            eneat0=eneat0+0.5d0*cp(1,indspe(l0))*cp(1,indspe(l1))*vlo0
          end if
        end do
        end do
        end do
      end do
    end do
  end select
  !$omp critical
  eneat=eneat+eneat0
  !$omp end critical
  !$omp barrier
  !$omp single
    call mpi_allreduce(eneat,eneat0,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
    eneat=eneat0
  !$omp end single
  !$omp barrier

end subroutine totalenergy_space_03

subroutine totalenergy_space_04( &
 jelcalc,nperi,natom,num_spe,indspe,ncpx,ncpy,ncpz,nmesh,nzmax,key_jel_calc, & ! <
 xmax,ymax,zmax,dx,dy,dz,cp,chrjel,endjel,strjel,                            & ! <
 atx,aty,atz,                                                                & ! <
 enejel)                                                                       ! >
use mod_mpi,       only: myrx,myry,myrz
implicit none
integer,intent(in)  :: jelcalc,nperi
integer,intent(in)  :: natom,num_spe,indspe(natom),ncpx,ncpy,ncpz,nmesh,nzmax,key_jel_calc
real*8, intent(in)  :: xmax,ymax,zmax,dx,dy,dz,cp(8,num_spe),chrjel,endjel,strjel
real*8, intent(in)  :: atx(natom),aty(natom),atz(natom)
real*8, intent(out) :: enejel
integer :: na,kz,ix,iy,iz
integer :: ncpx_d,ncpy_d,ncpz_d
real*8  :: pi,ddx,ddy,ddz,omega
real*8  :: cp1,cp2,cp3,cp4,cp5,vkz,vk2,x,y,z,zz,r,enejel0

  pi=dacos(-1.0d0)
  ncpx_d=ncpx*nmesh
  ncpy_d=ncpy*nmesh
  ncpz_d=ncpz*nmesh
  !$omp single
  enejel= 0.0d0
  !$omp end single
  !$omp barrier

  if (jelcalc==key_jel_calc) then
  if (nperi==3) then
    ddx= dx/dble(nmesh)
    ddy= dy/dble(nmesh)
    ddz= dz/dble(nmesh)
    omega= 8.0d0*xmax*ymax*zmax
    do na=1,natom
      cp1=cp(1,indspe(na))
      cp2=cp(2,indspe(na))
      cp3=cp(3,indspe(na))
      cp4=cp(4,indspe(na))
      cp5=cp(5,indspe(na))
      do kz=-nzmax,nzmax
        if (kz**2 .ne. 0) then
          vkz=pi/zmax*kz
          vk2=vkz*vkz
          enejel0=0.0d0
          !$omp do
          do iz=1,ncpz_d
          do iy=1,ncpy_d
          do ix=1,ncpx_d
            x=dble(ix+myrx*ncpx_d)*ddx-0.5d0*ddx-xmax-atx(na)
            y=dble(iy+myry*ncpy_d)*ddy-0.5d0*ddy-ymax-aty(na)
            z=dble(iz+myrz*ncpz_d)*ddz-0.5d0*ddz-zmax-atz(na)
            zz=dble(iz+myrz*ncpz_d)*ddz-0.5d0*ddz-zmax
            if (x .gt.  xmax) x=x-2.0d0*xmax
            if (x .le. -xmax) x=x+2.0d0*xmax
            if (y .gt.  ymax) y=y-2.0d0*ymax
            if (y .le. -ymax) y=y+2.0d0*ymax
            if (z .gt.  zmax) z=z-2.0d0*zmax
            if (z .le. -zmax) z=z+2.0d0*zmax
            r=dsqrt(x*x+y*y+z*z)
            enejel0=enejel0+4.0d0*pi/omega*chrjel/(strjel-endjel)/vk2 &
              *(dsin(vkz*(strjel-zz))-dsin(vkz*(endjel-zz)))/vkz &
              *cp1*(cp4*(cp2/pi)**1.5d0*dexp(-cp2*r*r)+cp5*(cp3/pi)**1.5d0*dexp(-cp3*r*r))*ddx*ddy*ddz
          end do
          end do
          end do
          !$omp critical
          enejel=enejel+enejel0
          !$omp end critical
          !$omp barrier
        end if
      end do
    end do
  end if ! (nperi==3)
  end if

end subroutine totalenergy_space_04


subroutine totalenergy_space_05( &
 natom,num_spe,nperi,npoint,num_atcell,nradmx, & ! <
 key_pp_paw,key_natpri_in,                     & ! <
 tnumele,xmax,ymax,zmax,                       & ! <
 natpri,natpri_inf,indspe,ntyppp,nradct,       & ! <
 wt,radial,dradial,vhsmtr,vhaugr,vhtrur,       & ! <
 eneof)                                          ! >
implicit none
integer,intent(in)::natom,num_spe,nperi,npoint,num_atcell,nradmx
integer,intent(in)::key_pp_paw,key_natpri_in
real*8, intent(in)::tnumele,xmax,ymax,zmax
integer,intent(in) ::natpri(natom),natpri_inf(natom),indspe(natom),ntyppp(num_spe),nradct(num_spe)
real*8, intent(in) ::wt(npoint),radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8, intent(in) ::vhtrur(nradmx,npoint,num_atcell)
real*8, intent(in) ::vhsmtr(nradmx,npoint,num_atcell)
real*8, intent(in) ::vhaugr(nradmx,npoint,num_atcell)
real*8, intent(out)::eneof
integer il,ir,ipri,na
real*8 omega,omegain,r,dr,eneof0

  !$omp single
  eneof=0.0d0
  !$omp end single
  !$omp barrier

  eneof0=0.0d0
  if (nperi .eq. 3) then
    omega=8.0d0*xmax*ymax*zmax
    omegain=1.0d0/omega
    do na=1,natom
      if ((ntyppp(indspe(na))==key_pp_paw) .and. (natpri(na)==key_natpri_in)) then
        ipri=natpri_inf(na)
        !$omp do
        do il=1,npoint
          do ir=2,nradct(indspe(na))-1
            r=radial(ir,indspe(na))
            dr=dradial(ir,indspe(na))
            eneof0=eneof0+(vhsmtr(ir,il,ipri)+vhaugr(ir,il,ipri)-vhtrur(ir,il,ipri))*r*r*dr*wt(il)
          end do
        end do
      end if
    end do
    !$omp critical
    eneof=eneof+0.5d0*eneof0*omegain*tnumele
    !$omp end critical
    !$omp barrier
  end if

end subroutine totalenergy_space_05


subroutine totalenergy_space_06( &
 key_polcon_atoms,key_polcon_asa,key_polcon2_size,key_natpri_in, & ! <
 natom,nspv,npolcon,npolcon2,natpri,                             & ! <
 spinpol,polconb,                                                & ! < 
 enebc)                                                            ! >
implicit none
integer,intent(in) :: key_polcon_atoms,key_polcon_asa,key_polcon2_size,key_natpri_in
integer,intent(in) :: natom,nspv
integer,intent(in) :: npolcon,npolcon2(0:natom),natpri(natom)
real*8, intent(in) :: spinpol(nspv-1,0:natom),polconb(3,natom)
real*8, intent(out):: enebc
integer:: na,ns 

  enebc= 0.0d0 
  if ((npolcon==key_polcon_atoms).or.(npolcon==key_polcon_asa)) then
    do na= 1,natom
      if (natpri(na)==key_natpri_in) then
      if (npolcon2(na)==key_polcon2_size) then
        do ns= 1,nspv-1
          enebc= enebc -spinpol(ns,na)*polconb(ns,na)
        enddo
      endif
      endif
    enddo
  endif 

end subroutine totalenergy_space_06


end module
