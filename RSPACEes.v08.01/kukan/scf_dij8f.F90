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
! **********  scf_dij8f.F90 11/19/2013-01  **********

module mod_scf_dij
implicit none
real*8,allocatable::vtmpm(:,:,:)
contains

subroutine scf_dij( &
 natom,npolcondim,num_spe,num_atcell,nradmx,nprjmx,nums,ncol,nspv,npoint, & ! <
 lrhomx,num_list_d,num_ppcell_d,lmx,nprmx,                                & ! <
 ncpx_d,ncpy_d,ncpz_d,                                                    & ! <
 npolcon,                                                                 & ! <
 key_natpri_in,key_natpri_inps,key_pp_paw,                                & ! <
 key_polcon_atoms,key_polcon_asa,key_polcon2_none,                        & ! <
 ddx,ddy,ddz,                                                             & ! <
 npolcon2,ntyppp,nradct,nprj,npr,indspe,natpri,natpri_inf,                & ! <
 nwexp,natinfd,napsd,natprid,lpmx,                                        & ! <
 nlind,noind,lstvecd2,                                                    & ! <
 lstdx,lstdy,lstdz,ndatx,ndaty,ndatz,                                     & ! <
 awf,pwf,wail,radial,dradial,rfac,                                        & ! <
 atx,aty,atz,                                                             & ! <
 vh_dense,vloc_dense,                                                     & ! <
 vhtrur,vhsmtr,vhaugr,                                                    & ! <
 vxctru,vxcsmt,vcorer_all,vcorer,qijl,akv,polconb,yylm,wt,                & ! <
 dij)                                                                       ! >
use mod_mpi
use mod_stopp
implicit none
integer,intent(in) ::natom,npolcondim,num_spe,num_atcell,nradmx,nprjmx,nums,ncol,nspv,npoint
integer,intent(in) ::lrhomx,num_list_d,num_ppcell_d,lmx,nprmx
integer,intent(in) ::ncpx_d,ncpy_d,ncpz_d
integer,intent(in) ::npolcon
integer,intent(in) ::key_natpri_in,key_natpri_inps,key_pp_paw,key_polcon_atoms,key_polcon_asa,key_polcon2_none
integer,intent(in) ::npolcon2(0:npolcondim),ntyppp(num_spe),nradct(num_spe),nprj(num_spe),npr(0:lmx-1,num_spe),indspe(natom)
integer,intent(in) ::natpri(natom),natpri_inf(natom),nwexp(num_spe),natinfd(natom),napsd(natom),natprid(natom),lpmx(num_spe)
integer,intent(in) ::nlind(nprjmx,num_spe),noind(nprjmx,num_spe),lstvecd2(num_list_d,num_ppcell_d)
integer,intent(in) ::lstdx(num_list_d,num_ppcell_d),lstdy(num_list_d,num_ppcell_d),lstdz(num_list_d,num_ppcell_d)
integer,intent(in) ::ndatx(natom),ndaty(natom),ndatz(natom)
real*8, intent(in) ::ddx,ddy,ddz
real*8, intent(in) ::awf(nradmx,nprmx*lmx,num_spe),pwf(nradmx,nprmx*lmx,num_spe),wail((nprmx**2+nprmx)/2,lmx,num_spe)
real*8, intent(in) ::radial(nradmx,num_spe),dradial(nradmx,num_spe),rfac(0:2*(lmx-1),num_spe)
real*8, intent(in) ::atx(natom),aty(natom),atz(natom)
real*8, intent(in) ::vh_dense(ncpx_d,ncpy_d,ncpz_d)
real*8, intent(in) ::vloc_dense(ncpx_d,ncpy_d,ncpz_d)
real*8, intent(in) ::vhtrur(nradmx,npoint,num_atcell)
real*8, intent(in) ::vhsmtr(nradmx,npoint,num_atcell)
real*8, intent(in) ::vhaugr(nradmx,npoint,num_atcell)
real*8, intent(in) ::vxctru(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::vxcsmt(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::vcorer_all(nradmx,npoint,num_atcell)
real*8, intent(in) ::vcorer(nradmx,npoint,num_atcell)
real*8, intent(in) ::qijl(nprjmx,nprjmx,lrhomx,natom)
real*8, intent(in) ::akv(nprjmx,nprjmx,num_spe)
real*8, intent(in) ::polconb(3,npolcondim)
real*8, intent(in) ::yylm(npoint,lrhomx),wt(npoint)
real*8, intent(out)::dij(nprjmx,nprjmx,nums*ncol,natom)
integer na,ipri,ispe,ns,i,j,ij,ii,jj,lla,ipr,imm
real*8  tmpr
real*8,allocatable ::dijtmp(:,:,:)
real*8,allocatable ::dtilij_int(:,:)
real*8,allocatable ::vcom(:,:),vhm(:,:),vxcm(:,:,:)

  allocate( dijtmp(nprjmx,nprjmx,nums*ncol) )
! matrix elements \tilda{Dij^1} and \hat{Dij} [see PRB59 1758 (1999)],
! moments of one center Coulomb potential of nucleus, Hartree, and exchange correlation potentials
  allocate(  &
   dtilij_int((nprjmx**2+nprjmx)/2,num_atcell),  &
   vcom((nprjmx**2+nprjmx)/2,num_atcell), vhm((nprjmx**2+nprjmx)/2,num_atcell),  &
   vxcm((nprjmx**2+nprjmx)/2,nspv,num_atcell) )

! **********  <\psi_ae|v_core+v_h+v_xc|\psi_ae>-<\psi_ps|v_core+v_h+v_xc|\psi_ps>  **********
  allocate(vtmpm(nradmx,3+2*nums*ncol,num_atcell))
!$omp parallel default(shared)
  call scf_dij_01( &
   key_natpri_in,key_pp_paw,natom,num_spe,num_atcell,lrhomx,lmx,nprmx,nprjmx,nradmx,npoint,nspv, & ! <
   indspe,natpri,natpri_inf,ntyppp,nprj,nradct,nlind,noind,                                      & ! <
   yylm,wt,dradial,vhtrur,vhsmtr,vhaugr,vxctru,vxcsmt,vcorer,awf,pwf,                            & ! <
   vhm,vcom,vxcm)                                                                                  ! >
!$omp end parallel
  deallocate(vtmpm)

  do na=1,natom
    if (natpri(na) .eq. key_natpri_in) then
      ipri=natpri_inf(na)
      ij= 0
      do j=1,nprj(indspe(na))
        do i=1,j
          ij= ij + 1
          vcom(ij,ipri)=vcom(ij,ipri)+akv(i,j,indspe(na))
        end do
      end do
    end if
  end do
! *******************************************************************************************

! **********  \int (v_core+v_h) Q_{ij}^L dr [dtilij_int], \int (v_core+v_h) Q_{ij}^L dxdydz [dhatij]  **********
  if (lrhomx > 25) call stopp ('scf_dij: lrhomx must be <= 25!')
!$omp parallel default(shared)
  call scf_dij_02( &
   key_natpri_in,key_natpri_inps,key_pp_paw,                                                                        & ! <
   natom,num_spe,num_atcell,nums*ncol,nradmx,npoint,lmx,lrhomx,nprjmx,ncpx_d,ncpy_d,ncpz_d,num_list_d,num_ppcell_d, & ! <
   indspe,natpri,natpri_inf,natinfd,ndatx,ndaty,ndatz,napsd,natprid,ntyppp,nwexp,lpmx,nprj,nradct,                  & ! <
   lstdx,lstdy,lstdz,lstvecd2,                                                                                      & ! <
   ddx,ddy,ddz,radial,dradial,rfac,vhsmtr,vhaugr,vh_dense,vcorer_all,vloc_dense,qijl,atx,aty,atz,yylm,wt,           & ! <
   dtilij_int,dij)                                                                                                    ! >
!$omp end parallel
! now [ dij(:,:,1,:)=dhatij , dij(:,:,ns,:)=0.0d0 if ns>1 ]

  do na= 1,natom
    if (natprid(na) .eq. key_natpri_inps) then
      ispe= indspe(na)
      ij=0
      do j=1,nprj(ispe)
        do i=1,j
          ij=ij+1
          if (natpri(na) .eq. key_natpri_in) then
            ipri=natpri_inf(na)
            dij(i,j,1,na)= dij(i,j,1,na) -dtilij_int(ij,ipri)+vcom(ij,ipri)+vhm(ij,ipri)
          end if
          if (nums>1) dij(i,j,nums,na)= dij(i,j,1,na)
          if (natpri(na) .eq. key_natpri_in) then
            do ns= 1,nspv
              dij(i,j,ns,na)= dij(i,j,ns,na) + vxcm(ij,ns,ipri)
            end do
          end if
          if (nspv<ncol) dij(i,j,nums,na)= dij(i,j,1,na)
        end do
      end do
!     ********* <\psi_ae|B_con|\psi_ae> *********
      if ( ((npolcon==key_polcon_atoms).or.(npolcon==key_polcon_asa)) .and. (natpri(na)==key_natpri_in) ) then
        if (npolcon2(na)/=key_polcon2_none) then
          do ns= 1,nspv
            tmpr= polconb(max(1,ns-1),na)
            if (ns==2) tmpr= -tmpr
            ipr= 0
            do lla= 0,lpmx(ispe)
              ij= 0
              do j= 1,npr(lla,ispe)
                do i= 1,j
                  ij= ij +1
                  do imm= 1,2*lla+1
                    ii= ipr+(i-1)*(2*lla+1)+imm
                    jj= ipr+(j-1)*(2*lla+1)+imm
                    dij(ii,jj,ns,na)= dij(ii,jj,ns,na) +tmpr*wail(ij,lla+1,ispe)
                  enddo
                enddo
              enddo
              ipr= ipr+npr(lla,ispe)*(2*lla+1)
            enddo ! lla
          enddo ! ns
        endif
      endif ! ( ((npolcon==key_polcon_atoms).or.(npolcon==key_polcon_asa)) .and. (natpri(na)==key_natpri_in) )
!     *******************************************
!     fill upper triangle of D_ij matrix
      do ns=1,nums*ncol
        do j=1,nprj(indspe(na))
          do i=1,j-1
            dij(j,i,ns,na)= dij(i,j,ns,na)
          end do
        end do
      end do
    end if ! (natprid(na) .eq. key_natpri_inps)
  end do ! na
! ****************************************************************************************************
  do na= 1,natom
    call mpi_allreduce(dij(1,1,1,na),dijtmp(1,1,1),nprjmx*nprjmx*nums*ncol  &
                      ,mpi_double_precision,mpi_sum,mpicom_space,mpij)
    dij(:,:,:,na)= dijtmp(:,:,:)
  end do

  deallocate( dijtmp )
  deallocate( dtilij_int, vcom, vhm, vxcm )

  return

end subroutine scf_dij


subroutine scf_dij_01( &
 key_natpri_in,key_pp_paw,natom,num_spe,num_atcell,lrhomx,lmx,nprmx,nprjmx,nradmx,npoint,nspv, & ! <
 indspe,natpri,natpri_inf,ntyppp,nprj,nradct,nlind,noind,                                      & ! <
 yylm,wt,dradial,vhtrur,vhsmtr,vhaugr,vxctru,vxcsmt,vcorer,awf,pwf,                            & ! <
 vhm,vcom,vxcm)                                                                                  ! >
implicit none
integer,intent(in) ::key_natpri_in,key_pp_paw
integer,intent(in) ::natom,num_spe,num_atcell,lrhomx,lmx,nprmx,nprjmx,nradmx,npoint,nspv
integer,intent(in) ::indspe(natom),natpri(natom),natpri_inf(natom)
integer,intent(in) ::ntyppp(num_spe),nprj(num_spe),nradct(num_spe),nlind(nprjmx,num_spe),noind(nprjmx,num_spe)
real*8, intent(in) ::yylm(npoint,lrhomx),wt(npoint),dradial(nradmx,num_spe)
real*8, intent(in) ::vhtrur(nradmx,npoint,num_atcell),vhsmtr(nradmx,npoint,num_atcell),vhaugr(nradmx,npoint,num_atcell)
real*8, intent(in) ::vxctru(nradmx,npoint,nspv,num_atcell),vxcsmt(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::vcorer(nradmx,npoint,num_atcell)
real*8, intent(in) ::awf(nradmx,nprmx*lmx,num_spe),pwf(nradmx,nprmx*lmx,num_spe)
real*8, intent(out)::vhm((nprjmx**2+nprjmx)/2,num_atcell)
real*8, intent(out)::vcom((nprjmx**2+nprjmx)/2,num_atcell)
real*8, intent(out)::vxcm((nprjmx**2+nprjmx)/2,nspv,num_atcell)
integer na,ipri,ispe,j,j1,j2,i,i1,i2,ij,il,ir,ns
real*8, allocatable::sum(:)
real*8  tmp,dr,tmp0,tmp1

  allocate(sum(3+2*nspv))

!$omp do
  do ir=1,(nprjmx**2+nprjmx)/2*num_atcell
    vhm(ir,1)=0.0d0
    vcom(ir,1)=0.0d0
  end do
!$omp do
  do ir=1,(nprjmx**2+nprjmx)/2*nspv*num_atcell
    vxcm(ir,1,1)=0.0d0
  end do

  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
      ipri=natpri_inf(na)
      ispe=indspe(na)
      ij= 0
      do j=1,nprj(ispe)
      j1=nlind(j,ispe)
      j2=noind(j,ispe)
      do i=1,j
        ij= ij + 1
        i1=nlind(i,ispe)
        i2=noind(i,ispe)
!$omp do
        do ir=1,nradmx*(3+2*nspv)
          vtmpm(ir,1,ipri)=0.0d0
        end do
        do il=1,npoint
          tmp=yylm(il,j1)*yylm(il,i1)*wt(il)
!$omp do
          do ir=2,nradct(ispe)-1
            vtmpm(ir,1,ipri)=vtmpm(ir,1,ipri)+tmp*vhtrur(ir,il,ipri)
            vtmpm(ir,2,ipri)=vtmpm(ir,2,ipri)+tmp*(vhsmtr(ir,il,ipri)+vhaugr(ir,il,ipri))
            vtmpm(ir,3,ipri)=vtmpm(ir,3,ipri)+tmp*vcorer(ir,il,ipri)
          end do
          do ns= 1,nspv
!$omp do
            do ir=2,nradct(ispe)-1
              vtmpm(ir,2+2*ns,ipri)=vtmpm(ir,2+2*ns,ipri)+tmp*vxctru(ir,il,ns,ipri)
              vtmpm(ir,3+2*ns,ipri)=vtmpm(ir,3+2*ns,ipri)+tmp*vxcsmt(ir,il,ns,ipri)
            end do
          end do
        end do
        sum(:)= 0.0d0
!$omp do
        do ir=2,nradct(ispe)-1
          dr=dradial(ir,ispe)
          tmp0=awf(ir,j2,ispe)*awf(ir,i2,ispe)*dr
          tmp1=pwf(ir,j2,ispe)*pwf(ir,i2,ispe)*dr
          sum(1)=sum(1)+vtmpm(ir,1,ipri)*tmp0
          sum(2)=sum(2)+vtmpm(ir,2,ipri)*tmp1
          sum(3)=sum(3)+vtmpm(ir,3,ipri)*(tmp0-tmp1)
          do ns= 1,nspv
            sum(2+2*ns)=sum(2+2*ns)+vtmpm(ir,2+2*ns,ipri)*tmp0
            sum(3+2*ns)=sum(3+2*ns)+vtmpm(ir,3+2*ns,ipri)*tmp1
          end do
        end do
!$omp critical
        vhm(ij,ipri)=vhm(ij,ipri)+sum(1)-sum(2)
        vcom(ij,ipri)=vcom(ij,ipri)+sum(3)
        do ns=1,nspv
          vxcm(ij,ns,ipri)=vxcm(ij,ns,ipri)+sum(2+2*ns)-sum(3+2*ns)
        end do
!$omp end critical
!$omp barrier
      end do
      end do
    end if
  end do

  deallocate(sum)

  return

end subroutine scf_dij_01


subroutine scf_dij_02( &
 key_natpri_in,key_natpri_inps,key_pp_paw,                                                                       & ! <
 natom,num_spe,num_atcell,numsncol,nradmx,npoint,lmx,lrhomx,nprjmx,ncpx_d,ncpy_d,ncpz_d,num_list_d,num_ppcell_d, & ! <
 indspe,natpri,natpri_inf,natinfd,ndatx,ndaty,ndatz,napsd,natprid,ntyppp,nwexp,lpmx,nprj,nradct,                 & ! <
 lstdx,lstdy,lstdz,lstvecd2,                                                                                     & ! <
 ddx,ddy,ddz,radial,dradial,rfac,vhsmtr,vhaugr,vh_dense,vcorer_all,vloc_dense,qijl,atx,aty,atz,yylm,wt,          & ! <
 dtilij_int,dhatij)                                                                                                ! >
implicit none
integer,intent(in) ::key_natpri_in,key_natpri_inps,key_pp_paw
integer,intent(in) ::natom,num_spe,num_atcell,numsncol,nradmx,npoint,lmx,lrhomx,nprjmx
integer,intent(in) ::ncpx_d,ncpy_d,ncpz_d,num_list_d,num_ppcell_d
integer,intent(in) ::indspe(natom),natpri(natom),natpri_inf(natom),natinfd(natom)
integer,intent(in) ::ndatx(natom),ndaty(natom),ndatz(natom),napsd(natom),natprid(natom)
integer,intent(in) ::ntyppp(num_spe),nwexp(num_spe),lpmx(num_spe),nprj(num_spe),nradct(num_spe)
integer,intent(in) ::lstdx(num_list_d,num_ppcell_d),lstdy(num_list_d,num_ppcell_d),lstdz(num_list_d,num_ppcell_d)
integer,intent(in) ::lstvecd2(num_list_d,num_ppcell_d)
real*8, intent(in) ::ddx,ddy,ddz
real*8, intent(in) ::radial(nradmx,num_spe),dradial(nradmx,num_spe),rfac(0:2*(lmx-1),num_spe)
real*8, intent(in) ::vhsmtr(nradmx,npoint,num_atcell),vhaugr(nradmx,npoint,num_atcell),vh_dense(ncpx_d,ncpy_d,ncpz_d)
real*8, intent(in) ::vcorer_all(nradmx,npoint,num_atcell),vloc_dense(ncpx_d,ncpy_d,ncpz_d)
real*8, intent(in) ::qijl(nprjmx,nprjmx,lrhomx,natom)
real*8, intent(in) ::atx(natom),aty(natom),atz(natom)
real*8, intent(in) ::yylm(npoint,lrhomx)
real*8, intent(in) ::wt(npoint)
real*8, intent(out)::dtilij_int((nprjmx**2+nprjmx)/2,num_atcell)
real*8, intent(out)::dhatij(nprjmx,nprjmx,numsncol,natom) ! dhatij doesn't need a spin index, but dij does
integer na,ispe,ipri,i,j,ij,il,ir,iapsd
real*8 pi,sqfourpi,sq3fourpi,sq1516pi,sq0516pi,sq1504pi,sq0716pi,sq2132pi,sq10516pi,sq3532pi,sq10504pi &
      ,sq09256pi,sq4532pi,sq4564pi,sq31532pi,sq315256pi,sq4516pi,sq31516pi
real*8 rcut,rcutin,rcut2in,rcut3in,rcut4in,r,r2,r3,r4,x,y,z,x2,y2,z2,sum0,tmp,dr,aug0,aug1,aug2,aug3,aug4 &
      ,augbaser0,augbaser1,augbaser2,augbaser3,augbaser4 &
      ,rin,r2in,r3in,r4in,xy,yz,zx,xyz,ddxyz,tmppot,tmpo0,tmpo1,tmpo2,tmpo3,tmpo4
real*8 tmp01,tmp02,tmp03,tmp04,tmp05,tmp06,tmp07,tmp08,tmp09,tmp10 &
      ,tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,tmp20 &
      ,tmp21,tmp22,tmp23,tmp24,tmp25
real*8 sum01,sum02,sum03,sum04,sum05,sum06,sum07,sum08,sum09,sum10 &
      ,sum11,sum12,sum13,sum14,sum15,sum16,sum17,sum18,sum19,sum20 &
      ,sum21,sum22,sum23,sum24,sum25
real*8 ylm01,ylm02,ylm03,ylm04,ylm05,ylm06,ylm07,ylm08,ylm09,ylm10 &
      ,ylm11,ylm12,ylm13,ylm14,ylm15,ylm16,ylm17,ylm18,ylm19,ylm20 &
      ,ylm21,ylm22,ylm23,ylm24,ylm25

  pi= 4.0d0*datan(1.0d0)
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
  do i=1,(nprjmx**2+nprjmx)/2*num_atcell
    dtilij_int(i,1)=0.0d0
  end do
!$omp do
  do i=1,nprjmx*nprjmx*numsncol*natom
    dhatij(i,1,1,1)=0.0d0
  end do

  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
      ipri=natpri_inf(na)
      ispe=indspe(na)
      rcut=radial(nradct(ispe),ispe)
      rcutin=1.0d0/rcut
      rcut2in=rcutin*rcutin
      if (lpmx(ispe)>1) then
        rcut3in=rcut2in*rcutin
        rcut4in=rcut3in*rcutin
      end if
      ij= 0
      do j=1,nprj(ispe)
      do i=1,j
        ij= ij + 1
        sum0=0.0d0
        do il=1,npoint
          tmp01=yylm(il, 1)*wt(il)*qijl(i,j, 1,na)
          if (lpmx(ispe)>0) then
            tmp02=yylm(il, 2)*wt(il)*qijl(i,j, 2,na)
            tmp03=yylm(il, 3)*wt(il)*qijl(i,j, 3,na)
            tmp04=yylm(il, 4)*wt(il)*qijl(i,j, 4,na)
            tmp05=yylm(il, 5)*wt(il)*qijl(i,j, 5,na)
            tmp06=yylm(il, 6)*wt(il)*qijl(i,j, 6,na)
            tmp07=yylm(il, 7)*wt(il)*qijl(i,j, 7,na)
            tmp08=yylm(il, 8)*wt(il)*qijl(i,j, 8,na)
            tmp09=yylm(il, 9)*wt(il)*qijl(i,j, 9,na)
          end if
          if (lpmx(ispe)>1) then
            tmp10=yylm(il,10)*wt(il)*qijl(i,j,10,na)
            tmp11=yylm(il,11)*wt(il)*qijl(i,j,11,na)
            tmp12=yylm(il,12)*wt(il)*qijl(i,j,12,na)
            tmp13=yylm(il,13)*wt(il)*qijl(i,j,13,na)
            tmp14=yylm(il,14)*wt(il)*qijl(i,j,14,na)
            tmp15=yylm(il,15)*wt(il)*qijl(i,j,15,na)
            tmp16=yylm(il,16)*wt(il)*qijl(i,j,16,na)
            tmp17=yylm(il,17)*wt(il)*qijl(i,j,17,na)
            tmp18=yylm(il,18)*wt(il)*qijl(i,j,18,na)
            tmp19=yylm(il,19)*wt(il)*qijl(i,j,19,na)
            tmp20=yylm(il,20)*wt(il)*qijl(i,j,20,na)
            tmp21=yylm(il,21)*wt(il)*qijl(i,j,21,na)
            tmp22=yylm(il,22)*wt(il)*qijl(i,j,22,na)
            tmp23=yylm(il,23)*wt(il)*qijl(i,j,23,na)
            tmp24=yylm(il,24)*wt(il)*qijl(i,j,24,na)
            tmp25=yylm(il,25)*wt(il)*qijl(i,j,25,na)
          end if
          select case (lpmx(ispe))
          case (0)
!$omp do
            do ir=2,nradct(ispe)-1
              r=radial(ir,ispe)
              r2=r*r
              dr=dradial(ir,ispe)
              aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
              aug3=aug4*(1.0d0-r2*rcut2in)
              aug2=aug3*(1.0d0-r2*rcut2in)
              aug1=aug2*(1.0d0-r2*rcut2in)
              aug0=aug1*(1.0d0-r2*rcut2in)
              augbaser0=rfac(0,ispe)*aug0
              tmp= augbaser0*tmp01
              sum0= sum0 + (vcorer_all(ir,il,ipri)+vhsmtr(ir,il,ipri)+vhaugr(ir,il,ipri)) *tmp*r*r*dr
            end do !ir
          case (1)
!$omp do
            do ir=2,nradct(ispe)-1
              r=radial(ir,ispe)
              r2=r*r
              dr=dradial(ir,ispe)
              aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
              aug3=aug4*(1.0d0-r2*rcut2in)
              aug2=aug3*(1.0d0-r2*rcut2in)
              aug1=aug2*(1.0d0-r2*rcut2in)
              aug0=aug1*(1.0d0-r2*rcut2in)
              augbaser0=rfac(0,ispe)*aug0
              augbaser1=rfac(1,ispe)*r *rcutin *aug1
              augbaser2=rfac(2,ispe)*r2*rcut2in*aug2
              tmp= augbaser0*tmp01  &
                  +augbaser1*(tmp02+tmp03+tmp04) &
                  +augbaser2*(tmp05+tmp06+tmp07+tmp08+tmp09)
              sum0= sum0 + (vcorer_all(ir,il,ipri)+vhsmtr(ir,il,ipri)+vhaugr(ir,il,ipri)) *tmp*r*r*dr
            end do !ir
          case (2)
!$omp do
            do ir=2,nradct(ispe)-1
              r=radial(ir,ispe)
              r2=r*r
              r3=r2*r
              r4=r3*r
              dr=dradial(ir,ispe)
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
              tmp= augbaser0*tmp01  &
                  +augbaser1*(tmp02+tmp03+tmp04) &
                  +augbaser2*(tmp05+tmp06+tmp07+tmp08+tmp09) &
                  +augbaser3*(tmp10+tmp11+tmp12+tmp13+tmp14+tmp15+tmp16) &
                  +augbaser4*(tmp17+tmp18+tmp19+tmp20+tmp21+tmp22+tmp23+tmp24+tmp25)
              sum0= sum0 + (vcorer_all(ir,il,ipri)+vhsmtr(ir,il,ipri)+vhaugr(ir,il,ipri)) *tmp*r*r*dr
            end do !ir
          end select
        end do !il
!$omp critical
        dtilij_int(ij,ipri)=dtilij_int(ij,ipri)+sum0
!$omp end critical
!$omp barrier
      end do !i
      end do !j
    end if
  end do

  ddxyz=ddx*ddy*ddz
  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natprid(na) .eq. key_natpri_inps)) then
    iapsd=napsd(na)
    ispe=indspe(na)
    rcut=radial(nradct(ispe),ispe)
    rcutin=1.0d0/rcut
    rcut2in=rcutin*rcutin

    sum01=0.0d0
    if (lpmx(ispe)>0) then
      sum02=0.0d0
      sum03=0.0d0
      sum04=0.0d0
      sum05=0.0d0
      sum06=0.0d0
      sum07=0.0d0
      sum08=0.0d0
      sum09=0.0d0
    end if
    if (lpmx(ispe)>1) then
      rcut3in=rcut2in*rcutin
      rcut4in=rcut3in*rcutin
      sum10=0.0d0
      sum11=0.0d0
      sum12=0.0d0
      sum13=0.0d0
      sum14=0.0d0
      sum15=0.0d0
      sum16=0.0d0
      sum17=0.0d0
      sum18=0.0d0
      sum19=0.0d0
      sum20=0.0d0
      sum21=0.0d0
      sum22=0.0d0
      sum23=0.0d0
      sum24=0.0d0
      sum25=0.0d0
    end if

    select case(lpmx(ispe))
    case(0)
!$omp do
      do i=1,natinfd(na)
        x=lstdx(i,iapsd)*ddx-atx(na)+(ndatx(na)*ddx-0.5d0*ddx)
        y=lstdy(i,iapsd)*ddy-aty(na)+(ndaty(na)*ddy-0.5d0*ddy)
        z=lstdz(i,iapsd)*ddz-atz(na)+(ndatz(na)*ddz-0.5d0*ddz)
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
          tmppot=vloc_dense(lstvecd2(i,iapsd),1,1)+vh_dense(lstvecd2(i,iapsd),1,1)
          ylm01= sqfourpi
          augbaser0=rfac(0,ispe)*aug0
          tmpo0=augbaser0*tmppot*ddxyz
          sum01=sum01+ylm01*tmpo0
!        end if
      end do
    case(1)
!$omp do
      do i=1,natinfd(na)
        x=lstdx(i,iapsd)*ddx-atx(na)+(ndatx(na)*ddx-0.5d0*ddx)
        y=lstdy(i,iapsd)*ddy-aty(na)+(ndaty(na)*ddy-0.5d0*ddy)
        z=lstdz(i,iapsd)*ddz-atz(na)+(ndatz(na)*ddz-0.5d0*ddz)
        x2=x*x
        y2=y*y
        z2=z*z
        xy=x*y
        yz=y*z
        zx=z*x
        r2=x2+y2+z2
        r=dsqrt(r2)
        rin=1.0d0/r
        r2in=rin*rin
!        if (r .le. radial(nradct(ispe),ispe)) then
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)
          tmppot=vloc_dense(lstvecd2(i,iapsd),1,1)+vh_dense(lstvecd2(i,iapsd),1,1)
          ylm01= sqfourpi
          ylm02= sq3fourpi*x*rin
          ylm03= sq3fourpi*y*rin
          ylm04= sq3fourpi*z*rin
          ylm05= sq1516pi*(x2-y2)*r2in
          ylm06=-sq1504pi*zx*r2in
          ylm07= sq0516pi*(3.0d0*z2-r2)*r2in
          ylm08=-sq1504pi*yz*r2in
          ylm09= sq1504pi*xy*r2in
          augbaser0=rfac(0,ispe)*aug0
          augbaser1=rfac(1,ispe)*r *rcutin *aug1
          augbaser2=rfac(2,ispe)*r2*rcut2in*aug2
          tmpo0=augbaser0*tmppot*ddxyz
          tmpo1=augbaser1*tmppot*ddxyz
          tmpo2=augbaser2*tmppot*ddxyz
          sum01=sum01+ylm01*tmpo0
          sum02=sum02+ylm02*tmpo1
          sum03=sum03+ylm03*tmpo1
          sum04=sum04+ylm04*tmpo1
          sum05=sum05+ylm05*tmpo2
          sum06=sum06+ylm06*tmpo2
          sum07=sum07+ylm07*tmpo2
          sum08=sum08+ylm08*tmpo2
          sum09=sum09+ylm09*tmpo2
!        end if
      end do
    case(2)
!$omp do
      do i=1,natinfd(na)
        x=lstdx(i,iapsd)*ddx-atx(na)+(ndatx(na)*ddx-0.5d0*ddx)
        y=lstdy(i,iapsd)*ddy-aty(na)+(ndaty(na)*ddy-0.5d0*ddy)
        z=lstdz(i,iapsd)*ddz-atz(na)+(ndatz(na)*ddz-0.5d0*ddz)
        x2=x*x
        y2=y*y
        z2=z*z
        xy=x*y
        yz=y*z
        zx=z*x
        xyz=xy*z
        r2=x2+y2+z2
        r=dsqrt(r2)
        r3=r2*r
        r4=r3*r
        rin=1.0d0/r
        r2in=rin*rin
        r3in=r2in*rin
        r4in=r3in*rin
!        if (r .le. radial(nradct(ispe),ispe)) then
          aug4=(1.0d0-r2*rcut2in)**(nwexp(ispe)-4)
          aug3=aug4*(1.0d0-r2*rcut2in)
          aug2=aug3*(1.0d0-r2*rcut2in)
          aug1=aug2*(1.0d0-r2*rcut2in)
          aug0=aug1*(1.0d0-r2*rcut2in)
          tmppot=vloc_dense(lstvecd2(i,iapsd),1,1)+vh_dense(lstvecd2(i,iapsd),1,1)
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
          augbaser0=rfac(0,ispe)*aug0
          augbaser1=rfac(1,ispe)*r *rcutin *aug1
          augbaser2=rfac(2,ispe)*r2*rcut2in*aug2
          augbaser3=rfac(3,ispe)*r3*rcut3in*aug3
          augbaser4=rfac(4,ispe)*r4*rcut4in*aug4
          tmpo0=augbaser0*tmppot*ddxyz
          tmpo1=augbaser1*tmppot*ddxyz
          tmpo2=augbaser2*tmppot*ddxyz
          tmpo3=augbaser3*tmppot*ddxyz
          tmpo4=augbaser4*tmppot*ddxyz
          sum01=sum01+ylm01*tmpo0
          sum02=sum02+ylm02*tmpo1
          sum03=sum03+ylm03*tmpo1
          sum04=sum04+ylm04*tmpo1
          sum05=sum05+ylm05*tmpo2
          sum06=sum06+ylm06*tmpo2
          sum07=sum07+ylm07*tmpo2
          sum08=sum08+ylm08*tmpo2
          sum09=sum09+ylm09*tmpo2
          sum10=sum10+ylm10*tmpo3
          sum11=sum11+ylm11*tmpo3
          sum12=sum12+ylm12*tmpo3
          sum13=sum13+ylm13*tmpo3
          sum14=sum14+ylm14*tmpo3
          sum15=sum15+ylm15*tmpo3
          sum16=sum16+ylm16*tmpo3
          sum17=sum17+ylm17*tmpo4
          sum18=sum18+ylm18*tmpo4
          sum19=sum19+ylm19*tmpo4
          sum20=sum20+ylm20*tmpo4
          sum21=sum21+ylm21*tmpo4
          sum22=sum22+ylm22*tmpo4
          sum23=sum23+ylm23*tmpo4
          sum24=sum24+ylm24*tmpo4
          sum25=sum25+ylm25*tmpo4
!        end if
      end do
    end select

!$omp critical
    do j=1,nprj(ispe)
      do i=1,j
        dhatij(i,j,1,na)=dhatij(i,j,1,na) &
                   +qijl(i,j, 1,na)*sum01
        if (lpmx(ispe)>0) dhatij(i,j,1,na)=dhatij(i,j,1,na) &
                   +qijl(i,j, 2,na)*sum02+qijl(i,j, 3,na)*sum03 &
                   +qijl(i,j, 4,na)*sum04+qijl(i,j, 5,na)*sum05+qijl(i,j, 6,na)*sum06 &
                   +qijl(i,j, 7,na)*sum07+qijl(i,j, 8,na)*sum08+qijl(i,j, 9,na)*sum09
        if (lpmx(ispe)>1) dhatij(i,j,1,na)=dhatij(i,j,1,na) &
                   +qijl(i,j,10,na)*sum10+qijl(i,j,11,na)*sum11+qijl(i,j,12,na)*sum12 &
                   +qijl(i,j,13,na)*sum13+qijl(i,j,14,na)*sum14+qijl(i,j,15,na)*sum15 &
                   +qijl(i,j,16,na)*sum16+qijl(i,j,17,na)*sum17+qijl(i,j,18,na)*sum18 &
                   +qijl(i,j,19,na)*sum19+qijl(i,j,20,na)*sum20+qijl(i,j,21,na)*sum21 &
                   +qijl(i,j,22,na)*sum22+qijl(i,j,23,na)*sum23+qijl(i,j,24,na)*sum24 &
                   +qijl(i,j,25,na)*sum25
      end do
    end do
!$omp end critical
!$omp barrier
    end if
  end do

  return
end subroutine scf_dij_02


subroutine scf_dij_soc( &
 natom,num_atcell,num_spe,nprmx,nprjmx,lmx,nradmx,npoint,nspv,            & ! < 
 key_natpri_in,key_pp_paw,key_soc_calc,                                   & ! <
 natpri,natpri_inf,indspe,natsoc,mspe,ntyppp,nradct,nprj,npr,nlind,noind, & ! <
 socang,radial,dradial,wt,awf,aeeig,aepot,                                & ! <
 rhotrucorer,vxctru,                                                      & ! <
 dijsoc)                                                                    ! > 
use mod_mpi
use mod_stopp
implicit none
integer,intent(in) ::natom,num_atcell,num_spe,nprmx,nprjmx,lmx,nradmx,npoint,nspv
integer,intent(in) ::key_natpri_in,key_pp_paw,key_soc_calc
integer,intent(in) ::natpri(natom),natpri_inf(natom),indspe(natom),natsoc(natom,0:lmx-1)
integer,intent(in) ::mspe(num_spe),ntyppp(num_spe),nradct(num_spe)
integer,intent(in) ::nprj(num_spe),npr(0:lmx-1,num_spe),nlind(nprjmx,num_spe),noind(nprjmx,num_spe)
real*8, intent(in) ::socang(3)
real*8, intent(in) ::radial(nradmx,num_spe),dradial(nradmx,num_spe),wt(npoint)
real*8, intent(in) ::awf(nradmx,nprmx*lmx,num_spe),aeeig(nprmx*lmx,num_spe),aepot(nradmx,num_spe)
real*8, intent(in) ::rhotrucorer(nradmx,npoint,nspv,num_atcell),vxctru(nradmx,npoint,nspv,num_atcell)
real*8, intent(out)::dijsoc(nprjmx,nprjmx,3,natom)
real*8 :: pi
integer:: j,na
integer:: ir0,ispe,ipri,nrad
integer,allocatable:: iatsoc(:)
real*8, allocatable:: crot(:),srot(:),vxcr(:),dvxcr(:,:),radrho(:)
real*8, allocatable:: zeta(:,:),dijtmp(:,:,:)

  allocate( iatsoc(lmx-1) )
  allocate( crot(3), srot(3) )
  allocate( vxcr(nradmx), dvxcr(2,nradmx+1), radrho(nradmx) )
  allocate( zeta((nprmx*(nprmx+1))/2,lmx-1) )

  pi= 4.0d0*datan(1.0d0) 
  do j= 1,3
    crot(j)= dcos(pi*socang(j))
    srot(j)= dsin(pi*socang(j))
  enddo 

  !$omp parallel default(shared) private(j)  
  !$omp do
  do j= 1,nprjmx*nprjmx*3*natom
    dijsoc(j,1,1,1)= 0.0d0
  enddo
  !$omp end parallel  

  do na=1,natom
  if ((natpri(na)==key_natpri_in).and.(natsoc(na,0)==key_soc_calc)) then
    ispe= indspe(na)
    ipri= natpri_inf(na)
    nrad= nradct(ispe)
    do j= 1,lmx-1
      iatsoc(j)= natsoc(na,j)
    enddo  

    if (ntyppp(ispe)/=key_pp_paw) call stopp('scf_dij_soc: Spin-orbit coupling can only be calculated for PAW atoms.') 

    !$omp parallel default(shared) 

    call scf_dij_soc_radial( &
     num_atcell,num_spe,nprmx,lmx,nradmx,npoint,nspv,ispe,ipri,nrad, & ! <
     npr,mspe,wt,radial,dradial,awf,aeeig,aepot,rhotrucorer,vxctru,  & ! <
     zeta,                                                           & ! > 
     ir0,vxcr,dvxcr,radrho)                                            ! W

    call scf_dij_soc_angular( &
     nprmx,nprjmx,lmx,num_spe,ispe,          & ! <
     key_soc_calc,                           & ! <
     iatsoc,nprj,nlind,noind,crot,srot,zeta, & ! <
     dijsoc(1,1,1,na))                         ! >

    !$omp end parallel

  endif
  enddo ! na

  deallocate( iatsoc,crot,srot,vxcr,dvxcr,radrho,zeta ) 

  allocate( dijtmp(nprjmx,nprjmx,3) )
  do na= 1,natom
    call mpi_allreduce(dijsoc(1,1,1,na),dijtmp(1,1,1),nprjmx*nprjmx*3,mpi_double_precision,mpi_sum,mpicom_space,mpij)
    dijsoc(:,:,:,na)= dijtmp(:,:,:)
  end do
  deallocate( dijtmp ) 

end subroutine scf_dij_soc


subroutine scf_dij_soc_radial( &
 num_atcell,num_spe,nprmx,lmx,nradmx,npoint,nspv,ispe,ipri,nrad, & ! <
 npr,mspe,wt,radial,dradial,awf,aeeig,aepot,rhotrucorer,vxctru,  & ! <
 zeta,                                                           & ! >
 ir0,vxcr,dvxcr,radrho)                                            ! W
use mod_stopp
implicit none
real*8,parameter:: epsrad = 1.0d-7         ! minimal radius that is not regarded as zero 
real*8,parameter:: epsdrad= 1.0d-20        ! minimal radius difference that is not regarded as zero 
real*8,parameter:: fineinv= 137.0359895d0  ! 1 / (fine structure constant)
integer,intent(in) ::num_atcell,num_spe,nprmx,lmx,nradmx,npoint,nspv,ispe,ipri,nrad
integer,intent(in) ::npr(0:lmx-1,num_spe),mspe(num_spe)
real*8, intent(in) ::wt(npoint),radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8, intent(in) ::awf(nradmx,nprmx*lmx,num_spe),aeeig(0:lmx-1,num_spe),aepot(nradmx,num_spe)
real*8, intent(in) ::rhotrucorer(nradmx,npoint,nspv,num_atcell),vxctru(nradmx,npoint,nspv,num_atcell)
real*8, intent(out)::zeta((nprmx*(nprmx+1))/2,lmx-1) 
integer,intent(out)::ir0 
! real*8, intent(out)::vxcr(nradmx),dvxcr(2,nradmx+1),radrho(nradmx)
real*8, intent(out)::vxcr(nradmx),dvxcr(2,nradmx),radrho(nradmx)
integer:: mspv,ns,ir,il,i2,j2,ij2,ill
real*8 :: pi,c2,dvxc,vso,tmp  
real*8, allocatable:: zeta_priv(:,:)

  allocate( zeta_priv((nprmx*(nprmx+1))/2,lmx-1) )

  pi= 4.0d0*datan(1.0d0) 
  c2= fineinv**2
  mspv= min(2,nspv) 

  !$omp single
  ir0= 1
  do ir= 1,2 
    if (( radial(ir,ispe)<epsrad ).or.( radial(ir+1,ispe)-radial(ir,ispe)<epsdrad) ) ir0= ir+1
  enddo 
  if (ir0>2) call stopp('scf_dij_soc_radial: first radial value(s) too small or too dense')
  !$omp end single 

  !$omp do
  do i2= 1,((nprmx*(nprmx+1))/2)*(lmx-1)
    zeta(i2,1)     = 0.0d0
  enddo 
  zeta_priv(:,:)= 0.0d0 

  !$omp do 
  do ir= ir0,nrad
    vxcr(ir)= 0.0d0
    radrho(ir)= 0.0d0 
    do ns= 1,mspv 
      do il= 1,npoint 
        vxcr(  ir)= vxcr(  ir) +vxctru(     ir,il,ns,ipri)*wt(il)
        radrho(ir)= radrho(ir) +rhotrucorer(ir,il,ns,ipri)*wt(il)
      enddo 
    enddo 
    vxcr(  ir)= vxcr(  ir)/(4.0d0*pi*dble(mspv))
    radrho(ir)= -radrho(ir)*radial(ir,ispe)**2*dradial(ir,ispe)
  enddo 

  !$omp do 
  do ir= ir0+1,nrad-2 ! ?? 
  ! do ir= ir0+1,nrad-1
    dvxcr(1,ir)= (vxcr(ir  )-vxcr(ir-1)) / (radial(ir  ,ispe)-radial(ir-1,ispe))
    dvxcr(2,ir)= (vxcr(ir+1)-vxcr(ir-1)) / (radial(ir+1,ispe)-radial(ir-1,ispe))
  enddo 
  !$omp single
  dvxcr(1,ir0   )= 0.0d0
  dvxcr(2,ir0   )= 0.0d0
  dvxcr(1,nrad-1)= (vxcr(nrad-1)-vxcr(nrad-2)) / (radial(nrad-1,ispe)-radial(nrad-2,ispe))
  dvxcr(2,nrad-1)= 0.0d0
  dvxcr(1,nrad  )= 0.0d0
  ! dvxcr(1,nrad  )= (vxcr(nrad)-vxcr(nrad-1)) / (radial(nrad,ispe)-radial(nrad-1,ispe))
  ! dvxcr(2,nrad  )= 0.0d0
  ! dvxcr(1,nrad+1)= 0.0d0

  do ir= ir0+1,nrad
    radrho(ir)= radrho(ir-1) +radrho(ir)
  enddo
  !$omp end single 

  !$omp do 
  do ir= ir0,nrad-1
  ! do ir= ir0,nrad
    dvxc= dvxcr(1,ir)+dvxcr(1,ir+1)-dvxcr(2,ir) 
    tmp= ( dvxc +(dble(mspe(ispe))+radrho(ir))/radial(ir,ispe)**2 ) * c2 * (dradial(ir,ispe)/radial(ir,ispe))
    do ill= 1,lmx-1 
      ij2= 0 
      do j2= ill*nprmx+1,ill*nprmx+npr(ill,ispe)
        do i2= ill*nprmx+1,j2
          ij2= ij2+1 
          vso= tmp / ( 2.0d0*c2 +(aeeig(i2,ispe)+aeeig(j2,ispe))/2.0d0 -aepot(ir,ispe)/radial(ir,ispe) )**2
          zeta_priv(ij2,ill)= zeta_priv(ij2,ill) + awf(ir,i2,ispe)*awf(ir,j2,ispe)*vso
        enddo
      enddo
    enddo   
  enddo ! ir
  !$omp critical 
  zeta(:,:)= zeta(:,:) +zeta_priv(:,:) 
  !$omp end critical 

  deallocate( zeta_priv )

end subroutine scf_dij_soc_radial 


subroutine scf_dij_soc_angular( &
 nprmx,nprjmx,lmx,num_spe,ispe,          & ! <
 key_soc_calc,                           & ! <
 iatsoc,nprj,nlind,noind,crot,srot,zeta, & ! < 
 dijsoc)                                   ! > 
use mod_stopp
implicit none
! data table for angular matrix elements
integer,dimension(3,13), parameter :: myly   = reshape( (/ &
 +1,0,0, 0,+1,0, 0,0,+1, 0,-1,0, 0,0,0 , 0,-3,0, 0,0,+1, +1,0,0, 0,0,+3, 16,0,0, 0,0,+1, 0,0,0 , 0,+1,0 /),shape(myly))
 ! 3,2     4,2     4,3     6,5     7,5     7,6     8,5     8,6     8,7     9,5     9,6     9,7     9,8
integer,intent(in) ::nprmx,nprjmx,lmx,num_spe,ispe 
integer,intent(in) ::key_soc_calc
integer,intent(in) ::iatsoc(1:lmx-1),nprj(num_spe),nlind(nprjmx,num_spe),noind(nprjmx,num_spe)
real*8, intent(in) ::crot(3),srot(3),zeta((nprmx*(nprmx+1))/2,lmx-1)
real*8, intent(out)::dijsoc(nprjmx,nprjmx,3)
integer:: i,j,i1,j1,i2,j2,ij2,ill,jll,ipos  
real*8,allocatable:: tmp(:) 

  allocate( tmp(3) )

  !$omp do 
  do j= 1,nprj(ispe)
    do i= 1,nprj(ispe)
      dijsoc(i,j,1)= 0.0d0
      dijsoc(i,j,2)= 0.0d0
      dijsoc(i,j,3)= 0.0d0
    enddo 
  enddo 
  !$omp do 
  do j= 1,nprj(ispe)
    j1= nlind(j,ispe)
    if ((j1>1).and.(j<nprj(ispe))) then  
      jll= 1
      if (j1>4) jll= 2  
      if (j1>9) call stopp ('scf_dij_soc_angular: f-states not implemented!')
      j2= noind(j,ispe)
      do i= j+1,nprj(ispe)
        i1= nlind(i,ispe)
        ill= 0
        if (i1>1) ill= 1  
        if (i1>4) ill= 2
        if (i1>9) ill= 3
        if ( (ill==jll) .and. (i1/=j1) .and. (iatsoc(ill)==key_soc_calc) ) then
          ipos= max(i1,j1)-ill**2
          ipos= (4*ill**3-3*ill**2-ill)/6 +(ipos**2-3*ipos+2)/2 +min(i1,j1)-ill**2
          i2= noind(i,ispe)
          ij2= max(i2,j2)-ill*nprmx-1
          ij2= (ij2*(ij2+1))/2 +min(i2,j2)-ill*nprmx
          tmp(1)= dsqrt(dble(abs(myly(1,ipos))))*dble(sign(1,(i1-j1)*myly(1,ipos))) * zeta(ij2,ill)
          tmp(2)= dsqrt(dble(abs(myly(2,ipos))))*dble(sign(1,(i1-j1)*myly(2,ipos))) * zeta(ij2,ill)
          tmp(3)= dsqrt(dble(abs(myly(3,ipos))))*dble(sign(1,(i1-j1)*myly(3,ipos))) * zeta(ij2,ill)
          dijsoc(i,j,1)= +crot(1)*tmp(1) +srot(1)*(crot(2)*tmp(3)-srot(2)*tmp(2)) 
          dijsoc(i,j,2)= +srot(2)*tmp(3) +crot(2)*tmp(2) 
          dijsoc(i,j,3)= -srot(1)*tmp(1) +crot(1)*(crot(2)*tmp(3)-srot(2)*tmp(2)) 
          dijsoc(j,i,1)= -dijsoc(i,j,1)
          dijsoc(j,i,2)= -dijsoc(i,j,2)
          dijsoc(j,i,3)= -dijsoc(i,j,3)
        endif 
      enddo ! i 
    endif ! (j1>1) 
  enddo ! j 

  deallocate( tmp )

end subroutine scf_dij_soc_angular


end module

