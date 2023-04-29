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
! **********  dim_global8f.F90  06/20/2018-01  **********

module dim_arrays
contains

subroutine allocate_arrays
! unused arrays are allocated with length 1, as they are are passed to some general subroutines
use mod_mpi
use var_arrays
use var_read_input_kukan
use var_read_ppfile,      only: nradmx,nprjmx,lrhomx
implicit none

  if (myrank_glbl<nprocw) then

  allocate(dij(nprjmx,nprjmx,nums*ncol,natom)) ! ?? it would be better to define dij on an array that's localized in space
  if (nso==0) then
    allocate(dijsoc(1,1,1,1))
  else
    allocate(dijsoc(nprjmx,nprjmx,3,natom))
  endif
  allocate(fatx(natom),faty(natom),fatz(natom))

  allocate(rhosmt_i(ncpx,ncpy,ncpz,nspv))
  allocate(rhosmt_o(ncpx,ncpy,ncpz,nspv))
  allocate(atocc((nprjmx**2+nprjmx)/2,nspv,num_atcell))

  if (nrc==0) then
    allocate(rspsep(nprjmx,num_ppcell,neigmx,nums,numk))
    allocate(svecre(ncpx,ncpy,ncpz,neigmx,nums,numk))
    allocate(ssvre(ncpx,ncpy,ncpz,neigmx,nums,numk))
    allocate(hsvre(ncpx,ncpy,ncpz,neigmx))
    allocate(cspsep(1,1,1,1,1),sveccm(1,1,1,1,1,1),ssvcm(1,1,1,1,1,1),hsvcm(1,1,1,1,1))
  else
    allocate(cspsep(nprjmx,num_ppcell,neigmx,nums,numk))
    allocate(sveccm(ncpx,ncpy,ncpz,neigmx,nums,numk))
    allocate(ssvcm(ncpx,ncpy,ncpz,neigmx,nums,numk))
    allocate(hsvcm(ncpx,ncpy,ncpz,neigmx,ncol))
    allocate(rspsep(1,1,1,1,1),svecre(1,1,1,1,1,1),ssvre(1,1,1,1,1,1),hsvre(1,1,1,1))
  endif

  allocate(sval(neigmx,nums+1-ncol,numk),fnele(neigmx,nums+1-ncol,numk))
  allocate(residual_states(neigmx,nums+1-ncol,numk))

  allocate(veff(ncpx,ncpy,ncpz,nspv))

  if ((nso==1).and.(ncol==1)) then
    allocate(svecsoc(ncpx,ncpy,ncpz,2*neigmx,2,numk))
    allocate(cspsepsoc(nprjmx,num_ppcell,2*neigmx,2,numk))
    allocate(svalsoc(2*neigmx,1,numk),fnelesoc(2*neigmx,1,numk))
  else
    allocate(svecsoc(1,1,1,1,1,1))
    allocate(cspsepsoc(1,1,1,1,1))
    allocate(svalsoc(1,1,1),fnelesoc(1,1,1))
  endif

  if (myr_kpt==0) then

    allocate(rhotrur_i(nradmx,npoint,nspv,num_atcell),rhosmtr_i(nradmx,npoint,nspv,num_atcell))
    allocate(rhotrur_o(nradmx,npoint,nspv,num_atcell),rhosmtr_o(nradmx,npoint,nspv,num_atcell))

    if ((nso==0).or.(ncol==2)) then
      allocate(sval_wfc(neigmx,nums+1-ncol,numkmx),fnele_wfc(neigmx,nums+1-ncol,numkmx))
    else
      allocate(sval_wfc(2*neigmx,1,numkmx),fnele_wfc(2*neigmx,1,numkmx))
    endif
    allocate(residual_states_wfc(neigmx,nums+1-ncol,numkmx))

    allocate(atmpole(10,natom))
    allocate( &
     vh_coarse(ncpx,ncpy,ncpz), &
     vx(ncpx,ncpy,ncpz,nspv), &
     rho_coarse(ncpx,ncpy,ncpz), &
     pwei(ncpx,ncpy,ncpz,natom))

    allocate( &
     vh_dense(ncpx_d,ncpy_d,ncpz_d), &
     vx_dense(ncpx_d,ncpy_d,ncpz_d,nspv), &
     ex_dense(ncpx_d,ncpy_d,ncpz_d), &
     rho_aug_dense(ncpx_d,ncpy_d,ncpz_d), &
     rhosmt_pcc_dense(ncpx_d,ncpy_d,ncpz_d,nspv) )

    allocate( rhoaugr(     nradmx,npoint,     num_atcell) )
    allocate( rhotrucorer( nradmx,npoint,nspv,num_atcell) )
    allocate( rhosmt_pccr( nradmx,npoint,nspv,num_atcell) )
    allocate( vxctru(nradmx,npoint,nspv,num_atcell) )
    allocate( exctru(nradmx,npoint,     num_atcell) )
    allocate( vxcsmt(nradmx,npoint,nspv,num_atcell) )
    allocate( excsmt(nradmx,npoint,     num_atcell) )
    allocate( vhtrur(nradmx,npoint,     num_atcell) )
    allocate( vhsmtr(nradmx,npoint,     num_atcell) )
    allocate( vhaugr(nradmx,npoint,     num_atcell) )
    allocate(vboundx(-(nfh-1):nfh,ncpy_d,ncpz_d,9*(3-nperi)/3+1,(natom-1)*(3-nperi)/3+1))
    allocate(vboundy(ncpx_d,-(nfh-1):nfh,ncpz_d,9*(3-nperi)/2+1,(natom-1)*(3-nperi)/2+1))
    allocate(vboundz(ncpx_d,ncpy_d,-(nfh-1):nfh,9*(1-nperi/3)+1,(natom-1)*(1-nperi/3)+1))

    allocate(point(npoint,3),wt(npoint) &
            ,         yylm(npoint,lrhomx),      dylm_dtheta(npoint,lrhomx) &
            ,d2ylm_dtheta2(npoint,lrhomx),        dylm_dphi(npoint,lrhomx) &
            ,  d2ylm_dphi2(npoint,lrhomx),d2ylm_dtheta_dphi(npoint,lrhomx))

    if (nspv>1) then
      allocate(spinpol(nspv-1,0:natom))
    else
      allocate(spinpol(1,0:0))
    endif

    allocate(natcell_inf(num_atcell))

  else

    allocate(rhotrur_i(1,1,1,1),rhosmtr_i(1,1,1,1))
    allocate(rhotrur_o(1,1,1,1),rhosmtr_o(1,1,1,1))

    allocate(sval_wfc(1,1,1),fnele_wfc(1,1,1))
    allocate(residual_states_wfc(1,1,1))

    allocate(spinpol(1,0:0))

  endif ! (myr_kpt==0) else

  endif ! (myrank_glbl<nprocw)

end subroutine allocate_arrays


subroutine allocate_atom(natom,npolcon,key_polcon_atoms,key_polcon_asa)
use var_read_ppfile,only:lmx  
use var_read_input_kukan, only:natsoc,numz,nmdx,nmdy,nmdz,atx,aty,atz,watom,npolcon2,polconb,polconmag,polconeta
use mod_mpi
implicit none
integer, intent(in)::natom,npolcon,key_polcon_atoms,key_polcon_asa
integer::npolcondim

  if (myrank_glbl<nprocw) then

  allocate(natsoc(natom,0:lmx-1))
  allocate(numz(natom))
  allocate(nmdx(natom),nmdy(natom),nmdz(natom))
  allocate(atx(natom),aty(natom),atz(natom),watom(natom))

  if (myr_kpt==0) then
    if ((npolcon==key_polcon_atoms).or.(npolcon==key_polcon_asa)) then
      npolcondim= natom
    else
      npolcondim= 1
    endif
    allocate(npolcon2(0:npolcondim))
    allocate(polconb(3,npolcondim),polconmag(3,npolcondim),polconeta(2,0:npolcondim))
  endif

  endif ! (myrank_glbl<nprocw)

end subroutine allocate_atom


subroutine allocate_skpxyz(numkmx)
use var_read_input_kukan, only:skpxyz,nwskp
implicit none
integer, intent(in)::numkmx
  allocate(nwskp(numkmx))
  allocate(skpxyz(numkmx,3))
  return
end subroutine allocate_skpxyz


subroutine allocate_skp(numk)
use var_read_input_kukan, only:skpxx,skpyy,skpzz,nwskk
implicit none
integer, intent(in)::numk
  allocate(nwskk(numk))
  allocate(skpxx(numk),skpyy(numk),skpzz(numk))
  return
end subroutine allocate_skp


end module dim_arrays

! ========================================================================================================

module dim_read_ppfile
contains

subroutine allocate_read_ppfile
use mod_mpi
use var_read_input_kukan, only: num_spe,nso 
use var_read_ppfile
implicit none

  if (myrank_glbl<nprocw) then

  ! fixed values
  allocate(ntyppp(num_spe))
  allocate(nprj(num_spe))

  if (myr_kpt==0) then
    ! fixed values
    allocate(grwein(15))
    allocate(nrprj(num_spe),lpmx(num_spe))
    allocate(npr(0:lmx-1,num_spe))
    allocate(sss0(nprmx*lmx,nprmx*lmx,num_spe),akv0(nprmx*lmx,nprmx*lmx,num_spe),prj(nradmx,nprmx*lmx,num_spe))
    allocate(gmax(num_spe))
    allocate(nradps(num_spe))
    allocate(nradct(num_spe))
    allocate(cp(8,num_spe))
    allocate(radial(nradmx,num_spe),dradial(nradmx,num_spe))
    allocate(potc(nradmx,num_spe),awf(nradmx,nprmx*lmx,num_spe),pwf(nradmx,nprmx*lmx,num_spe))
    allocate(rhocore(nradmx,num_spe))
    allocate(rhopcc(nradmx,num_spe))
    allocate(aeeig(nprmx*lmx,num_spe),aepot(nradmx,num_spe))
  endif ! (myr_kpt==0)

  endif ! (myrank_glbl<nprocw)

end subroutine allocate_read_ppfile

end module dim_read_ppfile

! ========================================================================================================

module dim_pseudopotentials
contains

subroutine allocate_pseudopotentials
use mod_mpi
use var_pseudopotentials
use var_read_ppfile, only:nradmx,nprjmx,nprmx,lmx,lrhomx
use var_read_input_kukan, only: &
 natom,num_spe,num_ppcell,num_ppcell_d,ncpx,ncpy,ncpz, &
 ncpx_d,ncpy_d,ncpz_d,npoint,num_atcell,nqmx
implicit none

  if (myrank_glbl<nprocw) then

  ! fixed values
  allocate(indspe(natom))
  allocate(nlind(nprjmx,num_spe),noind(nprjmx,num_spe))
  allocate(wail((nprmx**2+nprmx)/2,lmx,num_spe),wpil((nprmx**2+nprmx)/2,lmx,num_spe))
  allocate(mspe(num_spe))

  ! geometry-dependend values
  allocate(natpri(natom))
  allocate(natpri_inf(natom))
  allocate(naps(natom))
  allocate(natinf(natom))
  allocate(natx(natom),naty(natom),natz(natom))
  allocate(latom(natom))
  allocate(sss(nprjmx,nprjmx,num_spe))

  if (myr_kpt==0) then

    ! fixed values
    allocate(nwexp(num_spe))
    allocate(nqct(num_spe),nqctpcc(num_spe))
    allocate(coef(0:nqmx,0:7,num_spe))
    allocate(rfac(0:2*(lmx-1),num_spe))
    allocate(akv(nprjmx,nprjmx,num_spe))
    allocate(vjell(ncpz))

    ! geometry-dependend values
    allocate(ndatx(natom),ndaty(natom),ndatz(natom))
    allocate(natprid(natom))
    allocate(napsd(natom))
    allocate(natinfd(natom),natinfd_vloc(natom))
    allocate(vloc_coarse(ncpx,ncpy,ncpz))
    allocate(vloc_dense(ncpx_d,ncpy_d,ncpz_d))
    allocate(   vloc_scw(ncpx,ncpy,ncpz,natom))
    allocate(dvlocdx_scw(ncpx,ncpy,ncpz,natom))
    allocate(dvlocdy_scw(ncpx,ncpy,ncpz,natom))
    allocate(dvlocdz_scw(ncpx,ncpy,ncpz,natom))
    allocate(    vcorer(nradmx,npoint,num_atcell))
    allocate(vcorer_all(nradmx,npoint,num_atcell))
    allocate(qijl(nprjmx,nprjmx,lrhomx,natom))
    allocate(rhopcc_dense(ncpx_d,ncpy_d,ncpz_d))
    allocate(rhopccr(nradmx,npoint,num_atcell))

  endif ! (myr_kpt==0)

  endif ! (myrank_glbl<nprocw)

end subroutine allocate_pseudopotentials


subroutine deallocate_pseudopotentials
use var_pseudopotentials
use mod_mpi
implicit none

  if (myrank_glbl<nprocw) then

  ! fixed values
  deallocate(indspe)
  deallocate(wail,wpil)
  deallocate(mspe)

  ! geometry-dependend values
  deallocate(nlind,noind)
  deallocate(natpri)
  deallocate(natpri_inf)
  deallocate(naps)
  deallocate(natinf)
  deallocate(natx,naty,natz)
  deallocate(sss)

  if (myr_kpt==0) then

    ! fixed values
    deallocate(nwexp)
    deallocate(nqct,nqctpcc)
    deallocate(coef)
    deallocate(rfac)
    deallocate(akv)
    deallocate(vjell)

    ! geometry-dependend values
    deallocate(ndatx,ndaty,ndatz)
    deallocate(natprid)
    deallocate(napsd)
    deallocate(natinfd,natinfd_vloc)
    deallocate(vloc_coarse)
    deallocate(vloc_dense)
    deallocate(   vloc_scw)
    deallocate(dvlocdx_scw)
    deallocate(dvlocdy_scw)
    deallocate(dvlocdz_scw)
    deallocate(    vcorer)
    deallocate(vcorer_all)
    deallocate(qijl)
    deallocate(rhopcc_dense)
    deallocate(rhopccr)

  endif ! (myr_kpt==0)

  endif ! (myrank_glbl<nprocw)

end subroutine deallocate_pseudopotentials

end module dim_pseudopotentials

! ========================================================================================================

module dim_listvec
contains

subroutine allocate_listvec
use mod_mpi
use var_listvec
use var_read_input_kukan ,only: num_ppcell,num_ppcell_d
use var_read_ppfile,only: nprjmx

  if (myrank_glbl<nprocw) then

  allocate(lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell))
  allocate(lstvec2(num_list,num_ppcell))
  allocate(vnlocp(num_list,nprjmx,num_ppcell))
  if (myr_kpt==0) then
    allocate(lstdx(num_list_d,num_ppcell_d),lstdy(num_list_d,num_ppcell_d),lstdz(num_list_d,num_ppcell_d))
    allocate(lstvecd2(num_list_d,num_ppcell_d))
    allocate(   vloc_hdp(num_list_d,num_ppcell_d))
    allocate(dvlocdx_hdp(num_list_d,num_ppcell_d))
    allocate(dvlocdy_hdp(num_list_d,num_ppcell_d))
    allocate(dvlocdz_hdp(num_list_d,num_ppcell_d))
    allocate(   rhoaug3d(num_list_d,num_ppcell_d))
    allocate(drhoaug3ddx(num_list_d,num_ppcell_d))
    allocate(drhoaug3ddy(num_list_d,num_ppcell_d))
    allocate(drhoaug3ddz(num_list_d,num_ppcell_d))
  endif

  endif ! (myrank_glbl<nprocw)

end subroutine allocate_listvec

end module dim_listvec

! ========================================================================================================

module dim_kslaplacian
contains

subroutine allocate_kslaplacian
use mod_mpi
use var_kslaplacian, only: acfd
implicit none

  if (myrank_glbl<nprocw) then

    allocate(acfd(0:8))

  endif ! (myrank_glbl<nprocw)

end subroutine allocate_kslaplacian


subroutine deallocate_kslaplacian
use mod_mpi
use var_kslaplacian, only: acfd
implicit none

  if (myrank_glbl<nprocw) then

    deallocate(acfd)

  endif ! (myrank_glbl<nprocw)

end subroutine deallocate_kslaplacian

end module dim_kslaplacian

! ========================================================================================================

