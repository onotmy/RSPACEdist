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
! **********  kukan8f.F90 04/27/2023-01  **********

use mod_mpi
use mod_stopp
use mod_disptime

use var_keys
use var_scalars 
use var_arrays
use var_read_input_kukan
use var_read_ppfile
use var_pseudopotentials
use var_listvec
use var_kslaplacian

use dim_arrays
use dim_read_ppfile
use dim_pseudopotentials
use dim_listvec

use mod_output,                 only: output_compcond, output_countdat, output_energy, output_bands, output_molecdyn, &
                                      output_force, output_polcon, output_xcrysden
use mod_read_input_kukan,       only: read_input_readinp, &
                                      read_input_atomxyz, read_input_atommag, read_input_mpisetup, read_input_kpt, &
                                      read_input_kptsnd, read_input_kptrcv
use mod_read_ppfile,            only: read_ppfile
use mod_rayleighritz,           only: rayleighritz_blacs_ini,rayleighritz_blacs_fin
use mod_pseudocalc,             only: pseudocalc, pseudocalc_broadcast
use mod_kslaplacian,            only: kslaplacian_initialize
use mod_scf_charge,             only: scf_charge
use mod_scf_chargemixing,       only: scf_chargemixing,scf_chargemixing_initialize,scf_chargemixing_finalize
use mod_force,                  only: force_eig_r,force_eig_c,force
use mod_scf_diag,               only: scf_diag
use mod_orthogonalization,      only: orthogonalization_r,orthogonalization_c
use mod_setup_pseudopot,        only: setup_pseudopot,setup_pseudopot_broadcast
use mod_filterpseudopotentials, only: filterpseudopotentials
use mod_fuzzycell,              only: fuzzycell
use mod_totalenergy,            only: totalenergy_space,totalenergy_eig
use mod_setup_initialwave,      only: setup_initialwave
use mod_rw_wfpot,               only: read_wf,read_rho,read_vht,write_wf,write_rhovht,get_wfpotdir
use mod_scf_potentials,         only: scf_potentials
use mod_scf_augcharge,          only: scf_augcharge
use mod_scf_dij,                only: scf_dij,scf_dij_soc
use mod_scf_occupation,         only: scf_occupation
use mod_scf_diffcharge,         only: scf_diffcharge
use mod_scf_polcon,             only: scf_polcon,scf_polcon_magproj,scf_polcon_initialcheck
use mod_trans,                  only: trans_d2c_smtcharge
use mod_tools,                  only: tools_sphglpts,tools_genyylm,tools_listvecdim,tools_indspecies,tools_jelliumpotential, &
                                      tools_shiftatoms,tools_calenespr,tools_moveatoms,tools_countelectron,tools_scffile, &
                                      tools_nlpcomsplit,tools_nlpcomsplit_nnn,tools_potbroadcast,tools_rodr_sum_nlp, &
                                      tools_maginitrot,tools_deri_yylm
use mod_proccount,              only: proccount

implicit none

integer:: i,n,iscf,imd,jfile
integer:: lbrydn
logical:: l_newdens 
logical:: scfctr(3),scfctr2,l_scf
logical:: lspec
logical:: ldirinout
integer:: ix,iy,iz
real*8 :: stime,stimesub,stime0
real*8 :: pi
character:: chdirinp*200,chdirout*200,fname*200

  call mpi_init(mpij)
  call mpi_comm_rank(mpi_comm_world,myrank_glbl,mpij)

  if (myrank_glbl == 0) then
    inquire(file='iodir.txt',exist=ldirinout)
    if (ldirinout) then
      read(5,*) chdirinp
      read(5,*) chdirout
    else
      chdirinp=''
      chdirout=''
    end if
  end if
  call mpi_bcast(chdirinp,len(chdirinp),mpi_character,0,mpi_comm_world,mpij)
  call mpi_bcast(chdirout,len(chdirout),mpi_character,0,mpi_comm_world,mpij)

  comtime=0.0d0
  comcount=0.0d0

  nndisp=66
  if (myrank_glbl==0) then
    fname='output_kukan.dat'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(nndisp,file=fname)
    fname='mdresult.dat'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(12,file=fname)
    write(6,*) 'computed by kukan8f 04/27/2023-01'
    write(nndisp,*) 'computed by kukan8f 04/27/2023-01'
    write(12,*)    'computed by kukan8f 04/27/2023-01'
  end if !myrank_glbl==0
  call stopp_initialize (nndisp)
  call proccount(6)
  call proccount(nndisp)
  call proccount(12)
  call read_input_readinp( &
   chdirinp,chdirout,ndisp,                                                     & ! <
   key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp,key_jel_nocalc,key_jel_calc, & ! <
   key_polcon_none,key_polcon_occ,key_ksym_inv,                                 & ! <
   nxmax,nymax,nzmax,npxmax,npymax,npzmax,                                      & ! >
   nso,nsym,neigmx,natom,num_atcell,num_ppcell,num_ppcell_d,nperi,npopshow,     & ! >
   numkmx,numkx,numky,numkz,ksym,kband,nums,ncol,nspv,ncgmin,ncgmax,            & ! >
   ncgres,nprecon_cg,nprecon_diis,nsdmax,nkscg,ndiismax,ncgscf,nretcg,          & ! >
   nrrz,nchange,                                                                & ! >
   looplimit,nbrydn,nmd_start,nmd_end,ngdiis,npolcon,                           & ! >
   npre,nevhist,northo,nradmx,nprjmx,lsphel,lrhomx,                             & ! >
   nfiltyp,nqmx,nmesh,new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz,          & ! >
   nint1dmax,nf,nfdg,npmesh,nfh,npoint,nrc,ncpx,ncpy,ncpz,                      & ! >
   ncpx_d,ncpy_d,ncpz_d,                                                        & ! >
   nwskptot,jelcalc,lmx,nprmx,                                                  & ! >
   xmax,ymax,zmax,socang,gmaxps,epsvh,epssd,ratio_diis,                         & ! >
   eps_scf,eta,etamag,tmstep,sconst,                                            & ! >
   biasx,biasy,biasz,tf,tfmin,tfmax,chrgd,polconocc,endjel,chrjel,              & ! >
   fcut,eps,eps_eig_diis,alambda_diis,alambda_min,alambda_max,psctoff,          & ! >
   psftrad,psctrat,psext,filpp,rctpcc,veta,zs_pre,pol_pre,                      & ! >
   dx,dy,dz,ddx,ddy,ddz,strjel,                                                 & ! >
   cexco,catmfn,lveffout,lcalcpdos,                                             & ! >
   .true.)                                                                        ! <
   if (ndisp == 6) call stopp_initialize (ndisp)
   if ((ndisp /= 6) .and. (ndisp /= 66)) ndisp=nndisp
  call read_input_mpisetup( &
   nxmax,nymax,nzmax,numkmx, & ! <
   numk)                       ! >

! **********  read & distribute atom coordinates, etc.  **********
  call allocate_atom( &
   natom,npolcon,key_polcon_atoms,key_polcon_asa) ! <
  call read_input_atomxyz( &
   chdirinp,                                                    & ! <
   ndisp,catmfn,natom,nsym,lmx,                                 & ! <
   key_sym_any,key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp, & ! <
   key_soc_nocalc,key_soc_calc,                                 & ! <
   xmax,ymax,zmax,                                              & ! <
   num_spe,                                                     & ! >
   numz,nmdx,nmdy,nmdz,natsoc,                                  & ! >
   atx,aty,atz,watom)                                             ! >
  if (myr_kpt==0) then
    if ((npolcon==key_polcon_atoms).or.(npolcon==key_polcon_asa)) then
      call read_input_atommag( &
       chdirinp,                             & ! <
       ndisp,natom,key_polcon2_dir,          & ! <
       npolcon2,polconb,polconmag,polconeta)   ! >
    else
      npolcon2(:)= key_polcon2_none
    end if
  end if !myr_kpt==0
  call mpi_bcast(num_spe,1,mpi_integer,0,mpi_comm_world,mpij)
  if (myrank_glbl<nprocw) call mpi_bcast(natsoc(1,0),natom,mpi_integer,0,mpicom_kpt,mpij)
! **********************************************************

! **********  read & distribute k-point  **********
  call allocate_skp(numk) ! <
  if (myr_kpt==0) then
    call mpi_bcast(numkmx,1,mpi_integer,0,mpicom_space,mpij)
    call allocate_skpxyz(numkmx) ! <
    call read_input_kpt(                  &
     chdirout,numkmx,                     & ! <
     nwskp,skpxyz)                          ! >
    call read_input_kptsnd(numkmx,numk,skpxyz,nwskp, & ! <
                           skpxx,skpyy,skpzz,nwskk)    ! >
  else
    call read_input_kptrcv(numk,                  & ! <
                           skpxx,skpyy,skpzz,nwskk) ! >
  end if !myr_kpt==0
! *************************************************

  call rayleighritz_blacs_ini(neigmx,nrc)

  if (myrank_glbl<nprocw) then

  allocate(mpicom_atom(natom))
  call allocate_arrays
  call allocate_pseudopotentials
  call allocate_read_ppfile
  if (myr_kpt==0) call scf_chargemixing_initialize( &
   ncpx,ncpy,ncpz,nradmx,nspv,nbrydn,npoint,num_atcell) ! <
  ksconv=.false.
  ksitmax=0

! **********  find the species of the atoms  **********
  if (myr_kpt==0) call tools_indspecies( &
   natom,num_spe,numz, & ! <
   mspe,indspe)          ! >
  call mpi_bcast(indspe,natom,mpi_integer,0,mpicom_kpt,mpij)
  call mpi_bcast(mspe,num_spe,mpi_integer,0,mpicom_kpt,mpij)
! ***************************************************

! **********  read pseudopotential database  **********
  if (myr_kpt==0) then
    pi=dacos(-1.0d0)
!    gridmax=dmax1(dx,dy,dz,pi/gmaxps)
    gridmax=pi/gmaxps
    call read_ppfile( &  
     chdirinp,cexco,                              & ! <
     ndisp,num_spe,nradmx,nprmx,lmx,              & ! <
     key_pp_ncps,key_pp_paw,                      & ! <
     mspe,                                        & ! <
     gmaxqp,                                      & ! >
     ntyppp,nradps,nradct,npr,lpmx,nrprj,         & ! >
     cp,radial,dradial,                           & ! >
     potc,awf,pwf,rhocore,rhopcc,grwein,          & ! >
     sss0,akv0,prj,gmax,aeeig,aepot)                ! >
    call setup_pseudopot( &
     num_spe,nradmx,nprmx,nprjmx,lmx,lpmx,npr,nrprj,nradct, & ! <
     radial,dradial,grwein,gmaxqp,gmax,sss0,akv0,awf,pwf,   & ! <
     nwexp,nprj,nlind,noind,sss,akv,rfac,wail,wpil)           ! >
    call filterpseudopotentials( &
     chdirout,                                                                                 & ! <
     ndisp,nfiltyp,num_spe,nqmx,nmesh,nradmx,nprmx,lmx,nradps,nradct,nrprj,mspe,               & ! <
     psftrad,psctrat,psext,psctoff,gridmax,filpp,rctpcc,radial,dradial,rhopcc,potc,prj,cp,     & ! <
     nqct,nqctpcc,coef)                                                                          ! >
  end if ! (myr_kpt==0)
  call setup_pseudopot_broadcast( &
   num_spe,nprmx,nprjmx,lmx,              & ! <
   ntyppp,nprj,nlind,noind,sss,wail,wpil)   ! X
! *****************************************************

! **********  allocate list-vector dependent quantities  **********
  if ((npxmax<0) .or. (npymax<0) .or. (npzmax<0)) then
    if (myrank_glbl==0) write(ndisp,*) 'npxmax, npymax, and npzmax are updated below.'
    if (myrank_glbl==0) write(   12,*) 'npxmax, npymax, and npzmax are updated below.'
  end if
  if (myrank_glbl==0) call tools_listvecdim( &
   natom,num_spe,nradmx,nxmax,nymax,nzmax,ncpx,ncpy,ncpz,nmesh,nfdg, & ! <
   indspe,nradct,dx,dy,dz,psctoff,radial,                            & ! <
   num_list,num_list_d,npxmax,npymax,npzmax)                           ! >
  call mpi_bcast(num_list,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(npxmax,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(npymax,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(npzmax,1,mpi_integer,0,mpi_comm_world,mpij)
  if (myr_kpt==0) call mpi_bcast(num_list_d,1,mpi_integer,0,mpicom_space,mpij)
  call allocate_listvec
  if (myrank_glbl==0) then
    do i=1,2
      if (i==1) jfile= ndisp
      if (i==2) jfile= 12
      write(jfile,*) '==========  data computed by listvecdim  =========='
      write(jfile,*) 'npxmax      :',npxmax
      write(jfile,*) 'npymax      :',npymax
      write(jfile,*) 'npzmax      :',npzmax
      write(jfile,*) '===================================================='
    end do ! i 
  end if
! *****************************************************************

! **********  compute the number of electrons  **********
  if (myr_kpt==0) call tools_countelectron( &
   natom,num_spe,chrgd,indspe,cp, & ! <
   tnumele)                         ! >
! *******************************************************

! **********  check the size of the constrained magnetization  **********
  if ((npolcon/=key_polcon_none).and.(myrank_glbl==0)) call scf_polcon_initialcheck( &
   natom,num_spe,ncol,npolcon,npolcon2,indspe,ntyppp,                                           & ! <
   key_polcon_occ,key_polcon_atoms,key_polcon_asa,key_polcon2_none,key_polcon2_size,key_pp_paw, & ! <
   polconocc,polconmag,tnumele)                                                                   ! <
! ***********************************************************************

! **********  compute weight of spherical integration  **********
  if (myr_kpt==0) then
    call tools_sphglpts(lsphel,npoint,point,wt)
    call tools_genyylm(npoint,lrhomx,point,yylm,ndisp)
    call tools_genyylm(npoint,lrhomx,point,yylm,ndisp)
    call tools_deri_yylm(npoint,lrhomx,point,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2,d2ylm_dtheta_dphi)
  end if ! (myr_kpt==0)
! ***************************************************************

! **********  read or prepare initial wavefunctions and electron density  **********
  l_newdens= (nevhist/=-1)
  call get_wfpotdir( &
   ndisp,                     & ! <
   scfctr2,chdirinp,chdirout)   ! >
  if (npre==0) then
    call read_wf( &
     -1,-1,trim(chdirinp),nrc,nums,ncol,neigmx,numk,ncpx,ncpy,ncpz,ksym, & ! <
     key_ksym_inv2,                                                      & ! <
     svecre,sveccm,sval)                                                   ! >
    if (myr_kpt==0) then
      call read_vht( &
       -1,trim(chdirinp),nevhist,ncpx,ncpy,ncpz,nmesh, & ! <
       vh_dense)                                         ! >
      call read_rho( &
       -1,trim(chdirinp),nevhist,nspv,ncpx,ncpy,ncpz,nradmx,npoint,num_atcell, & ! <
       lbrydn,n,natcell_inf,rhosmt_i,rhotrur_i,rhosmtr_i)                        ! >
      if ((npolcon2(0)==key_polcon2_rot).or.(npolcon==key_polcon_asa)) then
        natpri(:)= key_natpri_out ! here, it might also be key_natpri_inps
        do i= 1,n
          natpri(natcell_inf(i))= key_natpri_in
          natpri_inf(natcell_inf(i))= i
        end do
      end if 
      if (npolcon2(0)==key_polcon2_rot) call tools_maginitrot( &
       key_natpri_in,nspv,natom,num_spe,num_atcell,nradmx,npoint,ncpx,ncpy,ncpz, & ! <
       natpri,natpri_inf,indspe,nradct,radial,dradial,wt,                        & ! <
       polconmag,                                                                & ! <
       rhosmt_i,rhotrur_i,rhosmtr_i)                                               ! X
      if (npolcon==key_polcon_asa) then
        call scf_polcon_magproj( &
         natom,num_spe,num_atcell,nspv,nradmx,npoint,natpri,natpri_inf,indspe,nradct,npolcon2, & ! <
         key_natpri_in,key_polcon2_none,                                                       & ! <
         polconmag,                                                                            & ! <
         rhotrur_i)                                                                              ! X
      end if
    end if 
    if ((nso==1).and.(ncol==1).and.((nevhist==1).or.(nevhist==3))) then 
      if (myr_kpt==0) then 
        call read_rho( &
         -1,trim(chdirinp),4,nspv,ncpx,ncpy,ncpz,nradmx,npoint,num_atcell, & ! <
         i,n,natcell_inf,rhosmt_o,rhotrur_o,rhosmtr_o)                       ! >
         l_newdens= (i==0).and.(lbrydn==0)
         if ((i==0).and.(lbrydn==1)) then
           rhosmt_o = rhosmt_i
           rhotrur_o= rhotrur_i
           rhosmtr_o= rhosmtr_i
         end if 
      end if 
      call mpi_bcast(l_newdens,1,mpi_logical,0,mpi_comm_world,mpij)
    end if
  else
    if (myr_kpt/=0) allocate(cp(8,num_spe))
    call setup_initialwave( &
      nrc,natom,num_spe,nums,ncol,numk,neigmx,ncpx,ncpy,ncpz,zs_pre,pol_pre, & ! <
      xmax,ymax,zmax,dx,dy,dz,                                               & ! <
      mspe,indspe,                                                           & ! <
      cp,atx,aty,atz,                                                        & ! X
      svecre,sveccm,sval)                                                      ! >
    if (myr_kpt/=0) deallocate(cp)
    if (npre .eq. 2) then
      call write_wf( &
       -1,-1,trim(chdirout),nrc,nums,ncol,neigmx,numk,ncpx,ncpy,ncpz,-1,-1,-1,ksym, & ! <
       svecre,sveccm,sval)                                                            ! <
      goto 2000
    end if 
    if (myr_kpt==0) then
      lbrydn= 0
      vh_dense(:,:,:)=0.0d0
    end if ! (myr_kpt==0)
  end if
  residual_states(:,:,:)= 0.0d0 
! **********************************************************************************

! **********  display computational conditions  **********
  if (myrank_glbl==0) call output_compcond( &
   ndisp,12,nperi,numkmx,nmesh,nums,nspv,npre,                            & ! < 
   xmax,ymax,zmax,dx,dy,dz,skpxyz,nwskp,                                  & ! < 
   tnumele,chrgd,biasx,biasy,biasz,tmstep,tf,fcut,eta,etamag)   ! <
! ********************************************************

! **********  display input for B_con  **********
  if ( (myrank_glbl==0) .and. ((npolcon==key_polcon_atoms).or.(npolcon==key_polcon_asa)) ) call output_polcon( &
   1,ndisp,12,key_polcon2_new,key_polcon2_rot,key_polcon2_none,key_polcon2_dir,key_polcon2_size,key_polcon2_fix, & ! <
   natom,npolcon2,polconb,polconmag,polconeta)                                                                     ! <
! ***********************************************
  
! **********  update npolcon2(0) after output_polcon  **********
  if (myr_kpt==0) then 
    if (npolcon2(0)==key_polcon2_rot) npolcon2(0)= key_polcon2_new  
  end if 
! **************************************************************

  call starttime(stime0)

! **********  compute jellium potential on real-space grid  **********
  if (myr_kpt==0) then
    call tools_jelliumpotential( &
     key_jel_calc,jelcalc,ncpz,nzmax,xmax,ymax,zmax,chrjel,endjel,strjel, & ! <
     vjell)                                                                 ! >
  end if ! (myr_kpt==0)
! ********************************************************************

! **********  compute occupation number  **********
  call scf_occupation( &
   ndim_ke,npre,npolcon,ndisp,0,nso,nums,ncol,neigmx,numkmx,numk,nwskp,nwskptot, & ! <
   key_polcon_occ,                                                               & ! <
   tf,tfmax,tfmin,tnumele,polconocc,sval,residual_states,                        & ! <
   ksconv,ksitmax,                                                               & ! X
   ferm,sval_wfc,residual_states_wfc,fnele_wfc,fnele)                              ! >
! *************************************************

! __________  start of structural optimization  __________
  imd=nmd_start
1040 continue

  imd=imd+1

  if (myr_kpt==0) then

    call tools_shiftatoms( &
     nperi,natom,xmax,ymax,zmax,dx,dy,dz, & ! <
     atx,aty,atz)                           ! X

  call starttime(stime)

! **********  compute fuzzy cell weight  **********
    if ((nperi .lt. 3) .or. (npopshow .eq. 1)) call fuzzycell( &
     natom,nperi,ncpx,ncpy,ncpz,xmax,ymax,zmax,dx,dy,dz,atx,aty,atz, & ! <
     pwei)                                                             ! >
! *************************************************

  end if ! (myr_kpt==0)

! **********  compute pseudopotentials on grid points  **********
  if (myr_kpt==0) &
    call pseudocalc(natom,num_spe,nperi,nfdg,npmesh,nmesh,nf,nfh,ndisp,lrhomx,npoint,lsphel,nfiltyp,nqmx,nprjmx,nradmx,nprmx,lmx, & ! <
                    jelcalc,nint1dmax,                                                                                     & ! <
                    ncpx,ncpy,ncpz,ncpx_d,ncpy_d,ncpz_d,npxmax,npymax,npzmax,nxmax,nymax,nzmax,                            & ! <
                    new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz,                                                       & ! <
                    num_atcell,num_ppcell,num_ppcell_d,num_list,num_list_d,                                                & ! <
                    key_pp_paw,key_natpri_in,key_natpri_inps,key_natpri_out,key_jel_calc,                                  & ! <
                    eps,veta,psctoff,psftrad,filpp,rctpcc,chrjel,strjel,endjel,                                            & ! <
                    xmax,ymax,zmax,biasx,biasy,biasz,                                                                      & ! <
                    indspe,nlind,noind,nqct,nqctpcc,ntyppp,nprj,lpmx,nradct,                                               & ! <
                    coef,cp,radial,dradial,potc,awf,pwf,                                                                   & ! <
                    yylm,point,wt,                                                                                         & ! <
                    atx,aty,atz,                                                                                           & ! <
                    natpri,natprid,natpri_inf,natinf,natinfd,natinfd_vloc,naps,napsd,lstvec2,lstvecd2,                     & ! >
                    natx,naty,natz,ndatx,ndaty,ndatz,lstx,lsty,lstz,lstdx,lstdy,lstdz,                                     & ! >
                    qijl,vloc_coarse,vloc_dense,vcorer_all,vcorer,vnlocp,rhopcc_dense,rhopccr,                             & ! >
                    vloc_scw,vloc_hdp,dvlocdx_scw,dvlocdy_scw,dvlocdz_scw,dvlocdx_hdp,dvlocdy_hdp,dvlocdz_hdp,             & ! >
                    vboundx,vboundy,vboundz,vjell)                                                                           ! >
  call pseudocalc_broadcast( &
   nprjmx,natom,num_list,num_ppcell,                                           & ! <
   natpri,natpri_inf,naps,natinf,natx,naty,natz,lstx,lsty,lstz,lstvec2,vnlocp)   ! X
! ***************************************************************

! **********  construct communicators for integration of NLP of pp  **********
  lspec=.false.
  do iz=1,2
  do iy=1,2
  do ix=1,2
    if ((nprocx>=4) .and. (nprocy>=4) .and. (nprocz>=4) .and. &
        (2*npxmax<=ix*ncpx) .and. (2*npymax<=iy*ncpy) .and. (2*npzmax<=iz*ncpz) .and. &
        (mod(nprocx,2*ix)==0) .and. (mod(nprocy,2*iy)==0) .and. (mod(nprocz,2*iz)==0) .and. (.not. lspec)) then
      if (myrank_glbl == 0) write(ndisp,'(a,3i1,a)') 'special routine for mpi_communicator ',2*ix,2*iy,2*iz,' is used. See kukan.'
      call tools_nlpcomsplit_nnn(.true.,natom,natx,naty,natz,ncpx,ncpy,ncpz,nxmax,nymax,nzmax,2*ix,2*iy,2*iz)
      lspec=.true.
    end if
  end do
  end do
  end do
  if (.not. lspec) then
    if (myrank_glbl == 0) write(ndisp,'(a,3i1,a)') 'mpi_communicator is assigned to each atom. See kukan.'
    call tools_nlpcomsplit(.true.,natom,key_natpri_out,natpri)
  end if
! ****************************************************************************
! **********  construct index array for integration of NLP of pp  **********
  call tools_rodr_sum_nlp(natom,key_natpri_out,natpri,latom)
! **************************************************************************
  if (myrank_glbl == 0) write (ndisp,*) '==========  nonlocal pp integration order  =========='
  do i=1,natom
     if (myrank_glbl == 0) write (ndisp,*) i,latom(i)
  end do
  if (myrank_glbl == 0) write (ndisp,*) '====================================================='
  call endtime(ndisp,stime,'[TI] atomic potentials :')

  call starttime(stime)

! **********  orthonormalize wave function (svec)  **********
  if (nrc==0) then
    call orthogonalization_r(natom,neigmx,nprjmx,nums,numk,num_ppcell,num_spe,northo,num_list, & ! <
                             key_ortho_cmpt_innerproduct,ndisp,                                & ! <
                             ncpx,ncpy,ncpz,                                                   & ! <
                             key_ortho_cmpt_innerproduct,                                      & ! <
                             key_natpri_in,key_natpri_inps,                                    & ! <
                             dx,dy,dz,                                                         & ! <
                             nprj,indspe,natpri,naps,natinf,lstvec2,latom,ntyppp,              & ! <
                             vnlocp,sss,                                                       & ! <
                             svecre,ssvre,rspsep)                                                ! X
  else
    call orthogonalization_c(natom,neigmx,nprjmx,nums,ncol,numk,num_ppcell,num_spe,northo,num_list, & ! <
                             key_ortho_cmpt_innerproduct,ndisp,                                     & ! <
                             ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                                   & ! <
                             key_ortho_cmpt_innerproduct,                                           & ! <
                             key_natpri_in,key_natpri_inps,                                         & ! <
                             dx,dy,dz,                                                              & ! <
                             nprj,indspe,natpri,naps,natinf,lstvec2,latom,ntyppp,                   & ! <
                             lstx,lsty,lstz,natx,naty,natz,                                         & ! <
                             skpxx,skpyy,skpzz,                                                     & ! <
                             vnlocp,sss,                                                            & ! <
                             sveccm,ssvcm,cspsep)                                                     ! X
  end if
! ***********************************************************
  call endtime(ndisp,stime,'[TI] GS for ini. wf. :')

  call starttime(stime)
! **********  compute charge density and spin polarization from input wf  *******************************
! rspsep and cspsep are computed above in the routine for orthonormalization of WFs
  if (l_newdens) then
    call scf_charge( &
     key_natpri_in,key_pp_paw,key_sym_any,key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp,                                   & ! <
     ndim_ke,nrc,nsym,numk,neigmx,nums,ncol,nspv,ncpx,ncpy,ncpz,num_spe,natom,num_atcell,num_ppcell,nprjmx,nprmx,lmx,lrhomx, & ! <
     nradmx,npoint,natpri,naps,                                                                                              & ! <
     natpri_inf,indspe,ntyppp,nprj,nlind,noind,nradct,dx,dy,dz,yylm,wt,awf,pwf,radial,dradial,                               & ! <
     fnele,   svecre,sveccm, rspsep,cspsep,                                                                                  & ! <
     atocc,rhosmt_o,rhotrur_o,rhosmtr_o,spinpol)                                                                               ! >
    if ((myr_kpt==0).and.(npolcon==key_polcon_asa)) call scf_polcon_magproj( &
     natom,num_spe,num_atcell,nspv,nradmx,npoint,natpri,natpri_inf,indspe,nradct,npolcon2, & ! <
     key_natpri_in,key_polcon2_none,                                                       & ! <
     polconmag,                                                                            & ! <
     rhotrur_o)                                                                              ! X
  end if
  l_newdens= (nevhist/=-1)
! *******************************************************************************************************
  call endtime(ndisp,stime,'[TI] setup charge :')

! __________  start of s.c.f  __________
  iscf=0
1020 fname='count.dat'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  if (myrank_glbl==0) open (11,file=fname)
  iscf=iscf+1

  call starttime(stime)
  if (myr_kpt==0) then

! **********  mix the charge densities  **********
    call scf_chargemixing( &
     (nevhist==-1),                                              & ! <
     key_natpri_in,key_pp_paw,                                   & ! <
     nbrydn,ndisp,                                               & ! <
     natom,num_spe,num_atcell,nradmx,nspv,npoint,ncpx,ncpy,ncpz, & ! <
     ntyppp,nradct,indspe,natpri,natpri_inf,                     & ! <
     eta,etamag,dx,dy,dz,radial,dradial,wt,                      & ! <
     rhosmt_o,rhotrur_o,rhosmtr_o,                               & ! <
     lbrydn,                                                     & ! X
     rhosmt_i,rhotrur_i,rhosmtr_i)                                 ! X
! ************************************************

! **********  compute charge moments and potential **********
    call scf_potentials( &
     key_natpri_in,key_natpri_inps,key_pp_paw,                                        & ! <
     nmesh,ncpx,ncpy,ncpz,natom,num_atcell,num_spe,nradmx,npoint,lrhomx,lmx,          & ! <
     num_list_d,num_ppcell_d,nspv,                                                    & ! <
     npopshow,nevhist,ndisp,nperi,nf,nfh,ncgres,ncgmin,ncgmax,                        & ! <
     indspe,natpri,natprid,natpri_inf,napsd,natinfd,ndatx,ndaty,ndatz,                & ! <
     ntyppp,nradct,lpmx,nwexp,lstvecd2,lstdx,lstdy,lstdz,cexco,                       & ! <
     epsvh,xmax,ymax,zmax,dx,dy,dz,tnumele,biasx,biasy,biasz,                         & ! <
     yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2,d2ylm_dtheta_dphi,          & ! <
     wt,point,radial,dradial,cp,rfac,atx,aty,atz,pwei,                                & ! <
     rhocore,rhopccr,rhotrur_i,rhosmtr_i,rhopcc_dense,rhosmt_i,vloc_coarse,           & ! <
     vboundx,vboundy,vboundz,                                                         & ! <
     atmpole,rhotrucorer,rhosmt_pcc_dense,rhosmt_pccr,rhoaugr,rhoaug3d,rho_aug_dense, & ! >
     veff(1,1,1,1),vh_coarse,vhtrur,vhsmtr,vhaugr,vx,vx_dense,vxctru,vxcsmt,          & ! >
     ex_dense,exctru,excsmt,                                                          & ! >
     vh_dense)                                                                          ! X
! **********************************************************

! **********  compute constraining field **********
    if ((npolcon==key_polcon_atoms).or.(npolcon==key_polcon_asa)) then
      if (npolcon2(0)==key_polcon2_new) call scf_polcon( &
       natom,nspv,natpri,npolcon2,key_natpri_in,key_polcon2_dir,key_polcon2_size, & ! <
       polconeta,polconmag,spinpol,                                               & ! <
       polconb)                                                                     ! X
    end if
! *************************************************

! **********  compute D_{ij}  **********
    i= 1
    if ((npolcon==key_polcon_atoms).or.(npolcon==key_polcon_asa)) i= natom
    call scf_dij( &
     natom,i,num_spe,num_atcell,nradmx,nprjmx,nums,ncol,nspv,npoint, & ! <
     lrhomx,num_list_d,num_ppcell_d,lmx,nprmx,                       & ! <
     ncpx_d,ncpy_d,ncpz_d,                                           & ! <
     npolcon,                                                        & ! <
     key_natpri_in,key_natpri_inps,key_pp_paw,                       & ! <
     key_polcon_atoms,key_polcon_asa,key_polcon2_none,               & ! <
     ddx,ddy,ddz,                                                    & ! <
     npolcon2,ntyppp,nradct,nprj,npr,indspe,natpri,natpri_inf,       & ! <
     nwexp,natinfd,napsd,natprid,lpmx,                               & ! <
     nlind,noind,lstvecd2,                                           & ! <
     lstdx,lstdy,lstdz,ndatx,ndaty,ndatz,                            & ! <
     awf,pwf,wail,radial,dradial,rfac,                               & ! <
     atx,aty,atz,                                                    & ! <
     vh_dense,vloc_dense,                                            & ! <
     vhtrur,vhsmtr,vhaugr,                                           & ! <
     vxctru,vxcsmt,vcorer_all,vcorer,qijl,akv,polconb,yylm,wt,       & ! <
     dij)                                               ! >
    if (nso==1) call scf_dij_soc( &
     natom,num_atcell,num_spe,nprmx,nprjmx,lmx,nradmx,npoint,nspv,            & ! <
     key_natpri_in,key_pp_paw,key_soc_calc,                                   & ! <
     natpri,natpri_inf,indspe,natsoc,mspe,ntyppp,nradct,nprj,npr,nlind,noind, & ! <
     socang,radial,dradial,wt,awf,aeeig,aepot,                                & ! <
     rhotrucorer,vxctru,                                                      & ! <
     dijsoc)                                                                    ! >
! **************************************

  end if ! (myr_kpt==0)

! **********  broadcast V_eff and D_{ij} within mpicom_kpt  ***
  call tools_potbroadcast( &
   nso,nums,ncol,nspv,ncpx,ncpy,ncpz,natom,nprjmx,                  & ! <
   veff,dij,dijsoc)                                                   ! X
! ************************************************************
  call endtime(ndisp,stime,'[TI] potentials & Dij :')


  call starttime(stime)
! **********  solve eigen problem  **********
  call scf_diag( &
   nrc,nso,nprjmx,num_ppcell,num_spe,iscf,ncgscf,nretcg,                 & ! < 
   natom,neigmx,nums,ncol,nspv,numk,nperi,northo,nf,ndisp,num_list,      & ! <
   ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                                  & ! <
   nchange,nrrz,                                                         & ! <
   key_ortho_nocmpt_innerproduct,key_ortho_cmpt_innerproduct,            & ! <
   key_natpri_in,key_natpri_inps,key_natpri_out,key_soc_calc,key_pp_paw, & ! <
   skpxx,skpyy,skpzz,dx,dy,dz,xmax,ymax,zmax,                            & ! <
   dij,dijsoc,veff,vnlocp,sss,                                           & ! <
   nkscg,nsdmax,nprecon_cg,                                              & ! <  ! only cg
   epssd,                                                                & ! <  !   "
   ndiismax,nprecon_diis,                                                & ! <  ! only diis
   eps_eig_diis,alambda_diis,ratio_diis,alambda_min,alambda_max,         & ! <  !   "
   indspe,natpri,naps,nprj,natinf,lstvec2,latom,natsoc(1,0),ntyppp,      & ! <
   lstx,lsty,lstz,natx,naty,natz,                                        & ! <
   rspsep,cspsep,svecre,sveccm,hsvre,hsvcm,ssvre,ssvcm,                  & ! X
   sval,svalsoc,cspsepsoc,svecsoc,                                       & ! >
   residual_states,ksconv,ksitmax)                                         ! >
! *******************************************
  call endtime(ndisp,stime,'[TI] eigenproblem :')

  call starttime(stime)
! **********  compute occupation number  **********
  if ((nso==0).or.(ncol==2)) then
    call scf_occupation( &
     ndim_ke,0,npolcon,ndisp,0,nso,nums,ncol,neigmx,numkmx,numk,nwskp,nwskptot, & ! <
     key_polcon_occ,                                                            & ! <
     tf,tfmax,tfmin,tnumele,polconocc,sval   ,residual_states,                  & ! <
     ksconv,ksitmax,                                                            & ! X
     ferm,sval_wfc,residual_states_wfc,fnele_wfc,fnele   )                        ! >
  else
    call scf_occupation( &
     ndim_ke,0,npolcon,ndisp,1,nso,nums,ncol,neigmx,numkmx,numk,nwskp,nwskptot, & ! <
     key_polcon_occ,                                                            & ! <
     tf,tfmax,tfmin,tnumele,polconocc,svalsoc,residual_states,                  & ! <
     ksconv,ksitmax,                                                            & ! X
     ferm,sval_wfc,residual_states_wfc,fnele_wfc,fnelesoc)                        ! >
  end if
! *************************************************

! **********  compute new charge density  ***************************************************************
  if ((nso==0).or.(ncol==2)) then
    call scf_charge( &
     key_natpri_in,key_pp_paw,key_sym_any,key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp,                                     & ! <
     ndim_ke,nrc,nsym,numk,  neigmx,nums,ncol,nspv,ncpx,ncpy,ncpz,num_spe,natom,num_atcell,num_ppcell,nprjmx,nprmx,lmx,lrhomx, & ! <
     nradmx,npoint,natpri,naps,                                                                                                & ! <
     natpri_inf,indspe,ntyppp,nprj,nlind,noind,nradct,dx,dy,dz,yylm,wt,awf,pwf,radial,dradial,                                 & ! <
     fnele,   svecre,sveccm, rspsep,   cspsep,                                                                                 & ! <
     atocc,rhosmt_o,rhotrur_o,rhosmtr_o,spinpol)                                                                                 ! >
  else
    call scf_charge( &
     key_natpri_in,key_pp_paw,key_sym_any,key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp,                                     & ! <
     ndim_ke,1,  nsym,numk,2*neigmx,   2,   2,nspv,ncpx,ncpy,ncpz,num_spe,natom,num_atcell,num_ppcell,nprjmx,nprmx,lmx,lrhomx, & ! <
     nradmx,npoint,natpri,naps,                                                                                                & ! <
     natpri_inf,indspe,ntyppp,nprj,nlind,noind,nradct,dx,dy,dz,yylm,wt,awf,pwf,radial,dradial,                                 & ! <
     fnelesoc,svecre,svecsoc,rspsep,cspsepsoc,                                                                                 & ! <
     atocc,rhosmt_o,rhotrur_o,rhosmtr_o,spinpol)                                                                                 ! >
  end if 
  if ((myr_kpt==0).and.(npolcon==key_polcon_asa)) call scf_polcon_magproj( &
   natom,num_spe,num_atcell,nspv,nradmx,npoint,natpri,natpri_inf,indspe,nradct,npolcon2, & ! <
   key_natpri_in,key_polcon2_none,                                                       & ! <
   polconmag,                                                                            & ! <
   rhotrur_o)                                                                              ! X
! *******************************************************************************************************

! **********  compute difference in charge densities for SCF iteration  **********
  if (myr_kpt==0) then
    call scf_diffcharge(ncpx,ncpy,ncpz,nspv,rhosmt_o,rhosmt_i,diff_charge)
    diff_charge=dsqrt(diff_charge*dx*dy*dz/(8.0d0*xmax*ymax*zmax))
  end if ! (myr_kpt==0)
! ********************************************************************************
  call endtime(ndisp,stime,'[TI] occ. #, new charge & diff charge :')

  if (myrank_glbl==0) l_scf= (diff_charge < eps_scf) .or. (iscf>=looplimit)
  call mpi_bcast(l_scf,1,mpi_logical,0,mpi_comm_world,mpij)

! =========================================

! ****** display s.c.f. output ***************************************

  if (myrank_glbl==0) then
! **********  write 'count.dat'  **********
    i= 1
    if ((npolcon==key_polcon_atoms).or.(npolcon==key_polcon_asa)) i= natom
    call output_countdat( &
     ndisp,11,iscf,nperi,npopshow,nso,nspv,nums,ncol,natom,i,numkmx,neigmx,                              & ! <
     npolcon,npolcon2,key_polcon_atoms,key_polcon_asa,key_polcon2_none,                                      & ! <
     nwskp,nwskptot,                                                                                         & ! <
     ksconv,ksitmax,tnumele,ferm,spinpol,polconb,atmpole,diff_charge,residual_states_wfc,sval_wfc,fnele_wfc)   ! <
    close (11)
    call flush (ndisp)
  end if ! (myrank_glbl==0)
! *****************************************

  scfctr(1)= l_scf .and. (imd >= nmd_end)
  scfctr(3)= l_scf
  call tools_scffile(chdirout,scfctr)

  if ( (scfctr(1) .and. scfctr2) .or. scfctr(2) ) then
! **********  output of wf&pot files, B_con  **********
    call write_wf( &
     -1,-1,trim(chdirout),nrc,nums,ncol,neigmx,numk,ncpx,ncpy,ncpz,-1,-1,-1,ksym, & ! <
     svecre,sveccm,sval)                                                            ! <
    if (myr_kpt==0) then
      i= 1
      if ((nso==1).and.(ncol==1)) i= 2 
      call write_rhovht( &
       -1,trim(chdirout),i,0,npre,nevhist,nspv,ncpx,ncpy,ncpz,nmesh,nradmx,npoint,num_atcell,-1,-1,-1, & ! <
       natom,key_natpri_in,natpri,                                                                     & ! <
       rhosmt_i,rhotrur_i,rhosmtr_i,rhosmt_o,rhotrur_o,rhosmtr_o,vh_dense)                               ! <
    end if 
    if ( (myrank_glbl==0) .and. ((npolcon==key_polcon_atoms).or.(npolcon==key_polcon_asa)) ) call output_polcon( &
     3,ndisp,10,key_polcon2_new,key_polcon2_rot,key_polcon2_none,key_polcon2_dir,key_polcon2_size,key_polcon2_fix, & ! <
     natom,npolcon2,polconb,polconmag,polconeta)                                                                     ! <
! *****************************************************
  end if

  if (scfctr(3)) then

  call starttime(stime)
! **********  compute total energy  **********
    if (myr_kpt==0) then
      i= 1
      if ((npolcon==key_polcon_atoms).or.(npolcon==key_polcon_asa)) i= natom
      call totalenergy_space( &
       key_natpri_in,key_pp_paw,key_jel_calc,                                    & ! <
       key_polcon_atoms,key_polcon_asa,key_polcon2_size,                         & ! <
       i,natom,num_spe,nperi,nmesh,nspv,npoint,ncpx,ncpy,ncpz,num_atcell,nradmx, & ! <
       jelcalc,nint1dmax,nzmax,new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz,  & ! <
       natpri,natpri_inf,indspe,ntyppp,nradct,npolcon,npolcon2,                  & ! <
       wt,radial,dradial,cp,veta,tnumele,dx,dy,dz,xmax,ymax,zmax,                & ! <
       biasx,biasy,biasz,chrjel,endjel,strjel,atx,aty,atz,                       & ! <
       rhosmt_i,rho_aug_dense,rhosmt_pcc_dense,                                  & ! <
       rhotrur_i,rhosmtr_i,rhoaugr,rhosmt_pccr,rhotrucorer,spinpol,polconb,      & ! <
       vh_dense,vx,vhtrur,vhsmtr,vhaugr,vxctru,vxcsmt,ex_dense,exctru,excsmt,    & ! <
       eneco,eneha,eneex,eneat,eneof,enejel,enebc)                                 ! >
    end if 
    if (myr_space==0) then
      if ((nso==0).or.(ncol==2)) then
        call totalenergy_eig( &
         nums,ncol,numk,  neigmx,nwskk,nwskptot,tf,   sval,   fnele, & ! <
         eneel,eneth)                                                  ! >
      else
        call totalenergy_eig( &
         2,2      ,numk,2*neigmx,nwskk,nwskptot,tf,svalsoc,fnelesoc, & ! <
         eneel,eneth)                                                  ! >
      end if 
    end if 
! ********************************************
  call endtime(ndisp,stime,'[TI] total energy :')

! **********  write energy  **********
    if (myrank_glbl==0 .and. sconst>eps) call tools_calenespr(chdirinp,natom,xmax,ymax,zmax,atx,aty,atz,sconst,enespr)
    if (myrank_glbl==0) call output_energy( &
     chdirout,                                                                                   & ! <
     ndisp,12,iscf,nperi,npopshow,nso,nspv,nums,ncol,natom,numkmx,neigmx,nwskp,nwskptot,kband,   & ! <
     tnumele,ferm,spinpol,atmpole,diff_charge,residual_states_wfc,sval_wfc,fnele_wfc,sconst,eps, & ! <
     eneel,eneex,eneha,eneat,eneco,eneth,eneof,enejel,enebc,enespr)                                ! <
! ************************************

! **********  output KS-effective potential  **********
  if (lveffout) then
  if (myr_kpt == 0) call potout(chdirout,ndisp,nspv,nums,ncpx,ncpy,ncpz,myrank_glbl,nprocx,nprocy,nprocz, & !<
                                mpicom_space,mpi_status_size,mpi_double_precision,mpistat,  & !<
                                ferm,veff)                                                    !<
  end if
! *****************************************************

  call starttime(stime)
! **********  force routines  **********

    if (ncol==2) then
      ! ?? To do!
      if (myrank_glbl==0) write(ndisp,*) 'Forces not implemented for noncollinear magnetism !'
      goto 1011
    end if

  call starttime(stimesub)
! **********  compute charge moments and potential from output charge **********
    if (myr_kpt==0) call scf_potentials( &
     key_natpri_in,key_natpri_inps,key_pp_paw,                                        & ! <
     nmesh,ncpx,ncpy,ncpz,natom,num_atcell,num_spe,nradmx,npoint,lrhomx,lmx,          & ! <
     num_list_d,num_ppcell_d,nspv,                                                    & ! <
     npopshow,nevhist,ndisp,nperi,nf,nfh,ncgres,ncgmin,ncgmax,                        & ! <
     indspe,natpri,natprid,natpri_inf,napsd,natinfd,ndatx,ndaty,ndatz,                & ! <
     ntyppp,nradct,lpmx,nwexp,lstvecd2,lstdx,lstdy,lstdz,cexco,                       & ! <
     epsvh,xmax,ymax,zmax,dx,dy,dz,tnumele,biasx,biasy,biasz,                         & ! <
     yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2,d2ylm_dtheta_dphi,          & ! <
     wt,point,radial,dradial,cp,rfac,atx,aty,atz,pwei,                                & ! <
     rhocore,rhopccr,rhotrur_o,rhosmtr_o,rhopcc_dense,rhosmt_o,vloc_coarse,           & ! <
     vboundx,vboundy,vboundz,                                                         & ! <
     atmpole,rhotrucorer,rhosmt_pcc_dense,rhosmt_pccr,rhoaugr,rhoaug3d,rho_aug_dense, & ! >
     veff(1,1,1,1),vh_coarse,vhtrur,vhsmtr,vhaugr,vx,vx_dense,vxctru,vxcsmt,          & ! >
     ex_dense,exctru,excsmt,                                                          & ! >
     vh_dense)                                                                          ! X
! ******************************************************************************
  call endtime(ndisp,stimesub,'[TI_sub] scf_potentials :')

  call starttime(stimesub)
! **********  broadcast V_eff and D_{ij} within mpicom_kpt  ***
    call tools_potbroadcast( &
     nso,nums,ncol,nspv,ncpx,ncpy,ncpz,natom,nprjmx,                  & ! <
     veff,dij,dijsoc)                                                   ! X
! ************************************************************
  call endtime(ndisp,stimesub,'[TI_sub] tools_potbroadcast :')

! **********  compute force acting on atoms  **********
  call starttime(stimesub)
    if (nrc==0) then
      call force_eig_r( &
       nperi,numk,nums,neigmx,nf,ncpx,ncpy,ncpz,num_list,natom,num_spe,num_ppcell,nprjmx, & ! <
       key_natpri_in,key_natpri_inps,                                                     & ! <
       lstvec2,latom,lstx,lsty,lstz,nprj,natpri,naps,indspe,natinf,                       & ! <
       xmax,ymax,zmax,sss,vnlocp,dij,                                                     & ! <
       svecre,sval,fnele,                                                                 & ! <
       fatx,faty,fatz)                                                                      ! >
    else
      call force_eig_c( &
       nperi,numk,nums,ncol,neigmx,nf,ncpx,ncpy,ncpz,num_list,natom,num_spe,num_ppcell,nprjmx, & ! <
       key_natpri_in,key_natpri_inps,                                                          & ! <
       lstvec2,latom,lstx,lsty,lstz,nprj,natpri,naps,indspe,natinf,                            & ! <
       xmax,ymax,zmax,sss,vnlocp,dij,                                                          & ! <
       natx,naty,natz,skpxx,skpyy,skpzz,                                                       & ! <
       sveccm,sval,fnele,                                                                      & ! <
       fatx,faty,fatz)                                                                           ! >
    end if
  call endtime(ndisp,stimesub,'[TI_sub] force_eig :')
    if (myr_kpt==0) then
  call starttime(stimesub)
      call scf_augcharge( &
       1,key_natpri_in,key_natpri_inps,key_pp_paw,                       & ! <
       nspv,nradmx,npoint,lrhomx,lmx,natom,num_atcell,num_spe,           & ! <
       num_list_d,num_ppcell_d,                                          & ! <
       indspe,natpri,natprid,natpri_inf,napsd,natinfd,ndatx,ndaty,ndatz, & ! <
       ntyppp,nradct,lpmx,nwexp,lstdx,lstdy,lstdz,                       & ! <
       ddx,ddy,ddz,yylm,wt,radial,dradial,rfac,                          & ! <
       atx,aty,atz,rhotrur_o,rhosmtr_o,                                  & ! <
       rhoaugr,rhoaug3d, drhoaug3ddx,drhoaug3ddy,drhoaug3ddz)              ! >
  call endtime(ndisp,stimesub,'[TI_sub] scf_augcharge :')
  call starttime(stimesub)
      call trans_d2c_smtcharge( &
       nmesh,nperi,ndisp,nf,ncpx,ncpy,ncpz, & ! <
       rho_aug_dense,                       & ! <
       rho_coarse)                            ! >
  call endtime(ndisp,stimesub,'[TI_sub] trans_d2c_smtcharge :')
  call starttime(stimesub)
      call force(natom,num_spe,num_atcell,num_ppcell_d,num_list_d,nqmx,nradmx,npoint,                 & ! <
                 nmesh,nspv,nperi,nint1dmax,nzmax,jelcalc,                                            & ! <
                 ncpx,ncpy,ncpz,ncpx_d,ncpy_d,ncpz_d,                                                 & ! <
                 new_pwx,new_pwy,new_pwz,new_rsx,new_rsy,new_rsz,                                     & ! <
                 key_jel_calc,key_natpri_in,key_natpri_inps,                                          & ! <
                 indspe,natpri,natpri_inf,nradct,natprid,napsd,natinfd,natinfd_vloc,lstvecd2,nqctpcc, & ! <
                 lstdx,lstdy,lstdz,                                                                   & ! <
                 veta,psctoff,psftrad,radial,dradial,point,wt,xmax,ymax,zmax,biasx,biasy,biasz,       & ! <
                 cp,coef,                                                                             & ! <
                 atx,aty,atz,                                                                         & ! <
                 dvlocdx_scw,dvlocdy_scw,dvlocdz_scw,                                                 & ! <
                 dvlocdx_hdp,dvlocdy_hdp,dvlocdz_hdp,                                                 & ! <
                 rhotrur_o,rhosmtr_o,rhoaugr,drhoaug3ddx,drhoaug3ddy,drhoaug3ddz,                     & ! <
                 rho_coarse,rho_aug_dense,vloc_dense,vh_dense,vx_dense,                               & ! <
                 chrjel,strjel,endjel,                                                                & ! <
                 fatx,faty,fatz,                                                                      & ! X
                 fjel)                                                                                  ! >
    end if
  call endtime(ndisp,stimesub,'[TI_sub] force :')
    if (myrank_glbl==0) call output_force(ndisp,12,natom,fatx,faty,fatz)
    1011 continue
! *****************************************************
  call endtime(ndisp,stime,'[TI] force :')

! **************************************

  end if ! (scfctr(3))
! ************************************

  if ( .not. (l_scf .or. scfctr(1)) )  goto 1020
  ! __________  end of s.c.f  __________

! **********  write old atom coordinate  **********
  if (myrank_glbl==0) call output_molecdyn( &
   chdirout,                                                                           & ! <
   1,ndisp,12,catmfn,natom,num_spe,lmx,indspe,numz,nmdx,nmdy,nmdz,natsoc,ntyppp,watom, & ! <
   xmax,ymax,zmax,imd,atx,aty,atz,fatx,faty,fatz)                                        ! <
! *******************************************************

! **********  move atoms  **********
  if ((myr_kpt==0) .and. (l_scf)) call tools_moveatoms(chdirinp,chdirout,natom,ngdiis,ndisp,watom,nmdx,nmdy,nmdz,xmax,ymax,zmax,tmstep,fcut,sconst,eps, & ! <
                                                             atx,aty,atz,fatx,faty,fatz)                                                ! X
! **********************************

! **********  write new atom coordinate  **********
  if (myrank_glbl==0 .and. sconst>eps) call output_force(ndisp,12,natom,fatx,faty,fatz)
  if (myrank_glbl==0) call output_molecdyn( &
   chdirout,                                                                           & ! <
   2,ndisp,12,catmfn,natom,num_spe,lmx,indspe,numz,nmdx,nmdy,nmdz,natsoc,ntyppp,watom, & ! <
   xmax,ymax,zmax,imd,atx,aty,atz,fatx,faty,fatz)                                        ! <
! *******************************************************

! **********  destruct communicators  **********
  lspec=.false.
  do iz=1,2
  do iy=1,2
  do ix=1,2
    if ((nprocx>=4) .and. (nprocy>=4) .and. (nprocz>=4) .and. &
        (2*npxmax<=ix*ncpx) .and. (2*npymax<=iy*ncpy) .and. (2*npzmax<=iz*ncpz) .and. &
        (mod(nprocx,2*ix)==0) .and. (mod(nprocy,2*iy)==0) .and. (mod(nprocz,2*iz)==0) .and. (.not. lspec)) then
      call tools_nlpcomsplit_nnn(.false.,natom,natx,naty,natz,ncpx,ncpy,ncpz,nxmax,nymax,nzmax,2*ix,2*iy,2*iz)
      lspec=.true.
    end if
  end do
  end do
  end do
  if (.not. lspec) call tools_nlpcomsplit(.false.,natom,key_natpri_out,natpri)
! **********************************************

  if (myrank_glbl==0) call flush (12)

  if ( (imd < nmd_end) .and. (.not. scfctr(1)) ) goto 1040
! __________  end of structural optimization  __________

  call endtime(ndisp,stime0,'real time')
  call endtime(12,stime0,'real time')
  if (myrank_glbl==0) then
    write (ndisp,*,err=9999) 'com. time',comtime,'(sec)'
    write (12   ,*,err=9999) 'com. time',comtime,'(sec)'
    write (ndisp,*,err=9999) 'com. vol.',comcount,'(byte)'
    write (12   ,*,err=9999) 'com. vol.',comcount,'(byte)'
    close(12)
  end if

! **********  write data to atom.xyz  **********
  if (myrank_glbl==0) call output_molecdyn( &
   chdirout,                                                                           & ! <
   3,ndisp,10,catmfn,natom,num_spe,lmx,indspe,numz,nmdx,nmdy,nmdz,natsoc,ntyppp,watom, & ! <
   xmax,ymax,zmax,imd,atx,aty,atz,fatx,faty,fatz)                                        ! <
  if (myrank_glbl==0) call output_xcrysden( &
   chdirout,natom,ndisp,10,atx,aty,atz,fatx,faty,fatz,numz) ! <
! **********************************************

! **********  write pdos  **********
  if (lcalcpdos) &
  call atomicocc_pdos( &
    chdirout,                                                     & ! <
    key_natpri_in,key_pp_paw,nrc,natom,nradmx,num_spe,num_ppcell, & ! <
    nprmx,lmx,nprjmx,nums,ncol,nspv,numk,neigmx,                  & ! <
    natpri,naps,indspe,ntyppp,nprj,nradct,nlind,noind,            & ! <
    awf,pwf,dradial,akv,                                          & ! <
    rspsep,cspsep)                                                  ! <
! **********************************

! **********  write band structure files  **********
  if (kband==key_band_yes) call output_bands( & 
   chdirout,                              & ! <
   nso,nums,ncol,neigmx,numk,numkmx,      & ! <
   sval,svalsoc)                            ! <
! **************************************************

!! **********  output KS-effective potential  **********
!  if (lveffout) then
!  if (myr_kpt == 0) call potout(chdirout,ndisp,nspv,nums,ncpx,ncpy,ncpz,myrank_glbl,nprocx,nprocy,nprocz, & !<
!                                mpicom_space,mpi_status_size,mpi_double_precision,mpistat,  & !<
!                                ferm,veff)                                                    !<
!  end if
!! *****************************************************

 2000 continue

  call deallocate_pseudopotentials
  if (myr_kpt==0) call scf_chargemixing_finalize

  end if ! (myrank_glbl<nprocw)

  if ((myrank_glbl==0).and.(ndisp/=6)) close(nndisp)

  call mpi_barrier(mpi_comm_world,mpij)
  call mpi_comm_free(mpicom_space,mpij)
  if (nprock>1) call mpi_comm_free(mpicom_kpt,mpij)
  call rayleighritz_blacs_fin
  call mpi_barrier(mpi_comm_world,mpij)
  call mpi_finalize(mpij)

stop
  9999 call stopp('error writing output')
end

