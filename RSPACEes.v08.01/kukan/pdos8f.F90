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
! **********  pdos8f.f90 04/18/2023-01  **********

use mod_mpi!, only:nprocx,nprocy,nprocz
use mod_stopp

use var_keys
use var_arrays
use var_read_input_kukan
use var_read_ppfile
use var_pseudopotentials

use dim_arrays

use mod_read_input_kukan
use mod_rw_wfpot

implicit none
integer :: i,nk,ns,ix,iy,iz,ii,ie,l,jz,ktot,na,nsumk
integer :: npxmax,npymax,npzmax
real*8  :: x,y,z
real*8, parameter :: dmin=1.0e-20
real*8 ene_min,ene_max,ferm,alph,de,e,weight,pi,afact
integer nene,numl
logical:: ldirinout
character:: chdirinp*200,chdirout*200,fname*200
!integer mpij,mpi_comm_world,myrank_glbl,nprocs,nprocw,nprock,myr_kpt,myr_space,myrz,myry,myrx

! ** for pdos **
character*6 :: corb(9),ctmpp
integer ij,i1,i2,j1,j2,itmpp
!integer indspe(natom)
!integer,allocatable:: nprj(:),nradct(:),nlind(:,:),noind(:,:),orb_nrm(:,:) ! already defined
real*8,allocatable:: orb_dos(:,:,:,:)
real*8,allocatable:: pdos(:,:,:)

pi=dacos(-1.0d0)

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

  nndisp=66
  fname='output_pdos.dat'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  open(nndisp,file=fname)
  write(6,*) 'computed by pdos8f 11/03/2022-01'
  write(nndisp,*) 'computed by pdos8f 11/03/2022-01'
  call stopp_initialize (nndisp)
  call mpi_comm_size(mpi_comm_world,i,mpij)
  if (i /= 1) call stopp ('error. # of processes defined in script file must be 1.')

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
  nprocs= nprocx*nprocy*nprocz
  nprocw= nprocs*nprock
  ncpx=2*nxmax/nprocx
  ncpy=2*nymax/nprocy
  ncpz=2*nzmax/nprocz

  call read_input_dos(chdirinp,chdirout,ene_min,ene_max,ferm,alph,nene,ndisp)

!**************************************************************************************************

!**************************************************************************************************
! **** routine for PDOS of an atom at the top in atom.xyz ****
! ** input file "input_dos.txt" should be prepared **
! reuse sval, ferm
  corb(1)='   s  '
  corb(2)='  px  '
  corb(3)='  py  '
  corb(4)='  pz  '
  corb(5)='dxx-yy'
  corb(6)=' dzx  '
  corb(7)='dzz-rr'
  corb(8)=' dyz  '
  corb(9)=' dxy  '

  fname='occ_pdos.txt'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  open(31,file=fname)
  i1=0
  numl=-1
  do while (i1 == 0)
    read(31,*,iostat=i1) i
    numl=numl+1
  end do
  numl=numl/(neigmx*numkmx*max(nums,nspv))
  write(ndisp,*) '# of spherical harmonics is', numl
  close(31)

  allocate(orb_dos(numl,neigmx,nums,numkmx))
  allocate(pdos(numl,0:nene,max(nums,nspv)))

  na=1
  orb_dos=0.0d0
  fname='occ_pdos.txt'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  open(31,file=fname)
  do nk=1,numkmx
    do l= 1,neigmx
      do ns= 1,max(nums,nspv)
        do i=1,numl
          read(31,*) itmpp,itmpp,itmpp,ctmpp,orb_dos(i,l,ns,nk)
        end do
      end do ! ns
    end do ! l
  end do ! nk
  close(31)
!**************************************************************************************************
!**************************************************************************************************

  ktot=0
  fname='.kpt_tmp.txt'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  open(10,file=fname)
    do nk=1,numkmx
      read(10,*) x,y,z,ii
      ktot=ktot+ii
    end do
  close(10)

  nsumk=0
  pdos=0.0d0
  do i=0,nprocw-1
    numk= numkmx/nprock
    myr_kpt=i/nprocs
    if (myr_kpt >= nprock-mod(numkmx,nprock) ) numk=numk+1
    myr_space=mod(i,nprocs)
    myrz=myr_space/(nprocx*nprocy)
    myry=(myr_space-myrz*nprocx*nprocy)/nprocx
    myrx=myr_space-myrz*nprocx*nprocy-myry*nprocx
    allocate(svecre(ncpx*(1-nrc)+nrc,ncpy*(1-nrc)+nrc,ncpz*(1-nrc)+nrc,neigmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk*(1-nrc)+nrc))
    allocate(sveccm(ncpx*nrc-nrc+1,ncpy*nrc-nrc+1,ncpz*nrc-nrc+1,neigmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1))
    allocate(sval(neigmx,nums+1-ncol,numk))
    allocate(nwskp(numk))

    fname='.kpt_tmp.txt'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(10,file=fname)
      do nk=1,nsumk
      read(10,*) x,y,z,ii
      end do
      do nk=1,numk
      read(10,*) x,y,z,nwskp(nk)
      end do
    close(10)
    if (myr_space == nprocs-1) nsumk=nsumk+numk
      call read_wf(  &
       myr_space,myr_kpt,'',nrc,nums,ncol,neigmx,numk,ncpx,ncpy,ncpz,ksym, & ! <
       key_ksym_inv2,                                                 & ! <
       svecre,sveccm,sval)                                              ! >

      write(ndisp,*) myr_space,myr_kpt
      de=(ene_max-ene_min)/nene
      afact=1.0d0
      if (nums == 1) afact=2.0d0
!$omp parallel default(shared),private(ie,ns,nk,l,e,weight,ix,iy,iz,jz,ii,ij,i1,j1,i2,j2) 
!$omp do
        do ie=0,nene
          do ns=1,nums
            do nk=1,numk
              do l=1,neigmx
                e=sval(l,ns,nk)-(ferm+(ene_min+ie*de))
                if (alph*e*e > 40.0d0) then
                  weight=0.0d0
                else
                  weight=afact*dsqrt(alph/pi)*dexp(-alph*e*e)*nwskp(nk)/ktot
                end if
                do i1=1,numl
                  pdos(i1,ie,ns)=pdos(i1,ie,ns)+orb_dos(i1,l,ns,nk)*weight
                end do
              end do
            end do
          end do
        end do
!$omp end parallel

    deallocate(sval,svecre,sveccm)
    deallocate(nwskp)
  end do
  fname='.kpt_tmp.txt'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  open(10,file=fname)
  close(10,status='delete')

  if(nums.eq.1) then
    fname='pdos.dat'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(10,file=fname)
    write(10,fmt='(10a18)') "E-EF(eV)",(corb(i1),i1=1,numl)
    do ie=0,nene
      e=ene_min+ie*de
      write(10,fmt='(10E18.8)') e*27.212d0,(pdos(i1,ie,1),i1=1,numl)
    end do
    close(10)
  else if(nums.eq.2) then
    fname='pdos.dat'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(10,file=fname)
    write(10,fmt='(20a18)') "E-EF(eV)",(corb(i1),i1=1,numl),(corb(i1),i1=1,numl)
    do ie=0,nene
      e=ene_min+ie*de
      write(10,fmt='(20E18.8)') e*27.212d0,( pdos(i1,ie,1),i1=1,numl),(-pdos(i1,ie,2),i1=1,numl)
    end do
    close(10)
  end if

  deallocate(pdos,orb_dos)

  close(nndisp)

  call mpi_barrier(mpi_comm_world,mpij)
  call mpi_finalize(mpij)

stop
 2000 call stopp ('error found when reading kpt_tmp.txt')
 2010 call stopp ('error found when removing kpt_tmp.txt')
end
