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
! **********  orbcharge8f.f90 04/18/2023-01  **********

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
integer :: i,nk,ns,na,ix,iy,iz,ii,ie,l,jx,jy,jz,ktot,nk0,ns0,neig0,nk1
integer :: npxmax,npymax,npzmax,nsumk
real*8  :: x,y,z
real*8,allocatable::ele(:,:,:)
real*8 ene_min,ene_max,ferm,alph,de,e,weight,pi
integer nene
logical:: ldirinout
character:: chdirinp*200,chdirout*200,fname*200
!integer mpij,mpi_comm_world,myrank_glbl,nprocs,nprocw,nprock,myr_kpt,myr_space,myrz,myry,myrx
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
  fname='output_orbcharge.dat'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  open(nndisp,file=fname)
  write(6,*) 'computed by orbcharge8f  04/18/2023-01'
  write(nndisp,*) 'computed by orbcharge8f  04/18/2023-01'
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

  ktot=0
  fname='.kpt_tmp.txt'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  open(10,file=fname)
    do nk=1,numkmx
      read(10,*) x,y,z,ii
      ktot=ktot+ii
    end do
  close(10)

  call read_input_orbcharge(chdirinp,chdirout,ene_min,ene_max,ferm,neig0,ns0,nk0,ndisp)

  nsumk=0
  allocate(ele(2*nxmax,2*nymax,2*nzmax))
  allocate(atx(natom),aty(natom),atz(natom),fatx(natom),faty(natom),fatz(natom))
  allocate(numz(natom))
  ele=0.0d0
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

    if (myr_space == 0) nsumk=nsumk+numk

    fname='.kpt_tmp.txt'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(10,file=fname)
      do nk=1,nsumk-numk
      read(10,*) x,y,z,ii
      end do
      do nk=1,numk
      read(10,*) x,y,z,nwskp(nk)
      end do
    close(10)

    if ((nk0 < 1) .or. (nk0 .ge. nsumk-numk+1) .and. (nk0 .le. nsumk)) then
      call read_wf(  &
       myr_space,myr_kpt,'',nrc,nums,ncol,neigmx,numk,ncpx,ncpy,ncpz,ksym, & ! <
       key_ksym_inv2,                                                 & ! <
       svecre,sveccm,sval)                                              ! >
      write(ndisp,*) myr_space,myr_kpt
      if (neig0 < 1) then
        nk1=nk0-(nsumk-numk)
        do nk=1,numk
        do ns=1,nums
        if (((nk0 < 1) .or. (nk1 == nk)) .and. ((ns0 < 1) .or. (ns0 == ns)))then
        select case(nrc)
        case(0)
          do l=1,neigmx
            if ((sval(l,ns0,nk)-ferm .lt. ene_max) .and. (sval(l,ns0,nk)-ferm .gt. ene_min)) then
              do iz=1,ncpz
              do iy=1,ncpy
              do ix=1,ncpx
                jz=iz+ncpz*myrz
                jy=iy+ncpy*myry
                jx=ix+ncpx*myrx
                ele(jx,jy,jz)=ele(jx,jy,jz)+svecre(ix,iy,iz,l,ns0,nk)*svecre(ix,iy,iz,l,ns0,nk)*dble(nwskp(nk))/dble(ktot)
              end do
              end do
              end do
            end if
          end do
        case(1)
          do l=1,neigmx
            if ((sval(l,ns0,nk)-ferm .lt. ene_max) .and. (sval(l,ns0,nk)-ferm .gt. ene_min)) then
              do iz=1,ncpz
              do iy=1,ncpy
              do ix=1,ncpx
                jz=iz+ncpz*myrz
                jy=iy+ncpy*myry
                jx=ix+ncpx*myrx
                ele(jx,jy,jz)=ele(jx,jy,jz)+real(dconjg(sveccm(ix,iy,iz,l,ns0,nk))*sveccm(ix,iy,iz,l,ns0,nk))*dble(nwskp(nk))/dble(ktot)
              end do
              end do
              end do
            end if
          end do
        end select
        end if
        end do
        end do
      else
        nk=nk0-(nsumk-numk)
        ns=ns0
        select case(nrc)
        case(0)
          do iz=1,ncpz
          do iy=1,ncpy
          do ix=1,ncpx
            jz=iz+ncpz*myrz
            jy=iy+ncpy*myry
            jx=ix+ncpx*myrx
            ele(jx,jy,jz)=ele(jx,jy,jz)+svecre(ix,iy,iz,neig0,ns,nk)*svecre(ix,iy,iz,neig0,ns,nk)
          end do
          end do
          end do
        case(1)
          do iz=1,ncpz
          do iy=1,ncpy
          do ix=1,ncpx
            jz=iz+ncpz*myrz
            jy=iy+ncpy*myry
            jx=ix+ncpx*myrx
            ele(jx,jy,jz)=ele(jx,jy,jz)+real(dconjg(sveccm(ix,iy,iz,neig0,ns,nk))*sveccm(ix,iy,iz,neig0,ns,nk))
          end do
          end do
          end do
        end select
      end if
    end if
    deallocate(svecre,sveccm,sval)
    deallocate(nwskp)
  end do

  fname='.kpt_tmp.txt'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  open(10,file=fname)
  close(10,status='delete')

  fname='ele.dat'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  open(10,file=fname)
    do iz=1,2*nzmax
    do iy=1,2*nymax
    do ix=1,2*nxmax
      x=(ix-0.5d0)*xmax/nxmax-xmax
      y=(iy-0.5d0)*ymax/nymax-ymax
      z=(iz-0.5d0)*zmax/nzmax-zmax
      write(10,'(5e20.10)') x,y,z,ele(ix,iy,iz)
    end do
    end do
    end do
  close(10)

  fname='atom.xyz'
  if (len_trim(chdirinp) > 0) fname=trim(chdirinp)//'/'//fname
  open(10,file=fname)
  read (10,fmt='(1x)')
  do na=1,natom
    read(10,*) atx(na),aty(na),atz(na),numz(na)
  end do
  close(10)

  fatx=0.0d0
  faty=0.0d0
  fatz=0.0d0
  fname='ele.xsf'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  open(10,file=fname)
  write(10,*) 'ATOMS'
  do na=1,natom
    write(10,'(i3,3x,6f13.8)') numz(na),atx(na)*0.52918d0,aty(na)*0.52918d0,atz(na)*0.52918d0,fatx(na),faty(na),fatz(na)
  end do
  write(10,*) 'BEGIN_BLOCK_DATAGRID_3D'
  write(10,*) ' my_first_example_of_3D_datagrid'
  write(10,*) ' BEGIN_DATAGRID_3D_this_is_3Dgrid#1'
  write(10,'(3i6)') 2*nxmax,2*nymax,2*nzmax
  write(10,'(3f8.3)') (-xmax+0.5d0*dx)*0.52918d0,(-ymax+0.5d0*dy)*0.52918d0,(-zmax+0.5d0*dz)*0.52918d0
  write(10,'(3f8.3)') (2*nxmax*dx)*0.52918d0,0.0d0,0.0d0
  write(10,'(3f8.3)') 0.0d0,(2*nymax*dy)*0.52918d0,0.0d0
  write(10,'(3f8.3)') 0.0d0,0.0d0,(2*nzmax*dz)*0.52918d0
  do iz=1,2*nzmax
  do iy=1,2*nymax
  write(10,'(1000f15.10)') (ele(ix,iy,iz),ix=1,2*nxmax)
  end do
  end do
  write(10,*) ' END_DATAGRID_3D'
  write(10,*) 'END_BLOCK_DATAGRID_3D'
  close(10)

  close(nndisp)

  call mpi_barrier(mpi_comm_world,mpij)
  call mpi_finalize(mpij)

stop
 2000 call stopp ('error found when reading kpt_tmp.txt')
 2010 call stopp ('error found when removing kpt_tmp.txt')
end
