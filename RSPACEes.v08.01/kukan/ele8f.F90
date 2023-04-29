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
! **********  ele8f.f90 04/18/2023-01  **********

use mod_mpi
use mod_stopp

use var_keys
use var_scalars
use var_arrays
use var_read_input_kukan
use var_read_ppfile
use var_pseudopotentials
use var_listvec

use dim_arrays

use mod_read_input_kukan
use mod_rw_wfpot

implicit none
integer :: i,ns,ix,iy,iz,lbrydn,n,na
real*8  :: x,y,z
real*8,allocatable::rho(:,:,:,:)
logical:: ldirinout
character:: chdirinp*200,chdirout*200,fname*200

  call mpi_init(mpij)
  call mpi_comm_rank(mpi_comm_world,myrank_glbl,mpij)
  mpicom_space=mpi_comm_world

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
  fname='output_ele.dat'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  open(nndisp,file=fname)
  write(6,*) 'computed by ele8f 04/18/2023-01'
  write(nndisp,*) 'computed by ele8f 04/18/2023-01'
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
  fname='.kpt_tmp.txt'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  open(10,file=fname)
  close(10,status='delete')
  allocate(rho(2*nxmax,2*nymax,2*nzmax,nums))
  allocate(atx(natom),aty(natom),atz(natom),fatx(natom),faty(natom),fatz(natom))
  allocate(numz(natom))
  allocate(rhosmt_i(ncpx,ncpy,ncpz,nspv))
  allocate(rhotrur_i(nradmx,npoint,nspv,num_atcell),rhosmtr_i(nradmx,npoint,nspv,num_atcell))
  allocate(natcell_inf(num_atcell))

  do i=0,nprocw-1

    numk= numkmx/nprock
    myr_kpt=i/nprocs
    if (myr_kpt >= nprock-mod(numkmx,nprock) ) numk=numk+1
    myr_space=mod(i,nprocs)
    myrz=myr_space/(nprocx*nprocy)
    myry=(myr_space-myrz*nprocx*nprocy)/nprocx
    myrx=myr_space-myrz*nprocx*nprocy-myry*nprocx

    if (myr_kpt==0) then
      call read_rho( &
       i,'',nevhist,nspv,ncpx,ncpy,ncpz,nradmx,npoint,num_atcell, & ! <
       lbrydn,n,natcell_inf,rhosmt_i,rhotrur_i,rhosmtr_i)               ! >

      do ns=1,nums
      do iz=1,ncpz
      do iy=1,ncpy
      do ix=1,ncpx
        rho(ix+myrx*ncpx,iy+myry*ncpy,iz+myrz*ncpz,ns)=rhosmt_i(ix,iy,iz,ns)
      end do
      end do
      end do
      end do
    end if

  end do

  fname='ele.dat'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  open(10,file=fname)
    do iz=1,2*nzmax
    do iy=1,2*nymax
    do ix=1,2*nxmax
      x=(ix-0.5d0)*xmax/nxmax-xmax
      y=(iy-0.5d0)*ymax/nymax-ymax
      z=(iz-0.5d0)*zmax/nzmax-zmax
      write(10,'(5e20.10)') x,y,z,(rho(ix,iy,iz,ns),ns=1,nums)
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
  if (nums==1) then
    fname='ele.xsf'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(10,file=fname)
    write(10,'(a)') "CRYSTAL"
    write(10,'(a)') "PRIMVEC"
    write(10,'(3(f15.9))') 2*xmax*0.52918d0, 0.0d0, 0.0d0
    write(10,'(3(f15.9))') 0.0d0, 2*ymax*0.52918d0, 0.0d0
    write(10,'(3(f15.9))') 0.0d0, 0.0d0, 2*zmax*0.52918d0
    write(10,'(a)') "PRIMCOORD"
    write(10,'(i9,i9)') natom, 1
    do na=1,natom
      write(10,'(i3,3x,6f13.8)') numz(na),atx(na)*0.52918d0,aty(na)*0.52918d0,atz(na)*0.52918d0,fatx(na),faty(na),fatz(na)
    end do
    write(10,*) 'BEGIN_BLOCK_DATAGRID_3D'
    write(10,*) ' 3D_PWSCF'
    write(10,*) ' BEGIN_DATAGRID_3D_UNKNOWN'
    write(10,'(3i6)') 2*nxmax,2*nymax,2*nzmax
!    write(10,'(3f8.3)') (-xmax+0.5d0*dx)*0.52918d0,(-ymax+0.5d0*dy)*0.52918d0,(-zmax+0.5d0*dz)*0.52918d0
    write(10,'(3f8.3)') -xmax*0.52918d0,-ymax*0.52918d0,-zmax*0.52918d0
    write(10,'(3f8.3)') (2*nxmax*dx)*0.52918d0,0.0d0,0.0d0
    write(10,'(3f8.3)') 0.0d0,(2*nymax*dy)*0.52918d0,0.0d0
    write(10,'(3f8.3)') 0.0d0,0.0d0,(2*nzmax*dz)*0.52918d0
    do iz=1,2*nzmax
    do iy=1,2*nymax
    write(10,'(1000f8.3)') (rho(ix,iy,iz,1),ix=1,2*nxmax)
    end do
    end do
    write(10,*) ' END_DATAGRID_3D'
    write(10,*) 'END_BLOCK_DATAGRID_3D'
    close(10)

    fname='ele.cube'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(10,file=fname)
    write(10,'(a)') "GAUSSIAN CUBE FILE"
    write(10,'(a)') "Generated by RSPACE"
    write(10,'(i7,3(f12.6))') natom, -xmax, -ymax, -zmax
    write(10,'(i7,3(f12.6))') 2*nxmax, dx, 0.0d0, 0.0d0
    write(10,'(i7,3(f12.6))') 2*nymax, 0.0d0, dy, 0.0d0
    write(10,'(i7,3(f12.6))') 2*nzmax, 0.0d0, 0.0d0, dz
    do na=1, natom
        write(10,'(i7,4(f12.6))')  numz(na), 0.0d0, atx(na), aty(na), atz(na)
    end do
    i = 0
    do ix=1,2*nxmax
        do iy=1,2*nymax
            do iz=1,2*nzmax
                write(10, '(e12.4)', advance='no') rho(ix,iy,iz,1)
                i = i + 1
                if (mod(i, 6) == 0) write(10, '(a)') ""
            end do
        end do
    end do
    write(10, '(a)') ""
    close(10)

  else
    fname='eleup.xsf'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(10,file=fname)
    write(10,'(a)') "CRYSTAL"
    write(10,'(a)') "PRIMVEC"
    write(10,'(3(f15.9))') 2*xmax*0.52918d0, 0.0d0, 0.0d0
    write(10,'(3(f15.9))') 0.0d0, 2*ymax*0.52918d0, 0.0d0
    write(10,'(3(f15.9))') 0.0d0, 0.0d0, 2*zmax*0.52918d0
    write(10,'(a)') "PRIMCOORD"
    write(10,'(i9,i9)') natom, 1
    do na=1,natom
      write(10,'(i3,3x,6f13.8)') numz(na),atx(na)*0.52918d0,aty(na)*0.52918d0,atz(na)*0.52918d0,fatx(na),faty(na),fatz(na)
    end do
    write(10,*) 'BEGIN_BLOCK_DATAGRID_3D'
    write(10,*) ' 3D_PWSCF'
    write(10,*) ' BEGIN_DATAGRID_3D_UNKNOWN'
    write(10,'(3i6)') 2*nxmax,2*nymax,2*nzmax
!    write(10,'(3f8.3)') (-xmax+0.5d0*dx)*0.52918d0,(-ymax+0.5d0*dy)*0.52918d0,(-zmax+0.5d0*dz)*0.52918d0
    write(10,'(3f8.3)') -xmax*0.52918d0,-ymax*0.52918d0,-zmax*0.52918d0
    write(10,'(3f8.3)') (2*nxmax*dx)*0.52918d0,0.0d0,0.0d0
    write(10,'(3f8.3)') 0.0d0,(2*nymax*dy)*0.52918d0,0.0d0
    write(10,'(3f8.3)') 0.0d0,0.0d0,(2*nzmax*dz)*0.52918d0
    do iz=1,2*nzmax
    do iy=1,2*nymax
    write(10,'(1000f8.3)') (rho(ix,iy,iz,1),ix=1,2*nxmax)
    end do
    end do
    write(10,*) ' END_DATAGRID_3D'
    write(10,*) 'END_BLOCK_DATAGRID_3D'
    close(10)
    fname='eledn.xsf'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(10,file=fname)
    write(10,'(a)') "CRYSTAL"
    write(10,'(a)') "PRIMVEC"
    write(10,'(3(f15.9))') 2*xmax*0.52918d0, 0.0d0, 0.0d0
    write(10,'(3(f15.9))') 0.0d0, 2*ymax*0.52918d0, 0.0d0
    write(10,'(3(f15.9))') 0.0d0, 0.0d0, 2*zmax*0.52918d0
    write(10,'(a)') "PRIMCOORD"
    write(10,'(i9,i9)') natom, 1
    do na=1,natom
      write(10,'(i3,3x,6f13.8)') numz(na),atx(na)*0.52918d0,aty(na)*0.52918d0,atz(na)*0.52918d0,fatx(na),faty(na),fatz(na)
    end do
    write(10,*) 'BEGIN_BLOCK_DATAGRID_3D'
    write(10,*) ' 3D_PWSCF'
    write(10,*) ' BEGIN_DATAGRID_3D_UNKNOWN'
    write(10,'(3i6)') 2*nxmax,2*nymax,2*nzmax
!    write(10,'(3f8.3)') (-xmax+0.5d0*dx)*0.52918d0,(-ymax+0.5d0*dy)*0.52918d0,(-zmax+0.5d0*dz)*0.52918d0
    write(10,'(3f8.3)') -xmax*0.52918d0,-ymax*0.52918d0,-zmax*0.52918d0
    write(10,'(3f8.3)') (2*nxmax*dx)*0.52918d0,0.0d0,0.0d0
    write(10,'(3f8.3)') 0.0d0,(2*nymax*dy)*0.52918d0,0.0d0
    write(10,'(3f8.3)') 0.0d0,0.0d0,(2*nzmax*dz)*0.52918d0
    do iz=1,2*nzmax
    do iy=1,2*nymax
    write(10,'(1000f8.3)') (rho(ix,iy,iz,nums),ix=1,2*nxmax)
    end do
    end do
    write(10,*) ' END_DATAGRID_3D'
    write(10,*) 'END_BLOCK_DATAGRID_3D'
    close(10)
    fname='elediff.xsf'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(10,file=fname)
    write(10,'(a)') "CRYSTAL"
    write(10,'(a)') "PRIMVEC"
    write(10,'(3(f15.9))') 2*xmax*0.52918d0, 0.0d0, 0.0d0
    write(10,'(3(f15.9))') 0.0d0, 2*ymax*0.52918d0, 0.0d0
    write(10,'(3(f15.9))') 0.0d0, 0.0d0, 2*zmax*0.52918d0
    write(10,'(a)') "PRIMCOORD"
    write(10,'(i9,i9)') natom, 1
    do na=1,natom
      write(10,'(i3,3x,6f13.8)') numz(na),atx(na)*0.52918d0,aty(na)*0.52918d0,atz(na)*0.52918d0,fatx(na),faty(na),fatz(na)
    end do
    write(10,*) 'BEGIN_BLOCK_DATAGRID_3D'
    write(10,*) ' 3D_PWSCF'
    write(10,*) ' BEGIN_DATAGRID_3D_UNKNOWN'
    write(10,'(3i6)') 2*nxmax,2*nymax,2*nzmax
!    write(10,'(3f8.3)') (-xmax+0.5d0*dx)*0.52918d0,(-ymax+0.5d0*dy)*0.52918d0,(-zmax+0.5d0*dz)*0.52918d0
    write(10,'(3f8.3)') -xmax*0.52918d0,-ymax*0.52918d0,-zmax*0.52918d0
    write(10,'(3f8.3)') (2*nxmax*dx)*0.52918d0,0.0d0,0.0d0
    write(10,'(3f8.3)') 0.0d0,(2*nymax*dy)*0.52918d0,0.0d0
    write(10,'(3f8.3)') 0.0d0,0.0d0,(2*nzmax*dz)*0.52918d0
    do iz=1,2*nzmax
    do iy=1,2*nymax
    write(10,'(1000f11.6)') (rho(ix,iy,iz,1)-rho(ix,iy,iz,nums),ix=1,2*nxmax)
    end do
    end do
    write(10,*) ' END_DATAGRID_3D'
    write(10,*) 'END_BLOCK_DATAGRID_3D'
    close(10)


    fname='eleup.cube'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(10,file=fname)
    write(10,'(a)') "GAUSSIAN CUBE FILE"
    write(10,'(a)') "Generated by RSPACE"
    write(10,'(i7,3(f12.6))') natom, -xmax, -ymax, -zmax
    write(10,'(i7,3(f12.6))') 2*nxmax, dx, 0.0d0, 0.0d0
    write(10,'(i7,3(f12.6))') 2*nymax, 0.0d0, dy, 0.0d0
    write(10,'(i7,3(f12.6))') 2*nzmax, 0.0d0, 0.0d0, dz
    do na=1, natom
        write(10,'(i7,4(f12.6))')  numz(na), 0.0d0, atx(na), aty(na), atz(na)
    end do
    i = 0
    do ix=1,2*nxmax
        do iy=1,2*nymax
            do iz=1,2*nzmax
                write(10, '(e12.4)', advance='no') rho(ix,iy,iz,1)
                i = i + 1
                if (mod(i, 6) == 0) write(10, '(a)') ""
            end do
        end do
    end do
    write(10, '(a)') ""
    close(10)

    fname='eledn.cube'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(10,file=fname)
    write(10,'(a)') "GAUSSIAN CUBE FILE"
    write(10,'(a)') "Generated by RSPACE"
    write(10,'(i7,3(f12.6))') natom, -xmax, -ymax, -zmax
    write(10,'(i7,3(f12.6))') 2*nxmax, dx, 0.0d0, 0.0d0
    write(10,'(i7,3(f12.6))') 2*nymax, 0.0d0, dy, 0.0d0
    write(10,'(i7,3(f12.6))') 2*nzmax, 0.0d0, 0.0d0, dz
    do na=1, natom
        write(10,'(i7,4(f12.6))')  numz(na), 0.0d0, atx(na), aty(na), atz(na)
    end do
    i = 0
    do ix=1,2*nxmax
        do iy=1,2*nymax
            do iz=1,2*nzmax
                write(10, '(e12.4)', advance='no') rho(ix,iy,iz,nums)
                i = i + 1
                if (mod(i, 6) == 0) write(10, '(a)') ""
            end do
        end do
    end do
    write(10, '(a)') ""
    close(10)

    fname='elediff.cube'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(10,file=fname)
    write(10,'(a)') "GAUSSIAN CUBE FILE"
    write(10,'(a)') "Generated by RSPACE"
    write(10,'(i7,3(f12.6))') natom, -xmax, -ymax, -zmax
    write(10,'(i7,3(f12.6))') 2*nxmax, dx, 0.0d0, 0.0d0
    write(10,'(i7,3(f12.6))') 2*nymax, 0.0d0, dy, 0.0d0
    write(10,'(i7,3(f12.6))') 2*nzmax, 0.0d0, 0.0d0, dz
    do na=1, natom
        write(10,'(i7,4(f12.6))')  numz(na), 0.0d0, atx(na), aty(na), atz(na)
    end do
    i = 0
    do ix=1,2*nxmax
        do iy=1,2*nymax
            do iz=1,2*nzmax
                write(10, '(e12.4)', advance='no') rho(ix,iy,iz,1)-rho(ix,iy,iz,nums)
                i = i + 1
                if (mod(i, 6) == 0) write(10, '(a)') ""
            end do
        end do
    end do
    write(10, '(a)') ""
    close(10)
    



  end if

  deallocate(rho)
  deallocate(atx,aty,atz,fatx,faty,fatz,numz)
  deallocate(rhosmt_i)
  deallocate(rhotrur_i,rhosmtr_i)

  close(nndisp)

  call mpi_barrier(mpi_comm_world,mpij)
  call mpi_finalize(mpij)

stop
 2000 call stopp ('error found when reading kpt_tmp.txt')
 2010 call stopp ('error found when removing kpt_tmp.txt')
end
