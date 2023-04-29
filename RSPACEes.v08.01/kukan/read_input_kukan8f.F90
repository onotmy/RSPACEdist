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
! **********  read_input_kukan8f.F90 04/18/2023-01  **********

module mod_read_input_kukan
implicit none
contains


subroutine read_input_kpt( &
 chdir,numkmx,             & ! <
 nwskp,                    & ! >
 skpxyz)                     ! >
use mod_mpi
use mod_stopp
implicit none
character,intent(in)::chdir*200
integer,intent(in)::numkmx
integer,intent(out)::nwskp(numkmx)
real*8 ,intent(out)::skpxyz(numkmx,3)
integer i
character::fname*200
  if (myrank_glbl==0) then
    fname='.kpt_tmp.txt'
    if (len_trim(chdir) > 0) fname=trim(chdir)//'/'//fname
    open(10,file=fname,form='formatted',err=2000)
    do i=1,numkmx
      read(10,*) skpxyz(i,1),skpxyz(i,2),skpxyz(i,3),nwskp(i)
    end do
    close(10,status='delete',err=2010)
  end if
  call mpi_bcast(skpxyz,numkmx*3,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(nwskp ,numkmx  ,mpi_integer         ,0,mpicom_space,mpij)

  return
 2000 call stopp ('error found when reading .kpt_tmp.txt')
 2010 call stopp ('error found when removing .kpt_tmp.txt')
end subroutine read_input_kpt


subroutine read_input_readinp( &
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
 lshow)                                                                         ! <
use mod_mpi
use mod_stopp
use mod_kmesh
implicit none
character,intent(in)::chdirinp*200,chdirout*200
integer,intent(inout)::ndisp
integer,intent(in)::key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp,key_jel_nocalc,key_jel_calc
integer,intent(in)::key_polcon_none,key_polcon_occ,key_ksym_inv
real*8, intent(out)::xmax,ymax,zmax
integer,intent(out)::nxmax,nymax,nzmax,npxmax,npymax,npzmax
integer,intent(out)::nso
real*8, intent(out)::socang(3) 
integer,intent(out)::nsym
integer,intent(out)::neigmx
integer,intent(out)::natom
integer,intent(out)::num_atcell
integer,intent(out)::num_ppcell
integer,intent(out)::num_ppcell_d
real*8, intent(out)::gmaxps
integer,intent(out)::nperi,npopshow
integer,intent(out)::numkmx,numkx,numky,numkz,ksym,kband
character(len=7),intent(out)::cexco
integer,intent(out)::nums
integer,intent(out)::ncol,nspv
real*8 ,intent(out)::epsvh,epssd,ratio_diis,eps_scf
integer,intent(out)::ncgmin,ncgmax,ncgres
integer,intent(out)::nprecon_cg,nprecon_diis
integer,intent(out)::nsdmax,nkscg
integer,intent(out)::ndiismax,ncgscf,nretcg
integer,intent(out)::nrrz,nchange
real*8 ,intent(out)::eta,etamag(2)
integer,intent(out)::looplimit
integer,intent(out)::nbrydn
real*8 ,intent(out)::tmstep,sconst
integer,intent(out)::nmd_start,nmd_end,ngdiis
real*8 ,intent(out)::biasx,biasy,biasz
real*8 ,intent(out)::tf,tfmin,tfmax
real*8 ,intent(out)::chrgd,polconocc
real*8 ,intent(out)::endjel,chrjel
real*8 ,intent(out)::fcut
integer,intent(out)::npolcon
integer,intent(out)::npre,nevhist
integer,intent(out)::northo
real*8 ,intent(out)::eps,eps_eig_diis,alambda_diis,alambda_min,alambda_max
integer,intent(out)::nradmx,nprjmx
integer,intent(out)::lsphel,lrhomx
integer,intent(out)::nfiltyp
real*8 ,intent(out)::psctoff
integer,intent(out)::nqmx
real*8 ,intent(out)::psftrad
real*8 ,intent(out)::psctrat,psext
real*8 ,intent(out)::filpp,rctpcc
integer,intent(out)::nmesh
real*8 ,intent(out)::veta
integer,intent(out)::new_pwx,new_pwy,new_pwz
integer,intent(out)::new_rsx,new_rsy,new_rsz
integer,intent(out)::nint1dmax
integer,intent(out)::nf,nfdg,npmesh,nfh
real*8 ,intent(out)::zs_pre,pol_pre
integer,intent(out)::npoint,nrc
real*8 ,intent(out)::dx,dy,dz
real*8 ,intent(out)::ddx,ddy,ddz
integer,intent(out)::ncpx,ncpy,ncpz
integer,intent(out)::ncpx_d,ncpy_d,ncpz_d
integer,intent(out)::nwskptot
integer,intent(out)::jelcalc
real*8 ,intent(out)::strjel
integer,intent(out)::lmx,nprmx
integer,allocatable::nwkp(:)
real*8 ,allocatable::skpx(:),skpy(:),skpz(:)
integer :: kmeshgen
real*8  pi
integer i,j,ios,loop,jfile
integer nk,nrprjmx,nlmax
character catmfn*50
logical lveffout,lcalcpdos,lopen,lshow
character fname*200
  namelist /nml_inp_prm_kukan/ &
   ndisp       , catmfn      , nprocx      , nprocy      , nprocz      , nprock      , &
   xmax        , ymax        , zmax        , nxmax       , nymax       , nzmax       , &
   npxmax      , npymax      , npzmax      , nsym        , neigmx      , natom       , &
   num_atcell  , num_ppcell  , num_ppcell_d, nperi       , npopshow    , numkx       , &
   numky       , numkz       , ksym        , kband       , skpx        , skpy        , &
   skpz        , nwkp        , cexco       , nspv        , epsvh       , epssd       , &
   ratio_diis  , eps_scf     , ncgmin      , ncgmax      , ncgres      , nprecon_cg  , &
   nprecon_diis, nsdmax      , nkscg       , ndiismax    , ncgscf      , nretcg      , &
   nrrz        , nchange     , eta         , looplimit   , nbrydn      , etamag      , &
   tmstep      , nmd_start   , nmd_end     , ngdiis      , sconst      , biasx       , &
   biasy       , biasz       , tf          , tfmin       , tfmax       , chrgd       , &
   npolcon     , polconocc   , endjel      , chrjel      , fcut        , npre        , &
   nevhist     , northo      , lveffout    , lcalcpdos   , nso         , socang      , &
   eps         , eps_eig_diis, alambda_diis, alambda_min , alambda_max , nradmx      , &
   nrprjmx     , nprjmx      , lsphel      , nlmax       , lrhomx      , nfiltyp     , &
   gmaxps      , psctoff     , nqmx        , psftrad     , psctrat     , psext       , &
   filpp       , rctpcc      , veta        , new_pwx     , new_pwy     , new_pwz     , &
   new_rsx     , new_rsy     , new_rsz     , nint1dmax   , nf          , nfdg        , &
   nfh         , nmesh       , npmesh      , zs_pre      , pol_pre     , kmeshgen

  pi= 4.0d0*datan(1.0d0)
  if (myrank_glbl == 0) then
    allocate(nwkp(100000),skpx(100000),skpy(100000),skpz(100000))
    ndisp       =     66   ! outputfile (6: display, other: output.dat)
    catmfn      ='atom.xyz'! filename of atomic coordinate
    nprocx      =      1   ! # of processes (x)
    nprocy      =      1   ! # of processes (y)
    nprocz      =      1   ! # of processes (z)
    nprock      =      1   ! # of processes (k)
    xmax        =  2.5d0   ! length of supercell (x in bohr, total length is 2*xmax)
    ymax        =  2.5d0   ! length of supercell (y in bohr, total length is 2*ymax)
    zmax        =  8.0d0   ! length of supercell (z in bohr, total length is 2*zmax)
    nxmax       =      8   ! # of grid points (x, total number is 2*nxmax)
    nymax       =      8   ! # of grid points (y, total number is 2*nymax)
    nzmax       =     26   ! # of grid points (z, total number is 2*nzmax)
    npxmax      =     -1   ! # of grid points in augmented sphere in (x, total number is 2*npxmax)
    npymax      =     -1   ! (if one of them is less than zero, they are automatically determined.)
    npzmax      =     -1
    nsym        =      0   ! symmetric operation. 0: non, 1: BCC, 2: FCC, 3: DIA, 4: HCP
    neigmx      =     18   ! # of states per k point
    natom       =      2   ! # of atoms
    num_atcell  =     -1   ! # of atoms per sub-domain (if less than 0, it corresponds to natom)
    num_ppcell  =     -1   ! # of non-localparts of p.p. per sub-domain on coarse grid (if less than 0, it corresponds to natom)
    num_ppcell_d=     -1   ! # of non-localparts of p.p. per sub-domain on dense grid (if less than 0, it corresponds to natom)
    nperi       =      3   ! switchs for periodic boundary conditions (0; isolated, 3; periodic)
    npopshow    =      0   ! switchs for display of atomic population (0;noshow 1;show)
    numkx       =      1   ! # of sampling k points (x)
    numky       =      1   ! # of sampling k points (y)
    numkz       =      1   ! # of sampling k points (z)
    ksym        =      0   ! symmetry of kpoints
    kband       =      0   ! swich to output band map data. 0: not output, 1: output
    kmeshgen    =      0   ! k grid generation (0: manual, 1: auto without symmetry, 2: auto)
    skpx(1)     =  0.0d0   ! kx
    skpy(1)     =  0.0d0   ! ky
    skpz(1)     =  0.0d0   ! kz
    nwkp(1)     =      1   ! weight for the point
    cexco       =  'vwn'   ! type of exchange correlation functional vwn,pz,pbe,pw91
    nspv        =      2   ! spin (1; degenerate, 2; free_collinear, 4; free_noncollinear)
    epsvh       = 1.0d-12  ! criteria of the convergency for CG (P. eq.)
    epssd       = 1.0d-6   ! criteria of the convergency for CG (KS eq.)
    ratio_diis  =  0.3d0   ! criteria of the convergency for DIIS (KS eq.)
    eps_scf     = 1.0d-6   ! criteria of the convergency for SCF
    ncgmin      =      1   ! min. # of its. for CG (P. eq.)
    ncgmax      =    800   ! max. # of its. for CG (P. eq.)
    ncgres      =    801   ! restart for CG (P. eq.):N (1/N its.)
    nprecon_cg  =      1   ! switch for preconditioning CG (KS eq.) (0;no 1;yes)
    nprecon_diis=      1   ! switch for preconditioning  DIIS (KS eq.) (0;no 1;yes)
    nsdmax      =      4   ! max. # of its. for CG (KS eq.)
    nkscg       =      1   ! switch for CG (KS eq.) (0;no 1;yes)
    ndiismax    =      4   ! max. # of its. for DIIS (KS eq.)
    ncgscf      =    200   ! min. # of SCF its. using CG before DIIS
    nretcg      =   2000   ! retry of CG:N (1/N scf)
    nrrz        =      1   ! rayleigh-ritz or orthogonalization of WFs in the case of DIIS:N (1/N scf)
    nchange     =      1   ! re-order eigenstates:N (1/N scf)
    eta         = 0.02d0   ! mixing ratio of charge density
    looplimit   =   1000   ! max. # of its. for SCF
    nbrydn      =     20   ! # of steps for Broyden mixing
    etamag(1)   = 0.02d0   ! mixing ratio for magnetic moment
    etamag(2)   =  2.0d0   ! increment factor for mixing of magnetic moment
    tmstep      =  0.0d0   ! time step of Str. Opt. (a.u.)
    nmd_start   =      0   ! its. # of SO
    nmd_end     =      1   ! its. # of SO (nmd_end-nmd_start is the total # of it. of SO)
    ngdiis      =      1   ! # of steps for GDIIS
    sconst      =  0.0d0   ! spring constant for NEB
    biasx       =  0.0d0   ! electric field in x (a.u.)
    biasy       =  0.0d0   ! electric field in y (a.u.)
    biasz       =  0.0d0   ! electric field in z (a.u.)
    tf          = 3.0d-3   ! Temp. for Fermi dist. (KT) (a.u.)
    tfmin       = -3.0d2   ! min. of expected Fermi level (a.u.)
    tfmax       =  3.0d2   ! max. of expected Fermi level (a.u.)
    chrgd       =  0.0d0   ! (tot.neg.charge)-(tot.pos.charge)
    npolcon     =      0   ! spin pol.(0; free, 1; cf. atom.mag, 2; constr.tot.)
    polconocc   =  2.7d0   ! tot. spin pol.
    endjel      =  0.0d0   ! edge of jellium (if you do not use jellium, set chrjel to be zero.)
    chrjel      =  0.0d0   ! charge of jellium (if you do not use jellium, set chrjel to be zero.)
    fcut        = 1.0d-5   ! cutoff of force (a.u.)
    npre        =      1   ! swich for computation of initial wavefunctions & charge density (0; read from files, 1; generate, 2; generate only)
    nevhist     =      3   ! swich for charge potential & density (0: wf only, 1: wf+rho, 2: wf+Htr pot., 3: wf+rho+Htr pot., -1: non SCF calc.)
    northo      =      1   ! # of its. for orthogonalization
    lveffout    = .false.  ! swich for outout KS effective potential Potential.txt
    lcalcpdos   = .false.  ! swich for output pdos parameters
    nso         =      0   ! spin-orbit parameters
    socang(1)   =  0.0d0   ! spin-orbit parameters
    socang(2)   =  0.0d0   ! spin-orbit parameters
    socang(3)   =  0.0d0   ! spin-orbit parameters
    eps         =1.0d-16   ! computer epsilon
    eps_eig_diis=1.0d-14   ! parameters for DIIS
    alambda_diis=  0.5d0   ! parameters for DIIS
    alambda_min =  0.1d0   ! parameters for DIIS
    alambda_max =  1.0d0   ! parameters for DIIS

! parameters for pseudopotential
    nradmx      =   1502   ! max. # of radial grids in pseusopotential
    nrprjmx     =      6   ! max. # of radial functions in pseudopotential. e.g., s*2 p*2 d*2 =6
    nprjmx      =     18   ! max. # of projectors in pseudopotential. e.g., s*2 p*6 d*10 =18

! parameters for numerical integration in augmented sphere
    lsphel      =      8   ! order of spherical integration
    nlmax       =      4   ! order of spherical harminics
    lrhomx      =     25   ! the number of spherical harminics

! parameters for filtering of pseudopotential
    nfiltyp     =      1   ! filtering type of pseudopotential (0; King-Smith type [KS], 1; Fermi distribution [FD])
    gmaxps      = 20.0d0   ! cutoff of pseusopotential (in sqrt(Ry))
    psctoff     = 1.05d0   ! cutoff radius radius of pp
    nqmx        =    400   ! the number of waves to expand projectors [KS]
    psftrad     = 10.0d0   ! period of the waves [KS]
    psctrat     =  1.0d0   ! cutoff ratio of the waves to expand projectors [KS]
    psext       =  2.0d0   ! cutoff ratio of the waves to vanish the projectors outside the augmented sphere [KS]
    filpp       = 0.02d0   ! filtering parameter in the case of nfiltyp=1 [FD]
    rctpcc      =  2.0d0   ! cutoff radius ratio of pcc charge in the case of nfiltyp=1 [FD]

! parameters for Ewald summation
    veta        =  0.2d0   ! separating parameter
    new_pwx     =     11   ! parameters for Ewald summation (x, plane wave part)
    new_pwy     =     11   ! parameters for Ewald summation (y, plane wave part)
    new_pwz     =     11   ! parameters for Ewald summation (z, plane wave part)
    new_rsx     =     11   ! parameters for Ewald summation (x, real part)
    new_rsy     =     11   ! parameters for Ewald summation (y, real part)
    new_rsz     =     11   ! parameters for Ewald summation (z, real part)
    nint1dmax   =   1000   ! # of grid points for numerical integration for 1D periodic boundary condition

! parameters for finite-difference & double grid
    nf          =      4   ! order of finite difference for KS Eq.
    nfdg        =      4   ! order of Lagrange interpolation of DG tech.
    nfh         =      4   ! order of finite difference for P. Eq.
    nmesh       =      2   ! the number of grids for double-grid technique (potential)
    npmesh      =      2   ! the number of grids for double-grid technique (nonlocal pseudopotential)

! parameters for initial guess for wave func. and charge dens.
    zs_pre      =  0.5d0   ! initial spread of wave functions
    pol_pre     =  1.5d0   ! initial spin polarization of charge density
    fname='parameters.inp'
    if (len_trim(chdirinp) > 0) fname=trim(chdirinp)//'/'//fname
    open(10,file=fname,form="formatted",status="old" )
    read(10, nml=nml_inp_prm_kukan, iostat=ios)
      if ( ios/=0 ) then
        write(ndisp,*) 'error in nml_inp_prm_kukan!! iostat=',ios
      end if
    close(10)
  end if

  call mpi_bcast(ios,1,mpi_integer,0,mpi_comm_world,mpij)
  if ( ios/=0 ) then
    call mpi_barrier(mpi_comm_world,mpij)
    call mpi_finalize(mpij)
    stop
  end if

  if (myrank_glbl == 0) then
    if (num_atcell < 0) num_atcell=natom
    if (num_ppcell < 0) num_ppcell=natom
    if (num_ppcell_d < 0) num_ppcell_d=natom
    if (kmeshgen > 0) then
      call generate_kmesh(numkx, numky, numkz, &
        kmeshgen, size(skpx), skpx, skpy, skpz, nwkp, numkmx)
    else
      numkmx=numkx*numky*numkz
      if (numkmx<1) goto 2130
    end if
    cexco= adjustl(cexco)
    do j= 1,len(cexco)
      if ((cexco(j:j)>='A').and.(cexco(j:j)<='Z')) &
       cexco(j:j)=char(ichar(cexco(j:j))+ichar('a')-ichar('A'))
    end do
    if (nperi<3) new_pwz=0
    if (nperi<2) new_pwy=0
    if (nperi<1) new_pwx=0
    if (nperi<3) new_rsz=0
    if (nperi<2) new_rsy=0
    if (nperi<1) new_rsx=0

    if ((nsym .ge. 1) .and. (nsym .le. 3)) then
      ymax=xmax
      zmax=xmax
      nymax=nxmax
      nzmax=nxmax
    end if

    select case (nspv)
    case (1)
      nums= 1
      ncol= 1
    case (2)
      nums= 2
      ncol= 1
    case (4)
      nums= 2
      ncol= 2
    case default
      nspv= -1
    end select 

    if (nspv==1) npolcon= key_polcon_none
    if ((npolcon==key_polcon_occ).and.((ncol==2).or.(nso==1))) npolcon= key_polcon_none

    do nk=1,numkmx
      skpx(nk)=pi/xmax*skpx(nk)
      skpy(nk)=pi/ymax*skpy(nk)
      skpz(nk)=pi/zmax*skpz(nk)
    end do

    lmx= nlmax/2 + 1
    nprmx= nrprjmx/max(1,lmx) ! sign of lmx will be checked below

    ! real or complex KS equation ?
    if ((nperi==0).and.(ncol==1)) then
      nrc= 0
    else
      nrc= 1
    end if

    ! positiveness of n?max,nproc? will be checked below
    dx=xmax/dble(max(nxmax,1))
    dy=ymax/dble(max(nymax,1))
    dz=zmax/dble(max(nzmax,1))
    ncpx=(2*nxmax)/max(nprocx,1)
    ncpy=(2*nymax)/max(nprocy,1)
    ncpz=(2*nzmax)/max(nprocz,1)

    ddx=dx/nmesh
    ddy=dy/nmesh
    ddz=dz/nmesh
    ncpx_d=ncpx*nmesh
    ncpy_d=ncpy*nmesh
    ncpz_d=ncpz*nmesh
    npoint=((lsphel+1)/2)*lsphel*4

    nwskptot=0
    do nk=1,numkmx
      nwskptot=nwskptot+nwkp(nk)
    end do

    strjel=-zmax-(endjel+zmax)
    if (abs(chrjel) .gt. 1.0d-15) then
      jelcalc= key_jel_calc
    else
      jelcalc= key_jel_nocalc
    end if

    if (nspv==-1) then
      write(ndisp,*) 'error!!'
      call stopp('nspv must be in {1,2,4}.')
    end if
    if ((nprocx<1).or.(nprocy<1).or.(nprocz<1)) then
      write(ndisp,*) 'error!!'
      call stopp('number of subdomains must be positive in each direction.')
    end if
    if ((nperi/=0).and.(2*min(nxmax,nymax,nzmax)<nf)) then
      write(ndisp,*) 'error!!'
      call stopp('not enough gridpoints for the order of finite difference for KS eq. (problem with per. bound. cond.)')
    end if
    if ((nperi/=0).and.(2*nmesh*min(nxmax,nymax,nzmax)<nfh)) then
      write(ndisp,*) 'error!!'
      call stopp('not enough gridpoints for the order of finite difference for Poisson eq. (problem with per. bound. cond.)')
    end if
    if ((max(nprocx,nprocy,nprocz)>1).and.(min(ncpx,ncpy,ncpz)<nf)) then
      write(ndisp,*) 'error!!'
      call stopp('too many real-space decompositions for the number of gridpoints (overlap problem at KS eq.) ')
    end if
    if ((max(nprocx,nprocy,nprocz)>1).and.(nmesh*min(ncpx,ncpy,ncpz)<nfh)) then
      write(ndisp,*) 'error!!'
      call stopp('too many real-space decompositions for the number of gridpoints (overlap problem at Poisson eq.) ')
    end if
    if ((nso<0).or.(nso>1)) then
      write(ndisp,*) 'error!!'
      call stopp('spin-orbit switch nso must be in {0,1}.')
    end if
    if (((nsym .eq. 1) .and. (natom .ne. 2)) .and. ((nsym .eq. 1) .and. (natom .ne. 16))) then
      write(ndisp,*) 'error!!'
      call stopp('# of atom (natom) should be 2 or 16 when you use BCC type symmetric operations.')
    end if
    if ((nsym .eq. 2) .and. (natom .ne. 4)) then
      write(ndisp,*) 'error!!'
      call stopp('# of atom (natom) should be 4 when you use FCC type symmetric operations.')
    end if
    if ((nsym .eq. 3) .and. (natom .ne. 8)) then
      write(ndisp,*) 'error!!'
      call stopp('# of atom (natom) should be 8 when you use Diamond type symmetric operations.')
    end if
    if ((nsym .eq. 4) .and. (natom .ne. 8)) then
      write(ndisp,*) 'error!!'
      call stopp('# of atom (natom) should be 8 when you use HCP type symmetric operations.')
    end if
    if ((nsym .eq. 1) .and. (natom .ne. 16) .and. (mod(nxmax,2) .ne. 0)) then
      write(ndisp,*) 'error!!'
      call stopp('nxmax should be even when you use large BCC type symmetric operations.')
    end if
    if ((nsym .eq. 3) .and. (mod(nxmax,2) .ne. 0)) then
      write(ndisp,*) 'error!!'
      call stopp('nxmax should be even when you use Diamond type symmetric operations.')
    end if
    if ((nsym .gt. 0) .and. (nprocx*nprocy*nprocz .gt. 1)) then
      write(ndisp,*) 'error!!'
      call stopp('# of Proc. should be 1 when you use symmetric operations.')
    end if
    if ((ksym==key_ksym_inv).and.(ncol==2)) then
      write(ndisp,*) 'error!!'
      call stopp('noncollinear wavefunctions break inversion symmetry in k-space') 
    end if
    if (mod(lsphel,2) .eq. 1) then
      write(ndisp,*) 'error!!'
      call stopp('lsphel should be even number')
    end if
    if (mod(nlmax,2)/=0) then
      write(ndisp,*) 'error!!'
      call stopp('nlmax should be even')
    end if
    if ((nlmax+1)**2/=lrhomx) then
      write(ndisp,*) 'error!!'
      call stopp('nlmax and lrhomx are not consistent')
    end if
    if ((lmx>3).or.(lmx<1)) then
      write(ndisp,*) 'error!!'
      call stopp('only s,p,d-electrons implemented')
    end if
    if ((lmx==1).and.(nso==1)) then
      write(ndisp,*) 'error!!'
      call stopp('no spin-orbit coupling is s-electron systems')
    end if
    if (mod(nrprjmx,lmx)/=0) then
      write(ndisp,*) 'error!!'
      call stopp('nrprjmx should be multiple of lmx')
    end if
    if (nprmx>2) then
      write(ndisp,*) 'error!!'
      call stopp('implementation allows only 2 radial projetors per l')
    end if
    if ((jelcalc==key_jel_calc).and.(nperi/=3)) then
      write(ndisp,*) 'error!!'
      call stopp('jellium potential cannot be used for nperi/=3')
    end if
    if ((nmesh<1).or.(npmesh<1)) then
      write(ndisp,*) 'error!!'
      call stopp('double grid parameters nmesh, npmesh must be positive')
    endif
    if (nbrydn<0) then
      write(ndisp,*) 'error!!'
      call stopp('nbrydn cannot be nagative')
    endif
    if (eta<0.0d0) then
      write(ndisp,*) 'error!!'
      call stopp('charge mixing parameter eta must be positive')
    endif 
    if (etamag(1)<0.0d0) then
      write(ndisp,*) 'error!!'
      call stopp('magnetization mixing parameter etamag(1) must be positive')
    endif 
    if (nchange<0) then
      write(ndisp,*) 'error!!'
      call stopp('nchange must be positive (or 0, if reorderstates is not used)')
    endif
   
    do loop=1,2
      if (loop==1) jfile= ndisp
      if (loop==2) jfile= 12
      inquire(unit=jfile,opened=lopen)
      if (lopen) then
        if (nsym .eq. key_sym_bcc) then
          write(jfile,*) 'WARNING!!'
          write(jfile,*) 'BCC type symmetric operation is implemented.'
        end if
        if (nsym .eq. key_sym_fcc) then
          write(jfile,*) 'WARNING!!'
          write(jfile,*) 'FCC type symmetric operation is implemented.'
        end if
        if (nsym .eq. key_sym_dia) then
          write(jfile,*) 'WARNING!!'
          write(jfile,*) 'Diamond type symmetric operation is implemented.'
        end if
        if (nsym .eq. key_sym_hcp) then
          write(jfile,*) 'WARNING!!'
          write(jfile,*) 'HCP type symmetric operation is implemented.'
        end if
        if (lshow) then
          write(jfile,*) '==========  data from nml_inp_prm_kukan  =========='
          write(jfile,*) 'ndisp       :',ndisp
          write(jfile,*) 'catmfn      :',catmfn
          write(jfile,*) 'nprocx      :',nprocx
          write(jfile,*) 'nprocy      :',nprocy
          write(jfile,*) 'nprocz      :',nprocz
          write(jfile,*) 'nprock      :',nprock
          write(jfile,*) 'xmax        :',xmax
          write(jfile,*) 'ymax        :',ymax
          write(jfile,*) 'zmax        :',zmax
          write(jfile,*) 'nxmax       :',nxmax
          write(jfile,*) 'nymax       :',nymax
          write(jfile,*) 'nzmax       :',nzmax
          write(jfile,*) 'npxmax      :',npxmax
          write(jfile,*) 'npymax      :',npymax
          write(jfile,*) 'npzmax      :',npzmax
          write(jfile,*) 'nso         :',nso 
          write(jfile,fmt='(1x,a,2x,3(1x,e10.3))') 'socang      :',socang(:)
          write(jfile,*) 'nsym        :',nsym
          write(jfile,*) 'neigmx      :',neigmx
          write(jfile,*) 'natom       :',natom
          write(jfile,*) 'num_atcell  :',num_atcell
          write(jfile,*) 'num_ppcell  :',num_ppcell
          write(jfile,*) 'num_ppcell_d:',num_ppcell_d
          write(jfile,*) 'gmaxps      :',gmaxps
          write(jfile,*) 'nperi       :',nperi
          write(jfile,*) 'npopshow    :',npopshow
          write(jfile,*) 'numkx       :',numkx
          write(jfile,*) 'numky       :',numky
          write(jfile,*) 'numkz       :',numkz
          write(jfile,*) 'ksym        :',ksym
          write(jfile,*) 'kband       :',kband
          write(jfile,*) 'kmeshgen    :',kmeshgen
          write(jfile,*) 'cexco       :',adjustr(cexco)
          write(jfile,*) 'nspv        :',nspv
          write(jfile,*) 'epsvh       :',epsvh
          write(jfile,*) 'epssd       :',epssd
          write(jfile,*) 'ratio_diis  :',ratio_diis
          write(jfile,*) 'eps_scf     :',eps_scf
          write(jfile,*) 'ncgmin      :',ncgmin
          write(jfile,*) 'ncgmax      :',ncgmax
          write(jfile,*) 'ncgres      :',ncgres
          write(jfile,*) 'nprecon_cg  :',nprecon_cg
          write(jfile,*) 'nprecon_diis:',nprecon_diis
          write(jfile,*) 'nsdmax      :',nsdmax
          write(jfile,*) 'nkscg       :',nkscg
          write(jfile,*) 'ndiismax    :',ndiismax
          write(jfile,*) 'ncgscf      :',ncgscf
          write(jfile,*) 'nretcg      :',nretcg
          write(jfile,*) 'nrrz        :',nrrz
          write(jfile,*) 'nchange     :',nchange
          write(jfile,*) 'eta         :',eta
          write(jfile,fmt='(1x,a,4(1x,e11.4))') 'etamag      :',(etamag(j),j=1,2)
          write(jfile,*) 'looplimit   :',looplimit
          write(jfile,*) 'nbrydn      :',nbrydn
          write(jfile,*) 'tmstep      :',tmstep
          write(jfile,*) 'nmd_start   :',nmd_start
          write(jfile,*) 'nmd_end     :',nmd_end
          write(jfile,*) 'ngdiis      :',ngdiis
          write(jfile,*) 'sconst      :',sconst
          write(jfile,*) 'biasx       :',biasx
          write(jfile,*) 'biasy       :',biasy
          write(jfile,*) 'biasz       :',biasz
          write(jfile,*) 'tf          :',tf
          write(jfile,*) 'tfmin       :',tfmin
          write(jfile,*) 'tfmax       :',tfmax
          write(jfile,*) 'chrgd       :',chrgd
          write(jfile,*) 'npolcon     :',npolcon
          write(jfile,*) 'polconocc   :',polconocc
          write(jfile,*) 'endjel      :',endjel
          write(jfile,*) 'chrjel      :',chrjel
          write(jfile,*) 'fcut        :',fcut
          write(jfile,*) 'npre        :',npre
          write(jfile,*) 'nevhist     :',nevhist
          write(jfile,*) 'northo      :',northo
          write(jfile,*) 'lveffout    :',lveffout
          write(jfile,*) 'lcalcpdos   :',lcalcpdos
          write(jfile,*) 'eps         :',eps
          write(jfile,*) 'eps_eig_diis:',eps_eig_diis
          write(jfile,*) 'alambda_diis:',alambda_diis
          write(jfile,*) 'alambda_min :',alambda_min
          write(jfile,*) 'alambda_max :',alambda_max
          write(jfile,*) 'nradmx      :',nradmx
          write(jfile,*) 'nrprjmx     :',nrprjmx
          write(jfile,*) 'nprjmx      :',nprjmx
          write(jfile,*) 'lsphel      :',lsphel
          write(jfile,*) 'nlmax       :',nlmax
          write(jfile,*) 'lrhomx      :',lrhomx
          write(jfile,*) 'nfiltyp     :',nfiltyp
          write(jfile,*) 'psctoff     :',psctoff
          write(jfile,*) 'nqmx        :',nqmx
          write(jfile,*) 'psftrad     :',psftrad
          write(jfile,*) 'psctrat     :',psctrat
          write(jfile,*) 'psext       :',psext
          write(jfile,*) 'filpp       :',filpp
          write(jfile,*) 'rctpcc      :',rctpcc
          write(jfile,*) 'nmesh       :',nmesh
          write(jfile,*) 'npmesh      :',npmesh
          write(jfile,*) 'veta        :',veta
          write(jfile,*) 'new_pwx     :',new_pwx
          write(jfile,*) 'new_pwy     :',new_pwy
          write(jfile,*) 'new_pwz     :',new_pwz
          write(jfile,*) 'new_rsx     :',new_rsx
          write(jfile,*) 'new_rsy     :',new_rsy
          write(jfile,*) 'new_rsz     :',new_rsz
          write(jfile,*) 'nint1dmax   :',nint1dmax
          write(jfile,*) 'nf          :',nf
          write(jfile,*) 'nfdg        :',nfdg
          write(jfile,*) 'nfh         :',nfh
          write(jfile,*) 'zs_pre      :',zs_pre
          write(jfile,*) 'pol_pre     :',pol_pre
          write(jfile,*) '==================================================='
        end if !lshow
      end if !lopen
    end do ! loop
   
    fname='.kpt_tmp.txt'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(10,file=fname,form='formatted')
    do i=1,numkmx
      write(10,'(3e25.15,i10)') skpx(i),skpy(i),skpz(i),nwkp(i)
    end do
    close(10)
    deallocate(nwkp,skpx,skpy,skpz)

  end if !myrank_glbl==0

  call mpi_bcast(ndisp,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nprocx,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nprocy,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nprocz,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nprock,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nxmax,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nymax,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nzmax,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(npxmax,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(npymax,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(npzmax,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nso,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nsym,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(neigmx,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(natom,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(num_atcell,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(num_ppcell,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(num_ppcell_d,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nperi,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(npopshow,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(numkmx,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(numkx,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(numky,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(numkz,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ksym,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(kband,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nums,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ncol,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nspv,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ncgmin,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ncgmax,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ncgres,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nprecon_cg,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nprecon_diis,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nsdmax,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nkscg,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ndiismax,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ncgscf,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nretcg,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nrrz,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nchange,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(looplimit,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nbrydn,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nmd_start,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nmd_end,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ngdiis,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(sconst,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(npolcon,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(npre,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nevhist,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(northo,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nradmx,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nprjmx,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(lsphel,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(lrhomx,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nfiltyp,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nqmx,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nmesh,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(new_pwx,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(new_pwy,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(new_pwz,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(new_rsx,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(new_rsy,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(new_rsz,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nint1dmax,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nf,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nfdg,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(npmesh,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nfh,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(npoint,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nrc,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ncpx,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ncpy,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ncpz,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ncpx_d,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ncpy_d,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ncpz_d,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nwskptot,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(jelcalc,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(lmx,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nprmx,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(xmax,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(ymax,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(zmax,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(socang,3,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(gmaxps,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(epsvh,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(epssd,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(ratio_diis,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(eps_scf,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(eta,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(etamag,2,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(tmstep,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(biasx,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(biasy,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(biasz,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(tf,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(tfmin,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(tfmax,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(chrgd,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(polconocc,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(endjel,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(chrjel,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(fcut,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(eps,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(eps_eig_diis,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(alambda_diis,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(alambda_min,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(alambda_max,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(psctoff,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(psftrad,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(psctrat,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(psext,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(filpp,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(rctpcc,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(veta,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(zs_pre,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(pol_pre,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(dx,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(dy,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(dz,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(ddx,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(ddy,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(ddz,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(strjel,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(catmfn,50,mpi_character,0,mpi_comm_world,mpij)
  call mpi_bcast(cexco,7,mpi_character,0,mpi_comm_world,mpij)
  call mpi_bcast(lveffout,1,mpi_logical,0,mpi_comm_world,mpij)
  call mpi_bcast(lcalcpdos,1,mpi_logical,0,mpi_comm_world,mpij)

  return
 2130 call stopp ('error found in read.inp when reading numkx,numky,numkz,ksym,kband')
end subroutine read_input_readinp


subroutine read_input_atomxyz( &
 chdir,                                                       & ! <
 ndisp,catmfn,natom,nsym,lmx,                                 & ! <
 key_sym_any,key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp, & ! <
 key_soc_nocalc,key_soc_calc,                                 & ! <
 xmax,ymax,zmax,                                              & ! <
 num_spe,                                                     & ! >
 numz,nmdx,nmdy,nmdz,natsoc,                                  & ! >
 atx,aty,atz,watom)                                             ! >
use mod_mpi
use mod_stopp
use mod_tools, only:tools_symmetrizeatom
implicit none
character, intent(in) :: chdir*200
integer,intent(in)::ndisp,natom,nsym,lmx
integer,intent(in)::key_sym_any,key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp
integer,intent(in)::key_soc_nocalc,key_soc_calc 
real*8, intent(in)::xmax,ymax,zmax
integer,intent(out)::num_spe
integer,intent(out)::numz(natom),nmdx(natom),nmdy(natom),nmdz(natom),natsoc(natom,0:lmx-1)
real*8 ,intent(out)::atx(natom),aty(natom),atz(natom),watom(natom)
character, intent(in)::catmfn*50
integer  :: j
integer  :: na,ispe
character:: chargs*200,fname*200
integer, allocatable::mspe1(:)

  if (myrank_glbl==0) then
    allocate(mspe1(natom))
    fname=catmfn
    if (len_trim(chdir) > 0) fname=trim(chdir)//'/'//fname
    open (10,file=fname)
    read (10,fmt='(1x)',err=9998)
    do na=1,natom
      read (10,fmt='(a)',err=9998) chargs
      backspace(10)
      call read_input_argnum(chargs,j)
      if (j<9) then
        read (10,*,err=9998) atx(na),aty(na),atz(na),numz(na),nmdx(na),nmdy(na),nmdz(na),watom(na)
        natsoc(na,0:lmx-1)= key_soc_nocalc
      else
        read (10,*,err=9998) atx(na),aty(na),atz(na),numz(na),nmdx(na),nmdy(na),nmdz(na),watom(na),chargs(1:max(1,lmx-1))
        do j= 1,lmx-1
          read(chargs(j:j),fmt='(i1)',err=9998) natsoc(na,j)
        enddo  
      endif
      natsoc(na,0)= key_soc_nocalc
      do j= 1,lmx-1 
        if (natsoc(na,j)==key_soc_calc) then
          natsoc(na,0)= key_soc_calc 
        else 
          if (natsoc(na,j)/=key_soc_nocalc) call stopp ('error found in atom.xyz: natsoc contains unknown key')
        endif
      enddo 
    end do
    close(10)
    write(ndisp,*) 'atom.xyz has been read.'
    if (nsym .ne. key_sym_any) &
      call tools_symmetrizeatom(natom,nsym,key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp,xmax,ymax,zmax,atx,aty,atz)
    num_spe= 0
    do na= 1,natom
      ispe= 0
      do while (ispe<num_spe)
        ispe= ispe + 1
        if (mspe1(ispe)==numz(na)) ispe= num_spe+1
      end do
      if (ispe==num_spe) then
        num_spe= num_spe + 1
        mspe1(num_spe)= numz(na)
      end if
    end do
    deallocate(mspe1)
  end if !myrank_glbl==0

  call mpi_bcast(atx,natom,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(aty,natom,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(atz,natom,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(numz,natom,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nmdx,natom,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nmdy,natom,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nmdz,natom,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(watom,natom,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(natsoc,natom*lmx,mpi_integer,0,mpi_comm_world,mpij)

return

 9998 call stopp('error found in atom.xyz')
end subroutine read_input_atomxyz


subroutine read_input_atommag( &
 chdir,                                & ! <
 ndisp,natom,key_polcon2_dir,          & ! <
 npolcon2,polconb,polconmag,polconeta)   ! >
use mod_mpi
use mod_stopp
implicit none
real*8, parameter:: eps_p=1.0d-10
character,intent(in) ::chdir*200
integer,intent(in) ::ndisp,natom
integer,intent(in) ::key_polcon2_dir
integer,intent(out)::npolcon2(0:natom)
real*8 ,intent(out)::polconb(3,natom),polconmag(3,natom),polconeta(2,0:natom)
character ::fname*200
integer:: na,ns
logical:: isfile
real*8 :: sprod,pol(3)

  if (myrank_glbl==0) then
    fname='atom.mag'
    if (len_trim(chdir) > 0) fname=trim(chdir)//'/'//fname
    inquire(file=fname,exist=isfile)
    if (.not. isfile) call stopp('atom.mag does not exist, although npolcon=key_polcon_[atoms/asa]')
    open(10,file=fname,form='formatted',action='read')
    read(10,*,err=9999)
    read(10,*,err=9999) polconeta(:,0), npolcon2(0)
    read(10,*,err=9999)
    read(10,*,err=9999)
    do na= 1,natom
      read(10,*,err=9999) polconb(:,na), npolcon2(na), polconmag(:,na), polconeta(:,na)
      if (npolcon2(na)==key_polcon2_dir) then
        sprod= 0.0d0
        do ns= 1,3
          sprod= sprod +polconmag(ns,na)**2
        enddo
        sprod= dsqrt(sprod)
        if (sprod>eps_p) then
          do ns= 1,3
            pol(ns)= polconmag(ns,na)/sprod
          enddo
        else
          pol(:)= 0.0d0
        endif
        sprod= 0.0d0
        do ns= 1,3
          sprod= sprod +pol(ns)*polconb(ns,na)
        enddo
        do ns= 1,3
          polconb(ns,na)= polconb(ns,na) -pol(ns)*sprod
        enddo
      endif
    enddo
    close(10)
    write(ndisp,*) 'atom.mag has been read.'
  end if !myrank_glbl==0
  call mpi_bcast(npolcon2,natom+1,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(polconb,3*natom,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(polconmag,3*natom,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(polconeta,2*(natom+1),mpi_double_precision,0,mpicom_space,mpij)

  return

 9999 call stopp('error found in atom.mag')
end subroutine read_input_atommag


subroutine read_input_mpisetup( &
 nxmax,nymax,nzmax,numkmx, & ! <
 numk)                       ! >
use mod_mpi
use mod_stopp
implicit none
integer, intent(in)  :: nxmax,nymax,nzmax,numkmx
integer, intent(out) :: numk
integer :: irank, rankidx

  call mpi_comm_size(mpi_comm_world,nprocworld,mpij)

  nprocs= nprocx*nprocy*nprocz
  nprocw= nprocs*nprock
  if (myrank_glbl==0) then
    if (nprocw > nprocworld) call stopp('more processors required than supplied')
    if (nprocw < nprocworld) call stopp('less processors required than supplied')
    if ( (nprocx<1).or.((nprocx>1).and.(mod(nprocx,2)/=0)) &
     .or.(nprocy<1).or.((nprocy>1).and.(mod(nprocy,2)/=0)) &
     .or.(nprocz<1).or.((nprocz>1).and.(mod(nprocz,2)/=0)) ) &
    call stopp('nprox{x,y,z} must be 1 or a positive even number')
    if (mod(2*nxmax,nprocx) .ne. 0) call stopp('error mod(2*nxmax,nprocx) ne 0')
    if (mod(2*nymax,nprocy) .ne. 0) call stopp('error mod(2*nymax,nprocy) ne 0')
    if (mod(2*nzmax,nprocz) .ne. 0) call stopp('error mod(2*nzmax,nprocz) ne 0')
    if (nprock <= 0) call stopp('nprock should be positive integer')
  end if !myrank_glbl==0

  rankidx= mod(myrank_glbl,nprocs)
  call mpi_comm_split(mpi_comm_world,rankidx,myrank_glbl,mpicom_kpt,mpij)
  call mpi_comm_rank(mpicom_kpt,myr_kpt,mpij)
  call mpi_comm_split(mpi_comm_world,myr_kpt,myrank_glbl,mpicom_space,mpij)
  call mpi_comm_rank(mpicom_space,myr_space,mpij)

  myrz=myr_space/(nprocx*nprocy)
  myry=(myr_space-myrz*nprocx*nprocy)/nprocx
  myrx=myr_space-myrz*nprocx*nprocy-myry*nprocx

  allocate( numkproc(0:nprock-1) )
  do irank= 0,nprock-1
    numkproc(irank)= numkmx/nprock
    if ( irank >= nprock-mod(numkmx,nprock) ) numkproc(irank)= numkproc(irank)+1
  end do
  numk= numkproc(myr_kpt)
  if (myr_kpt/=0) deallocate(numkproc)

  if (myr_kpt==0) then
    ndim_ke= 1
  else
    ndim_ke= 0
  end if

end subroutine read_input_mpisetup


subroutine read_input_kptsnd(numkmx,numk,skpxyz,nwskp, & ! <
                             skpxx,skpyy,skpzz,nwskk)    ! >
use mod_mpi
use mod_stopp
implicit none
integer,intent(in) ::numkmx,numk
integer,intent(in) ::nwskp(numkmx)
real*8, intent(in) ::skpxyz(numkmx,3)
integer,intent(out)::nwskk(numk)
real*8, intent(out)::skpxx(numk),skpyy(numk),skpzz(numk)
integer irank,nk,nmk

  skpxx(:)= skpxyz(1:numk,1)
  skpyy(:)= skpxyz(1:numk,2)
  skpzz(:)= skpxyz(1:numk,3)
  nwskk(:)= nwskp(  :numk  )
  nk= 1
  do irank= 1,nprock-1
    nk= nk + numkproc(irank-1)
    nmk= max(numkproc(irank),1)
    call mpi_send(skpxyz(nk,1),nmk,mpi_double_precision,irank,0,mpicom_kpt,mpij)
    call mpi_send(skpxyz(nk,2),nmk,mpi_double_precision,irank,0,mpicom_kpt,mpij)
    call mpi_send(skpxyz(nk,3),nmk,mpi_double_precision,irank,0,mpicom_kpt,mpij)
    call mpi_send(nwskp( nk  ),nmk,mpi_integer         ,irank,0,mpicom_kpt,mpij)
  end do

end subroutine read_input_kptsnd


subroutine read_input_kptrcv(numk,                   & ! <
                            skpxx,skpyy,skpzz,nwskk)   ! >
use mod_mpi
use mod_stopp
implicit none
integer,intent(in) ::numk
integer,intent(out)::nwskk(numk)
real*8, intent(out)::skpxx(numk),skpyy(numk),skpzz(numk)

      call mpi_recv(skpxx,numk,mpi_double_precision,0,0,mpicom_kpt,mpistat,mpij)
  call mpi_recv(skpyy,numk,mpi_double_precision,0,0,mpicom_kpt,mpistat,mpij)
  call mpi_recv(skpzz,numk,mpi_double_precision,0,0,mpicom_kpt,mpistat,mpij)
  call mpi_recv(nwskk,numk,mpi_integer         ,0,0,mpicom_kpt,mpistat,mpij)

  end subroutine read_input_kptrcv


  subroutine read_input_argnum( &
   chargs, & ! <
   argnum)   ! >
  implicit none
  character(*), intent(in) :: chargs
  integer,      intent(out):: argnum
  integer :: pos,pos0

  argnum= 0
  pos0= len_trim(chargs)-len_trim(adjustl(chargs))
  pos= index(chargs,'!')-1
  if (pos<0) pos= len(chargs)
  pos= len_trim(chargs(1:pos))
  do while (pos>pos0)
    argnum= argnum +1
    pos= index(trim(chargs(1:pos)),' ',.true.) -1
  end do

end subroutine read_input_argnum


subroutine read_input_dos(chdirinp,chdirout,ene_min,ene_max,ferm,alph,nene,ndisp)
use mod_mpi
implicit none
character :: chdirinp*200,chdirout*200
real*8    :: ene_min,ene_max,ferm,alph
integer   :: nene,ndisp,ios,jfile,loop
logical lopen
character :: fname*200
  namelist /nml_inp_ldos/ ene_min, ene_max, ferm, alph, nene

  ene_min=-1.0d0/27.212d0
  ene_max=-ene_min
  alph=10000.0d0
  nene=100
  ferm=0.0d0
  fname='fermilevel.dat'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  inquire(file=fname,exist=lopen)
  if (lopen) then
    open(10,file=fname,form='formatted',action='read')
    read(10,*) ferm
    close(10)
  end if
  fname='input_dos.txt'
  if (len_trim(chdirinp) > 0) fname=trim(chdirinp)//'/'//fname
  open(10,file=fname,form='formatted',status='old' )
    read(10, nml=nml_inp_ldos, iostat=ios)
    if ( ios/=0 ) then
      write(ndisp,*) 'error in nml_inp_ldos!! iostat=',ios
    end if
  close(10)
  do loop=1,2
    if (loop==1) jfile= ndisp
    if (loop==2) jfile= 12
    inquire(unit=jfile,opened=lopen)
    if (lopen) then
      write(jfile,*) '==========  data from nml_inp_dos  =========='
      write(jfile,*) 'ene_min     :',ene_min
      write(jfile,*) 'ene_max     :',ene_max
      write(jfile,*) 'nene        :',nene
      write(jfile,*) 'ferm        :',ferm
      write(jfile,*) 'alph        :',alph
      write(jfile,*) '============================================='
    end if !lopen
  end do ! loop
  call mpi_bcast(ene_min,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(ene_max,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(nene,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ferm,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(alph,1,mpi_double_precision,0,mpi_comm_world,mpij)
  return
end subroutine read_input_dos


subroutine read_input_orbcharge(chdirinp,chdirout,ene_min,ene_max,ferm,neig,ns,nk,ndisp)
use mod_mpi
implicit none
character :: chdirinp*200,chdirout*200
real*8    :: ene_min,ene_max,ferm
integer :: neig,ns,nk,ndisp,ios,loop,jfile
logical lopen,lerr
character :: fname*200
  namelist /nml_inp_orbcharge/ ene_min,ene_max,ferm,neig,ns,nk

  ene_min=0.0d0
  ene_max=0.0d0
  neig=1
  ns=1
  nk=1
  ferm=0.0d0
  fname='fermilevel.dat'
  if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
  inquire(file=fname,exist=lopen)
  if (lopen) then
    open(10,file=fname,form='formatted',action='read')
    read(10,*) ferm
    close(10)
  end if
  fname='input_orbcharge.txt'
  if (len_trim(chdirinp) > 0) fname=trim(chdirinp)//'/'//fname
  open(10,file=fname,form='formatted',status='old' )
    read(10, nml=nml_inp_orbcharge, iostat=ios)
    if ( ios/=0 ) then
      write(ndisp,*) 'error in nml_inp_ldos!! iostat=',ios
    end if
  close(10)
  lerr=.false.
  do loop=1,2
    if (loop==1) jfile= ndisp
    if (loop==2) jfile= 12
    inquire(unit=jfile,opened=lopen)
    if (lopen) then
      write(jfile,*) '==========  data from nml_inp_orbchage  =========='
      write(jfile,*) 'ene_min     :',ene_min
      write(jfile,*) 'ene_max     :',ene_max
      write(jfile,*) 'ene_max     :',ene_max
      write(jfile,*) 'ferm        :',ferm
      write(jfile,*) 'ns          :',ns
      write(jfile,*) 'nk          :',nk
      write(jfile,*) '=================================================='
      if (neig < 1) then
        write(jfile,*) 'charge is computed by ENERGY'
        if (nk < 1) then
          write(jfile,*) 'charge is integrated in BZ'
        end if
        if (ns < 1) then
          write(jfile,*) 'charge is integrated in SPIN'
        end if
      else
        write(jfile,*) 'charge is computed by INDEX'
        if (nk < 1) then
          write(jfile,*) 'ERROR! nk shoud be larger than 1 if neig > 0.'
          lerr=.true.
        end if
        if (ns < 1) then
          write(jfile,*) 'ERROR! ns shoud be larger than 1 if neig > 0.'
          lerr=.true.
        end if
      end if
    end if !lopen
  end do ! loop
  if (lerr) then
    call mpi_barrier(mpi_comm_world,mpij)
    call mpi_finalize(mpij)
    stop
  end if
  call mpi_bcast(ene_min,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(ene_max,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(neig,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(ferm,1,mpi_double_precision,0,mpi_comm_world,mpij)
  call mpi_bcast(ns,1,mpi_integer,0,mpi_comm_world,mpij)
  call mpi_bcast(nk,1,mpi_integer,0,mpi_comm_world,mpij)
  return
end subroutine read_input_orbcharge


end module mod_read_input_kukan
