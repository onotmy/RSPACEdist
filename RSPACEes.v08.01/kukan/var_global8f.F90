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
! **********  var_global8f.F90  12/05/2022-01  **********

module var_scalars
implicit none 

! *****  convergence of cg/diis   
logical:: ksconv
! *****  max. number of iterations needed by cg/diis   
integer:: ksitmax
! *****  total electronic charge 
real*8 :: tnumele
! *****  distance of in charge and magnetization densities  
real*8 :: diff_charge
! *****  grid point cutoff for pseudopotential filtering  
real*8 :: gridmax
! *****  Fermi energy 
real*8 :: ferm
! *****  components of the total energy
real*8 :: eneco,eneel,eneha,eneex,eneat,eneth,eneof,enejel,enebc,enespr
! *****  jellium term of the forces 
real*8 :: fjel

end module var_scalars 

! ========================================================================================================

module var_arrays
implicit none

! *****  sum of matrix elements \hat{Dij}+Dij^1-\tilda{Dij^1} [see PRB59 1758 (1999)]
real*8,allocatable::dij(:,:,:,:)
! *****  matrix elements of true partial waves with spin-orbit operator
real*8,allocatable::dijsoc(:,:,:,:)
! *****  force acting on atoms
real*8,allocatable::fatx(:),faty(:),fatz(:)
! *****  occupation of atomic wave functions
real*8,allocatable::atocc(:,:,:)
! *****  smooth charge: input,output (coarse)
real*8,allocatable::rhosmt_i(:,:,:,:),rhosmt_o(:,:,:,:)
! *****  true, and smooth charges: input,output (radial grid)
real*8,allocatable::rhotrur_i(:,:,:,:),rhosmtr_i(:,:,:,:),rhotrur_o(:,:,:,:),rhosmtr_o(:,:,:,:)
! *****  wave functions
real*8,    allocatable::svecre(:,:,:,:,:,:)
complex*16,allocatable::sveccm(:,:,:,:,:,:)
! *****  wave functions, used for spin-orbit coupling in second variation
complex*16,allocatable::svecsoc(:,:,:,:,:,:)
! *****  product between H and \Psi, H is Hamiltonian matrix
real*8,    allocatable::hsvre(:,:,:,:)
complex*16,allocatable::hsvcm(:,:,:,:,:)
! *****  product between S and \Psi, S is overlap matrix of PAW
real*8,    allocatable::ssvre(:,:,:,:,:,:)
complex*16,allocatable::ssvcm(:,:,:,:,:,:)
! *****  inner product between projector and \Psi
real*8,    allocatable::rspsep(:,:,:,:,:)
complex*16,allocatable::cspsep(:,:,:,:,:)
! *****  inner product between projector and \Psi, used for spin-orbit coupling in second variation
complex*16,allocatable::cspsepsoc(:,:,:,:,:)
! *****  eigen state energy, occupation number (only wfc on local cpu)
real*8,allocatable::sval(:,:,:),fnele(:,:,:)
! *****  eigen state energy, occupation number for spin-orbit coupling in second variation (only wfc on local cpu)
real*8,allocatable::svalsoc(:,:,:),fnelesoc(:,:,:)
! *****  residual norm of eigen states (only wf of local cpu)
real*8,allocatable::residual_states(:,:,:)
! *****  temporaly arrayt for product between S and \Psi, S is overlap matrix of PAW
!real*8,allocatable::asvre(:,:,:,:),asvim(:,:,:,:)
! *****  effective potential on coarse grid
real*8,allocatable::veff(:,:,:,:)

! *****  eigen state energy, occupation number (all wfc)
real*8,allocatable::sval_wfc(:,:,:),fnele_wfc(:,:,:)
! *****  residual norm of eigen states (all wfc)
real*8,allocatable::residual_states_wfc(:,:,:)
! *****  atomic monopole computed by fuzzy cell method
real*8,allocatable::atmpole(:,:)
! *****  hartree potential (coarse), exchange correlation potential (coarse)
real*8,allocatable::vh_coarse(:,:,:),vx(:,:,:,:)
! *****  smooth + hugmented charge (coarse, inversely interpolated)
real*8,allocatable::rho_coarse(:,:,:)
! *****  fuzzy cell weight (coarse)
real*8,allocatable::pwei(:,:,:,:)
! *****  hartree potential (dense), exchange correlation potential (dense), exchange correlation energy (dense)
real*8,allocatable::vh_dense(:,:,:),vx_dense(:,:,:,:),ex_dense(:,:,:)
! *****  augmented charge (dense), smooth + pcc charge (dense)
real*8,allocatable::rho_aug_dense(:,:,:),rhosmt_pcc_dense(:,:,:,:)
! *****  boundary value of Hartree potential
real*8,allocatable::vboundx(:,:,:,:,:),vboundy(:,:,:,:,:),vboundz(:,:,:,:,:)
! *****  augmented charges (radial grid)
real*8,allocatable::rhoaugr(:,:,:)
! *****  true + core and smooth + pcc charges (radial grid)
real*8,allocatable::rhotrucorer(:,:,:,:),rhosmt_pccr(:,:,:,:)
! *****  exchange correlation potentials and energy of true + core charge (radial grid)
real*8,allocatable::vxctru(:,:,:,:),exctru(:,:,:)
! *****  exchange correlation potentials and energy of smooth + pcc charge (radial grid)
real*8,allocatable::vxcsmt(:,:,:,:),excsmt(:,:,:)
! *****  Hartree potentials of true, smooth, and augmented charges (radial grid)
real*8,allocatable::vhtrur(:,:,:),vhsmtr(:,:,:),vhaugr(:,:,:)
! ! *****  moments of true, smooth, and augmented charges
! real*8,allocatable::rhotrum(:,:,:),rhosmtm(:,:,:),rhoaugm(:,:,:)
! *****  point and weight for spherical integration, and spherical harmonics
real*8,allocatable::point(:,:),wt(:),yylm(:,:),dylm_dtheta(:,:),d2ylm_dtheta2(:,:),dylm_dphi(:,:),d2ylm_dphi2(:,:),d2ylm_dtheta_dphi(:,:)
! ***** spin polarization
real*8,allocatable::spinpol(:,:)

! *****  density atom information (currently only dummy)
integer,allocatable::natcell_inf(:)

! ! *****  smooth + augmented charge (dense)
! real*8,allocatable::rho_dense(:,:,:,:)

end module var_arrays

! ========================================================================================================

module var_read_input_kukan
implicit none

! ==========  parameters read from 'read.inp'  ==========
! *****  outputfile handle 
integer ndisp,nndisp
! *****  length of supercell
real*8 xmax,ymax,zmax
! *****  the number of grid points
integer nxmax,nymax,nzmax
! *****  spin-orbit switch 
integer nso 
! *****  Euler angles to rotate spin-orbit coordinate system 
real*8  socang(3)
! *****  symmetry
integer nsym
! *****  the number of states per k point
integer neigmx
! *****  the number of atoms
integer natom
! *****  the number of atoms per sub-domain
integer num_atcell
! *****  the number of non-localparts of p.p. per sub-domain on coarse grid
integer num_ppcell
! *****  the number of non-localparts of p.p. per sub-domain on dense grid
integer num_ppcell_d
! *****  cutoff of pseusopotential
real*8 gmaxps
! *****  switchs for periodic boundary conditions and atomic population
integer nperi,npopshow
! *****  the number of sampling k-points
integer numkx,numky,numkz,numkmx
! *****  symmetry of sampling k-points
integer ksym
! *****  band structure switch  
integer kband 
! *****  weight of sampling k-point
integer,allocatable::nwskp(:)
! *****  sampling k-point
real*8,allocatable::skpxyz(:,:)
! *****  type of exchange correlation functional
character cexco*7
! *****  filename of atomic coordinate
character catmfn*50
! *****  the number of spins
integer nums
! *****  independent spin channels (2 for collinear magnetism, 1 else), spin components of potential (1,2,4)
integer ncol, nspv
! *****  criterias of the convergencies
real*8 epsvh,epssd,ratio_diis,eps_scf
! *****  min. # of its. for CG (P. eq.), max. # of its. for CG (P. eq.), restart for CG (P. eq.)
integer ncgmin,ncgmax,ncgres
! *****  switches for preconditioning CG and DIIS (KS eq.)
integer nprecon_cg,nprecon_diis
! *****  max. # of its. for CG (KS eq.), switch for CG (KS eq.)
integer nsdmax,nkscg
! *****  max. # of its. for DIIS (KS eq.), min. # of SCF its. using CG before DIIS, retry of CG
integer ndiismax,ncgscf,nretcg
! *****  rayleigh-ritz or orthogonalization of WFs in the case of DIIS, re-order eigenstates
integer nrrz,nchange
! *****  mixing ratio of charge density
real*8 eta
! *****  parameters for mixing of magnetization (alpha_mag, magnorm, beta, theta_max) 
real*8 etamag(2)
! *****  max. # of its. for SCF
integer looplimit
! *****  # of steps for Broyden mixing
integer nbrydn
! *****  time step of Str. Opt., spring const. for NEB
real*8 tmstep,sconst
! *****  its. # of SO, its. # of SO, # of points for NEB
integer nmd_start,nmd_end,ngdiis
! *****  electric field
real*8 biasx,biasy,biasz
! *****  Temp. for Fermi dist., min. and max. of expected Fermi level
real*8 tf,tfmin,tfmax
! *****  (total negative charge) - (total positive charge)
real*8 chrgd
! *****  constrained total magnetization (constrained via occupation numbers)
real*8 polconocc
! *****  the edge and charge of jellium
real*8 endjel,chrjel
! *****  cutoff of force
real*8 fcut
! *****  swich for constrained magnetization
integer npolcon
! *****  switch for computation of initial wavefunctions, initial density and potential
integer npre,nevhist
! *****  # of its. for orthogonalization
integer northo
! ----------------------------------------------------------------
! ***** computer epsilon
real*8 eps
! ***** parameters for DIIS
real*8 eps_eig_diis,alambda_diis,alambda_min,alambda_max
! ! ***** order of spherical integration, order of spherical harmonics, the number of spherical harmonics
! integer lsphel,lmx,lrhomx
! ***** filtering type of pseudopotential
integer nfiltyp
! ***** cutoff radius ratio of pp
real*8 psctoff
! ***** the number of waves to expand projectors
integer nqmx
! ***** period of the waves
real*8 psftrad
! ***** cutoff ratio of the waves to expand projectors, cutoff ratio of the waves to vanish the projectors outside the augmented sphere
real*8 psctrat,psext
! ***** filtering parameter for the Fermi distribution
real*8 filpp,rctpcc
! ***** the number of grids for double-grid technique
integer nmesh,npmesh
! ***** parameter for Ewald summation
real*8 veta
! ***** parameters for Ewald summation (plane wave part)
integer new_pwx,new_pwy,new_pwz
! ***** parameters for Ewald summation (real part)
integer new_rsx,new_rsy,new_rsz
! ***** parameter for 1D periodic boundary condition
integer nint1dmax
! ***** order of finite difference for KS Eq., order of Lagrange interpolation of DG tech., order of finite difference for P. Eq.
integer nf,nfdg,nfh
! ***** initial spread of wave functions, initial spin polarization of charge density
real*8 zs_pre,pol_pre
! ==========  parameters read from 'atom.xyz'  ==========
! ***** the number of species
integer num_spe
! *****  atomic numbers
integer,allocatable::numz(:)
! *****  switch for structural optimization
integer,allocatable::nmdx(:),nmdy(:),nmdz(:)
! *****  atom coordinate
real*8,allocatable::atx(:),aty(:),atz(:)
! *****  atom mass
real*8,allocatable::watom(:)
! *****  spin-orbit switch  
integer,allocatable::natsoc(:,:) 
! ==========  parameters read from 'atom.mag'  ==========
! ***** switch for constraining magnetic moment
integer,allocatable:: npolcon2(:)
! ***** constraining B-field
real*8, allocatable:: polconb(:,:)
! ***** magnetic moment or direction
real*8, allocatable:: polconmag(:,:)
! ***** proportionality factor of B_con and -m
real*8, allocatable:: polconeta(:,:)
! ==========  parameters calculated in read_input_readinp  ==========
! *****  the number of points for spherical integration, switch for real/complex wave functions
integer :: npoint,nrc
! *****  coarse grid spacing
real*8  :: dx,dy,dz
! *****  dense grid spacing
real*8  :: ddx,ddy,ddz
! *****  the number of coarse grid points in the subdomain per one dimension
integer :: ncpx,ncpy,ncpz
! *****  the number of dense grid points in the subdomain per one dimension
integer :: ncpx_d,ncpy_d,ncpz_d
! *****  total weight of sampling k points
integer :: nwskptot
! *****  swich for jellium calculation
integer :: jelcalc
! *****  starting edge of jellium electrode
real*8  :: strjel

! *****  swich for output KS effective potential Potential.txt
logical  :: lveffout
! *****  swich for output pdos parameters
logical  :: lcalcpdos

! ==========  cpu-specific parameters  ==========
! *****  the number of eigen value per process, the number of k point per process
integer :: numk
! *****  k points
real*8, allocatable :: skpxx(:),skpyy(:),skpzz(:)
! *****  the weight for k points
integer,allocatable :: nwskk(:)

end module var_read_input_kukan

! ========================================================================================================

module var_read_ppfile
implicit none

! fixed values
! ***** max. number of radial grid points
integer nradmx
! ***** max. number of projectors, max. number of projector func. per l
integer nprjmx,nprmx
! ***** order of spherical integration, order of spherical harmonics, the number of spherical harmonics
integer lsphel,lmx,lrhomx
! ***** type of pseudopotential
integer,allocatable::ntyppp(:)
! ***** index of the cutoff radius
integer,allocatable::nradct(:)
! ***** the number of projectors
integer,allocatable::nprj(:)
! ***** the number of valence electrons, parameters for core potential
real*8,allocatable::cp(:,:)
! ***** radial grid points, derivative of radial grid points
real*8,allocatable::radial(:,:),dradial(:,:)
! ***** core potential, all-electron wave function, pseudized wave function, core electron, pcc charge
real*8,allocatable::potc(:,:),awf(:,:,:),pwf(:,:,:),rhocore(:,:),rhopcc(:,:)
! ***** all-electron eigenvalues, all-electron potential 
real*8,allocatable::aeeig(:,:),aepot(:,:) 
! ***** parameters for Weinert type augmented charge
real*8,allocatable::grwein(:)
! ***** overlap matrix, product between projectors and kinetic operator + core potential, projectors
real*8,allocatable::sss0(:,:,:),akv0(:,:,:),prj(:,:,:)
! ***** cutoff of pseudpotential
real*8,allocatable::gmax(:)
! ***** the number of projector per one azimuthal quantum number (1 or 2), the max. azimuthal quantum numbers of projectors (starts from 0)
integer,allocatable::npr(:,:),lpmx(:)
! ***** total azimuthal quantum numbers of projectors, the number of grid points on radial grid
integer,allocatable::nrprj(:),nradps(:)
! ***** parameter for Weinert type augmented charge, should be same with pseusopotential generation code
real*8 gmaxqp

end module var_read_ppfile

! ========================================================================================================

module var_pseudopotentials
implicit none

! fixed values
! ***** atom information in the subdomain, index number of species, atomic number in spieces labelling
integer,allocatable::indspe(:),mspe(:)
! ***** parameter for Weinert type augmented charge
integer,allocatable::nwexp(:)
! ***** the index of quantum number l and m for projectors, the index of quantum number l for projectors (starts from 1)
integer,allocatable::nlind(:,:),noind(:,:)
! ***** pseudopotential parameters, the number of Bessel functions for expansion, the number of Bessel functions for expansion of PCC charge
integer,allocatable::nqct(:),nqctpcc(:)
! ***** pseudopotential parameters for nonlocal and hard local parts and pcc charge
real*8,allocatable::coef(:,:,:)
! ***** denominator for augmented charge, see J. Math. Phys. 22, 2433 (1981)
real*8,allocatable::rfac(:,:)
! ***** <\phi|\phi> and <\phi|-1/2\nabla^2+v_loc|\phi> for PAW
real*8,allocatable::sss(:,:,:),akv(:,:,:)
! ***** \phi_smt * \phi_tru
real*8,allocatable::wail(:,:,:),wpil(:,:,:)
! ***** jellium potential
real*8,allocatable::vjell(:)

! geometry-dependend values
! ***** atomic position (coarse grid)
integer,allocatable::natx(:),naty(:),natz(:)
! ***** atomic position (dense grid)
integer,allocatable::ndatx(:),ndaty(:),ndatz(:)
! ***** atom information in the subdomain, inside subdomain or outside (coarse/dense grid)
integer,allocatable::natpri(:),natprid(:)
! ***** atom information in the subdomain, index number in the subdomain from atomic position
integer,allocatable::natpri_inf(:)
! ***** atom information in the subdomain, index number in the subdomain from cutoff radius of pseudopotential (coarse/dense grid)
integer,allocatable::naps(:),napsd(:)
! ***** atom information in the subdomain, index number in the subdomain from non-local part of pseudopotential (coarse/dense grid)
integer,allocatable::natinf(:),natinfd(:),natinfd_vloc(:)
! ***** index array for the order of the integration in MPI for the NLP
integer,allocatable::latom(:)
! ***** local potentials on coarse grid
real*8,allocatable::vloc_coarse(:,:,:)
! ***** local potentials on dense grid
real*8,allocatable::vloc_dense(:,:,:)
! ***** local soft potentials for respective atom on coarse grid
real*8,allocatable::vloc_scw(:,:,:,:)
! ***** derivative of vloc_scw
real*8,allocatable::dvlocdx_scw(:,:,:,:)
real*8,allocatable::dvlocdy_scw(:,:,:,:)
real*8,allocatable::dvlocdz_scw(:,:,:,:)
! ***** local potential on radial grid EXCEPT that of owner atom of augmented sphere
real*8,allocatable::vcorer(:,:,:)
! ***** local potential on radial grid INCLUDING that of owner atom of augmented sphere
real*8,allocatable::vcorer_all(:,:,:)
! ***** moments of charges, \int Qij |r|^l Y(\hat{r}) dr in Eq. (26) of PRB59 1758 (1999)
real*8,allocatable::qijl(:,:,:,:)
! ***** PCC charge on dense cartesian grid
real*8,allocatable::rhopcc_dense(:,:,:)
! ***** PCC charge on radial grid
real*8,allocatable::rhopccr(:,:,:)

end module var_pseudopotentials

! ========================================================================================================

module var_listvec
implicit none

! *****  the number of grid points in core region
integer:: npxmax,npymax,npzmax
! *****  dimension of list vector (coarse grid), lstvec2, and lst(x,y,z), those of dense grid
integer:: num_list,num_list_d
! ***** list vector for (x,y,z) index (coarse grid)
integer,allocatable:: lstx(:,:),lsty(:,:),lstz(:,:)
! ***** list vector for (x,y,z) index (dense grid)
integer,allocatable:: lstdx(:,:),lstdy(:,:),lstdz(:,:)
! ***** list vector for the array defined on whole region (coarse grid)
integer,allocatable:: lstvec2(:,:)
! ***** list vector for the array defined on whole region (dense grid)
integer,allocatable:: lstvecd2(:,:)
! ***** nonlocal part of pseudopotential
real*8, allocatable:: vnlocp(:,:,:)
! ***** local hard potentials for respective atom on dense grid
real*8, allocatable:: vloc_hdp(:,:)
! ***** derivative of vloc_hdp
real*8, allocatable:: dvlocdx_hdp(:,:),dvlocdy_hdp(:,:),dvlocdz_hdp(:,:)
! ***** augmented charge (dense)
real*8, allocatable:: rhoaug3d(:,:)
! ***** derivative of augmented charge (dense)
real*8, allocatable:: drhoaug3ddx(:,:),drhoaug3ddy(:,:),drhoaug3ddz(:,:)

end module var_listvec

! ========================================================================================================

module var_kslaplacian
implicit none

! ***** coefficients for finite difference
real*8,allocatable::acfd(:)

end module var_kslaplacian

! ========================================================================================================

