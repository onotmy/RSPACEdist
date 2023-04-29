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
!     **********  keys8f.F90 03/15/2013  **********

      module var_keys

! *** type of pseudpotential
!   * pseudopotential is NCPS
      integer, parameter, public :: key_pp_ncps = 0
!   * pseudopotential is PAW
      integer, parameter, public :: key_pp_paw  = 1

! *** atom's relation to real-space decomposition
!   * nucleus is inside sub-domain
      integer, parameter, public :: key_natpri_in   = 2
!   * nonlocal region of pp is overlaped in sub-domain
      integer, parameter, public :: key_natpri_inps = 1
!   * nonlocal region of pp is not overlaped in sub-domain
      integer, parameter, public :: key_natpri_out  = 0

! *** symmetry
!   * symmetricity any
      integer, parameter, public :: key_sym_any = 0
!   * symmetricity bcc
      integer, parameter, public :: key_sym_bcc = 1
!   * symmetricity fcc
      integer, parameter, public :: key_sym_fcc = 2
!   * symmetricity dia
      integer, parameter, public :: key_sym_dia = 3
!   * symmetricity hcp
      integer, parameter, public :: key_sym_hcp = 4

! *** symmetry of k-points
!   * no symmetry specified
      integer, parameter, public :: key_ksym_none = 0
!   * half Brillouin zone filled 
      integer, parameter, public :: key_ksym_inv  = 1
!   * order of k-points is +k1,+k2,+k3,..,+kn,-kn,-k(n-1),...,-k2,-k1 
      integer, parameter, public :: key_ksym_inv2 = 2

! *** jellium potential
!   * don't use jellium potential
      integer, parameter, public :: key_jel_nocalc = 0
!   * use jellium potential
      integer, parameter, public :: key_jel_calc   = 1

! *** constrained magnetization
!   * don't constrain
      integer, parameter, public :: key_polcon_none  = 0
!   * constrain atoms' magnetic moments via constraining B-fields
      integer, parameter, public :: key_polcon_atoms = 1
!   * as prev. line, additionally project mag. in the spheres on constr. direction
      integer, parameter, public :: key_polcon_asa   = 2
!   * constrain total spin polarization via occupation numbers
      integer, parameter, public :: key_polcon_occ   = 3
!   * update constraining fields
      integer, parameter, public :: key_polcon2_new  = 0
!   * update constraining fields, rotate initial mag. moments before 1st iteration
      integer, parameter, public :: key_polcon2_rot  = 1
!   * specific atom: don't constrain
      integer, parameter, public :: key_polcon2_none = 0
!   * specific atom: constrain magnetization direction
      integer, parameter, public :: key_polcon2_dir  = 1
!   * specific atom: constrain magnetization direction and magnetic moment
      integer, parameter, public :: key_polcon2_size = 2
!   * fixed constraining field
      integer, parameter, public :: key_polcon2_fix  = 3

! *** spin-orbit coupling
!   * do not calculate this atom
      integer, parameter, public :: key_soc_nocalc = 0
!   * calculate this atom
      integer, parameter, public :: key_soc_calc   = 1

! *** band structure files
!   * do not create "band*.dat" files  
      integer, parameter, public :: key_band_no  = 0
!   * write single-particle eigenvalues (sval) vs. k-point indices (nk) in "band*.dat" files  
      integer, parameter, public :: key_band_yes = 1

! *** other
! ***
!   * augcharge compute moment, charge, and derivative of charge
      integer, parameter, public :: key_augcharge_cmpt_all        = 0
!   * augcharge compute moment and charge
      integer, parameter, public :: key_augcharge_cmpt_mom_charge = 1
! ***      
      integer, parameter, public :: key_ortho_nocmpt_innerproduct = 0
      integer, parameter, public :: key_ortho_cmpt_innerproduct   = 1


      end module var_keys
