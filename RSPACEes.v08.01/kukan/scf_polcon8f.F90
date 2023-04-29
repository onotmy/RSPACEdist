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
! **********  scf_polcon8f.F90 2013/06/28-01 **********

module mod_scf_polcon
implicit none
contains


subroutine scf_polcon_initialcheck( &
 natom,num_spe,ncol,npolcon,npolcon2,indspe,ntyppp,                                           & ! <
 key_polcon_occ,key_polcon_atoms,key_polcon_asa,key_polcon2_none,key_polcon2_size,key_pp_paw, & ! <
 polconocc,polconmag,tnumele)                                                                   ! <
use mod_stopp
implicit none
integer, intent(in)::natom,num_spe,ncol,npolcon,npolcon2(0:natom),indspe(natom),ntyppp(num_spe)
integer, intent(in)::key_polcon_occ,key_polcon_atoms,key_polcon_asa,key_polcon2_none,key_polcon2_size,key_pp_paw
real*8,  intent(in)::polconocc,polconmag(3,natom),tnumele
integer na
real*8  poltot

  if (npolcon==key_polcon_occ) poltot= dabs(polconocc)
  if ((npolcon==key_polcon_atoms).or.(npolcon==key_polcon_asa)) then
    poltot= 0.0d0
    do na= 1,natom
      if (npolcon2(na)==key_polcon2_size) then
        if (ncol==2) then
          poltot= poltot + dabs(polconmag(1,na))
        else
          poltot= poltot + dsqrt(polconmag(1,na)**2+polconmag(2,na)**2+polconmag(3,na)**2)
        endif
      endif
      if ((npolcon2(na)/=key_polcon2_none).and.(ntyppp(indspe(na))/=key_pp_paw)) &
       call stopp('scf_polcon_initialcheck: local constraining B-field is implemented only for PAW atoms!')
    enddo
  endif
  if (poltot>tnumele) call stopp('scf_polcon_initialcheck: total spin polarization exceeds the number of electrons!')

end subroutine scf_polcon_initialcheck


subroutine scf_polcon( &
 natom,nspv,natpri,npolcon2,key_natpri_in,key_polcon2_dir,key_polcon2_size, & ! <
 polconeta,polconmag,spinpol,                                               & ! <
 polconb)                                                                     ! X
use mod_mpi
implicit none
real*8, parameter:: eps_p2=1.0d-10
integer,intent(in)   :: natom,nspv,natpri(natom),npolcon2(0:natom)
integer,intent(in)   :: key_natpri_in,key_polcon2_dir,key_polcon2_size
real*8, intent(in)   :: polconeta(2,0:natom),polconmag(3,natom),spinpol(max(1,nspv-1),0:natom)
real*8, intent(inout):: polconb(3,natom)
logical:: somemag
integer:: na,ns
real*8 :: polmul1,polmul,pol,bnew(3)

  do na= 1,natom
    if ((npolcon2(na)==key_polcon2_dir).or.(npolcon2(na)==key_polcon2_size)) then

      if (natpri(na)==key_natpri_in) then
        polmul1= 0.0d0
        do ns= 1,nspv-1
          polmul1= polmul1 +polconmag(ns,na)**2
        enddo
        somemag= (polmul1>eps_p2)
        if (somemag) then
          polmul= 0.0d0
          do ns= 1,nspv-1
            polmul= polmul +polconmag(ns,na)*spinpol(ns,na)/polmul1
          enddo
        endif
        do ns= 1,nspv-1
          if (somemag) then
            pol= polconmag(ns,na)*polmul
          else
            pol= spinpol(ns,na)
          endif
          bnew(ns)= polconb(ns,na) +polconeta(1,0)*polconeta(1,na)*(spinpol(ns,na)-pol)
          if (npolcon2(na)==key_polcon2_size) bnew(ns)= bnew(ns) +polconeta(2,0)*polconeta(2,na)*(pol-polconmag(ns,na))
        enddo
      endif

      if (natpri(na)/=key_natpri_in) bnew(:)= 0.0d0
      call mpi_allreduce(bnew(1),polconb(1,na),max(1,nspv-1),mpi_double_precision,mpi_sum,mpicom_space,mpij)

    endif
  enddo

end subroutine scf_polcon


subroutine scf_polcon_magproj( &
 natom,num_spe,num_atcell,nspv,nradmx,npoint,natpri,natpri_inf,indspe,nradct,npolcon2, & ! <
 key_natpri_in,key_polcon2_none,                                                       & ! <
 polconmag,                                                                            & ! <
 rhotrur)                                                                                ! X
use mod_stopp
implicit none
real*8, parameter:: eps_p=1.0d-10
integer,intent(in)   :: natom,num_spe,num_atcell,nspv,nradmx,npoint
integer,intent(in)   :: natpri(natom),natpri_inf(natom),indspe(natom),nradct(num_spe)
integer,intent(in)   :: npolcon2(0:natom),key_natpri_in,key_polcon2_none
real*8, intent(in)   :: polconmag(3,natom)
real*8, intent(inout):: rhotrur(nradmx,npoint,nspv,num_atcell)
integer:: na,ipri,nrad,il,ir,ns
real*8 :: sprod,pol(3),rho(4)

  if (nspv/=4) call stopp('scf_polcon_magproj: nspv /= 4 !')

  do na= 1,natom
    if ((natpri(na)==key_natpri_in).and.(npolcon2(na)/=key_polcon2_none)) then
      ipri= natpri_inf(na)
      nrad= nradct(indspe(na))
      sprod= 0.0d0
      do ns= 1,nspv-1
        sprod= sprod +polconmag(ns,na)**2
      enddo
      sprod= dsqrt(sprod)
      if (sprod>eps_p) then
        do ns= 1,nspv-1
          pol(ns)= polconmag(ns,na)/sprod
        enddo
      else
        pol(:)= 0.0d0
      endif
      do il= 1,npoint
      do ir= 1,nrad
        do ns= 1,nspv
          rho(ns)= rhotrur(ir,il,ns,ipri)
        enddo
        rho(1)= rho(1) +rho(2)
        rho(2)= rho(1)/2.0d0 -rho(2)
        ! now: rho(1)=n , rho(2)=m_z/2 , rho(3)=m_x/2 , rho(4)=m_y/2
        sprod= 0.0d0
        do ns= 1,3
          sprod= sprod +pol(ns)*rho(ns+1)
        enddo
        do ns= 1,3
          rho(ns+1)= pol(ns)*sprod
        enddo
        rho(1)= rho(1)/2.0d0+rho(2)
        rho(2)= rho(1)-rho(2)*2.0d0
        do ns= 1,nspv
          rhotrur(ir,il,ns,ipri)= rho(ns)
        enddo
      enddo
      enddo
    endif
  enddo

end subroutine scf_polcon_magproj


end module

