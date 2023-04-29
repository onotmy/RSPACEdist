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
! **********  kmesh.F90 12/14/2022-05  **********

module mod_kmesh
  implicit none
  contains
  
subroutine generate_kmesh(numkx, numky, numkz, isym, nlen, skpx, skpy, skpz, nwkp, numk)
  implicit none
  integer, intent(in) :: numkx, numky, numkz, isym, nlen
  real(8), intent(out) :: skpx(nlen), skpy(nlen), skpz(nlen)
  integer, intent(out) :: nwkp(nlen)
  integer, intent(out) :: numk
  
  call generate_uniform_grid()
  if (isym > 1) then
    ! Impose time-reversal symmetry
    call reduce_k_symmetry(0)
  end if

  contains

  ! Generate homogenious k distribution
  subroutine generate_uniform_grid()
    implicit none
    integer :: ik_count
    real(8) :: px, py, pz
    integer :: ix, iy, iz

    ik_count = 0
    do iz = 1, numkz
    do iy = 1, numky
        do ix = 1, numkx
        px = dble(ix-1) / numkx
        py = dble(iy-1) / numky
        pz = dble(iz-1) / numkz
        ! Reduce coordinate into (0.5:0.5]
        if (px > 0.5d0) px = px - 1.0d0
        if (py > 0.5d0) py = py - 1.0d0
        if (pz > 0.5d0) pz = pz - 1.0d0
        ! Calculate coordinate
        ik_count = ik_count + 1
        skpx(ik_count) = px
        skpy(ik_count) = py
        skpz(ik_count) = pz
        nwkp(ik_count) = 1
        end do
    end do
    end do
    numk = ik_count
  end subroutine generate_uniform_grid

  ! Exclude equivalent k-points on given symmetry operation
  subroutine reduce_k_symmetry(itype)
    implicit none
    integer, intent(in) :: itype
    real(8) :: skpx_tmp(nlen)
    real(8) :: skpy_tmp(nlen)
    real(8) :: skpz_tmp(nlen)
    integer :: nwkp_tmp(nlen)
    real(8) :: px, py, pz
    real(8) :: qx, qy, qz
    integer :: ik, jk, ik_sym, ik_count
    integer :: idx, idy, idz
    logical :: flag_sym

    ik_count = 0
    do ik = 1, numk
    px = skpx(ik)
    py = skpy(ik)
    pz = skpz(ik)
    ! Check equivalent point is already exists:
    flag_sym = .false.
    do jk = 1, ik_count
        ! Transformed coordinate of k point:
        qx = skpx_tmp(jk)
        qy = skpy_tmp(jk)
        qz = skpz_tmp(jk)
        select case (itype)
        case(0) ! test time reversal symmetry
        qx=-qx; qy=-qy; qz=-qz;
        ! case(1) ! test x symmetry
        !   qx=-qx; qy=+qy; qz=+qz;
        ! case(2) ! test y symmetry
        !   qx=+qx; qy=-qy; qz=+qz;
        ! case(3) ! test z symmetry
        !   qx=+qx; qy=+qy; qz=-qz;
        end select
        ! Distance between p and q:
        idx = mod(nint(numkx*(px-qx)), numkx)
        idy = mod(nint(numky*(py-qy)), numky)
        idz = mod(nint(numkz*(pz-qz)), numkz)
        if ((idx==0) .and. (idy==0) .and. (idz==0)) then
        flag_sym = .true.
        ik_sym = jk
        exit
        end if
    end do

    if (flag_sym) then
        nwkp_tmp(ik_sym) = nwkp_tmp(ik_sym) + nwkp(ik)
    else
        ik_count = ik_count + 1
        skpx_tmp(ik_count) = skpx(ik)
        skpy_tmp(ik_count) = skpy(ik)
        skpz_tmp(ik_count) = skpz(ik)
        nwkp_tmp(ik_count) = nwkp(ik)
    end if
    end do
    numk = ik_count
    skpx(1:numk) = skpx_tmp(1:numk)
    skpy(1:numk) = skpy_tmp(1:numk)
    skpz(1:numk) = skpz_tmp(1:numk)
    nwkp(1:numk) = nwkp_tmp(1:numk)
  end subroutine reduce_k_symmetry

end subroutine generate_kmesh
  
  
  
  
end module mod_kmesh
