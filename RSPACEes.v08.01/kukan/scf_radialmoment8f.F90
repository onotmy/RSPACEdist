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
! **********  scf_radialmoment8e.f90 02/24/2012-01  **********

module mod_scf_radialmoment
implicit none

contains

subroutine scf_radialmoment( &
 key_natpri_in,key_pp_paw,                                 & ! <
 natom,num_spe,num_atcell,nradmx,npoint,lrhomx,nspv,nperi, & ! <
 indspe,natpri,natpri_inf,ntyppp,nradct,lpmx,              & ! <
 xmax,ymax,zmax,tnumele,                                   & ! <
 yylm,wt,                                                  & ! <
 rhotrur,rhosmtr,rhoaugr,                                  & ! <
 rhotrum,rhosmtm,rhoaugm)                                    ! >
implicit none
integer,intent(in) ::key_natpri_in,key_pp_paw
integer,intent(in) ::natom,num_spe,num_atcell,nradmx,npoint,lrhomx,nspv,nperi
integer,intent(in) ::indspe(natom),natpri(natom),natpri_inf(natom)
integer,intent(in) ::ntyppp(num_spe),nradct(num_spe),lpmx(num_spe)
real*8, intent(in) ::xmax,ymax,zmax
real*8, intent(in) ::tnumele
real*8, intent(in) ::yylm(npoint,lrhomx), wt(npoint)
real*8, intent(in) ::rhotrur(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::rhosmtr(nradmx,npoint,nspv,num_atcell)
real*8, intent(in) ::rhoaugr(nradmx,npoint,num_atcell)
real*8, intent(out)::rhotrum(nradmx,lrhomx,num_atcell)
real*8, intent(out)::rhosmtm(nradmx,lrhomx,num_atcell)
real*8, intent(out)::rhoaugm(nradmx,lrhomx,num_atcell)
integer:: na,ipri,ns,la,il,ir, ispe,lrho
real*8 :: omega,tmp,tmpyw

  tmp=0.0d0
  if (nperi==3) then
    omega=8.0d0*xmax*ymax*zmax
    tmp=-tnumele/omega
  end if

  do na=1,natom
    if ((ntyppp(indspe(na)) .eq. key_pp_paw) .and. (natpri(na) .eq. key_natpri_in)) then
      ipri=natpri_inf(na)
      ispe=indspe(na)
      lrho=(2*lpmx(ispe)+1)**2
!$omp do
      do ir=1,nradmx*lrhomx
        rhosmtm(ir,1,ipri)=0.0d0
        rhotrum(ir,1,ipri)=0.0d0
        rhoaugm(ir,1,ipri)=0.0d0
      end do
      do ns=1,min(2,nspv)
        if (ns==1) then
!$omp do
          do la=1,lrho
            do il=1,npoint
              tmpyw= yylm(il,la)*wt(il)
              do ir=2,nradct(ispe)-1
                rhoaugm(ir,la,ipri)= rhoaugm(ir,la,ipri) + rhoaugr(ir,il,ipri)*tmpyw
                rhosmtm(ir,la,ipri)= rhosmtm(ir,la,ipri) + rhosmtr(ir,il,ns,ipri)*tmpyw + tmp*tmpyw
                rhotrum(ir,la,ipri)= rhotrum(ir,la,ipri) + rhotrur(ir,il,ns,ipri)*tmpyw + tmp*tmpyw
              end do
            end do
          end do
        else
!$omp do
          do la=1,lrho
            do il=1,npoint
              tmpyw= yylm(il,la)*wt(il)
              do ir=2,nradct(ispe)-1
                rhosmtm(ir,la,ipri)= rhosmtm(ir,la,ipri) + rhosmtr(ir,il,ns,ipri)*tmpyw
                rhotrum(ir,la,ipri)= rhotrum(ir,la,ipri) + rhotrur(ir,il,ns,ipri)*tmpyw
              end do
            end do
          end do
        end if
      end do
    end if
  end do
  return

end subroutine


end module
