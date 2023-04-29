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
! **********  ggartp8f.F90 06/18/2014-01  **********

module mod_vxpot_ggartp
implicit none
contains

subroutine vxpot_ggartp(num_spe,nradmx,nrad,npoint,lrhomx,ispe,yylm,dylm_dtheta,d2ylm_dtheta2,dylm_dphi,d2ylm_dphi2 & ! <
                       ,d2ylm_dtheta_dphi,wt,point,radial,dradial,rhortp                                            & ! <
                       ,agr,gggr,g2r)                                                                                 ! >
use mod_tools, only:tools_deri_grdch
integer,intent(in) :: num_spe,nradmx,nrad,npoint,lrhomx,ispe
real*8, intent(in) :: yylm(npoint,lrhomx),dylm_dtheta(npoint,lrhomx)
real*8, intent(in) :: d2ylm_dtheta2(npoint,lrhomx),dylm_dphi(npoint,lrhomx)
real*8, intent(in) :: d2ylm_dphi2(npoint,lrhomx),d2ylm_dtheta_dphi(npoint,lrhomx)
real*8, intent(in) :: wt(npoint),point(npoint,3)
real*8, intent(in) :: radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8, intent(in) :: rhortp(nradmx,npoint)
real*8, intent(out):: agr(nradmx,npoint),gggr(nradmx,npoint),g2r(nradmx,npoint)

integer :: la,il,ir
real*8  :: theta,rho,rhodr,rhodt,rhodf,rhodrr,rhodtt,rhodff,rhodtf,rhodrt,rhodrf
real*8  :: r2,r3,pi
real*8  :: sint1,sint2,tant1,grr,grt,grf,dagrr,dagrt,dagrf
real*8,allocatable::rho_m(:,:),rholm(:,:),rholmdr(:,:),rholmdrr(:,:)

  pi=dacos(-1.0d0)

  allocate(rho_m(nradmx,lrhomx),rholm(nradmx,lrhomx),rholmdr(nradmx,lrhomx),rholmdrr(nradmx,lrhomx))

  rho_m(:,:) =0.d0

!$omp parallel default(shared),private(ir)
  do la = 1,lrhomx
  do il = 1,npoint
!$omp do
  do ir = 1,nradmx
    rho_m(ir,la) = rho_m(ir,la) + rhortp(ir,il)*yylm(il,la)*wt(il)
  end do
  end do
  end do
!$omp end parallel

  rholmdr (:,:)=0.d0
  rholmdrr(:,:)=0.d0

  do la=1,lrhomx
    call tools_deri_grdch(nradmx,radial(1,ispe),dradial(1,ispe),rho_m(1,la),rholmdr(1,la),rholmdrr(1,la))
  end do

!$omp parallel default(firstprivate),shared(point,radial,agr,g2r,gggr,yylm,dylm_dtheta,d2ylm_dtheta2 &
!$omp  ,dylm_dphi,d2ylm_dphi2,d2ylm_dtheta_dphi,rho_m,rholmdr,rholmdrr)
!$omp do
  do ir=1,nradmx
    do il=1,npoint
    rho   =0.d0
    rhodr =0.d0
    rhodt =0.d0
    rhodf =0.d0
    rhodrr=0.d0
    rhodtt=0.d0
    rhodff=0.d0
    rhodtf=0.d0
    rhodrt=0.d0
    rhodrf=0.d0

      do la=1,lrhomx
           rho=     rho+             yylm(il,la)*   rho_m(ir,la)
         rhodr=   rhodr+             yylm(il,la)* rholmdr(ir,la)
         rhodt=   rhodt+      dylm_dtheta(il,la)*   rho_m(ir,la)
         rhodf=   rhodf+        dylm_dphi(il,la)*   rho_m(ir,la)
        rhodrr=  rhodrr+             yylm(il,la)*rholmdrr(ir,la)
        rhodtt=  rhodtt+    d2ylm_dtheta2(il,la)*   rho_m(ir,la)
        rhodff=  rhodff+      d2ylm_dphi2(il,la)*   rho_m(ir,la)
        rhodtf=  rhodtf+d2ylm_dtheta_dphi(il,la)*   rho_m(ir,la)
        rhodrt=  rhodrt+      dylm_dtheta(il,la)* rholmdr(ir,la)
        rhodrf=  rhodrf+        dylm_dphi(il,la)* rholmdr(ir,la)
      end do

      theta=dacos(point(il,3)/dsqrt(point(il,1)*point(il,1)+point(il,2)*point(il,2)+point(il,3)*point(il,3)))

      r2 = radial(ir,ispe)*radial(ir,ispe)
      r3 = r2*radial(ir,ispe)

      sint1 =dsin(theta)
      sint2 =sint1*sint1
      tant1 =dtan(theta)

      grr = rhodr
      grt = rhodt/radial(ir,ispe)
      grf = rhodf/(radial(ir,ispe)*sint1)
      agr(ir,il)=dsqrt(grr*grr+grt*grt+grf*grf)
      dagrr = (rhodr*rhodrr*r3+rhodt* (rhodrt*radial(ir,ispe)-rhodt)+ &
               rhodf* (rhodrf*radial(ir,ispe)-rhodf)/sint2)/agr(ir,il)/r3
      dagrt = (rhodr*rhodrt*r2+rhodt*rhodtt+ &
               rhodf* (-rhodf/tant1+rhodtf)/sint2)/ (agr(ir,il)*r3)
      dagrf = (rhodr*rhodrf*r2+rhodt*rhodtf+rhodf*rhodff/sint2)/ &
               (agr(ir,il)*r3*sint1)
      g2r(ir,il)=rhodrr+2.d0*rhodr/radial(ir,ispe)+(rhodtt+rhodt/tant1+rhodff/sint2)/r2
      gggr(ir,il)=grr*dagrr+grt*dagrt+grf*dagrf

    end do
  end do
!$omp end parallel

  deallocate(rho_m,rholm,rholmdr,rholmdrr)

  return
end subroutine vxpot_ggartp


end module
