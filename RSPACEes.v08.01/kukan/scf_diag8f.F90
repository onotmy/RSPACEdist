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
! **********  scf_diag8f.F90 06/20/2018-01  **********

module mod_scf_diag
contains

subroutine scf_diag( &
 nrc,nso,nprjmx,num_ppcell,num_spe,iscf,ncgscf,nretcg,                 & ! <
 natom,neigmx,nums,ncol,nspv,numk,nperi,northo,nf,ndisp,num_list,      & ! <
 ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                                  & ! <
 nchange,nrrz,                                                         & ! <
 key_ortho_nocmpt_innerproduct,key_ortho_cmpt_innerproduct,            & ! <
 key_natpri_in,key_natpri_inps,key_natpri_out,key_soc_calc,key_pp_paw, & ! <
 skpxx,skpyy,skpzz,dx,dy,dz,xmax,ymax,zmax,                            & ! <
 dij,dijsoc,veff,vnlocp,sss,                                           & ! <
 nkscg,nsdmax,nprecon_cg,                                              & ! <  ! only cg
 epssd,                                                                & ! <  !   "
 ndiismax,nprecon_diis,                                                & ! <  ! only diis
 eps_eig_diis,alambda_diis,ratio_diis,alambda_min,alambda_max,         & ! <  !   "
 indspe,natpri,naps,nprj,natinf,lstvec2,latom,natsoc,ntyppp,           & ! <
 lstx,lsty,lstz,natx,naty,natz,                                        & ! <
 rspsep,cspsep,svecre,sveccm,hsvre,hsvcm,ssvre,ssvcm,                  & ! X
 sval,svalsoc,cspsepsoc,svecsoc,                                       & ! >
 residual_states,ksconv,ksitmax)                                         ! >
use mod_mpi, only: myrank_glbl
use mod_stopp
use mod_kscg,         only:kscg_r,kscg_c
use mod_ksdiis,       only:ksdiis_r,ksdiis_c
use mod_reorderstates
use mod_rayleighritz, only:rayleighritz_r,rayleighritz_c,rayleighritzsoc_r,rayleighritzsoc_c
implicit none
integer,    intent(in)   ::nrc,nso,nprjmx,num_ppcell,num_spe,iscf,ncgscf,nretcg
integer,    intent(in)   ::natom,neigmx,nums,ncol,nspv,numk,nperi,northo,nf,ndisp,num_list
integer,    intent(in)   ::ncpx,ncpy,ncpz,npxmax,npymax,npzmax
integer,    intent(in)   ::nchange,nrrz
integer,    intent(in)   ::key_ortho_nocmpt_innerproduct,key_ortho_cmpt_innerproduct
integer,    intent(in)   ::key_natpri_in,key_natpri_inps,key_natpri_out,key_soc_calc,key_pp_paw
real*8,     intent(in)   ::skpxx(numk),skpyy(numk),skpzz(numk),dx,dy,dz,xmax,ymax,zmax
real*8,     intent(in)   ::dij(nprjmx,nprjmx,nums*ncol,natom)
real*8,     intent(in)   ::dijsoc(nprjmx*nso-nso+1,nprjmx*nso-nso+1,3*nso-nso+1,natom*nso-nso+1)
real*8,     intent(in)   ::veff(ncpx,ncpy,ncpz,nspv)
real*8,     intent(in)   ::vnlocp(num_list,nprjmx,num_ppcell),sss(nprjmx,nprjmx,num_spe)
integer,    intent(in)   ::nkscg,nsdmax,nprecon_cg
real*8,     intent(in)   ::epssd
integer,    intent(in)   ::ndiismax,nprecon_diis
real*8,     intent(in)   ::eps_eig_diis,alambda_diis,ratio_diis,alambda_min,alambda_max
integer,    intent(in)   ::indspe(natom),natpri(natom),naps(natom),nprj(num_spe)
integer,    intent(in)   ::natinf(natom),lstvec2(num_list,num_ppcell),latom(natom),natsoc(natom),ntyppp(num_spe)
integer,    intent(in)   ::lstx(num_list,num_ppcell),lsty(num_list,num_ppcell),lstz(num_list,num_ppcell)
integer,    intent(in)   ::natx(natom),naty(natom),natz(natom)
real*8,     intent(inout)::rspsep(nprjmx*(1-nrc)+nrc,num_ppcell*(1-nrc)+nrc, &
                                  neigmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk*(1-nrc)+nrc)
complex*16, intent(inout)::cspsep(nprjmx*nrc-nrc+1,num_ppcell*nrc-nrc+1, &
                                  neigmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1)
real*8,     intent(inout)::svecre(ncpx*(1-nrc)+nrc,ncpy*(1-nrc)+nrc,ncpz*(1-nrc)+nrc, &
                                  neigmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk*(1-nrc)+nrc)
complex*16, intent(inout)::sveccm(ncpx*nrc-nrc+1,ncpy*nrc-nrc+1,ncpz*nrc-nrc+1, &
                                  neigmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1)
real*8,     intent(inout)::hsvre(ncpx*(1-nrc)+nrc,ncpy*(1-nrc)+nrc,ncpz*(1-nrc)+nrc, &
                                 neigmx*(1-nrc)+nrc)
complex*16, intent(inout)::hsvcm(ncpx*nrc-nrc+1,ncpy*nrc-nrc+1,ncpz*nrc-nrc+1, &
                                 neigmx*nrc-nrc+1,ncol*nrc-nrc+1)
real*8,     intent(inout)::ssvre(ncpx*(1-nrc)+nrc,ncpy*(1-nrc)+nrc,ncpz*(1-nrc)+nrc, &
                                 neigmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk*(1-nrc)+nrc)
complex*16, intent(inout)::ssvcm(ncpx*nrc-nrc+1,ncpy*nrc-nrc+1,ncpz*nrc-nrc+1, &
                                 neigmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1)
real*8,     intent(out)  ::sval(neigmx,nums+1-ncol,numk)
real*8,     intent(out)  ::svalsoc((2*neigmx-1)*nso*(2-ncol)+1,1,(numk-1)*nso*(2-ncol)+1)
complex*16, intent(out)  ::cspsepsoc((nprjmx-1)*nso*(2-ncol)+1,(num_ppcell-1)*nso*(2-ncol)+1, &
                                     (2*neigmx-1)*nso*(2-ncol)+1,nso*(2-ncol)+1,(numk-1)*nso*(2-ncol)+1)
complex*16, intent(out)  ::svecsoc((ncpx-1)*nso*(2-ncol)+1,(ncpy-1)*nso*(2-ncol)+1,(ncpz-1)*nso*(2-ncol)+1, &
                                   (2*neigmx-1)*nso*(2-ncol)+1,nso*(2-ncol)+1,(numk-1)*nso*(2-ncol)+1)
real*8,     intent(out)  ::residual_states(neigmx,nums+1-ncol,numk)
logical,    intent(out)  ::ksconv
integer,    intent(out)  ::ksitmax
logical :: l_rrb,l_rra,l_reo,l_soc2nd,l_hsv
integer :: ns,nk

  l_soc2nd= (nso==1) .and. (ncol==1)                                           ! soc in second variation 
  if (nrrz/=0) then
    l_rrb= (mod(iscf,abs(nrrz))==0) .and. (nrrz>0)                             ! rayleigh-ritz before cg/diis 
    l_rra= (mod(iscf,abs(nrrz))==0) .and. (nrrz<0)                             ! rayleigh-ritz after cg/diis
    l_hsv= l_rra .or. l_soc2nd .or. ( (mod(iscf+1,nrrz)==0) .and. (nrrz>0) )   ! H*Psi output from cg/diis is used afterwards 
  else 
    l_rrb= .false.
    l_rra= .false.
    l_hsv= l_soc2nd
  endif 
  if (nchange>0) then 
    l_reo= (mod(iscf,nchange)==0) .and. (.not. l_rra) .and. (.not. l_soc2nd)   ! reorder states after cg/diis 
  else
    l_reo= .false.
  endif

  ksconv = .true.
  ksitmax= 0

  do nk=1,numk

    if (nrc==0) then 

      do ns= 1,nums

        if (l_rrb) then
!         **********  rayleigh-ritz method  *******************
          call rayleighritz_r( &
           natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,nspv, & ! <
           num_list,ncpx,ncpy,ncpz,                               & ! <
           nperi,nk,ns,nf,ndisp,                                  & ! <
           key_natpri_in,key_natpri_inps,                         & ! <
           dx,dy,dz,                                              & ! <
           nprj,natinf,lstvec2,indspe,natpri,naps,ntyppp,         & ! <
           dij,veff,vnlocp,                                       & ! <
           l_rrb,                                                 & ! <
           rspsep,svecre,ssvre,                                   & ! X
           hsvre,                                                 & ! K
           sval)                                                    ! >
!         *****************************************************
        endif ! (l_rrb)
  
!       **********  find eigenstates  **********
        if ((iscf .le. ncgscf) .or. (mod(iscf,nretcg) .eq. 0)) then
          if (myrank_glbl .eq. 0) write(ndisp,*,err=9999) 'solver for KS Eq. ==> CG'
          if (myrank_glbl .eq. 0) write(   11,*,err=9999) 'solver for KS Eq. ==> CG'
          call kscg_r( &
           natom,num_spe,neigmx,nums,nspv,numk,nperi,northo,ns,nk,nf,nkscg,nsdmax, & ! <
           nprecon_cg,nprjmx,num_list,num_ppcell,ndisp,                            & ! <
           ncpx,ncpy,ncpz,                                                         & ! <
           key_ortho_nocmpt_innerproduct,key_ortho_cmpt_innerproduct,              & ! <
           key_natpri_in,key_natpri_inps,key_natpri_out,                           & ! <
           dx,dy,dz,epssd,                                                         & ! <
           indspe,natpri,naps,nprj,natinf,lstvec2,latom,ntyppp,                    & ! <
           dij,veff,vnlocp,sss,                                                    & ! <
           ksconv,ksitmax,                                                         & ! X
           sval,rspsep,svecre,ssvre,                                               & ! X
           residual_states(1,ns,nk),                                               & ! >
           hsvre)                                                                    ! >
        else
          if (myrank_glbl .eq. 0) write(ndisp,*,err=9999) 'solver for KS Eq. ==> DIIS'
          if (myrank_glbl .eq. 0) write(   11,*,err=9999) 'solver for KS Eq. ==> DIIS'
          call ksdiis_r( &
           l_hsv,natom,num_spe,neigmx,nums,nspv,numk,nperi,ndiismax,      & ! <
           nprecon_diis,nk,ns,nf,nprjmx,num_list,num_ppcell,ndisp,northo, & ! <
           ncpx,ncpy,ncpz,                                                & ! <
           key_ortho_cmpt_innerproduct,                                   & ! <
           key_natpri_in,key_natpri_inps,                                 & ! <
           eps_eig_diis,alambda_diis,                                     & ! <
           ratio_diis,alambda_min,alambda_max,epssd,                      & ! <
           xmax,ymax,zmax,                                                & ! <
           indspe,natpri,naps,nprj,natinf,lstvec2,latom,ntyppp,           & ! <
           dij,veff,vnlocp,sss,                                           & ! <
           ksconv,ksitmax,                                                & ! X
           rspsep,sval,svecre,ssvre,                                      & ! X
           residual_states(1,ns,nk),                                      & ! X
           hsvre)                                                           ! >
        end if
!       ****************************************

        if (l_rra) then
!         **********  rayleigh-ritz method  *******************
          call rayleighritz_r( &
           natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,nspv, & ! <
           num_list,ncpx,ncpy,ncpz,                               & ! <
           nperi,nk,ns,nf,ndisp,                                  & ! <
           key_natpri_in,key_natpri_inps,                         & ! <
           dx,dy,dz,                                              & ! <
           nprj,natinf,lstvec2,indspe,natpri,naps,ntyppp,         & ! <
           dij,veff,vnlocp,                                       & ! <
           l_rrb,                                                 & ! <
           rspsep,svecre,ssvre,                                   & ! X
           hsvre,                                                 & ! K
           sval)                                                    ! >
!         *****************************************************
        endif ! (l_rra)

      enddo ! ns

      if (l_soc2nd) then 
!     **********  rayleigh-ritz method for soc in 2nd variation  **********
        call rayleighritzsoc_r( &
         natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,nspv, & ! <
         num_list,ncpx,ncpy,ncpz,                               & ! <
         nperi,nk,nf,ndisp,                                     & ! <
         key_natpri_in,key_natpri_inps,key_soc_calc,            & ! <
         dx,dy,dz,                                              & ! <
         nprj,natinf,lstvec2,indspe,natpri,naps,natsoc,ntyppp,  & ! <
         dij,dijsoc,veff,vnlocp,                                & ! <
         rspsep,cspsepsoc,svecre,svecsoc,                       & ! X
         hsvre,                                                 & ! K
         svalsoc)                                                 ! >
!     *********************************************************************
      endif

      if (l_reo) then 
!       **********  re-order eigenvalue and vectors  ********
        do ns=1,nums
          call reorderstates(nrc,num_ppcell,nprjmx,neigmx,nums,ncol,numk,ncpx,ncpy,ncpz,ns,nk,        & ! <
                            sval,residual_states,rspsep,cspsep,svecre,sveccm,ssvre,ssvcm,hsvre,hsvcm)   ! X
        end do
!       *****************************************************
      end if

    else ! nrc

      do ns= 1,nums+1-ncol

        if (l_rrb) then
!         **********  rayleigh-ritz method  *******************
          call rayleighritz_c( &
           nso,natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,ncol,nspv, & ! <
           num_list,ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                   & ! <
           nperi,nk,ns,nf,ndisp,                                           & ! <
           key_natpri_in,key_natpri_inps,key_soc_calc,                     & ! <
           dx,dy,dz,xmax,ymax,zmax,skpxx,skpyy,skpzz,                      & ! <
           nprj,natinf,lstvec2,indspe,natpri,naps,natsoc,ntyppp,           & ! <
           lstx,lsty,lstz,natx,naty,natz,                                  & ! <
           dij,dijsoc,veff,vnlocp,                                         & ! <
           l_rrb,                                                          & ! <
           cspsep,sveccm,ssvcm,                                            & ! X
           hsvcm,                                                          & ! K
           sval)                                                             ! >
!         *****************************************************
        endif ! (l_rrb) 

!       **********  find eigenstates  **********
        if ((iscf .le. ncgscf) .or. (mod(iscf,nretcg) .eq. 0)) then
          if (myrank_glbl .eq. 0) write(ndisp,*,err=9999) 'solver for KS Eq. ==> CG'
          if (myrank_glbl .eq. 0) write(   11,*,err=9999) 'solver for KS Eq. ==> CG'
          call kscg_c( &
           nso,natom,num_spe,neigmx,nums,ncol,nspv,numk,nperi,northo,ns,nk,nf,nkscg,nsdmax, & ! <
           nprecon_cg,nprjmx,num_list,num_ppcell,ndisp,                                     & ! <
           ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                                             & ! <
           key_ortho_nocmpt_innerproduct,key_ortho_cmpt_innerproduct,                       & ! <
           key_natpri_in,key_natpri_inps,key_natpri_out,key_soc_calc,                       & ! <
           dx,dy,dz,skpxx(nk),skpyy(nk),skpzz(nk),epssd,                                    & ! <
           indspe,natpri,naps,nprj,natinf,lstvec2,latom,natsoc,ntyppp,                      & ! <
           lstx,lsty,lstz,natx,naty,natz,                                                   & ! <
           dij,dijsoc,veff,vnlocp,sss,                                                      & ! <
           ksconv,ksitmax,                                                                  & ! X
           sval,cspsep,sveccm,ssvcm,                                                        & ! X
           residual_states(1,ns,nk),                                                        & ! >
           hsvcm)                                                             ! >
        else
          if (myrank_glbl .eq. 0) write(ndisp,*,err=9999) 'solver for KS Eq. ==> DIIS'
          if (myrank_glbl .eq. 0) write(   11,*,err=9999) 'solver for KS Eq. ==> DIIS'
          call ksdiis_c( &
           l_hsv,nso,natom,num_spe,neigmx,nums,ncol,nspv,numk,nperi,ndiismax, & ! <
           nprecon_diis,nk,ns,nf,nprjmx,num_list,num_ppcell,ndisp,northo,     & ! <
           ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                               & ! <
           key_ortho_cmpt_innerproduct,                                       & ! <
           key_natpri_in,key_natpri_inps,key_soc_calc,                        & ! <
           eps_eig_diis,alambda_diis,                                         & ! <
           ratio_diis,alambda_min,alambda_max,epssd,                          & ! <
           xmax,ymax,zmax,skpxx(nk),skpyy(nk),skpzz(nk),                      & ! <
           indspe,natpri,naps,nprj,natinf,lstvec2,latom,natsoc,ntyppp,        & ! <
           lstx,lsty,lstz,natx,naty,natz,                                     & ! <
           dij,dijsoc,veff,vnlocp,sss,                                        & ! <
           ksconv,ksitmax,                                                    & ! X
           cspsep,sval,sveccm,ssvcm,                                          & ! X
           residual_states(1,ns,nk),                                          & ! X
           hsvcm)                                                               ! >
        end if
!       ****************************************

        if (l_rra) then
!         **********  rayleigh-ritz method  *******************
          call rayleighritz_c( &
           nso,natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,ncol,nspv, & ! <
           num_list,ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                   & ! <
           nperi,nk,ns,nf,ndisp,                                           & ! <
           key_natpri_in,key_natpri_inps,key_soc_calc,                     & ! <
           dx,dy,dz,xmax,ymax,zmax,skpxx,skpyy,skpzz,                      & ! <
           nprj,natinf,lstvec2,indspe,natpri,naps,natsoc,ntyppp,           & ! <
           lstx,lsty,lstz,natx,naty,natz,                                  & ! <
           dij,dijsoc,veff,vnlocp,                                         & ! <
           l_rrb,                                                          & ! <
           cspsep,sveccm,ssvcm,                                            & ! X
           hsvcm,                                                          & ! K
           sval)                                                             ! >
!         *****************************************************
        end if ! (l_rra)

      enddo ! ns

      if (l_soc2nd) then 
!       **********  rayleigh-ritz method  *******************
        call rayleighritzsoc_c( &
         nso,natom,num_ppcell,num_spe,nprjmx,neigmx,numk,nums,ncol,nspv, & ! <
         num_list,ncpx,ncpy,ncpz,npxmax,npymax,npzmax,                   & ! <
         nperi,nk,nf,ndisp,                                              & ! <
         key_natpri_in,key_natpri_inps,key_soc_calc,                     & ! <
         dx,dy,dz,xmax,ymax,zmax,skpxx,skpyy,skpzz,                      & ! <
         nprj,natinf,lstvec2,indspe,natpri,naps,natsoc,ntyppp,           & ! <
         lstx,lsty,lstz,natx,naty,natz,                                  & ! <
         dij,dijsoc,veff,vnlocp,                                         & ! <
         cspsep,cspsepsoc,sveccm,svecsoc,                                & ! X
         hsvcm,                                                          & ! K
         svalsoc)                                                          ! >
!       *****************************************************
      endif

      if (l_reo) then 
!       **********  re-order eigenvalue and vectors  ********
        do ns=1,nums-ncol+1
          call reorderstates(nrc,num_ppcell,nprjmx,neigmx,nums,ncol,numk,ncpx,ncpy,ncpz,ns,nk,        & ! <
                            sval,residual_states,rspsep,cspsep,svecre,sveccm,ssvre,ssvcm,hsvre,hsvcm)   ! X
        end do
!       *****************************************************
      end if

    endif ! (nrc==0) else 

  enddo ! nk

return
  9999 call stopp ('scf_diag: error writing comments to file')
end subroutine scf_diag

end module mod_scf_diag

