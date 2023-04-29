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
! **********  filterpseudopotentials8f.F90 12/05/2022-01  **********

      module mod_filterpseudopotentials
      implicit none
      contains


subroutine filterpseudopotentials( &
 chdir,                                                                                & ! <
 ndisp,nfiltyp,num_spe,nqmx,nmesh,nradmx,nprmx,lmx,nradps,nradct,nrprj,mspe,           & ! <
 psftrad,psctrat,psext,psctoff,gridmax,filpp,rctpcc,radial,dradial,rhopcc,potc,prj,cp, & ! <
 nqct,nqctpcc,coef)                                                                      ! >
use mod_mpi
use mod_stopp
use mod_mathfunctions, only: fermidis
implicit none
character,intent(in)::chdir*200
integer,intent(in) ::ndisp,nfiltyp,num_spe,nqmx,nmesh,nradmx,nprmx,lmx
integer,intent(in) ::nradps(num_spe),nradct(num_spe),nrprj(num_spe),mspe(num_spe)
real*8, intent(in) ::psftrad,psctrat,psext,psctoff,gridmax,filpp,rctpcc
real*8, intent(in) ::radial(nradmx,num_spe),dradial(nradmx,num_spe)
real*8, intent(in) ::rhopcc(nradmx,num_spe),potc(nradmx,num_spe),prj(nradmx,nprmx*lmx,num_spe),cp(8,num_spe)
integer,intent(out)::nqct(num_spe),nqctpcc(num_spe)
real*8, intent(out)::coef(0:nqmx,0:7,num_spe)

integer::ierr
integer::ispe,i,ir,ircut,nmat,nmatpcc,iq,jq,i1,i2,i3,l,ll
real*8 ::pi,twopicbin,fourpi,r,dr
real*8 ::rcuts,qqq,dqq,qqi,qqj,z
real*8 ::besselj0,besselj1,besselj2,besselj0i,besselj1i,besselj2i,besselj0j,besselj1j,besselj2j
real*8 ::tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7
real*8 ::weight0,weight1
character cfilename*7,fname*200
integer,allocatable::ipiv(:)
real*8, allocatable::fctaa(:,:),fcta(:,:)
real*8, allocatable::fctb(:),worklu(:)
real*8, allocatable::prjps(:,:),vlocr(:)
logical :: lopen

real*8::eps ; parameter( eps=1.0d-16 )


  if (myrank_glbl .eq. 0) then

  allocate(fctaa(0:nqmx,0:nqmx),fcta(nqmx,nqmx))
  allocate(fctb(nqmx),worklu(nqmx))
  allocate(vlocr(nradmx))
  allocate(prjps(nradmx,nprmx*lmx))
  allocate(ipiv(nqmx))

  pi=dacos(-1.0d0)
  twopicbin=1.0d0/(2.0d0*pi)**3
  fourpi=4.0d0*pi

  do 1020 ispe=1,num_spe

    do i= 1,nprmx*lmx
      prjps(1:nradct(ispe),i)= prj(1:nradct(ispe),i,ispe)
      if (nradps(ispe)>nradct(ispe)) prjps(nradct(ispe)+1:nradps(ispe),i)=0.0d0
    enddo 

! ==========  correct the local potential  ==========
    do ir=2,nradct(ispe)
      r=radial(ir,ispe)
      vlocr(ir)=potc(ir,ispe)-(cp(7,ispe)*dcos(cp(6,ispe)*r)+cp(8,ispe))
    end do
    do ir=nradct(ispe)+1,nradps(ispe)
      vlocr(ir)=0.0d0
    end do
! ===================================================

! ==========  compute cutoff of the pseudopotential  ==========
    inquire(unit=12,opened=lopen)
    if (nfiltyp .eq. 0) then
      rcuts=radial(nradct(ispe),ispe)*psftrad
      nqct(ispe)=int(psctrat*rcuts/gridmax)
      nmat=int(psext*rcuts/gridmax)-nqct(ispe)
      nqctpcc(ispe)=int(psctrat*rcuts/gridmax)*nmesh
      nmatpcc=int(psext*rcuts/gridmax)*nmesh-nqctpcc(ispe)
      write (ndisp,*) 'Filtering of pseudopotentials is implemented by King-Smith method.'
      if (lopen) write (12,*) 'Filtering of pseudopotentials is implemented by King-Smith method.'
    else
      rcuts=radial(nradct(ispe),ispe)*psftrad
      nqct(ispe)=int(rcuts/gridmax)
      nmat=0
      nqctpcc(ispe)=int(rcuts/gridmax)*nmesh
      nmatpcc=0
      if (nfiltyp .eq. 1) then
        write (ndisp,*) 'Filtering of pseudopotentials is implemented by Fermi distribution.'
        if (lopen) write (12,*) 'Filtering of pseudopotentials is implemented by Fermi distribution.'
      else
        write (ndisp,*) 'Mask function for filtered pseudopotentials is not used.'
        if (lopen) write (12,*) 'Mask function for filtered pseudopotentials is not used.'
      endif
    end if
    if (nmatpcc+nqctpcc(ispe) .gt. nqmx) then
      write (ndisp,*) 'error nqmx is too small. nqmx should be >',nmatpcc+nqctpcc(ispe)
      stop
    end if
    if ((nfiltyp .eq. 0) .and. (nmat .le. 1)) then
      write (ndisp,*) 'Fourier filtering for atom number',mspe(ispe),'is switched off.'
      if (lopen) write (12,*) 'Fourier filtering for atom number',mspe(ispe),'is switched off.'
    end if
! =============================================================

! ==========  compute coefficients of Bessel function  ==========
  coef(:,:,ispe)=0.0d0
  if (nrprj(ispe) .eq. 0) then
    iq=0
    qqq=pi/rcuts*iq
    tmp0=0.0d0
    do ir=2,nradps(ispe)
      r=radial(ir,ispe)
      dr=dradial(ir,ispe)
      z=qqq*r
      besselj0=1.0d0
      tmp0=tmp0+fourpi*r*r*besselj0*vlocr(ir)*dr
    end do
    coef(iq,0,ispe)=tmp0
    do iq=1,nqct(ispe)
      qqq=pi/rcuts*iq
      tmp0=0.0d0
      do ir=2,nradps(ispe)
        r=radial(ir,ispe)
        dr=dradial(ir,ispe)
        z=qqq*r
        besselj0=dsin(z)/(z)
        tmp0=tmp0+fourpi*r*r*besselj0*vlocr(ir)*dr
      end do
      coef(iq,0,ispe)=tmp0
    end do
  end if
  if (nrprj(ispe) .eq. 2) then
    iq=0
    qqq=pi/rcuts*iq
    tmp0=0.0d0
    tmp1=0.0d0
    tmp2=0.0d0
    do ir=2,nradps(ispe)
      r=radial(ir,ispe)
      dr=dradial(ir,ispe)
      z=qqq*r
      besselj0=1.0d0
      tmp0=tmp0+fourpi*r*r*besselj0*vlocr(ir)*dr
      tmp1=tmp1+fourpi*r*besselj0*prjps(ir,1)*dr
      tmp2=tmp2+fourpi*r*besselj0*prjps(ir,2)*dr
    end do
    coef(iq,0,ispe)=tmp0
    coef(iq,1,ispe)=tmp1
    coef(iq,2,ispe)=tmp2
    do iq=1,nqct(ispe)
      qqq=pi/rcuts*iq
      tmp0=0.0d0
      tmp1=0.0d0
      tmp2=0.0d0
      do ir=2,nradps(ispe)
        r=radial(ir,ispe)
        dr=dradial(ir,ispe)
        z=qqq*r
        besselj0=dsin(z)/(z)
        tmp0=tmp0+fourpi*r*r*besselj0*vlocr(ir)*dr
        tmp1=tmp1+fourpi*r*besselj0*prjps(ir,1)*dr
        tmp2=tmp2+fourpi*r*besselj0*prjps(ir,2)*dr
      end do
      coef(iq,0,ispe)=tmp0
      coef(iq,1,ispe)=tmp1
      coef(iq,2,ispe)=tmp2
    end do
  end if
  if (nrprj(ispe) .eq. 4) then
    iq=0
    qqq=pi/rcuts*iq
    tmp0=0.0d0
    tmp1=0.0d0
    tmp2=0.0d0
    do ir=2,nradps(ispe)
      r=radial(ir,ispe)
      dr=dradial(ir,ispe)
      z=qqq*r
      besselj0=1.0d0
      tmp0=tmp0+fourpi*r*r*besselj0*vlocr(ir)*dr
      tmp1=tmp1+fourpi*r*besselj0*prjps(ir,1)*dr
      tmp2=tmp2+fourpi*r*besselj0*prjps(ir,2)*dr
    end do
    coef(iq,0,ispe)=tmp0
    coef(iq,1,ispe)=tmp1
    coef(iq,2,ispe)=tmp2
    do iq=1,nqct(ispe)
      qqq=pi/rcuts*iq
      tmp0=0.0d0
      tmp1=0.0d0
      tmp2=0.0d0
      tmp3=0.0d0
      tmp4=0.0d0
      do ir=2,nradps(ispe)
        r=radial(ir,ispe)
        dr=dradial(ir,ispe)
        z=qqq*r
        besselj0=dsin(z)/(z)
        besselj1=(dsin(z)-z*dcos(z))/(z*z)
        tmp0=tmp0+fourpi*r*r*besselj0*vlocr(ir)*dr
        tmp1=tmp1+fourpi*r*besselj0*prjps(ir,1)*dr
        tmp2=tmp2+fourpi*r*besselj0*prjps(ir,2)*dr
        tmp3=tmp3+fourpi*r*besselj1*prjps(ir,3)*dr
        tmp4=tmp4+fourpi*r*besselj1*prjps(ir,4)*dr
      end do
      coef(iq,0,ispe)=tmp0
      coef(iq,1,ispe)=tmp1
      coef(iq,2,ispe)=tmp2
      coef(iq,3,ispe)=tmp3
      coef(iq,4,ispe)=tmp4
    end do
  end if
  if (nrprj(ispe) .eq. 6) then
    iq=0
    qqq=pi/rcuts*iq
    tmp0=0.0d0
    tmp1=0.0d0
    tmp2=0.0d0
    do ir=2,nradps(ispe)
      r=radial(ir,ispe)
      dr=dradial(ir,ispe)
      z=qqq*r
      besselj0=1.0d0
      tmp0=tmp0+fourpi*r*r*besselj0*vlocr(ir)*dr
      tmp1=tmp1+fourpi*r*besselj0*prjps(ir,1)*dr
      tmp2=tmp2+fourpi*r*besselj0*prjps(ir,2)*dr
    end do
    coef(iq,0,ispe)=tmp0
    coef(iq,1,ispe)=tmp1
    coef(iq,2,ispe)=tmp2
    do iq=1,nqct(ispe)
      qqq=pi/rcuts*iq
      tmp0=0.0d0
      tmp1=0.0d0
      tmp2=0.0d0
      tmp3=0.0d0
      tmp4=0.0d0
      tmp5=0.0d0
      tmp6=0.0d0
      do ir=2,nradps(ispe)
        r=radial(ir,ispe)
        dr=dradial(ir,ispe)
        z=qqq*r
        besselj0=dsin(z)/(z)
        besselj1=(dsin(z)-z*dcos(z))/(z*z)
        besselj2=((3.0d0-z*z)*dsin(z)-3.0d0*z*dcos(z))/(z*z*z)
        tmp0=tmp0+fourpi*r*r*besselj0*vlocr(ir)*dr
        tmp1=tmp1+fourpi*r*besselj0*prjps(ir,1)*dr
        tmp2=tmp2+fourpi*r*besselj0*prjps(ir,2)*dr
        tmp3=tmp3+fourpi*r*besselj1*prjps(ir,3)*dr
        tmp4=tmp4+fourpi*r*besselj1*prjps(ir,4)*dr
        tmp5=tmp5+fourpi*r*besselj2*prjps(ir,5)*dr
        tmp6=tmp6+fourpi*r*besselj2*prjps(ir,6)*dr
      end do
      coef(iq,0,ispe)=tmp0
      coef(iq,1,ispe)=tmp1
      coef(iq,2,ispe)=tmp2
      coef(iq,3,ispe)=tmp3
      coef(iq,4,ispe)=tmp4
      coef(iq,5,ispe)=tmp5
      coef(iq,6,ispe)=tmp6
    end do
  end if
  if ((nrprj(ispe) .ne. 0) .and. (nrprj(ispe) .ne. 2) &
    .and. (nrprj(ispe) .ne. 4) .and. (nrprj(ispe) .ne. 6))  then
    call stopp('Error detected in filtering pseudopotential')
  end if
  iq=0
  qqq=pi/rcuts*iq
  tmp0=0.0d0
  do ir=2,nradps(ispe)
    r=radial(ir,ispe)
    dr=dradial(ir,ispe)
    z=qqq*r
    besselj0=1.0d0
    tmp0=tmp0+fourpi*r*r*besselj0*rhopcc(ir,ispe)*dr
  end do
  coef(iq,7,ispe)=tmp0
  do iq=1,nqctpcc(ispe)
    qqq=pi/rcuts*iq
    tmp0=0.0d0
    do ir=2,nradps(ispe)
      r=radial(ir,ispe)
      dr=dradial(ir,ispe)
      z=qqq*r
      besselj0=dsin(z)/(z)
      tmp0=tmp0+fourpi*r*r*besselj0*rhopcc(ir,ispe)*dr
    end do
    coef(iq,7,ispe)=tmp0
  end do
! ===============================================================

  if (nfiltyp .eq. 0) then
! ==========  vanish the value of PP outside of cutoff radius  ==========
  ircut=0
  do ir=2,nradps(ispe)
    if (radial(ir,ispe) .le. radial(nradct(ispe),ispe)*psctoff) ircut=ir
  end do
  if (ircut .eq. 0) call stopp('Pseudopotential data is wrong!')
  do 1010 l=0,6
    ll=-1
    if (l .eq. 0) ll=0
    if ((nrprj(ispe) .ge. 2) .and. (l .eq. 1)) ll=0
    if ((nrprj(ispe) .ge. 2) .and. (l .eq. 2)) ll=0
    if ((nrprj(ispe) .ge. 4) .and. (l .eq. 3)) ll=1
    if ((nrprj(ispe) .ge. 4) .and. (l .eq. 4)) ll=1
    if ((nrprj(ispe) .ge. 6) .and. (l .eq. 5)) ll=2
    if ((nrprj(ispe) .ge. 6) .and. (l .eq. 6)) ll=2
    if (l .eq. 7) ll=0
    if (ll .eq. -1) goto 1010
    fctaa=0.0d0
    fcta=0.0d0
    fctb=0.0d0
    dqq=pi/rcuts
    do iq=1,nmat
      qqi=dqq*(iq+nqct(ispe))
      fcta(iq+(iq-1)*nmat,1)=pi/2.0d0*qqi*qqi
    end do
    if (ll .eq. 0) then
    do jq=1,nqct(ispe)+nmat
    do iq=1,nqct(ispe)+nmat
      qqj=dqq*jq
      qqi=dqq*iq
      do ir=2,ircut
        r=radial(ir,ispe)
        dr=dradial(ir,ispe)
        z=qqi*r
        besselj0i=dsin(z)/(z)
        z=qqj*r
        besselj0j=dsin(z)/(z)
        fctaa(iq,jq)=fctaa(iq,jq)+qqi**2*qqj**2*besselj0i*besselj0j*r*r*dr
      end do
    end do
    end do
    end if
    if (ll .eq. 1) then
    do jq=1,nqct(ispe)+nmat
    do iq=1,nqct(ispe)+nmat
      qqj=dqq*jq
      qqi=dqq*iq
      do ir=2,ircut
        r=radial(ir,ispe)
        dr=dradial(ir,ispe)
        z=qqi*r
        besselj1i=(dsin(z)-z*dcos(z))/(z*z)
        z=qqj*r
        besselj1j=(dsin(z)-z*dcos(z))/(z*z)
        fctaa(iq,jq)=fctaa(iq,jq)+qqi**2*qqj**2*besselj1i*besselj1j*r*r*dr
      end do
    end do
    end do
    end if
    if (ll .eq. 2) then
    do jq=1,nqct(ispe)+nmat
    do iq=1,nqct(ispe)+nmat
      qqj=dqq*jq
      qqi=dqq*iq
      do ir=2,ircut
        r=radial(ir,ispe)
        dr=dradial(ir,ispe)
        z=qqi*r
        besselj2i=((3.0d0-z*z)*dsin(z)-3.0d0*z*dcos(z))/(z*z*z)
        z=qqj*r
        besselj2j=((3.0d0-z*z)*dsin(z)-3.0d0*z*dcos(z))/(z*z*z)
        fctaa(iq,jq)=fctaa(iq,jq)+qqi**2*qqj**2*besselj2i*besselj2j*r*r*dr
      end do
    end do
    end do
    end if
    do jq=1,nmat
    do iq=1,nmat
      fcta(iq+(jq-1)*nmat,1)=fcta(iq+(jq-1)*nmat,1)-fctaa(iq+nqct(ispe),jq+nqct(ispe))*dqq
    end do
    end do
    do jq=1,nmat
    do iq=0,nqct(ispe)
      qqj=dqq*(jq+nqct(ispe))
      qqi=dqq*iq
      fctb(jq)=fctb(jq)+fctaa(iq,jq+nqct(ispe))*coef(iq,l,ispe)*dqq
    end do
    end do
    if (nmat .gt. 1) then
      call dgetri(nmat,fcta,nmat,ipiv,worklu,nqmx,ierr)
      do jq=1,nmat
      do iq=1,nmat
        coef(jq+nqct(ispe),l,ispe)=coef(jq+nqct(ispe),l,ispe)+fctb(iq)*fcta(iq+(jq-1)*nmat,1)
      end do
      end do
    end if
 1010 continue
  fctaa=0.0d0
  fcta=0.0d0
  fctb=0.0d0
  dqq=pi/rcuts
  do iq=1,nmatpcc
    qqi=dqq*(iq+nqctpcc(ispe))
    fcta(iq+(iq-1)*nmatpcc,1)=pi/2.0d0*qqi*qqi
  end do
  do jq=1,nqctpcc(ispe)+nmatpcc
  do iq=1,nqctpcc(ispe)+nmatpcc
    qqj=dqq*jq
    qqi=dqq*iq
    do ir=2,ircut
      r=radial(ir,ispe)
      dr=dradial(ir,ispe)
      z=qqi*r
      besselj0i=dsin(z)/(z)
      z=qqj*r
      besselj0j=dsin(z)/(z)
      fctaa(iq,jq)=fctaa(iq,jq)+qqi**2*qqj**2*besselj0i*besselj0j*r*r*dr
    end do
  end do
  end do
  do jq=1,nmatpcc
  do iq=1,nmatpcc
    fcta(iq+(jq-1)*nmatpcc,1)=fcta(iq+(jq-1)*nmatpcc,1)-fctaa(iq+nqctpcc(ispe),jq+nqctpcc(ispe))*dqq
  end do
  end do
  do jq=1,nmatpcc
  do iq=0,nqctpcc(ispe)
    qqj=dqq*(jq+nqctpcc(ispe))
    qqi=dqq*iq
    fctb(jq)=fctb(jq)+fctaa(iq,jq+nqctpcc(ispe))*coef(iq,7,ispe)*dqq
  end do
  end do
  if (nmatpcc .gt. 1) then
    call dgetri(nmatpcc,fcta,nmatpcc,ipiv,worklu,nqmx,ierr)
    do jq=1,nmatpcc
    do iq=1,nmatpcc
      coef(jq+nqctpcc(ispe),7,ispe)=coef(jq+nqctpcc(ispe),7,ispe)+fctb(iq)*fcta(iq+(jq-1)*nmatpcc,1)
    end do
    end do
  end if
! =======================================================================

  nqct(ispe)=nqct(ispe)+nmat
  nqctpcc(ispe)=nqctpcc(ispe)+nmatpcc

  end if ! (nfiltyp .eq. 0)

! ==========  replicate the pseudopotentials  ==========
  cfilename(1:4)='psr.'
  i1=mspe(ispe)/100
  i2=(mspe(ispe)-i1*100)/10
  i3=(mspe(ispe)-i1*100-i2*10)
  cfilename(5:5)=char(i1+48)
  cfilename(6:6)=char(i2+48)
  cfilename(7:7)=char(i3+48)
  fname=cfilename
  if (len_trim(chdir) > 0) fname=trim(chdir)//'/'//fname
  open(11,file=fname)
  write(11,*) 'radial data'
  write(11,*) 'radius,local,s1,s2,p1,p2,d1,d2,pcc'
  do ir=2,nradps(ispe)
    r=radial(ir,ispe)
    dr=dradial(ir,ispe)
    if (nfiltyp .eq. 1) then
      call fermidis(1.0d0,r,filpp,radial(nradct(ispe),ispe), weight0)
      call fermidis(1.0d0,r,filpp,rctpcc*radial(nradct(ispe),ispe), weight1)
    else
      weight0=1.0d0
      weight1=1.0d0
    end if

    tmp0=0.0d0
    tmp1=0.0d0
    tmp2=0.0d0
    tmp3=0.0d0
    tmp4=0.0d0
    tmp5=0.0d0
    tmp6=0.0d0
    tmp7=0.0d0
    do iq=1,nqct(ispe)
      dqq=pi/rcuts
      qqq=dqq*iq
      z=qqq*r
      besselj0=dsin(z)/(z)
      besselj1=(dsin(z)-z*dcos(z))/(z*z)
      besselj2=((3.0d0-z*z)*dsin(z)-3.0d0*z*dcos(z))/(z*z*z)
      tmp0=tmp0+fourpi*coef(iq,0,ispe)*besselj0*qqq*qqq*dqq*twopicbin*weight0
      tmp1=tmp1+fourpi*coef(iq,1,ispe)*besselj0*qqq*qqq*dqq*twopicbin*weight0
      tmp2=tmp2+fourpi*coef(iq,2,ispe)*besselj0*qqq*qqq*dqq*twopicbin*weight0
      tmp3=tmp3+fourpi*coef(iq,3,ispe)*besselj1*qqq*qqq*dqq*twopicbin*weight0
      tmp4=tmp4+fourpi*coef(iq,4,ispe)*besselj1*qqq*qqq*dqq*twopicbin*weight0
      tmp5=tmp5+fourpi*coef(iq,5,ispe)*besselj2*qqq*qqq*dqq*twopicbin*weight0
      tmp6=tmp6+fourpi*coef(iq,6,ispe)*besselj2*qqq*qqq*dqq*twopicbin*weight0
    end do
    do iq=1,nqctpcc(ispe)
      dqq=pi/rcuts
      qqq=dqq*iq
      z=qqq*r
      besselj0=dsin(z)/(z)
      tmp7=tmp7+fourpi*coef(iq,7,ispe)*besselj0*qqq*qqq*dqq*twopicbin*weight1
    end do
    if (r .lt. rcuts) then
      write(11,1111) r,tmp0,tmp1*r,tmp2*r,tmp3*r,tmp4*r,tmp5*r,tmp6*r,tmp7, &
       vlocr(ir),prjps(ir,1),prjps(ir,2),prjps(ir,3),prjps(ir,4),prjps(ir,5),prjps(ir,6),rhopcc(ir,ispe)
    end if
  end do
  write(11,*) 'coefficients'
  write(11,*) 'wave vector,local,s1,s2,p1,p2,d1,d2,pcc'
  do iq=0,nqctpcc(ispe)
    dqq=pi/rcuts
    qqq=dqq*iq
    write(11,1111) qqq,(coef(iq,i,ispe),i=0,7)
  end do
! ======================================================
  write(11,*) 'cutoff(high)=',dqq*nqct(ispe),'(a.u.-1)','cutoff(low)=',dqq*(nqct(ispe)-nmat),'(a.u.-1)'
  close(11)
! ==========  correct the local potential and pcc charge  ==========
!  do ir=2,nradps(ispe)
!    r=radial(ir,ispe)
!    tmp0=0.0d0
!    tmp7=0.0d0
!    do iq=1,nqct(ispe)
!      dqq=pi/rcuts
!      qqq=dqq*iq
!      z=qqq*r
!      besselj0=dsin(z)/(z)
!      tmp0=tmp0+fourpi*coef(iq,0,ispe)*besselj0*qqq*qqq*dqq*twopicbin
!    end do
!    do iq=1,nqctpcc(ispe)
!      dqq=pi/rcuts
!      qqq=dqq*iq
!      z=qqq*r
!      besselj0=dsin(z)/(z)
!      tmp7=tmp7+fourpi*coef(iq,7,ispe)*besselj0*qqq*qqq*dqq*twopicbin
!    end do
!    vlocr(ir,ispe)=tmp0
!    rhopcc(ir,ispe)=tmp7
!  end do
! ==================================================================

 1020 continue ! ispe=1,num_spe
 1111 format(20e16.7)

  deallocate(fctaa,fcta)
  deallocate(fctb,worklu)
  deallocate(vlocr,prjps)
  deallocate(ipiv)

  endif ! (myrank_glbl .eq. 0)

  call mpi_bcast(coef,(nqmx+1)*8*num_spe,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(nqct,num_spe,mpi_integer,0,mpicom_space,mpij)
  call mpi_bcast(nqctpcc,num_spe,mpi_integer,0,mpicom_space,mpij)

end subroutine filterpseudopotentials

end module
