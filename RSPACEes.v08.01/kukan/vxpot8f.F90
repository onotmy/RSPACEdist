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
! **********  scf_vxpot8f.F90 07/06/2017-01  **********

module mod_vxpot
implicit none
real*8,parameter::rhocut=1.0d-12
contains

! ===================================================================
subroutine vxpot_lda_vwn(n,ppp,vxc,exc)
implicit real*8 (a-h,o-z)
parameter (tol_15=1.0d-15)
integer n
integer i
real*8 exc(n),vxc(n)
real*8 ppp(n)

data ap,xp0,bp,cp,qp,cp1,cp2,cp3/0.0621814d0,-0.10498d0,3.72744d0, &
     12.9352d0,6.1519908d0,1.2117833d0,1.1435257d0,-0.031167608d0/
data three/3.0d0/

  fpi    = 16.0d0 * atan(1.0d0)
  thfpi  = three / fpi
  onthrd = 1.0d0/3.0d0

!ocl norecurrence(vxc,exc)
!$omp do
  do i=1,n
! **********  for LDA  **********
    rho0=ppp(i)
    rho0 = max(tol_15,rho0)
    rs = (thfpi/rho0)**onthrd
    x = sqrt(rs)
    xpx = x*x + bp*x + cp
    beta = 1.d0/ (2.74208d0+3.182d0*x+0.09873d0*x*x+0.18268d0*x**3)
    dbeta = - (0.27402d0*x+0.09873d0+1.591d0/x)*beta**2
    atnp = atan(qp/ (2.d0*x+bp))
    ecp = ap* (log(x*x/xpx)+cp1*atnp-cp3* (log((x-xp0)**2/xpx)+cp2*atnp))
    ec = ecp
    ex = -0.9163306d0/rs
    exc(i) = ec + ex
    tp1    = (x*x+bp*x)/xpx
    ucp    = ecp - ap/3.d0* (1.d0-tp1-cp3* (x/ (x-xp0)-tp1-xp0*x/xpx))
    vxc(i) = ucp - 1.221774d0/rs
    vxc(i) = vxc(i)*0.5d0
    exc(i) = exc(i)*0.5d0
! *******************************
  end do

  return
end subroutine vxpot_lda_vwn

! ===================================================================

subroutine vxpot_lda_vwns(n,ppp,vxc,exc)
implicit real*8 (a-h,o-z)
parameter (tol_15=1.0d-15)
integer n
integer i
real*8 exc(n),vxc(n,2)
real*8 ppp(n,2)

  data ap,xp0,bp,cp,qp,cp1,cp2,cp3/0.0621814d0,-0.10498d0,3.72744d0, &
       12.9352d0,6.1519908d0,1.2117833d0,1.1435257d0,-0.031167608d0/
  data af,xf0,bf,cf,qf,cf1,cf2,cf3/0.0310907d0,-0.32500d0,7.06042d0, &
       18.0578d0,4.7309269d0,2.9847935d0,2.7100059d0,-0.1446006d0/
  data one,three/1.0d0,3.0d0/

  fpi    = 16.0d0 * atan(1.0d0)
  thfpi  = three / fpi
  onthrd = 1.0d0/3.0d0
  tothrd = onthrd + onthrd
  fothrd = tothrd + tothrd

!ocl norecurrence(vxc,exc)
!$omp do
  do i=1,n
    rho=ppp(i,1)+ppp(i,2)
    ! rhomag=dabs(ppp(i,1)-ppp(i,2))
    rhomag=-ppp(i,1)+ppp(i,2) ! revised on 17.Aug.2007 by Ono
    rho = max(tol_15,rho)
    smag = sign(one,rhomag)
    rhomag = smag*min(rho-tol_15,abs(rhomag))
    rs = (thfpi/rho)**onthrd
    s = rhomag/rho
    x = sqrt(rs)
    xpx = x*x + bp*x + cp
    xfx = x*x + bf*x + cf
    s4 = s**4 - 1.d0
    cbrt1 = (1.d0-s) ** onthrd
    cbrt2 = (1.d0+s) ** onthrd
    fs = ((1.d0+s) ** fothrd + (1.d0-s) ** fothrd -2.d0)/(2.d0 ** fothrd -2.d0)
    beta = 1.d0/ (2.74208d0+3.182d0*x+0.09873d0*x*x+0.18268d0*x**3)
    dfs = fothrd * (cbrt1-cbrt2)/ (2.d0 ** fothrd - 2.d0)
    dbeta = - (0.27402d0*x+0.09873d0+1.591d0/x)*beta**2
    atnp = atan(qp/ (2.d0*x+bp))
    atnf = atan(qf/ (2.d0*x+bf))
    ecp = ap* (log(x*x/xpx)+cp1*atnp-cp3* (log((x-xp0)**2/xpx)+cp2*atnp))
    ecf = af* (log(x*x/xfx)+cf1*atnf-cf3* (log((x-xf0)**2/xfx)+cf2*atnf))
    ec = ecp + fs* (ecf-ecp)* (1.d0+s4*beta)
    ex = -0.9163306d0/rs* (1.d0 + 0.25992105d0 * fs)
    exc(i) = ec + ex
    tp1 = (x*x+bp*x)/xpx
    tf1 = (x*x+bf*x)/xfx
    ucp = ecp - ap/3.d0* (1.d0-tp1-cp3* (x/ (x-xp0)-tp1-xp0*x/xpx))
    ucf = ecf - af/3.d0* (1.d0-tf1-cf3* (x/ (x-xf0)-tf1-xf0*x/xfx))
    uc0 = ucp + (ucf-ucp)*fs
    uc10 = uc0 + (ecf-ecp)* (s+1.d0)*dfs
    uc20 = uc0 + (ecf-ecp)* (s-1.d0)*dfs
    duc = (ucf-ucp)*beta*s4*fs + (ecf-ecp)* (-rs/3.d0)*dbeta*s4*fs
    duc1 = duc + (ecf-ecp)*beta* (s+1.d0)* (4.d0*s**3*fs+s4*dfs)
    duc2 = duc + (ecf-ecp)*beta* (s-1.d0)* (4.d0*s**3*fs+s4*dfs)
    uc1 = uc10 + duc1
    uc2 = uc20 + duc2
    vxc(i,1) = uc1 - 1.221774d0/rs*cbrt1
    vxc(i,2) = uc2 - 1.221774d0/rs*cbrt2
    vxc(i,1) = vxc(i,1)*0.5d0
    vxc(i,2) = vxc(i,2)*0.5d0
    exc(i)   = exc(i)*0.5d0
  end do

  return
end subroutine vxpot_lda_vwns

! ===================================================================

subroutine vxpot_lda_pz(n,rho,vx,ex)
implicit real*8 (a-h,o-z)
parameter (tol_15=1.0d-15)
integer n,i
real*8 rho(n),vx(n),ex(n)

  third=1.0d0/3.0d0
  third4=4.0d0/3.0d0
  pi=dacos(-1.0d0)
!$omp do
  do i=1,n
    if (rho(i) .gt. tol_15) then
      rs=(3.0d0/(4.0d0*pi*(rho(i))))**(1.0d0/3.0d0)
      if (rs.ge.1.0d0) then
        exc=-0.1423d0/(1.0d0+1.0529d0*sqrt(rs)+0.3334d0*rs)
        vx(i)=exc-rs/3.0d0*(0.1423d0*(0.3334d0+0.52645d0/sqrt(rs))/(1.0d0+1.0529d0*sqrt(rs)+0.3334d0*rs)**2.0d0)
      else
        exc=-0.048d0+0.0311d0*log(rs)-0.0116d0*rs+0.002d0*rs*log(rs)
        vx(i)=exc-rs/3.0d0*(0.0311d0/rs-0.0096d0+0.002d0*log(rs))
      end if
      exc=exc-0.75d0*(2.25d0/pi**2)**third/rs
      ex(i)=exc
      vx(i)=vx(i)-0.75d0*(2.25d0/pi**2)**third/rs-rs/3.0d0*(0.75d0*(2.25d0/pi**2)**third/rs**2.0d0)
    else
      vx(i)=0.0d0
      ex(i)=0.0d0
    end if
  end do
  return 
end subroutine vxpot_lda_pz


subroutine vxpot_lda_pzs(n,rho,vx,ex)
implicit real*8 (a-h,o-z)
parameter (tol_15=1.0d-15)
integer n,i
real*8:: rho(n,2),vx(n,2),ex(n)

  third=1.0d0/3.0d0
  third4=4.0d0/3.0d0
  pi=dacos(-1.0d0)
!$omp do
  do i=1,n
    ppp=rho(i,1)+rho(i,2)
    if ((ppp.gt.tol_15).and.(rho(i,1).gt.tol_15).and.(rho(i,2).gt.tol_15)) then
      xi=(rho(i,1)-rho(i,2))/ppp
      fxi= ((1.0d0+xi)**third4+(1.0d0-xi)**third4-2.0d0)/(2.0d0*(2.0d0**third-1.0d0))
      dfxi=(2.0d0*(1.0d0+xi)**third-2.0d0*(1.0d0-xi)**third)/(2.0d0**third-1.0d0)/3.0d0
      rs=(3.0d0/(4.0d0*pi*ppp))**third
      drs=-(0.75d0/pi)**third/3.0d0*ppp**(-third4)
      eexp=-0.75d0*(2.25d0/pi**2)**third/rs
      eexf=-0.75d0*(4.50d0/pi**2)**third/rs
      deexp=0.75d0*(2.25d0/pi**2)**third/rs**2
      deexf=0.75d0*(4.50d0/pi**2)**third/rs**2
      if (rs .ge. 1.0d0) then
        eecp=-0.1423d0/(1.0d0+1.0529d0*dsqrt(rs)+0.3334d0*rs)
        eecf=-0.0843d0/(1.0d0+1.3981d0*dsqrt(rs)+0.2611d0*rs)
        deecp=0.1423d0*(0.3334d0+0.52645d0/dsqrt(rs))/(1.0d0+1.0529d0*dsqrt(rs)+0.3334d0*rs)**2
        deecf=0.0843d0*(0.2611d0+0.69905d0/dsqrt(rs))/(1.0d0+1.3981d0*dsqrt(rs)+0.2611d0*rs)**2
      else
        eecp=0.03110d0*log(rs)-0.0480d0+0.0020d0*rs*log(rs)-0.0116d0*rs
        eecf=0.01555d0*log(rs)-0.0269d0+0.0007d0*rs*log(rs)-0.0048d0*rs
        deecp=0.03110d0/rs+0.0020d0*log(rs)-0.0096d0
        deecf=0.01555d0/rs+0.0007d0*log(rs)-0.0041d0
      end if
      eexcp=eexp+eecp
      eexcf=eexf+eecf
      deexcp=deexp+deecp
      deexcf=deexf+deecf
      deexc=deexcp+(deexcf-deexcp)*fxi
      dexcxi=(eexcf-eexcp)*dfxi
      ex(i)=eexcp+(eexcf-eexcp)*fxi
      vx(i,1)=ex(i)+ppp*deexc*drs+dexcxi*(1.0d0-xi)
      vx(i,2)=ex(i)+ppp*deexc*drs-dexcxi*(1.0d0+xi)
    else
      vx(i,1)=0.0d0
      vx(i,2)=0.0d0
      ex(i)=0.0d0
    end if
  end do
  return
end subroutine vxpot_lda_pzs

! ===================================================================

subroutine vxpot_gga_pw91(n,rho,abgr,ggabg,ggr,vx,exu)
implicit real*8 (a-h,o-z)
integer,intent(in) :: n
real*8 ,intent(in) :: rho(n),abgr(n),ggabg(n),ggr(n)
real*8 ,intent(out):: vx(n),exu(n)
integer :: i

!$omp do
  do i=1,n
    up=rho(i)*0.5d0
    dn=rho(i)*0.5d0
    agrup=abgr(i)*0.5d0
    agrdn=abgr(i)*0.5d0
    delgrup=ggabg(i)*0.25d0
    delgrdn=ggabg(i)*0.25d0
    uplap=ggr(i)*0.5d0
    dnlap=ggr(i)*0.5d0
    agr=abgr(i)
    delgr=ggabg(i)
    thrd=1.d0/3.d0;thrd2=2.d0*thrd
    pi32=29.608813203268075856503472999628d0
    pi=3.1415926535897932384626433832795d0
    alpha=1.91915829267751300662482032624669d0
!**********  EXCHANGE PART  **********
    rho2=2.d0*up
    if(rho2.gt.rhocut)then
      fk=(pi32*rho2)**thrd
      s=2.d0*agrup/(2.d0*fk*rho2)
      u=4.d0*delgrup/(rho2*rho2*(2.d0*fk)**3)
      v=2.d0*uplap/(rho2*(2.d0*fk)**2)
!----------  subroutine exchpw91  ----------
      a1=0.19645D0;a2=0.27430D0;a3=0.15084D0;a4=100.d0
      ax=-0.7385588D0;a=7.7956D0;b1=0.004d0
      thrd4=1.33333333333D0
! for Becke exchange, set a3=b1=0
      FAC = AX*RHO2**THRD
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.D0/DSQRT(1.D0+A*A*S2)
      P1 = DLOG(A*S+1.D0/P0)
      P2 = DEXP(-A4*S2)
      P3 = 1.D0/(1.D0+A1*S*P1+B1*S4)
      P4 = 1.D0+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      exuppw91 = FAC*F
!  ENERGY DONE. NOW THE POTENTIAL:
      P5 = B1*S2-(A2-A3*P2)
      P6 = A1*S*(P1+A*S*P0)
      P7 = 2.D0*(A2-A3*P2)+2.D0*A3*A4*S2*P2-4.D0*B1*S2*F
      FS = P3*(P3*P5*P6+P7)
      P8 = 2.D0*S*(B1-A3*A4*P2)
      P9 = A1*P1+A*A1*S*P0*(3.D0-A*A*S2*P0*P0)
      P10 = 4.D0*A3*A4*S*P2*(2.D0-A4*S2)-8.D0*B1*S*F-4.D0*B1*S3*FS
      P11 = -P3*P3*(A1*P1+A*A1*S*P0+4.D0*B1*S3)
      FSS = P3*P3*(P5*P9+P6*P8)+2.D0*P3*P5*P6*P11+P3*P10+P7*P11
      vxuppw91 = FAC*(THRD4*F-(U-THRD4*S3)*FSS-V*FS)
!-------------------------------------------
    else
      exuppw91=0.d0
      vxuppw91=0.d0
    endif
! construct total density and contribution to ex
    rhot=up+dn
    exdnpw91=exuppw91
    vxdnpw91=vxuppw91
    expw91=(exuppw91*up+exdnpw91*dn)/rhot
!*************************************
!**********  CORRELATION PART  **********
    if(rhot.gt.rhocut) then
      zet=(up-dn)/rhot
      g=((1.d0+zet)**thrd2+(1.d0-zet)**thrd2)/2.d0
      fk=(pi32*rhot)**thrd
      rs=alpha/fk
      sk=sqrt(4.d0*fk/pi)
      twoksg=2.d0*sk*g
      t=agr/(twoksg*rhot)
      uu=delgr/(rhot*rhot*twoksg**3)
      rholap=uplap+dnlap
      vv=rholap/(rhot*twoksg**2)
      ww=(agrup**2-agrdn**2-zet*agr**2)/(rhot*rhot*twoksg**2)
!==========  subroutine corlsd  ==========
      gam=0.5198421D0;fzz=1.709921D0
      thrd4=1.333333333333D0
      F = ((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
!----------  subroutine corlsd  ----------
      A=0.0310907D0;A1=0.21370D0;B1=7.5957D0;B2=3.5876D0;B3=1.6382D0;B4=0.49294D0;P=1.00D0
      P1 = P + 1.D0
      Q0 = -2.D0*A*(1.D0+A1*RS)
      RS12 = DSQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      eu=gg;eurs=ggrs
!-----------------------------------------
!----------  subroutine corlsd  ----------
      A=0.01554535D0;A1=0.20548D0;B1=14.1189D0;B2=6.1977D0;B3=3.3662D0;B4=0.62517D0;P=1.00D0
      P1 = P + 1.D0
      Q0 = -2.D0*A*(1.D0+A1*RS)
      RS12 = DSQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      ep=gg;eprs=ggrs
!-----------------------------------------
!----------  subroutine corlsd  ----------
      A=0.0168869D0;A1=0.11125D0;B1=10.357D0;B2=3.6231D0;B3=0.88026D0;B4=0.49671D0;P=1.00D0
      P1 = P + 1.D0
      Q0 = -2.D0*A*(1.D0+A1*RS)
      RS12 = DSQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      alfm=gg;alfrsm=ggrs
!-----------------------------------------
!  ALFM IS MINUS THE SPIN STIFFNESS ALFC
      ALFC = -ALFM
      Z4 = ZET**4
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
!  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
              -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
!=========================================
!----------  subroutine CORpw91  ----------
      xnu=15.75592D0;cc0=0.004235D0;cx=-0.001667212D0
      alf=0.09D0
      c1=0.002568D0;c2=0.023266D0;c3=7.389D-6;c4=8.723D0;c5=0.472D0;c6=7.389D-2;a4=100.D0
      thrdm=-0.333333333333D0
      BET = XNU*CC0
      DELT = 2.D0*ALF/BET
      G3 = G**3
      G4 = G3*G
      PON = -DELT*EC/(G3*BET)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      Q6 = C1+C2*RS+C3*RS2
      Q7 = 1.D0+C4*RS+C5*RS2+C6*RS3
      CC = -CX + Q6/Q7
      R0 = 0.663436444d0*rs
      R1 = A4*R0*G4
      COEFF = CC-CC0-3.D0*CX/7.D0
      R2 = XNU*COEFF*G3
      R3 = DEXP(-R1*T2)
      H0 = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0+H1
      CCRS = (C2+2.*C3*RS)/Q7 - Q6*(C4+2.*C5*RS+3.*C6*RS2)/Q7**2
      RSTHRD = RS/3.D0
      R4 = RSTHRD*CCRS/COEFF
      GZ = ((1.D0+ZET)**THRDM - (1.D0-ZET)**THRDM)/3.D0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      H0B = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      H0RS = -RSTHRD*H0B*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      H0BT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      H0RST = RSTHRD*T2*H0BT*BEC*ECRS
      H0Z = 3.D0*GZ*H0/G + H0B*(BG*GZ+BEC*ECZET)
      H0T = 2.*BET*G3*Q9/Q8
      H0ZT = 3.D0*GZ*H0T/G+H0BT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      H0TT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      H1RS = R3*R2*T2*(-R4+R1*T2/3.D0)
      FACT4 = 2.D0-R1*T2
      H1RST = R3*R2*T2*(2.D0*R4*(1.D0-R1*T2)-THRD2*R1*T2*FACT4)
      H1Z = GZ*R3*R2*T2*(3.D0-4.D0*R1*T2)/G
      H1T = 2.D0*R3*R2*(1.D0-R1*T2)
      H1ZT = 2.D0*GZ*R3*R2*(3.D0-11.D0*R1*T2+4.D0*R1*R1*T4)/G
      H1TT = 4.D0*R3*R2*R1*T*(-2.D0+R1*T2)
      HRS = H0RS+H1RS
      HRST = H0RST+H1RST
      HT = H0T+H1T
      HTT = H0TT+H1TT
      HZ = H0Z+H1Z
      HZT = H0ZT+H1ZT
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
!------------------------------------------
      ecpw91=ec+h
      vcuppw91=vcup+dvcup
    else
      ecpw91=0.0d0
      vcuppw91=0.0d0
    end if
!****************************************
    exu(i)=expw91+ecpw91
    vx(i)=vxuppw91+vcuppw91
  end do
end subroutine vxpot_gga_pw91


subroutine vxpot_gga_pw91s(n,rho,abgr,ggabg,ggr,vx,exu)
implicit real*8 (a-h,o-z)
integer,intent(in) :: n
real*8 ,intent(in) :: rho(n,2),abgr(n,3),ggabg(n,3),ggr(n,3)
real*8 ,intent(out):: vx(n,2),exu(n)
integer :: i

!$omp do
  do i=1,n
    up=rho(i,1)
    dn=rho(i,2)
    agrup=abgr(i,1)
    agrdn=abgr(i,2)
    delgrup=ggabg(i,1)
    delgrdn=ggabg(i,2)
    uplap=ggr(i,1)
    dnlap=ggr(i,2)
    agr=abgr(i,3)
    delgr=ggabg(i,3)
    thrd=1.d0/3.d0;thrd2=2.d0*thrd
    pi32=29.608813203268075856503472999628d0
    pi=3.1415926535897932384626433832795d0
    alpha=1.91915829267751300662482032624669d0
!**********  EXCHANGE PART  **********
    rho2=2.d0*up
    if(rho2.gt.rhocut)then
      fk=(pi32*rho2)**thrd
      s=2.d0*agrup/(2.d0*fk*rho2)
      u=4.d0*delgrup/(rho2*rho2*(2.d0*fk)**3)
      v=2.d0*uplap/(rho2*(2.d0*fk)**2)
!----------  subroutine exchpw91  ----------
      a1=0.19645D0;a2=0.27430D0;a3=0.15084D0;a4=100.d0
      ax=-0.7385588D0;a=7.7956D0;b1=0.004d0
      thrd4=1.33333333333D0
! for Becke exchange, set a3=b1=0
      FAC = AX*RHO2**THRD
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.D0/DSQRT(1.D0+A*A*S2)
      P1 = DLOG(A*S+1.D0/P0)
      P2 = DEXP(-A4*S2)
      P3 = 1.D0/(1.D0+A1*S*P1+B1*S4)
      P4 = 1.D0+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      exuppw91 = FAC*F
!  ENERGY DONE. NOW THE POTENTIAL:
      P5 = B1*S2-(A2-A3*P2)
      P6 = A1*S*(P1+A*S*P0)
      P7 = 2.D0*(A2-A3*P2)+2.D0*A3*A4*S2*P2-4.D0*B1*S2*F
      FS = P3*(P3*P5*P6+P7)
      P8 = 2.D0*S*(B1-A3*A4*P2)
      P9 = A1*P1+A*A1*S*P0*(3.D0-A*A*S2*P0*P0)
      P10 = 4.D0*A3*A4*S*P2*(2.D0-A4*S2)-8.D0*B1*S*F-4.D0*B1*S3*FS
      P11 = -P3*P3*(A1*P1+A*A1*S*P0+4.D0*B1*S3)
      FSS = P3*P3*(P5*P9+P6*P8)+2.D0*P3*P5*P6*P11+P3*P10+P7*P11
      vxuppw91 = FAC*(THRD4*F-(U-THRD4*S3)*FSS-V*FS)
!-------------------------------------------
    else
      exuppw91=0.d0
      vxuppw91=0.d0
    endif
! repeat for down
    rho2=2.d0*dn
    if(rho2.gt.rhocut)then
      fk=(pi32*rho2)**thrd
      s=2.d0*agrdn/(2.d0*fk*rho2)
      u=4.d0*delgrdn/(rho2*rho2*(2.d0*fk)**3)
      v=2.d0*dnlap/(rho2*(2.d0*fk)**2)
!----------  subroutine exchpw91  ----------
      a1=0.19645D0;a2=0.27430D0;a3=0.15084D0;a4=100.d0
      ax=-0.7385588D0;a=7.7956D0;b1=0.004d0
      thrd4=1.33333333333D0
! for Becke exchange, set a3=b1=0
      FAC = AX*RHO2**THRD
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.D0/DSQRT(1.D0+A*A*S2)
      P1 = DLOG(A*S+1.D0/P0)
      P2 = DEXP(-A4*S2)
      P3 = 1.D0/(1.D0+A1*S*P1+B1*S4)
      P4 = 1.D0+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      exdnpw91 = FAC*F
!  ENERGY DONE. NOW THE POTENTIAL:
      P5 = B1*S2-(A2-A3*P2)
      P6 = A1*S*(P1+A*S*P0)
      P7 = 2.D0*(A2-A3*P2)+2.D0*A3*A4*S2*P2-4.D0*B1*S2*F
      FS = P3*(P3*P5*P6+P7)
      P8 = 2.D0*S*(B1-A3*A4*P2)
      P9 = A1*P1+A*A1*S*P0*(3.D0-A*A*S2*P0*P0)
      P10 = 4.D0*A3*A4*S*P2*(2.D0-A4*S2)-8.D0*B1*S*F-4.D0*B1*S3*FS
      P11 = -P3*P3*(A1*P1+A*A1*S*P0+4.D0*B1*S3)
      FSS = P3*P3*(P5*P9+P6*P8)+2.D0*P3*P5*P6*P11+P3*P10+P7*P11
      vxdnpw91 = FAC*(THRD4*F-(U-THRD4*S3)*FSS-V*FS)
!-------------------------------------------
    else
      exdnpw91=0.d0
      vxdnpw91=0.d0
    endif
! construct total density and contribution to ex
    rhot=up+dn
    expw91=(exuppw91*up+exdnpw91*dn)/rhot
!*************************************
!**********  CORRELATION PART  **********
    if(rhot.gt.rhocut) then
      zet=(up-dn)/rhot
      if (zet > 1.0d0-rhocut) zet=1.0d0-rhocut
      if (zet < rhocut-1.0d0) zet=rhocut-1.0d0
      g=((1.d0+zet)**thrd2+(1.d0-zet)**thrd2)/2.d0
      fk=(pi32*rhot)**thrd
      rs=alpha/fk
      sk=sqrt(4.d0*fk/pi)
      twoksg=2.d0*sk*g
      t=agr/(twoksg*rhot)
      uu=delgr/(rhot*rhot*twoksg**3)
      rholap=uplap+dnlap
      vv=rholap/(rhot*twoksg**2)
      ww=(agrup**2-agrdn**2-zet*agr**2)/(rhot*rhot*twoksg**2)
!==========  subroutine corlsd  ==========
      gam=0.5198421D0;fzz=1.709921D0
      thrd4=1.333333333333D0
      F = ((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
!----------  subroutine corlsd  ----------
      A=0.0310907D0;A1=0.21370D0;B1=7.5957D0;B2=3.5876D0;B3=1.6382D0;B4=0.49294D0;P=1.00D0
      P1 = P + 1.D0
      Q0 = -2.D0*A*(1.D0+A1*RS)
      RS12 = DSQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      eu=gg;eurs=ggrs
!-----------------------------------------
!----------  subroutine corlsd  ----------
      A=0.01554535D0;A1=0.20548D0;B1=14.1189D0;B2=6.1977D0;B3=3.3662D0;B4=0.62517D0;P=1.00D0
      P1 = P + 1.D0
      Q0 = -2.D0*A*(1.D0+A1*RS)
      RS12 = DSQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      ep=gg;eprs=ggrs
!-----------------------------------------
!----------  subroutine corlsd  ----------
      A=0.0168869D0;A1=0.11125D0;B1=10.357D0;B2=3.6231D0;B3=0.88026D0;B4=0.49671D0;P=1.00D0
      P1 = P + 1.D0
      Q0 = -2.D0*A*(1.D0+A1*RS)
      RS12 = DSQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      alfm=gg;alfrsm=ggrs
!-----------------------------------------
!  ALFM IS MINUS THE SPIN STIFFNESS ALFC
      ALFC = -ALFM
      Z4 = ZET**4
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
!  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
              -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
!=========================================
!----------  subroutine CORpw91  ----------
      xnu=15.75592D0;cc0=0.004235D0;cx=-0.001667212D0
      alf=0.09D0
      c1=0.002568D0;c2=0.023266D0;c3=7.389D-6;c4=8.723D0;c5=0.472D0;c6=7.389D-2;a4=100.D0
      thrdm=-0.333333333333D0
      BET = XNU*CC0
      DELT = 2.D0*ALF/BET
      G3 = G**3
      G4 = G3*G
      PON = -DELT*EC/(G3*BET)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      Q6 = C1+C2*RS+C3*RS2
      Q7 = 1.D0+C4*RS+C5*RS2+C6*RS3
      CC = -CX + Q6/Q7
      R0 = 0.663436444d0*rs
      R1 = A4*R0*G4
      COEFF = CC-CC0-3.D0*CX/7.D0
      R2 = XNU*COEFF*G3
      R3 = DEXP(-R1*T2)
      H0 = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0+H1
      CCRS = (C2+2.*C3*RS)/Q7 - Q6*(C4+2.*C5*RS+3.*C6*RS2)/Q7**2
      RSTHRD = RS/3.D0
      R4 = RSTHRD*CCRS/COEFF
      GZ = ((1.D0+ZET)**THRDM - (1.D0-ZET)**THRDM)/3.D0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      H0B = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      H0RS = -RSTHRD*H0B*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      H0BT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      H0RST = RSTHRD*T2*H0BT*BEC*ECRS
      H0Z = 3.D0*GZ*H0/G + H0B*(BG*GZ+BEC*ECZET)
      H0T = 2.*BET*G3*Q9/Q8
      H0ZT = 3.D0*GZ*H0T/G+H0BT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      H0TT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      H1RS = R3*R2*T2*(-R4+R1*T2/3.D0)
      FACT4 = 2.D0-R1*T2
      H1RST = R3*R2*T2*(2.D0*R4*(1.D0-R1*T2)-THRD2*R1*T2*FACT4)
      H1Z = GZ*R3*R2*T2*(3.D0-4.D0*R1*T2)/G
      H1T = 2.D0*R3*R2*(1.D0-R1*T2)
      H1ZT = 2.D0*GZ*R3*R2*(3.D0-11.D0*R1*T2+4.D0*R1*R1*T4)/G
      H1TT = 4.D0*R3*R2*R1*T*(-2.D0+R1*T2)
      HRS = H0RS+H1RS
      HRST = H0RST+H1RST
      HT = H0T+H1T
      HTT = H0TT+H1TT
      HZ = H0Z+H1Z
      HZT = H0ZT+H1ZT
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
!------------------------------------------
      ecpw91=ec+h
      vcuppw91=vcup+dvcup
      vcdnpw91=vcdn+dvcdn
    else
      ecpw91=0.0d0
      vcuppw91=0.0d0
      vcdnpw91=0.0d0
    end if
!****************************************
    exu(i)=expw91+ecpw91
    vx(i,1)=vxuppw91+vcuppw91
    vx(i,2)=vxdnpw91+vcdnpw91
  end do
  return
end subroutine vxpot_gga_pw91s

! ===================================================================

subroutine vxpot_gga_pbe(n,rho,abgr,ggabg,ggr,vx,exu)
implicit real*8 (a-h,o-z)
integer,intent(in) :: n
real*8 ,intent(in) :: rho(n),abgr(n),ggabg(n),ggr(n)
real*8 ,intent(out):: vx(n),exu(n)
integer :: i

!$omp do
  do i=1,n
    up=rho(i)*0.5d0
    dn=rho(i)*0.5d0
    agrup=abgr(i)*0.5d0
    agrdn=abgr(i)*0.5d0
    delgrup=ggabg(i)*0.25d0
    delgrdn=ggabg(i)*0.25d0
    uplap=ggr(i)*0.5d0
    dnlap=ggr(i)*0.5d0
    agr=abgr(i)
    delgr=ggabg(i)
    thrd=1.d0/3.d0;thrd2=2.d0*thrd
    pi32=29.608813203268075856503472999628d0
    pi=3.1415926535897932384626433832795d0
    alpha=1.91915829267751300662482032624669d0
!**********  EXCHANGE PART  **********
    rho2=2.d0*up
    if(rho2.gt.rhocut)then
      fk=(pi32*rho2)**thrd
      s=2.d0*agrup/(2.d0*fk*rho2)
      u=4.d0*delgrup/(rho2*rho2*(2.d0*fk)**3)
      v=2.d0*uplap/(rho2*(2.d0*fk)**2)
!----------  subroutine exchpbe  ----------
      thrd4=4.d0/3.d0
      ax=-0.738558766382022405884230032680836d0
      um=0.2195149727645171d0;uk=0.8040d0;ul=um/uk
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct LDA exchange energy density
      exunif = AX*rho2**THRD
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct PBE enhancement factor
      S2 = S*S
      P0=1.d0+ul*S2
      FxPBE = 1d0+uk-uk/P0
      exuppbe = exunif*FxPBE
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first and second derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
!  Fss=d Fs/ds
      Fs=2.d0*uk*ul/(P0*P0)
      Fss=-4.d0*ul*S*Fs/P0
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! calculate potential from [b](24) 
      vxuppbe = exunif*(THRD4*FxPBE-(U-THRD4*S2*s)*FSS-V*FS)
!------------------------------------------
    else
      exuppbe=0.d0
      vxuppbe=0.d0
    endif
! construct total density and contribution to ex
    rhot=up+dn
    exdnpbe=exuppbe
    vxdnpbe=vxuppbe
    expbe=(exuppbe*up+exdnpbe*dn)/rhot
!*************************************
!**********  CORRELATION PART  **********
    if(rhot.gt.rhocut) then
      zet=(up-dn)/rhot
      g=((1.d0+zet)**thrd2+(1.d0-zet)**thrd2)/2.d0
      fk=(pi32*rhot)**thrd
      rs=alpha/fk
      sk=sqrt(4.d0*fk/pi)
      twoksg=2.d0*sk*g
      t=agr/(twoksg*rhot)
      uu=delgr/(rhot*rhot*twoksg**3)
      rholap=uplap+dnlap
      vv=rholap/(rhot*twoksg**2)
      ww=(agrup**2-agrdn**2-zet*agr**2)/(rhot*rhot*twoksg**2)
!==========  subroutine corpbe  ==========
      thrdm=-thrd
      sixthm=thrdm/2.d0
      thrd4=4.d0*thrd
      GAM=0.5198420997897463295344212145565d0
      fzz=8.d0/(9.d0*GAM)
      gamma=0.03109069086965489503494086371273d0
      bet=0.06672455060314922d0;delt=bet/gamma
      eta=1.d-12
      rtrs=dsqrt(rs)
!----------  subroutine gcor2  ----------
      A=0.0310907D0;A1=0.21370D0;B1=7.5957D0;B2=3.5876D0;B3=1.6382D0;B4=0.49294D0
      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
      eu=gg;eurs=ggrs
!----------------------------------------
!----------  subroutine gcor2  ----------
      A=0.01554535D0;A1=0.20548D0;B1=14.1189D0;B2=6.1977D0;B3=3.3662D0;B4=0.62517D0
      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
      ep=gg;eprs=ggrs
!----------------------------------------
!----------  subroutine gcor2  ----------
      A=0.0168869D0;A1=0.11125D0;B1=10.357D0;B2=3.6231D0;B3=0.88026D0;B4=0.49671D0
      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
      alfm=gg;alfrsm=ggrs
!----------------------------------------
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
              -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      G=((1.d0+ZET)**thrd2+(1.d0-ZET)**thrd2)/2.d0
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      H = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3.D0
      GZ=(((1.d0+zet)**2+eta)**sixthm- &
      ((1.d0-zet)**2+eta)**sixthm)/3.d0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      hB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      hBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      hRST = RSTHRD*T2*hBT*BEC*ECRS
      hZ = 3.D0*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2.d0*BET*G3*Q9/Q8
      hZT = 3.D0*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      hTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
!=========================================
      ecpbe=ec+h
      vcuppbe=vcup+dvcup
      vcdnpbe=vcdn+dvcdn
    else
      ecpbe=0.0d0
      vcuppbe=0.0d0
      vcdnpbe=0.0d0
    end if
!****************************************
    exu(i)=expbe+ecpbe
    vx(i)=vxuppbe+vcuppbe
  end do

  return
end subroutine vxpot_gga_pbe


subroutine vxpot_gga_pbes(n,rho,abgr,ggabg,ggr,vx,exu)
implicit real*8 (a-h,o-z)
integer,intent(in) :: n
real*8 ,intent(in) :: rho(n,2),abgr(n,3),ggabg(n,3),ggr(n,3)
real*8 ,intent(out):: vx(n,2),exu(n)
integer :: i

!$omp do
  do i=1,n
    up=rho(i,1)
    dn=rho(i,2)
    agrup=abgr(i,1)
    agrdn=abgr(i,2)
    delgrup=ggabg(i,1)
    delgrdn=ggabg(i,2)
    uplap=ggr(i,1)
    dnlap=ggr(i,2)
    agr=abgr(i,3)
    delgr=ggabg(i,3)
    thrd=1.d0/3.d0;thrd2=2.d0*thrd
    pi32=29.608813203268075856503472999628d0
    pi=3.1415926535897932384626433832795d0
    alpha=1.91915829267751300662482032624669d0
!**********  EXCHANGE PART  **********
    rho2=2.d0*up
    if(rho2.gt.rhocut)then
      fk=(pi32*rho2)**thrd
      s=2.d0*agrup/(2.d0*fk*rho2)
      u=4.d0*delgrup/(rho2*rho2*(2.d0*fk)**3)
      v=2.d0*uplap/(rho2*(2.d0*fk)**2)
!----------  subroutine exchpbe  ----------
      thrd4=4.d0/3.d0
      ax=-0.738558766382022405884230032680836d0
      um=0.2195149727645171d0;uk=0.8040d0;ul=um/uk
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct LDA exchange energy density
      exunif = AX*rho2**THRD
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct PBE enhancement factor
      S2 = S*S
      P0=1.d0+ul*S2
      FxPBE = 1d0+uk-uk/P0
      exuppbe = exunif*FxPBE
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first and second derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
!  Fss=d Fs/ds
      Fs=2.d0*uk*ul/(P0*P0)
      Fss=-4.d0*ul*S*Fs/P0
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! calculate potential from [b](24) 
      vxuppbe = exunif*(THRD4*FxPBE-(U-THRD4*S2*s)*FSS-V*FS)
!------------------------------------------
    else
      exuppbe=0.d0
      vxuppbe=0.d0
    endif
! repeat for down
    rho2=2.d0*dn
    if(rho2.gt.rhocut)then
      fk=(pi32*rho2)**thrd
      s=2.d0*agrdn/(2.d0*fk*rho2)
      u=4.d0*delgrdn/(rho2*rho2*(2.d0*fk)**3)
      v=2.d0*dnlap/(rho2*(2.d0*fk)**2)
!----------  subroutine exchpbe  ----------
      thrd4=4.d0/3.d0
      ax=-0.738558766382022405884230032680836d0
      um=0.2195149727645171d0;uk=0.8040d0;ul=um/uk
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct LDA exchange energy density
      exunif = AX*rho2**THRD
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct PBE enhancement factor
      S2 = S*S
      P0=1.d0+ul*S2
      FxPBE = 1d0+uk-uk/P0
      exdnpbe = exunif*FxPBE
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first and second derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
!  Fss=d Fs/ds
      Fs=2.d0*uk*ul/(P0*P0)
      Fss=-4.d0*ul*S*Fs/P0
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! calculate potential from [b](24) 
      vxdnpbe = exunif*(THRD4*FxPBE-(U-THRD4*S2*s)*FSS-V*FS)
!------------------------------------------
    else
      exdnpbe=0.d0
      vxdnpbe=0.d0
    endif
! construct total density and contribution to ex
    rhot=up+dn
    expbe=(exuppbe*up+exdnpbe*dn)/rhot
!*************************************
!**********  CORRELATION PART  **********
    if(rhot.gt.rhocut) then
      zet=(up-dn)/rhot
      if (zet > 1.0d0-rhocut) zet=1.0d0-rhocut
      if (zet < rhocut-1.0d0) zet=rhocut-1.0d0
      g=((1.d0+zet)**thrd2+(1.d0-zet)**thrd2)/2.d0
      fk=(pi32*rhot)**thrd
      rs=alpha/fk
      sk=sqrt(4.d0*fk/pi)
      twoksg=2.d0*sk*g
      t=agr/(twoksg*rhot)
      uu=delgr/(rhot*rhot*twoksg**3)
      rholap=uplap+dnlap
      vv=rholap/(rhot*twoksg**2)
      ww=(agrup**2-agrdn**2-zet*agr**2)/(rhot*rhot*twoksg**2)
!==========  subroutine corpbe  ==========
      thrdm=-thrd
      sixthm=thrdm/2.d0
      thrd4=4.d0*thrd
      GAM=0.5198420997897463295344212145565d0
      fzz=8.d0/(9.d0*GAM)
      gamma=0.03109069086965489503494086371273d0
      bet=0.06672455060314922d0;delt=bet/gamma
      eta=1.d-12
      rtrs=dsqrt(rs)
!----------  subroutine gcor2  ----------
      A=0.0310907D0;A1=0.21370D0;B1=7.5957D0;B2=3.5876D0;B3=1.6382D0;B4=0.49294D0
      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
      eu=gg;eurs=ggrs
!----------------------------------------
!----------  subroutine gcor2  ----------
      A=0.01554535D0;A1=0.20548D0;B1=14.1189D0;B2=6.1977D0;B3=3.3662D0;B4=0.62517D0
      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
      ep=gg;eprs=ggrs
!----------------------------------------
!----------  subroutine gcor2  ----------
      A=0.0168869D0;A1=0.11125D0;B1=10.357D0;B2=3.6231D0;B3=0.88026D0;B4=0.49671D0
      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
      alfm=gg;alfrsm=ggrs
!----------------------------------------
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
              -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      G=((1.d0+ZET)**thrd2+(1.d0-ZET)**thrd2)/2.d0
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      H = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3.D0
      GZ=(((1.d0+zet)**2+eta)**sixthm- &
      ((1.d0-zet)**2+eta)**sixthm)/3.d0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      hB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      hBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      hRST = RSTHRD*T2*hBT*BEC*ECRS
      hZ = 3.D0*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2.d0*BET*G3*Q9/Q8
      hZT = 3.D0*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      hTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
!=========================================
      ecpbe=ec+h
      vcuppbe=vcup+dvcup
      vcdnpbe=vcdn+dvcdn
    else
      ecpbe=0.0d0
      vcuppbe=0.0d0
      vcdnpbe=0.0d0
    end if
!****************************************
    exu(i)=expbe+ecpbe
    vx(i,1)=vxuppbe+vcuppbe
    vx(i,2)=vxdnpbe+vcdnpbe
  end do

  return
end subroutine vxpot_gga_pbes

end module
