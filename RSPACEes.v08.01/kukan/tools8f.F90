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
! **********  tools8f.F90 12/05/2022-01  **********

module mod_tools
implicit none

contains

subroutine tools_sphglpts(lsphel,msp,p,wt)
implicit none
integer lsphel,msp
integer n,ij,m,i,j
real*8 p(msp,3),wt(msp)
real*8 w(lsphel),x(lsphel)
real*8 pi,x1,x2,x3,wj,hw,rs,rc
  pi=dacos(-1.0d0)

  n=lsphel
  call tools_sphglpts_grule(n,x,w)
  wj=pi/dfloat(n)
  hw=wj/2.0d0
  ij=0
  m=(n+1)/2
  do j=1,2*n
     rc=dcos(dfloat(2*j-1)*hw)
     rs=dsin(dfloat(2*j-1)*hw)
     do i=1,m
        x1=dsqrt(1.d0-x(i)*x(i))* rs
        x2=dsqrt(1.d0-x(i)*x(i))* rc
        x3=x(i)
        ij=ij + 1
        p(ij,1)= x1
        p(ij,2)= x2
        p(ij,3)= x3
        wt(ij)=w(i)*wj
        ij=ij+1
        p(ij,1)= x1
        p(ij,2)= x2
        p(ij,3)=-x3
        wt(ij)=w(i)*wj
     end do
  end do
return
end subroutine tools_sphglpts


subroutine tools_sphglpts_grule(n,x,w)
implicit none
integer n
real*8 w(n),x(n)
real*8 d1,d2pn,d3pn,d4pn,den,dp,dpn,e1,fx,h,pi,pk,pkm1,pkp1,t,t1,u,v,x0
real*8 p
integer i,it,k,m

  pi = 4.d0*datan(1.d0)
  m = (n+1)/2
  e1 = n* (n+1)
  do i = 1,m
     t = (4*i-1)*pi/ (4*n+2)
     x0 = (1.d0- (1.d0-1.d0/n)/ (8.d0*n*n))*dcos(t)
     do it = 1,2
        pkm1 = 1.d0
        pk = x0
        do k = 2,n
           t1 = x0*pk
           pkp1 = t1 - pkm1 - (t1-pkm1)/k + t1
           pkm1 = pk
           pk = pkp1
        end do
        den = 1.d0 - x0*x0
        d1 = n* (pkm1-x0*pk)
        dpn = d1/den
        d2pn = (2.d0*x0*dpn-e1*pk)/den
        d3pn = (4.d0*x0*d2pn+ (2.d0-e1)*dpn)/den
        d4pn = (6.d0*x0*d3pn+ (6.d0-e1)*d2pn)/den
        u = pk/dpn
        v = d2pn/dpn
        h = -u* (1.d0+0.5d0*u* (v+u* (v*v-u*d3pn/ (3.d0*dpn))))
        p = pk + h* (dpn+0.5d0*h* (d2pn+h/3.d0* (d3pn+0.25d0*h*d4pn)))
        dp = dpn + h* (d2pn+0.5d0*h* (d3pn+h*d4pn/3.d0))
        h = h - p/dp
        x0 = x0 + h
     end do
     x(i) = x0
     fx = d1 - h*e1* (pk+0.5d0*h* (dpn+h/3.d0* (d2pn+0.25d0*h* (d3pn+.2d0*h*d4pn))))
     w(i) = 2.d0* (1.d0-x(i)*x(i))/ (fx*fx)
  end do
  if (m+m.gt.n) x(m) = 0.d0
return
end subroutine tools_sphglpts_grule


subroutine tools_genyylm(npoint,lrhomx,point,yylm,ndisp)
implicit none
integer, intent(in)  :: npoint,lrhomx,ndisp
real*8,  intent(in)  :: point(npoint,3)
real*8,  intent(out) :: yylm(npoint,lrhomx)
integer j
real*8  pi,r,x,y,z
real*8  dacos
  pi=dacos(-1.0d0)
  if ((lrhomx/=1).and.(lrhomx/=9).and.(lrhomx/=25)) then
    write(ndisp,*) 'error in genyylm. lrhomx should be (2*l+1)**2, with l in {0,1,2}'
    stop
  end if
  do j=1,npoint
    r=1.0d0
    x=point(j,1)*r
    y=point(j,2)*r
    z=point(j,3)*r
    yylm(j, 1)= 1.0d0/dsqrt(4.0d0*pi)
    if (lrhomx>1) then
      yylm(j, 2)= dsqrt(3.0d0/(4.0d0*pi))*x/r
      yylm(j, 3)= dsqrt(3.0d0/(4.0d0*pi))*y/r
      yylm(j, 4)= dsqrt(3.0d0/(4.0d0*pi))*z/r
      yylm(j, 5)= dsqrt(15.0d0/16.0d0/pi)*(x*x-y*y)/r/r
      yylm(j, 6)=-dsqrt(15.0d0/4.0d0/pi)*z*x/r/r
      yylm(j, 7)= dsqrt(5.0d0/16.0d0/pi)*(3.0d0*z*z-r*r)/r/r
      yylm(j, 8)=-dsqrt(15.0d0/4.0d0/pi)*y*z/r/r
      yylm(j, 9)= dsqrt(15.0d0/4.0d0/pi)*x*y/r/r
    end if
    if (lrhomx>9) then
      yylm(j,10)= dsqrt(7.0d0/16.0d0/pi)*z*(5.0d0*z*z-3.0d0*r*r)/(r*r*r)
      yylm(j,11)=-dsqrt(21.0d0/32.0d0/pi)*x*(5.0d0*z*z-r*r)/(r*r*r)
      yylm(j,12)= dsqrt(105.0d0/16.0d0/pi)*z*(x*x-y*y)/(r*r*r)
      yylm(j,13)=-dsqrt(35.0d0/32.0d0/pi)*x*(x*x-3.0d0*y*y)/(r*r*r)
      yylm(j,14)=-dsqrt(21.0d0/32.0d0/pi)*y*(5.0d0*z*z-r*r)/(r*r*r)
      yylm(j,15)= dsqrt(105.0d0/4.0d0/pi)*x*y*z/(r*r*r)
      yylm(j,16)=-dsqrt(35.0d0/32.0d0/pi)*y*(3.0d0*x*x-y*y)/(r*r*r)
      yylm(j,17)= dsqrt(9.0d0/256.0d0/pi)*(35.0d0*z*z*z*z-30.0d0*z*z*r*r+3.0d0*r*r*r*r)/(r*r*r*r)
      yylm(j,18)= dsqrt(45.0d0/32.0d0/pi)*x*z*(7.0d0*z*z-3.0d0*r*r)/(r*r*r*r)
      yylm(j,19)= dsqrt(45.0d0/64.0d0/pi)*(x*x-y*y)*(7.0d0*z*z-1.0d0*r*r)/(r*r*r*r)
      yylm(j,20)= dsqrt(315.0d0/32.0d0/pi)*z*x*(x*x-3.0d0*y*y)/(r*r*r*r)
      yylm(j,21)= dsqrt(315.0d0/256.0d0/pi)*(x*x*x*x-6.0d0*x*x*y*y+y*y*y*y)/(r*r*r*r)
      yylm(j,22)= dsqrt(45.0d0/32.0d0/pi)*y*z*(7.0d0*z*z-3.0d0*r*r)/(r*r*r*r)
      yylm(j,23)= dsqrt(45.0d0/16.0d0/pi)*x*y*(7.0d0*z*z-1.0d0*r*r)/(r*r*r*r)
      yylm(j,24)= dsqrt(315.0d0/32.0d0/pi)*z*y*(y*y-3.0d0*x*x)/(r*r*r*r)
      yylm(j,25)= dsqrt(315.0d0/16.0d0/pi)*x*y*(x*x-y*y)/(r*r*r*r)
    end if
  end do
return
end subroutine tools_genyylm


subroutine tools_definefdcoef(nf,acfd0,acfd1,acfd2,acfd3,acfd4,acfd5,acfd6,acfd7,acfd8)
implicit none
integer nf
real*8 acfd0,acfd1,acfd2,acfd3,acfd4,acfd5,acfd6,acfd7,acfd8
  if (nf .eq. 1) then
    acfd0=-2.0d0
    acfd1=1.0d0
    acfd2=0.0d0
    acfd3=0.0d0
    acfd4=0.0d0
    acfd5=0.0d0
    acfd6=0.0d0
    acfd7=0.0d0
    acfd8=0.0d0
  end if
  if (nf .eq. 2) then
    acfd0=-5.0d0/2.0d0
    acfd1=4.0d0/3.0d0
    acfd2=-1.0d0/12.0d0
    acfd3=0.0d0
    acfd4=0.0d0
    acfd5=0.0d0
    acfd6=0.0d0
    acfd7=0.0d0
    acfd8=0.0d0
  end if
  if (nf .eq. 3) then
    acfd0=-49.0d0/18.0d0
    acfd1=3.0d0/2.0d0
    acfd2=-3.0d0/20.0d0
    acfd3=1.0d0/90.0d0
    acfd4=0.0d0
    acfd5=0.0d0
    acfd6=0.0d0
    acfd7=0.0d0
    acfd8=0.0d0
  end if
  if (nf .eq. 4) then
    acfd0=-205.0d0/72.0d0
    acfd1=8.0d0/5.0d0
    acfd2=-1.0d0/5.0d0
    acfd3=8.0d0/315.0d0
    acfd4=-1.0d0/560.0d0
    acfd5=0.0d0
    acfd6=0.0d0
    acfd7=0.0d0
    acfd8=0.0d0
  end if
  if (nf .eq. 5) then
    acfd0=-36883.0d0/12600.0d0
    acfd1=5.0d0/3.0d0
    acfd2=-5.0d0/21.0d0
    acfd3=5.0d0/126.0d0
    acfd4=-5.0d0/1008.0d0
    acfd5=1.0d0/3150.0d0
    acfd6=0.0d0
    acfd7=0.0d0
    acfd8=0.0d0
  end if
  if (nf .eq. 6) then
    acfd0=-5369.0d0/1800.0d0
    acfd1=12.0d0/7.0d0
    acfd2=-15.0d0/56.0d0
    acfd3=10.0d0/189.0d0
    acfd4=-1.0d0/112.0d0
    acfd5=2.0d0/1925.0d0
    acfd6=-1.0d0/16632.0d0
    acfd7=0.0d0
    acfd8=0.0d0
  end if
  if (nf .eq. 7) then
    acfd0=-3.0235941043083900d0
    acfd1= 1.7500000000000000d0
    acfd2=-0.2916666666666667d0
    acfd3= 6.4814814814814815d-2
    acfd4=-1.3257575757575758d-2
    acfd5= 2.1212121212121212d-3
    acfd6=-2.2662522662522663d-4
    acfd7= 1.1892869035726179d-5
    acfd8= 0.0d0
  end if
  if (nf .eq. 8) then
    acfd0=-3.0548441043083900d0
    acfd1= 1.7777777777777778d0
    acfd2=-0.3111111111111111d0
    acfd3= 7.5420875420875421d-2
    acfd4=-1.7676767676767677d-2
    acfd5= 3.4809634809634810d-3
    acfd6=-5.1800051800051800d-4
    acfd7= 5.0742907885765029d-5
    acfd8=-2.4281274281274281d-6
  end if
return
end subroutine tools_definefdcoef


subroutine tools_listvecdim( &
 natom,num_spe,nradmx,nxmax,nymax,nzmax,ncpx,ncpy,ncpz,nmesh,nfdg, & ! <
 indspe,nradct,dx,dy,dz,psctoff,radial,                            & ! <
 num_list,num_list_d,                                              & ! >
 npxmax,npymax,npzmax)                                               ! X
implicit none
integer,intent(in) :: natom,num_spe,nradmx,nxmax,nymax,nzmax,ncpx,ncpy,ncpz,nmesh,nfdg
integer,intent(in) :: indspe(natom),nradct(num_spe)
real*8, intent(in) :: dx,dy,dz,psctoff
real*8, intent(in) :: radial(nradmx,num_spe)
integer,intent(inout):: num_list,num_list_d,npxmax,npymax,npzmax
integer:: loop,nftmp,na,npxmx,npymx,npzmx,nlistx,nlisty,nlistz
real*8 :: rcut
logical lauto

  lauto=.false.
  if ((npxmax<0) .or. (npymax<0) .or. (npzmax<0)) lauto=.true.

  do loop=1,2

    if (loop==1) nftmp= nfdg
    if (loop==2) nftmp= 0

    if (lauto) then
      npxmx= 1
      npymx= 1
      npzmx= 1
      do na=1,natom
        rcut= radial(nradct(indspe(na)),indspe(na))*psctoff
        npxmx= max(npxmx, int(rcut/dx)+1+nftmp )
        npymx= max(npymx, int(rcut/dy)+1+nftmp )
        npzmx= max(npzmx, int(rcut/dz)+1+nftmp )
      end do
      if (loop==1) then
        npxmax    = npxmx
        npymax    = npymx
        npzmax    = npzmx
      end if
    else
      npxmx= npxmax
      npymx= npymax
      npzmx= npzmax
    end if

    if (2*nxmax .lt. 2*npxmx) then
      nlistx=2*npxmx
    else
      if (ncpx .gt. 2*npxmx) then
        nlistx=2*npxmx
      else
        nlistx=ncpx
      end if
    end if
    if (2*nymax .lt. 2*npymx) then
      nlisty=2*npymx
    else
      if (ncpy .gt. 2*npymx) then
        nlisty=2*npymx
      else
        nlisty=ncpy
      end if
    end if
    if (2*nzmax .lt. 2*npzmx) then
      nlistz=2*npzmx
    else
      if (ncpz .gt. 2*npzmx) then
        nlistz=2*npzmx
      else
        nlistz=ncpz
      end if
    end if
    if (loop==1) then
      num_list  = nlistx*nlisty*nlistz
    else
      num_list_d= nlistx*nlisty*nlistz*nmesh**3
    end if

  end do ! loop

end subroutine tools_listvecdim


subroutine tools_indspecies(natom,num_spe,numz,mspe,indspe)
use mod_mpi
use mod_stopp
implicit none
integer, intent(in)::natom
integer, intent(in)::num_spe
integer, intent(in)::numz(natom)
integer, intent(out)::mspe(num_spe)
integer, intent(out)::indspe(natom)
integer :: na,ispe,nspe

  nspe= 0
  mspe(:)=0
  do na= 1,natom
    ispe= 0
    do while (ispe<nspe)
      ispe= ispe + 1
      if (mspe(ispe)==numz(na)) then
        indspe(na)= ispe
        ispe= nspe+1
      end if
    end do
    if (ispe==nspe) then
      nspe= nspe + 1
      mspe(nspe)= numz(na)
      indspe(na)= nspe
    end if
  end do

end subroutine tools_indspecies


subroutine tools_correctpotential(ncpx_d,ncpy_d,ncpz_d,xmax,ymax,zmax,ddx,ddy,ddz,vh_dense)
use mod_mpi
implicit none
integer, intent(in)::ncpx_d,ncpy_d,ncpz_d
real*8, intent(in)::xmax,ymax,zmax,ddx,ddy,ddz
real*8, intent(inout)::vh_dense(ncpx_d,ncpy_d,ncpz_d)
integer ix,iy,iz
real*8 vh_ave,omega,omegain,tmp

  omega=8.0d0*xmax*ymax*zmax
  omegain=1.0d0/omega
  vh_ave=0.0d0
  do iz=1,ncpz_d
  do iy=1,ncpy_d
  do ix=1,ncpx_d
    vh_ave=vh_ave+vh_dense(ix,iy,iz)
  end do
  end do
  end do
  call mpi_allreduce(vh_ave,tmp,1,mpi_double_precision,mpi_sum,mpicom_space,mpij)
  vh_ave=tmp*ddx*ddy*ddz*omegain
  vh_dense=vh_dense-vh_ave

return
end subroutine tools_correctpotential


subroutine tools_jelliumpotential( &
 key_jel_calc,jelcalc,ncpz,nzmax,xmax,ymax,zmax,chrjel,endjel,strjel, & ! <
 vjell)                                                                 ! >
use mod_mpi,       only:myrz
implicit none
integer, intent(in) :: key_jel_calc
integer, intent(in) :: jelcalc,ncpz,nzmax
real*8,  intent(in) :: xmax,ymax,zmax,chrjel,endjel,strjel
real*8,  intent(out):: vjell(ncpz)
integer :: kz,iz
real*8  :: pi,omega,vkz,vk2,dz,z

  pi= 4.0d0*datan(1.0d0)

  vjell(:)=0.0d0
  if (jelcalc==key_jel_calc) then
    omega=xmax*ymax*zmax*8.0d0
    dz= zmax/dble(nzmax)
    do kz=-nzmax,nzmax
      if (kz**2 .ne. 0) then
        vkz=pi/zmax*kz
        vk2=vkz*vkz
        do iz=1,ncpz
          z=dble(myrz*ncpz+iz)*dz-zmax-0.5d0*dz
          vjell(iz)= vjell(iz) &
           -4.0d0*pi/omega*chrjel/(strjel-endjel)/vk2*(dsin(vkz*(strjel-z))-dsin(vkz*(endjel-z)))/vkz
        end do
      end if
    end do
  end if

end subroutine tools_jelliumpotential


subroutine tools_shiftatoms( &
 nperi,natom,xmax,ymax,zmax,dx,dy,dz, & ! <
 atx,aty,atz)                           ! X
use mod_stopp
implicit none
integer,intent(in)::nperi,natom
real*8, intent(in)::xmax,ymax,zmax,dx,dy,dz
real*8, intent(inout)::atx(natom),aty(natom),atz(natom)
integer na

  do na=1,natom

    if (nperi==0) then

      if ( (atx(na) .gt. xmax+0.5d0*dx).or.(atx(na) .lt.-xmax+0.5d0*dx) &
       .or.(aty(na) .gt. ymax+0.5d0*dy).or.(aty(na) .lt.-ymax+0.5d0*dy) &
       .or.(atz(na) .gt. zmax+0.5d0*dz).or.(atz(na) .lt.-zmax+0.5d0*dz) ) then
        call stopp ('tools_shiftatoms: atom position out of bounds')
      end if

    else

      do while (atx(na) .gt. xmax+0.5d0*dx)
        atx(na)=atx(na)-2.0d0*xmax
      end do
      do while (atx(na) .lt.-xmax+0.5d0*dx)
        atx(na)=atx(na)+2.0d0*xmax
      end do
      do while (aty(na) .gt. ymax+0.5d0*dy)
        aty(na)=aty(na)-2.0d0*ymax
      end do
      do while (aty(na) .lt.-ymax+0.5d0*dy)
        aty(na)=aty(na)+2.0d0*ymax
      end do
      do while (atz(na) .gt. zmax+0.5d0*dz)
        atz(na)=atz(na)-2.0d0*zmax
      end do
      do while (atz(na) .lt.-zmax+0.5d0*dz)
        atz(na)=atz(na)+2.0d0*zmax
      end do

    end if

  end do
return
end subroutine tools_shiftatoms


subroutine tools_calenespr(chdirinp,natom,xmax,ymax,zmax,atx,aty,atz,sconst,enespr)
use mod_mpi
implicit none
character,intent(in)               :: chdirinp*200
integer,intent(in)                 :: natom
real*8,intent(in)                  :: xmax,ymax,zmax,sconst
real*8,dimension(natom),intent(in) :: atx,aty,atz
real*8,intent(out)                 :: enespr
real*8,allocatable,dimension(:)    :: atxm,atym,atzm,atxp,atyp,atzp
integer :: i,ierr
character :: chara*1,fname*200

  enespr=0.0d0
  if (myrank_glbl==0) then
    allocate(atxm(natom),atym(natom),atzm(natom),atxp(natom),atyp(natom),atzp(natom))
    fname='atomm.xyz'
    if (len_trim(chdirinp) > 0) fname=trim(chdirinp)//'/'//fname
    open(25,file=fname,iostat=ierr)
      read(25,*,iostat=ierr) chara
      do i=1,natom
        read(25,*,iostat=ierr) atxm(i),atym(i),atzm(i)
        if (atx(i)-atxm(i) < -xmax) atxm(i)=atxm(i)-2.0d0*xmax
        if (atx(i)-atxm(i) >  xmax) atxm(i)=atxm(i)+2.0d0*xmax
        if (aty(i)-atym(i) < -ymax) atym(i)=atym(i)-2.0d0*ymax
        if (aty(i)-atym(i) >  ymax) atym(i)=atym(i)+2.0d0*ymax
        if (atz(i)-atzm(i) < -zmax) atzm(i)=atzm(i)-2.0d0*zmax
        if (atz(i)-atzm(i) >  zmax) atzm(i)=atzm(i)+2.0d0*zmax
      end do
    close(25)
    fname='atomp.xyz'
    if (len_trim(chdirinp) > 0) fname=trim(chdirinp)//'/'//fname
    open(25,file=fname,iostat=ierr)
      read(25,*,iostat=ierr) chara
      do i=1,natom
        read(25,*,iostat=ierr) atxp(i),atyp(i),atzp(i)
        if (atx(i)-atxp(i) < -xmax) atxp(i)=atxp(i)-2.0d0*xmax
        if (atx(i)-atxp(i) >  xmax) atxp(i)=atxp(i)+2.0d0*xmax
        if (aty(i)-atyp(i) < -ymax) atyp(i)=atyp(i)-2.0d0*ymax
        if (aty(i)-atyp(i) >  ymax) atyp(i)=atyp(i)+2.0d0*ymax
        if (atz(i)-atzp(i) < -zmax) atzp(i)=atzp(i)-2.0d0*zmax
        if (atz(i)-atzp(i) >  zmax) atzp(i)=atzp(i)+2.0d0*zmax
      end do
    close(25)
    do i=1,natom
      enespr=enespr+0.5d0*sconst*((atxp(i)-atx(i))**2+(atyp(i)-aty(i))**2+(atzp(i)-atz(i))**2 &
                                   +(atxm(i)-atx(i))**2+(atym(i)-aty(i))**2+(atzm(i)-atz(i))**2)
    end do
    deallocate(atxm,atym,atzm,atxp,atyp,atzp)
  end if

end subroutine tools_calenespr


subroutine tools_moveatoms(chdirinp,chdirout,natom,ngdiis,ndisp,watom,nmdx,nmdy,nmdz,xmax,ymax,zmax,tmstep,fcut,sconst,eps, &
                           atx,aty,atz,fatx,faty,fatz)
use mod_mpi
implicit none
character,intent(in)                  :: chdirinp*200,chdirout*200
integer,intent(in)                    :: natom,ngdiis,ndisp
real*8,intent(in)                     :: xmax,ymax,zmax,tmstep,fcut,sconst,eps

real*8,dimension(natom),intent(inout) :: atx,aty,atz
real*8,dimension(natom),intent(in) :: watom
real*8,dimension(natom),intent(inout) :: fatx,faty,fatz
integer,dimension(natom),intent(in)   :: nmdx,nmdy,nmdz

integer,allocatable,dimension(:)  :: ipiv
real*8,allocatable,dimension(:,:) :: dpos,dfor
real*8,allocatable,dimension(:)   :: psta,dste
real*8,allocatable,dimension(:,:) :: dmat
real*8,allocatable,dimension(:)   :: dvec
real*8,allocatable,dimension(:)   :: atxm,atym,atzm,atxp,atyp,atzp
integer :: mdstep
integer :: i,j,k,na,ndim,ierr
real*8 :: sum,cosphi,fphi,pi
real*8 :: taux,tauy,tauz,fsx,fsy,fsz,fsx1,fsy1,fsz1,fatx1,faty1,fatz1
character :: chara*1,fname*200
  pi=dacos(-1.0d0)

  if (myrank_glbl==0) then
!   ==========  for NEB  ==========
    if (sconst > eps) then
      allocate(atxm(natom),atym(natom),atzm(natom),atxp(natom),atyp(natom),atzp(natom))
      fname='atomm.xyz'
      if (len_trim(chdirinp) > 0) fname=trim(chdirinp)//'/'//fname
      open(25,file=fname,iostat=ierr)
        read(25,*,iostat=ierr) chara
        do i=1,natom
          read(25,*,iostat=ierr) atxm(i),atym(i),atzm(i)
          if (atx(i)-atxm(i) < -xmax) atxm(i)=atxm(i)-2.0d0*xmax
          if (atx(i)-atxm(i) >  xmax) atxm(i)=atxm(i)+2.0d0*xmax
          if (aty(i)-atym(i) < -ymax) atym(i)=atym(i)-2.0d0*ymax
          if (aty(i)-atym(i) >  ymax) atym(i)=atym(i)+2.0d0*ymax
          if (atz(i)-atzm(i) < -zmax) atzm(i)=atzm(i)-2.0d0*zmax
          if (atz(i)-atzm(i) >  zmax) atzm(i)=atzm(i)+2.0d0*zmax
        end do
      close(25)
      fname='atomp.xyz'
      if (len_trim(chdirinp) > 0) fname=trim(chdirinp)//'/'//fname
      open(25,file=fname,iostat=ierr)
        read(25,*,iostat=ierr) chara
        do i=1,natom
          read(25,*,iostat=ierr) atxp(i),atyp(i),atzp(i)
          if (atx(i)-atxp(i) < -xmax) atxp(i)=atxp(i)-2.0d0*xmax
          if (atx(i)-atxp(i) >  xmax) atxp(i)=atxp(i)+2.0d0*xmax
          if (aty(i)-atyp(i) < -ymax) atyp(i)=atyp(i)-2.0d0*ymax
          if (aty(i)-atyp(i) >  ymax) atyp(i)=atyp(i)+2.0d0*ymax
          if (atz(i)-atzp(i) < -zmax) atzp(i)=atzp(i)-2.0d0*zmax
          if (atz(i)-atzp(i) >  zmax) atzp(i)=atzp(i)+2.0d0*zmax
        end do
      close(25)
      do i=1,natom
        taux=atxp(i)-atxm(i)
        tauy=atyp(i)-atym(i)
        tauz=atzp(i)-atzm(i)
        sum=dsqrt(taux*taux+tauy*tauy+tauz*tauz)
        if (sum > 1.0d-16) then
          taux=taux/sum
          tauy=tauy/sum
          tauz=tauz/sum
          fatx1=fatx(i)-(fatx(i)*taux+faty(i)*tauy+fatz(i)*tauz)*taux
          faty1=faty(i)-(fatx(i)*taux+faty(i)*tauy+fatz(i)*tauz)*tauy
          fatz1=fatz(i)-(fatx(i)*taux+faty(i)*tauy+fatz(i)*tauz)*tauz
          fsx=sconst*(atxp(i)-2.0d0*atx(i)+atxm(i))
          fsy=sconst*(atyp(i)-2.0d0*aty(i)+atym(i))
          fsz=sconst*(atzp(i)-2.0d0*atz(i)+atzm(i))
          fsx1=(fsx*taux+fsy*tauy+fsz*tauz)*taux
          fsy1=(fsx*taux+fsy*tauy+fsz*tauz)*tauy
          fsz1=(fsx*taux+fsy*tauy+fsz*tauz)*tauz
          cosphi=((atxp(i)-atx(i))*(atx(i)-atxm(i))+(atyp(i)-aty(i))*(aty(i)-atym(i))+(atzp(i)-atz(i))*(atz(i)-atzm(i))) &
              /(dsqrt((atxp(i)-atx(i))**2+(atyp(i)-aty(i))**2+(atzp(i)-atz(i))**2) &
               *dsqrt((atx(i)-atxm(i))**2+(aty(i)-atym(i))**2+(atz(i)-atzm(i))**2))
          if ((cosphi >= 0.0d0) .and. (cosphi <= 1.0d0)) then
            fphi=0.5d0*(1.0d0+dcos(pi*cosphi))
          else
            fphi=1.0d0
          end if
          fatx(i)=(fatx1+fsx1)+fphi*(fsx-fsx1)
          faty(i)=(faty1+fsy1)+fphi*(fsy-fsy1)
          fatz(i)=(fatz1+fsz1)+fphi*(fsz-fsz1)
        end if
      end do
      deallocate(atxm,atym,atzm,atxp,atyp,atzp)
    end if
!   ===============================

!   ==========  for GDIIS  ==========
    allocate(dpos(3*natom,ngdiis))
    allocate(dfor(3*natom,ngdiis))
    allocate(dste(3*natom))
    allocate(psta(3*natom))
    fname='gdiis_hist.txt'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(25,file=fname,status='old',iostat=ierr)
    if (ierr==0) then
       read(25,*,iostat=ierr) mdstep
       if (mdstep<=ngdiis) then
         if (ierr==0) read(25,*,iostat=ierr) ((dpos(i,j),i=1,3*natom),j=1,mdstep)
         if (ierr==0) read(25,*,iostat=ierr) ((dfor(i,j),i=1,3*natom),j=1,mdstep)
       else
         ierr=1
       end if
    end if
    close(25)
    if (ierr/=0) then
      mdstep=0
      dpos=0.0d0
      dfor=0.0d0
    end if
    if (mdstep>=1) then
      do i=1,natom
        if (atx(i)-dpos(i        ,mdstep) < -xmax) then
          do j=1,mdstep
            dpos(i        ,j)=dpos(i        ,j)-2.0d0*xmax
          end do
        end if
        if (atx(i)-dpos(i        ,mdstep) >  xmax) then
          do j=1,mdstep
            dpos(i        ,j)=dpos(i        ,j)+2.0d0*xmax
          end do
        end if
        if (aty(i)-dpos(i+  natom,mdstep) < -ymax) then
          do j=1,mdstep
            dpos(i+  natom,j)=dpos(i+  natom,j)-2.0d0*ymax
          end do
        end if
        if (aty(i)-dpos(i+  natom,mdstep) >  ymax) then
          do j=1,mdstep
            dpos(i+  natom,j)=dpos(i+  natom,j)+2.0d0*ymax
          end do
        end if
        if (atz(i)-dpos(i+2*natom,mdstep) < -zmax) then
          do j=1,mdstep
            dpos(i+2*natom,j)=dpos(i+2*natom,j)-2.0d0*zmax
          end do
        end if
        if (atz(i)-dpos(i+2*natom,mdstep) >  zmax) then
          do j=1,mdstep
            dpos(i+2*natom,j)=dpos(i+2*natom,j)+2.0d0*zmax
          end do
        end if
      end do
    end if
 
    mdstep=mdstep+1
    ndim=mdstep+1
 
    allocate(ipiv(ndim),dmat(ndim,ndim),dvec(ndim))
 
    do i=1,natom
      dfor(i        ,mdstep)=fatx(i)*nmdx(i)
      dfor(i+natom  ,mdstep)=faty(i)*nmdy(i)
      dfor(i+2*natom,mdstep)=fatz(i)*nmdz(i)
      dpos(i        ,mdstep)=atx(i)
      dpos(i+natom  ,mdstep)=aty(i)
      dpos(i+2*natom,mdstep)=atz(i)
    end do
 
    if (mdstep>1) then
      dmat=0.0d0
      do i=1,ndim-1
        do j=1,ndim-1
          do k=1,3*natom
            dmat(i,j)=dmat(i,j)+dfor(k,i)*dfor(k,j)
          end do
        end do
      end do
      dmat(ndim,:)=1.0d0
      dmat(:,ndim)=1.0d0
      dmat(ndim,ndim)=0.0d0
      dvec=0.0d0
      dvec(ndim)=1.0d0
 
      call dgesv(ndim,1,dmat,ndim,ipiv,dvec,ndim,ierr)
 
      if (ierr==0) then
        dste=0.0d0
        psta=0.0d0
        do j=1,ndim-1
          do i=1,3*natom
            dste(i)=dste(i)+dvec(j)*dfor(i,j)
            psta(i)=psta(i)+dvec(j)*dpos(i,j)
          end do
        end do
        do i=1,natom
          atx(i)=psta(i        )+0.5d0*dste(i        )/watom(i)*tmstep**2*nmdx(i)
          aty(i)=psta(i+natom  )+0.5d0*dste(i+natom  )/watom(i)*tmstep**2*nmdy(i)
          atz(i)=psta(i+2*natom)+0.5d0*dste(i+2*natom)/watom(i)*tmstep**2*nmdz(i)
        end do
      else
        write(ndisp,*) '=====  GDIIS history is cleared.  ====='
        write(12,*) '=====  GDIIS history is cleared.  =====.'
        dfor(:,1)=dfor(:,mdstep)
        dpos(:,1)=dpos(:,mdstep)
        mdstep=1
      end if
 
    end if
 
    fname='gdiis_hist.txt'
    if (len_trim(chdirout) > 0) fname=trim(chdirout)//'/'//fname
    open(25,file=fname)
    if (mdstep==ngdiis) then
      write(25,*) mdstep-1
      write(25,*) ((dpos(i,j),i=1,3*natom),j=2,mdstep)
      write(25,*) ((dfor(i,j),i=1,3*natom),j=2,mdstep)
    else
      write(25,*) mdstep
      write(25,*) ((dpos(i,j),i=1,3*natom),j=1,mdstep)
      write(25,*) ((dfor(i,j),i=1,3*natom),j=1,mdstep)
    end if
    close(25)
    deallocate(ipiv,dmat,dvec)
    deallocate(dpos,dfor,dste,psta)
!   =================================

!   ==========  for simple str. opt.  ==========
    if (mdstep<=1) then
      do i=1,natom
        atx(i)=atx(i)+0.5d0*fatx(i)/watom(i)*tmstep**2*nmdx(i)
        aty(i)=aty(i)+0.5d0*faty(i)/watom(i)*tmstep**2*nmdy(i)
        atz(i)=atz(i)+0.5d0*fatz(i)/watom(i)*tmstep**2*nmdz(i)
      end do
    end if
!   ============================================

  end if

  call mpi_bcast(atx,natom,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(aty,natom,mpi_double_precision,0,mpicom_space,mpij)
  call mpi_bcast(atz,natom,mpi_double_precision,0,mpicom_space,mpij)

end subroutine tools_moveatoms


  subroutine tools_countelectron(natom,num_spe,chrgd,indspe,cp, tnumele)
  implicit none
  integer, intent(in) ::natom,num_spe
  real*8,  intent(in) ::chrgd
  integer, intent(in) ::indspe(natom)
  real*8,  intent(in) ::cp(8,num_spe)
  real*8,  intent(out)::tnumele
  integer na
  tnumele=chrgd
  do na=1,natom
     tnumele=tnumele+cp(1,indspe(na))
  end do
return
end subroutine tools_countelectron


!     **********  symmetric7g.f90 04/08/2008-01  **********
!     07/14/2009 module natpri_parameters was renamed to mod_keys
!     04/08/2008 temporary array in symmetric_02 were changed to rho_tmp to vre2
!     03/16/2008 routine for atomic coordinates is added.
!     03/15/2008 compile options for SX series are added.

!     This file contains
!      + symmetric_01,02,04
!     (symmetric_03 is not used.)

subroutine tools_symmetrizeatom(natom,nsym,key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp, &
                                xmax,ymax,zmax,atx,aty,atz)
implicit none
integer, intent(in)::natom,nsym
integer, intent(in)::key_sym_bcc,key_sym_fcc,key_sym_dia,key_sym_hcp
real*8, intent(in)::xmax,ymax,zmax
real*8, intent(inout)::atx(natom),aty(natom),atz(natom)

  if (nsym .eq. key_sym_bcc) then
    if (natom .eq. 2) then
      atx(1)=0.0d0
      aty(1)=0.0d0
      atz(1)=0.0d0
      atx(2)=-xmax
      aty(2)=-ymax
      atz(2)=-zmax
    else
      atx( 1)=-xmax
      aty( 1)=-ymax
      atz( 1)=-zmax
      atx( 2)=-xmax*0.5d0
      aty( 2)=-ymax*0.5d0
      atz( 2)=-zmax*0.5d0
      atx( 3)=0.0d0
      aty( 3)=-ymax
      atz( 3)=-zmax
      atx( 4)= xmax*0.5d0
      aty( 4)=-ymax*0.5d0
      atz( 4)=-zmax*0.5d0
      atx( 5)=-xmax
      aty( 5)=0.0d0
      atz( 5)=-zmax
      atx( 6)=-xmax*0.5d0
      aty( 6)= ymax*0.5d0
      atz( 6)=-zmax*0.5d0
      atx( 7)=0.0d0
      aty( 7)=0.0d0
      atz( 7)=-zmax
      atx( 8)= xmax*0.5d0
      aty( 8)= ymax*0.5d0
      atz( 8)=-zmax*0.5d0
      atx( 9)=-xmax
      aty( 9)=-ymax
      atz( 9)=0.0d0
      atx(10)=-xmax*0.5d0
      aty(10)=-ymax*0.5d0
      atz(10)= zmax*0.5d0
      atx(11)=0.0d0
      aty(11)=-ymax
      atz(11)=0.0d0
      atx(12)= xmax*0.5d0
      aty(12)=-ymax*0.5d0
      atz(12)= zmax*0.5d0
      atx(13)=-xmax
      aty(13)=0.0d0
      atz(13)=0.0d0
      atx(14)=-xmax*0.5d0
      aty(14)= ymax*0.5d0
      atz(14)= zmax*0.5d0
      atx(15)=0.0d0
      aty(15)=0.0d0
      atz(15)=0.0d0
      atx(16)= xmax*0.5d0
      aty(16)= ymax*0.5d0
      atz(16)= zmax*0.5d0
    end if
  end if
  if (nsym .eq. key_sym_fcc) then
    atx(1)=-xmax
    aty(1)=-ymax
    atz(1)=-zmax
    atx(2)=0.0d0
    aty(2)=0.0d0
    atz(2)=-zmax
    atx(3)=-xmax
    aty(3)=0.0d0
    atz(3)=0.0d0
    atx(4)=0.0d0
    aty(4)=-ymax
    atz(4)=0.0d0
  end if
  if (nsym .eq. key_sym_dia) then
    atx(1)=-xmax
    aty(1)=-ymax
    atz(1)=-zmax
    atx(2)=0.0d0
    aty(2)=0.0d0
    atz(2)=-zmax
    atx(3)=-xmax*0.5d0
    aty(3)=-ymax*0.5d0
    atz(3)=-zmax*0.5d0
    atx(4)= xmax*0.5d0
    aty(4)= ymax*0.5d0
    atz(4)=-zmax*0.5d0
    atx(5)=0.0d0
    aty(5)=-ymax
    atz(5)=0.0d0
    atx(6)=-xmax
    aty(6)=0.0d0
    atz(6)=0.0d0
    atx(7)=-xmax*0.5d0
    aty(7)= ymax*0.5d0
    atz(7)= zmax*0.5d0
    atx(8)= xmax*0.5d0
    aty(8)=-ymax*0.5d0
    atz(8)= zmax*0.5d0
  end if
  if (nsym .eq. key_sym_hcp) then
    atx(1)=-xmax
    aty(1)=0.0d0
    atz(1)=-zmax
    atx(2)=0.0d0
    aty(2)=0.0d0
    atz(2)=-zmax
    atx(3)=-xmax*0.5d0
    aty(3)=-ymax
    atz(3)=-zmax
    atx(4)= xmax*0.5d0
    aty(4)=-ymax
    atz(4)=-zmax
    atx(5)=-xmax
    aty(5)=-ymax*2.0d0/3.0d0
    atz(5)=0.0d0
    atx(6)=0.0d0
    aty(6)=-ymax*2.0d0/3.0d0
    atz(6)=0.0d0
    atx(7)=-xmax*0.5d0
    aty(7)=ymax/3.0d0
    atz(7)=0.0d0
    atx(8)= xmax*0.5d0
    aty(8)=ymax/3.0d0
    atz(8)=0.0d0
  end if

return
end subroutine tools_symmetrizeatom


subroutine tools_bisec(z,r,eps,a,b,c,ierr)
implicit none
real*8, intent(in):: z,r,eps
real*8, intent(out):: a,b,c
real*8 a0,a2,fr,f0,f1,f2
integer, intent(out):: ierr
integer l
  a0=dacos(-1.0d0)*0.5d0/r+eps
  a2=2.0d0*dacos(-1.0d0)*0.5d0/r-eps
  fr=-r/2.0d0
  f0=dtan(a0*r)/a0-fr
  f2=dtan(a2*r)/a2-fr
  l=0
5 l=l+1
  a=(a0+a2)*0.5d0
  f0=dtan(a0*r)/a0-fr
  f1=dtan(a*r)/a-fr
  if (f0*f1 .lt. 0.0d0) then
    a2=a
    else
    a0=a
  end if
  if (l .gt. 1000000) then
    ierr=l
    write(6,*) 'ERROR OCCURS in ppbisec'
    return
  end if
  if (a2-a0 .gt. eps) goto 5
  b=-z/r/r/a/dsin(a*r)
  c=-z/r-b*dcos(a*r)
  return
 end subroutine tools_bisec


subroutine tools_scffile( &
 chdir,scfctr) ! X
! Exit s.c.f. loop and write output files,
! if a file 'scfctr' exists and contains a 'T'.
use mod_mpi
implicit none
character, intent(in)  :: chdir*200
logical, intent(inout) :: scfctr(3)
logical :: scfctr1,scfctr3
integer :: j,fileok
character :: fname*200

  if (myrank_glbl==0) then

    scfctr1= scfctr(1)
    scfctr3= scfctr(3)

    fname='scfctr'
    if (len_trim(chdir) > 0) fname=trim(chdir)//'/'//fname
    inquire(file=fname,exist=scfctr(1))
    if (scfctr(1)) then
      open(13,file=fname,form='formatted',action='read')
      fileok= 0
      do j=1,3
        if (fileok==0) read(13,fmt=*,iostat=fileok) scfctr(j)
        if (fileok/=0) scfctr(j)=.false.
      end do
      close(13)
      open(13,file=fname,form='formatted',status='replace')
      write(13,fmt='(l1,2x,A)') .false.,'write files and stop'
      write(13,fmt='(l1,2x,A)') .false.,'write files'
      write(13,fmt='(l1,2x,A)') .false.,'write energy and forces'
      close(13)
    else
      scfctr(:)= .false.
    end if

    scfctr(1)= scfctr(1) .or. scfctr1
    scfctr(3)= scfctr(1) .or. scfctr(3) .or. scfctr3

  end if ! (myrank_glbl==0)
  call mpi_bcast(scfctr,3,mpi_logical,0,mpi_comm_world,mpij)

end subroutine tools_scffile


subroutine tools_nlpcomsplit (splfre,natom,key_natpri_out,natpri)
use mod_mpi
implicit none
logical, intent(in) :: splfre
integer, intent(in) :: natom,key_natpri_out
integer, intent(in) :: natpri(natom)
integer :: newcomm
integer :: na
  do na=1,natom
    if (splfre) then
      if (natpri(na)/=key_natpri_out) then
        call mpi_comm_split(mpicom_space,myr_kpt,myr_space,newcomm,mpij)
      else
        call mpi_comm_split(mpicom_space,mpi_undefined,myr_space,newcomm,mpij)
      end if
      mpicom_atom(na)=newcomm
    else
      if (natpri(na)/=key_natpri_out) call mpi_comm_free(mpicom_atom(na),mpij)
    end if
  end do
end subroutine tools_nlpcomsplit


subroutine tools_nlpcomsplit_nnn (splfre,natom,natx,naty,natz,ncpx,ncpy,ncpz,nxmax,nymax,nzmax,nx,ny,nz)
use mod_mpi
implicit none
logical, intent(in) :: splfre
integer, intent(in) :: natom,ncpx,ncpy,ncpz,nxmax,nymax,nzmax,nx,ny,nz
integer, intent(in) :: natx(natom),naty(natom),natz(natom)
integer :: na,i,ix,iy,iz,jx,jy,jz,kx,ky,kz,lx,ly,lz
  if (splfre) then
    allocate(ntilecomm(0:7))
    lx=nx/2
    ly=ny/2
    lz=nz/2
    do iz=0,1
    do iy=0,1
    do ix=0,1
      jx=mod(myrx/lx+ix,2)
      jy=mod(myry/ly+iy,2)
      jz=mod(myrz/lz+iz,2)
      kx=myrx/lx-jx
      ky=myry/ly-jy
      kz=myrz/lz-jz
      if (kx<0) kx=kx+nprocx/lx
      if (ky<0) ky=ky+nprocy/ly
      if (kz<0) kz=kz+nprocz/lz
      i=kz*nprocx*nprocy+ky*nprocx+kx
      call mpi_comm_split(mpicom_space,(myr_kpt+1)*(i+1),myr_space,ntilecomm(iz*4+iy*2+ix),mpij)
    end do
    end do
    end do
    do na=1,natom
      lx=4/nx
      ly=4/ny
      lz=4/nz
      ix= (natx(na)+nxmax-1-ncpx/lx)
      iy= (naty(na)+nymax-1-ncpy/ly)
      iz= (natz(na)+nzmax-1-ncpz/lz)
      do while (ix>2*nxmax)
        ix=ix-2*nxmax
      end do
      do while (ix<1)
        ix=ix+2*nxmax
      end do
      do while (iy>2*nymax)
        iy=iy-2*nymax
      end do
      do while (iy<1)
        iy=iy+2*nymax
      end do
      do while (iz>2*nzmax)
        iz=iz-2*nzmax
      end do
      do while (iz<1)
        iz=iz+2*nzmax
      end do
      ix=(ix/(nx/2*ncpx))
      iy=(iy/(ny/2*ncpy))
      iz=(iz/(nz/2*ncpz))
      kx=mod(ix,2)
      ky=mod(iy,2)
      kz=mod(iz,2)
      mpicom_atom(na)=ntilecomm(kz*4+ky*2+kx)
    end do
  else
    do i=0,7
      call mpi_comm_free(ntilecomm(i),mpij)
    end do
    deallocate(ntilecomm)
  end if
end subroutine tools_nlpcomsplit_nnn


subroutine tools_potbroadcast( &
 nso,nums,ncol,nspv,ncpx,ncpy,ncpz,natom,nprjmx, & ! <
 veff,dij,dijsoc)                                  ! X
use mod_mpi
implicit none
integer, intent(in)    :: nso,nums,ncol,nspv,ncpx,ncpy,ncpz,natom,nprjmx
real*8,  intent(inout) :: veff(ncpx,ncpy,ncpz,nspv)
real*8,  intent(inout) :: dij(nprjmx,nprjmx,nums*ncol,natom)
real*8,  intent(inout) :: dijsoc(nprjmx*nso-nso+1,nprjmx*nso-nso+1,3*nso-nso+1,natom*nso-nso+1)

  call mpi_bcast(veff,ncpx*ncpy*ncpz*nspv,mpi_double_precision,0,mpicom_kpt,mpij)
  call mpi_bcast(dij,nprjmx*nprjmx*nums*ncol*natom,mpi_double_precision,0,mpicom_kpt,mpij)
  if (nso==1) call mpi_bcast(dijsoc,nprjmx*nprjmx*3*natom,mpi_double_precision,0,mpicom_kpt,mpij)

end subroutine tools_potbroadcast


subroutine tools_rodr_sum_nlp(natom,key_natpri_out,natpri,latom)
use mod_mpi
implicit none
integer,intent(in)::natom,key_natpri_out
integer,intent(in)::natpri(natom)
integer,intent(out)::latom(natom)
integer,allocatable::ifl(:)
integer i,j,k,iend,lcut,lcut0,nstep

  nstep=nprocs/27
  if (nstep < 1) nstep=1
  if (natom >= nstep) then
    allocate(ifl(natom))
    ifl=0
    latom=0
    k=1
    do while (k <= natom)
      j=1
      do while (j <= natom)
        if (ifl(j) == 0) then
          i=k-1
          iend=k-nstep
          if (iend .lt. 1) iend=1
          lcut=0
          do while (i >= iend)
            if (natpri(latom(i))/=key_natpri_out) lcut=1
            i=i-1
          end do
          if (natpri(j)==key_natpri_out) lcut=0
          lcut0=0
          call mpi_allreduce(lcut,lcut0,1,mpi_integer,mpi_sum,mpicom_space,mpij)
          lcut=lcut0
          if (lcut == 0) then
            latom(k)=j
            ifl(j)=1
            j=natom+2
          end if
        end if
        j=j+1
      end do
      if (j == natom+1) then
        j=1
        do while (j <= natom)
          if (ifl(j) == 0) then
            latom(k)=j
            ifl(j)=1
            j=natom
          end if
          j=j+1
        end do
      end if
      k=k+1
    end do
    deallocate(ifl)
  else
    do i=1,natom
      latom(i)=i
    end do
  end if
  return
end subroutine tools_rodr_sum_nlp


subroutine tools_maginitrot( &
 key_natpri_in,nspv,natom,num_spe,num_atcell,nradmx,npoint,ncpx,ncpy,ncpz, & ! <
 natpri,natpri_inf,indspe,nradct,radial,dradial,wt,                        & ! <
 polconmag,                                                                & ! <
 rhosmt,rhotrur,rhosmtr)                                                     ! X
implicit none
real*8, parameter:: eps_p=1.0d-10
integer,intent(in)   :: key_natpri_in
integer,intent(in)   :: nspv,natom,num_spe,num_atcell,nradmx,npoint,ncpx,ncpy,ncpz 
integer,intent(in)   :: natpri(natom),natpri_inf(natom),indspe(natom),nradct(num_spe)
real*8, intent(in)   :: radial(nradmx,num_spe),dradial(nradmx,num_spe),wt(npoint)
real*8, intent(in)   :: polconmag(3,natom)
real*8, intent(inout):: rhosmt(ncpx,ncpy,ncpz,nspv) 
real*8, intent(inout):: rhotrur(nradmx,npoint,nspv,num_atcell),rhosmtr(nradmx,npoint,nspv,num_atcell)
integer:: na,ipri,ispe,nrad,il,ir,ns,ns2,ix,iy,iz 
real*8 :: pscl,chrg,rfac  
real*8,allocatable:: pol(:),spol1(:),spol2(:) 

  allocate( pol(3),spol1(3),spol2(3) )

  do na= 1,natom
  if (natpri(na)==key_natpri_in) then
    ipri= natpri_inf(na)
    ispe= indspe(na)
    nrad= nradct(indspe(na))

    pscl= dsqrt(dabs(polconmag(1,na)**2+polconmag(2,na)**2+polconmag(3,na)**2))
    if (pscl>eps_p) then 

      spol2(1)= polconmag(1,na)/pscl 
      spol2(2)= polconmag(2,na)/pscl 
      spol2(3)= polconmag(3,na)/pscl 

      spol1(:)= 0.0d0
      do ns= 1,4
        ns2= max(1,ns-1)
        do il= 1,npoint
        do ir= 2,nrad-1
          rfac= wt(il)*dradial(ir,ispe)*radial(ir,ispe)**2
          spol1(ns2)= spol1(ns2)+rfac*rhotrur(ir,il,ns,ipri)
        enddo
        enddo
        if (ns<3) then
          spol1(ns2)= -spol1(ns2)
        else 
          spol1(ns2)= 2.0d0*spol1(ns2) 
        endif
      enddo
      pscl= dsqrt(dabs(spol1(1)**2+spol1(2)**2+spol1(3)**3))

      if (pscl>eps_p) then
        spol1(1)= spol1(1)/pscl 
        spol1(2)= spol1(2)/pscl 
        spol1(3)= spol1(3)/pscl 
      else
        spol1(:)= 0.0d0
      endif 

      do il= 1,npoint 
      do ir= 2,nrad
        chrg  = rhotrur(ir,il,1,ipri)+rhotrur(ir,il,nspv-2,ipri)
        pol(1)= rhotrur(ir,il,1,ipri)-rhotrur(ir,il,nspv-2,ipri)
        pol(2)= 2.0d0*rhotrur(ir,il,nspv-1,ipri)
        pol(3)= 2.0d0*rhotrur(ir,il,nspv-0,ipri)
        pscl= pol(1)*spol1(1)+pol(2)*spol1(2)+pol(3)*spol1(3)
        pol(1)= pscl*spol2(1)
        pol(2)= pscl*spol2(2)
        pol(3)= pscl*spol2(3)
        rhotrur(ir,il,nspv-3,ipri)= (chrg+pol(1))/2.0d0  
        rhotrur(ir,il,nspv-2,ipri)= (chrg-pol(1))/2.0d0  
        rhotrur(ir,il,nspv-1,ipri)= pol(2)/2.0d0   
        rhotrur(ir,il,nspv-0,ipri)= pol(3)/2.0d0   
      enddo
      enddo

    endif 

    do il= 1,npoint  
    do ir= 2,nrad
      chrg= rhosmtr(ir,il,1,ipri)+rhosmtr(ir,il,nspv-2,ipri) 
      rhosmtr(ir,il,nspv-3,ipri)= chrg/2.0d0 
      rhosmtr(ir,il,nspv-2,ipri)= chrg/2.0d0 
      rhosmtr(ir,il,nspv-1,ipri)= 0.0d0 
      rhosmtr(ir,il,nspv-0,ipri)= 0.0d0 
    enddo
    enddo

  endif 
  enddo ! na  

  do iz= 1,ncpz
  do iy= 1,ncpy
  do ix= 1,ncpx
    chrg= rhosmt(ix,iy,iz,1)+rhosmt(ix,iy,iz,nspv-2)
    rhosmt(ix,iy,iz,nspv-3)= chrg/2.0d0 
    rhosmt(ix,iy,iz,nspv-2)= chrg/2.0d0 
    rhosmt(ix,iy,iz,nspv-1)= 0.0d0 
    rhosmt(ix,iy,iz,nspv-0)= 0.0d0 
  enddo
  enddo
  enddo

  deallocate( pol,spol1,spol2 )

end subroutine tools_maginitrot


subroutine tools_spinflip_svec( &
 nflip,nrc,nums,ncol,ncpx,ncpy,ncpz,neigmx,numk, & ! <
 svecre,sveccm,sval)                             ! X
use mod_stopp
use mod_mpi ! test 
implicit none
integer,   intent(in)   ::nflip,nrc,nums,ncol,ncpx,ncpy,ncpz,neigmx,numk
real*8,    intent(inout)::svecre(ncpx*(1-nrc)+nrc,ncpy*(1-nrc)+nrc,ncpz*(1-nrc)+nrc, &
                                 neigmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk*(1-nrc)+nrc)
complex*16,intent(inout)::sveccm(ncpx*nrc-nrc+1,ncpy*nrc-nrc+1,ncpz*nrc-nrc+1, &
                                 neigmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1)
real*8,    intent(inout)::sval(neigmx,nums+1-ncol,numk)
integer::nk
real*8,    allocatable::rbuf(:,:,:,:)
complex*16,allocatable::cbuf(:,:,:,:)

  if (nflip==1) then
    if ((nums/=2).or.(ncol/=1)) call stopp('tools_spinflip_svec: (nums/=2) or (ncol/=1)')
    if (nrc==0) then
      allocate(rbuf(ncpx,ncpy,ncpz,neigmx))
      do nk= 1,numk
        rbuf(:,:,:,:)          = svecre(:,:,:,:,1   ,nk)
        svecre(:,:,:,:,1   ,nk)= svecre(:,:,:,:,nums,nk)
        svecre(:,:,:,:,nums,nk)= rbuf(:,:,:,:)
      enddo
      deallocate(rbuf)
    else
      allocate(cbuf(ncpx,ncpy,ncpz,neigmx))
      do nk= 1,numk
        cbuf(:,:,:,:)          = sveccm(:,:,:,:,1   ,nk)
        sveccm(:,:,:,:,1   ,nk)= sveccm(:,:,:,:,nums,nk)
        sveccm(:,:,:,:,nums,nk)= cbuf(:,:,:,:)
      enddo
      deallocate(cbuf)
    endif
    allocate(rbuf(neigmx,1,1,1))
    do nk= 1,numk
      rbuf(:,1,1,1)  = sval(:,1,   nk) 
      sval(:,1   ,nk)= sval(:,nums,nk) 
      sval(:,nums,nk)= rbuf(:,1,1,1) 
    enddo
    deallocate(rbuf)
  endif

end subroutine tools_spinflip_svec


subroutine tools_spinflip_rho( &
 nflip,nspv,ncpx,ncpy,ncpz,nradmx,npoint,num_atcell, & ! <
 rhosmt,rhotrur,rhosmtr)                               ! X
use mod_stopp
implicit none
integer,   intent(in)   ::nflip,nspv,ncpx,ncpy,ncpz,nradmx,npoint,num_atcell
real*8,    intent(inout)::rhosmt(ncpx,ncpy,ncpz,nspv)
real*8,    intent(inout)::rhotrur(nradmx,npoint,nspv,num_atcell),rhosmtr(nradmx,npoint,nspv,num_atcell)
integer::na 
real*8,    allocatable::rbuf(:,:,:)

  if (nflip==1) then 
    if (nspv/=2) call stopp('tools_spinflip_rho: nspv /= 2')
    allocate(rbuf(ncpx,ncpy,ncpz))
    rbuf(:,:,:)       = rhosmt(:,:,:,1   )
    rhosmt(:,:,:,1   )= rhosmt(:,:,:,nspv)
    rhosmt(:,:,:,nspv)= rbuf(:,:,:)  
    deallocate(rbuf)
    allocate(rbuf(nradmx,npoint,1))
    do na= 1,num_atcell
      rbuf(:,:,1)         = rhotrur(:,:,1   ,na)
      rhotrur(:,:,1   ,na)= rhotrur(:,:,nspv,na)
      rhotrur(:,:,nspv,na)= rbuf(:,:,1)  
      rbuf(:,:,1)         = rhosmtr(:,:,1   ,na)
      rhosmtr(:,:,1   ,na)= rhosmtr(:,:,nspv,na)
      rhosmtr(:,:,nspv,na)= rbuf(:,:,1)  
    enddo
    deallocate(rbuf)
  endif 

end subroutine tools_spinflip_rho 


subroutine tools_deri_grdch(nradmx,radial,drdi,ro,drr,ddrr)
implicit none
integer,intent(in) :: nradmx
real*8, intent(in) :: radial(nradmx)
real*8, intent(in) :: drdi(nradmx)
real*8, intent(in) :: ro(nradmx)
real*8, intent(out):: drr(nradmx),ddrr(nradmx)

! .. Local Scholar ..
real*8  :: dfx,d2fx,ddrdi,rin
integer :: ir

!$omp parallel default(shared),private(ir,dfx,d2fx,ddrdi,rin)
!$omp do
  do ir=2,nradmx-1
   dfx  =   (ro(ir+1) - ro(ir-1))*0.5d0
   d2fx =    ro(ir-1) -2.d0*    ro(ir) +     ro(ir+1)
   ddrdi=radial(ir-1) -2.d0*radial(ir) + radial(ir+1)
   rin=1.0d0/drdi(ir)
    drr(ir) = dfx*rin
   ddrr(ir) = d2fx*rin*rin - dfx*ddrdi*rin*rin*rin
  end do
!$omp end parallel
  drr(1)=drr(2)
  drr(nradmx)=drr(nradmx-1)
  ddrr(1)=drr(2)
  ddrr(nradmx)=drr(nradmx-1)

  return
end subroutine tools_deri_grdch


subroutine tools_deri_yylm(npoint,lrhomx,point,dylm_dtheta,d2ylm_dtheta2, &
                     dylm_dphi,d2ylm_dphi2,d2ylm_dtheta_dphi)
implicit none
integer,intent(in)    :: npoint,lrhomx
real*8, intent(in)    :: point(npoint,3) 
real*8, intent(out)   :: dylm_dtheta(npoint,25), d2ylm_dtheta2(npoint,25), &
                         dylm_dphi(npoint,25), d2ylm_dphi2(npoint,25),     &
                         d2ylm_dtheta_dphi(npoint,25)
integer :: j
real*8  :: pi,r,x,y,z,x2,y2,z2,r2,x4,y4,z4,rin,rin2,rin3,rin4,ct,st,cf,sf
real*8  :: ct2,st2,cf2,sf2,ct3,st3,cf3,sf3,ct4,st4,cf4,sf4
  pi = 4.d0*atan(1.d0)

  do j=1,npoint
    r=1.d0
    x=point(j,1)*r   ! for comparing with Fleur
    y=point(j,2)*r   ! for comparing with Fleur
    z=point(j,3)*r   ! for comparing with Fleur
    x2=x*x
    y2=y*y
    z2=z*z
    r2=r*r
    x4=x2*x2
    y4=y2*y2
    z4=z2*z2
    rin=1.0d0/r
    rin2=rin*rin
    rin3=rin2*rin
    rin4=rin2*rin2
    ct=z/r
    st=dsqrt(x*x+y*y)/r
    cf=x/dsqrt(x*x+y*y)
    sf=y/dsqrt(x*x+y*y)
    ct2=ct*ct
    st2=st*st
    cf2=cf*cf
    sf2=sf*sf
    ct3=ct2*ct
    st3=st2*st
    cf3=cf2*cf
    sf3=sf2*sf
    ct4=ct2*ct2
    st4=st2*st2
    cf4=cf2*cf2
    sf4=sf2*sf2

!   start of L=0
     dylm_dtheta(j,1)     = 0.d0
    d2ylm_dtheta2(j,1)    = 0.d0
     dylm_dphi(j,1)       = 0.d0
    d2ylm_dphi2(j,1)      = 0.d0
    d2ylm_dtheta_dphi(j,1)= 0.d0
!   end   of L=0

!   start of L=1
     dylm_dtheta(j,2)     =  dsqrt(0.75d0/pi)*ct*cf
    d2ylm_dtheta2(j,2)    = -dsqrt(0.75d0/pi)*st*cf
     dylm_dphi(j,2)       = -dsqrt(0.75d0/pi)*st*sf
    d2ylm_dphi2(j,2)      = -dsqrt(0.75d0/pi)*st*cf
    d2ylm_dtheta_dphi(j,2)= -dsqrt(0.75d0/pi)*ct*sf

     dylm_dtheta(j,3)     =  dsqrt(0.75d0/pi)*ct*sf
    d2ylm_dtheta2(j,3)    = -dsqrt(0.75d0/pi)*st*sf
     dylm_dphi(j,3)       =  dsqrt(0.75d0/pi)*st*cf
    d2ylm_dphi2(j,3)      = -dsqrt(0.75d0/pi)*st*sf
    d2ylm_dtheta_dphi(j,3)=  dsqrt(0.75d0/pi)*ct*cf

     dylm_dtheta(j,4)     = -dsqrt(0.75d0/pi)*st
    d2ylm_dtheta2(j,4)    = -dsqrt(0.75d0/pi)*ct
     dylm_dphi(j,4)       =  0.d0
    d2ylm_dphi2(j,4)      =  0.d0
    d2ylm_dtheta_dphi(j,4)=  0.d0
!   end   of L=1

!   start of L=2
     dylm_dtheta(j,5)     =  dsqrt(15.0d0/(4.0d0*pi))*st*ct*(cf2-sf2)
    d2ylm_dtheta2(j,5)    =  dsqrt(15.0d0/(4.0d0*pi))*(ct2-st2)*(cf2-sf2)
     dylm_dphi(j,5)       = -dsqrt(15.d0/pi)*st2*sf*cf
    d2ylm_dphi2(j,5)      = -dsqrt(15.d0/pi)*st2*(cf2-sf2)
    d2ylm_dtheta_dphi(j,5)= -2.0d0*dsqrt(15.d0/pi)*st*ct*sf*cf

     dylm_dtheta(j,6)     = -0.5d0*dsqrt(15.d0/pi)*(ct2-st2)*cf
    d2ylm_dtheta2(j,6)    =  2.0d0*dsqrt(15.d0/pi)*st*ct*cf
     dylm_dphi(j,6)       =  0.5d0*dsqrt(15.d0/pi)*st*ct*sf
    d2ylm_dphi2(j,6)      =  0.5d0*dsqrt(15.d0/pi)*st*ct*cf
    d2ylm_dtheta_dphi(j,6)=  0.5d0*dsqrt(15.d0/pi)*(ct2-st2)*sf

     dylm_dtheta(j,7)     = -1.5d0*dsqrt(5.d0/pi)*ct*st
    d2ylm_dtheta2(j,7)    = -1.5d0*dsqrt(5.d0/pi)*(ct2-st2)
     dylm_dphi(j,7)       =  0.d0
    d2ylm_dphi2(j,7)      =  0.d0
    d2ylm_dtheta_dphi(j,7)=  0.d0

     dylm_dtheta(j,8)     = -0.5d0*dsqrt(15.0d0/pi)*(ct2-st2)*sf
    d2ylm_dtheta2(j,8)    =  2.0d0*dsqrt(15.0d0/pi)*st*ct*sf
     dylm_dphi(j,8)       = -0.5d0*dsqrt(15.0d0/pi)*st*ct*cf
    d2ylm_dphi2(j,8)      =  0.5d0*dsqrt(15.0d0/pi)*st*ct*sf
    d2ylm_dtheta_dphi(j,8)= -0.5d0*dsqrt(15.0d0/pi)*(ct2-st2)*cf

     dylm_dtheta(j,9)     =  dsqrt(15.0d0/pi)*st*ct*sf*cf
    d2ylm_dtheta2(j,9)    =  dsqrt(15.0d0/pi)*(ct2-st2)*sf*cf
     dylm_dphi(j,9)       =  0.5d0*dsqrt(15.0d0/pi)*st2*(cf2-sf2)
    d2ylm_dphi2(j,9)      = -2.0d0*dsqrt(15.0d0/pi)*st2*sf*cf
    d2ylm_dtheta_dphi(j,9)=  dsqrt(15.0d0/pi)*st*ct*(cf2-sf2)
!   end   of L=2

!   start of L=3
     dylm_dtheta(j,10)     = 0.75d0*dsqrt(7.d0/pi)*st*(1.0d0-5.0d0*ct2)
    d2ylm_dtheta2(j,10)    = 0.75d0*dsqrt(7.d0/pi)*(11.0d0*st2-4.0d0*ct2)*ct
     dylm_dphi(j,10)       = 0.d0
    d2ylm_dphi2(j,10)      = 0.d0
    d2ylm_dtheta_dphi(j,10)= 0.d0

     dylm_dtheta(j,11)     = -0.25d0*dsqrt(10.5d0/pi)*(4.0d0*ct2-11.0d0*st2)*ct*cf
    d2ylm_dtheta2(j,11)    = -0.25d0*dsqrt(10.5d0/pi)*(11.0d0*st2-34.0d0*ct2)*st*cf
     dylm_dphi(j,11)       =  0.25d0*dsqrt(10.5d0/pi)*st*sf*(5.0d0*ct2-1.0d0)
    d2ylm_dphi2(j,11)      =  0.25d0*dsqrt(10.5d0/pi)*st*cf*(5.0d0*ct2-1.0d0)
    d2ylm_dtheta_dphi(j,11)=  0.25d0*dsqrt(10.5d0/pi)*(4.0d0*ct2-11.0d0*st2)*ct*sf

     dylm_dtheta(j,12)     =  0.25d0*dsqrt(105.d0/pi)*(2.0d0*st*ct2-st3)*(cf2-sf2)
    d2ylm_dtheta2(j,12)    =  0.25d0*dsqrt(105.d0/pi)*(2.0d0*ct2-7.0d0*st2)*ct*(cf2-sf2)
     dylm_dphi(j,12)       = -dsqrt(105.d0/pi)*st2*ct*sf*cf
    d2ylm_dphi2(j,12)      = -dsqrt(105.d0/pi)*st2*ct*(cf2-sf2)
    d2ylm_dtheta_dphi(j,12)= -dsqrt(105.d0/pi)*(2.0d0*ct2-st2)*ct*sf*cf

     dylm_dtheta(j,13)     = -0.75d0*dsqrt(17.5d0/pi)*st2*ct*cf*(cf2-3.0d0*sf2)
    d2ylm_dtheta2(j,13)    = -0.75d0*dsqrt(17.5d0/pi)*(2.0d0*st*ct2-st3)*cf*(cf2-3.0d0*sf2)
     dylm_dphi(j,13)       = -0.75d0*dsqrt(17.5d0/pi)*st3*sf*(sf2-3.0d0*cf2)
    d2ylm_dphi2(j,13)      = -2.25d0*dsqrt(17.5d0/pi)*st3*cf*(3.0d0*sf2-cf2)
    d2ylm_dtheta_dphi(j,13)= -2.25d0*dsqrt(17.5d0/pi)*st2*ct*sf*(sf2-3.0d0*cf2)

     dylm_dtheta(j,14)     = -0.25d0*dsqrt(10.5d0/pi)*(4.0d0*ct2-11.0d0*st2)*ct*sf
    d2ylm_dtheta2(j,14)    = -0.25d0*dsqrt(10.5d0/pi)*(11.0d0*st2-34.0d0*ct2)*st*sf
     dylm_dphi(j,14)       = -0.25d0*dsqrt(10.5d0/pi)*st*(5.0d0*ct2-1.0d0)*cf
    d2ylm_dphi2(j,14)      =  0.25d0*dsqrt(10.5d0/pi)*st*(5.0d0*ct2-1.0d0)*sf
    d2ylm_dtheta_dphi(j,14)= -0.25d0*dsqrt(10.5d0/pi)*(4.0d0*ct2-11.0d0*st2)*ct*cf

     dylm_dtheta(j,15)     =  0.5d0*dsqrt(105.0d0/pi)*(2.0d0*ct2-st2)*st*sf*cf
    d2ylm_dtheta2(j,15)    =  0.5d0*dsqrt(105.0d0/pi)*(2.0d0*ct2-7.0d0*st2)*ct*sf*cf
     dylm_dphi(j,15)       =  0.5d0*dsqrt(105.0d0/pi)*st2*ct*(cf2-sf2)
    d2ylm_dphi2(j,15)      = -2.0d0*dsqrt(105.0d0/pi)*st2*ct*sf*cf
    d2ylm_dtheta_dphi(j,15)=  0.5d0*dsqrt(105.0d0/pi)*(2.0d0*ct2-st2)*st*(cf2-sf2)

     dylm_dtheta(j,16)     = -0.75d0*dsqrt(17.5d0/pi)*st2*ct*sf*(3.0d0*cf2-sf2)
    d2ylm_dtheta2(j,16)    = -0.75d0*dsqrt(17.5d0/pi)*(2.0d0*ct2-st2)*st*sf*(3.0d0*cf2-sf2)
     dylm_dphi(j,16)       = -0.75d0*dsqrt(17.5d0/pi)*st3*cf*(cf2-3.0d0*sf2)
    d2ylm_dphi2(j,16)      = -2.25d0*dsqrt(17.5d0/pi)*st3*sf*(sf2-3.0d0*cf2)
    d2ylm_dtheta_dphi(j,16)= -2.25d0*dsqrt(17.5d0/pi)*st2*ct*cf*(cf2-3.0d0*sf2)
!   end   of L=3

!   start of L=4
     dylm_dtheta(j,17)     =  15.d0/4.d0/dsqrt(pi)*(3.0d0-7.0d0*ct2)*ct*st
    d2ylm_dtheta2(j,17)    = -15.d0/4.d0/dsqrt(pi)*(3.d0*st4-21.0d0*st2*ct2+4.0d0*ct4)
     dylm_dphi(j,17)       =  0.d0
    d2ylm_dphi2(j,17)      =  0.d0
    d2ylm_dtheta_dphi(j,17)=  0.d0

     dylm_dtheta(j,18)     =  0.75d0*dsqrt(2.5d0/pi)*(4.0d0*ct4-21.0d0*st2*ct2+3.0d0*st4)*cf
    d2ylm_dtheta2(j,18)    =  1.5d0*dsqrt(2.5d0/pi)*(27.0d0*st3*ct-29.0d0*st*ct3)*cf
     dylm_dphi(j,18)       = -0.75d0*dsqrt(2.5d0/pi)*(7.0d0*ct3-3.0d0*ct)*st*sf
    d2ylm_dphi2(j,18)      = -0.75d0*dsqrt(2.5d0/pi)*(7.0d0*ct3-3.0d0*ct)*st*cf
    d2ylm_dtheta_dphi(j,18)= -0.75d0*dsqrt(2.5d0/pi)*(4.0d0*ct4-21.0d0*st2*ct2+3.0d0*st4)*sf

      dylm_dphi(j,19)       =  3.d0*dsqrt(5.d0)/2.d0/dsqrt(pi)*(r2-7.d0*z2)*x*y*rin4
     dylm_dtheta(j,19)     =  1.5d0*dsqrt(5.d0/pi)*(3.0d0*ct2-4.0d0*st2)*st*ct*(cf2-sf2)
    d2ylm_dtheta2(j,19)    =  1.50*dsqrt(5.0d0)/dsqrt(pi)*(4.0d0*st4-21.0d0*st2*ct2+3.0d0*ct4)*(cf2-sf2)
     dylm_dphi(j,19)       = -1.5d0*dsqrt(5.d0/pi)*(7.d0*ct2-1.0d0)*st2*sf*cf
    d2ylm_dphi2(j,19)      = -1.5d0*dsqrt(5.d0/pi)*(7.d0*ct2-1.0d0)*st2*(cf2-sf2)
    d2ylm_dtheta_dphi(j,19)= -6.0d0*dsqrt(5.d0/pi)*(3.d0*ct2-4.0d0*st2)*st*ct*sf*cf

     dylm_dtheta(j,20)     =  0.75d0*dsqrt(17.5d0/pi)*(3.0d0*st2*ct2-st4)*(4.0d0*cf3-3.0d0*cf)
    d2ylm_dtheta2(j,20)    =  1.5d0*dsqrt(17.5d0/pi)*(3.0d0*st*ct3-5.0d0*st3*ct)*(4.0d0*cf3-3.0d0*cf)
     dylm_dphi(j,20)       =  2.25d0*dsqrt(17.5d0/pi)*st3*ct*sf*(sf2-3.0d0*cf2)
    d2ylm_dphi2(j,20)      =  6.75d0*dsqrt(17.5d0/pi)*st3*ct*cf*(3.0d0*sf2-cf2)
    d2ylm_dtheta_dphi(j,20)=  2.25d0*dsqrt(17.5d0/pi)*(3.0d0*st2*ct2-st4)*sf*(sf2-3.0d0*cf2)

     dylm_dtheta(j,21)     = 12.0d0/16.0d0*dsqrt(35.0d0/pi)*st**3*ct*(1.0d0-8.0d0*cf*cf*sf*sf)
    d2ylm_dtheta2(j,21)    = 0.75d0*dsqrt(35.0d0/pi)*(3.0d0*st*st*ct*ct-st**4)*(1.0d0-8.0d0*cf*cf*sf*sf)
     dylm_dphi(j,21)       = 3.0d0*dsqrt(35.0d0/pi)*st**4*(cf*sf**3-cf**3*sf)
    d2ylm_dphi2(j,21)      =-3.0d0*dsqrt(35.0d0/pi)*st**4*(cf**4-6.0d0*cf*cf*sf*sf+sf**4)
    d2ylm_dtheta_dphi(j,21)= 12.0d0*dsqrt(35.0d0/pi)*st**3*ct*(cf*sf**3-cf**3*sf)

     dylm_dtheta(j,22)     =  0.75d0*dsqrt(2.5d0/pi)*(4.0d0*ct4-21.0d0*st2*ct2+3.0d0*st4)*sf
    d2ylm_dtheta2(j,22)    =  1.5d0*dsqrt(2.5d0/pi)*(27.0d0*st3*ct-29.0d0*st*ct3)*sf
     dylm_dphi(j,22)       =  0.75d0*dsqrt(2.5d0/pi)*(4.0d0*ct2-3.0d0*st2)*ct*st*cf
    d2ylm_dphi2(j,22)      = -0.75d0*dsqrt(2.5d0/pi)*(4.0d0*ct2-3.0d0*st2)*ct*st*sf
    d2ylm_dtheta_dphi(j,22)=  0.75d0*dsqrt(2.5d0/pi)*(4.0d0*ct4-21.0d0*st2*ct2+3.0d0*st4)*cf

     dylm_dtheta(j,23)     =  3.0d0*dsqrt(5.0d0/pi)*(3.0d0*st*ct3-4.0d0*st3*ct)*cf*sf
    d2ylm_dtheta2(j,23)    =  3.0d0*dsqrt(5.0d0/pi)*(3.0d0*ct4-21.0d0*st2*ct2+4.0d0*st4)*cf*sf
     dylm_dphi(j,23)       =  0.75d0*dsqrt(5.0d0/pi)*(6.0d0*ct2-st2)*st2*(cf2-sf2)
    d2ylm_dphi2(j,23)      = -3.0d0*dsqrt(5.0d0/pi)*(6.0d0*ct2-st2)*st2*sf*cf
    d2ylm_dtheta_dphi(j,23)=  3.0d0*dsqrt(5.0d0/pi)*(3.0d0*st*ct3-4.0d0*st3*ct)*(cf2-sf2)

     dylm_dtheta(j,24)     = -0.75d0*dsqrt(17.5d0/pi)*(3.0d0*st2*ct2-st4)*(3.0d0*sf-4.0d0*sf3)
    d2ylm_dtheta2(j,24)    = -1.5d0*dsqrt(17.5d0/pi)*(3.0d0*st*ct3-5.0d0*st3*ct)*(3.0d0*sf-4.0d0*sf3)
     dylm_dphi(j,24)       = -2.25d0*dsqrt(17.5d0/pi)*st3*ct*cf*(cf2-3.0d0*sf2)
    d2ylm_dphi2(j,24)      = -6.75d0*dsqrt(17.5d0/pi)*st3*ct*sf*(sf2-3.0d0*cf2)
    d2ylm_dtheta_dphi(j,24)= -2.25d0*dsqrt(17.5d0/pi)*(3.0d0*st*ct3-st4)*cf*(cf2-3.0d0*sf2)

     dylm_dtheta(j,25)     = 3.0d0*dsqrt(35.0d0/pi)*st**3*ct*sf*cf*(cf*cf-sf*sf)
    d2ylm_dtheta2(j,25)    = 3.0d0*dsqrt(35.0d0/pi)*(3.0d0*st*st*ct*ct-st**4)*sf*cf*(cf*cf-sf*sf)
     dylm_dphi(j,25)       = 0.75d0*dsqrt(35.0d0/pi)*st**4*(cf**4-6.0d0*sf*sf*cf*cf+sf**4)
    d2ylm_dphi2(j,25)      = 12.0d0*dsqrt(35.0d0/pi)*st**4*(sf**3*cf-sf*cf**3)
    d2ylm_dtheta_dphi(j,25)= 3.0d0*dsqrt(35.0d0/pi)*st**3*ct*(cf**4-6.0d0*sf*sf*cf*cf+sf**4)
! end of  L=4

  end do

  return
end subroutine tools_deri_yylm


end module
