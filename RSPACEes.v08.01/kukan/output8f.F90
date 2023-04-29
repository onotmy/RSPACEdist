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
! **********  output8f.F90 12/05/2022-01  **********

module mod_output
contains

subroutine output_compcond( &
 ndisp,njfile,nperi,numkmx,nmesh,nums,nspv,npre,                        & ! <
 xmax,ymax,zmax,dx,dy,dz,skpxyz,nwskp,                                  & ! <
 tnumele,chrgd,biasx,biasy,biasz,tmstep,tf,fcut,eta,etamag)   ! <
use mod_stopp
implicit none
integer, intent(in) :: ndisp,njfile,nperi,numkmx,nmesh,nums,nspv,npre
integer, intent(in) :: nwskp(numkmx)
real*8,  intent(in) :: &
 xmax,ymax,zmax,dx,dy,dz,skpxyz(numkmx,3),tnumele,chrgd, &
 biasx,biasy,biasz,tmstep,tf,fcut,eta,etamag(2)
integer :: loop,jfile,nk
real*8  :: pi

  pi= 4.0d0*datan(1.0d0)

  do loop=1,2
    if (loop==1) jfile= ndisp
    if (loop==2) jfile= njfile

    write (jfile,*,err=9999) 'cell size [lx,ly,lz (a.u.)]'
    write (jfile,fmt='(3d18.9,i6)',err=9999) 2.0d0*xmax,2.0d0*ymax,2.0d0*zmax
    write (jfile,*,err=9999) 'coarse grid spacing [hx,hy,hz (a.u.)]'
    write (jfile,fmt='(3d18.9,i6)',err=9999) dx,dy,dz
    write (jfile,*,err=9999) 'dense grid spacing [hx,hy,hz (a.u.)]'
    write (jfile,fmt='(3d18.9,i6)',err=9999) dx/nmesh,dy/nmesh,dz/nmesh
    write (jfile,*,err=9999) 'cutoff energy of wave functions [gmax_coarse,gmax_dense (Ry)]'
    write (jfile,fmt='(3d18.9,i6)',err=9999) (pi/(dmax1(dx,dy,dz)))**2,(nmesh*pi/(dmax1(dx,dy,dz)))**2
    if (nperi .gt. 0) then
      write (jfile,*,err=9999) 'sample k point [skx,sky,skz (2pi/l)]  weight'
      do nk=1,numkmx
        write (jfile,fmt='(3d18.9,i6)',err=9999) &
         skpxyz(nk,1)/(pi/xmax),skpxyz(nk,2)/(pi/ymax),skpxyz(nk,3)/(pi/zmax),nwskp(nk)
      end do
    end if
    write (jfile,*,err=9999) 'number of electron       charge         spin'
    write (jfile,fmt='(3d18.9,i6)',err=9999) tnumele,chrgd,dble(nums)
    write (jfile,*,err=9999) 'electric field [ex,ey,ez (a.u./a.u.)]'
    write (jfile,fmt='(3d18.9,i6)',err=9999) biasx,biasy,biasz
    write (jfile,*,err=9999) 'time step (a.u.)  (fs)'
    write (jfile,fmt='(3d18.9,i6)',err=9999) tmstep,tmstep*0.0241888d0
    write (jfile,*,err=9999) 'thermo (k)    kt (a.u.)'
    write (jfile,fmt='(3d18.9,i6)',err=9999) tf/3.1670d-6,tf
    write (jfile,*,err=9999) 'cufoff of force acting on atoms (a.u.)'
    write (jfile,fmt='(3d18.9,i6)',err=9999) fcut
    write (jfile,*,err=9999) 'charge mixing ratio of scf'
    write (jfile,fmt='(1d18.9,i6)',err=9999) eta
    if (nspv>1) then
      write (jfile,*,err=9999) 'mag. mix. rat., weight of mag. in F_broyd'
      write (jfile,fmt='(2d18.9,i6)',err=9999) etamag(1:2)
    endif
    if (npre .eq. 0) then
      write (jfile,*,err=9999) 'Wavefunctions & electron density are read from files.'
    else
      write (jfile,*,err=9999) 'Wavefunctions & electron density are generated automatically.'
    end if
  enddo ! loop

return
  9999 call stopp('output_compcond: error writing output')
end subroutine output_compcond


subroutine output_countdat( &
 ndisp,njfile,iscf,nperi,npopshow,nso,nspv,nums,ncol,natom,npolcondim,numkmx,neigmx,                 & ! <
 npolcon,npolcon2,key_polcon_atoms,key_polcon_asa,key_polcon2_none,                                      & ! <
 nwskp,nwskptot,                                                                                         & ! <
 ksconv,ksitmax,tnumele,ferm,spinpol,polconb,atmpole,diff_charge,residual_states_wfc,sval_wfc,fnele_wfc)   ! < 
use mod_stopp
implicit none
integer, intent(in) :: ndisp,njfile,iscf,nperi,npopshow,nso,nspv,nums,ncol,natom,npolcondim,numkmx,neigmx
integer, intent(in) :: npolcon,npolcon2(0:npolcondim),key_polcon_atoms,key_polcon_asa,key_polcon2_none
integer, intent(in) :: nwskp(numkmx),nwskptot
logical, intent(in) :: ksconv
integer, intent(in) :: ksitmax
real*8,  intent(in) :: tnumele,ferm
real*8,  intent(in) :: spinpol(max(1,nspv-1),0:natom*min(1,nspv-1)),polconb(3,npolcondim)
real*8,  intent(in) :: atmpole(10,natom)
real*8,  intent(in) :: diff_charge
real*8,  intent(in) :: residual_states_wfc(neigmx,nums+1-ncol,numkmx)
real*8,  intent(in) :: sval_wfc(           neigmx*(1+nso*(2-ncol)),1+(nums-1)*(2-ncol)*(1-nso),numkmx)
real*8,  intent(in) :: fnele_wfc(          neigmx*(1+nso*(2-ncol)),1+(nums-1)*(2-ncol)*(1-nso),numkmx)

integer,parameter:: nout= 1 ! after how any iterations should the long eigenvalue list be written to ndisp ? 

integer  :: loop,jfile,ns,na,nk,l 
real*8   :: wsqresmx,svalmin,fwghtinv
character:: chnspv*1
logical  :: l_outp

  write(chnspv,fmt='(i1)') nspv-1

  do loop=1,2
    if (loop==1) then
      jfile= ndisp
      if (nout>0) then
        l_outp= (mod(iscf-1,nout)==0)
      else
        l_outp= .false.
      endif 
    else
      jfile= njfile
      l_outp= .true.
    endif 

    write (jfile,*,err=9999) 'total charge',tnumele,'(electrons)'
    write (jfile,*,err=9999) 'fermi level',ferm,'(hartree)'
    if (nspv>1) then
      if (nspv==2) then
        write (jfile,*,err=9999) 'spin polarization',spinpol(1,0),'(\mu_B)'
      else
        write (jfile,fmt='(A,3(1x,e11.4))',err=9999) '(z,x,y)-polarization (in \mu_B): ', &
         (spinpol(ns,0),ns=1,3)
      endif
      do na=1,natom
        if ( ((npolcon==key_polcon_atoms).or.(npolcon==key_polcon_asa)) .and. (npolcon2(na)/=key_polcon2_none) ) then
          write(jfile,fmt='(3x,"atom #",i4,":  spinpol=",'//chnspv//'(1x,e11.4)," , B_con=",'//chnspv//'(1x,e10.3))',err=9999) &
           na, (spinpol(ns,na),ns=1,nspv-1), (polconb(ns,na),ns=1,nspv-1)
        else
          write(jfile,fmt='(3x,"atom #",i4,":  spinpol=",'//chnspv//'(1x,e11.4))',err=9999) &
           na, (spinpol(ns,na),ns=1,nspv-1)
        endif
      enddo
    end if
    if ((nperi .lt. 3) .or. (npopshow .eq. 1)) then
      write (jfile,*,err=9999) 'atomic charge','(electrons)'
      do na=1,natom
         write (jfile,*,err=9999) na,atmpole(1,na)
      end do
    end if
    if (l_outp) call output_evlist( &
     jfile,nso,nums,ncol,numkmx,neigmx,                     & ! <
     nwskp,nwskptot,residual_states_wfc,sval_wfc,fnele_wfc)   ! <
    wsqresmx= 0.0d0
    svalmin= sval_wfc(1,1,1)
    do nk=1,numkmx
      if (nwskp(nk)/=0) then
        fwghtinv= dble(nwskptot*nums)/dble(nwskp(nk)*2)
      else
        fwghtinv= 0.0d0
      endif
      if ((nso==0).or.(ncol==2)) then
        do ns=1,nums+1-ncol
          do l=1,neigmx
            wsqresmx= dmax1(wsqresmx,dsqrt(residual_states_wfc(l,ns,nk))*fnele_wfc(l,ns,nk)*fwghtinv)
            svalmin= dmin1(svalmin,sval_wfc(l,ns,nk))
          end do
        end do
      else
        do ns=1,nums
          do l=1,neigmx
            wsqresmx= dmax1(wsqresmx,dsqrt(residual_states_wfc(l,ns,nk)))
            svalmin= dmin1(svalmin,sval_wfc(l,ns,nk))
          enddo
        enddo
      endif
    end do
    write (jfile,fmt='(4hscf=,i4,3x,3hdp=,d20.10)',err=9999) iscf,diff_charge
    if (ksconv) then
      write (jfile,fmt='(a,i4,a)') 'max # of iterations (KS equation) =',ksitmax,' , all states converged'
    else
      write (jfile,fmt='(a,i4,a)') 'max # of iterations (KS equation) =',ksitmax,' , not all states converged'
    endif
    if ((nso==0).or.(ncol==2)) then
      write (jfile,fmt='(a,e20.10)') 'max (residual_norm * occupation) =',wsqresmx
    else
      write (jfile,fmt='(a,e20.10)') 'max residual_norm =',wsqresmx
    endif
    write (jfile,fmt='(a,e20.10)') 'smallest eigenvalue =',svalmin 
  enddo ! loop
  write (ndisp,*,err=9999) '------------------------------------------------------------------------'

return
  9999 call stopp('output_countdat: error writing output')
end subroutine output_countdat


subroutine output_energy( &
 chdir,                                                                                        & ! <
 ndisp,njfile,iscf,nperi,npopshow,nso,nspv,nums,ncol,natom,numkmx,neigmx,nwskp,nwskptot,kband, & ! <
 tnumele,ferm,spinpol,atmpole,diff_charge,residual_states_wfc,sval_wfc,fnele_wfc,sconst,eps,   & ! <
 eneel,eneex,eneha,eneat,eneco,eneth,eneof,enejel,enebc,enespr)                                  ! <
use mod_stopp
use var_keys, only: key_band_no
implicit none
character, intent(in) :: chdir*200
integer,   intent(in) :: ndisp,njfile,iscf,nperi,npopshow,nso,nspv,nums,ncol,natom,numkmx,neigmx
integer,   intent(in) :: nwskp(numkmx),nwskptot,kband
real*8,    intent(in) :: tnumele,ferm,sconst,eps
real*8,    intent(in) :: spinpol(max(1,nspv-1),0:natom*min(1,nspv-1))
real*8,    intent(in) :: atmpole(10,natom)
real*8,    intent(in) :: diff_charge
real*8,    intent(in) :: residual_states_wfc(neigmx,nums+1-ncol,numkmx)
real*8,    intent(in) :: sval_wfc(           neigmx,nums+1-ncol,numkmx)
real*8,    intent(in) :: fnele_wfc(          neigmx,nums+1-ncol,numkmx)
real*8,    intent(in) :: eneel,eneex,eneha,eneat,eneco,eneth,eneof,enejel,enebc,enespr

integer :: loop,jfile,ns,na
character :: fname*200

  do loop=1,2
    if (loop==1) jfile=ndisp
    if (loop==2) jfile=njfile

    write (jfile,*,err=9999) 'total charge',tnumele,'(electrons)'
    write (jfile,*,err=9999) 'fermi level',ferm,'(hartree)'
    if (nspv>1) then
      if (nspv==2) then
        write (jfile,*,err=9999) 'spin polarization',spinpol(1,0),'(\mu_B)'
      else
        write (jfile,fmt='(A,3(1x,e11.4))',err=9999) '(z,x,y)-polarization (in \mu_B): ', &
         (spinpol(ns,0),ns=1,3)
      endif
    end if
    if ((nperi .lt. 3) .or. (npopshow .eq. 1)) then
      write (jfile,*,err=9999) 'atomic charge','(electrons)'
      do na=1,natom
        write (jfile,*,err=9999) na,atmpole(1,na)
      end do
    end if
    call output_evlist( &
     jfile,nso,nums,ncol,numkmx,neigmx,                     & ! <
     nwskp,nwskptot,residual_states_wfc,sval_wfc,fnele_wfc)   ! <
    write (jfile,fmt='(4hScf=,i4,3x,3hdp=,d20.10)',err=9999) iscf,diff_charge
    write (jfile,*,err=9999) '   one electron energy     ',eneel,'(hartree)'
    write (jfile,*,err=9999) 'exchange correlation energy',eneex,'(hartree)'
    write (jfile,*,err=9999) '   hartree       energy    ',eneha,'(hartree)'
    write (jfile,*,err=9999) '   ewald         energy    ',eneat,'(hartree)'
    write (jfile,*,err=9999) '   field correction        ',eneco,'(hartree)'
    write (jfile,*,err=9999) '   B_con correction        ',enebc,'(hartree)'
    if (nperi .eq. 3) then
      write (jfile,*,err=9999) ' energy offset of ps. pot. ',eneof,'(hartree)'
      write (jfile,*,err=9999) ' energy jellium/ions       ',enejel,'(hartree)'
    end if
    write (jfile,*,err=9999) '  helmholtz free energy    ',eneel+eneex+eneha+eneat+eneco+eneth+eneof+enejel+enebc,'(hartree)'
    if (sconst>eps) write (jfile,*,err=9999) 'helmholtz free energy < spr',eneel+eneex+eneha+eneat+eneco+eneth+eneof+enejel+enebc+enespr,'(hartree)'
    write (jfile,*,err=9999) '   Total         energy    ',eneel+eneex+eneha+eneat+eneco+eneof+enejel+enebc,'(hartree)'
    if (sconst>eps) write (jfile,*,err=9999) '   Total    energy < spr   ',eneel+eneex+eneha+eneat+eneco+eneof+enejel+enebc+enespr,'(hartree)'
  enddo ! loop
  if (kband==key_band_no) then
    fname='fermilevel.dat'
    if (len_trim(chdir) > 0) fname=trim(chdir)//'/'//fname
    open(101,file=fname,form='formatted')
    write(101,*) ferm
    close(101)
  end if

return
  9999 call stopp('output_energy: error writing output')
end subroutine output_energy


subroutine output_bands( & 
 chdir,                            & ! <
 nso,nums,ncol,neigmx,numk,numkmx, & ! <
 sval,svalsoc)                       ! < 
use mod_mpi
use mod_stopp
implicit none
character, intent(in) :: chdir*200
integer, intent(in) :: nso,nums,ncol,neigmx,numk,numkmx
real*8,  intent(in) :: sval(neigmx,nums+1-ncol,numk),svalsoc(2*neigmx*nso-nso+1,1,numk*nso-nso+1)
logical  :: isfile,lopen
integer  :: nfile,irank,nso1,nums1,neigmx1,ns,l,nk,nk2 
character:: chns*3,chneigmx*3,fname*200
real*8 :: ferm
real*8,allocatable:: svalall(:,:,:) 

  if (myr_space==0) then

  ferm=0.0d0
  fname='fermilevel.dat'
  if (len_trim(chdir) > 0) fname=trim(chdir)//'/'//fname
  inquire(file=fname,exist=lopen)
  if (lopen) then
    open(101,file=fname,form='formatted',err=9999)
    read(101,*) ferm
    close(101)
  end if
    
  do nso1= 0,min(nso,2-ncol)

    if (myr_kpt/=0) then

      do nk= 1,numk
        if (nso1==0) then
          do ns= 1,nums+1-ncol
            call mpi_send(sval(1,ns,nk),neigmx,mpi_double_precision,0,0,mpicom_kpt,mpij)
          enddo
        else
          call mpi_send(svalsoc(1,1,nk),2*neigmx,mpi_double_precision,0,0,mpicom_kpt,mpij)
        endif
      enddo

    else 

      if (nso1==0) then
        neigmx1= neigmx
        nums1= nums+1-ncol 
      else
        neigmx1= 2*neigmx
        nums1= 1
      endif  

      allocate(svalall(neigmx1,nums1,numkmx))

      do nk= 1,numk
        if (nso1==0) then
          do ns= 1,nums1
            svalall(1:neigmx,ns,nk)= sval(1:neigmx,ns,nk)
          enddo
        else
          svalall(1:2*neigmx,1,nk)= svalsoc(1:2*neigmx,1,nk)
        endif
      enddo  

      nk= 1
      do irank= 1,nprock-1
        nk= nk + numkproc(irank-1)
        do nk2= nk,nk-1+max(numkproc(irank),1)
          if (nso1==0) then
            do ns= 1,nums1 
              call mpi_recv(svalall(1,ns,nk2),neigmx,mpi_double_precision,irank,0,mpicom_kpt,mpistat,mpij)
            enddo
          else
            call mpi_recv(svalall(1,1,nk2),2*neigmx,mpi_double_precision,irank,0,mpicom_kpt,mpistat,mpij)
          endif
        enddo
      enddo

      isfile= .true.
      nfile = 10
      do while (isfile)
        nfile= nfile+1
        inquire(unit=nfile,opened=isfile)
      enddo
      if (neigmx1>999) call stopp('output_bands: too many bands') 
      write(chneigmx,fmt='(i3)') neigmx1 
      do ns= 1,nums1 
        chns= ''
        if (nums1>1) write(chns,fmt='(i1)') ns
        if (nso1==1) chns= 'Soc'
        open(nfile,file='bands'//trim(chns)//'.dat',form='formatted',status='replace')
        do nk= 1,numkmx
          write(nfile,fmt='(i5,'//chneigmx//'(1x,e11.5))') nk,((svalall(l,ns,nk)-ferm)*27.212d0,l=1,neigmx1) 
        enddo
        close(nfile)
      enddo

      deallocate(svalall)

    endif ! (myr_kpt/=0) else 

  enddo ! nso1 

  endif ! (myr_space==0)

  return

9999 call stopp('output_bands: error writing bands')

end subroutine output_bands 


subroutine output_molecdyn( &
 chdir,                                                                                      & ! <
 ioput,ndisp,njfile,catmfn,natom,num_spe,lmx,indspe,numz,nmdx,nmdy,nmdz,natsoc,ntyppp,watom, & ! <
 xmax,ymax,zmax,imd,atx,aty,atz,fatx,faty,fatz)                                                ! <
use mod_stopp
implicit none
character, intent(in) :: chdir*200
integer,   intent(in) :: ioput,ndisp,njfile
integer,   intent(in) :: natom,num_spe,lmx 
integer,   intent(in) :: indspe(natom),numz(natom),nmdx(natom),nmdy(natom),nmdz(natom),natsoc(natom,0:lmx-1),ntyppp(num_spe)
real*8,    intent(in) :: watom(natom),xmax,ymax,zmax
integer,   intent(in) :: imd
real*8,    intent(in) :: atx(natom),aty(natom),atz(natom),fatx(natom),faty(natom),fatz(natom)
character,  intent(in):: catmfn*50
integer  :: jfile,loop,loopmx
integer  :: na,lmn
character:: chlmx,fname*200

  lmn= min(1,lmx-1)
  write(chlmx,fmt='(i1)') lmx-lmn  

  loopmx= 2
  if (ioput==3) loopmx= 1

  do loop=1,loopmx
    if (loop==1) jfile= njfile
    if (loop==2) jfile= ndisp

    select case (ioput)
    case (1)
      write (jfile,*,err=9999) 'mdstep=',imd
      write (jfile,*,err=9999) 'atom coordinate (input)'
    case (2)
      write (jfile,*,err=9999) 'atom coordinate (output)'
    case (3)
      fname=catmfn
      if (len_trim(chdir) > 0) fname=trim(chdir)//'/'//fname
      open (jfile,file=fname)
      write(jfile,fmt='(a)') '! [x], [y], [z], [atom number], switch [x], [y], [z], [weight], switches [soc], [pp], [na]'
    end select
    do na=1,natom
      write (jfile,fmt='(3d25.16,4i4,f10.2,3x,'//chlmx//'i1,i3,1x,i4,"a")',err=9999) &
       atx(na),aty(na),atz(na),numz(na),nmdx(na),nmdy(na),nmdz(na),watom(na),natsoc(na,lmn:lmx-1),ntyppp(indspe(na)),na
    enddo
    select case (ioput)
    case (2)
      ! call output_force(-1,jfile,natom,fx,fy,fz)
    case (3)
      write(jfile,fmt='(1x)')
      write(jfile,fmt='(a)') '! Bravais matrix (output only!)'
      write(jfile,fmt='(3d25.16)') 2.0d0*xmax,0.0d0,0.0d0
      write(jfile,fmt='(3d25.16)') 0.0d0,2.0d0*ymax,0.0d0
      write(jfile,fmt='(3d25.16)') 0.0d0,0.0d0,2.0d0*zmax
      close(jfile)
    end select
  enddo ! loop

return
  9999 call stopp('output_molecdyn: error writing output')
end subroutine output_molecdyn


subroutine output_force( &
 ndisp,njfile,natom, & ! <
 fatx,faty,fatz)       ! <
use mod_stopp
implicit none
integer, intent(in) :: ndisp,njfile
integer, intent(in) :: natom
real*8,  intent(in) :: fatx(natom),faty(natom),fatz(natom)
integer :: loop,jfile,na

  do loop=1,2
    if (loop==1) jfile=ndisp
    if (loop==2) jfile=njfile
    if (jfile>0) then

      write(jfile,*,err=9999) 'atomic force'
      do na=1,natom
        write(jfile,fmt='(i5,3d25.16)',err=9999) na,fatx(na),faty(na),fatz(na)
      enddo

    endif
  enddo ! loop

return
  9999 call stopp('output_force: error writing output')
end subroutine output_force


subroutine output_polcon( &
 ioput,ndisp,njfile,                                                                                & ! <
 key_polcon2_new,key_polcon2_rot,key_polcon2_none,key_polcon2_dir,key_polcon2_size,key_polcon2_fix, & ! <
 natom,npolcon2,polconb,polconmag,polconeta)                                                          ! <
use mod_stopp
implicit none
integer, intent(in) :: ioput,ndisp,njfile
integer, intent(in) :: key_polcon2_new,key_polcon2_rot,key_polcon2_none,key_polcon2_dir,key_polcon2_size,key_polcon2_fix
integer, intent(in) :: natom,npolcon2(0:natom)
real*8,  intent(in) :: polconb(3,natom),polconmag(3,natom),polconeta(2,0:natom)
integer :: jfile,loop,loopmx
integer :: na

  loopmx= 2
  if (ioput==3) loopmx= 1

  do loop=1,loopmx
    if (loop==1) jfile= njfile
    if (loop==2) jfile= ndisp

    select case (ioput)
    case (1)
      write(jfile,fmt='(1x)',err=9999)
      write (jfile,*,err=9999) 'B_con (input):'
      write(jfile,fmt='(1x)',err=9999)
    case (2)
      write(jfile,fmt='(1x)',err=9999)
      write (jfile,*,err=9999) 'B_con (output):'
      write(jfile,fmt='(1x)',err=9999)
    case (3)
      open (jfile,file='atom.mag',form='formatted',status='replace')
    end select
    write(jfile,fmt='(a)',err=9999) ' eta_perp    eta_para    switch0'
    write(jfile,fmt='(2(e10.3,2x),3x,i2)',err=9999) polconeta(:,0), npolcon2(0)
    write(jfile,fmt='(1x)',err=9999)
    write(jfile,fmt='(a)',err=9999) &
     '    B_z            B_x            B_y       switchN     M_z          M_x          M_y         rel._perp   rel._para     na'
    do na= 1,natom
      write(jfile,fmt='(3(e13.6,2x),2x,i1,4x,3(e11.4,2x),2(2x,e10.3),3x,i4,"a")',err=9999) &
       polconb(:,na), npolcon2(na), polconmag(:,na), polconeta(:,na), na
    enddo
    write(jfile,fmt='(1x)',err=9999)

    if (ioput==3) then
      write(jfile,fmt='(1x)',err=9999)
      write(jfile,fmt='("switch0:  (",i1,") update B_con")',err=9999) key_polcon2_new
      write(jfile,fmt='("          (",i1,") update B_con, rotate mag. moments before 1st iteration")',err=9999) key_polcon2_rot
      write(jfile,fmt='("          (",i1,") keep B_con fixed")',err=9999) key_polcon2_fix
      write(jfile,fmt='("switchN:  (",i1,") don''t constrain")',err=9999) key_polcon2_none
      write(jfile,fmt='("          (",i1,") constrain mag. direction")',err=9999) key_polcon2_dir
      write(jfile,fmt='("          (",i1,") constrain |m| and mag. direction")',err=9999) key_polcon2_size
      write(jfile,fmt='("          (",i1,") keep B_con fixed")',err=9999) key_polcon2_fix
      close(jfile)
    endif

  enddo

return
  9999 call stopp('output_polcon: error writing output')
end subroutine output_polcon


subroutine output_xcrysden( &
 chdir,natom,ndisp,njfile,atx,aty,atz,fatx,faty,fatz,numz) ! <
use mod_stopp
implicit none
character, intent(in) :: chdir*200
integer,   intent(in) :: natom,ndisp,njfile
real*8,    intent(in) :: atx(natom),aty(natom),atz(natom)
real*8,    intent(in) :: fatx(natom),faty(natom),fatz(natom)
integer,   intent(in) :: numz(natom)
integer :: na
character fname*200
  fname='xcrysden.xsf'
  if (len_trim(chdir) > 0) fname=trim(chdir)//'/'//fname
  open(njfile,file=fname,err=9999)
  write(njfile,*,err=9999) 'ATOMS'
  do na=1,natom
    write(njfile,'(i3,3x,6f13.8)',err=9999) numz(na),atx(na)*0.52918d0,aty(na)*0.52918d0,atz(na)*0.52918d0,fatx(na),faty(na),fatz(na)
  end do
  close(njfile,err=9999)
return
  9999 call stopp('output_xcrysden: error writing output')
end subroutine


subroutine output_evlist( &
 jfile,nso,nums,ncol,numkmx,neigmx,                     & ! <
 nwskp,nwskptot,residual_states_wfc,sval_wfc,fnele_wfc)   ! <
use mod_stopp
implicit none
integer, intent(in):: jfile,nso,nums,ncol,numkmx,neigmx
integer, intent(in):: nwskp(numkmx),nwskptot
real*8,  intent(in):: residual_states_wfc(neigmx,nums+1-ncol,numkmx)
real*8,  intent(in):: sval_wfc( neigmx*(1+nso*(2-ncol)),1+(nums-1)*(2-ncol)*(1-nso),numkmx)
real*8,  intent(in):: fnele_wfc(neigmx*(1+nso*(2-ncol)),1+(nums-1)*(2-ncol)*(1-nso),numkmx)
integer:: nk,ns,l,l2
real*8 :: fwghtinv 

  write(jfile,fmt='(1x)',err=9999)
  if ((nso==0).or.(ncol==2)) then
    write(jfile,fmt='(a)') '   k,spin,  band ,     eigen value    ,    residual norm  ,    occupation'  
  else
    write(jfile,fmt='(a)') '   k , band ,   eigen value    ,   occupation     | spin,band,   residual norm'
  endif 
  do nk=1,numkmx
    if (nwskp(nk)/=0) then
      fwghtinv= dble(nwskptot*nums)/dble(nwskp(nk)*2)
    else
      fwghtinv= 0.0d0
    endif
    if ((nso==0).or.(ncol==2)) then
      do ns=1,nums+1-ncol
      do l=1,neigmx
        write (jfile,fmt='(i4,1x,i4,1x,i6,1x,3e20.10)',err=9999) &
         nk,ns,l,sval_wfc(l,ns,nk),dsqrt(residual_states_wfc(l,ns,nk)),fnele_wfc(l,ns,nk)*fwghtinv
      end do
      end do
    else
      do ns=1,2
        do l2=1,neigmx
          l= (ns-1)*neigmx+l2
          if (ns<=nums) then
            write (jfile,fmt='(i4,1x,i6,2e19.10," | ",i1,1x,i6,e19.10)',err=9999) &
             nk,l,sval_wfc(l,1,nk),fnele_wfc(l,1,nk)*fwghtinv, ns,l2,dsqrt(residual_states_wfc(l2,ns,nk))
          else
            write (jfile,fmt='(i4,1x,i6,2e19.10," | ")',err=9999) &
             nk,l,sval_wfc(l,1,nk),fnele_wfc(l,1,nk)*fwghtinv
          endif
        enddo
      enddo
    endif
  end do

return
  9999 call stopp('output_evlist: error writing output')
end subroutine output_evlist 


end module mod_output

