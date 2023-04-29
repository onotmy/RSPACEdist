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
! **********  scf_occupation8f.F90 06/20/2018-01  **********

module mod_scf_occupation
implicit none
contains

subroutine scf_occupation( &
 nd,mpre,npolcon,ndisp,nso1,nso,nums,ncol,neigmx,numkmx,numk,nwskp,nwskptot, & ! <
 key_polcon_occ,                                                             & ! <
 tf,tfmax,tfmin,tnumele,polconocc,sval,residual_states,                      & ! <
 ksconv,ksitmax,                                                             & ! X
 ferm,sval_wfc,residual_states_wfc,fnele_wfc,fnele)                            ! >
use mod_mpi
use mod_mathfunctions, only: fermidis
implicit none
real*8, parameter :: eps_occ=1.0d-8
integer, intent(in)   :: nd,mpre,npolcon,ndisp
integer, intent(in)   :: nso1,nso,nums,ncol,neigmx,numkmx,numk
integer, intent(in)   :: nwskp(numkmx),nwskptot
integer, intent(in)   :: key_polcon_occ
real*8,  intent(in)   :: tf,tfmax,tfmin,tnumele,polconocc
real*8,  intent(in)   :: sval( neigmx*(1+nso1*nso*(2-ncol)),1+(nums-1)*(2-ncol)*(1-nso1*nso),numk)
real*8,  intent(in)   :: residual_states(neigmx,nums+1-ncol,numk)
logical, intent(inout):: ksconv
integer, intent(inout):: ksitmax
real*8,  intent(out)  :: ferm
real*8,  intent(out)  :: sval_wfc( neigmx*(1+nso*(2-ncol))*nd-nd+1,(1+(nums-1)*(2-ncol)*(1-nso))*nd-nd+1,numkmx*nd-nd+1)
real*8,  intent(out)  :: residual_states_wfc(neigmx*nd-nd+1,(nums+1-ncol)*nd-nd+1,numkmx*nd-nd+1)
real*8,  intent(out)  :: fnele_wfc(neigmx*(1+nso*(2-ncol))*nd-nd+1,(1+(nums-1)*(2-ncol)*(1-nso))*nd-nd+1,numkmx*nd-nd+1)
real*8,  intent(out)  :: fnele(neigmx*(1+nso1*nso*(2-ncol)),1+(nums-1)*(2-ncol)*(1-nso1*nso),numk)
logical :: l_help
integer :: neigmxdim,neigmx2,nums2,ncol2 
integer :: ns,irank,l,nk,nk2,i_help

  if (myr_space==0) then

    if (myr_kpt==0) then
      nk= 1
      do irank= 1,nprock-1
        nk= nk + numkproc(irank-1)
        do nk2= nk,nk-1+max(numkproc(irank),1)
          do ns= 1,nums+1-ncol
            call mpi_recv(residual_states_wfc(1,ns,nk2),neigmx,mpi_double_precision,irank,0,mpicom_kpt,mpistat,mpij)
          end do
          if ((nso==0).or.(ncol==2)) then 
            do ns= 1,nums+1-ncol
              call mpi_recv(sval_wfc(1,ns,nk2),neigmx,mpi_double_precision,irank,0,mpicom_kpt,mpistat,mpij)
            end do
          else
            if (nso1==0) then 
              do ns= 1,nums
                call mpi_recv(sval_wfc(1+(ns-1)*neigmx,1,nk2),neigmx,mpi_double_precision,irank,0,mpicom_kpt,mpistat,mpij)
              end do
            else 
              call mpi_recv(sval_wfc(1,1,nk2),2*neigmx,mpi_double_precision,irank,0,mpicom_kpt,mpistat,mpij)
            end if 
          end if 
        end do
      end do
    else
      do nk= 1,numk
        do ns= 1,nums+1-ncol
          call mpi_send(residual_states(1,ns,nk),neigmx,mpi_double_precision,0,0,mpicom_kpt,mpij)
        end do
        if ((nso==0).or.(ncol==2)) then
          do ns= 1,nums+1-ncol
            call mpi_send(sval(1,ns,nk),neigmx,mpi_double_precision,0,0,mpicom_kpt,mpij)
          end do
        else
          if (nso1==0) then 
            do ns= 1,nums
              call mpi_send(sval(1,ns,nk),neigmx,mpi_double_precision,0,0,mpicom_kpt,mpij)
            end do
          else 
            call mpi_send(sval(1,1,nk),2*neigmx,mpi_double_precision,0,0,mpicom_kpt,mpij)
          end if 
        end if 
      end do
    end if

    l_help= ksconv
    call mpi_reduce(l_help,ksconv,1,mpi_logical,mpi_land,0,mpicom_kpt,mpij)
    i_help= ksitmax
    call mpi_reduce(i_help,ksitmax,1,mpi_integer,mpi_max,0,mpicom_kpt,mpij)

  end if ! (myr_space==0)

  if (myrank_glbl==0) then

    do nk= 1,numk
      do ns= 1,nums+1-ncol
        do l=1,neigmx
          residual_states_wfc(l,ns,nk)= residual_states(l,ns,nk)
        end do
      end do
      if ((nso==0).or.(ncol==2)) then 
        do ns= 1,nums+1-ncol
          do l=1,neigmx
            sval_wfc(l,ns,nk)= sval(l,ns,nk)
          end do
        end do
      else
        if (nso1==0) then 
          do ns= 1,nums
            do l=1,neigmx
              sval_wfc(l+(ns-1)*neigmx,1,nk)= sval(l,ns,nk)
            end do
          end do 
        else 
          do l=1,2*neigmx
            sval_wfc(l,1,nk)= sval(l,1,nk)
          end do
        end if
      end if
    end do 

    if ((nso==0).or.(ncol==2)) then
      neigmxdim= neigmx
      nums2  = nums 
      ncol2  = ncol 
      neigmx2  = neigmx
    else
      neigmxdim= 2*neigmx
      if (nso1==0) then 
        nums2= nums  
        ncol2= nums
        neigmx2= nums*neigmx
      else 
        nums2= 2 
        ncol2= 2
        neigmx2= 2*neigmx
      end if
    end if   

    call scf_occupation_fnele( &
     mpre,npolcon,ndisp,nums2,ncol2,neigmxdim,neigmx2,numkmx,nwskp,nwskptot, & ! <
     key_polcon_occ,                                                     & ! <
     tf,tfmax,tfmin,tnumele,polconocc,sval_wfc,                          & ! <
     ferm,fnele_wfc)                                                       ! >

    do nk= 1,numk
      if ((nso==0).or.(ncol==2)) then
        do ns= 1,nums+1-ncol
          do l=1,neigmx
            fnele(l,ns,nk)= fnele_wfc(l,ns,nk)
          end do
        end do
      else
        if (nso1==0) then
          do ns= 1,nums
            do l=1,neigmx
              fnele(l,ns,nk)= fnele_wfc(l+(ns-1)*neigmx,1,nk)
            end do
          end do
        else
          do l=1,2*neigmx
            fnele(l,1,nk)= fnele_wfc(l,1,nk)
          end do
        end if
      end if
    end do

  else

    ferm= 0.0d0
    sval_wfc(:,:,:)= 0.0d0
    residual_states_wfc(:,:,:)= 0.0d0
    fnele_wfc(:,:,:)= 0.0d0

  end if ! (myrank_glbl==0) else

  if (myr_space==0) then
    if (myr_kpt==0) then

      nk= 1
      do irank= 1,nprock-1
        nk= nk + numkproc(irank-1)
        do nk2= nk,nk-1+max(numkproc(irank),1)
          if ((nso==0).or.(ncol==2)) then 
            do ns= 1,nums+1-ncol
              call mpi_send(fnele_wfc(1,ns,nk2),neigmx,mpi_double_precision,irank,0,mpicom_kpt,mpij)
            end do
          else
            if (nso1==0) then
              do ns= 1,nums
                call mpi_send(fnele_wfc(1+(ns-1)*neigmx,1,nk2),neigmx,mpi_double_precision,irank,0,mpicom_kpt,mpij)
              end do
            else
              call mpi_send(fnele_wfc(1,1,nk2),2*neigmx,mpi_double_precision,irank,0,mpicom_kpt,mpij)
            end if
          end if
        end do
      end do
    else

      do nk= 1,numk
        if ((nso==0).or.(ncol==2)) then
          do ns= 1,nums+1-ncol
            call mpi_recv(fnele(1,ns,nk),neigmx,mpi_double_precision,0,0,mpicom_kpt,mpistat,mpij)
          enddo
        else
          if (nso1==0) then
            do ns= 1,nums
              call mpi_recv(fnele(1,ns,nk),neigmx,mpi_double_precision,0,0,mpicom_kpt,mpistat,mpij)
            enddo
          else
            call mpi_recv(fnele(1,1,nk),2*neigmx,mpi_double_precision,0,0,mpicom_kpt,mpistat,mpij)
          end if
        end if
      end do

    end if
  end if
  if ((nso==0).or.(ncol==2)) then 
    call mpi_bcast(fnele,neigmx*(nums+1-ncol)*numk,mpi_double_precision,0,mpicom_space,mpij)
  else
    if (nso1==0) then
      call mpi_bcast(fnele,neigmx*nums*numk,mpi_double_precision,0,mpicom_space,mpij)
    else
      call mpi_bcast(fnele,2*neigmx*numk,mpi_double_precision,0,mpicom_space,mpij)
    endif
  end if

end subroutine scf_occupation


subroutine scf_occupation_fnele( &
 mpre,npolcon,ndisp,nums,ncol,neigmxdim,neigmx,numkmx,nwskp,nwskptot, & ! <
 key_polcon_occ,                                                    & ! <
 tf,tfmax,tfmin,tnumele,polconocc,sval_wfc,                         & ! <
 ferm,fnele_wfc)                                                      ! >
use mod_stopp
use mod_mathfunctions, only: fermidis
implicit none
real*8, parameter :: eps_occ=1.0d-8
integer, intent(in)   :: mpre,npolcon,ndisp
integer, intent(in)   :: nums,ncol,neigmxdim,neigmx,numkmx
integer, intent(in)   :: nwskp(numkmx),nwskptot
integer, intent(in)   :: key_polcon_occ
real*8,  intent(in)   :: tf,tfmax,tfmin,tnumele,polconocc
real*8,  intent(in)   :: sval_wfc( neigmxdim,nums+1-ncol,numkmx)
real*8,  intent(out)  :: ferm
real*8,  intent(out)  :: fnele_wfc(neigmxdim,nums+1-ncol,numkmx)
integer :: ns,ns1,ns2,nsp,nsp1,nsp2,l,nk,nk2,lbisec
real*8  :: fwght,ferm1,ferm3,fermp,bisec1,bisec3,bisec,tnumele0,tnumele1,fnelehigh

  if (mpre/=0) then

    do nk=1,numkmx
      fwght= dble(nwskp(nk)*2)/dble(nwskptot*nums)
      do ns=1,nums+1-ncol
        do l=1,neigmx
          fnele_wfc(l,ns,nk)=0.0d0
        end do
        do l=1,(ncol*nint(tnumele))/2
          fnele_wfc(l,ns,nk)=fwght
        end do
        if (mod(ncol*nint(tnumele),2)==1) fnele_wfc(nint(tnumele)/2+1,ns,nk)= fwght/2.0d0
      end do
    end do

  else

    nsp1= 1
    nsp2= 1
    ns1 = 1
    ns2 = nums+1-ncol
    tnumele1= tnumele

    if (npolcon==key_polcon_occ) then
      if ((nums/=2).or.(ncol/=1)) call stopp( &
       'scf_occupation_fnele: only in the col.mag. non-soc case, occ.numbers can be modified to constr. the spin pol.')
      nsp1= 1
      nsp2= 2
    end if

    do nsp= nsp1,nsp2

      if (npolcon==key_polcon_occ) then
        ns1= nsp
        ns2= nsp
        tnumele1= tnumele/2.0d0 +polconocc*dble(3-2*nsp)/2.0d0
      endif

      ferm1=tfmax
      ferm3=tfmin
      ferm=ferm1
      tnumele0=0.0d0

      do nk=1,numkmx
        fwght= dble(nwskp(nk)*2)/dble(nwskptot*nums)
        do ns=ns1,ns2
          do l=1,neigmx
            call fermidis(fwght,sval_wfc(l,ns,nk),tf,ferm, fnele_wfc(l,ns,nk))
            tnumele0=tnumele0+fnele_wfc(l,ns,nk)
          end do
        end do
      end do
      bisec1=tnumele0-tnumele1
      ferm=ferm3
      tnumele0=0.0d0
      do nk=1,numkmx
        fwght= dble(nwskp(nk)*2)/dble(nwskptot*nums)
        do ns=ns1,ns2
          do l=1,neigmx
            call fermidis(fwght,sval_wfc(l,ns,nk),tf,ferm, fnele_wfc(l,ns,nk))
            tnumele0=tnumele0+fnele_wfc(l,ns,nk)
          end do
        end do
      end do
      bisec3=tnumele0-tnumele1
      lbisec=0
      do while (lbisec < 100)
        lbisec=lbisec+1
        ferm=(ferm1+ferm3)*0.5d0
        tnumele0=0.0d0
        do nk=1,numkmx
          fwght= dble(nwskp(nk)*2)/dble(nwskptot*nums)
          do ns=ns1,ns2
            do l=1,neigmx
              call fermidis(fwght,sval_wfc(l,ns,nk),tf,ferm, fnele_wfc(l,ns,nk))
              tnumele0=tnumele0+fnele_wfc(l,ns,nk)
            end do
          end do
        end do
        bisec=tnumele0-tnumele1
        if (dabs(bisec) .lt. 1.0d-12) lbisec=100
        if (bisec1*bisec .gt. 0.0d0) then
          ferm1=ferm
          bisec1=bisec
        end if
        if (bisec*bisec3 .gt. 0.0d0) then
          ferm3=ferm
          bisec3=bisec
        end if
      end do ! (lbisec < 100)

      if (npolcon==key_polcon_occ) then
        if (nsp==1) then
          fermp= ferm
        else
          write(ndisp,fmt='("Delta (Fermi energy) = ",e11.4," htr")') fermp-ferm
          ferm= (fermp+ferm)/2.0d0
        end if
      end if

    end do ! nsp

    fnelehigh= 0.0d0
    ns2= 1
    nk2= 1
    do nk= 1,numkmx
      if (nwskp(nk)>0) then
        fwght= dble(nwskptot*nums)/dble(nwskp(nk)*2)
      else
        fwght= 0.0d0
      endif
      do ns= 1,nums+1-ncol
        if (fnelehigh<fwght*fnele_wfc(neigmx,ns,nk)) then
          fnelehigh= fwght*fnele_wfc(neigmx,ns,nk)
          ns2= ns
          nk2= nk
        end if
      end do
    end do
    if (fnelehigh>eps_occ) then
      write(ndisp,fmt='("WARNING!! At (nk=",i4,", ns=",i2,") all states are at least ",e11.4," times occupied !")') &
       nk2,ns2,fnelehigh
      write(ndisp,fmt='("(maybe increase neigmx)")')
    end if

  end if ! (mpre/=0) else

end subroutine scf_occupation_fnele


end module
