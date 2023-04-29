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
! **********  atomicocc_pdos8f.F90 12/05/2022-01  **********

subroutine atomicocc_pdos( &
 chdir,                                                        & ! <
 key_natpri_in,key_pp_paw,nrc,natom,nradmx,num_spe,num_ppcell, & ! <
 nprmx,lmx,nprjmx,nums,ncol,nspv,numk,neigmx,                  & ! <
 natpri,naps,indspe,ntyppp,nprj,nradct,nlind,noind,            & ! <
 awf,pwf,dradial,akv,                                          & ! <
 rspsep,cspsep)                                                  ! <
use mod_mpi
implicit none
character,  intent(in)  :: chdir*200
integer,    intent(in)  :: key_natpri_in,key_pp_paw
integer,    intent(in)  :: nrc,natom,nradmx,num_spe,num_ppcell,nprmx,lmx,nprjmx,nums,ncol,nspv,numk,neigmx
integer,    intent(in)  :: natpri(natom),naps(natom),indspe(natom)
integer,    intent(in)  :: ntyppp(num_spe),nprj(num_spe),nradct(num_spe),nlind(nprjmx,num_spe),noind(nprjmx,num_spe)
real*8,     intent(in)  :: awf(nradmx,nprmx*lmx,num_spe),pwf(nradmx,nprmx*lmx,num_spe)
real*8,     intent(in)  :: rspsep(nprjmx*(1-nrc)+nrc,num_ppcell*(1-nrc)+nrc,neigmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk*(1-nrc)+nrc)
real*8,     intent(in)  :: dradial(nradmx,num_spe),akv(nprjmx,nprjmx,num_spe)
complex*16, intent(in)  :: cspsep(nprjmx*nrc-nrc+1,num_ppcell*nrc-nrc+1,neigmx*nrc-nrc+1,nums*nrc-nrc+1,numk*nrc-nrc+1)
real*8,     allocatable :: rspsep1(:,:,:,:,:)
real*8,     allocatable :: orb_nrm(:,:,:),orb_dos(:)
complex*16, allocatable :: cspsep1(:,:,:,:,:)
integer     :: na,iaps,nk,l,ns,i,ir,i1,i2,j,j1,j2,irankk,numk1,nnk
real*8     :: tmpr,sum_tru
character*6 :: corb(9)
character :: fname*200
  corb(1)='   s  '
  corb(2)='  px  '
  corb(3)='  py  '
  corb(4)='  pz  '
  corb(5)='dxx-yy'
  corb(6)=' dzx  '
  corb(7)='dzz-rr'
  corb(8)=' dyz  '
  corb(9)=' dxy  '

  if (myr_kpt == 0) then
    allocate(orb_nrm(nprjmx,nprjmx,num_spe),orb_dos(lmx*lmx))
    na=1
    do j=1,nprj(indspe(na))
      j1=nlind(j,indspe(na))
      j2=noind(j,indspe(na))
      do i=1,nprj(indspe(na))
        i1=nlind(i,indspe(na))
        i2=noind(i,indspe(na))
        if (i1==j1) then
          if (ntyppp(indspe(na)) .eq. key_pp_paw) then
            sum_tru=0.0d0
            do ir=2,nradct(indspe(na))
              sum_tru=sum_tru+(awf(ir,i2,indspe(na))*awf(ir,j2,indspe(na)))*dradial(ir,indspe(na))
            end do
             orb_nrm(i,j,indspe(na))=sum_tru
          else
             orb_nrm(i,j,indspe(na))=dabs(akv(i,j,indspe(na)))
          end if
        end if
      end do
    end do
  end if

  nnk=0
  do irankk=0,nprock-1
    if (myr_kpt == 0) then
      if (irankk == 0) then
        numk1=numk
      else
        call mpi_recv(numk1,1,mpi_integer,irankk,0,mpicom_kpt,mpistat,mpij)
      end if
      allocate(rspsep1(nprjmx*(1-nrc)+nrc,num_ppcell*(1-nrc)+nrc,neigmx*(1-nrc)+nrc,nums*(1-nrc)+nrc,numk1*(1-nrc)+nrc))
      allocate(cspsep1(nprjmx*nrc-nrc+1,num_ppcell*nrc-nrc+1,neigmx*nrc-nrc+1,nums*nrc-nrc+1,numk1*nrc-nrc+1))
      if (irankk == 0) then
        if (nrc==0) then
          rspsep1=rspsep
        else
          cspsep1=cspsep
        end if
      else
        if (nrc==0) then
          call mpi_recv(rspsep1,nprjmx*num_ppcell*neigmx*nums*numk1,mpi_double_precision,irankk,0,mpicom_kpt,mpistat,mpij)
        else
          call mpi_recv(cspsep1,nprjmx*num_ppcell*neigmx*nums*numk1,mpi_double_complex,irankk,0,mpicom_kpt,mpistat,mpij)
        end if
      end if

      na=1
      if (natpri(na)==key_natpri_in) then
        if (irankk == 0) then
          fname='occ_pdos.txt'
          if (len_trim(chdir) > 0) fname=trim(chdir)//'/'//fname
          open(101,file=fname)
        end if
        iaps= naps(na)
  
!        if (irankk == 0) write(101,*) '    k     s   ieig      index1    index2   sqrt(Re(c_1^+*c_2))'
        do nk=1,numk1
          nnk=nnk+1
          do l= 1,neigmx
            do ns= 1,max(nums,nspv)
              orb_dos=0.0d0
              do j=1,nprj(indspe(na))
                j1=nlind(j,indspe(na))
                j2=noind(j,indspe(na))
                do i=1,nprj(indspe(na))
                  i1=nlind(i,indspe(na))
                  i2=noind(i,indspe(na))
                  if (i1==j1) then
                    if (nrc==0) then
                      tmpr = rspsep1(i,iaps,l,ns,nk)*rspsep1(j,iaps,l,ns,nk)
                    else
                      tmpr = dreal(dconjg(cspsep1(i,iaps,l,ns,nk))*cspsep1(j,iaps,l,ns,nk))
                    end if
                    if (ntyppp(indspe(na)) .ne. key_pp_paw) then
                      orb_dos(i1) = orb_dos(i1) + tmpr*orb_nrm(i,j,indspe(na))
                    else
                      orb_dos(i1) = orb_dos(i1) + tmpr*orb_nrm(i,j,indspe(na))**2
                    end if
                  end if
                end do ! i
              end do ! j
              do i=1,nprj(indspe(na))
                i1=nlind(i,indspe(na))
                i2=noind(i,indspe(na))
                if (mod(i2+1,2)+1 ==1) write(101,'(3i6,a8,e25.15)') nnk,ns,l,corb(i1),orb_dos(i1)
              end do
            end do ! ns
          end do ! l
        end do ! nk
  
      if (irankk == nprock-1) close(101)
      end if

    else

      if (irankk == myr_kpt) then
        call mpi_send(numk,1,mpi_integer,0,0,mpicom_kpt,mpij)
        if (nrc==0) then
          call mpi_send(rspsep,nprjmx*num_ppcell*neigmx*nums*numk,mpi_double_precision,0,0,mpicom_kpt,mpij)
        else
          call mpi_send(cspsep,nprjmx*num_ppcell*neigmx*nums*numk,mpi_double_complex,0,0,mpicom_kpt,mpij)
        end if
      end if

    end if

    call mpi_barrier(mpi_comm_world,mpij)
    if (myr_kpt == 0) deallocate(rspsep1,cspsep1)
  end do

  if (myr_kpt == 0) deallocate(orb_nrm,orb_dos)

  return

end subroutine atomicocc_pdos
