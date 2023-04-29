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
! **********  mathfunctions8f.F90 06/15/2014-01  **********

module mod_mathfunctions
implicit none
contains


subroutine fermidis( &
 fwght,r,t,rc, & ! <
 fdis)           ! >
implicit none
real*8,intent(in) ::fwght,r,t,rc
real*8,intent(out)::fdis
real*8 acoef
  acoef=(r-rc)/t
  if (dabs(acoef) .lt. 1.0d2) then
    fdis=1.0d0/(dexp(acoef)+1.0d0)
  else
    if (acoef .gt. 0.0d0) fdis=0.0d0
    if (acoef .lt. 0.0d0) fdis=1.0d0
  end if
  fdis= fdis*fwght
end subroutine fermidis


function expint1(x,l0)
implicit none
real*8,  intent(in) :: x
integer, intent(inout) :: l0
integer, parameter :: imax=100
real*8,  parameter :: eps=1.0d-15,fpmin=1.0d-30,euler=0.5772156649015329d0
integer i
real*8 a,b,c,d,del,af,h,expint1

  if(x <= 1.0d0) then
    expint1=-log(x)-euler
    af=1.0d0
    do i=1,imax
      af=-af*x/i
      del=-af/i
      expint1=expint1+del
      if(abs(del) < abs(expint1)*eps) return
    end do
    l0=l0+1
    return
  else
    b=x+1.0d0
    c=1.0d0/fpmin
    d=1.0d0/b
    h=d
    do i=1,imax
      a=-i*i
      b=b+2.0d0
      d=1.0d0/(a*d+b)
      c=b+a/c
      del=c*d
      h=h*del
      if (abs(del-1.0d0) < eps) then
        expint1=h*exp(-x)
        return
      end if
    end do
    l0=l0+1
    return
  end if
  return
end function expint1


end module
