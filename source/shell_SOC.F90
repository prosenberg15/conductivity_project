implicit none

integer :: ndim, nk, nkx, nky, nkz, nkp, ikx, iky, ikz, ik
integer :: l, ir, i, j, inshell, nk_sim, nptcl
real(kind=8), allocatable :: ksq(:)
real(kind=8) :: lambda, a, k, scaleBZ, pi
character(len=1) :: ndim_c
character(len=20) :: lambda_c, nk_c, nptcl_c

pi = acos(-1.d0)

call get_command_argument(1,ndim_c)
call get_command_argument(2,nk_c)
call get_command_argument(3,nptcl_c)
call get_command_argument(4,lambda_c)

read(ndim_c,*) ndim
read(lambda_c,*) lambda
read(nk_c,*) nk_sim
read(nptcl_c,*) nptcl

!lambda = lambda * sqrt(2*pi*nptcl)/nk_sim
scaleBZ = (2*pi/nk_sim)**2

nk  = 10
nkx = nk
nky = nk
if (ndim.eq.3) then
   nkz = nk
else
   nkz = 0
endif

nkp = 0
do ikz = -nkz,nkz
   do iky = -nky,nky
      do ikx = -nkx,nkx
         k=ikz**2+iky**2+ikx**2
         if (k .le. nk**2) nkp = nkp+1
      enddo
   enddo
enddo

nkp = 2*nkp
allocate(ksq(nkp))
ik = 0
do ikz = -nkz,nkz
   do iky = -nky,nky
      do ikx = -nkx,nkx
         k=ikz**2+iky**2+ikx**2
         if (k .le. nk**2) then
            k=scaleBZ*k
            ik = ik + 1
            ksq(ik) = k - lambda *sqrt(k)
            ik = ik + 1
            ksq(ik) = k + lambda *sqrt(k)
         endif
      enddo
   enddo
enddo

! Sort with heapsort
l  = nkp/2+1
ir = nkp
do
   if (l .gt. 1) then
      l  = l-1
      a  = ksq(l)
   else
      a = ksq(ir)
      ksq(ir) = ksq(1)
      ir = ir-1
      if (ir .eq. 1) then
         ksq(1) = a
         exit
      end if
   end if
   i = l
   j = l+l
   if (j .le. ir) then
      do
         if (j .lt. ir) then
            if (ksq(j) .lt. ksq(j+1)) j = j+1
         endif
         if (a .lt. ksq(j)) then
            ksq(i) = ksq(j)
            i = j
            j = j+j
         else
            j=ir+1
         end if
         if (j .gt. ir) exit
      enddo
   endif
   ksq(i) = a
enddo


inshell=1
do ik = 2, nkp
   if (ksq(ik) .eq. ksq(ik-1)) then
      inshell = inshell + 1
   else
      if (ik-1 .ge. nptcl)then
         write(*,*) inshell, (ik-1), ksq(ik-1)
         exit
      endif
      inshell = 1
   endif
enddo

end
