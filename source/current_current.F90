subroutine current_current(nk,nsites,U,V,ener)

integer,intent(IN)::nk(3)
integer,intent(IN)::nsites
real(kind=8),intent(IN)::ener(2*nsites)
complex(kind=8),intent(IN)::U(2*nsites,2*nsites),V(2*nsites,2*nsites)

integer::i,j,k,ix,iy,iz,ir,sir,irs,is,sirp,isp
integer::plus(2*nsites)
real(kind=8)::r(3,2*nsites),eig(2*nsites)
real(kind=8)::q,omega,pi
complex(kind=8),parameter::Zero=cmplx(0.d0,0.d0)
complex(kind=8),parameter::One =cmplx(1.d0,0.d0)
complex(kind=8),parameter::Xi  =cmplx(0.d0,1.d0)
complex(kind=8)::Lambdar(nsites),Ekin,Lambdaq
complex(kind=8)::Ar,Br,Arp,Brp,E,Om

!print U and V
open(2,file='U_matrix',status='unknown')
do i=1,2*nsites
  do j=1,2*nsites
    write(2,*)j,i,dble(U(j,i)),aimag(U(j,i)),ener(i)
  enddo
enddo
close(2)

open(2,file='V_matrix',status='unknown')
do i=1,2*nsites
  do j=1,2*nsites
    write(2,*)j,i,dble(V(j,i)),aimag(V(j,i)),ener(i)
  enddo
enddo
close(2)



!assume to receive in input negative eigenvalues of HFB hamiltonian
eig(:)=-ener(:)

!build the correspondence site --> site + x, for current in x direction
k=0
do iz=1,nk(3)
  do iy=1,nk(2)
    do ix=1,nk(1)
      k=k+1
      r(1,k)=dble(ix-1)
      r(2,k)=dble(iy-1)
      r(3,k)=dble(iz-1)
      if(mod(k,nk(1)).ne.0)then
        plus(k)       =k+1
        plus(k+nsites)=k+1+nsites
      else
        plus(k)       =k-nk(1)+1
        plus(k+nsites)=k-nk(1)+1+nsites
      endif
    enddo
  enddo
enddo

!compute < - k_x >

Ekin=Zero

do ir=1,nsites
  do sir=0,1
    irs=ir+sir*nsites
    do i=1,2*nsites
      Ekin=Ekin+(conjg(V(plus(irs),i))*V(irs,i)+conjg(V(irs,i))*V(plus(irs),i))/dble(nsites)
    enddo
  enddo
enddo

open(2,file='x_kinetic',status='unknown')
write(2,*)dble(Ekin),aimag(Ekin)
close(2)



!compute retarded correlation Dret(J(r,t) J(0, 0))
Lambdar(:)=Zero
omega=0.d0

do ir=1,nsites
  do sir=0,1
    irs=ir+sir*nsites

    is=1
    do sirp=0,1
      isp=is+sirp*nsites

      do i=1,2*nsites
        do j=1,2*nsites

          Ar=conjg( U(plus(irs),i)*V(irs,j) - U(irs,i)*V(plus(irs),j) )
          Br=       V(plus(irs),i)*U(irs,j) - V(irs,i)*U(plus(irs),j)

!          Arp=V(plus(isp),i)*U(isp,j)-V(isp,i)*U(plus(isp),j)-V(plus(isp),j)*U(isp,i)+V(isp,j)*U(plus(isp),i)
          Arp=V(plus(isp),j)*U(isp,i)-V(isp,j)*U(plus(isp),i)-V(plus(isp),i)*U(isp,j)+V(isp,i)*U(plus(isp),j)

!          Brp=conjg( U(plus(isp),i)*V(isp,j)-U(isp,i)*V(plus(isp),j)-U(plus(isp),j)*V(isp,i)+U(isp,j)*V(plus(isp),i) )
          Brp=conjg(U(plus(isp),j)*V(isp,i)-U(isp,j)*V(plus(isp),i)-U(plus(isp),i)*V(isp,j)+U(isp,i)*V(plus(isp),j) )
   
!DEBUG
!          write(*,*)
!          write(*,*)'ir,sir,is,sirp '
!          write(*,*)ir,sir,is,sirp
!          write(*,*)
!          write(*,*)'plus(irs) = ',plus(irs)
!          write(*,*)'plus(isp) = ',plus(isp) 
!          write(*,*)
!          write(*,*)'i,j'
!          write(*,*)i,j
!          write(*,*)'Energy = ',eig(i)+eig(j)
!          write(*,*)
!          write(*,*)'Ar = ',Ar
!          write(*,*)'Br = ',Br
!          write(*,*)'Arp = ',Arp
!          write(*,*)'Brp = ',Brp
!          write(*,*)

 
          E=cmplx(eig(i)+eig(j),0.d0)
          Om=cmplx(omega,0.d0)

         ! Lambdar(ir)=Lambdar(ir)+((Ar*Arp)/(Om+E))-((Br*Brp)/(Om-E))
          Lambdar(ir)=Lambdar(ir)-((Br*Brp)/(Om-E))
!          write(*,*)'Accumulation ',ir,Lambdar(ir)

        enddo
      enddo
    enddo
  enddo
enddo

open(2,file='lambda_q',status='unknown')
pi=dacos(-1.d0)
do i=1,nk(2)
  q=dble(i)*2.d0*pi/dble(nk(2))
  Lambdaq=Zero
  do ir=1,nsites
    Lambdaq=Lambdaq+Lambdar(ir)*exp(-Xi*q*r(2,ir)) !/dble(Nsites)
  enddo
  write(2,*)q,dble(Lambdaq),aimag(Lambdaq)
enddo
close(2)


end subroutine current_current
