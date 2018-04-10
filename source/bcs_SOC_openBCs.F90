program bcs_openBCs

  use mpiwrap
  implicit none

  integer                        :: nk(3), nsites, nmesh, msite, getCN, BCs, IVext 
  integer                        :: Nthetax, Nthetay, thetax, thetay, IH0
  real(kind=8)                   :: thop, lambda, U, mu, hz, temp, alpha, beta, tol
  real(kind=8)                   :: theta_twist(4), kpoint_x, kpoint_y, kpoint_z
  complex(kind=8),allocatable    :: Hzero(:,:)
  complex(kind=8),allocatable    :: phi(:,:,:)           
  real(kind=8),allocatable       :: Vext(:)     
  logical                        :: openBCs, calc_CN
  character*75                   :: filename
  character*10                   :: Lxstr, Lystr, Ustr, mustr, hzstr

  integer                        :: i,j
  real(kind=8)                   :: x,y
  real(kind=8),parameter         :: pi=3.14159265358979323846264338


  open(10, file="hfb_params")
  read(10,*) nk(1:3)    !size of the lattice
  read(10,*) BCs        !0 for PBC/TBC, 1 for open boundary conditions
  read(10,*) kpoint_x
  read(10,*) kpoint_y
  read(10,*) kpoint_z
  read(10,*) lambda
  read(10,*) thop
  read(10,*) U
  read(10,*) mu
  read(10,*) hz
  read(10,*) temp      !temperature
  read(10,*) alpha     !newgap(j) = alpha*newgap(j) + beta*gap(j)
  read(10,*) beta
  read(10,*) IVext     !1 if an external field is present, 0 otherwise
  read(10,*) IH0       !1 if read hopping matrix from file, 0 otherwise
  close(10)

  if(temp.gt.1.d-6)then
    write(*,*)
    write(*,*)'Finite temperature is not available yet, I will run GS HFB '
    write(*,*)
    temp=0.d0
  endif


  if(BCs.gt.0) then
     openBCs = .true.       !OBC
  else 
     openBCs = .false.      !PBC
  endif
  calc_CN = .false.

  nsites = product(nk)

  allocate(Vext(nsites))
  Vext=0.d0  
  if(IVext.gt.0)then
    open(10,file='external_random_field.input',status='old')
    do i=1,nsites
      read(10,*)x,y,Vext(i)
    enddo
    close(10)
  endif

  allocate(Hzero(2*nsites,2*nsites))
  Hzero=cmplx(0.d0,0.d0)
  if(IH0.gt.0)then
    open(10,file='free_hamiltonian.input',status='old')
    do i=1,2*nsites
      do j=1,2*nsites
        read(10,*)Hzero(i,j)
      enddo
    enddo
    close(10)
  endif

  tol   = 1.d-6

  write(Lxstr, '(i4)') nk(1)
  write(Lystr, '(i4)') nk(2)
  write(Ustr,'(f5.1)') U
  write(mustr, '(f5.2)') mu
  write(hzstr, '(f5.2)') hz
  Ustr  = trim(adjustl(Ustr))
  mustr = trim(adjustl(mustr))
  hzstr = trim(adjustl(hzstr))

 filename=trim(adjustl("HFB_results_Lx"))//trim(adjustl(Lxstr))//trim(adjustl('_Ly')) &
      //trim(adjustl(Lystr))//trim(adjustl('_U'))//trim(adjustl(Ustr))//trim(adjustl('_mu')) &
      //trim(adjustl(mustr))//trim(adjustl('_hz'))//trim(adjustl(hzstr))//trim(adjustl('.dat'))
 open(50,file=filename,status='unknown')


 theta_twist(1) = kpoint_x             !thetax_up
 theta_twist(2) = -1.d0*theta_twist(1) !thetax_dn
 theta_twist(3) = kpoint_y             !thetay_up
 theta_twist(4) = -1.d0*theta_twist(3) !thetay_dn


 theta_twist(1) = 2*pi*theta_twist(1)
 theta_twist(2) = 2*pi*theta_twist(2)
 theta_twist(3) = 2*pi*theta_twist(3)
 theta_twist(4) = 2*pi*theta_twist(4)

 call scf_loop(nk,nsites,iH0,Hzero,thop,msite,nmesh,phi,theta_twist, lambda,                   &
     &         U, mu, hz, temp, Vext, alpha, beta, tol, openBCs,calc_CN)
 close(50)

end program bcs_openBCs

!------------------------------------------------------------------!

! SCF loop
subroutine scf_loop(nk,nsites,iH0,Hzero,thop,msite,nmesh,phi,theta_twist,lambda,                   &
    &               U,mu,hz,temp,Vext,alpha,beta, tol,openBCs,calc_CN)

  use mpiwrap
  implicit none

  integer, intent(in)            :: nk(3), nsites, msite, nmesh,iH0
  logical, intent(in)            :: openBCs, calc_CN
  real(kind=8), intent(in)       :: thop, lambda, U, mu, hz, temp, alpha, beta, tol, theta_twist(4)
  real(kind=8), intent(in)       :: Vext(nsites)
  complex(kind=8), intent(in)    :: Hzero(2*nsites,2*nsites)
  complex(kind=8), intent(inout) :: phi(nmesh,4*nsites,2*nsites)
  real(kind=8)                   :: old_energy, new_energy, energy_diff, gap_sum
  real(kind=8)                   :: nptcl, spinpol
  complex(kind=8), allocatable   :: gap(:), newgap(:)
  integer                        :: i, j
   
  allocate(gap(nsites))
  allocate(newgap(nsites))

  old_energy=0.d0
  gap(:) = 0.01d0

  open(22,file='run.information',status='unknown')
  do

     call get_gap(nk, nsites, iH0, Hzero, thop, theta_twist, lambda, U, mu, hz, temp, Vext, gap, newgap, &
          new_energy, alpha, beta, openBCs)

     ! update gap
     gap(:) = newgap(:)

     gap_sum=0.d0
     do j = 1, nsites, 1
        gap_sum = gap_sum + dsqrt(dimag(newgap(j))**2+dreal(newgap(j))**2)
     end do
     ! print gap
     write(22,*) 'new gap (sum/sites) :', gap_sum/nsites
     flush(22)

     ! check if energy has converged
     energy_diff = new_energy - old_energy
     write(22,*) 'E_diff : ', energy_diff, new_energy
     flush(22)

     ! update energy
     old_energy = new_energy

     ! if energy has converged, exit
     if( abs(energy_diff) .lt. tol ) then

        ! print gap
        open(2,file='order_parameter.out')
         do j = 1, nsites, 1
           write(2,*) j, dble(gap(j)), aimag(gap(j))
         end do
        close(2)

        ! print properties
        call get_dens_pol(nk, nsites, iH0, Hzero, thop, theta_twist, lambda, U, mu, hz, Vext, gap, openBCs)

        exit
     end if
 
  end do
  close(22)

  deallocate(gap)
  deallocate(newgap)

end subroutine scf_loop

!------------------------------------------------------------------!

! Get gap for fixed mu and h 
subroutine get_gap(nk, nsites, iH0, Hzero, thop, theta_twist, lambda, U, mu, hz, temp, Vext, gap, newgap, &
     new_energy, alpha, beta, openBCs) 

  use mpiwrap
  implicit none

  integer, intent(in)            :: nk(3), nsites, iH0
  logical, intent(in)            :: openBCs
  real(kind=8), intent(in)       :: thop, lambda, U, mu, hz, temp, theta_twist(4)
  real(kind=8), intent(in)       :: Vext(nsites)
  complex(kind=8), intent(in)    :: Hzero(2*nsites,2*nsites)
  real(kind=8), intent(inout)    :: new_energy
  real(kind=8), intent(in)       :: alpha, beta !mixing parameters
  complex(kind=8), intent(in)    :: gap(nsites)
  complex(kind=8), intent(inout) :: newgap(nsites)

  real(kind=8)                   :: rwork(50000)
  real(kind=8), allocatable      :: ener(:)
  integer, parameter             :: lwork=50000
  complex(kind=8)                :: work(lwork)

  integer                        :: Hdim, i, j
  real(kind=8)                   :: energy_sum, f
  complex(kind=8), allocatable   :: H(:,:)
  complex(kind=8)                :: pmat(2*nsites,2*nsites)
  real(kind=8), parameter        :: one=1.d0
  real(kind=8), parameter        :: zero=0.d0
  
  complex(kind=8)                :: sum1, sum2
  integer                        :: norb

  energy_sum = 0.d0
  Hdim       = 4*nsites

  allocate(H(Hdim,Hdim))
  allocate(ener(Hdim))

  newgap(:)  = 0.d0
  H(:,:)     = 0.d0
  ener(:)    = 0.d0

  ! build H
  if(iH0.eq.0)then
    if(openBCs) then
       call buildHopen(nsites, Hdim, H, nk, thop, lambda, gap, mu, hz, Vext)
    else
       call buildH(nsites, Hdim, H, nk, thop, theta_twist, lambda, gap, mu, hz, Vext)
    end if
  else
    call buildHext(nsites, Hdim, H, nk, Hzero, gap, mu, hz, Vext)
  endif

  ! check to make sure Hamiltonian is hermitian
  call check_Hermite_c(H,Hdim)

  ! diagonalize H
  call zheev("V","L",Hdim,H,Hdim,ener,work,lwork,rwork,i)

  norb = 0
     do i = 1, 4*nsites
        f = merge(1.0, 0.0, ener(i) .lt. 0) !1.0 / (exp(beta*ener(i)) + 1.d0)
           energy_sum = energy_sum + ener(i)*f
           norb = norb + int(f)
     enddo


  pmat(:,:) = 0.d0

!    H   (  V*  U )  =   (  V*  U )  ( -E    0  )
!        (  U*  V )      (  U*  V )  (  0    E  )

!ordering eigenvalues, first the negative ones

  ! calculate pairing matrix
  call zgemm('N','C',2*nsites,2*nsites,2*nsites,one,H(1:2*nsites,1:2*nsites), &               
       2*nsites,H(1+2*nsites:4*nsites,1:2*nsites),2*nsites,zero,pmat,2*nsites)

  !call mpiwrap_sum(newgap)

  do j = 1, nsites
     newgap(j) = U*pmat(j,j+nsites)

     ! mix old gap and new gap --> gap for next iteration
     newgap(j) = alpha*newgap(j) + beta*gap(j)
  end do

  new_energy = energy_sum

end subroutine get_gap

!------------------------------------------------------------------!

! Get particle number and polarization 
subroutine get_dens_pol(nk, nsites, iH0, Hzero, thop, theta_twist, lambda, U, mu, hz, Vext, gap, openBCs)

  use mpiwrap
  implicit none

  integer, intent(in)            :: nk(3), nsites, iH0
  logical, intent(in)            :: openBCs
  real(kind=8), intent(in)       :: thop, lambda, U, mu, hz, theta_twist(4)
  real(kind=8), intent(in)       :: Vext(nsites)
  complex(kind=8), intent(in)    :: Hzero(2*nsites,2*nsites)
  complex(kind=8), intent(in)    :: gap(nsites)

  real(kind=8)                   :: rwork(50000)
  real(kind=8), allocatable      :: ener(:)
  integer, parameter             :: lwork=50000
  complex(kind=8)                :: work(lwork)

  integer                        :: Hdim, Udim, i, j, k, ix, iy, Dimen,idim, jdim
  real(kind=8)                   :: f, ptcl_sum, pol_sum
  complex(kind=8), allocatable   :: H(:,:)
  complex(kind=8), allocatable   :: Ubcs(:,:),Vbcs(:,:)
  complex(kind=8)                :: dmat(2*nsites,2*nsites),dmath(2*nsites,2*nsites)
  complex(kind=8)                :: pmat(2*nsites,2*nsites),pmatbar(2*nsites,2*nsites)
  complex(kind=8)                :: Q(2*nsites,2*nsites), tmp(2*nsites,2*nsites)
  complex(kind=8), parameter     :: xi=dcmplx(0d0,1d0)
  real(kind=8), parameter        :: one=1.d0
  real(kind=8), parameter        :: zero=0.d0
  complex(kind=8)                :: sumup, sumdn, sumtot, en_gs
  real(kind=8)                   :: hop(nsites,nsites), KE, MUE, MUEH
  real(kind=8)                   :: rvec(3,2*nsites),rav(3),metric(3,3)

  ptcl_sum    = 0.d0
  pol_sum     = 0.d0
  KE          = 0.d0
  MUE         = 0.d0
  MUEH        = 0.d0
  Hdim        = 4*nsites
  Udim        = 2*nsites

  allocate(H(Hdim,Hdim))
  allocate(ener(Hdim))
  allocate(Ubcs(Udim,Udim))
  allocate(Vbcs(Udim,Udim))

  H(:,:)      = 0.d0
  ener(:)     = 0.d0

  ! build H
  if(iH0.eq.0)then
    if(openBCs) then
       call buildHopen(nsites, Hdim, H, nk, thop, lambda, gap, mu, hz, Vext)
    else
       call buildH(nsites, Hdim, H, nk, thop, theta_twist, lambda, gap, mu, hz, Vext)
    end if
  else
    call buildHext(nsites, Hdim, H, nk, Hzero, gap, mu, hz, Vext)
  endif



  open(2,file='nonzero_elements_ham.out',status='unknown')
  do i = 1, 4*nsites
     do j = 1, 4*nsites
        if(abs(H(i,j)).gt.1.d-12) then
          write(2,*)i,j,dble(H(i,j)),aimag(H(i,j))
        endif
     enddo
  enddo
  close(2)

  do i = 1, nsites
    do j = 1, nsites 
      if((i.ne.j).and.(abs(H(i,j)).gt.0))then
        hop(i,j) = 1
      else
        hop(i,j) = 0
      endif
    enddo
  enddo

  call check_Hermite_c(H,Hdim)

  ! diagonalize H
  call zheev("V","L",Hdim,H,Hdim,ener,work,lwork,rwork,i)

  open(2,file='eigenvalues_hfb_hamiltonian.out',status='unknown')
  do i = 1, 4*nsites
    write(2,*)i,ener(i)
  enddo
  close(2)

  ! get matrices U and V
!
!    H    (  V*  U  )  =   (  V*  U )  (  -E    0  )
!         (  U*  V  )      (  U*  V )  (   0    E  )


  Vbcs(1:2*nsites,1:2*nsites)=conjg(H(1:2*nsites         ,1:2*nsites))
  Ubcs(1:2*nsites,1:2*nsites)=conjg(H(1+2*nsites:4*nsites,1:2*nsites))


  ! get density matrix and pairing tensor

  dmat(:,:) = 0.d0
  pmat(:,:) = 0.d0
  pmatbar(:,:) = 0.d0


  ! dmat_{ij} =  (V* VT)_{i,j} = <c+_j c_i>
  call zgemm('N','C',2*nsites,2*nsites,2*nsites,one,H(1:2*nsites,1:2*nsites), &
       2*nsites,H(1:2*nsites,1:2*nsites),2*nsites,zero,dmat,2*nsites)

  !  pmat_{ij}    = (V* UT)_{ij} = <c_jc_i> 
  call zgemm('N','C',2*nsites,2*nsites,2*nsites,one,H(1:2*nsites,1:2*nsites), &
       2*nsites,H(1+2*nsites:4*nsites,1:2*nsites),2*nsites,zero,pmat,2*nsites)

  !  pmatbar_{ij} = (U* VT)_{ij} = <cdagger_i cdagger_j>
  call zgemm('N','C',2*nsites,2*nsites,2*nsites,one,H(1+2*nsites:4*nsites,1:2*nsites), &
       2*nsites,H(1:2*nsites,1:2*nsites),2*nsites,zero,pmatbar,2*nsites)



  open(2,file='density_matrix.out',status='unknown')
  do i = 1, 2*nsites
    do j = 1, 2*nsites
      write(2,*)i,j,dmat(i,j)
    enddo
  enddo
  close(2)
  open(2,file='pairing_matrix.out',status='unknown')
  do i = 1, 2*nsites
    do j = 1, 2*nsites
      write(2,*)i,j,pmat(i,j)
    enddo
  enddo
  write(2,*)
  write(2,*)
  do i = 1, 2*nsites
    do j = 1, 2*nsites
      write(2,*)i,j,pmatbar(i,j)
    enddo
  enddo
  close(2)

  sumup = 0.d0
  sumdn = 0.d0
  do i = 1, nsites
     sumup = sumup + dmat(i,i)
     sumdn = sumdn + dmat(i+nsites,i+nsites)
  enddo

   write(50,*) 'N_up: ', sumup, 'N_dn: ', sumdn, 'P: ', sumup - sumdn
   write(*,*) 'N_up: ', sumup, 'N_dn: ', sumdn, 'P: ', sumup - sumdn


! Calculate current-current correlation function

   call current_current(nk,nsites,Ubcs,Vbcs,ener(1:2*nsites))


! Calculate quantum metric tensor (2D)

   i=0
   do iy=1,nk(2)   
     do ix=1,nk(1)
       i=i+1

       rvec(1,i)=dble(ix-1)-nk(1)/2
       rvec(1,i+nsites)=rvec(1,i)
 
       rvec(2,i)=dble(iy-1)-nk(2)/2
       rvec(2,i+nsites)=rvec(2,i)
     enddo
   enddo

   rav=0.d0
   do i = 1, 2*nsites
     rav(1)=rav(1)+rvec(1,i)*dmat(i,i)
     rav(2)=rav(2)+rvec(2,i)*dmat(i,i)
   enddo

   metric=0.d0
   do idim=1,2
     do jdim=1,2
       do i = 1, 2*nsites
         do j = 1, 2*nsites
           metric(idim,jdim)=metric(idim,jdim)+rvec(idim,i)*rvec(jdim,j)*      &
                      (dmat(i,i)*dmat(j,j)-pmatbar(i,j)*pmat(i,j)  &
                      -dmat(i,j)*dmat(j,i))
           if(i.eq.j)then
             metric(idim,jdim)=metric(idim,jdim)+rvec(idim,i)*rvec(jdim,j)*dmat(i,j)
           endif
         enddo
       enddo
       metric(idim,jdim)=metric(idim,jdim)-rav(idim)*rav(jdim)
     enddo
   enddo    

   sumtot=sumup+sumdn
 
   open(2,file='localization_length')
   write(2,*)dble((metric(1,1))/sumtot)
   close(2)


   open(2,file='quantum_metric_tensor.out')
   do idim=1,2
     do jdim=1,2
       write(2,*)idim,jdim,metric(idim,jdim)/sumtot
     enddo
   enddo
   write(2,*)
   write(2,*)
   write(2,*)'<r_alpha>'
   do idim=1,2
     write(2,*)idim,rav(idim)
   enddo
   write(2,*)
   write(2,*)
   write(2,*)'<N> = ',sumtot
   write(2,*)
   write(2,*)
   write(2,*)'Lattice '
   write(2,*)
   do i=1,2*nsites
     write(2,*)rvec(1,i),rvec(2,i)
   enddo
   close(2)


   ! calcualte E0 from Green's functions
   en_gs = 0.d0
   do i = 1, nsites
      do j = 1, nsites
         if(i.eq.j) then
            !en_gs = en_gs - (mu+hz)*(0.5d0-dmat(i,i)) - (mu-hz)*(0.5d0-dmat(i+nsites,i+nsites))
            en_gs = en_gs + ((mu-Vext(i))+hz)*(0.5d0-dmat(i,i)) + ((mu-Vext(i))-hz)*(0.5d0-dmat(i+nsites,i+nsites))
            en_gs = en_gs + conjg(gap(i))*pmat(i,i+nsites) - gap(i)*pmatbar(i,i+nsites)
            MUE = MUE + dmat(i,i) + dmat(i+nsites,i+nsites)!(1.d0-dmat(i,i)) + (1.d0-dmat(i+nsites,i+nsites))
            MUEH = MUEH + (1.d0 - dmat(i,i) - dmat(i+nsites,i+nsites))
            !write(*,*) i, MUE, MUEH, (1.d0-dmat(i,i)), (1.d0-dmat(i+nsites,i+nsites))
         endif
         if(hop(i,j).gt.0) then
            en_gs = en_gs - thop*dmat(j,i) - thop*dmat(j+nsites,i+nsites)
            KE = KE - thop*dmat(j,i) - thop*dmat(j+nsites,i+nsites)
         endif
         ! spin-flipping hopping along -x-direction
         if ((i.eq.(j-1)) .and. .not. ((mod(i,nk(1)).eq.0).and.(j.gt.i).and.(mod(j,nk(1)).ne.0)) &
              .and. .not. ((mod(j,nk(1)).eq.0).and.(i.gt.j).and.(mod(i,nk(1)).ne.0))) then
            en_gs = en_gs - lambda*dmat(j,i+nsites)
            en_gs = en_gs + lambda*dmat(j+nsites,i)
         endif
         ! spin-flipping hopping along +x-direction
         if ((i.eq.(j+1)) .and. .not. ((mod(i,nk(1)).eq.0).and.(j.gt.i).and.(mod(j,nk(1)).ne.0)) &
              .and. .not. ((mod(j,nk(1)).eq.0).and.(i.gt.j).and.(mod(i,nk(1)).ne.0))) then
            en_gs = en_gs + lambda*dmat(j,i+nsites)
            en_gs = en_gs - lambda*dmat(j+nsites,i)
         endif
         ! spin-flipping hopping along -y-direction
         if ((i.eq.(j-nk(2))) .and. .not. ((mod(i,nk(1)).eq.0).and.(j.gt.i).and.(mod(j,nk(1)).ne.0)) &
              .and. .not. ((mod(j,nk(1)).eq.0).and.(i.gt.j).and.(mod(i,nk(1)).ne.0))) then             
            en_gs = en_gs - xi*lambda*dmat(j,i+nsites)
            en_gs = en_gs - xi*lambda*dmat(j+nsites,i)
         endif
         ! spin-flipping hopping along +y-direction
         if ((i.eq.(j+nk(2))) .and. .not. ((mod(i,nk(1)).eq.0).and.(j.gt.i).and.(mod(j,nk(1)).ne.0)) &
              .and. .not. ((mod(j,nk(1)).eq.0).and.(i.gt.j).and.(mod(i,nk(1)).ne.0))) then 
            en_gs = en_gs + xi*lambda*dmat(j,i+nsites)
            en_gs = en_gs + xi*lambda*dmat(j+nsites,i)
         endif

      enddo
   enddo

end subroutine get_dens_pol

!------------------------------------------------------------------!

! Get wfn labeled by thetax, thetay
subroutine get_wf(phi,msite,nmesh,nk,nsites,thop,theta_twist,lambda,U,mu,hz,Vext,gap)

  use mpiwrap
  implicit none

  integer, intent(in)            :: nk(3), nsites, nmesh, msite
  real(kind=8), intent(in)       :: thop, lambda, U, mu, hz, theta_twist(4)
  real(kind=8), intent(in)       :: Vext(nsites)
  complex(kind=8), intent(inout) :: phi(nmesh,4*nsites,4*nsites)
  complex(kind=8), intent(in)    :: gap(nsites)

  real(kind=8)                   :: rwork(50000)
  real(kind=8), allocatable      :: ener(:)
  integer, parameter             :: lwork=50000
  complex(kind=8)                :: work(lwork)

  integer                        :: Hdim, i, j, k
  real(kind=8)                   :: f, ptcl_sum, pol_sum
  complex(kind=8), allocatable   :: H(:,:)
  complex(kind=8), parameter     :: xi=dcmplx(0d0,1d0)
  real(kind=8), parameter        :: one=1.d0
  real(kind=8), parameter        :: zero=0.d0

  Hdim        = 4*nsites

  allocate(H(Hdim,Hdim))
  allocate(ener(Hdim))

  H(:,:)      = 0.d0
  ener(:)     = 0.d0

  ! build H
  call buildH(nsites, Hdim, H, nk, thop, theta_twist, lambda, gap, mu, hz, Vext)
  ! call buildHk

  ! check to make sure Hamiltonian is hermitian
  !call check_Hermite_c(H,Hdim)

  ! diagonalize H
  call zheev("V","L",Hdim,H,Hdim,ener,work,lwork,rwork,i)

  !do i = 1, 4*nsites
  !   if(ener(i).lt.0) then
  !      nocc(msite) = nocc(msite) + 1
  !   endif
  !enddo

  phi(msite,:,1:2*nsites) = H(:,1:2*nsites)
  
end subroutine get_wf

!------------------------------------------------------------------!

!build hamiltonian using free part from external file
subroutine buildHext(nsites, Hdim, Hopen, nk, Hzero, gap, mu, hz, Vext)

implicit none

real(kind=8), intent(in)       :: mu, hz
real(kind=8), intent(in)       :: Vext(nsites)
complex(kind=8), intent(in)    :: Hzero(2*nsites,2*nsites)
complex(kind=8), intent(in)    :: gap(nsites)
complex(kind=8), intent(inout) :: Hopen(Hdim,Hdim)

integer                        :: i,j,nk(3),nsites
integer                        :: Hdim

Hopen=cmplx(0.d0,0.d0)

do i=1,nsites,1
   do j=1,nsites,1
     if(i.eq.j)then
       Hopen(i,j)                   =  0.5d0*(Hzero(i,j)              -(mu-Vext(i))-hz)         !spin up
       Hopen(i+nsites,j+nsites)     =  0.5d0*(Hzero(i+nsites,j+nsites)-(mu-Vext(i))+hz)         !spin dn
       Hopen(i+2*nsites,j+2*nsites) = -Hopen(i,j)
       Hopen(i+3*nsites,j+3*nsites) = -Hopen(i+nsites,j+nsites)
       Hopen(i,j+3*nsites)          =  0.5*gap(i)
       Hopen(i+nsites,j+2*nsites)   = -0.5*gap(i)
       Hopen(i+2*nsites,j+nsites)   = -0.5*conjg(gap(i))
       Hopen(i+3*nsites,j)          =  0.5*conjg(gap(i))
     else
       Hopen(i,j)                   =  0.5d0*Hzero(i,j) 
       Hopen(i+nsites,j+nsites)     =  0.5d0*Hzero(i+nsites,j+nsites)
       Hopen(j+2*nsites,i+2*nsites) = -Hopen(i,j)
       Hopen(j+3*nsites,i+3*nsites) = -Hopen(i+nsites,j+nsites)
     endif
   enddo
enddo


end subroutine buildHext




! Build Hamiltonian matrix for system with open BCs
subroutine buildHopen(nsites,Hdim,Hopen,nk,thop,lambda,gap,mu,hz,Vext)

implicit none

real(kind=8), intent(in)       :: thop, mu, hz, lambda
real(kind=8), intent(in)       :: Vext(nsites)
complex(kind=8), intent(in)    :: gap(nsites)
complex(kind=8), intent(inout) :: Hopen(Hdim,Hdim)

integer                        :: i,j,nk(3),nsites
integer                        :: Hdim

Hopen=cmplx(0.d0,0.d0)

do i=1,nsites,1
   do j=1,nsites,1
      ! on-site terms
      if(i.eq.j) then
         Hopen(i,j)=0.5*(-1.0*(mu-Vext(i))-hz)
         Hopen(i+nsites,j+nsites)=0.5*(-1.0*(mu-Vext(i))+hz)
         Hopen(i+2*nsites,j+2*nsites)=0.5*((mu-Vext(i))+hz)
         Hopen(i+3*nsites,j+3*nsites)=0.5*((mu-Vext(i))-hz)
         Hopen(i,j+3*nsites)=0.5*gap(i)
         Hopen(i+nsites,j+2*nsites)=-0.5*gap(i)
         Hopen(i+2*nsites,j+nsites)=-0.5*conjg(gap(i))
         Hopen(i+3*nsites,j)=0.5*conjg(gap(i))
      endif

      ! spin-conserved hopping -x direction
      if ( (i.eq.(j-1))  .and. (mod(i,nk(1)).ne.0) )   then

         Hopen(i,j)=-0.5*thop
         Hopen(i+nsites,j+nsites)=-0.5*thop
         Hopen(j+2*nsites,i+2*nsites)=0.5*thop
         Hopen(j+3*nsites,i+3*nsites)=0.5*thop

      endif

      ! spin-conserved hopping +x direction
      if ( (i.eq.(j+1))  .and. (mod(j,nk(1)).ne.0) )   then

         Hopen(i,j)=-0.5*thop
         Hopen(i+nsites,j+nsites)=-0.5*thop
         Hopen(j+2*nsites,i+2*nsites)=0.5*thop
         Hopen(j+3*nsites,i+3*nsites)=0.5*thop

      endif

      ! spin-conserved hopping -y direction
      if ( (i.eq.(j-nk(1))) )   then

         Hopen(i,j)=-0.5*thop
         Hopen(i+nsites,j+nsites)=-0.5*thop
         Hopen(j+2*nsites,i+2*nsites)=0.5*thop
         Hopen(j+3*nsites,i+3*nsites)=0.5*thop

      endif

      ! spin-conserved hopping +y direction
      if ( (i.eq.(j+nk(1))) )   then

         Hopen(i,j)=-0.5*thop
         Hopen(i+nsites,j+nsites)=-0.5*thop
         Hopen(j+2*nsites,i+2*nsites)=0.5*thop
         Hopen(j+3*nsites,i+3*nsites)=0.5*thop

      endif


      ! spin-conserved hopping b/w nearest neighbors 
      !if (((i.eq.(j+1)).or.(i.eq.(j-1)).or.(i.eq.(j-nk(2))).or.(i.eq.(j+nk(2)))).and. &
      !   .not. ((mod(i,nk(1)).eq.0).and.(j.gt.i).and.(mod(j,nk(1)).ne.0)) .and. &
      !   .not. ((mod(j,nk(1)).eq.0).and.(i.gt.j).and.(mod(i,nk(1)).ne.0))) then
      !   Hopen(i,j)=-0.5*thop
      !   Hopen(i+nsites,j+nsites)=-0.5*thop
      !   !Hopen(i+2*nsites,j+2*nsites)=0.5*thop
      !   !Hopen(i+3*nsites,j+3*nsites)=0.5*thop
      !   Hopen(j+2*nsites,i+2*nsites)=0.5*thop
      !   Hopen(j+3*nsites,i+3*nsites)=0.5*thop
      !endif
      
      ! spin-flipping hopping along -x-direction
   !   if ((i.eq.(j-1)) .and. .not. ((mod(i,nk(1)).eq.0).and.(j.gt.i).and.(mod(j,nk(1)).ne.0)) &
   !        .and. .not. ((mod(j,nk(1)).eq.0).and.(i.gt.j).and.(mod(i,nk(1)).ne.0))) then
       if ( (i.eq.(j-1))  .and. (mod(i,nk(1)).ne.0) )   then
         Hopen(i+nsites,j)=-0.5*lambda
         Hopen(i,j+nsites)=0.5*lambda
         !Hopen(i+3*nsites,j+2*nsites)=0.5*lambda
         !Hopen(i+2*nsites,j+3*nsites)=-0.5*lambda
         Hopen(j+3*nsites,i+2*nsites)=-0.5*lambda
         Hopen(j+2*nsites,i+3*nsites)=0.5*lambda
      endif
      ! spin-flipping hopping along +x-direction
   !   if ((i.eq.(j+1)) .and. .not. ((mod(i,nk(1)).eq.0).and.(j.gt.i).and.(mod(j,nk(1)).ne.0)) &
   !        .and. .not. ((mod(j,nk(1)).eq.0).and.(i.gt.j).and.(mod(i,nk(1)).ne.0))) then
      if ( (i.eq.(j+1))  .and. (mod(j,nk(1)).ne.0) )   then
         Hopen(i+nsites,j)=0.5*lambda
         Hopen(i,j+nsites)=-0.5*lambda
         !Hopen(i+3*nsites,j+2*nsites)=-0.5*lambda
         !Hopen(i+2*nsites,j+3*nsites)=0.5*lambda
         Hopen(j+3*nsites,i+2*nsites)=0.5*lambda
         Hopen(j+2*nsites,i+3*nsites)=-0.5*lambda
      endif
      ! spin-flipping hopping along -y-direction
     ! if ((i.eq.(j-nk(2))) .and. .not. ((mod(i,nk(1)).eq.0).and.(j.gt.i).and.(mod(j,nk(1)).ne.0)) &
     !      .and. .not. ((mod(j,nk(1)).eq.0).and.(i.gt.j).and.(mod(i,nk(1)).ne.0))) then 
      if ( (i.eq.(j-nk(1))) )   then
         Hopen(i+nsites,j)=cmplx(0.d0,-0.5*lambda,kind=8)
         Hopen(i,j+nsites)=cmplx(0.d0,-0.5*lambda,kind=8)
         !Hopen(i+3*nsites,j+2*nsites)=cmplx(0.d0,0.5*lambda,kind=8)
         !Hopen(i+2*nsites,j+3*nsites)=cmplx(0.d0,0.5*lambda,kind=8)
         Hopen(j+3*nsites,i+2*nsites)=cmplx(0.d0,0.5*lambda,kind=8)
         Hopen(j+2*nsites,i+3*nsites)=cmplx(0.d0,0.5*lambda,kind=8)
      endif
      ! spin-flipping hopping along +y-direction
  !    if ((i.eq.(j+nk(2))) .and. .not. ((mod(i,nk(1)).eq.0).and.(j.gt.i).and.(mod(j,nk(1)).ne.0)) &
  !         .and. .not. ((mod(j,nk(1)).eq.0).and.(i.gt.j).and.(mod(i,nk(1)).ne.0))) then 
      if ( (i.eq.(j+nk(1))) )   then
         Hopen(i+nsites,j)=cmplx(0.d0,0.5*lambda,kind=8)
         Hopen(i,j+nsites)=cmplx(0.d0,0.5*lambda,kind=8)
         !Hopen(i+3*nsites,j+2*nsites)=cmplx(0.d0,-0.5*lambda,kind=8)
         !Hopen(i+2*nsites,j+3*nsites)=cmplx(0.d0,-0.5*lambda,kind=8)
         Hopen(j+3*nsites,i+2*nsites)=cmplx(0.d0,-0.5*lambda,kind=8)
         Hopen(j+2*nsites,i+3*nsites)=cmplx(0.d0,-0.5*lambda,kind=8)
      endif
   enddo
enddo

end subroutine buildHopen

!------------------------------------------------------------------!

! Build Hamiltonian matrix for system with periodic BCs
subroutine buildH(nsites,Hdim,Hopen,nk,thop,theta_twist,lambda,gap,mu,hz,Vext)

implicit none

real(kind=8), intent(in)       :: thop, mu, hz, lambda, theta_twist(4)
real(kind=8), intent(in)       :: Vext(nsites)
complex(kind=8), intent(in)    :: gap(nsites)
complex(kind=8), intent(inout) :: Hopen(Hdim,Hdim)

integer                        :: i,j,nk(3),nsites
integer                        :: Hdim
complex(kind=8), parameter     :: xi=dcmplx(0d0,1d0)

do i=1,nsites,1
   do j=1,nsites,1
      ! on-site terms
      if(i.eq.j) then
         Hopen(i,j)=0.5*(-1.0*(mu-Vext(i))-hz)
         Hopen(i+nsites,j+nsites)=0.5*(-1.0*(mu-Vext(i))+hz)
         Hopen(i+2*nsites,j+2*nsites)=0.5*((mu-Vext(i))+hz)
         Hopen(i+3*nsites,j+3*nsites)=0.5*((mu-Vext(i))-hz)
         Hopen(i,j+3*nsites)=0.5*gap(i)
         Hopen(i+nsites,j+2*nsites)=-0.5*gap(i)
         Hopen(i+2*nsites,j+nsites)=-0.5*conjg(gap(i))
         Hopen(i+3*nsites,j)=0.5*conjg(gap(i))
      endif

      ! spin-conserved hopping b/w nearest neighbors
      !if (((i.eq.(j+1)).and.(mod(j,nk(1)).ne.0)).or. &
      !    ((i.eq.(j-1)).and.(mod(i,nk(1)).ne.0)).or. &
      !    (i.eq.(j-nk(2))).or.(i.eq.(j+nk(2))) &
      !     .or.(i.eq.(j+nk(1)*(nk(2)-1))) &
      !     .or.(i.eq.(j-nk(1)*(nk(2)-1))) &
      !     .or.((mod(i,nk(1)).eq.0).and.(i.eq.(j+nk(1)-1))) &
      !     .or.((mod(j,nk(1)).eq.0).and.(i.eq.(j-nk(1)+1)))) then
      !   Hopen(i,j)=-0.5*thop
      !   Hopen(i+nsites,j+nsites)=-0.5*thop
      !   !Hopen(i+2*nsites,j+2*nsites)=0.5*thop
      !   !Hopen(i+3*nsites,j+3*nsites)=0.5*thop
      !   Hopen(j+2*nsites,i+2*nsites)=0.5*thop
      !   Hopen(j+3*nsites,i+3*nsites)=0.5*thop
      !endif
      
       ! -x-direction
       if (((i.eq.(j-1)).and.(mod(i,nk(1)).ne.0)).or. &
           ((mod(i,nk(1)).eq.0).and.(i.eq.(j+nk(1)-1)))) then
          Hopen(i,j)=-0.5*thop*exp(-1.d0*xi*theta_twist(1)/nk(1))
          Hopen(i+nsites,j+nsites)=-0.5*thop*exp(-1.d0*xi*theta_twist(2)/nk(1))
          !Hopen(i+2*nsites,j+2*nsites)=0.5*thop
          !Hopen(i+3*nsites,j+3*nsites)=0.5*thop
!OKKIO
!          Hopen(j+2*nsites,i+2*nsites)=0.5*thop*exp(1.d0*xi*theta_twist(1)/nk(1))
!          Hopen(j+3*nsites,i+3*nsites)=0.5*thop*exp(1.d0*xi*theta_twist(2)/nk(1))
          Hopen(j+2*nsites,i+2*nsites)=-Hopen(i,j)
          Hopen(j+3*nsites,i+3*nsites)=-Hopen(i+nsites,j+nsites)

          Hopen(i+nsites,j)=-0.5*lambda*exp(1.d0*xi*(theta_twist(2)-theta_twist(1))/nk(1))
          Hopen(i,j+nsites)=0.5*lambda*exp(-1.d0*xi*(theta_twist(2)-theta_twist(1))/nk(1))
      !    !Hopen(i+3*nsites,j+2*nsites)=0.5*lambda
      !    !Hopen(i+2*nsites,j+3*nsites)=-0.5*lambda
          Hopen(j+3*nsites,i+2*nsites)=-0.5*lambda*exp(-1.d0*xi*(theta_twist(2)-theta_twist(1))/nk(1))
          Hopen(j+2*nsites,i+3*nsites)=0.5*lambda*exp(1.d0*xi*(theta_twist(2)-theta_twist(1))/nk(1))
       endif
       ! +x-direction
       if (((i.eq.(j+1)).and.(mod(j,nk(1)).ne.0)).or. &
            ((mod(j,nk(1)).eq.0).and.(i.eq.(j-nk(1)+1)))) then
          Hopen(i,j)=-0.5*thop*exp(1.d0*xi*theta_twist(1)/nk(1))
          Hopen(i+nsites,j+nsites)=-0.5*thop*exp(1.d0*xi*theta_twist(2)/nk(1))
          !Hopen(i+2*nsites,j+2*nsites)=0.5*thop
          !Hopen(i+3*nsites,j+3*nsites)=0.5*thop
!OKKIO
!          Hopen(j+2*nsites,i+2*nsites)=0.5*thop*exp(-1.d0*xi*theta_twist(1)/nk(1))
!          Hopen(j+3*nsites,i+3*nsites)=0.5*thop*exp(-1.d0*xi*theta_twist(2)/nk(1))
          Hopen(j+2*nsites,i+2*nsites)=-Hopen(i,j)
          Hopen(j+3*nsites,i+3*nsites)=-Hopen(i+nsites,j+nsites)

          Hopen(i+nsites,j)=0.5*lambda*exp(-1.d0*xi*(theta_twist(2)-theta_twist(1))/nk(1))
          Hopen(i,j+nsites)=-0.5*lambda*exp(1.d0*xi*(theta_twist(2)-theta_twist(1))/nk(1))
      !    !Hopen(i+3*nsites,j+2*nsites)=-0.5*lambda
      !    !Hopen(i+2*nsites,j+3*nsites)=0.5*lambda
          Hopen(j+3*nsites,i+2*nsites)=0.5*lambda*exp(1.d0*xi*(theta_twist(2)-theta_twist(1))/nk(1))
          Hopen(j+2*nsites,i+3*nsites)=-0.5*lambda*exp(-1.d0*xi*(theta_twist(2)-theta_twist(1))/nk(1))
       endif
       ! -y-direction
       if ((i.eq.(j-nk(1))).or.(i.eq.(j+nk(1)*(nk(2)-1)))) then
          Hopen(i,j)=-0.5*thop*exp(-1.d0*xi*theta_twist(3)/nk(2))
          Hopen(i+nsites,j+nsites)=-0.5*thop*exp(-1.d0*xi*theta_twist(4)/nk(2))
          !Hopen(i+2*nsites,j+2*nsites)=0.5*thop
          !Hopen(i+3*nsites,j+3*nsites)=0.5*thop
!OKKIO
!          Hopen(j+2*nsites,i+2*nsites)=0.5*thop*exp(1.d0*xi*theta_twist(3)/nk(2))
!          Hopen(j+3*nsites,i+3*nsites)=0.5*thop*exp(1.d0*xi*theta_twist(4)/nk(2))
          Hopen(j+2*nsites,i+2*nsites)=-Hopen(i,j)
          Hopen(j+3*nsites,i+3*nsites)=-Hopen(i+nsites,j+nsites)

          Hopen(i+nsites,j)=cmplx(0.d0,-0.5*lambda,kind=8)*exp(1.d0*xi*(theta_twist(4)-theta_twist(3))/nk(2))
          Hopen(i,j+nsites)=cmplx(0.d0,-0.5*lambda,kind=8)*exp(-1.d0*xi*(theta_twist(4)-theta_twist(3))/nk(2))
      !    !Hopen(i+3*nsites,j+2*nsites)=cmplx(0.d0,0.5*lambda,kind=8)
      !    !Hopen(i+2*nsites,j+3*nsites)=cmplx(0.d0,0.5*lambda,kind=8)
          Hopen(j+3*nsites,i+2*nsites)=cmplx(0.d0,0.5*lambda,kind=8)*exp(-1.d0*xi*(theta_twist(4)-theta_twist(3))/nk(2))
          Hopen(j+2*nsites,i+3*nsites)=cmplx(0.d0,0.5*lambda,kind=8)*exp(1.d0*xi*(theta_twist(4)-theta_twist(3))/nk(2))
       endif
       ! +y-direction
       if ((i.eq.(j+nk(1))).or.(i.eq.(j-nk(1)*(nk(2)-1)))) then
          Hopen(i,j)=-0.5*thop*exp(1.d0*xi*theta_twist(3)/nk(2))
          Hopen(i+nsites,j+nsites)=-0.5*thop*exp(1.d0*xi*theta_twist(4)/nk(2))
          !Hopen(i+2*nsites,j+2*nsites)=0.5*thop
          !Hopen(i+3*nsites,j+3*nsites)=0.5*thop
!OKKIO
!          Hopen(j+2*nsites,i+2*nsites)=0.5*thop*exp(-1.d0*xi*theta_twist(3)/nk(2))
!         Hopen(j+3*nsites,i+3*nsites)=0.5*thop*exp(-1.d0*xi*theta_twist(4)/nk(2))
          Hopen(j+2*nsites,i+2*nsites)=-Hopen(i,j)
          Hopen(j+3*nsites,i+3*nsites)=-Hopen(i+nsites,j+nsites)

          Hopen(i+nsites,j)=cmplx(0.d0,0.5*lambda,kind=8)*exp(-1.d0*xi*(theta_twist(4)-theta_twist(3))/nk(2))
          Hopen(i,j+nsites)=cmplx(0.d0,0.5*lambda,kind=8)*exp(1.d0*xi*(theta_twist(4)-theta_twist(3))/nk(2))
      !    !Hopen(i+3*nsites,j+2*nsites)=cmplx(0.d0,-0.5*lambda,kind=8)
      !    !Hopen(i+2*nsites,j+3*nsites)=cmplx(0.d0,-0.5*lambda,kind=8)
          Hopen(j+3*nsites,i+2*nsites)=cmplx(0.d0,-0.5*lambda,kind=8)*exp(1.d0*xi*(theta_twist(4)-theta_twist(3))/nk(2))
          Hopen(j+2*nsites,i+3*nsites)=cmplx(0.d0,-0.5*lambda,kind=8)*exp(-1.d0*xi*(theta_twist(4)-theta_twist(3))/nk(2))
       endif

   enddo
enddo

end subroutine buildH

!------------------------------------------------------------------!

! Build Hamiltonian matrix for system with periodic BCs
subroutine buildHk(nsites,Hdim,H,nk,thop,lambda,gap,mu,hz)

implicit none

real(kind=8), intent(in)       :: thop, mu, hz, lambda
complex(kind=8), intent(in)    :: gap
complex(kind=8), intent(inout) :: H(Hdim,Hdim)

integer                        :: i,j,nk(3),nsites
integer                        :: Hdim
complex(kind=8), parameter     :: xi=dcmplx(0d0,1d0)

!eku = -2*thop*(cos(kx)+cos(ky))-mu-hz
!ekd = -2*thop*(cos(kx)+cos(ky))-mu+hz
!soc = 2*lambda*(sin(ky)+xi*sin(kx))

!H(1,1) = eku
!H(2,1) = soc  ; H(2,2) = ekd 
!H(3,1) = 0.d0 ; H(3,2) = -gap ; H(3,3) = -eku
!H(4,1) = gap  ; H(4,2) = 0.d0 ; H(4,3) = conjg(soc) ; H(4,4) = -ekd

end subroutine buildHk

!------------------------------------------------------------------!

!check the Hermition of the matrix H.
subroutine check_Hermite_c(H,n)
implicit none
integer, intent(IN):: n
complex(kind=8), intent(IN):: H(n,n)
real(kind=8):: error
integer:: i, j

error=0d0
do i=1, n
   do j=i+1,n
      !if(abs(H(i,j)-conjg(H(j,i))).gt.1d-8) then
      !   write(*,*) i, j, H(i,j), H(j,i)
      !endif
      error=error+ abs(H(i,j)-conjg(H(j,i)))
   enddo 
enddo 

if (error>1d-8) then 
   print *, 'H not Hermitian', error
   stop
endif 

end subroutine check_Hermite_c

!Evaluate the Inverse of one complex Matrix.
subroutine inverse(a,n)
implicit none
integer,intent(IN)::n
complex(kind=8),intent(INOUT)::a(n,n)
complex(kind=8),allocatable::work(:)
complex(kind=8)::work_test(1)
integer::lwork
integer::ipiv(n),info


call zgetrf(n,n,a,n,ipiv,info)
if (info < 0) then
  write(*,"(' zgetrf error info=',I4)") info
  write(*,"(' The',I4,'-th parameter is illegal.')") info
  stop
else if (info > 0) then
  write(*,"(' zgetrf error info=',I4)") info
  write(*,"(' Factorization is done. But the ',I4,'-th diagonal element is zero.')") info
  stop
! we can calculate determinant if everything is fine
end if


!get the best lwork
call zgetri(n,a,n,ipiv,work_test,-1,info)
lwork=work_test(1)!;lwork=10*n
allocate(work(lwork))


!get the inverse
call zgetri(n,a,n,ipiv,work,lwork,info)


!deallocate arrays
deallocate(work)


!Check the point
if (info < 0) then
  write(*,"(' zgetri error info=',I4)") info
  write(*,"(' The',I4,'-th parameter is illegal.')") info
  stop
else if (info > 0) then
  write(*,"(' zgetri error info=',I4)") info
  write(*,"(' The ',I4,'-th diagonal element is zero. Inversion cannot be performed.')") info
  stop
end if

end subroutine inverse

subroutine calc_ChernN(phi,nk,nsites,Nthetax,Nthetay,nmesh)
implicit none

integer, intent(in)           :: nk(3), nsites, nmesh, Nthetax, Nthetay
complex(kind=8), intent(in)   :: phi(nmesh,4*nsites,2*nsites)
complex(kind=8)               :: Ulink(nmesh,4)
complex(kind=8)               :: omega
complex(kind=8)               :: ChernN
integer                       :: j, thetax, thetay, s1, s2, s3, s4

complex(kind=8), parameter    :: xi=dcmplx(0d0,1d0)
real(kind=8),parameter        :: pi=3.14159265358979323846264338

write(50,*) 'calculating Chern Number...'

! calculate Chern number for system with PBC
ChernN = 0.d0
Ulink(:,:) = 0.d0

j = 1
do thetay = 0, Nthetax-1
   do thetax = 0, Nthetay-1

      ! get site numbers for overlap
      s1 = (thetax+1) + nk(2)*thetay
      s2 = (mod(thetax+1,Nthetax)+1) + nk(2)*thetay
      s3 = (mod(thetax+1,Nthetax)+1) + nk(2)*mod(thetay+1,Nthetay)
      s4 = (thetax+1) + nk(2)*mod(thetay+1,Nthetay)

      ! calculate overlap on each 4 site square within the theta mesh grid
      call deter_overlap(4*nsites,2*nsites,phi(s1,:,1:2*nsites),phi(s2,:,1:2*nsites),Ulink(j,1))
      call deter_overlap(4*nsites,2*nsites,phi(s2,:,1:2*nsites),phi(s3,:,1:2*nsites),Ulink(j,2))
      call deter_overlap(4*nsites,2*nsites,phi(s3,:,1:2*nsites),phi(s4,:,1:2*nsites),Ulink(j,3))
      call deter_overlap(4*nsites,2*nsites,phi(s4,:,1:2*nsites),phi(s1,:,1:2*nsites),Ulink(j,4))

      Ulink(j,1)=Ulink(j,1)/abs(Ulink(j,1))
      Ulink(j,2)=Ulink(j,2)/abs(Ulink(j,2))
      Ulink(j,3)=Ulink(j,3)/abs(Ulink(j,3))
      Ulink(j,4)=Ulink(j,4)/abs(Ulink(j,4))

      omega=product(Ulink(j,:))
      omega=datan2(dimag(omega),dble(omega))

      ChernN = ChernN + (1.d0/(2*pi))*omega
      
      j = j + 1

   enddo
enddo

write(50,*) 'Chern Number : ', dreal(ChernN), dimag(ChernN)

end subroutine calc_ChernN

subroutine deter_overlap(n,m,phi_l,phi_r,imp)
use caldet_module
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(1,n,m),phi_r(1,n,m)
complex(kind=8)           :: a(n,m), b(n,m)
complex(kind=8),intent(INOUT)::imp
complex(kind=8):: ovlpinv_temp(m,m)
complex(kind=8)::one=dcmplx(1.d0,0.d0)
complex(kind=8)::zero=dcmplx(0.d0,0.d0)

a(:,:) = phi_l(1,:,:)
b(:,:) = phi_r(1,:,:)

 call zgemm('C','N',m,m,n,one,a,n,b,n,zero,ovlpinv_temp,m)
 call caldet(m,ovlpinv_temp(1:m,1:m),imp)
end subroutine deter_overlap
