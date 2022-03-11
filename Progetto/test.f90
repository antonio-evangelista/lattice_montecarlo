program test
  use precision    
  implicit none

  integer :: N, j, i, t , W, ii
  real(dp), allocatable :: obs(:) ! osservabile
  real(dp), allocatable :: B(:)
  real(dp), allocatable :: auto_corr(:) ! funzione di autocorrelazione normalizzata
  real(dp), allocatable :: tau_int(:), delta_tau_int(:) ! autocorrelazione integrata
  character(100) :: arg 
  real(dp) :: obs_bar, norm


  N=1200
  allocate(obs(1:N))
  
  open(9, file='solution.dat')
  open(10, file='auto_corr.dat')
  open(11, file='tau_int.dat')
  obs_bar=0

  do ii=1,N
     read(9,*) j, obs(j)
  end do
  
  allocate(B(1:N))
  allocate(auto_corr(0:N-1))
  allocate(tau_int(1:N-1))
  allocate(delta_tau_int(1:N-1))

  obs_bar=sum(obs)/real(N,8)
  print*, obs_bar
  
  do i=1, N
     B(i)=obs(i)-obs_bar
  end do

  do i=1,N
     auto_corr(0)=auto_corr(0) + B(i)**2
  end do
  norm=auto_corr(0)/real(N,8)
  write(10,*) 0, 1.0_dp
     
  do t= 1, N-1
     auto_corr(t)=0
     do i=1,N-t
        auto_corr(t)=auto_corr(t) + B(i)*B(i+t)
     end do
     auto_corr(t)=auto_corr(t)*(1.0_dp/real(N-t))
     auto_corr(t)=auto_corr(t)/norm
     print*, t, norm
     write(10, *) t , auto_corr(t)
  end do  

  do W= 1 , N-1
     tau_int(W)=0
     do t= 1, W
        tau_int(W)=tau_int(W) + auto_corr(t)
     end do
     tau_int(W)= 0.5_dp + tau_int(W)
     delta_tau_int(W)=sqrt((real(4*W + 2)/real(N)))*abs(tau_int(W))
     write(11, *) W , tau_int(W) , delta_tau_int(W)
  end do

  deallocate(B)
  deallocate(auto_corr)
  deallocate(tau_int)
  deallocate(delta_tau_int)

  close(9)
  close(10)
  close(11)

  end program
  
  

  
  
  

  

  

  
  
  
