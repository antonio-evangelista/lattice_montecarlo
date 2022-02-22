program test
  use precision    
  implicit none

  integer :: N, j, i, t , W
  real(dp), allocatable :: obs(:) ! osservabile
  real(dp), allocatable :: B(:)
  real(dp), allocatable :: auto_corr(:) ! funzione di autocorrelazione normalizzata
  real(dp), allocatable :: tau_int(:), delta_tau_int(:) ! autocorrelazione integrata
  character(100) :: arg 
  real(dp) :: obs_bar, norm

  ! call get_command_argument(1,arg)
  ! read(arg,*) N
  N=2000
  allocate(obs(0:N))
  
  open(9, file='solution.dat')
  open(10, file='auto_corr.dat')
  open(11, file='tau_int.dat')
  
  read(9,*) j, obs(j)
  
  allocate(B(0:N))
  allocate(auto_corr(0:N-1))
  allocate(tau_int(1:N))
  allocate(delta_tau_int(1:N))

  obs_bar=real(sum(obs)/N,8)
  print*, obs_bar
  
  do i=0, N
     B(i)=obs(i)-obs_bar
  end do
  
  do t= 0, N-1
     auto_corr(t)=0
     do i=1,N-t
        auto_corr(t)=auto_corr(t) + B(i)*B(i+t)
     end do
     auto_corr(t)=auto_corr(t)*(1.0_dp/real(N-t))
     norm=auto_corr(0)
     auto_corr(t)=auto_corr(t)/norm
     write(10, *) t , auto_corr(t)
  end do  

  do W= 1 , N-1
     tau_int(W)=0
     do t= 1, W
        tau_int(W)=tau_int(W) + auto_corr(t)
     end do
     tau_int(W)= 0.5_dp + tau_int(W)
     delta_tau_int(W)=sqrt((real(4*W + 2)/real(N)))*tau_int(W)
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
  
  

  
  
  

  

  

  
  
  
