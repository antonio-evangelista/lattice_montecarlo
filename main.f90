program main
  use solver
  use functions

  implicit none

  real(dp), allocatable :: phi(:,:,:) ! campo scalare reale in 2+1 dimensioni
  real(dp), allocatable :: pi(:,:,:)  ! impulso fittizio per la dinamica del montecarlo
  real(dp) :: lambda     ! coupling dimensionale dell'interazione
  real(dp) :: m2      ! massa quadra dello scalare
  real(dp) :: k       ! parametro di massa adimensionale
  real(dp) :: g       ! coupling adimensionale
  integer  :: L       ! parametro di discretizzazione del reticolo
  integer  :: Nstep   ! numero di step della catena di Markov (indicativamente credo sui 10^5)
  real(dp) :: a       ! passo reticolare lo prendo uguale sia nella direzione spaziale che temporale
  integer  :: i
  character(50), allocatable :: args(:)
  real(dp), allocatable :: param(:)

  
  if (command_argument_count() < 5) then
      write(*,*) "Nstep L a m2 lambda"
      stop
  endif

  num_args = command_argument_count()
  allocate(args(num_args))
  allocate(param(num_args))
    
  do i = 1, num_args
     call get_command_argument(i, args(i))   
     read(args(i),*) param(i)
   end do

   Nstep=int(param(1))
   L=(param(2))
   a=param(3)
   m2=param(4)
   lambda=param(5)
  
  allocate(phi(L,L,L))
  allocate(pi(L,L,L))

  call random_number(phi)
 
  do i=1, N
     call(init_pi(L,pi)) ! inizializzo l'impulso,forse non serve allocare pi all'interno di init_pi
     call(metropolis(phi,pi))

     ! per calcolare il correlatore devo capire quando arriva a termalizzazione
     ! una volta raggiunta la termalizzazione sommo: sum_x phi(x)phi(0)esp(-S(phi))
  end do
  
  
  


end program main
  
