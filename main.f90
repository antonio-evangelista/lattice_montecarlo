program main
  use solver
  use functions

  implicit none

  real(dp), allocatable :: phi(:,:,:) ! campo scalare reale in 2+1 dimensioni
  real(dp), allocatable :: pi(:,:,:)  ! impulso fittizio per la dinamica del montecarlo
  real(dp), allocatable :: correlatore(:,:,:) ! correlatore a due punti
  real(dp) :: lambda     ! coupling dimensionale dell'interazione
  real(dp) :: m2      ! massa quadra dello scalare
  real(dp) :: k       ! parametro di massa adimensionale
  real(dp) :: g       ! coupling adimensionale
  real(dp) :: obs_prova  ! un buon parametro di prova può essere somma sui punti del reticolo del campo
                         ! infatti per la simmetria della teoria questo dovrebbe essere nullo.
  
  integer  :: L       ! parametro di discretizzazione del reticolo
  real(dp) :: a       ! passo reticolare lo prendo uguale sia nella direzione spaziale che temporale
  integer  :: Nstep   ! numero di step della catena di Markov (indicativamente credo sui 10^5)
  integer  :: i, j, m_sep
  character(50), allocatable :: args(:)
  real(dp), allocatable :: param(:)

  
  if (command_argument_count() < 6) then
      write(*,*) "Nstep L a m2 lambda m_sep"
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
   L=int(param(2))
   a=param(3)
   m2=param(4)
   lambda=param(5)
   m_sep=int(param(6))
   
   g=lambda*a
   k=a*a*m2+ 6_dp ! in generale sarebbe m2*a^2 + 2*D dove D è la dimensione D=2+1 in questo caso
  
  allocate(phi(L,L,L))
  allocate(pi(L,L,L))
  allocate(correlatore(L,L,L))

  call random_number(phi) ! inizializzo randomicamente il campo
  obs_prova=sum(phi)
  j=0 ! contatore della misura

  open(newunit=funit, file="solution.dat")
  write(funit,*) "j" , "obs_prova" ! , "correlatore(x0)" , "correlatore(x1)" , "correlatore(x2)"
  write(funit,*)  j  ,  obs_prova  ! , correlatore(:,1,1) , correlatore(1,:,1) , correlatore(1,1,:)  
   
  do i=1, Nstep
     call(init_pi(L,pi)) ! inizializzo l'impulso,forse non serve allocare pi all'interno di init_pi
     phi(:,:,:)= call(metropolis(phi,pi)
     
     ! per calcolare il correlatore devo capire quando arriva a termalizzazione
     ! una volta raggiunta la termalizzazione sommo: sum_phi phi(x)phi(0)
     ! per ridurre la correlazione non devo misurare ad ogni step
     if mod(i,m_sep)
      j=j+1
      obs_prova=sum(phi)
      write(funit,*) j , obs_prova
     end if
       
     ! correlatore(:,:,:)+=phi(:,:,:)*phi(1,1,1)
     
  end do

  close(funit)
  
end program main
  
