program main
  use solver
  use function
  use precision

  implicit none
 
  real(dp), allocatable :: phi(:,:,:) ! campo scalare reale in 2+1 dimensioni
  real(dp), allocatable :: pi(:,:,:)  ! impulso fittizio per la dinamica del montecarlo
  !real(dp), allocatable :: correlatore(:,:,:) ! correlatore a due punti
  real(dp) :: lambda     ! coupling dimensionale dell'interazione
  real(dp) :: m2      ! massa quadra dello scalare
  real(dp) :: k_p       ! parametro di massa adimensionale
  real(dp) :: g_p       ! coupling adimensionale
  real(dp) :: obs_prova  ! un buon parametro di prova può essere somma sui punti del reticolo del campo
                         ! infatti per la simmetria della teoria questo dovrebbe essere nullo.
  real(dp) :: split_en, cum_split_en  
  integer  :: L       ! parametro di discretizzazione del reticolo
  real(dp) :: a       ! passo reticolare lo prendo uguale sia nella direzione spaziale che temporale
  real(dp) :: tau_p
  integer  :: n_p
  integer  :: Nstep   ! numero di step della catena di Markov (indicativamente credo sui 10^5)
  integer  :: i, j, m_sep, num_args, num, num_int
  integer :: ii,jj,kk
  
  character(100), allocatable :: args(:)
  real(dp), allocatable :: param(:)

  
  if (command_argument_count() < 8) then
      write(*,*) "Nstep L a m2 lambda m_sep tau n"
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
   tau=param(7)
   n=int(param(8))
   
   g=lambda*a
   k=a*a*m2+ 6.0_dp ! in generale sarebbe m2*a^2 + 2*D dove D è la dimensione D=2+1 in questo caso
  
  allocate(phi(0:L-1,0:L-1,0:L-1))
  allocate(pi(0:L-1,0:L-1,0:L-1))
  ! allocate(correlatore(0:L-1,0:L-1,0:L-1))

  call random_number(phi) ! inizializzo randomicamente il campo
  obs_prova=sum(phi)
  j=0 ! contatore della misura

   open(9, file='solution.dat')
   write(9,*)  j  ,  obs_prova  

   cum_split_en=0
  
  do i=1, Nstep
     ! inizializzo l'impulso
     call init_pi(L,pi)
  
     call metropolis(phi,pi, split_en)
     split_en=abs(split_en)
     
     ! per calcolare il correlatore devo capire quando arriva a termalizzazione
     ! una volta raggiunta la termalizzazione sommo: sum_phi phi(x)phi(0)
     ! per ridurre la correlazione non devo misurare ad ogni step
     
     if (mod(i,m_sep)==0) then
      j=j+1
      obs_prova=sum(phi)
      cum_split_en=cum_split_en + split_en
      print *, i , split_en
      write(9,*) j , obs_prova !, split_en , cum_split_en
   end if
   split_en=0       
     ! correlatore(:,:,:)+=phi(:,:,:)*phi(1,1,1)  
  end do

  j=real(j,dp)
  cum_split_en=cum_split_en/j
  print*, j, cum_split_en, num_int

    close(9)
     
end program main
  
