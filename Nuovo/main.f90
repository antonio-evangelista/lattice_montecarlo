program main
  use solver
  use function
  use precision

  implicit none
 
  real(dp), allocatable :: phi(:,:,:) ! campo scalare reale in 2+1 dimensioni
  real(dp), allocatable :: pi(:,:,:)  ! impulso fittizio per la dinamica del montecarlo
  real(dp), allocatable :: correlatore(:) ! correlatore a due punti
  real(dp) :: lambda     ! coupling dimensionale dell'interazione
  real(dp) :: m2      ! massa quadra dello scalare
  real(dp) :: k       ! parametro di massa adimensionale
  real(dp) :: g       ! coupling adimensionale
  real(dp) :: obs_prova  ! un buon parametro di prova può essere somma sui punti del reticolo del campo
                         ! infatti per la simmetria della teoria questo dovrebbe essere nullo.
  real(dp) :: prova_pi
  integer  :: L       ! parametro di discretizzazione del reticolo
  real(dp) :: a       ! passo reticolare lo prendo uguale sia nella direzione spaziale che temporale
  real(dp) :: tau
  integer  :: n
  integer  :: Nstep  ! numero di step della catena di Markov (indicativamente credo sui 10^5)
  integer  ::  num_args, num_cut
  integer  :: i, j, m_sep, ii, jj, kk
  real(dp) :: split_en, cum_split, prova

  character(100), allocatable :: args(:)
  real(dp), allocatable :: param(:)

  
  if (command_argument_count() < 9) then
      write(*,*) "Nstep L a m2 lambda m_sep tau n num_cut"
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
   num_cut=int(param(9))
   
   g=lambda*a
   k=a*a*m2+ 6_dp ! in generale sarebbe m2*a^2 + 2*D dove D è la dimensione D=2+1 in questo caso
  
  allocate(phi(0:L-1,0:L-1,0:L-1))
  allocate(pi(0:L-1,0:L-1,0:L-1))
  allocate(correlatore(0:L-1))
  
  call random_number(phi)          ! inizializzo randomicamente il campo
  phi= phi - 0.5_dp
  j=0 ! contatore della misura
 
  open(9, file='solution.dat')
  !open(10, file='check.dat')
  !open(11, file='energy_err.dat', access='append')
  open(12, file='coorelatore.dat')
  obs_prova=0
  split_en=0
  cum_split=0
  correlatore=0
  
  do i=1, Nstep
     call init_pi(L,pi) ! inizializzo l'impulso
  !   do ii=0,L-1
  !      do jj=0, L-1
  !         do kk= 0, L-1
  !            write(10,*) pi(ii,jj,kk)
  !         end do
  !      end do
  !   end do
          
     call metropolis(phi, pi, split_en, k, g, tau, n)

     if (mod(i,m_sep)==0) then
        j=j+1
        if (j>=num_cut) then
           do ii= 0, L-1
              correlatore(ii)= correlatore(ii) + sum(phi(ii,:,:)*phi(0,0,0))
           end do
        end if
        !split_en=abs(split_en)
        !cum_split= cum_split + split_en
        obs_prova= obs_prova + sum(phi)
        write(9,*) j , sum(phi)
     end if 
  end do

  j=real(j,8)
  !cum_split=cum_split/j
  obs_prova=obs_prova/j
  print*, j, obs_prova, cum_split
  
  j=j-num_cut ! ridefinisco j per mediare il correlatore
  correlatore=correlatore/j
  do ii=0, L-1
     write(12,*) ii, correlatore(ii)
  end do
  !write(11,*) (tau/real(n,8))**2 , cum_split
  
  close(9)
  !close(10)
  !close(11)
  close(12)
  
end program main
  
