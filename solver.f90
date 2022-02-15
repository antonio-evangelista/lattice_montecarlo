module solver
  use function
  use precision
  implicit none

  public :: metropolis
  public :: leapfrog
  real(dp), public :: tau, k, g
  integer, public :: n

contains
  
  subroutine leapfrog(u, v, tau, n)  !devo capire come selzione u o v
    real(dp), intent(inout) :: u(:,:,:) !prende il ruolo del campo
    real(dp), intent(inout) :: v(:,:,:) !prende il ruolo dell'impulso associato al campo
    real(dp), intent(in) :: tau  ! inizialmente dovrebbe essere un parametro libero che però se faccio varie simulazioni
                                 ! in base a quanto miglioro l'autocorrelazione, dt=tau è un pò poco mentre già un paio di
                                 ! iterazioni dovrebbero bastare
    integer, intent(in) :: n
    real(dp), allocatable   :: v1(:,:,:) ! variabile di appoggi
    real(dp), allocatable   :: grad_s1(:,:,:)                  ! la derivata del potenziale hamiltoniano, cioè dell'azione
    real(dp) :: dt                       ! passo di integrazione nel tempo di evoluzione del microstato: n*dt=tau
    real(dp) :: molt                      ! fattori che entrano nell'operatore simplettico di ordine 4 
    integer  :: dim
    
    
    integer ::i,j, L

    dt=(tau/n)
    L=size(v,dim=1)
    allocate(v1(L,L,L))
    allocate(grad_s1(L,L,L))

    do i=1, n
       do j=1,3
          if (j==1 .or. j==3) then
             molt=1.351207191959658_dp
          else
             molt=1.702414383919315_dp
          end if

          call grad_s(u+0.5_dp*molt*dt*v, k, g, grad_s1) 
          v1=v 
          v=v - molt*dt*grad_s1
          u=u + 0.5_dp*molt*dt*(v1+v)
       end do
    end do

    deallocate(v1) !occupo solo memoria altrimenti
    deallocate(grad_s1)
    
  end subroutine leapfrog
  
  subroutine metropolis(u, v) !poichè chiamo una funzione che richiede anche altri parametri devo specificarli anche qui?
    real(dp), intent(inout)  :: u(:,:,:)
    real(dp), intent(inout)  :: v(:,:,:)
    real(dp), allocatable :: phi_f(:,:,:)
    real(dp), allocatable :: pi_f(:,:,:)
    real(dp) :: H_i
    real(dp) :: H_f
    real(dp) :: delta_H
    real(dp)  :: prob
    real(dp)  :: r
    integer   :: L

    L=size(u,dim=1)
    
    allocate(phi_f(L,L,L))
    allocate(pi_f(L,L,L))
    phi_f=u 
    pi_f=v

    call leapfrog(phi_f, pi_f, tau, n)
    
    H_i   = H(u,v)
    H_f   = H(phi_f,pi_f)
    delta_H =H_f-H_i

    prob=min(1.0_dp, exp(-delta_H))

  if (prob==1) then
    u=phi_f
  else 
    call random_number(r)
    if (r < prob) then
       u=phi_f
    end if
  end if

 deallocate(phi_f)
 deallocate(pi_f)

end subroutine metropolis

end module solver
