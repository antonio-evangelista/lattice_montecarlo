module solver
  use function
  use precision
  implicit none

  public :: metropolis
  public :: leapfrog
  !real(dp), public :: k, g, tau
  !integer, public :: n

contains
  
  subroutine leapfrog(u, v, k, g, tau, n) 
    real(dp), intent(inout) :: u(0:,0:,0:) !prende il ruolo del campo
    real(dp), intent(inout) :: v(0:,0:,0:) !prende il ruolo dell'impulso associato al campo
    real(dp), intent(in) :: k, g, tau
    integer, intent(in)  :: n
    real(dp), allocatable   :: v1(:,:,:) ! variabile di appoggi
    real(dp), allocatable   :: grad_s1(:,:,:)                  ! la derivata del potenziale hamiltoniano, cioè dell'azione
    real(dp) :: dt                       ! passo di integrazione nel tempo di evoluzione del microstato: n*dt=tau
    real(dp) :: molt  ! fattori che entrano nell'operatore simplettico di ordine 4 
    integer :: i, L 
   

    dt=tau/real(n,8)
    L=size(v,dim=1)
    allocate(v1(0:L-1, 0:L-1, 0:L-1))
    allocate(grad_s1(0:L-1,0:L-1,0:L-1))
    molt=1.0_dp
    
    do i=1, n
       !do j=1,3
       !   if (j==1 .or. j==3) then
       !      molt=1.351207191959658_dp
       !   else
       !      molt=-1.702414383919315_dp
       !   end if
       
          v1=v      
          call grad_s(u+0.5_dp*molt*dt*v, grad_s1, k, g) 
           v=v - molt*dt*grad_s1
          u=u + 0.5_dp*molt*dt*(v1+v)
       !end do
       end do
       
    deallocate(v1) 
    deallocate(grad_s1)
    
  end subroutine leapfrog
  
  subroutine metropolis(u, v, delta_H, k, g, tau, n) !poichè chiamo una funzione che richiede anche altri parametri devo specificarli anche qui?
    real(dp), intent(inout)  :: u(0:,0:,0:)
    real(dp), intent(inout)  :: v(0:,0:,0:)
    real(dp), intent(inout) :: delta_H
    real(dp), intent(in) :: k, g, tau
    integer, intent(in)  :: n
    real(dp), allocatable :: phi_f(:,:,:)
    real(dp), allocatable :: pi_f(:,:,:)
    real(dp) :: H_i
    real(dp) :: H_f
    real(dp)  :: prob
    real(dp)  :: r
    integer   :: L
   

    L=size(u,dim=1)
    
    allocate(phi_f(0:L-1,0:L-1,0:L-1))
    allocate(pi_f(0:L-1,0:L-1,0:L-1))
    phi_f=u 
    pi_f=v
    
    call leapfrog(phi_f, pi_f, k, g, tau, n)

    H_i   = H(u,v, k, g)
    H_f   = H(phi_f,pi_f, k, g)
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
