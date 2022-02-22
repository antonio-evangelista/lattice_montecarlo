module solver
  use function
  use precision
  implicit none

  public :: metropolis
  public :: leapfrog
  real(dp), public :: tau
  integer, public :: n

contains
  
  subroutine leapfrog(u, v, tau, n)  !devo capire come selzione u o v
    real(dp), intent(inout) :: u(0:,0:,0:) !prende il ruolo del campo
    real(dp), intent(inout) :: v(0:,0:,0:) !prende il ruolo dell'impulso associato al campo
    real(dp), intent(in) :: tau  
    integer, intent(in) :: n
    real(dp), allocatable   :: v1(:,:,:) ! variabile di appoggi
    real(dp), allocatable   :: grad_s1(:,:,:) ! la derivata del potenziale hamiltoniano, cioè dell'azione
    real(dp) :: dt       ! passo di integrazione nel tempo di evoluzione del microstato: n*dt=tau
    real(dp) :: molt     ! fattori che entrano nell'operatore simplettico di ordine 4 
    integer ::i,j, L

    dt=(tau/n)
    L=size(v,dim=1)
    allocate(v1(0:L-1,0:L-1,0:L-1))
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
          call grad_s(u+0.5_dp*molt*dt*v, grad_s1)
          v=v - molt*dt*grad_s1
          u=u + 0.5_dp*molt*dt*(v1+v)
       !end do
    end do

    deallocate(v1) 
    deallocate(grad_s1)
    
  end subroutine leapfrog
  
  subroutine metropolis(u, v, delta_H) 
    real(dp), intent(inout)  :: u(0:,0:,0:)
    real(dp), intent(in)  :: v(0:,0:,0:)
    real(dp), intent(inout) :: delta_H
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
    else
       u=u
    end if
  end if

 deallocate(phi_f)
 deallocate(pi_f)

end subroutine metropolis

end module solver
