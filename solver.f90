module solver
  use function
  implicit none

  public :: metropolis
  public :: leapfrog

contains

  subroutine simplex_2(u, v, dt, mul)
    real(dp), intent(in) :: u(:,:,:) !campo
    real(dp), intent(in) :: v(:,:,:) !impulso fittizio
    real(dp), allocatable :: u_i(:,:,:) 
    real(dp), intent(out) :: u_f(:,:,:)
    real(dp), allocatable :: v_i(:,:,:)
    real(dp), intent(out) :: v_f(:,:,:)
    real(dp) :: grad_s1 !la derivata del potenziale, quindi dell'azione(devo ancora definire sia azione che sua derivata)
    real(dp) :: grad_s2
    real(dp) :: dt
    real(dp) :: mul
    integer :: dim

    dim=size(u)
    allocate(u_i(dim))
    allocate(u_f(dim))
    allocate(v_i(dim))
    allocate(v_f(dim))  

    dt= mul*dt

    u_i=u
    v_i=v
    
    grad_s1=call(grad_s(u_i))
    u_f=u_i + dt*(v_f -0.5_dp*dt*grad_s1)
    grad_s2=call(grad_s(u_f))
    v_f=v_i - 0.5_dp*dt*(grad_s1+grad_s2)

  end subroutine simplex_2
  
  
  !l'idea è di implementare un integratore simplettico del quardo ordine
  subroutine leapfrog(u, v, tau, n)  !devo capire come selzione u o v
    real(dp), intent(inout) :: u(:) !prende il ruolo del campo
    real(dp), intent(inout) :: v(:) !prende il ruolo dell'impulso associato al campo
    real(dp) :: dt !passo di integrazione nel tempo di evoluzione del microstato: n*dt=tau
    real(dp) :: tau
    integer :: n,i,j

    dt=(tau/n)_dp
    do i=0, n-1
      do j=0,2
        u(:,:,:)=call(simplex_2(u,v,alpha))
        v(:,:,:)=call(simplex_2(u,v,alpha))
      end do
    end do
    
  end subroutine leapfrog
  
  subroutine metropolis(phi, pi)
    real(dp) :: H_i
    real(dp) :: H_f
    real(dp) :: delta_H
    real(dp), allocatable   :: phi_i(:,:,:)
    real(dp), allocatable   :: pi_i(:,:,:)
    real(dp), allocatable   :: phi_f(:,:,:)
    real(dp), allocatable   :: pi_f(:,:,:)
    real(dp), allocatable   :: phi(:,:,:)
    real(dp)  :: prob
    real(dp)  :: r
    integer   :: n

    n=size(phi)
    phi_i=phi
    pi_i=pi
    allocate(phi_f(n))
    allocate(pi_f(n))

    phi_f =call(leapfrog(phi_i, pi_i,1))
    pi_f  =call(leapfrog(phi_i, pi_i,2))
    H_i   =call(H(phi_i,pi_i))
    H_f   =call(H(phi_f,pi_f))
    delta_H =H_f-H_i

    prob=min(1, exp(-delta_H))

 case(prob==1)
    phi=phi_f
 case(prob== exp(-delta_H))
    call random_number(r)
    if r < prob
    phi=phi_f
 else if
    phi=phi_i
 end if

end subroutine metropolis

end module solver
