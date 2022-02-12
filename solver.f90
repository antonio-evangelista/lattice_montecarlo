module solver
  use function
  implicit none

  public :: metropolis
  public :: leapfrog

contains
  
  !l'idea è di implementare un integratore simplettico del quardo ordine
  subroutine leapfrog(u, v, tau, n)  !devo capire come selzione u o v
    real(dp), intent(inout) :: u(:,:,:) !prende il ruolo del campo
    real(dp), intent(inout) :: v(:,:,:) !prende il ruolo dell'impulso associato al campo
    real(dp), allocatable   :: v1(:,:,:)
    real(dp) :: grad_s1 !la derivata del potenziale hamiltoniano, cioè dell'azione
    real(dp) :: dt      ! passo di integrazione nel tempo di evoluzione del microstato: n*dt=tau
    real(dp) :: mul    ! fattori che entrano nell'operatore simplettico di ordine 4 
    integer  :: dim
    real(dp) :: tau  ! non ho capito in base a cosa lo decido
    integer  :: n,i,j

    dt=(tau/n)_dp
    allocate(v1(L,L,L)) !c'è da aggiustare qualcosa

    do i=1, n
       do j=1,3
          if j==(1.or.3) mul=!valore
          elseif mul=!valore
          end if

          grad_s1=call(grad_s(u+0.5_dp*mul*dt*v))
          v1=v 
          v=v - mul*dt*grad_s1
          u=u + 0.5_dp*mul*dt*(v1+v)
       end do
    end do
    
  end subroutine leapfrog
  
  subroutine metropolis(phi, pi)
    real(dp) :: H_i
    real(dp) :: H_f
    real(dp) :: delta_H
    real(dp), intent(inout)  :: phi(:,:,:)
    real(dp), allocatable :: phi_i(:,:,:)
    real(dp), allocatable :: pi_i(:,:,:)
    real(dp), allocatable :: phi_f(:,:,:)
    real(dp), allocatable :: pi_f(:,:,:)
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
