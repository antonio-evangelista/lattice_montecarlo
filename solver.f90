module solver
  use function
  implicit none

  public :: metropolis
  public :: leapfrog

contains
  
  subroutine leapfrog(u, v, tau, n)  !devo capire come selzione u o v
    real(dp), intent(inout) :: u(:,:,:) !prende il ruolo del campo
    real(dp), intent(inout) :: v(:,:,:) !prende il ruolo dell'impulso associato al campo
    real(dp), allocatable   :: v1(:,:,:)
    real(dp) :: grad_s1 !la derivata del potenziale hamiltoniano, cioè dell'azione
    real(dp) :: dt      ! passo di integrazione nel tempo di evoluzione del microstato: n*dt=tau
    real(dp) :: mul    ! fattori che entrano nell'operatore simplettico di ordine 4 
    integer  :: dim
    real(dp) :: tau  ! inizialmente dovrebbe essere un parametro libero che però se faccio varie simulazioni
                     ! in base a quanto miglioro l'autocorrelazione, dt=tau è un pò poco mentre già un paio di
                     ! iterazioni dovrebbero bastare
    integer  :: n,i,j, L

    dt=(tau/n)_dp
    L=size(v,dim=1)
    allocate(v1(L,L,L)) 

    do i=1, n
       do j=1,3
          if j==(1.or.3) mul=1.351207191959658_dp
          elseif mul=2.702414383919316_dp 
          end if

          grad_s1=call(grad_s(u+0.5_dp*mul*dt*v))
          v1=v 
          v=v - mul*dt*grad_s1
          u=u + 0.5_dp*mul*dt*(v1+v)
       end do
    end do

    deallocate(v1) !occupo solo memoria altrimenti
    
  end subroutine leapfrog
  
  subroutine metropolis(phi, pi) !poichè chiamo una funzione che richiede anche altri parametri devo specificarli anche qui?
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
    integer   :: L

    L=size(phi,dim=1)
    phi_i=phi !forse devo prima allocare
    pi_i=pi
    allocate(phi_f(L,L,L))
    allocate(pi_f(L,L,L))

    phi_f =call(leapfrog(phi_i, pi_i, tau, n))
    pi_f  =call(leapfrog(phi_i, pi_i, tau, n))
    H_i   =call(H(phi_i,pi_i))
    H_f   =call(H(phi_f,pi_f))
    delta_H =H_f-H_i

    prob=min(1, exp(-delta_H))

  if prob==1
    phi=phi_f
  else 
    call random_number(r)
    if r < prob
       phi=phi_f
    else
       phi=phi_i
    end if
  end if

 deallocate(phi_i)
 deallocate(phi_f)
 deallocate(pi_i)
 deallocate(pi_f)

end subroutine metropolis

end module solver
