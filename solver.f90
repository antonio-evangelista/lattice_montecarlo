module solver
  use function
  implicit none

  public :: metropolis
  public :: evolution

contains

  subroutine metropolis(phi_i, pi_i)
    real(dp), intent(inout) :: H_i
    real(dp), intent(inout) :: H_f
    real(dp)                :: delta_H
    real(dp), allocatable   :: phi_i(:)
    real(dp), allocatable   :: pi_i(:)
    real(dp), allocatable   :: phi_f(:)
    real(dp), allocatable   :: pi_f(:)
    real(dp), allocatable   :: phi(:)
    real(dp)  :: prob
    real(dp)  :: r
    integer   :: n

    n=size(phi)
    phi_i=phi
    allocate(pi_i)
    allocate(phi_f(n))
    allocate(pi_f(n))

    phi_f =call(evolution(phi_i, pi_i,1))
    pi_f  =call(evolution(phi_i, pi_i,2))
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
