module function
  implicit none

  public :: S
  public :: H

  real(dp), allocatable :: phi(:)
  real(dp), parameter :: pigreco= 3.14159 26535 89793

  !forse interface è inutile 
  interface
     function H(phi, pi) result(y)
       real(dp) , allocatable  :: phi(:)
       real(dp) , intent(in)   :: pi
       real(dp) , intent(out)  :: y(t)
     end function H
  end interface

  contains

  subroutine S(phi)
    real(dp), allocatable :: phi(:)
    real(dp) :: kin
    real(dp), intent(out) :: y
    integer :: i, d, n

    n=size(phi)
    
    do i=0, 2
       !kin+=
    end do
    
    do i=0, n
       y+= kin ! + blablabla
    end do
  end subroutine S

  subroutine init_pi(phi,pi)
    real(dp), intent(in) :: phi(:)
    real(dp), allocatable :: x(:)
    real(dp), intent(out) :: pi(:)
    integer :: n, i

    n=size(phi)
    allocate(x(n))
    allocate(pi(n))  
    call random_number(x)

    do i=1,n
       pi(i)=sqrt(-2*log(1-x(i)))*cos(2*pigreco*(1-x(i)))
    end do    
  end subroutine init_pi

  function H(phi,pi) return(y)
    real(dp), intent(in) :: action
    real(dp), intent(in) :: pi(:)
    real(dp), intent(in) :: phi(:)
    real(dp), intent(out) :: y
    integer :: i, n

    n=size(pi)

    action=call(S(phi))
    do i=1,n
       y+=pi(i)*pi(i)
    end do
    y+=action
  end function H
  
  
end module
