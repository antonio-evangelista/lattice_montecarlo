module function
  implicit none

  public :: S
  public :: H
  public :: jump

  real(dp), allocatable :: phi(:,:,:)
  real(dp), parameter :: pigreco= 3.14159 26535 89793

  interface
     function H(phi, pi) result(y)
       real(dp) , allocatable  :: phi(:,:,:)
       real(dp) , allocatable  :: pi(:,:,:)
       real(dp) , intent(out)  :: y
     end function H
  end interface

contains

  subroutine jump(phi, mu)
    real(dp), intent(inout) :: phi(:,:,:)
    integer :: mu,L

    case(mu==0)
       phi(i,:,:)=phi(i+1,:,:)
    case(mu==1)
       phi(:,i,:)=phi(:,i+1,:)
    case(mu==2)
       phi(:,:,i)=phi(:,:,i+1)
    case(mu==3)   !da 3 a 5 si intendono le direzioni negative
       phi(i,:,:)=phi(i+L-1,:,:)
    case(mu==4)
       phi(:,i,:)=phi(:,i+L-1,:)
    case(mu==5)
       phi(:,:,i)=phi(:,:,i+L-1)

  end subroutine jump
     

  ! S=int d3x -phi dmu dmu phi + 1/2m^2 phi^2 + lambda/4! phi^4
  ! il laplaciano l'ho discretizzato con la derivata centrale

  subroutine S(phi, k, g, y)
    real(dp), intent(in) :: phi(:,:,:) !compo scalare in 2+1 dimensioni
    real(dp), intent(in) :: k !costante adimonsionale che fa le veci della massa
    real(dp), intent(in) :: g !coupling dell'interazione
    real(dp)  allocatable :: kin(:,:,:) !termine cinetico
    real(dp), intent(out) :: y
    integer :: i, mu, nu, d, n

    n=size(phi)
    
    do i=0, n !c'è da ripensarlo un attimo in modo vettoriale
       
        kin=0
        do mu=0, 2
           nu=mu+3
           kin+=call(jump(phi(i), mu)) + call(jump(phi(i), nu)
        end do
    
        y+= -phi(i)*kin + 0.5_dp*k*phi(i)*phi(i)+ &
             0.0416666666666666_dp**phi(i)*phi(i)*phi(i)*phi(i)
    end do
  end subroutine S

  subroutine init_pi(phi,pi)
    real(dp), intent(in) :: phi(:,:,:)
    real(dp), allocatable :: x(:,:,:)
    real(dp), intent(out) :: pi(:,:,:)
    integer :: n, i

    n=size(phi)
    allocate(x(n))
    allocate(pi(n))  
    call random_number(x)

    do i=1,n
       pi(i,:,:)=sqrt(-2*log(1-x(i:,:)))*cos(2*pigreco*(1-x(i,:,:)))
       pi(:,i,:)=sqrt(-2*log(1-x(:,i,:)))*cos(2*pigreco*(1-x(:,i,:)))
       pi(:,:,i)=sqrt(-2*log(1-x(:,:,i)))*cos(2*pigreco*(1-x(:,:,i)))
    end do    
  end subroutine init_pi
  
  subroutine grad_s(phi, y)
    real(dp), intent(in) :: phi(:,:,:)
    real(dp)  :: kin
    real(dp), intent(out) :: y(:,:,:)
    
     kin=0
        do mu=0, 2
           nu=mu+3
           kin(:,:,:)=call(jump(phi, mu)) + call(jump(phi, nu)
        end do
    y=-kin+ k*phi + 0.166666666666666_dp*phi*phi*phi !devo stare attento perche voglio che le componenti siano elevate al cubo e non siano fatti prodotti scalari
  end subroutine grad_s

  function H(phi,pi) return(y)
    real(dp), intent(in) :: action
    real(dp),            :: pi(:,:,:)
    real(dp), intent(in) :: phi(:,:,:)
    real(dp), intent(out) :: y
    integer :: i, n

    n=size(pi)

    action=call(S(phi))
    y=action + pi*pi
    
  end function H
  
  
end module
