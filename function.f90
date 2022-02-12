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
       if i==L
          phi(i,:,:)=phi(1,:,:) ! devo stare attento a cose definisco gli array se partono da 0 o da 1, ho considerato da 1
       else
          phi(i,:,:)=phi(i+1,:,:)
       end if       
    case(mu==1)
       if i==L
          phi(:,i,:)=phi(:,1,:)
       else
          phi(:,i,:)=phi(:,i+1,:)
       end if
    case(mu==2)
       if i==L
          phi(:,:,i)=phi(:,:,1)
       else if
          phi(:,:,i)=phi(:,:,i+1)
       end if
    case(mu==3)   !da 3 a 5 si intendono le direzioni negative
       if i==1
          phi(i,:,:)=phi(L,:,:)
       else
          phi(i,:,:)=phi(i-1,:,:)
       end if   
    case(mu==4)
       if i==1
          phi(:,i,:)=phi(:,L,:)
       else
          phi(:,i,:)=phi(:,i-1,:)
       end if
    case(mu==5)
       if i==1
          phi(:,:,i)=phi(:,:,L)
       else
          phi(:,:,i)=phi(:,:,i-1)
       end if
  end subroutine jump
     

  ! S=int d3x -1/2 phi dmu dmu phi + 1/2m^2 phi^2 + lambda/4! phi^4
  ! il laplaciano l'ho discretizzato con la derivata centrale

  subroutine S(phi, k, g, y)
    real(dp), intent(in) :: phi(:,:,:) !compo scalare in 2+1 dimensioni
    real(dp), intent(in) :: k !costante adimonsionale che fa le veci della massa
    real(dp), intent(in) :: g !coupling dell'interazione
    real(dp)  allocatable :: kin(:,:,:) !termine cinetico
    real(dp), intent(out) :: y
    integer :: i, mu, nu, d, dim
    
    dim=size(phi)
    allocate(kin(dim,dim,dim))
    
       
        kin=0
        do mu=0, 2
           nu=mu+3
           kin=call(jump(phi, mu)) + call(jump(phi, nu)
        end do
    
        y=sum(-phi*kin + 0.5_dp*k*phi*phi+ & 
             0.0416666666666666_dp*phi*phi*phi*phi) ! sommo su tutte le posizioni, quindi su tutti gli elementi dell'array
        
  end subroutine S

  subroutine init_pi(L,pi)
    real(dp), allocatable :: x(:,:,:)
    real(dp), intent(out) :: pi(:,:,:)
    integer, intent(in) :: L
    integer :: i

    allocate(x(L,L,L))
    allocate(pi(L,L,L))  
    call random_number(x)

    do i=1,L
       pi(i,:,:)=sqrt(-2*log(1-x(i:,:)))*cos(2*pigreco*(1-x(i,:,:)))
       pi(:,i,:)=sqrt(-2*log(1-x(:,i,:)))*cos(2*pigreco*(1-x(:,i,:)))
       pi(:,:,i)=sqrt(-2*log(1-x(:,:,i)))*cos(2*pigreco*(1-x(:,:,i)))
    end do    
  end subroutine init_pi
  
  subroutine grad_s(phi, y)
    real(dp), intent(in)  :: phi(:,:,:)
    real(dp), allocatable :: kin(:,:,:)
    real(dp), intent(out) :: y(:,:,:)
    
    dim=size(phi)
    allocate(kin(dim,dim,dim))  
    allocate(y(dim,dim,dim))
    
    kin=0
        do mu=0, 2
           nu=mu+3
           kin(:,:,:)=call(jump(phi, mu)) + call(jump(phi, nu)
        end do
        
    y=-kin+ k*phi + 0.166666666666666_dp*phi*phi*phi ! voglio che ogni componente dell'array sia elevato al cubo
  end subroutine grad_s

  function H(phi,pi) return(y)
    real(dp) :: action
    real(dp) :: ham_kin 
    real(dp), intent(in) :: pi(:,:,:)
    real(dp), intent(in) :: phi(:,:,:)
    real(dp), intent(out) :: y
   
    action=call(S(phi))
    ham_kin=sum(pi*pi)  
    y=action +0.5_dp*ham_kin
    
  end function H
  
  
end module
