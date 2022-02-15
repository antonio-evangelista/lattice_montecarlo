module function
  use precision
  implicit none

  public :: S
  public :: H
  public :: jump
  public :: init_pi
  public :: func
  public :: kin
  
  real(dp), parameter :: pigreco = 3.141592653589793  

  interface
     function func(u, v) result (y)
       use precision
       real(dp) , intent(in)  :: u(:,:,:)
       real(dp) , intent(in)  :: v(:,:,:)
       real(dp) :: y
     end function 
  end interface

contains

  subroutine jump(u, mu)
    real(dp), intent(inout) :: u(:,:,:)
    integer, intent(in) :: mu
    integer :: L, i

    L=size(u,dim=1)

    if (mu==0) then
       if (i==L) then
          u(i,:,:)=u(1,:,:) ! devo stare attento a cose definisco gli array se partono da 0 o da 1, ho considerato da 1
       else
          u(i,:,:)=u(i+1,:,:)
       end if       
    elseif (mu==1) then
       if (i==L) then
          u(:,i,:)=u(:,1,:)
       else
          u(:,i,:)=u(:,i+1,:)
       end if
    elseif (mu==2) then
       if (i==L) then
          u(:,:,i)=u(:,:,1)
       else
          u(:,:,i)=u(:,:,i+1)
       end if
    elseif (mu==3) then   !da 3 a 5 si intendono le direzioni negative
       if (i==1) then
          u(i,:,:)=u(L,:,:)
       else
          u(i,:,:)=u(i-1,:,:)
       end if   
    elseif (mu==4) then
       if (i==1) then
          u(:,i,:)=u(:,L,:)
       else
          u(:,i,:)=u(:,i-1,:)
       end if
    elseif (mu==5) then
       if (i==1) then
          u(:,:,i)=u(:,:,L)
       else
          u(:,:,i)=u(:,:,i-1)
       end if
    end if
    
  end subroutine jump
     

  ! S=int d3x -1/2 phi dmu dmu phi + 1/2m^2 phi^2 + lambda/4! phi^4
  ! il laplaciano l'ho discretizzato con la derivata centrale

  subroutine kin(u,cin)
    real(dp), intent(in) :: u(:,:,:) !compo scalare in 2+1 dimensioni
    real(dp), allocatable, intent(out) :: cin(:,:,:) !termine cinetico
    real(dp), allocatable :: phi1(:,:,:)
    real(dp), allocatable :: phi2(:,:,:)
    integer ::  mu, nu, L

    L=size(u,dim=1)
    allocate(cin(L,L,L))
    allocate(phi1(L,L,L))
    allocate(phi2(L,L,L))

    phi1=u
    phi2=u
    cin=0
    
        do mu=0, 2
           nu=mu+3
           
           call jump(phi1, mu)
           call jump(phi2, nu)
           cin=cin + phi1+phi2
        end do
        
    deallocate(phi1)
    deallocate(phi2)
    
  end subroutine kin
  

  subroutine S(u, k, g, y)
    real(dp), intent(in) :: u(:,:,:) !compo scalare in 2+1 dimensioni
    real(dp), intent(in) :: k !costante adimonsionale che fa le veci della massa
    real(dp), intent(in) :: g !coupling dell'interazione
    real(dp), intent(out) :: y
    real(dp), allocatable :: cin(:,:,:)
    integer ::  L
    
    L=size(u,dim=1)
    allocate(cin(L,L,L))

    call kin(u,cin)    
    y=sum(-u*cin + 0.5_dp*k*u*u+ & 
             0.0416666666666666_dp*u*u*u*u) ! sommo su tutte le posizioni, quindi su tutti gli elementi dell'array

   deallocate(cin)
        
  end subroutine S

  subroutine init_pi(L,v)
    real(dp), allocatable :: x(:,:,:)
    real(dp), allocatable, intent(out) :: v(:,:,:)
    integer, intent(in) :: L
    integer :: i

    allocate(x(L,L,L))
    allocate(v(L,L,L))  
    call random_number(x)

    do i=1,L
       v(i,:,:)=sqrt(-2*log(1-x(i,:,:)))*cos(2*pigreco*(1-x(i,:,:)))
       v(:,i,:)=sqrt(-2*log(1-x(:,i,:)))*cos(2*pigreco*(1-x(:,i,:)))
       v(:,:,i)=sqrt(-2*log(1-x(:,:,i)))*cos(2*pigreco*(1-x(:,:,i)))
    end do

    deallocate(x)
    
  end subroutine init_pi
  
  subroutine grad_s(u, k, g, y)
    real(dp), intent(in)  :: u(:,:,:)
    real(dp), intent(in)  :: k, g    
    real(dp), allocatable, intent(out) :: y(:,:,:)
    real(dp), allocatable :: cin(:,:,:)
    integer :: L
    
    L=size(u)
    allocate(cin(L,L,L))  
    allocate(y(L,L,L))

    call kin(u,cin)
        
    y=-cin+ k*u + 0.166666666666666_dp*g*u*u*u ! voglio che ogni componente dell'array sia elevato al cubo

    deallocate(cin)
    
  end subroutine grad_s

  function H(u,v) result(y)
    real(dp) :: azione
    real(dp) :: ham_kin 
    real(dp), intent(in) :: v(:,:,:)
    real(dp), intent(in) :: u(:,:,:)
    real(dp) :: y ! , intent(out)
    real(dp) :: g
    real(dp) :: k
   
    call S(u, k, g, azione)
    ham_kin=sum(v*v)  
    y=azione +0.5_dp*ham_kin
    
  end function H
  
  
end module
