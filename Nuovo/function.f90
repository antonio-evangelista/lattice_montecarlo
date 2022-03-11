module function
  use precision
  implicit none

  public :: S
  public :: H
  public :: init_pi
  public :: func
  
  real(dp), parameter :: pigreco = 3.141592653589793
  !real(dp), public :: k,g

  interface
     function func(u, v) result (y)
       use precision
       real(dp) , intent(in)  :: u(:,:,:)
       real(dp) , intent(in)  :: v(:,:,:)
       real(dp) :: y
     end function 
  end interface

contains


  subroutine S(u, y, k, g)
    real(dp), intent(in) :: u(0:,0:,0:) !compo scalare in 2+1 dimensioni
    real(dp), intent(out) :: y
    real(dp), intent(in) :: k, g
    integer ::  L, ii, jj, kk, ii_m, ii_p, jj_m, jj_p, kk_m, kk_p
    
    L=size(u,dim=1)

    ! sommo su tutte le posizioni, quindi su tutti gli elementi dell'array
    y=sum(0.5_dp*k*u*u+ 0.0416666666666666_dp*g*u*u*u*u) 

    do ii = 0, L-1
      ii_m = mod(ii-1+L,L)
      ii_p = mod(ii+1+L,L)
      do jj = 0, L-1
        jj_m = mod(jj-1+L,L)
        jj_p = mod(jj+1+L,L)
        do kk = 0, L-1
          kk_m = mod(kk-1+L,L)
          kk_p = mod(kk+1+L,L)
 
          y = y - 0.5_dp*u(ii,jj,kk)*(u(ii_p,jj,kk)+u(ii_m,jj,kk) + &
                              &       u(ii,jj_p,kk)+u(ii,jj_m,kk) + & 
                              &       u(ii,jj,kk_p)+u(ii,jj,kk_m)) 
           
        end do
      end do  
   end do

  end subroutine S

  subroutine init_pi(L,v)
    integer, intent(in) :: L
    real(dp), intent(inout) :: v(0:,0:,0:)
    real(dp), allocatable :: x1(:,:,:)
    real(dp), allocatable :: x2(:,:,:)

    allocate(x1(0:L-1,0:L-1,0:L-1))
    allocate(x2(0:L-1,0:L-1,0:L-1))
    call random_number(x1)
    call random_number(x2)

    v = sqrt(-2.0_dp*log(1.0_dp-x1))*cos(2.0_dp*pigreco*(1.0_dp-x2))

    deallocate(x1)
    deallocate(x2)
    
  end subroutine init_pi
  
  subroutine grad_s(u, y, k, g)
    real(dp), intent(in)  :: u(0:,0:,0:)
    real(dp), intent(inout) :: y(0:,0:,0:)
    real(dp), intent(in) :: k, g
    integer ::  L, ii, jj, kk, ii_m, ii_p, jj_m, jj_p, kk_m, kk_p
    
    L=size(u, dim=1)

    y = k*u + 0.166666666666666_dp*g*u*u*u 
    ! voglio che ogni componente dell'array sia elevato al cubo

    do ii = 0, L-1
      ii_m = mod(ii-1+L,L)
      ii_p = mod(ii+1+L,L)
      do jj = 0, L-1
        jj_m = mod(jj-1+L,L)
        jj_p = mod(jj+1+L,L)
        do kk = 0, L-1
          kk_m = mod(kk-1+L,L)
          kk_p = mod(kk+1+L,L)

          y(ii,jj,kk) = y(ii,jj,kk) - (u(ii_p,jj,kk)+u(ii_m,jj,kk) + &
                                    &  u(ii,jj_p,kk)+u(ii,jj_m,kk) + & 
                                    &  u(ii,jj,kk_p)+u(ii,jj,kk_m)) 
          
        end do
      end do  
   end do
    
  end subroutine grad_s

  function H(u,v, k, g) result(y)
    real(dp), intent(in) :: v(0:,0:,0:)
    real(dp), intent(in) :: u(0:,0:,0:)
    real(dp), intent(in) :: k, g
    real(dp) :: y ! , intent(out)
    real(dp) :: azione
    real(dp) :: ham_kin
       
    call S(u, azione, k, g)
    ham_kin=sum(v*v)  
    y=azione +0.5_dp*ham_kin
    
  end function H
  
  
end module
