module math
  
implicit none

integer, parameter, public :: dp = kind(1.0D0)
integer                    :: i, j, k
real(dp)                   :: R, S

contains
  
  subroutine gaussElim(N, A, B)
    implicit none
    ! Define variables and arrays
    integer, intent(in)                     :: N
    real(dp), dimension(N,N), intent(inout) :: A  
    real(dp), dimension(N),   intent(inout) :: B
    real(dp), allocatable,    dimension(:)  :: x  
    ! Initialize variables
    allocate(x(N), source=0.0_dp)
    ! Doing Calculation
    write (*, *) "Using Gaussian Elimination (Standard) function"
    ! Forward Substitution
    k = 2
    do while (k <= N)
     do i = k, N
       R = A(i,k-1)/A(k-1,k-1)
       do j = k-1, N
         A(j,i) = A(j,i) - R * A(j,k-1)
       enddo
       B(i) = B(i) - R*B(k-1)
     enddo
     k = k + 1
    enddo
    ! Backward Substitution
    x(N) = B(N)/A(N,N)
    k = N-1
    do while (k > 0)
     S = 0.0_dp
     do i = k, N-1
       S = S + A(i+1,k)*x(i+1)
     enddo
     x(k) = ( B(k) - S ) / A(k,k)
     k = k - 1
    enddo
    B = x
    ! Deallocate additional memory used
    deallocate(x)    
  end subroutine gaussElim
  
  subroutine gaussElim_Sp(N, A, B)
    implicit none
    ! Define variables and arrays
    integer, intent(in)                     :: N
    real(dp), dimension(N,N), intent(inout) :: A  
    real(dp), dimension(N),   intent(inout) :: B
    real(dp), allocatable,    dimension(:)  :: D, x
    real(dp), allocatable,    dimension(:)  :: E, F  
    ! Initialize variables
    allocate(D(N),   source=0.0_dp)
    allocate(x(N),   source=0.0_dp)
    allocate(E(N-1), source=0.0_dp)
    allocate(F(N-1), source=0.0_dp)      
    ! Convert A(N,N) to 1D diagonal elements i.e D(N), E(N), F(N)
    do i = 1, N
      D(i) = A(i,i)
      if (i+1 <= N) then
        E(i) = A(i, i+1)
        F(i) = A(i+1, i)
      endif
    enddo
    ! Doing Calculation
    write (*, *) "Using Gaussian Elimination (Tridiagonal) function"
    ! Forward Substitution
    do i = 2, N
      R = ( F(i-1) / D(i-1) ) * E(i-1)
      D(i) = D(i) - R * E(i-1)
      B(i) = B(i) - R * B(i-1)
    enddo
    ! Backward Substitution
    x(N) = B(N)/D(N)
    k = N-1
    do while (k > 0)
      x(k) = ( B(k) - E(k)*x(k+1) ) / D(k)
      k = k - 1
    enddo
    B = x
    ! Deallocate additional memory used
    deallocate(D, x)
    deallocate(E, F) 
  end subroutine gaussElim_Sp
  
  subroutine LU_lapack(N, A, B)
    implicit none
    ! Define variables and arrays
    character(len=1)                        :: TRANS
    integer                                 :: LDA, LDB, NRHS, INFO
    real(dp), allocatable,    dimension(:)  :: x
    integer, intent(in)                     :: N 
    integer,  dimension(N)                  :: IPIV
    real(dp), dimension(N),   intent(inout) :: B
    real(dp), dimension(N,N), intent(inout) :: A
    ! Initialize lapack variables
    TRANS = 'N'
    LDA   =  N
    LDB   =  N
    NRHS  =  1
    INFO  =  0
    IPIV  =  0
    ! Doing Calculation
    write (*, *) "Using LU decomposition (Lapack) function"
    ! Call factorisation function LU function
    call dgetrf(N, N, A, LDA, IPIV, INFO)
    ! Call matrix solver from P*L*D matrix
    call dgetrs(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
  end subroutine LU_lapack

end module math