program main

   use math, only : &
      dp,           &
      gaussElim,    &
      gaussElim_Sp, &
      LU_lapack
   
   implicit none
      
   integer  :: N = 4, m
   real(dp), allocatable, dimension(:,:) :: A
   real(dp), allocatable, dimension(:)   :: B
   
   allocate(A(N,N),  source=0.0_dp)
   allocate(B(N),    source=0.0_dp)
   
   do m = 1, N
     A(m,m) = -2.0_dp
     if (m+1 <= N) then
       A(m+1, m) = 1.0_dp
       A(m, m+1) = 1.0_dp
     endif
   enddo
   
   B(1) = -20.0_dp
   B(N) = -70.0_dp

#if GE
   call gaussElim(N, A, B)
#elif GESp
   call gaussElim_Sp(N, A, B)
#else
   call LU_lapack(N, A, B)
#endif

   write (*, '(a15, (4f12.5, 1x))') "Solution :",B

   deallocate(A)
   deallocate(B)

end program main
