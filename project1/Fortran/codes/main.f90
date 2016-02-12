subroutine func(N, fi, ui, arg, h)
 
 use math, only : dp
 implicit none
 integer,         intent(in)  :: N
 real(kind = dp), intent(in)  :: h 
 real(kind = dp), intent(in),  dimension(N) :: arg
 real(kind = dp), intent(out), dimension(N) :: fi, ui

 fi = 100.0_dp * exp(-10.0_dp * arg) * h**2
 ui = 1.0_dp - ( 1.0_dp - exp(-10.0_dp) ) * arg - exp( -10.0_dp * arg)
 
end subroutine func

program main
   
   use gnufor
   use math, only : &
      dp,           &
      gaussElim,    &
      gaussElim_Sp, &
      LU_lapack
   
   implicit none
      
   integer  :: N , m
   character(len=30) :: fname
   character(len=3)  :: solver
   real(dp) :: ua, ub, h
   real(dp), allocatable, dimension(:,:) :: A
   real(dp), allocatable, dimension(:)   :: Bn, Ba, Xn

   print *, "Enter size of matrix :"
   read  *, N
   print *, "Enter solution regime :"
   read  *, ua, ub
   
   allocate(A(N-1,N-1),  source=0.0_dp)
   allocate(Bn(N-1),     source=0.0_dp)
   allocate(Ba(N-1),     source=0.0_dp)
   allocate(Xn(N-1),     source=0.0_dp)
   
   h = (ub - ua)/N
   print *,h
   do m = 1, N-1
     A(m,m) = 2.0_dp
     if (m+1 <= N-1) then
       A(m+1, m) = -1.0_dp
       A(m, m+1) = -1.0_dp
     endif
     Xn(m) = m * h
   enddo
 
   call func(N-2, Bn, Ba, Xn, h)
#if GE
   solver = "GE"
   call gaussElim(N-1, A, Bn)
#elif GESp
   solver = "GEs"
   call gaussElim_Sp(N-1, A, Bn)
#else
   solver = "LU"
   call LU_lapack(N-1, A, Bn)
#endif
   
   write (fname, "(A13,I5.5,A1)") "comp_ana_num_",N,"_"
   fname = trim(fname)//trim(solver)//".ps"
   print *,fname
   call plot(Xn, Ba, Xn, Bn, filename=fname, terminal='ps', &
             title="Comparison of analytical solution vs numerical solution", &
             xlabel="x", ylabel="U(x)", &
             key1="analytical solution", key2="numerical solution")
   call execute_command_line("ps2pdf14 "//fname, wait=.true.)
   
   deallocate(A)
   deallocate(Bn)
   deallocate(Ba)
   deallocate(Xn)
   
end program main
