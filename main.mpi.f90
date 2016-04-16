 program cg
 use conjugate_gradient_method
 use mpi
 implicit none

 character(len=32) :: arg, mode
 integer :: N = 10, sz,rk
 integer :: i,j,k,it,err, part_sz 
 real*8, dimension(:), allocatable :: A, rowsA
 real*8, dimension(:), allocatable :: b, x, Ap, solve
 real*8, dimension(:), allocatable :: r, p, rnew, pnew
 real*8 :: alpha, beta, rr, rr2, bnorm, pAp, norm
 real*8, parameter :: eps = 1.0E-9
 real*8 :: t1, t2, t3, eps2

 do
    call get_command_argument(i,arg)
    if (len_trim(arg) == 0) exit      
    if (i == 1) read(arg(:),*)N ! dimension of matrix
    if (i == 2) mode = trim(arg)
    i = i + 1
 enddo

 call MPI_INIT(err)
 call MPI_COMM_SIZE(MPI_COMM_WORLD,sz,err)
 call MPI_COMM_RANK(MPI_COMM_WORLD,rk,err)
 part_sz = N / sz

 if (rk == 0) then
    t1 = MPI_WTIME()
    if (mod(N,sz) /= 0)stop 1
    allocate(A(N*N))
    A = 0.0
    allocate(b(N))
    b = 0.0
    allocate(solve(N))
    ! check args
    if (trim(mode) == '-identity') mode = '-test'
    select case (trim(mode))
        case('-test')
            print *, "Test mode is enable: A=E, b=x."   
            do i = 1, N
                A((i-1) * N + i) = 1.0  
                b(i) = i * 1.0
            enddo
        case('-3diag')
            print *, "Test mode is enable: A is 3-diag matrix. x=5."   
            do i = 2, N-1
                A((i-1) * N + i - 1: (i-1) * N + i + 1) = (/1.07,-3.07,1.07/)  
                b(i) = -5.50 
            enddo
            A(1:2) = (/-3.50,1.50/); A(N*N-1:N*N) = (/1.50,-3.50/)
            b(1) = -150.0; b(N) = -160.0
    end select
    
    rr = dotprod(b,b)
    bnorm = sqrt(rr)
    norm = 1.0
 endif
 call MPI_BCAST(norm,1,MPI_REAL8,0,MPI_COMM_WORLD,err)

 allocate(rowsA(part_sz * N))
 allocate(r(part_sz))
 allocate(x(part_sz))
 allocate(p(N)) 
 allocate(Ap(part_sz))
 allocate(rnew(part_sz))
 allocate(pnew(part_sz)) 
 if (rk == 0) p = b
 call MPI_BCAST(p,N,MPI_REAL8,0,MPI_COMM_WORLD,err)
 call MPI_SCATTER(A, part_sz * N, MPI_REAL8, rowsA, part_sz * N, MPI_REAL8, 0, MPI_COMM_WORLD, err)
 call MPI_SCATTER(b, part_sz, MPI_REAL8, r, part_sz, MPI_REAL8, 0, MPI_COMM_WORLD, err)
 if (rk == 0) deallocate(A)

 do while (norm > eps)
    Ap = multmv(rowsA,p)
    pAp = dotprod(Ap,p(rk * part_sz + 1:(rk+1) * part_sz))    
    call MPI_REDUCE(pAp,alpha,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,err)
    if (rk == 0) alpha = rr / alpha
    call MPI_BCAST(alpha,1,MPI_REAL8,0,MPI_COMM_WORLD,err)
    
    rnew = r - alpha * Ap
    rr2 = dotprod(rnew, rnew)   
    call MPI_REDUCE(rr2,beta,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,err)
    if (rk == 0) beta = beta / rr
    call MPI_BCAST(beta,1,MPI_REAL8,0,MPI_COMM_WORLD,err)
    
    x = x + alpha * p(rk * part_sz + 1:(rk+1) * part_sz)  
    pnew = rnew + beta * p(rk * part_sz + 1:(rk+1) * part_sz) 
    r = rnew
    call MPI_ALLGATHER(pnew,part_sz,MPI_REAL8,p,part_sz,MPI_REAL8,MPI_COMM_WORLD,err)
    
    if (rk == 0) then
        rr = beta * rr
        norm = sqrt(rr) / bnorm
        it = it + 1
        !print *, "On ", it, " norm = ", norm
    endif
    call MPI_BCAST(norm,1,MPI_REAL8,0,MPI_COMM_WORLD,err)
 enddo
 call MPI_GATHER(x,part_sz,MPI_REAL8,solve,part_sz,MPI_REAL8,0,MPI_COMM_WORLD,err)

 if (rk == 0) then
    t2 = MPI_WTIME()
    print "('||r||/||b|| = ',ES10.2, ' ?? ',i5 ' ????????.')", sqrt(rr) / bnorm, it
    print '("????? ",f10.5,"? ??? ",i7,"-?????? ?? ", i3 " ?????????.")', t2 - t1, N, sz
    print '("CHECK of solve (method): ",L)',check_solve(a,solve,b,eps)
    print *,""
    if(N <= 10)print *,"x: ", solve
    deallocate(A)
 endif
 call MPI_FINALIZE(err)
 end program cg
