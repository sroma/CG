program cg
use conjugate_gradient_method
use omp_lib
implicit none

 character(len=32) :: arg, mode
integer :: N = 10, NUMPROC = 1
integer :: i,j,k,it 
real*8, dimension(:,:), allocatable :: A
real*8, dimension(:), allocatable :: b, x, Ap
real*8, dimension(:,:), allocatable :: r, p
real*8 :: alpha, beta, rr, rr2,  bnorm, pAp
real*8, parameter :: eps = 1.0E-9
real*8 :: t1, t2, t3

! command arguments
do
    call get_command_argument(i,arg)
    if (len_trim(arg) == 0) exit      
    if (i == 1) read(arg(:),*)N ! dimension of matrix
    if (i == 2) read(arg(:),*)NUMPROC  ! num of processes
    if (i == 3) mode = trim(arg)
    i = i + 1
enddo

print '("Dimension: ",i7)',N
print '("Num of processes: ",i3)', NUMPROC
print '("Accuracy: ", ES10.2)', eps

t1 = omp_get_wtime()

allocate(A(n,n))
A = 0.0
allocate(b(n))
allocate(x(n))
b = 0.0
allocate(r(n,2))
allocate(p(n,2))
allocate(Ap(n))

 call omp_set_num_threads(NUMPROC)
! check args
if (trim(mode) == '-identity') mode = '-test'
select case (trim(mode))
    case('-test')
       print *, "Test mode is enable: A=E, b=x."   
       !$OMP PARALLEL DO PRIVATE(i), SHARED(a,b,x)
       do i = 1, n
            A(i,i) = 1.0  
            b(i) = i * 1.0
       enddo
    case('-3diag')
       print *, "Test mode is enable: A is 3-diag matrix"   
       !$OMP PARALLEL DO PRIVATE(i), SHARED(a,b,x)
       do i = 2, n-1
            A(i-1:i+1,i) = (/1.0,-3.0,1.0/)  
            b(i) = -5.0 
         enddo
       A(1:2,1) = (/-3.0,1.0/); A(N-1:N,N) = (/1.0,-3.0/)
       b(1) = -10.0; b(n) = -10.0
end select

x = 0.0
r = 0.0
p = 0.0
r(:,1) = b
p(:,1) = r(:,1)
rr = dotprod(b,b)
bnorm = sqrt(rr)

t2 = omp_get_wtime()
! method
it = 0
do while (sqrt(rr) / bnorm > eps)
    Ap = multmv(A,p(:,1))
    pAp =  dotprod(Ap,p(:,1))
    alpha = rr / pAp   
    r(:,2) = r(:,1) - alpha * Ap 
    rr2 = dotprod(r(:,2),r(:,2)) 
    beta = rr2 / rr
    p(:,2) = r(:,2) + beta * p(:,1) 
    x = x + alpha * p(:,1)
    rr = rr2
    p(:,1) = p(:,2)
    r(:,1) = r(:,2)
    it = it + 1
enddo
t3 = omp_get_wtime()
print "('||r||/||b|| = ',ES10.2, ' on ',i5 ' iters.')", sqrt(rr) / bnorm, it
print '("Method time ",f10.5,"s for ",i7,"-size on ", i3 " threads.")', t3 - t2, N, NUMPROC
print '("Initial time, s: ",f10.5)',t2 - t1
print '("CHECK of solve (method): ",L)',check_solve(a,x,b,eps)
print *,""
if(N <= 10)print *,"x: ", x

deallocate(a)
deallocate(b)
deallocate(x)
deallocate(r)
deallocate(p)
deallocate(Ap)
end program cg
