import numpy as np

''' Solve linear system Ax=b,
    A is symmetric and positive-definite matrix,
    b is known column vector,
    x is unknown column vector,
    eps - accuracy.
    Return: 
        x - solution,
        it - number of iteration,
        rnorm * bnorm - accuracy.
'''
def CG_method(A, b, eps = 1.e-9):
    N = A.shape[0]
    assert (N == A.shape[1]), "Matrix A must be square"
    assert (N == b.shape[0]), "Dimension of vector b must be equal to dimension of matrix A"

    x = np.empty_like (b) 
    x[:] = 0.
    r = np.empty_like (b)
    p = np.empty_like (b)
    r[:] = b
    p[:] = b
    bnorm2 = float(np.dot(b.T,b))
    rnorm2 = bnorm2
    bnorm2 = 1. / bnorm2
    eps2 = eps * eps

    # main cycle
    it = 0 
    while (rnorm2 * bnorm2 > eps2):
        Ap = A * p
        alpha = rnorm2 / float(np.dot(p.T, Ap))
        x += alpha * p
        r -= alpha * Ap
        r2norm2 = float(np.dot(r.T,r))
        beta = r2norm2 / rnorm2
        p = r + beta * p
        rnorm2 = r2norm2
        it += 1
    return x, it, np.sqrt(rnorm2 * bnorm2)

'''
    Test run: Solve 3-diagonal matrix
'''
if __name__ == "__main__":
   
    N = 100
    A = np.zeros((N,N))
    b = np.zeros((N,1))
    b[0] = 6.6
    b[N-1] = 6.6
    A[0][0] = 5.
    A[0][1] = 1.
    A[N-1][N-2] = 1.
    A[N-1][N-1] = 5.
    for i in xrange(N):
        if (i != 0 and i != N-1): 
            A[i][i] = 5.
            A[i][i-1] = 1.
            A[i][i+1] = 1.
            b[i] = 7.7

    x,it,eps = CG_method(np.matrix(A), np.matrix(b))
    print "||r||/||b|| is {0} on {1:d} iteration".format(eps, it)
    if (x.size <= 10): print "x = {0}".format(x.T)
