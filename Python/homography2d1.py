import numpy as np
from scipy.optimize import least_squares
from normalise2dpts import normalise2dpts
from simon_H import simon_H

# HOMOGRAPHY2D - computes 2D homography
# Usage: H = homography2d(x1, x2)
#        H = homography2d(x)

def homography2d(*args):
    x1, x2 = checkargs(args)
    M = x1
    m = x2

    # Normalize each set of points
    x1, T1 = normalise2dpts(x1)
    x2, T2 = normalise2dpts(x2)

    # Estimation of the H between the model plane and its image
    Npts = x1.shape[1]
    A = np.zeros((3 * Npts, 9))
    O = np.array([0, 0, 0])

    for n in range(Npts):
        X = x1[:, n].T
        x = x2[0, n]
        y = x2[1, n]
        w = x2[2, n]
        A[3 * n] = np.concatenate([O, -w * X, y * X])
        A[3 * n + 1] = np.concatenate([w * X, O, -x * X])
        A[3 * n + 2] = np.concatenate([-y * X, x * X, O])

    #_, _, V = np.linalg.svd(A)
    U, D, Vt = np.linalg.svd(A)
    V = Vt.T

    # Extract homography
    #H1 = V[-1, :].reshape(3, 3).T
    H1 = V[:, 8].reshape(3, 3)

    # Denormalize
    H2 = np.linalg.inv(T2) @ H1 @ T1
    H = H2 / H2[2, 2]

    # Maximun likelihood estimation for the H
    options = {'method': 'lm'}
    res = least_squares(simon_H, H.flatten(), args=(m, M), method='lm')
    H = res.x.reshape(3, 3)
    H = H / H[2, 2]

    return H


def checkargs(arg):
    if len(arg) == 2:
        x1, x2 = arg
        if x1.shape != x2.shape:
            raise ValueError('x1 and x2 must have the same size')
        elif x1.shape[0] != 3:
            raise ValueError('x1 and x2 must be 3xN')
    elif len(arg) == 1:
        if arg[0].shape[0] != 6:
            raise ValueError('Single argument x must be 6xN')
        else:
            x1 = arg[0][:3, :]
            x2 = arg[0][3:, :]
    else:
        raise ValueError('Wrong number of arguments supplied')

    return x1, x2



