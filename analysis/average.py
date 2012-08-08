import numpy as np
import os
import MD.base.points as points


def average_position(M):
    A = np.zeros((M.shape[1],3))
    for i in range(M.shape[1]):
        count = 1
        p = np.array([0.0,0.0,0.0])
        p += M[0][i]
        for j in range(1,M.shape[0]):
            if points.dist_np(M[0][i],M[j][i])[0] <  5:
                p += M[j][i]
                count += 1
        p = p/count
        A[i] = p
    return A
