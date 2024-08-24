import numpy as np
import scipy


def analyticalBlochSolve(A, M_0, nbr_steps, dTau):
    A.astype(np.float64)
    M_0.astype(np.float64)
    propagator = scipy.linalg.expm(A*dTau)


    mTrajList = []
    mTrajList.append(M_0)
    
    for i_step in range(nbr_steps):
        M_0 = np.dot(propagator, M_0)
        mTrajList.append(M_0)

    return np.asarray(mTrajList, dtype=np.float64)