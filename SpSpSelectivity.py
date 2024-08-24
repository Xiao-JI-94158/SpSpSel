import numpy as np
import scipy


def constructBlochMatrix(deltaOmegaZ:float=0.0, omega1X:float=0.0, omega1Y:float=0.0, R_1:float=0.0, R_2:float=0.0):
    BlochMat = np.array([[-R_2, deltaOmegaZ, -omega1Y, 0],[-deltaOmegaZ, -R_2,omega1X,0],[omega1Y, -omega1X, R_1, R_1],[0,0,0,0]])
    return np.asarray(BlochMat, dtype=np.float64)

def solveAnalyticalBlochStepPropagator(A:np.array, Tau:float):
    A.astype(np.float64)
    return np.asarray(scipy.linalg.expm(A*Tau), dtype=np.float64)

def solveAnalyticalBlochTraj(A:np.array, M_0:np.array, nbr_steps:int, dTau:float):
    A.astype(np.float64)
    M_0.astype(np.float64)
    propagator = solveAnalyticalBlochStepPropagator(A, dTau)

    mTrajList = []
    mTrajList.append(M_0)
    
    for i_step in range(nbr_steps):
        M_0 = np.dot(propagator, M_0)
        mTrajList.append(M_0)

    return np.asarray(mTrajList, dtype=np.float64)