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

def solveBlochRfWaveformPropagatorList(rfwaveformX:np.array, rfwaveformY:np.array, Tau:float, offsetZ:float=0, R_1:float=0, R_2:float=0):
    nSteps = len(rfwaveformX)

    if not (nSteps==len(rfwaveformY)):
        raise ValueError(f"rfwaveformX ({len(rfwaveformX)}) and rfwaveformY ({len(rfwaveformY)}) are not of equal length")
    
    dTau = Tau / nSteps

    propagatorList = []

    for idx, rfVals in enumerate(zip(rfwaveformX,rfwaveformY)):
        BlochMat = np.array([[-R_2, offsetZ, -rfVals[1], 0],[-offsetZ, -R_2,rfVals[0],0],[rfVals[1], -rfVals[0], R_1, R_1],[0,0,0,0]])
        propagator = scipy.linalg.expm(BlochMat*dTau)
        propagatorList.append(propagator)


    return propagatorList

def solveBlochRfWaveformTraj(initialState:np.array, rfwaveformX:np.array, rfwaveformY:np.array, Tau:float, offsetZ:float=0, R_1:float=0, R_2:float=0):

    nSteps = len(rfwaveformX)

    if not (nSteps==len(rfwaveformY)):
        raise ValueError(f"rfwaveformX ({len(rfwaveformX)}) and rfwaveformY ({len(rfwaveformY)}) are not of equal length")
    
    dTau = Tau / nSteps

    mTrajList = []
    mTrajList.append(initialState)
    M_0 = initialState
    
    for idx, rfVals in enumerate(zip(rfwaveformX,rfwaveformY)):
        BlochMat = np.array([[-R_2, offsetZ, -rfVals[1], 0],[-offsetZ, -R_2,rfVals[0],0],[rfVals[1], -rfVals[0], R_1, R_1],[0,0,0,0]])
        propagator = scipy.linalg.expm(BlochMat*dTau)
        M_0 = np.dot(propagator, M_0)
        mTrajList.append(M_0)

    return np.asarray(mTrajList, dtype=np.float64)

def solveBlochRfWaveformOverallPropagator(rfwaveformX:np.array, rfwaveformY:np.array, Tau:float, offsetZ:float=0, R_1:float=0, R_2:float=0):
    nSteps = len(rfwaveformX)

    if not (nSteps==len(rfwaveformY)):
        raise ValueError(f"rfwaveformX ({len(rfwaveformX)}) and rfwaveformY ({len(rfwaveformY)}) are not of equal length")
    
    dTau = Tau / nSteps

    propagatorList = []

    for idx, rfVals in enumerate(zip(rfwaveformX,rfwaveformY)):
        BlochMat = np.array([[-R_2, offsetZ, -rfVals[1], 0],[-offsetZ, -R_2,rfVals[0],0],[rfVals[1], -rfVals[0], R_1, R_1],[0,0,0,0]])
        if idx==0:
            propagator = scipy.linalg.expm(BlochMat*dTau)
        else:
            propagator = np.dot(scipy.linalg.expm(BlochMat*dTau), propagator)

    return propagator