import numpy as np

def simon_HHH(params, m, M):
    # Unpack the params
    num = m.shape[2]
    R = []
    
    for i in range(num):
        R_new = params[(i * 6):(i * 6 + 6)]
        Q1 = R_new[0]
        Q2 = R_new[1]
        Q3 = R_new[2]
        TL = R_new[3:6]
        
        RL = np.array([
            [np.cos(Q2) * np.cos(Q1),  np.sin(Q2) * np.cos(Q1), -np.sin(Q1)],
            [-np.sin(Q2) * np.cos(Q3) + np.cos(Q2) * np.sin(Q1) * np.sin(Q3), np.cos(Q2) * np.cos(Q3) + np.sin(Q2) * np.sin(Q1) * np.sin(Q3), np.cos(Q1) * np.sin(Q3)],
            [np.sin(Q2) * np.sin(Q3) + np.cos(Q2) * np.sin(Q1) * np.cos(Q3), -np.cos(Q2) * np.sin(Q3) + np.sin(Q2) * np.sin(Q1) * np.cos(Q3), np.cos(Q1) * np.cos(Q3)]
        ])
        
        RT = np.hstack((RL[:, :2], TL.reshape(-1, 1)))
        R.append(RT)
    
    R = np.vstack(R)
    
    k1 = params[num * 6]
    k2 = params[num * 6 + 1]
    A = np.array([
        [params[num * 6 + 2], params[num * 6 + 3], params[num * 6 + 4]],
        [0, params[num * 6 + 5], params[num * 6 + 6]],
        [0, 0, 1]
    ])
    
    u0 = A[0, 2]
    v0 = A[1, 2]
    
    D = []
    d = []
    npts = m.shape[1]
    
    for flag in range(num):
        RT = R[(flag * 3):(flag * 3 + 3), :]
        XY = RT @ M
        UV = A @ XY
        UV = np.vstack((UV[0, :] / UV[2, :], UV[1, :] / UV[2, :], UV[2, :] / UV[2, :]))
        XY = np.vstack((XY[0, :] / XY[2, :], XY[1, :] / XY[2, :], XY[2, :] / XY[2, :]))
        
        radial_distortion = (XY[0, :] ** 2 + XY[1, :] ** 2)
        D.append(UV[0, :] + ((UV[0, :] - u0) * radial_distortion) * k1 + ((UV[0, :] - u0) * (radial_distortion ** 2)) * k2)
        D.append(UV[1, :] + ((UV[1, :] - v0) * radial_distortion) * k1 + ((UV[1, :] - v0) * (radial_distortion ** 2)) * k2)
        
        d.append(m[0, :, flag])
        d.append(m[1, :, flag])
    
    D = np.array(D)
    d = np.array(d)
    
    f = d - D
    return f.ravel()
