import numpy as np
from scipy.linalg import svd, eig
from scipy.optimize import least_squares
from homography2d1 import homography2d
from simon_HHH import simon_HHH

# Function for Zhang's Camera Calibration
# M: 2xN model plane points3 3
# m: 2xN corresponding image points
def Zhang(M, m):
    rows, npts = M.shape
    matrixone = np.ones((1, npts))  # 1 matrix
    M = np.vstack((M, matrixone))  # M = [X, Y, 1]^T
    num = m.shape[2]

    m_expanded = np.zeros((3, 256, 5))
    m_expanded[:2, :, :] = m  # 将原来的 m 填充到前两行

    # 将扩展后的矩阵的第三行赋值为 matrixone
    for i in range(m.shape[2]):
        m_expanded[2, :, i] = matrixone

    m = m_expanded

    # Estimate the homography matrix H
    H = np.zeros((3, 3, num))
    for i in range(num):
        H[:, :, i] = homography2d(M, m[:, :, i]).T  # Assuming homography2d1 is implemented

  
    # Solve the intrinsic parameters matrix A
    V = []
    for flag in range(num):
        v12 = np.array([
            H[0, 0, flag] * H[1, 0, flag],
            H[0, 0, flag] * H[1, 1, flag] + H[0, 1, flag] * H[1, 0, flag],
            H[0, 1, flag] * H[1, 1, flag],
            H[0, 2, flag] * H[1, 0, flag] + H[0, 0, flag] * H[1, 2, flag],
            H[0, 2, flag] * H[1, 1, flag] + H[0, 1, flag] * H[1, 2, flag],
            H[0, 2, flag] * H[1, 2, flag]
        ])
        v11 = np.array([
            H[0, 0, flag] ** 2,
            2 * H[0, 0, flag] * H[0, 1, flag],
            H[0, 1, flag] ** 2,
            2 * H[0, 2, flag] * H[0, 0, flag],
            2 * H[0, 2, flag] * H[0, 1, flag],
            H[0, 2, flag] ** 2
        ])
        v22 = np.array([
            H[1, 0, flag] ** 2,
            2 * H[1, 0, flag] * H[1, 1, flag],
            H[1, 1, flag] ** 2,
            2 * H[1, 2, flag] * H[1, 0, flag],
            2 * H[1, 2, flag] * H[1, 1, flag],
            H[1, 2, flag] ** 2
        ])
        V.append(v12)
        V.append(v11 - v22)

    V = np.array(V)
    k = V.T @ V
    u, d, v = svd(k)
    b = v[-1, :]  # Last row of V corresponds to smallest singular value

    # Compute intrinsic parameters
    v0 = (b[1] * b[3] - b[0] * b[4]) / (b[0] * b[2] - b[1] ** 2)
    s = b[5] - (b[3] ** 2 + v0 * (b[1] * b[3] - b[0] * b[4])) / b[0]
    alpha_u = np.sqrt(s / b[0])
    alpha_v = np.sqrt(s * b[0] / (b[0] * b[2] - b[1] ** 2))
    skewness = -b[1] * alpha_u ** 2 * alpha_v / s
    u0 = skewness * v0 / alpha_u - b[3] * alpha_u ** 2 / s

    A = np.array([
        [alpha_u, skewness, u0],
        [0, alpha_v, v0],
        [0, 0, 1]
    ])

    # Solve extrinsic parameters
    D = []
    d = []
    Rm = []
    for flag in range(num):
        s = (1 / np.linalg.norm(np.linalg.inv(A) @ H[0, :, flag]) +
             1 / np.linalg.norm(np.linalg.inv(A) @ H[1, :, flag])) / 2
        rl1 = s * np.linalg.inv(A) @ H[0, :, flag]
        rl2 = s * np.linalg.inv(A) @ H[1, :, flag]
        rl3 = np.cross(rl1, rl2)
        RL = np.column_stack((rl1, rl2, rl3))

        # Approximate RL to be a rotation matrix
        U, S, Vt = svd(RL)
        RL = U @ Vt

        TL = s * np.linalg.inv(A) @ H[2, :, flag]
        RT = np.column_stack((rl1, rl2, TL))
        XY = RT @ M
        UV = A @ XY
        UV = UV / UV[2, :]
        XY = XY / XY[2, :]

        for j in range(npts):
            D.append([
                (UV[0, j] - u0) * (XY[0, j] ** 2 + XY[1, j] ** 2),
                (UV[0, j] - u0) * (XY[0, j] ** 2 + XY[1, j] ** 2) ** 2
            ])
            D.append([
                (UV[1, j] - v0) * ((XY[0, j]) ** 2 + (XY[1, j]) ** 2),
                (UV[1, j] - v0) * ((XY[0, j]) ** 2 + (XY[1, j]) ** 2) ** 2
            ])
            d.append(m[0, j, flag] - UV[0, j])
            d.append(m[1, j, flag] - UV[1, j])
    
        r13 = RL[0, 2]
        r12 = RL[0, 1]
        r23 = RL[1, 2]

        Q1 = -np.arcsin(r13)
        Q2 = np.arcsin(r12 / np.cos(Q1))
        Q3 = np.arcsin(r23 / np.cos(Q1))

        result_matrix = np.array([
            [np.cos(Q2) * np.cos(Q1), np.sin(Q2) * np.cos(Q1), -np.sin(Q1)],
            [-np.sin(Q2) * np.cos(Q3) + np.cos(Q2) * np.sin(Q1) * np.sin(Q3), np.cos(Q2) * np.cos(Q3) + np.sin(Q2) * np.sin(Q1) * np.sin(Q3), np.cos(Q1) * np.sin(Q3)],
            [np.sin(Q2) * np.sin(Q3) + np.cos(Q2) * np.sin(Q1) * np.cos(Q3), -np.cos(Q2) * np.sin(Q3) + np.sin(Q2) * np.sin(Q1) * np.cos(Q3), np.cos(Q1) * np.cos(Q3)]
        ])

        R_new = [Q1, Q2, Q3, *TL]
        #Rm.append(R_new)
        Rm.extend(R_new)

    D = np.array(D)
    d = np.array(d)
    k = np.linalg.inv(D.T @ D) @ D.T @ d

    # Maximum Likelihood Estimation
    para = np.concatenate([Rm, [k[0], k[1], alpha_u, skewness, u0, alpha_v, v0]])
    result = least_squares(simon_HHH, para, args=(m, M), method='lm')

    x = result.x
    k1 = x[num * 6]
    k2 = x[num * 6 + 1]
    A = np.array([
        [x[num * 6 + 2], x[num * 6 + 3], x[num * 6 + 4]],
        [0, x[num * 6 + 5], x[num * 6 + 6]],
        [0, 0, 1]
    ])

    print(f"k1: {k1}, k2: {k2}")
    print(f"Intrinsic matrix A:\n{A}")


