import numpy as np

def simon_H(H, m, M):
    # 解包参数
    H = np.reshape(H, (3, 3))  # 将 H 重新塑造成 3x3 矩阵
    h = H
    
    # 解包 params_const
    m = np.vstack((m[:2, :], np.ones((1, m.shape[1]))))  # 将 m 扩展为 3 行
    M = np.vstack((M[:2, :], np.ones((1, M.shape[1]))))  # 将 M 扩展为 3 行
    
    # 计算 X
    X = h @ M  # 矩阵相乘
    X = np.vstack((X[0, :] / X[2, :], X[1, :] / X[2, :], X[2, :] / X[2, :]))  # 归一化 X，使第三行全为 1
    
    # 计算误差
    res = m - X
    req = np.concatenate((res[0, :], res[1, :]))  # 将前两行的误差展开成一个向量

    return req
