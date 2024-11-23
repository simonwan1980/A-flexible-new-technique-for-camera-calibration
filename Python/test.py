import numpy as np
import os
import logging
from Zhang import Zhang


# 设置日志记录器
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s')
logger = logging.getLogger(__name__)

# 读入数据并整理成可处理的形式
# M为模板数据
# m为提取出的数据

def load_data(filename):
    filepath = os.path.join(os.path.dirname(__file__), filename)
    if os.path.exists(filepath):
        return np.loadtxt(filepath)
    else:
        logger.warning(f"{filename} not found.")
        return np.zeros((8, 64))  # 使用默认值来避免文件未找到的错误

# 读取数据
M = load_data('Model.txt')  # M, m1, m2, m3, m4, m5初始得到都是8x64的一个矩阵
m1 = load_data('data1.txt')
m2 = load_data('data2.txt')
m3 = load_data('data3.txt')
m4 = load_data('data4.txt')
m5 = load_data('data5.txt')

# 重新整理数据
M = np.vstack([M[:, 0:2], M[:, 2:4], M[:, 4:6], M[:, 6:8]])
m1 = np.vstack([m1[:, 0:2], m1[:, 2:4], m1[:, 4:6], m1[:, 6:8]])
m2 = np.vstack([m2[:, 0:2], m2[:, 2:4], m2[:, 4:6], m2[:, 6:8]])
m3 = np.vstack([m3[:, 0:2], m3[:, 2:4], m3[:, 4:6], m3[:, 6:8]])
m4 = np.vstack([m4[:, 0:2], m4[:, 2:4], m4[:, 4:6], m4[:, 6:8]])
m5 = np.vstack([m5[:, 0:2], m5[:, 2:4], m5[:, 4:6], m5[:, 6:8]])

# 转置矩阵
M = M.T
m1 = m1.T
m2 = m2.T
m3 = m3.T
m4 = m4.T
m5 = m5.T

# 合并数据为三维数组
m = np.stack([m1, m2, m3, m4, m5], axis=2)

Zhang(M, m)

