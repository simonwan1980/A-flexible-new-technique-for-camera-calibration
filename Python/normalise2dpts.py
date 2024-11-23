import numpy as np

def normalise2dpts(pts, c=None):
    if pts.shape[0] != 3:
        raise ValueError("pts must be 3xN")

    # Find the indices of the points that are not at infinity
    finiteind = np.where(np.abs(pts[2, :]) > np.finfo(float).eps)[0]

    if len(finiteind) != pts.shape[1]:
        print("Warning: Some points are at infinity")

    # For the finite points ensure homogeneous coords have scale of 1
    pts[0, finiteind] = pts[0, finiteind] / pts[2, finiteind]
    pts[1, finiteind] = pts[1, finiteind] / pts[2, finiteind]
    pts[2, finiteind] = 1

    # Calculate the centroid if not provided
    if c is None:
        c = np.mean(pts[:2, finiteind], axis=1)

    # Shift origin to centroid
    newp = np.zeros_like(pts)
    newp[0, finiteind] = pts[0, finiteind] - c[0]
    newp[1, finiteind] = pts[1, finiteind] - c[1]

    meandist = np.mean(np.sqrt(newp[0, finiteind] ** 2 + newp[1, finiteind] ** 2))

    scale = np.sqrt(2) / meandist

    T = np.array([
        [scale, 0, -scale * c[0]],
        [0, scale, -scale * c[1]],
        [0, 0, 1]
    ])

    newpts = T @ pts

    return newpts, T
