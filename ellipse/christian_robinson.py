import numpy as np

def christian_robinson(s):

    R = 1737.4
    D = np.array([[1/R, 0, 0],[0, 1/R, 0], [0, 0, 1/R]])
    n,m = s.shape
    # print(m,n)
    H = np.zeros((m, 3))
    for i in range(m):
        H[i, :] = (np.dot(D, s[:, i])).T
        H[i, :] = H[i, :] / np.linalg.norm(H[i, :])

    # Calculate Dinv
    Dinv = np.zeros((3, 3))
    Dinv[0, 0] = 1 / D[0, 0]
    Dinv[1, 1] = 1 / D[1, 1]
    Dinv[2, 2] = 1 / D[2, 2]

    # Perform TLS to find n
    # m, n = H.shape
    ones_column = np.ones((m, 1))
    augmented_H = np.hstack((H, ones_column))
    augmented_H = np.array(augmented_H)
    # print(augmented_H)
    _, S, V = np.linalg.svd(augmented_H, full_matrices=False)
    print(V)
    n = -V[3, 0:3] / V[3, 3]
    # print(V)
    # print(n)
    if np.dot(H[0, :3], n) > 0:
        n = -n

    # Calculate position in planetary frame
    r = -1 / np.sqrt(np.dot(n, n) - 1) * np.dot(Dinv, n)
    return r