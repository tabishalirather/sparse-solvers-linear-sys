import matplotlib.pyplot as plt
import numpy as np
print("lock and load")


def tridiagonal_eigens(n):
    #     n = n
    #     i = np.arange(0, n)
    #     print("i: ", i)
    eigenvalues = np.zeros(n)
    for i in range(0, n):
        eigenvalues[i] = 4 * np.sin((i * np.pi) / (2*n))**2

    eigenvectors = np.zeros((n, n))
    for i in range(0, n):
        for j in range(0, n):
            eigenvectors[i, j] = np.sin((i * j * np.pi) / n)
    return eigenvalues, eigenvectors


d = 1
n = 11

eigenvalues, eigenvectors = tridiagonal_eigens(n)
print("eigenvalues: ", eigenvalues)
print("eigenvectors: ", eigenvectors)
index = np.arange(0, n)
print("index: ", index)

fig, axes = plt.subplots(1, 4, figsize=(16, 4))

data = [(1, 'o-', 'Eigenvector 1-smoothest'),
        (2, 's-', 'Eigenvector 2'),
        (5, 'd-', 'Eigenvector 5'),
        (10, 'x-', 'Eigenvector 10-loses smoothness')]

for ax, (i, style, label) in zip(axes, data):
    ax.plot(index, eigenvectors[i], style, label=label)
    ax.set_title(label)
    ax.set_xlabel('Index')
    ax.set_ylabel('Eigenvector Value')
    ax.grid()

plt.tight_layout()
plt.show()


# Cannot use this as this return eigens in random order, not partcular order
# tridiangular_matrix = np.diag(
#     2*np.ones(n-1), k=0) + np.diag(-1*np.ones(n-2), k=-1) + np.diag(-1*np.ones(n-2), k=1)
# print("tridiagonal matrix: \n", tridiangular_matrix)
