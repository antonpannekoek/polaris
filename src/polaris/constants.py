import numpy as np

# Required angles for demodulation process
ANGLES = {
    3: [0, 60, 120],
    4: [0, 45, 90, 135],
}


# Unnormalized demodulation matrix
# inner lists are the rows
MATRIX = {
    3: np.array([[1, 1, 1], [2, -1, -1], [0, -np.sqrt(3), np.sqrt(3)]]) / 3,
    4: np.array([[1, 1, 1, 1], [2, 0, -2, 0], [0, -2, 0, -2]]) / 4,
}
