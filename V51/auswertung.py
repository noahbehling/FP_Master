import numpy as np
import matplotlib.pyplot as plt

for i in range(3):
    f, U_a, delta_t = np.genfromtxt(
        f"data/linearverstärker_{i + 1}.txt", delimiter=", ", unpack=True
    )
    print(f)
