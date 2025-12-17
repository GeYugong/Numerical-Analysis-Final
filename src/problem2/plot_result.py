import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_temperature_field():
    data = pd.read_csv("temperature_field.csv", header=None).values
    Ny, Nx = data.shape
    x = np.linspace(1.0/(2*Nx), 1.0-1.0/(2*Nx), Nx)
    y = np.linspace(1.0/(2*Ny), 1.0-1.0/(2*Ny), Ny)
    X, Y = np.meshgrid(x, y)

    plt.figure(figsize=(10, 8))
    cp = plt.contourf(X, Y, data, 50, cmap='inferno')
    plt.colorbar(cp).set_label('Temperature ($^\circ$C)')
    cs = plt.contour(X, Y, data, 10, colors='white', alpha=0.5)
    plt.clabel(cs, inline=True, fontsize=8)
    plt.title('FVM Steady State Temperature Distribution')
    plt.show()