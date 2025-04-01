import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = np.loadtxt("PHI_dataUDS.csv", delimiter=",")

x = np.linspace(-1.0, 1.0, data.shape[1])
y = np.linspace(0.0, 1.0, data.shape[0])

X, Y = np.meshgrid(x, y)

plt.figure(figsize=(15, 7.5))
plt.pcolormesh(X, Y, data, shading="gouraud", cmap="jet")
plt.colorbar(label="ϕ")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")

plt.savefig('PHI.png')

plt.show()

#Plotting outlet
referenceData=pd.read_csv("convection_diffusion_results.csv")
xData=referenceData["Position x"]
gamma10=referenceData["gamma=10"]

results = pd.read_csv("PHI_outletUDS.csv", header=None, names=["x", "PHI"])

plt.figure(figsize=(10, 5))
plt.scatter(xData, gamma10, color='r', marker='o', label="Reference")
plt.plot(results["x"], results["PHI"], linestyle='-', color='b', label="Results")
plt.xlim([-1.0, 1.0])
plt.ylim(0.0)
plt.xlabel("X (m)")
plt.ylabel("ϕ")
plt.legend()
plt.show()