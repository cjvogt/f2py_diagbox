import numpy as np
import matplotlib.pyplot as plt

def getObs():
    #STANDARD RUN: A graph of the mld over the timeseries according to the observational output from the model 

    obs = np.loadtxt("observations.t")
    year = obs[:,0]
    mld = obs[:,3]

    plt.figure(figsize=(8, 6))
    plt.scatter(year, mld)
    plt.show()

if __name__ == "__main__":
    getObs()