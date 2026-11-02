import numpy as np
import matplotlib.pyplot as plt

def getObs():
    #STANDARD RUN: A graph of the mld over the timeseries according to the observational output from the model 

    obs = np.loadtxt("observations.t")
    year = obs[:,0]
    mld = obs[:,3]

    plt.figure(figsize=(8, 6))
    plt.scatter(year, mld)
    plt.gca().invert_yaxis()

    plt.title('Observational MLD over Timeseries from Standard Run', fontsize=14, fontweight='bold')
    plt.xlabel('Year', fontsize=12)
    plt.ylabel('MLD (m)', fontsize=12)
    plt.show()

if __name__ == "__main__":
    getObs()