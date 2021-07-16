IS_COLAB = False

if IS_COLAB:
    pass

# Python packages used in this chapter
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import pandas as pd
import urllib.request # use this library to download file from GitHub
from SBstoat.namedTimeseries import NamedTimeseries, TIME
from SBstoat.modelFitter import ModelFitter
import SBstoat
import tellurium as te

def getSharedCodes(moduleName):
  """
  Obtains common codes from the github repository.

  Parameters
  ----------
  moduleName: str
      name of the python module in the src directory
  """
  if IS_COLAB:
      url = "https://github.com/sys-bio/network-modeling-summer-school-2021/raw/main/src/%s.py" % moduleName
      local_python = "python.py"
      _, _ = urllib.request.urlretrieve(url=url, filename=local_python)
  else:
      local_python = "../../src/%s.py" % moduleName
  with open(local_python, "r") as fd:
    codeStr = "".join(fd.readlines())
  print(codeStr)
  exec(codeStr, globals())

# Acquire codes
getSharedCodes("util")

# TESTS
assert(isinstance(LINEAR_PATHWAY_DF, pd.DataFrame))

LINEAR_PATHWAY_ARR = LINEAR_PATHWAY_DF.to_numpy()

print(WOLF_MODEL)

print(LINEAR_PATHWAY_MODEL)

# True values of the parameters
PARAMETER_DCT = {"k1": 1, "k2": 2, "k3": 3, "k4":4}
# Create the Roadrunner object
rr = te.loada(LINEAR_PATHWAY_MODEL)
# Assign values to the kinetic constants
for name in PARAMETER_DCT.keys():
    rr[name] = PARAMETER_DCT[name]
# Run the simulation
modelArr = rr.simulate(1, 10, 100)
# Look at the result to ensure it makes sense
modelArr[0:10, :]

# Look at the observed data
LINEAR_PATHWAY_ARR[0:10, :]

# Calculate residuals
residualsArr = LINEAR_PATHWAY_ARR - modelArr
# Recover the time
residualsArr[:, 0] = modelArr[:, 0]
# Look at the data
residualsArr[0:10, :]

# Let's put this together into a single function that we can use later
def calcResiduals(modelStr, observedArr, parameterDct, startTime=1, endTime=10):
    """
    Calculates the residuals for parameter assignments for a model given the observed values.

    Parameters
    ----------
    modelStr: str
        Antimony model
    observedArr: np.array
        matrix of observed values; first column is time
    parameterDct: dict
        key: parameter name
        value: parameter value
    startTime: int
    endTime: int

    Returns
    -------
    np.array
        col 0: time
        col 1-: float
    """
    # Create the Roadrunner object
    rr = te.loada(modelStr)
    # Assign values to the kinetic constants
    for name in parameterDct.keys():
        rr[name] = parameterDct[name]
    # Run the simulation
    modelArr = rr.simulate(startTime, endTime, 10*endTime)
    # Calculate residuals
    residualsArr = observedArr - modelArr
    residualsArr[:, 0] = modelArr[:, 0] # recover time
    #
    return residualsArr

# TESTS
residualsArr = calcResiduals(LINEAR_PATHWAY_MODEL, LINEAR_PATHWAY_ARR, 
                             PARAMETER_DCT)
assert(len(residualsArr) > 0)
assert("NamedArray" in str(type(residualsArr)))

# Plot the result with separate plots for each species
fig = plt.figure(figsize=(12, 8))
numCol = 3
numRow = 2
# Use a customized plotting grid to avoid getting an empty plot
spec = gridspec.GridSpec(ncols=numCol, nrows=numRow, figure=fig)
# Do the plot
times = residualsArr[:, 0]
for idx, name in enumerate(modelArr.colnames[1:]):
    if idx < numCol:
        irow = 0
        icol = idx
    else:
        irow = 1
        icol = idx - numCol
    ax = fig.add_subplot(spec[irow, icol])
    ax.scatter(times, residualsArr[:, idx+1])
    title = name[1:-1] # Delete the square brackets around the species name
    ax.set_title(title)
    ax.set_ylim([-10, 10])

kVals = [n*0.2 for n in range(21)]
parameterName = "k1"
rsqDct = {} # Container for the simulation results
for k1Val in kVals:
    parameterDct = dict(PARAMETER_DCT)
    parameterDct[parameterName] = k1Val
    residualsArr = calcResiduals(LINEAR_PATHWAY_MODEL, LINEAR_PATHWAY_ARR,
                                 parameterDct)
    rsqDct[k1Val] = sum(sum(residualsArr**2))

# Look at the result
rsqDct

# Construct a plot
plt.plot(list(rsqDct.keys()), list(rsqDct.values()))
plt.xlabel("k1")
plt.ylabel("RSQ")
plt.title("Evaluation of k1 fits")

def functionToMinimize(parameterValues):
    """
    Parameters
    ----------
    parameterValues: list-float
        k1, k2, k3, k4
    """
    parameterDct = {"k%d" % (n + 1): parameterValues[n] for n in range(len(parameterValues))}
    residualsArr = calcResiduals(LINEAR_PATHWAY_MODEL, LINEAR_PATHWAY_ARR, parameterDct)
    return np.sum(residualsArr**2)

# TESTS
ssr1 = functionToMinimize([1, 2, 3, 4])
ssr2 = functionToMinimize([1, 1,  1, 1])
assert(ssr2 > ssr1)

from scipy.optimize import LinearConstraint

eye2Mat = np.identity(4) # Two dimensional identity matrix

constraints = LinearConstraint(eye2Mat, [0, 0, 0 , 0], [10, 10, 10, 10]),

# Optimization using gradient descent
from scipy.optimize import minimize

res = minimize(
    functionToMinimize,
    constraints=constraints,
    x0=(0, 0, 0, 0),
)

res.x

# Optimization using differential evolution
from scipy.optimize import differential_evolution

bounds = [(0, 10) for _ in range(4)]

res = differential_evolution(
    functionToMinimize,
    bounds=bounds,
)

res.x

print(WOLF_MODEL)

# Parameter ranges and initial value
upper = 1e5
WOLF_PARAMETERS = [
      SBstoat.Parameter("J1_k1", lower=0, value=1, upper=upper), #550
      SBstoat.Parameter("J1_Ki", lower=0, value=1, upper=upper), #1
      SBstoat.Parameter("J1_n", lower=0, value=1, upper=upper), #4
      SBstoat.Parameter("J2_k", lower=0, value=1, upper=upper), #9.8
      SBstoat.Parameter("J3_k", lower=0, value=1, upper=upper), #85.7
      SBstoat.Parameter("J4_kg", lower=0, value=1, upper=upper), #323.8
      SBstoat.Parameter("J4_kp", lower=0, value=1, upper=upper), #76411.1
      SBstoat.Parameter("J4_ka", lower=0, value=1, upper=upper), #57823.1
      SBstoat.Parameter("J4_kk", lower=0, value=1, upper=upper), #23.7
      SBstoat.Parameter("J5_k", lower=0, value=1, upper=upper), #80
      SBstoat.Parameter("J6_k", lower=0, value=1, upper=upper), #9.7
      SBstoat.Parameter("J7_k", lower=0, value=1, upper=upper), #2000
      SBstoat.Parameter("J8_k1", lower=0, value=1, upper=upper), #375
      SBstoat.Parameter("J8_k2", lower=0, value=1, upper=upper), #375
      SBstoat.Parameter("J9_k", lower=0, value=1, upper=upper), #28
      SBstoat.Parameter("J10_k", lower=0, value=1, upper=upper), #80                                                                                                                                                                          
]

def doFit(endTime=5, fitterMethods=["leastsq"], parametersToFit=WOLF_PARAMETERS,
          isTest=False):
    """
    Fit the Wolf Model.

    Parameters
    ----------
    endTime: int
    fitterMethods: list-str
    parametersToFit: list-SBstoat.Parameter

    Returns
    -------
    ModelFitter
    """
    model = te.loada(WOLF_MODEL)
    observedTS = NamedTimeseries(dataframe=WOLF_DF)
    # Find the last index to use
    lastIdx = len([t for t in observedTS[TIME] if t <= endTime])
    observedTS = observedTS[:lastIdx]
    # Construct the fitter and do the fit
    fitter = ModelFitter(model, observedTS, fitterMethods=fitterMethods, 
                         parametersToFit=parametersToFit)
    fitter.fitModel()
    if not isTest:
        print(fitter.reportFit())
        fitter.plotFitAll(figsize=(24,16), numCol=2, color=["red", "blue"])
    return fitter

# TESTS
result = doFit(isTest=True)
assert(isinstance(result, ModelFitter))

# Least squares fit
fitter = doFit(fitterMethods=["leastsq"])

# Differential evolution
doFit(fitterMethods=["differential_evolution"])

parametersToFit = [
      SBstoat.Parameter("J1_k1", lower=0, value=1, upper=1100), #550
      SBstoat.Parameter("J1_Ki", lower=0, value=1, upper=2), #1
      SBstoat.Parameter("J1_n", lower=0, value=1, upper=8), #4
      SBstoat.Parameter("J2_k", lower=0, value=1, upper=20), #9.8
      SBstoat.Parameter("J3_k", lower=0, value=1, upper=160), #85.7
      SBstoat.Parameter("J4_kg", lower=0, value=1, upper=660), #323.8
      SBstoat.Parameter("J4_kp", lower=0, value=1, upper=150000), #76411.1
      SBstoat.Parameter("J4_ka", lower=0, value=1, upper=110000), #57823.1
      SBstoat.Parameter("J4_kk", lower=0, value=1, upper=46), #23.7
      SBstoat.Parameter("J5_k", lower=0, value=1, upper=160), #80
      SBstoat.Parameter("J6_k", lower=0, value=1, upper=20), #9.7
      SBstoat.Parameter("J7_k", lower=0, value=1, upper=4000), #2000
      SBstoat.Parameter("J8_k1", lower=0, value=1, upper=750), #375
      SBstoat.Parameter("J8_k2", lower=0, value=1, upper=750), #375
      SBstoat.Parameter("J9_k", lower=0, value=1, upper=56), #28
      SBstoat.Parameter("J10_k", lower=0, value=1, upper=160), #80                                                                                                                                                                          
]

doFit(fitterMethods=["leastsqs"], parametersToFit=parametersToFit)

doFit(fitterMethods=["differential_evolution"], parametersToFit=parametersToFit)
