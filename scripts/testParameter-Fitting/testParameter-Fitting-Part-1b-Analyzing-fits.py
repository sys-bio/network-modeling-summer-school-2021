#!/usr/bin/env python
# coding: utf-8

# # PARAMETER FITTING PART 1B: ANALYZING FITS
# 

# This notebook provides explores the key elements of fitting parameters.

# # Preliminaries

# In[ ]:


IS_COLAB = False


# In[ ]:


if IS_COLAB:
    pass


# In[ ]:


# Python packages used in this chapter
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import urllib.request # use this library to download file from GitHub
import tellurium as te
from SBstoat.namedTimeseries import NamedTimeseries, TIME
from SBstoat.modelFitter import ModelFitter
import SBstoat


# In[ ]:


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


# # Linear Pathway Model & Data

# To make these ideas concrete, we use the linear pathway model.

# In[ ]:


print(LINEAR_PATHWAY_MODEL)


# Note that the values of the parameters are 0! This is not because their true value is 0.
# Rather, it's because we don't know the value of these parameters.

# It is very common to have parameters with unknown values in kinetics models. Typically, we know the range
# of the parameter values, but not a suitable value itself.
# 
# So, how do we run a model with known values for *all* of its parameters? The answer is that we estimate or **fit paramter values** based on observed data.
# 
# For this example, the observational data are in the file specified in ``LINEAR_PATHWAY_DATA``.

# In[ ]:


data = NamedTimeseries(dataframe=LINEAR_PATHWAY_DF)
fig, ax = plt.subplots(1)
for col in data.colnames:
    ax.scatter(data[TIME], data[col])
_ = plt.legend(data.colnames)


# The following sections provide details of parameter fitting using the ``PATHWAY_MODEL``.

# # Fitting Basics <a class="anchor" id="parameter-fitting-fitting-basics"></a>

# Some key concepts are:
# - A **parameter** is a constant whose value is unknown.
# - **Observed data** are data that are used for fitting parameter values.
# - **Fitted data** are data produced by a simulation for an assignment of values to parameters.
# - **Residuals** are the difference between observed data and simulated data.

# The basic algorithm for fitting is described below.
# The algorithm estimates ``parameters``, the vector of values of unknown parameters.
# The inputs to the algoritm are:
# - ``model`` - the simulation model
# - ``parameters`` parameters to estimate with their initial values
# - ``observedData`` observed data used for fitting
# 
# 1. **Initialization**
#    1. ``fittedData = simulate(model, parameters)``
#    1. ``residuals = observedData - fittedData``
# 1. **Loop**
#    1. ``parameters`` = adjust ``parameters`` based on ``residuals``
#    1. ``fittedData = simulation(parameters)``
#    1. ``residuals = observedData - fittedData``
#    1. if ``residuals`` are small
#        1. ``return parameters``

# The algorithm is fairly straight forward, with one exception.
# Step A - "adjust ``parameters`` based on ``residuals``" How are these adjustments made?
# 
# The adjustments are made based on an optimization that chooses parameter values that
# minimize the residuals.
# We give a brief discussion of optimization in the next section.

# # Exercise
# For the ``LINEAR_PATHWWAY_MODEL`` model and ``LINEAR_PATHWAY_DF`` data:
# 
# 1. Plot the residuals for all chemical species using a model with the true values of parameters. Use the values 1, 2, 3, 4 or ``k1``, ``k2``, ``k3``, ``k4``.
# Where are the residuals the largest? Why?
# 
# 1. Evaluate possible fits for ``k1``. You can quantify the "goodness" of a fit by calculating **RSQ**, the sum of squared residuals, for a set of parameter values. Run simulations with ``k1`` in the range $[0, 4]$ for 20 equally separated values. Plot RSQ. Where on the plot does the true value appear?
# 
# 1. Use the approach in (2) to propose an algorithm for parameter estimation. (No code is required; just a description in words.) What is the performance of this algorithm as the number of parameters increases?
