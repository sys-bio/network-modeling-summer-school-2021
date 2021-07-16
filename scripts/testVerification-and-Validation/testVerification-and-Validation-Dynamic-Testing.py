#!/usr/bin/env python
# coding: utf-8

# # VERIFICATION AND VALIDATION: DYNAMIC TESTING
# 

# # Preliminaries

# In[ ]:


IS_COLAB = True


# In[ ]:


if IS_COLAB:
    pass


# ## Imports

# In[ ]:


# Python packages used in this section
import tellurium as te
import matplotlib.pyplot as plt
import numpy as np
import urllib.request # use this library to download file from GitHub
from scipy import fftpack
from SBMLLint.tools.sbmllint import lint


# ## Helper Functions

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


# ## Constants

# In[ ]:


# Constants
ENDTIME = 10
NUMPOINTS = 100


# # Motivating Example

# Dynamic testing is about verifying the output produced by a simulation.
# 
# This is illustrated at the beginning of the chapter using the example
# of a linear pathway consisting of three species ``S1``, ``S2``, and ``S3``. That is, ``S1`` $\rightarrow$ ``S2`` $\rightarrow$ ``S3``, with ``S1`` initialized to 10. 
# 
# We showed that the following outcomes were expected from the simulation:
# 
# - **O1**: ``S1`` decreases monotonically from 10 to 0.
# - **O2**: ``S2`` increases initially, peaks, and then declises to 0.
# - **O3**: ``S3`` increases monotonically to 10. 
# 

# Let's consider an initial implementation that contains an error due to a typo in the kinetics law.

# In[ ]:


INCORRECT_MODEL = '''
R1: S1 -> S2; k1*S1  
R2: S2 -> S3; k2*S1

S1 = 10
k1 = 1; k2 = 1
'''
rr = te.loada(INCORRECT_MODEL)
data = rr.simulate(0, ENDTIME, NUMPOINTS)
rr.plot(data, title="First model of a linear pathway")


# We see that ``S1`` and ``S3`` generally follow the expected pattern. But ``S2`` does not.
# This is because of the error in the kinetics expression.
# We correct this below.

# In[ ]:


CORRECT_MODEL = '''
R1:  S1 -> S2; k1*S1  
R2a: S2 -> S3; k2*S2

S1 = 10
k1 = 1; k2 = 1;
'''
rr = te.loada(CORRECT_MODEL)
data = rr.simulate(1, 10)
rr.plot(data, title="Correct model of a linear pathway")


# # Testing Expected Outcomes

# The expected outcomes are the starting point for constructing verification tests.
# Implementing these tests requires two Boolean helper functions.
# A Boolean function returns either ``True`` or ``False``.
# The first helper function, ``isMonotone``, takes as input a list or array of values and
# the direction expected (+1 for monotone increasing and -1 for monotone decreasing).
# The function returns ``True`` if the data are monotone in the desired direction.

# In[ ]:


def isMonotone(values, direction=1):
    """
    Check if series is monotone in the desired direction.

    Parameters
    ----------
    values: np.array
    
       direction: 1 (increasing), -1 (decreasing)
    """
    arr = np.array(values)
    diff = direction*(arr[1:] - arr[:-1])
    return all([v > 0 for v in diff])


# ``isMonotone`` allows us to check the conditions O1 and O3. For O2, we need to check that the values initially increase,
# and then decrease. Curves with this shape are called *concave*. That is, there is some intermediate point at which the series changes from increasing to
# decreasing. 

# In[ ]:


def isConcave(values, direction=1):
    # Check if series is monotone in the desired direction
    #   direction: 1 (increasing), -1 (decreasing)
    arr = np.array(values)
    result = False
    for idx in range(1, len(values)-2):
        isIncreasing = isMonotone(arr[:idx], direction=1)
        isDecreasing = isMonotone(arr[idx:], direction=-1)
        if isIncreasing and isDecreasing:
            result = True
            break
    return result


# We illustrate these functions below.

# In[ ]:


# Illustration of isMonotone, isConcave
for values in [ [1, 2, 3, 4], [4, 3, 2, 1], [1, 4, 3, 2]]:
    print("\nAnalyzing the series: %s" % str(values))
    print("  Increasing? %s" % str(isMonotone(values, direction=1)))
    print("  Decreasing? %s" % str(isMonotone(values, direction=-1)))
    print("  Concave? %s" % str(isConcave(values)))


# # Test Infrastructure

# Now that we have a way to test the outcomes, we need to write codes that use these tests.
# Fortunately, python provides a convenient way to do this with its ``unittest`` infrastructure.
# Using the infrastructure requires the following:
# - Copy the test template. The template consists of:
#   - ``import unittest``
#   - ``class TestModel(unittest.TestCase)``
# - Write a setUp function within ``TestModel`` that loads the model and gets simulation data.
# - Write tests that check the appropriate values in the simulation output.
# 
# This is done below for the linear pathway.
# Note that the code expects that the variable ``MODEL`` is the antimony codes for the model being tested. We see that there
# is a separate test for each column.

# In[ ]:


import unittest

# Performs unittests on MODEL

# Define a class in which the tests will run
class TestModel(unittest.TestCase):
    
    def setUp(self):
        self.roadrunner = te.loada(MODEL)
        self.data = self.roadrunner.simulate(0, 10, 100)
        
    def testS1(self):
        self.assertTrue(isMonotone(self.data["[S1]"], direction=-1))
          
    def testS2(self):
        self.assertTrue(isConcave(self.data["[S2]"]))
                
    def testS3(self):
        self.assertTrue(isMonotone(self.data["[S3]"], direction=1))


# To run the model, you assign ``MODEL`` and then run the two lines that follow this assignment.

# In[ ]:


MODEL = CORRECT_MODEL
suite = unittest.TestLoader().loadTestsFromTestCase(TestModel)
_ = unittest.TextTestRunner().run(suite)


# These codes ran the 3 tests. All three passed for ``CORRECT_MODEL``. Now consider what happens with the model in ``INCORRECT_MODEL`` that fails on the test for ``S2``.

# In[ ]:


MODEL = INCORRECT_MODEL
suite = unittest.TestLoader().loadTestsFromTestCase(TestModel)
_ = unittest.TextTestRunner().run(suite)


# We have demonstrated these tests in the Jupyter notebook. However, in practice, the tests would be in
# a separate python modeule (a file ending in ``.py``). Ideally, these test files are run whenever there is a
# change to the model or its supporting code.
# Systems such as ``github`` provide a way to run these tests automatically whenever the model changes, a process referred to
# as *continuous integration*.

# # Exercise

# This exercise uses the Wolf model of glycolytic oscillations.
# 
# 1. What is the relationship between the concentrations of
# glucose, furctose_1_6_bisphosphate, and pyruvate after the initial transients.
# Does this relationship always hold?
# 
# 1. Propose a test that checks for that the order of the concentrations at steady state of glucose, furctose_1_6_bisphosphate, and pyruvate.
# Hint: Implement a function that checks that one vector
# is greater than a second vector after a specified index.
# 
# 1. Another characteristic of the Wolf model is the presence of oscillations at the same frequency.
# Describe in words how you would implement a test to
# checks for the presence of these oscillations.
