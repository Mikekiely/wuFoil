from skopt import gp_minimize
from skopt.callbacks import DeltaYStopper
from skopt.plots import plot_convergence
from wuFoil.airfoil import cst_Airfoil, Airfoil, get_cst_variables
from wuFoil.analysis import SU2_Analysis
import matplotlib.pyplot as plt
import os

"""
Script to optimize airfoil using scikit-optimize's Bayesian optimization function
Uses cst parameterization and SU2 RANS simulations
runs at a set cl of 0.7
chord_length = 1
"""

# Input parameters here
altitude = 35000
mach = 0.78
cl = .7
n_cst_variables = 9
base_airfoil = 'rae2822.dat'    # Base airfoil, just used to set bounds on cst variables
failed_case_cd = 0.3            # Value to assign CD when the case is failed
n_calls = 20                 # Maximum amount of function calls for optimization

def obj_fun(x):
    # Inputs x: a values of cst
    # Returns drag coefficient
    a = x
    af = cst_Airfoil(a)
    af.generate_mesh(show_graphics=False, hide_output=True)
    af.set_flight_conditions(altitude, mach, cl=cl)
    analysis = SU2_Analysis(af, hide_output=True)
    analysis.n_processes = 8
    # analysis.solver = 'Euler'
    analysis.convergence = -5
    analysis.run_analysis()
    if analysis.cd:
        print(f'Results: cd = {analysis.cd}')
        return analysis.cd
    else:
        print('Case Failed')
        return failed_case_cd

# get cst variables of rae2822
af = Airfoil(base_airfoil)
cst_vars = get_cst_variables(af, n_vars=n_cst_variables)

# set bounds - Realistically you should use a better method for this
bounds = []
range = .15
for a in cst_vars:
    bounds.append([a-range, a+range])
bounds[-1] = [0, .025]

# Run optimization
print(bounds)
res = gp_minimize(obj_fun,  # the function to minimize
                  bounds,  # the bounds on each dimension of x
                  acq_func="EI",  # the acquisition function
                  n_calls=n_calls,  # the number of evaluations of f
                  n_initial_points=15,  # the number of random initialization points
                  callback=DeltaYStopper(.0001, n_best=40),
                  initial_point_generator="lhs")

# Post Processing
print(res)
plot_convergence(res)
plt.rcParams["figure.figsize"] = (8, 14)
plt.show()

af = cst_Airfoil(res.x)
af.plot_airfoil(show_plot=True)


