import csv
import multiprocessing as mp
from wuFoil.airfoil import Airfoil, cst_Airfoil, get_cst_variables
from wuFoil.analysis import SU2_Analysis
from wuFoil.parallel_computation import analyze_batch
import numpy as np
import uuid

"""
Use this script to gather samples in parallel and output into csv file
This is useful for creating databases for machine learning algorithms

Randomizes variables based on normal distribution from base variables

DON'T USE THIS YET IT CRASHED MY COMPUTER IMMEDIATELY

"""


def main():
    # Initialize base airfoil
    # Likely you would want to use a list of predefined airfoils and randomly select one to randomize
    base_airfoil = 'rae2822.dat'

    base_altitude = 330000
    base_mach = .8
    base_cl = .65

    # get cst variables of rae2822
    af = Airfoil(base_airfoil)
    base_cst_vars = get_cst_variables(af, n_vars=8)

    batch_size = 20
    n_batches = 100
    header = ['a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'cl', 'ao', 're']
    with open('airfoil_results.csv', 'w', newline='') as file:
        # create csv file
        csv_writer = csv.writer(file)
        for _ in range(n_batches):

            # Create a batch of randomized airfoils
            airfoils = []
            for i in range(batch_size):
                # Create random set of cst variables and flight conditions
                a = []
                for a_base in base_cst_vars:
                    a.append(np.random.normal(a_base, .05, 1)[0])
                cl = np.random.normal(base_cl, .2, 1)[0]
                mach = np.random.normal(base_mach, .05, 1)[0]
                altitude = np.random.normal(base_altitude, 4000, 1)[0]

                airfoils.append(cst_Airfoil(a))
                airfoils[i].set_flight_conditions(altitude, mach, cl=cl)
            # run parallelized batch analysis
            cd, cl, aoa = analyze_batch(airfoils, analysis_method='xfoil')
            lines = [[a[i], [cd[i], cl[i], airfoils[i].flight_conditions.re]] for i in range(batch_size)]


if __name__ == '__main__':
    main()
