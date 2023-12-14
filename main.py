from wuFoil.parallel_computation import analyze_batch
from wuFoil.airfoil import Airfoil
from wuFoil.analysis import SU2_Analysis, xfoil_analysis
import numpy as np


def xfoil_Analysis(af, hide_output):
    pass


if __name__=='__main__':
    mach = .78
    cls = np.linspace(.5, .9, 10)
    afs = []
    for i, cl in enumerate(cls):
        afs.append(Airfoil('rae2822.dat'))
        # af = Airfoil('rae2822.dat')
        afs[i].set_flight_conditions(35000, mach, cl=cl)
        afs[i].interpolate_airfoil(100)

    # analyze_batch(afs, analysis_parameters={'solver': 'Euler', 'convergence': -5, 'hide_output': True})
    analyze_batch(afs, analysis_method='xfoil', analysis_parameters={'hide_output': True})
    # mach = .577777777777
    # af = Airfoil('rae2822.dat')
    # af.set_flight_conditions(35000, mach, cl=.7)
    # af.interpolate_airfoil(100)
    # af.generate_mesh(show_graphics=False, hide_output=True)
    # analysis = xfoil_analysis(af, hide_output=True)
    # analysis.hide_output = False
    # analysis.run_analysis()
    print('dont')

