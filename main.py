from wuFoil.optimization import random_batch_generator, pso_optimize
import wuFoil as wf
import numpy as np
from Scripts.PyAero_meshing import generate_PyAero_mesh
import multiprocessing as mp

if __name__ == '__main__':
    a =  {"a0": -0.19991011591078262, "a1": -0.08092302942484535, "a2": -0.3570295711640754, "a3": -0.09156345461325682,
          "a4": -0.11997658917931232, "a5": 0.10022769017451416, "a6": 0.235351398133694, "a7": 0.07046268596552027,
          "a8": 0.23262680749657774, "a9": 0.38934560619821973}

    a = {"a0": -0.19982976073067785, "a1": -0.08100217321129922, "a2": -0.3567980974235295, "a3": -0.09174725116055779,
         "a4": -0.1201764852185435, "a5": 0.10019887660886045, "a6": 0.235286545077955, "a7": 0.07026223383191912,
         "a8": 0.23240049154433798, "a9": 0.38950656071629347}
    ai = np.fromiter(a.values(), dtype=float)
    af = wf.cst_Airfoil(ai, chord_length=20.5, name='airfoil_2')
    af.set_flight_conditions(35000, .78, aoa=1)
    # af.plot_airfoil()
    generate_PyAero_mesh(af)
    # af.generate_mesh()
    analysis = wf.SU2_Analysis(af)
    analysis.n_processes = 4
    analysis.max_iter = 5000
    analysis.convergence = -4
    analysis.run_analysis()

