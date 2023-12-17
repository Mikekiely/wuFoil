import multiprocessing as mp
import logging
from wuFoil.airfoil import cst_Airfoil, Airfoil
from wuFoil.analysis import SU2_Analysis, xfoil_analysis
import uuid
import os

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def _test_sample(af, method, params, i):
    """
    Tests single airfoil for batch analysis
    Parameters:
        af: fully initialized airfoil
    """
    prefix = str(uuid.uuid4())
    af.name = prefix

    # Initiate analysis
    if method.startswith('SU2'):
        af.generate_mesh(show_graphics=False, hide_output=True)
        analysis = SU2_Analysis(af, hide_output=True)
    else:
        analysis = xfoil_analysis(af, hide_output=True)

    for name, value in params.items():
        if hasattr(analysis, name):
            setattr(analysis, name, value)
        else:
            logger.warning(f'Variable {name} non in variables for {method} analysis')
    analysis.run_analysis()

    # Delete created files
    files_to_delete = [file for file in os.listdir() if file.startswith(prefix)]
    for file in files_to_delete:
        try:
            os.remove((file))
        except:
            pass

    if analysis.cd:
        # unpack if cd is returned as a part of a list (xfoil analysis only)
        return_vars = [analysis.cd, analysis.cl, analysis.aoa]
        for i, var in enumerate(return_vars):
            if isinstance(var, list):
                print(var)
                return_vars[i] = var[0]
        print(f'Iteration {i}: Cl = {return_vars[0]}, Cd = {return_vars[1]}, aoa = {return_vars[2]}')
        return return_vars
    else:
        print('Case Failed')
        return [None, None, None]

def analyze_batch(airfoils: list[Airfoil], n_processes: int = None, analysis_method: str = 'SU2',
                  analysis_parameters: dict = {}, output_file: str = None):
    """
    Analyzes a group of airfoils using multiprocessing
    aoa, and cl can either be entered as a constant value for all airfoils or as a list of values,
    each corresponding to one airfoil
    Either aoa or cl must be input
    make sure to use airfoil.set_flight_conditions before using this

    Parameters:
    -   airfoils: <list(airfoil objects)> list of airfoils to be analyzed
    -   output_file: <str> csv file  to output values to, won't output to a file if no file is inputted
                    Add headers to file
    -   n_processes: number of parallel processes to run, defaults to number of available processors
    -   altitude: Altitude at which the airfoil is operating. Either a constant value for all airfoils or a list of values
    -   analysis_method: Type of analysis run. 'SU2' or 'xfoil'
    -   analysis_parameters: List of variables to change in analysis

    Returns:
    -   cl: list of cl values
    -   cd: list of cd values
    """

    # Check to make sure aoa or cl were input
    for af in airfoils:
        if not af.flight_conditions:
            logger.error('Please set airfoil flight conditions with either aoa or cl before using parallel computation')
        if not af.flight_conditions.aoa and not airfoils.flight_conditions.cl:
            logger.error('Angle of attack or Lift coefficient not set with the airfoil flight condition, please enter one or the other')

    # make sure valid analysis method was input
    valid_analysis_methods = ['SU2', 'xfoil']
    if analysis_method not in valid_analysis_methods:
        logger.error(f'Analysis method {analysis_method} not in valid methods. Please choose one of the following: SU2_RANS, SU2_EULER, xfoil')

    # Make sure flight conditions have been entered
    n = len(airfoils)

    if not n_processes:
        n_processes = mp.cpu_count()
    pool = mp.Pool(processes=n_processes)
    # Start parallel processes
    jobs = []

    for i, af in enumerate(airfoils):
        job = pool.apply_async(_test_sample, (af, analysis_method, analysis_parameters, i))
        jobs.append(job)

    pool.close()
    # Wait for results from all processes
    results = [job.get() for job in jobs]
    pool.join()
    cd = [var[0] for var in results]
    cl = [var[1] for var in results]
    aoa = [var[2] for var in results]
    return cd, cl, aoa
