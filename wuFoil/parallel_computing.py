import multiprocessing as mp
import logging
from wuFoil.airfoil import cst_Airfoil, Airfoil
from wuFoil.analysis import SU2_Analysis, xfoil_analysis
import uuid
import os
import csv

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def _test_sample(af, method, params, i, q=None, output_file=None, output_parameters=None):
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

    # Put output into queue to be printed to csv
    if output_file:
        output = []
        for param in output_parameters:
            if hasattr(af, param):
                output_value = getattr(af, param)
                if isinstance(output_value, list):
                    output.extend(output_value)
                else:
                    output.append(output_value)

        q.put(output)


    # Return cd, cl, and aoa
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

def _print_sample(q, output_file):
    # Listens for data printed to the queue and prints it to file
    with open(output_file, 'a', newline='') as f:
        writer = csv.writer(f)
        while 1:
            m = q.get()
            if m == 'kill':
                break
            writer.writerow(m)
            f.flush()

def analyze_batch(airfoils: list[Airfoil], n_processes: int = None, analysis_method: str = 'SU2',
                  analysis_parameters: dict = {}, output_file: str = None,
                  output_parameters: list = ['cst_variables', 'cd', 'cl', 'aoa', 'mach', 're']):
    """
    Analyzes a group of airfoils using multiprocessing
    aoa, and cl can either be entered as a constant value for all airfoils or as a list of values,
    each corresponding to one airfoil
    Either aoa or cl must be input
    make sure to use airfoil.set_flight_conditions before using this

    Parameters:
    -   airfoils: <list(airfoil objects)> list of airfoils to be analyzed. If they are cst airfoils make sure they have the same number of variables
    -   output_file: <str> csv file  to output values to, won't output to a file if no file is inputted
                    Add headers to file
    -   output_parameters: <list> list of variables to output to csv file, name must be the same as they appear in the
                            analysis class
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

    n = len(airfoils)

    # Set number of processes
    if not n_processes:
        n_processes = mp.cpu_count()
    pool = mp.Pool(processes=n_processes)
    # Start parallel processes
    jobs = []

    # handle output files

    # Create file if it doesn't exist
    if os.path.exists(output_file):
        # Output file exists, make sure its a csv
        if not output_file.lower().endswith('.csv'):
            logger.error('Output file must be formatted as a csv')

        # Check if cst variables is in output_parameters. If it is, add a0 through an to headers
        if 'cst_variables' in output_parameters:
            cst_index = output_parameters.index('cst_variables')
            output_parameters.pop(cst_index)
            output_parameters.extend([f'a{i}' for i in range(len(airfoils[0].cst_variables))])

        # Check headers of output file
        with open(output_file, 'r', newline='') as file:
            reader = csv.reader(file)
            headers = next(reader, [])

            # Check if headers in prexisting file match output parameters


    if output_file:
        manager = mp.Manager()
        q = manager.Queue()
        watcher = pool.apply_async(_print_sample, (q, output_file))
    else:
        q = None

    # Start parallel processes
    for i, af in enumerate(airfoils):
        job = pool.apply_async(_test_sample, (af, analysis_method, analysis_parameters, i, q, output_file, output_parameters))
        jobs.append(job)

    # Signal listener process to exit
    if output_file:
        q.put('kill')

    pool.close()
    # Wait for results from all processes
    results = [job.get() for job in jobs]
    pool.join()
    cd = [var[0] for var in results]
    cl = [var[1] for var in results]
    aoa = [var[2] for var in results]
    return cd, cl, aoa
