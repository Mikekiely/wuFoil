import os
import subprocess
import logging
import pandas as pd
from io import StringIO
from wuFoil.airfoil import Airfoil

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class SU2_Analysis:
    """
    Runs Analysis in SU2 with given airfoil

    Parameters:
    -   airfoil: airfoil to be tested
    -   aoa: <float or list of floats> angle of attack (degrees)
    -   cl: <float or list of floats> lift coefficient
    -   convergence: <float> su2 convergence criteria, log10 of density residual
    -   hide_output: <bool> Used to mute the data to screen
    -   altitude: <float> altitude in ft
    -   mach: <float> mach number
    -   Processes: <int> used for multi processing, number of processors used
    -   max_iter: <int> maximum number of iterations used in su2 simulation
    -   Timeout: <int> maximum time used for SU2 simulation in seconds
    """
    SU2_parameters = []
    flight_conditions = []
    airfoil = []
    prefix = ''
    history_file = 'history.csv'


    # Results
    cd = []
    cl = []

    def __init__(self, airfoil: Airfoil, hide_output: bool = False):
        """
        Set su2 solution parameters
        Note: either cl or aoa must be input
        """
        if not airfoil.flight_conditions:
            logger.error('Please set airfoil flight conditions before analysis')

        self.n_processes = 1
        self.airfoil = airfoil
        self.prefix = airfoil.name
        self.flight_conditions = airfoil.flight_conditions
        self.altitude = self.flight_conditions.altitude
        self.mach = self.flight_conditions.mach
        self.hide_output = hide_output
        self.convergence = -6
        self.solver = 'RANS'  # Solution Method
        self.turb_model = 'SA'  # Turbulence Model
        self.max_iter = 1000  # Maximum Iterations
        self.mesh_name = self.prefix + '.su2'
        self.chord_length = airfoil.chord_length  # Chord length (ft)
        self.re = self.flight_conditions.re

        self.aoa = self.flight_conditions.aoa
        self.cl = self.flight_conditions.cl

        dir = os.path.dirname(os.path.abspath(__file__))
        self.base_cfg = os.path.join(dir, 'base.cfg')

    def run_analysis(self):
        """
        Run analysis with preset parameters
        """

        # Check if airfoil mesh exists, create it if it doesn't
        # if not os.path.exists(f'{self.prefix}.su2'):
        #     logger.error(f'Mesh file {self.prefix}.su2 not found in folder')
        #     return

        # print(os.path.exists(f'{self.prefix}.su2'))
        mesh_file_path = self.airfoil.mesh_dir / self.airfoil.name
        mesh_file_path = str(mesh_file_path) + '.su2'
        if not os.path.exists(mesh_file_path):
            print('Please generate mesh first')

        self.write_cfg_file()

        # Create commands to run SU2
        commands = ''
        file_name = f'{self.prefix}.cfg'
        # commands += f"cd {self.airfoil.base_dir}\n"
        if self.n_processes > 1:
            commands += f"mpiexec -n {self.n_processes} "
        commands += f"SU2_CFD {str(file_name)}"
        print(commands)
        try:
            if self.hide_output:
                proc = subprocess.Popen(['cmd', '/c', commands], shell=False, env=os.environ, stdin=subprocess.PIPE,
                                        stdout=subprocess.DEVNULL, start_new_session=True, cwd=self.airfoil.base_dir)
            else:
                proc = subprocess.Popen(['cmd', '/c', commands], shell=False, stdin=subprocess.PIPE, env=os.environ,
                                        start_new_session=True,  cwd=self.airfoil.base_dir)

            try:
                stdout, stderr = proc.communicate()
            except subprocess.TimeoutExpired:
                proc.kill()
                out, errs = proc.communicate()
        except FileNotFoundError:
            logger.error(f'SU2_RUN not found in environment variables')
            return
        except Exception as e:
            logger.warning(f'Unknown error occured: {e}')
            return

        self.read_output()

    def write_cfg_file(self):
        """
        Writes cfg file by altering specific lines in the base cfg file
        """
        dir = os.path.dirname(os.path.abspath(__file__))
        base_file = os.path.join(dir, 'base.cfg')
        path = self.airfoil.base_dir
        with open(base_file, 'r') as f:
            lines = f.readlines()
        newlines = []
        # TODO fix this elif chain
        file_name = f'{self.prefix}.cfg'
        file_path = path / file_name
        with open(file_path, 'w') as f:
            for line in lines:
                if line.startswith('SOLVER'):
                    line = f'SOLVER= {self.solver.upper()}\n'
                elif line.startswith('KIND_TURB_MODEL'):
                    line = f'KIND_TURB_MODEL= {self.turb_model}\n'
                elif line.startswith('MACH_NUMBER'):
                    line = f'MACH_NUMBER= {self.mach}\n'
                elif line.startswith('AOA'):
                    line = f'AOA= {self.flight_conditions.aoa}\n'
                elif line.startswith('FREESTREAM_PRESSURE'):
                    line = f'FREESTREAM_PRESSURE= {self.flight_conditions.p}\n'
                elif line.startswith('FREESTREAM_TEMPERATURE'):
                    line = f'FREESTREAM_TEMPERATURE= {self.flight_conditions.t}\n'
                elif line.startswith('REYNOLDS_NUMBER'):
                    line = f'REYNOLDS_NUMBER= {self.re}\n'
                elif line.startswith('REYNOLDS_LENGTH'):
                    line = f'REYNOLDS_LENGTH= {self.flight_conditions.length}\n'
                elif line.startswith('REF_ORIGIN_MOMENT_X'):
                    line = f'REF_ORIGIN_MOMENT_X= {.25*self.flight_conditions.length}\n'
                elif line.startswith('REF_LENGTH'):
                    line = f'REF_LENGTH= {self.flight_conditions.length}\n'
                elif line.startswith('REF_AREA'):
                    line = f'REF_AREA= {self.flight_conditions.length}\n'
                elif line.startswith('ITER') and not line.startswith('ITER_'):
                    line = f'ITER= {self.max_iter}\n'
                elif line.startswith('CONV_RESIDUAL_MINVAL'):
                    line = f'CONV_RESIDUAL_MINVAL= {self.convergence}\n'
                elif line.startswith('CONV_FILENAME'):
                    line = f'CONV_FILENAME= {self.prefix}_history\n'
                # elif line.startswith('MESH_FILENAME'):
                #     line = f'MESH_FILENAME= {self.prefix}.su2\n'
                elif line.startswith('FIXED_CL_MODE') and self.cl:
                    line = f'FIXED_CL_MODE= YES\n'
                elif line.startswith('TARGET_CL') and self.cl:
                    line = f'TARGET_CL= {self.flight_conditions.cl}\n'
                elif line.startswith('DV_VALUE'):
                    line = f'DV_VALUE= {self.flight_conditions.length}\n'
                f.writelines(line)

    def read_output(self):
        """
        Load history file and set data variables
        """
        # Load history file
        path = self.airfoil.base_dir
        file_name = f'{self.prefix}_history.csv'
        file_path = path / file_name
        df = pd.read_csv(file_path, sep=r'\s*,\s*', engine='python')
        # Make sure convergence was reached
        self.cd = []
        self.cl = []
        if df['"rms[Rho]"'].iloc[-1] < self.convergence:
            self.cd = df['"CD"'].iloc[-1]
            self.cl = df['"CL"'].iloc[-1]
            self.aoa = df['"AoA"'].iloc[-1]
        else:
            print('convergence failed')


class xfoil_analysis:
    """
    Runs analysis in xfoil with given parameters

    Parameters:
    -   airfoil: airfoil to be tested
    -   aoa: <float or list of floats> angle of attack (degrees)
    -   cl: <float or list of floats> lift coefficient
    -   hide_output: <bool> Used to mute the data to screen
    -   altitude: <float> altitude in ft
    -   mach: <float> mach number
    -   max_iter: <int> maximum number of iterations used in su2 simulation
    """
    flight_conditions = []
    aoa = []
    cl = []
    cd = []
    airfoil = []
    output_file = 'data.dat'
    x_tr_top = []       # transition points
    x_tr_bottom = []

    def __init__(self, airfoil, aoa=None, cl=None, interp_points=None,
                 hide_output=False, max_iter: int = 40):
        """
        Initializes variables for xfoil analysis
        Can either define analysis by a set angle of attack or by a set lift coefficient, cl
        Additionally, you can analyze several angles of attack by inputting them in as a list

        Input cl and aoa values override those in airfoil.flight_conditions

        Note: must either specify aoa or cl
        """
        if not airfoil.flight_conditions:
            logger.error('Please set airfoil flight conditions before analysis')
        self.prefix = airfoil.name
        self.hide_output = hide_output
        self.output_file = self.prefix + '.dat'
        self.airfoil = airfoil
        self.flight_conditions = airfoil.flight_conditions
        self.interp_points = interp_points
        if aoa:
            if type(aoa) is not list:
                aoa = [aoa]
            self.aoa = aoa
        elif cl:
            if type(cl) is float or type(cl) is int:
                cl = []
            self.cl = cl
        else:
            self.aoa = [airfoil.flight_conditions.aoa]
            self.cl = [airfoil.flight_conditions.cl]
        self.max_iter = max_iter
        self.re = self.flight_conditions.re

    def run_analysis(self):
        """
        Creates keystrokes for xfoil analysis, runs xfoil
        """
        # Delete data file if it exists
        if os.path.exists(self.output_file):
            os.remove(self.output_file)

        self.create_airfoil_file()
        # Write commands to be input into xfoil
        keystrokes = (f"load\n"
                      f"{self.prefix}_xfoil.dat\n"
                      f"{self.airfoil.name}\n")
        # interpolate airfoil if required
        if self.interp_points:
            keystrokes += f'pane {self.interp_points}\n'
        keystrokes += f"oper\n" \
                      f"visc {self.flight_conditions.re}\n" \
                      f"M {self.flight_conditions.mach}\n" \
                      f"PACC\n" \
                      f"{self.output_file}\n" \
                      f"\n"


        keystrokes += f"iter {self.max_iter}\n"
        # Add all aoa and cl cases
        if self.cl:
            for c in self.cl:
                keystrokes += f"cl {c}\n"
        else:
            for alpha in self.aoa:
                keystrokes += f"alfa {alpha}\n"

        keystrokes += f"\nquit\n"

        # Run xfoil
        try:
            common_args = {'stdin': subprocess.PIPE,
                           'text': True,
                            'creationflags': subprocess.CREATE_NO_WINDOW,
                           }

            if self.hide_output:
                common_args['stdout'] = subprocess.DEVNULL

            process = subprocess.call(['xfoil', '/c', keystrokes], **common_args)
        except FileNotFoundError:
            logger.error('Xfoil.exe file not found, please add to working directory or add xfoil to environment variables')
            return

        os.remove(f"{self.prefix}_xfoil.dat")

        self.read_output()

    def read_output(self):
        with open(self.output_file, 'r') as file:
            lines = file.readlines()

        start_index = 0
        for i, line in enumerate(lines):
            if line.replace(" ", "").startswith('---'):
                titles_index = i-1
                start_index = i + 1
                break

        lines_df = StringIO("".join(lines[start_index:]))
        column_names = lines[titles_index].split()
        # Extract data starting from the line where numeric data starts
        df = pd.read_csv(lines_df, delim_whitespace=True, names=column_names)
        self.aoa = []
        self.cl = []
        self.x_tr_top = []
        self.x_tr_bottom = []
        for i in range(len(df['alpha'])):
            self.aoa.append(df['alpha'].iloc[i])
            self.cl.append(df['CL'].iloc[i])
            self.cd.append(df['CD'].iloc[i])
            self.x_tr_top.append(df['Top_Xtr'].iloc[i])
            self.x_tr_bottom.append(df['Bot_Xtr'].iloc[i])

    def create_airfoil_file(self):
        """
        Writes airfoil coordinates into a formatted file readable by xfoil
        """
        with open(f'{self.prefix}_xfoil.dat', 'w') as file:
            for point in self.airfoil.points:
                file.writelines(str(point[0]) + ' ' + str(point[1]) + '\n')





