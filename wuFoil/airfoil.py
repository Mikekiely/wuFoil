import logging
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splprep, splev
from scipy.optimize import minimize
from wuFoil.flight_conditions import FlightConditions
from wuFoil.meshing import mesh_parameters, generate_mesh

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class Airfoil:
    """
    A class containing all the methods and variables pertaining to regular airfoils
    To be used to generate and analyze airfoils from a dat file

    All coordinate points are input non-dimensionally
    Coordinates format: [[x0, y0], [x1, y1] ... [xn, yn]]
    Upper surface coordinates range from 1 to 0, lower surface coordinates range from 0 to 1
    Points variable combines these coordinates to make one closed loop from trailing edge to trailing edge

    All spline interpolation is handled in meshing/analysis, interpolation in this module is strictly for analysis

    Attributes:
    -   lower_surface: [[x0, y0], [x1, y1] ... [xn, yn]]
            X and Y coordinates of points on airfoil lower surface, formatted as numpy array
    -   upper_surface: [[x0, y0], [x1, y1] ... [xn, yn]]
            X and Y coordinates of points on airfoil upper surface, formatted as numpy array
    -   points: [[x0, y0], [x1, y1] ... [xn, yn]]
            x and y coordinates of all points on upper and lower surfaces in a closed loop
    -   chord_length: float
            Chord length, used for analysis
    -   max_camber: float
            value of maximum camber, nomalized by chord length
    -   x_mc: float
            x value at location of maximum camber
    -   max_thickness: float
            thickness normalized by chord length
    -   x_tc: float
            x value at location of maximum thickness
    -   mesh_parameters: dict
            Variables to be used in meshing the airfoil, all defined by chord lengths, see descriptions below
    """
    name = 'airfoil'

    upper_surface = []
    lower_surface = []
    points = []

    chord_length = 1

    max_camber = 0
    x_mc = 0
    max_thickness = 0
    x_tc = 0
    flight_conditions = []  # Only used for setting y+
    mesh_parameters = []

    def __init__(self, dat_file: str, chord_length=1, name='airfoil'):
        """
        Initialize airfoil coordinates from data file

        Parameters:
        -   input_params: data file containing airfoil coordinates. Accepted formats: Selig, Lednicer, should contain x
                          points for 0 and 1
        """
        self.name = name
        self.read_input(dat_file)
        self.mesh_parameters = mesh_parameters()
        self.chord_length = chord_length
        # self.set_tc_camber()

    def set_flight_conditions(self, altitude, mach, aoa = None, cl = None):
        self.flight_conditions = FlightConditions(altitude, mach, self.chord_length, aoa=aoa, cl=cl)

    def read_input(self, dat_file: str):
        """
        Read airfoil data file and generate lists of coordinates upper_surface and lower_surface

        Parameters:
        -   dat_file: <str> data file containing airfoil coordinates. Accepted formats: Selig, Lednicer
        """
        try:
            with open(dat_file, 'r') as f:
                # Read file and clean lines to readable format
                self.lower_surface = []
                self.upper_surface = []
                self.points = []
                lines = f.readlines()
                lower_surface_started = False
                airfoil_started = False
                selig = True
                lednicer = False
                for line in lines:
                    # Get data from lines, ignore titles and empty lines
                    try:
                        x, y = map(float, line.strip().split())
                    except:
                        continue

                    # skip titles and header, detect when the airfoil points start being listed
                    if not airfoil_started:
                        if x == 1:
                            selig = True
                            lednicer = False
                            airfoil_started = True
                        elif x == 0:
                            selig = False
                            lednicer = True
                            airfoil_started = True
                        else:
                            continue

                    # Handle selig airfoils
                    if selig:
                        if x == 0 and not lower_surface_started:
                            lower_surface_started = True
                            self.upper_surface.append([x, y])

                        if not lower_surface_started:
                            self.upper_surface.append([x, y])
                        else:
                            self.lower_surface.append([x, y])

                    if lednicer:
                        if x == 1.0 and not lower_surface_started:
                            lower_surface_started = True
                            self.upper_surface.append([x, y])
                            continue

                        if not lower_surface_started:
                            self.upper_surface.append([x, y])
                        else:
                            self.lower_surface.append([x, y])

                if lednicer:
                    self.upper_surface = self.upper_surface[::-1]
                self.points = self.upper_surface[:-1] + self.lower_surface
        except FileNotFoundError:
            logger.error(f"Airfoil File {dat_file} not found")
            return
        except Exception as e:
            logger.error(f"An unexpected error occured: {e}")
            return

    def plot_airfoil(self, show_plot: bool = True, line_type: str = '-'):
        """
        Create basic matplotlib plot of airfoil
        legend entry of airfoil defaults to the airfoil name

        Parameters:
            show_plot: <bool> Show plot or not. can be useful to set as False if you want to plot multiple airfoils
            line_type: <str> matplotlib argument of line type
        """
        x, y = zip(*self.points)
        plt.plot(x, y, line_type, label=self.name)
        plt.ylim([-.2, .2])
        if show_plot:
            plt.show()

    def set_tc_camber(self):
        """
        Sets value of maximum camber and thickness
        uses a spline interpolation of the top and bottom surfaces to ensure accuracy
        """
        n_points = 300  # number of points evaluated on spline

        x_u, y_u = zip(*self.upper_surface)
        x_l, y_l = zip(*self.lower_surface)
        # Calculate the spline representation of surfaces
        tck_u, _ = splprep([x_u, y_u], s=0)
        tck_l, _ = splprep([x_l, y_l], s=0)

        u_interp = np.linspace(0, 1, n_points)

        # Interpolate surface
        x_u, y_u = splev(u_interp, tck_u)
        x_l, y_l = splev(u_interp, tck_l)

        max_camber = 0
        max_thickness = 0
        x_mt = 0
        x_mc = 0
        x_u = x_u[::-1]
        y_u = y_u[::-1]
        for i in range(n_points):
            thickness = y_u[i] - y_l[i]
            camber = .5 * (y_u[i] + y_l[i])
            if thickness > max_thickness:
                max_thickness = thickness
                x_mt = x_u[i]
            if np.abs(camber) > max_camber:
                max_camber = camber
                x_mc = x_u[i]

        # plt.plot(x_l, y_l)
        # plt.plot(x_u, y_u)

        self.max_thickness = max_thickness
        self.x_tc = x_mt
        self.max_camber = max_camber
        self.x_mc = x_mc

    def interpolate_airfoil(self, n_points: int):
        """
        Uses cubic bspline to interpolate and refine airfoil

        Parameters:
        -   n_points: number of points used for interpolation on both upper and lower surfaces
                      will be lower than the resultant amount of points due to the leading edge refinement
        """
        x_interp_u, y_interp_u = zip(*self.upper_surface)
        x_interp_l, y_interp_l = zip(*self.lower_surface)
        tck_u, u_u = splprep([x_interp_u, y_interp_u], k=3, s=0)
        tck_l, u_l = splprep([x_interp_l, y_interp_l], k=3, s=0)
        u = np.linspace(0, 1, n_points)
        x_u, y_u = splev(u, tck_u)
        x_l, y_l = splev(u, tck_l)
        upper_surface = []
        lower_surface = []
        for x, y in zip(x_u, y_u):
            upper_surface.append([x, y])
        for x, y in zip(x_l, y_l):
            lower_surface.append([x, y])

        upper_surface[0][0] = 1
        lower_surface[0] = [0, 0]
        upper_surface[-1] = [0, 0]
        lower_surface[-1][0] = 1

        lower_surface[-1] = upper_surface[0]
        self.lower_surface = lower_surface
        self.upper_surface = upper_surface
        self.points = self.upper_surface[:-1] + self.lower_surface[::-1]

        self.refine_leading_edge(tck_u=tck_u, tck_l=tck_l, u=u)

    def refine_leading_edge(self, tol: float = 7, tck_u=None, tck_l=None, u=None, max_iter: int = 10):
        """
        Refines leading edge to ensure roundness using b-splines
        Uses recursive formula to ensure angle at leading edge is greater than 180 - tol
        Typical tolerance should be around 7 to 10

        Parameters:
        -   tol: Tolerance to determine convergence, converge when leading edge angle > 180 - tol (typically 8-10)
        -   tck_u, tck_l: Knot vectors from previously determined b splines
                          Useful when already interpolating airfoil, otherwise you will have to re evaluate the spline
                          Don't input these yourself unless you really know what you're doing
        -   max_iter: Maximum iterations of recursive loop
        """

        if tck_u is None:
            tck_u = []
        u = u.tolist()

        # Leading edge refinement
        upper_surface = self.upper_surface
        lower_surface = self.lower_surface
        for i in range(max_iter):
            # Extract the last two points from each list
            upper_point2, upper_point1 = np.array(upper_surface[-2:])
            lower_point1, lower_point2 = np.array(lower_surface[:2])

            # Calculate vectors representing the last two line segments
            vector_upper = upper_point2 - upper_point1
            vector_lower = lower_point2 - lower_point1

            # Calculate the cosine of the angle
            cosine_angle = np.dot(vector_upper, vector_lower) / (
                    np.linalg.norm(vector_upper) * np.linalg.norm(vector_lower))

            # Calculate the angle in radians
            le_angle = np.arccos(cosine_angle)

            # Convert the angle to degrees
            le_angle = np.degrees(le_angle)

            # check if convergence is reached, if not insert a new point halfway between last two points
            if np.abs(le_angle - 180) <= tol:
                break
            else:
                u_interp = .5 * (1 + u[-2])
                u.insert(-1, u_interp)

                y_u_interp = splev(u_interp, tck_u)
                y_l_interp = splev(1-u_interp, tck_l)
                upper_surface.insert(-1, [float(y_u_interp[0]), float(y_u_interp[1])])
                lower_surface.insert(1, [float(y_l_interp[0]), float(y_l_interp[1])])

        self.upper_surface = upper_surface
        self.lower_surface = lower_surface
        self.points = upper_surface[:-1] + lower_surface

    def set_desired_yplus(self, y_plus: float, altitude: float, mach: float):
        """
        Sets cell thickness at the first layer of the boundary layer to reach the desired y+ value
        Equations found here: https://www.cfd-online.com/Wiki/Y_plus_wall_distance_estimation

        Parameters:
        y_plus: <float> Desired y+ value at boundary layer, typically ~1.0
        chord_length: <float> Chord length of airfoil, used for reynolds number calculation
        altitude, mach: <float> flight conditions for determining atmospheric conditions

        Returns:
            first_cell_thickness: thickness of first cell
        """
        fc = FlightConditions(altitude, mach, self.chord_length)
        cf = (2 * np.log10(fc.re) - .65) ** -2.3    # skin friction coefficient
        tau_w = cf * .5 * fc.rho * fc.v ** 2        # Wall shear stress
        u_star = np.sqrt(tau_w / fc.rho)            # Friction velocity
        y = y_plus * fc.mu / (fc.rho * u_star)      # First cell thickness
        self.mesh_parameters.first_cell_thickness = y
        return y

    def generate_mesh(self, show_graphics: bool = True, output_format: str = '.su2',
                      hide_output: bool = False):
        # Just a way to generate the mesh through the airfoil class
        # Let me know if this is bad practice, it seems fine to me I don't want to write a multi hundred line code here
        generate_mesh(self, show_graphics=show_graphics, output_format=output_format,
                      hide_output=hide_output)


class cst_Airfoil(Airfoil):
    """
    A specific subclass of airfoil which uses CST parameterization to shape the airfoil instead of a data file
    Useful for optimization
    Still only works with sharp trailing edges
    For more info on CST Parameterization, see paper below

    Kulfan, B. M., “‘CST’ Universal Parametric Geometry Representation Method With Applications to Supersonic Aircraft”,
    Journal of Aircraft, Vol. 45, No. 1, Jan-Feb 2008, https://doi.org/10.2514/1.29958

    Parameters:
    -   a_l: cst variables used for lower surface
    -   a_u: cst variables used for upper surface
    -   yte: y value of airfoil at trailing edge
    """
    a_l = []
    a_u = []
    yte = 0

    def __init__(self, cst_variables: list[float], chord_length: float = 1, n_airfoil: int = 100,
                 leading_edge_refinement : bool = True, name='airfoil'):
        """
        Sets relevant airfoil parameters and calls cst shape function to set airfoil shap

        Parameters:
        -   cst_variables: <list of ints> Curvature coefficient weights for CST parameterization
                           Can use any amount of coefficients, 4-15 is the usual range
                           Even number of coefficients, y value at trailing edge = 0
                           Odd number of coefficients, y value at trailing edge = last coefficient
        -   n_airfoil: <int> number of points at which the airfoil will be evaluated on both top and bottom surface
                       Note - number of points on airfoil will be slightly higher than n_airfoil due to leading
                       edge refinement
                       Default = 100
        -   chord_length: <float> chord length of airfoil in feet, defualt = 1
        -   leading_edge_refinement: <bool> run leading edge refinement or not
        """
        # Set cst variables
        # Handle parameters with specified trailing edge height
        if len(cst_variables) % 2 == 1:
            self.yte = cst_variables[-1]
            a = cst_variables[:-1]
        else:
            a = cst_variables

        mid = len(a) // 2
        self.name = name
        self.chord_length = chord_length
        self.a_l = a[:mid]
        self.a_u = a[mid:]
        self.set_cst_shape(n_airfoil=n_airfoil, leading_edge_refinement=leading_edge_refinement)
        self.mesh_parameters = mesh_parameters()

    def set_cst_shape(self, n_airfoil: int = 100, leading_edge_refinement=True):
        """
        Sets shape of the airfoil according to the pre-defined cst parameters
        Refines leading edge by adding points until the resultant leading edge angle is close enough to 180

        Parameters:
        -   n_airfoil: <int> number of points interpolated at top or bottom of airfoil
        -   leading_edge_refinement: <bool> determines if the leading edge should be refined
        """

        # Set cst variables
        x = np.linspace(1, 0, n_airfoil)
        x_l = x
        y_l = []
        x_u = x
        y_u = []

        # evaluate cst equations
        for i in range(n_airfoil):
            y_u.append(self.evaluate_point(x_u[i], surface='upper'))
            y_l.append(self.evaluate_point(x_l[i], surface='lower'))

        upper_surface = []
        lower_surface = []
        for x, y in zip(x_u, y_u):
            upper_surface.append([x, y])
        for x, y in zip(x_l, y_l):
            lower_surface.append([x, y])

        # Leading edge refinement
        if leading_edge_refinement:
            max_iter = 10
            tol = 8    # degrees
            for i in range(max_iter):
                # Extract the last two points from each list
                upper_point1, upper_point2 = np.array(upper_surface[-2:])
                lower_point1, lower_point2 = np.array(lower_surface[-2:])

                # Calculate vectors representing the last two line segments
                vector_upper = upper_point2 - upper_point1
                vector_lower = lower_point2 - lower_point1

                # Calculate the cosine of the angle
                cosine_angle = np.dot(vector_upper, vector_lower) / (
                            np.linalg.norm(vector_upper) * np.linalg.norm(vector_lower))

                # Calculate the angle in radians
                le_angle = np.arccos(cosine_angle)

                # Convert the angle to degrees
                le_angle = np.degrees(le_angle)

                # check if convergence is reached, if not insert a new point halfway between last two points
                if np.abs(le_angle - 180) <= tol:
                    break
                else:
                    # Calculate new point
                    x = .5 * upper_point1[0]
                    # CST equations
                    y_u_interp = self.evaluate_point(x, surface='upper')
                    y_l_interp = self.evaluate_point(x, surface='lower')
                    upper_surface.insert(-1, [x, y_u_interp])
                    lower_surface.insert(-1, [x, y_l_interp])

        self.upper_surface = upper_surface
        self.lower_surface = lower_surface
        self.points = upper_surface[:-1] + lower_surface[::-1]

    def evaluate_point(self, x, surface='upper'):
        """
        Evaluate single point at defined x value

        Parameters:
        -   x: float, x value to be evaluated
        -   surface: 'upper' or 'lower', which surface of the airfoil to evaluate

        Returns:
        -   y: y value at given x
        """
        if surface == 'upper':
            a = self.a_u
        else:
            a = self.a_l

        n_vars = len(a)
        n = n_vars - 1

        s = 0
        for i in range(n_vars):
            s += a[i] * (np.math.factorial(n) / (np.math.factorial(i) * np.math.factorial(n - i))) * \
                  x ** i * (1 - x) ** (n - i)

        y = s * x ** .5 * (1 - x) ** 1 + x * self.yte / 2
        return y


def get_cst_variables(airfoil, tol: float = 1e-9, n_vars: int = 9,
                      a_guess = None, plot_airfoil: bool = False):
    """
    Uses SLSQP optimization to find cst variables for a given airfoil
    Optimizes based off mean least squares evaluated at all given points in dat file

    Parameters:
        airfoil: Object of airfoil class
        a_guess: first guess of cst variables, default for the rae 2822 airfoil
        tol: mean least squares tolerance for optimization function, default = 1e-9
        plot_airfoil: plot original points vs generated airfoil after

    returns:
        cst_variables: list
            optimum found cst variables
    """
    # Define objective function for optimization
    if not a_guess:
        a_guess = [0] * n_vars
    def objfun(a):
        a = a.tolist()
        af_cst = cst_Airfoil(a, leading_edge_refinement=False)
        error = 0
        n_points = len(airfoil.upper_surface) + len(airfoil.lower_surface)
        # Calculate mean least squares error over whole airfoil with given points
        for point in airfoil.upper_surface:
            x = point[0]
            y = point[1]
            y_cst = af_cst.evaluate_point(x, surface='upper')
            error += (y_cst - y) ** 2 / n_points

        for point in airfoil.lower_surface:
            x = point[0]
            y = point[1]
            y_cst = af_cst.evaluate_point(x, surface='lower')
            error += (y_cst - y) ** 2 / n_points
        return error
    # run optimization
    res = minimize(objfun, a_guess, method='SLSQP', tol=tol)
    cst_variables = res.x.tolist()

    # Plot airfoils
    if plot_airfoil:
        af_cst = cst_Airfoil(cst_variables)
        x, y = zip(*airfoil.points)
        xcst, ycst = zip(*af_cst.points)
        plt.plot(x, y, 'o')
        plt.plot(xcst, ycst)
        plt.ylim([-.2, .2])
        plt.show()

    return cst_variables
