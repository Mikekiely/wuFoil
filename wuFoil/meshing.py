import io
import os, sys
from contextlib import contextmanager
import gmsh
from scipy.interpolate import splprep, splev
import numpy as np


def generate_mesh(airfoil, show_graphics: bool = True, output_format: str = '.su2',
                  hide_output: bool = True):
    """
    Creates C-block mesh around airfoil using gmsh
    stores result in file airfoil.su2

    All units normalized by chord length

    Uses the dictionary mesh_parameters to model mesh
    Interpolates airfoil using cubic splines across whole closed loop for accuracy

    Currently only works with SU2 output format, I may add CGNS in the future
    Use meshio or something else to convert the mesh to another format if you want for now

    Naming conventions (can get very confusing, this is as much for me as it is for anyone else)
    -   Airfoil segments: Upper, Leading Edge, and Lower + _line
    -   Farfield Boundaries: Inlet, Top, Bottom, Top_wake, Bottom_wake, Outlet_top, Outlet_bottom + _line
    -   Points: af_te, af_top, af_bottom, inlet_top, inlet_bottom, bottom_te, top_te, wake_top, wake_bottom, wake_te
    -   Connecting lines: Named after points they connect, eg. afUpper_inletTop, afTe_bottomTe
    -   Sections: Inlet_section, Top_section, Bottom_section

    Parameters:
    -   Airfoil: <airfoil> airfoil or cst_airfoil object
    -   show_graphics: <bool> yes or no show gmsh gui
    -   output_format: <str> format of mesh to output (not currently used)
    """
    length_le = airfoil.mesh_parameters.leading_edge_length * airfoil.chord_length
    te_thickness = airfoil.mesh_parameters.trailing_edge_thickness  * airfoil.chord_length
    r_inlet = airfoil.mesh_parameters.inlet_radius * airfoil.chord_length
    downstream_distance = airfoil.mesh_parameters.downstream_distance * airfoil.chord_length
    first_layer_thickness = airfoil.mesh_parameters.first_cell_thickness * airfoil.chord_length
    n_airfoil = airfoil.mesh_parameters.n_airfoil
    n_volume = airfoil.mesh_parameters.n_volume
    n_wake = airfoil.mesh_parameters.n_wake
    n_leading_edge = airfoil.mesh_parameters.n_leading_edge
    center = (length_le, 0, 0)
    te = airfoil.upper_surface[0]
    points = [[x * airfoil.chord_length for x in coords] for coords in airfoil.points]

    # Split airfoil into 3 segments - upper, lower, and leading_edge. Create list of points in each
    upper_aft_started = True
    leading_edge_started = False
    lower_aft_started = False
    leading_edge_point = False
    af_le_pointdata = []
    af_upper_pointdata = []
    af_lower_pointdata = []
    for point in points:
        x = point[0]
        y = point[1]

        if x <= length_le and upper_aft_started:
            af_upper_pointdata.append([x, y, 0])
            leading_edge_started = True
            upper_aft_started = False

        if x >= length_le and leading_edge_started:
            af_le_pointdata.append([x, y, 0])
            leading_edge_started = False
            lower_aft_started = True

        if upper_aft_started:
            af_upper_pointdata.append([x, y, 0])
        elif leading_edge_started:
            if x == 0.0 and not leading_edge_point:
                af_le_pointdata.append([x, y, 0])
                leading_edge_point = True
            else:
                af_le_pointdata.append([x, y, 0])
        elif lower_aft_started:
            af_lower_pointdata.append([x, y, 0])

    gmsh.initialize()
    gmsh.model.add(airfoil.name)

    # hide output if specified
    if hide_output:
        gmsh.option.setNumber("General.Terminal", 0)

    # remove redundant points, will be added back later
    af_upper_pointdata.pop()
    af_lower_pointdata.pop()
    af_lower_pointdata.pop(0)

    model = gmsh.model.geo

    # Define points in airfoil in gmsh, create splines
    af_upper_points = [model.addPoint(x, y, z) for x, y, z in af_upper_pointdata]
    af_le_points = [model.addPoint(x, y, z) for x, y, z in af_le_pointdata]
    af_lower_points = [model.addPoint(x, y, z) for x, y, z in af_lower_pointdata]

    # Reinsert redundant points, set end points as equal to each other for meshing purposes
    af_top_point = af_le_points[0]
    af_bottom_point = af_le_points[-1]
    af_te_point = af_upper_points[0]
    af_upper_points.append(af_top_point)
    af_lower_points.insert(0, af_bottom_point)
    af_lower_points.append(af_te_point)

    # Create b splines for airfoil
    af_upper = model.addBSpline(af_upper_points)
    af_lower = model.addBSpline(af_lower_points)
    af_le = model.addBSpline(af_le_points)
    af = model.addCurveLoop([af_upper, af_le, af_lower])

    # Create Outer boundaries and surfaces

    # Inlet
    center_point = model.addPoint(center[0], center[1], 0)
    inlet_top_point = model.addPoint(center[0], r_inlet, 0)
    inlet_bottom_point = model.addPoint(center[0], -r_inlet, 0)
    front_line = model.addCircleArc(inlet_top_point, center_point, inlet_bottom_point)
    afTop_inletTop = model.addLine(af_top_point, inlet_top_point)
    inletBottom_afBottom = model.addLine(inlet_bottom_point, af_bottom_point)

    inlet_loop = model.addCurveLoop([front_line, inletBottom_afBottom, -af_le, afTop_inletTop])
    inlet_section = model.addPlaneSurface([inlet_loop])

    # Top Section
    top_te_point = model.addPoint(af_upper_pointdata[0][0], r_inlet, 0)
    top_line = model.addLine(top_te_point, inlet_top_point)
    topTe_afTe = model.addLine(top_te_point, af_te_point)

    top_loop = model.addCurveLoop([-af_upper, -afTop_inletTop, top_line, -topTe_afTe])
    top_section = model.addPlaneSurface([top_loop])

    # Bottom Section
    bottom_te_point = model.addPoint(af_upper_pointdata[0][0], -r_inlet, 0)
    bottom_line = model.addLine(inlet_bottom_point, bottom_te_point)
    afTe_bottomTe = model.addLine(af_te_point, bottom_te_point)

    bottom_loop = model.addCurveLoop([-inletBottom_afBottom, -af_lower, -afTe_bottomTe, bottom_line])
    bottom_section = model.addPlaneSurface([bottom_loop])

    # Top Wake Section
    top_wake_point = model.addPoint(downstream_distance, r_inlet, 0)
    center_wake_point = model.addPoint(downstream_distance, 0, 0)
    top_wake_line = model.addLine(top_wake_point, top_te_point)
    center_wake_line = model.addLine(af_te_point, center_wake_point)
    outlet_top = model.addLine(center_wake_point, top_wake_point)

    top_wake_loop = model.addCurveLoop([topTe_afTe, center_wake_line, outlet_top, top_wake_line])
    top_wake_section = model.addPlaneSurface([top_wake_loop])

    # Bottom Wake Section
    bottom_wake_point = model.addPoint(downstream_distance, -r_inlet, 0)
    outlet_bottom = model.addLine(bottom_wake_point, center_wake_point)
    bottom_wake_line = model.addLine(bottom_te_point, bottom_wake_point)

    bottom_wake_loop = model.addCurveLoop([-center_wake_line, outlet_bottom, bottom_wake_line, afTe_bottomTe])
    bottom_wake_section = model.addPlaneSurface([bottom_wake_loop])

    # flow_domain = model.addSurfaceLoop([inlet_section,bottom_section, bottom_wake_section, top_wake_section, bottom_wake_section])

    # Specify Meshing Parameters

    # Inlet Section
    airfoil.mesh_parameters.set_growth_factors()
    boundary_growth = airfoil.mesh_parameters.boundary_growth

    model.mesh.setTransfiniteCurve(front_line, n_leading_edge, "Bump", coef=-.1)
    model.mesh.setTransfiniteCurve(af_le, n_leading_edge)
    model.mesh.setTransfiniteCurve(afTop_inletTop, n_volume, "Progression", boundary_growth)
    model.mesh.setTransfiniteCurve(inletBottom_afBottom, n_volume, "Progression", -boundary_growth)
    model.mesh.setTransfiniteSurface(inlet_section)
    model.mesh.setRecombine(2, inlet_section)

    # Top Section
    te_growth = 1 / np.exp(np.log(te_thickness) / (n_airfoil - 1))
    te_growth = airfoil.mesh_parameters.te_growth
    model.mesh.setTransfiniteCurve(topTe_afTe, n_volume, "Progression", -boundary_growth)
    model.mesh.setTransfiniteCurve(af_upper, n_airfoil, "Progression", -te_growth)
    model.mesh.setTransfiniteCurve(top_line, n_airfoil, "Progression", -te_growth)
    model.mesh.setTransfiniteSurface(top_section)
    model.mesh.setRecombine(2, top_section)

    # Bottom Section
    model.mesh.setTransfiniteCurve(afTe_bottomTe, n_volume, "Progression", boundary_growth)
    model.mesh.setTransfiniteCurve(af_lower, n_airfoil, "Progression", te_growth)
    model.mesh.setTransfiniteCurve(bottom_line, n_airfoil, "Progression", -te_growth)
    model.mesh.setTransfiniteSurface(bottom_section)
    model.mesh.setRecombine(2, bottom_section)

    # Top Wake Section
    target_thickness = (1 - length_le) / n_airfoil
    target_thickness = airfoil.mesh_parameters.trailing_edge_thickness
    growth = 1 / np.exp(np.log(target_thickness) / (n_wake - 1))
    model.mesh.setTransfiniteCurve(top_wake_line, n_wake, "Progression", -growth)
    model.mesh.setTransfiniteCurve(center_wake_line, n_wake, "Progression", growth)
    model.mesh.setTransfiniteCurve(outlet_top, n_volume, "Progression", boundary_growth)
    model.mesh.setTransfiniteSurface(top_wake_section)
    model.mesh.setRecombine(2, top_wake_section)

    # Bottom Wake Section
    model.mesh.setTransfiniteCurve(bottom_wake_line, n_wake, "Progression", growth)
    model.mesh.setTransfiniteCurve(outlet_bottom, n_volume, "Progression", -boundary_growth)
    model.mesh.setTransfiniteSurface(bottom_wake_section)
    model.mesh.setRecombine(2, bottom_wake_section)

    model.addPhysicalGroup(1, [af_le, af_lower, af_upper], name='AIRFOIL')
    model.addPhysicalGroup(1, [outlet_top, outlet_bottom, top_line, top_wake_line, bottom_line, bottom_wake_line], name='OUTLET')
    model.addPhysicalGroup(1, [front_line], name='INLET')

    model.addPhysicalGroup(2, [inlet_section, top_section, top_wake_section, bottom_section, bottom_wake_section],
                           name='FLOWDOMAIN')

    # volume = model.extrude([(2, flow_domain)], 0, 0, 1)
    # gmsh.option.setNumber('Mesh.SaveAll', 1)
    model.synchronize()


    gmsh.model.mesh.generate(2)
    file_path = airfoil.name + '.su2'
    gmsh.write(file_path)
    gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
    if show_graphics:
        gmsh.fltk.run()
    gmsh.finalize()

    # Reformat the mesh file to be readable in SU2
    # GMSH forces you to add a named marker for the volume elements, use this to remove it
    # Create a temporary file to write the modified content
    temp_file_path = file_path + '.temp'
    with open(file_path, 'r') as source_file, open(temp_file_path, 'w') as temp_file:
        found_target_line = False
        for line in source_file:
            if line.startswith('NMARK'):
                line = 'NMARK= 3\n'
            if line.strip() == 'MARKER_TAG= FlowDomain':
                found_target_line = True
            if not found_target_line:
                temp_file.write(line)

    # Replace the original file with the temporary file
    os.replace(temp_file_path, file_path)


class mesh_parameters:
    """
    Stores variables used for meshing
    """
    def __init__(self):
        self.first_cell_thickness = 1e-5     # Desired thickness at first layer on boundary layer
        self.trailing_edge_thickness = 1e-2  # Desired cell thickness at trailing edge
        self.inlet_radius = 15               # Radius of the inlet, defines domain length to the front, top and bottom of airfoil
        self.downstream_distance = 25        # Domain length downstream of airfoil
        self.n_airfoil = 200               # Number of points to be used on the top and bottom surfaces of airfoil
        self.n_volume = 125                  # Number of points between airfoil and inlet
        self.n_wake = 150                    # Number of points behind the airfoil
        self.n_leading_edge = 100             # Number of points on the curved leading edge of airfoil
        self.leading_edge_length = .05        # Distance at which the leading edge of the airfoil is defined
        self.boundary_growth = 0
        self.te_growth = 0

    def set_growth_factors(self):
        self.boundary_growth = 1 / np.exp(np.log(self.first_cell_thickness) / (self.n_volume - 1))
        self.te_growth = (self.trailing_edge_thickness / .1) ** (1 / (self.n_airfoil - 1))
