import os
import subprocess
import shutil

def generate_PyAero_mesh(airfoil):
    write_airfoil_file(airfoil)

    fid = open('input_file.in', 'w')
    fid.write('cd PyAero\n')
    fid.write('python src/PyAero.py -no-gui data/Batch/batch_control.json\n')
    fid.close()

    mesh_file = 'airfoil.su2'

    try:
        subprocess.call('cmd.exe < input_file.in', shell=True, timeout=60, stdout=subprocess.DEVNULL)
    except:
        if os.path.exists(mesh_file):
            os.remove(mesh_file)
        return 0

    if os.path.exists(mesh_file):
        os.remove(mesh_file)

    current_path = os.getcwd()
    file_path = current_path + '/PyAero/data/OUTPUT/airfoil.su2'
    shutil.move(file_path, current_path)

    # reformat marker tags
    with open('airfoil.su2', 'r') as source_file, open('airfoil.temp', 'w') as temp_file:
        for line in source_file:
            if line.startswith('MARKER_TAG='):
                if line.strip() == 'MARKER_TAG= 1':
                    temp_file.write('MARKER_TAG= AIRFOIL\n')
                elif line.strip() == 'MARKER_TAG= 2':
                    temp_file.write('MARKER_TAG= INLET\n')
                elif line.strip() == 'MARKER_TAG= 3':
                    temp_file.write('MARKER_TAG= OUTLET\n')
            else:
                temp_file.write(line)
    os.replace('airfoil.temp', 'airfoil.su2')

    # Scale Mesh
    commands = f"SU2_DEF airfoil.cfg"
    proc = subprocess.Popen(['cmd', '/c', commands], shell=False, env=os.environ, stdin=subprocess.PIPE,
                            stdout=subprocess.DEVNULL, start_new_session=True)


def write_airfoil_file(airfoil):
    if os.path.exists('PyAero/data/Airfoils/airfoil/airfoil.dat'):
        os.remove('PyAero/data/Airfoils/airfoil/airfoil.dat')

    with open('PyAero/data/Airfoils/airfoil/airfoil.dat', 'w') as fid:
        fid = open('PyAero/data/Airfoils/airfoil/airfoil.dat', 'w')
        fid.write('#\n# Airfoil\n#\n')
        for point in airfoil.points:
            fid.write('{0}  {1}\n'.format(point[0], point[1]))

