import os
import subprocess
import shutil
import time

def generate_PyAero_mesh(airfoil):
    write_airfoil_file(airfoil)

    # fid = open('input_file.in', 'w')
    # fid.write('cd PyAero\n')
    # fid.write('python src/PyAero.py -no-gui data/Batch/batch_control.json\n')
    # fid.close()
    commands = 'python src/PyAero.py -no-gui data/Batch/batch_control.json\n'
    mesh_file = f'{airfoil.name}.su2'

    try:
        # subprocess.call('cmd.exe < input_file.in', shell=True, timeout=60, stdout=subprocess.DEVNULL, cwd='PyAero')
        proc = subprocess.Popen(['cmd', '/c', commands], shell=True, stdin=subprocess.PIPE, start_new_session=True,
                                cwd='PyAero', stdout=subprocess.DEVNULL)
        stdout, stderr = proc.communicate()
    except:
        if os.path.exists(mesh_file):
            os.remove(mesh_file)
        return 0

    if os.path.exists(mesh_file):
        os.remove(mesh_file)

    current_path = os.getcwd()
    file_path = current_path + '/PyAero/data/OUTPUT/airfoil.su2'

    file_name = f"{airfoil.name}.su2"
    destination_path = airfoil.mesh_dir / file_name

    temp_path = str(destination_path) + '.temp'

    shutil.move(file_path, destination_path)

    # Read from source and write to temp
    with open(destination_path, 'r') as source_file:
        lines = []
        found_target_line = False
        for line in source_file:
            if line.startswith('MARKER_TAG='):
                if line.strip() == 'MARKER_TAG= 1':
                    line = 'MARKER_TAG= AIRFOIL\n'
                elif line.strip() == 'MARKER_TAG= 2':
                    line = 'MARKER_TAG= INLET\n'
                elif line.strip() == 'MARKER_TAG= 3':
                    line = 'MARKER_TAG= OUTLET\n'

            lines.append(line)
    with open(temp_path, 'w') as temp_file:
        temp_file.writelines(lines)
    time.sleep(.05)
    os.replace(temp_path, destination_path)

    # Scale Mesh
    commands = f"SU2_DEF airfoil.cfg"
    proc = subprocess.Popen(['cmd', '/c', commands], shell=False, env=os.environ, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL,
                            start_new_session=True)



def write_airfoil_file(airfoil):
    if os.path.exists('PyAero/data/Airfoils/airfoil/airfoil.dat'):
        os.remove('PyAero/data/Airfoils/airfoil/airfoil.dat')

    with open('PyAero/data/Airfoils/airfoil/airfoil.dat', 'w') as fid:
        fid = open('PyAero/data/Airfoils/airfoil/airfoil.dat', 'w')
        fid.write('#\n# Airfoil\n#\n')
        for point in airfoil.points:
            fid.write('{0}  {1}\n'.format(point[0], point[1]))

