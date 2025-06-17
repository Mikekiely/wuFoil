from pathlib import Path

def setup_output_paths():
    base = Path("outputs")
    mesh = base / "meshes"
    results = base / "results"
    airfoils = base / "airfoils"

    # create the directories if they don't exist
    mesh.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    airfoils.mkdir(parents=True, exist_ok=True)

    return base, mesh, results, airfoils