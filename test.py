import sys, time, argparse, logging
import math
import numpy as np       # needed by molsym
import qcportal as ptl
from molmass import Formula   # this is used to obtain the mass of the molecule
import qcelemental as qcel # needed by molsym
import molsym # for the symmetry point

welcome_msg = """       
                  Welcome to the Pre-Exponential factor calculator! 

Description: The pre-ext uses the xxxxxx formula uses xxxxxxx

Author: b-unnit, namrata-rani10
            """


def calculation_msg(
    mol_col: str, mol: str, level_theory: str
) -> str:
    """
    Format a message for the start of the pre-exponential factor calculation.

    Args:
    - mol_col: Name of the molecule collection.
    - mol: Name of the target molecule.
    - level_theory: Method and basis at wich the molecule is optimized.

    Returns:
    - Formatted message string.
    """
    return f"""
-----------------------------------------------------------------------------------------
Starting the calculation of the pre-exponential factor of the molecule {mol} in the collection {mol_col} optimized at the {level_theory} level of theory
-----------------------------------------------------------------------------------------
    """


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments for the script.

    Returns:
    - Namespace containing the parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="""
A command line interface to calculate the pre-exponential factor of a given molecule inside of a collcetion.
This CLI is part of the Binding Energy Evaluation Platform (BEEP).
    """
    )

    parser.add_argument(
        "--client-address",
        default="localhost:7777",
        help="The URL address and port of the QCFractal server (default: localhost:7777)",
    )
    parser.add_argument(
        "--username",
        default=None,
        help="The username for the database client (Default = None)",
    )
    parser.add_argument(
        "--password",
        default=None,
        help="The password for the database client (Default = None)",
    )
    parser.add_argument(
        "--molecule",
        required=True,
        help="Molecule to be sampled (from a QCFractal OptimizationDataSet collection). Type 'all' to calculate all molecules in a collection",
    )
    parser.add_argument(
        "--molecule-collection",
        default="Small_molecules",
        help="The name of the collection containing molecules or radicals (default: Small_molecules)",
    )
    parser.add_argument(
        "--level-of-theory",
        default=["blyp_def2-svp"],
        help="The level of theory in which the molecule is optimized, in the format: method_basis (default: blyp_def2-svp)",
    )          
    parser.add_argument(
        "--temperature",
        default=["10"],
        help="Temperature for the calculation, in K (default: 10)",
    )  

  
    return parser.parse_args()


def check_collection_existence(
    client: FractalClient,
    collection: List,
    collection_type: str = "OptimizationDataset",
):
    """
    Check the existence of collections and raise DatasetNotFound error if not found.

    Args:
    - client: QCFractal client object
    - collection: QCFractal Datasets.
    - collection_type: type of Optimization Dataset

    Raises:
    - DatasetNotFound: If any of the specified collections do not exist.
    """
    try:
        client.get_collection(collection_type, collection)
    except KeyError:
        raise DatasetNotFound(
            f"Collection {collection} does not exist. Please create it first. Exiting..."
        )  


def check_optimized_molecule(
    ds: OptimizationDataset, opt_lot: str, mol: str
) -> None:
    """
    Check if all molecules are optimized at the requested level of theory.

    Args:
    - ds: OptimizationDataset containing the optimization records.
    - opt_lot: Level of theory string.
    - mol: Molecule name to check.

    Raises:
    - LevelOfTheoryNotFound: If the level of theory for a molecule or the entry itself does not exist.
    - ValueError: If optimizations are incomplete or encountered an error.
    """
    try:
        rr = ds.get_record(mol, opt_lot)
    except KeyError:
        raise LevelOfTheoryNotFound(
            f"{opt_lot} level of theory for {mol} or the entry itself does not exist in {ds.name} collection. Add the molecule and optimize it first\n"
        )
    if rr.status == "INCOMPLETE":
        raise ValueError(f" Optimization has status {rr.status} restart it or wait")
    elif rr.status == "ERROR":
        raise ValueError(f" Optimization has status {rr.status} restart it or wait")


def get_xyz(
    dataset: str, mol_name: str, level_theory: str, collection_type: str = "OptimizationDataset") -> str:
    """
    Extract the xyz of the molecule

    Args:
    - dataset: dataset containing the molecule.
    - mol_name: molecule name in the dataset.
    - level_theory: Level of theory at which the molecule is optimized.
    - collection_type: Type of optimization dataset (Default = OptimizationDataset.

    Returns:
    - XYZ file, excluding total number of atoms, charge, multiplicity and number of atom for each element present
    """
    #client = ptl.FractalClient(address=client_address, username = username, password = password, verify=False)
    ds_opt = client.get_collection(collection_type,dataset)
    rr = ds_opt.get_record(mol_name, level_theory)
    mol = rr.get_final_molecule()    
    geom = mol.to_string(dtype="xyz")
    xyz_list = geom.splitlines()[2:]     
    xyz = '\n'.join(xyz_list)
    return(xyz)


def get_mass(
    xyz: str) -> float:
    """
    Calculates the mass of a molecule. Utilizes the molmass package. Caps must be correctly used, so its extracted from the xyz.
    """
    lines = xyz.strip().splitlines()
    symbols = []
    for line in lines:
        parts = line.split()
        symbols.append(parts[0])
    mol_form = ''.join(symbols)
    f = Formula(mol_form)
    mass = f.mass * 1.66053906660e-27
    return(mass)


def sym_num(
    xyz: str):
    """
    Gives the symmetry number from the xyz (doesnt work with linear molecules)
    Symmetry numbers given by tables By P. W. ATKINS, M. S. CHILD, and C. S. G. PHILLIPS
    
    Args:
    - xyz: xyz file (only coordinates, no multiplicity, charge, etc)

    Returns:
    - Symmetry number
    """
    
    group_to_number = {
        "C1": 1, "Cs": 1, "Ch": 1, "Ci": 1, "S2": 1,
        "C2": 2, "C3": 3, "C4": 4, "C5": 5, "C6": 6, "C7": 7, "C8": 8,
        "D2": 4, "D3": 6, "D4": 8, "D5": 10, "D6": 12,
        "C2v": 2, "C3v": 3, "C4v": 4, "C5v": 5, "C6v": 6,
        "C2h": 2, "C3h": 3, "C4h": 4, "C5h": 5, "C6h": 6,
        "D2h": 4, "D3h": 6, "D4h": 8, "D5h": 10, "D6h": 12,
        "D2d": 4, "D3d": 6, "D4d": 8, "D5d": 10, "D6d": 12,
        "S4": 2, "S6": 3, "S8": 4,
        "T": 12, "Td": 12, "Th": 12, "O": 24, "Oh": 24,
        "I": 60, "Ih": 60,
    }

    schema = qcel.models.Molecule.from_data(xyz).dict()
    mol = molsym.Molecule.from_schema(schema)
    pg, (paxis, saxis) = molsym.find_point_group(mol)

    if pg not in group_to_number:
        raise Keyerror(f"Group {pg} not found in the dictionary... that shouldn't happen")

    return(group_to_number.get(pg))

def get_moments_of_inertia(xyz):
    """Calculate the principal moments of inertia from symbols and coordinates."""

    symbols = []
    coordinates = []

    # Split input into lines and extract symbol + coordinates from each line
    for line in xyz.strip().splitlines():
        parts = line.split()
        symbol = parts[0]
        x, y, z = map(float, parts[1:])
        symbols.append(symbol)
        coordinates.append([x, y, z])
    
    masses = np.array([get_atomic_mass(sym) for sym in symbols])
    coords = np.array(coordinates) * 1e-10

    # Calculate the center of mass
    total_mass = np.sum(masses)
    center_of_mass = np.sum(masses[:, np.newaxis] * coords, axis=0) / total_mass

    # Shift coordinates to the center of mass
    shifted_coords = coords - center_of_mass

    # Calculate the moment of inertia tensor
    I = np.zeros((3, 3))
    for m, r in zip(masses, shifted_coords):
        I[0, 0] += m * (r[1]**2 + r[2]**2)
        I[1, 1] += m * (r[0]**2 + r[2]**2)
        I[2, 2] += m * (r[0]**2 + r[1]**2)
        I[0, 1] -= m * r[0] * r[1]
        I[0, 2] -= m * r[0] * r[2]
        I[1, 2] -= m * r[1] * r[2]

    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]

    # Diagonalize the moment of inertia tensor to get principal moments
    eigenvalues, _ = np.linalg.eigh(I)
    Ia, Ib, Ic = np.sort(eigenvalues)  # Sort moments of inertia

    # Adjust for linear molecules
    if len(symbols) <= 3:  # Assuming linear for diatomic and triatomic
        Ia = 0  # Set the moment of inertia along the molecular axis to zero

    return Ia, Ib, Ic

def pre_exponential_factor(m, T, sigma, Ia, Ib, Ic):
    """Calculate the pre-exponential factor (v) for desorption."""
    kB = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34  # Planck's constant in J·s
    pi = math.pi

    # Translational contribution
    translational_part = ((2 * pi * m * kB * T) / h**2)**(3 / 2)

    # Rotational contribution (considering Ia = 0 for linear molecules)
    if Ia == 0:
        rotational_part = (8 * pi**2 * kB * T / h**2)**(3 / 2) * math.sqrt(Ib * Ic)
    else:
        rotational_part = (8 * pi**2 * kB * T / h**2)**(3 / 2) * math.sqrt(Ia * Ib * Ic)

    # Final pre-exponential factor
    v = ((kB * T) / h) * translational_part * (pi**0.5 / sigma) * rotational_part

    return v

def main():
    # Call the arguments
    args = parse_arguments()

    # Create a logger
    logger = logging.getLogger("preexpt_calc")
    logger.setLevel(logging.INFO)

    # File handler for logging to a file
    log_file = (
        "pre_exp_factor_" + args.molecule + "_" + args.level_of_theory + ".log"
    )
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(file_handler)

    logger.info(welcome_msg)

    client = ptl.FractalClient(  
        address=args.client_address,
        verify=False,
        username=args.username,
        password=args.password,
    )

    #Name of the molecule and level of theory    
    mol = args.molecule
    mol_col = args.molecule_collection
    mol_lot = args.level_of_theory
    temp = args.temperature
    
    #Check for collection existence
    check_collection_existence(client, mol_col)


    if mol == "all": 
        for molecule in mol_col?
        # Check if all the molecules are optimized at the requested level of theory
        check_optimized_molecule(mol_col, mol_lot, mol)
    
        logger.info(
            calculation_msg(mol_col, mol, mol_lot)
        )
                
        #WRITE THE FOR CYCLE
        
    else:   
        # Check if the molecule is optimized at the requested level of theory
        check_optimized_molecule(mol_col, mol_lot, mol)
    
        logger.info(
            calculation_msg(mol_col, mol, mol_lot)
        )    
    
        mol_xyz = get_xyz(mol_col,mol,mol_lot)  
        sym_num = sym_num(mol_xyz)
        mol_mass = get_mass(mol_xyz)
    
        Ia, Ib, Ic = get_moments_of_inertia(xyz)
        logger.info(f"Principal moments of inertia for {mol} (kg·m²): Ia={Ia:.3e}, Ib={Ib:.3e}, Ic={Ic:.3e}")
    
        v = pre_exponential_factor(mol_mass, temp, sym_num, Ia, Ib, Ic)
        logger.info(f"Pre-exponential factor for {mol} (v): {v:.3e} s⁻¹")


    if __name__ == "__main__":
    main()
