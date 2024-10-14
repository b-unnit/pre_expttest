import qcportal as ptl
from molmass import Formula   #this is used to obtain the mass of the molecule
from posym import SymmetryMolecule #used to obtain the symmetry of the molecule from the xyz

def sampling_model_msg(
    mol_collection: str, target_mol: str, level_theory: str
) -> str:
    """
    Format a message for the start of sampling of a surface model.

    Args:
    - mol_collection: Name of the molecule collection.
    - target_mol: Name of the target molecule.
    - level_theory: Method and basis at wich the molecule is optimized.

    Returns:
    - Formatted message string.
    """
    return f"""
-----------------------------------------------------------------------------------------
Starting the calculation of the pre-exponential factor of the molecule {target_mol} in the collection {mol_collection} optimized at the {level_theory} level of theory
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
        "--client_address",
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
        help="The name of the molecule to be sampled (from a QCFractal OptimizationDataSet collection)",
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
        "--tag",
        type=str,
        default="tera_opt",
        help="The tag to use to specify the qcfractal-manager for the calculation (default: tera_opt)",
    )
    parser.add_argument(
        "--Range of temperatures",
        type=str,
        default="tera_opt",
        help="The tag to use to specify the qcfractal-manager for the calculation (default: tera_opt)",
    )

    return parser.parse_args()


def check_collection_existence(
    client: FractalClient,
    *collections: List,
    collection_type: str = "OptimizationDataset",
):
    """
    Check the existence of collections and raise DatasetNotFound error if not found.

    Args:
    - client: QCFractal client object
    - *collections: List of QCFractal Datasets.
    - collection_type: type of Optimization Dataset

    Raises:
    - DatasetNotFound: If any of the specified collections do not exist.
    """
    for collection in collections:
        try:
            client.get_collection(collection_type, collection)
        except KeyError:
            raise DatasetNotFound(
                f"Collection {collection} does not exist. Please create it first. Exiting..."
            )  


def check_optimized_molecule(
    ds: OptimizationDataset, opt_lot: str, mol_names: List[str]
) -> None:
    """
    Check if all molecules are optimized at the requested level of theory.

    Args:
    - ds: OptimizationDataset containing the optimization records.
    - opt_lot: Level of theory string.
    - mol_names: List of molecule names to check.

    Raises:
    - LevelOfTheoryNotFound: If the level of theory for a molecule or the entry itself does not exist.
    - ValueError: If optimizations are incomplete or encountered an error.
    """
    for mol in list(mol_names):
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
    client_adress, username: str, password: str, dataset: str, mol_name: str, level_theory: str, collection_type: str = "OptimizationDataset") -> List[str]:
    """
    Extract the xyz of the molecule

    Args:
    - client_adress: 
    - username:
    - password:
    - dataset: dataset containing the molecule.
    - mol_name: molecule name in the dataset.
    - level_theory: Level of theory at which the molecule is optimized.
    - collection_type: Type of optimization dataset (Default = OptimizationDataset.

    Returns:
    - XYZ file, excluding total number of atoms, charge, multiplicity and number of atom for each element present
    """
    client = ptl.FractalClient(address=client_address, username = username, password = password, verify=False)
    ds_opt = client.get_collection(collection_type,dataset)
    rr = ds_opt.get_record(mol_name, level_theory)
    mol = rr.get_final_molecule()    
    geom = mol.to_string(dtype="xyz")
    xyz_list = geom.splitlines()[2:]     
    xyz = '\n'.join(xyz_list)
    return(xyz)


def get_mass(
    mol: str) -> float:
    """
    Calculates the mass of a molecule. Utilizes the molmass package.
    """
    f = Formula(mol)
    return(f.mass)


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
    schema = qcel.models.Molecule.from_data(xyz).dict()
    mol = molsym.Molecule.from_schema(schema)
    pg, (paxis, saxis) = molsym.find_point_group(mol)

    point_group_to_sym_number = {
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
    return(point_group_to_sym_number.get(pg, "Point group not found... that should not happen"))
    
    
