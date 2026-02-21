from __future__ import annotations

from typing import Callable, Any
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, PointGroupAnalyzer
from pymatgen.util.coord import coord_list_mapping_pbc, coord_list_mapping
from pymatgen.core.structure import Molecule, Structure
from pymatgen.core.operations import SymmOp

from bsym import SpaceGroup, SymmetryOperation, ConfigurationSpace, PointGroup
from bsym.configuration import load_configurations, save_configurations
from functools import partial
import numpy as np


def structure_cartesian_coordinates_mapping(
    structure: Structure,
    symmop: SymmOp
) -> np.ndarray:
    """
    Maps the coordinates of pymatgen ``Structure`` according to a ``SymmOp`` symmetry operation.

    Args:
        structure: The pymatgen ``Structure``.
        symmop: The pymatgen symmetry operation object.

    Returns:
        The mapped Cartesian coordinates.
    """
    return structure.lattice.get_cartesian_coords(symmop.operate_multi(structure.frac_coords))


def molecule_cartesian_coordinates_mapping(
    molecule: Molecule,
    symmop: SymmOp
) -> np.ndarray:
    """
    Maps the coordinates of pymatgen ``Molecule`` according to a ``SymmOp`` symmetry operation.

    Args:
        molecule: The pymatgen ``Molecule``.
        symmop: The pymatgen symmetry operation object.

    Returns:
        The mapped Cartesian coordinates.
    """
    return symmop.operate_multi(molecule.cart_coords)


def structure_mapping_list(
    new_structure: Structure,
    mapping_structure: Structure,
    atol: float
) -> list[int]:
    """
    Gives the index mapping between two pymatgen ``Structure`` objects.

    Args:
        new_structure: The new structure.
        mapping_structure: The structure to map to.
        atol: Absolute tolerance for coordinate matching.

    Returns:
        List of indices such that mapping_structure.sites[indices] == new_structure.sites
    """
    return coord_list_mapping_pbc(  # type: ignore[no-any-return]
        new_structure.frac_coords,
        mapping_structure.frac_coords,
        atol=atol
    )


def molecule_mapping_list(
    new_molecule: Molecule,
    mapping_molecule: Molecule,
    atol: float
) -> list[int]:
    """
    Gives the index mapping between two pymatgen ``Molecule`` objects.

    Args:
        new_molecule: The new molecule.
        mapping_molecule: The molecule to map to.
        atol: Absolute tolerance for coordinate matching.

    Returns:
        List of indices such that mapping_molecule.sites[indices] == new_molecule.sites
    """
    return coord_list_mapping(  # type: ignore[no-any-return]
        new_molecule.cart_coords,
        mapping_molecule.cart_coords,
        atol=atol
    )


def unique_symmetry_operations_as_vectors_from_structure(
    structure: Structure | Molecule,
    verbose: bool = False,
    subset: list[int] | None = None,
    atol: float = 1e-5
) -> list[list[int]]:
    """
    Uses `pymatgen`_ symmetry analysis to find the minimum complete set of symmetry operations 
    for the space group of a structure.

    Args:
        structure: Structure or Molecule to be analysed.
        verbose: Print verbose output including space/point group information.
        subset: List of atom indices to be used for generating the symmetry operations.
        atol: Tolerance factor for the ``pymatgen`` `coordinate mapping`_ under each 
            symmetry operation.

    Returns:
        A list of lists, containing the symmetry operations as vector mappings.

    .. _pymatgen:
        http://pymatgen.org
    .. _coordinate mapping:
        http://pymatgen.org/pymatgen.util.coord_utils.html#pymatgen.util.coord_utils.coord_list_mapping_pbc
    """
    if isinstance(structure, Structure):
        instantiate_structure: Callable[..., Structure | Molecule] = partial(
            Structure, lattice=structure.lattice, coords_are_cartesian=True
        )
        coord_mapping: Callable[[Structure | Molecule, SymmOp], np.ndarray] = structure_cartesian_coordinates_mapping  # type: ignore[assignment]
        mapping_list: Callable[[Structure | Molecule, Structure | Molecule, float], list[int]] = structure_mapping_list  # type: ignore[assignment]
        symmetry_analyzer: SpacegroupAnalyzer | PointGroupAnalyzer = SpacegroupAnalyzer(structure)
        if verbose:
            if isinstance(symmetry_analyzer, SpacegroupAnalyzer):
                print(f"The space group for this structure is {symmetry_analyzer.get_space_group_symbol()}")
    elif isinstance(structure, Molecule):
        instantiate_structure = Molecule
        coord_mapping = molecule_cartesian_coordinates_mapping  # type: ignore[assignment]
        mapping_list = molecule_mapping_list  # type: ignore[assignment]
        symmetry_analyzer = PointGroupAnalyzer(structure, tolerance=atol)
        if verbose:
            print(f"The point group for this structure is {symmetry_analyzer.sch_symbol}")
    else:
        raise ValueError('structure argument should be a Structure or Molecule object')
    
    symmetry_operations = symmetry_analyzer.get_symmetry_operations()
    seen: set[tuple[int, ...]] = set()
    mappings: list[list[int]] = []

    mapping_structure: Structure | Molecule
    if subset:
        species_subset = [spec for i, spec in enumerate(structure.species) if i in subset]
        cart_coords_subset = [coord for i, coord in enumerate(structure.cart_coords) if i in subset]
        mapping_structure = instantiate_structure(species=species_subset, coords=cart_coords_subset)
    else:
        mapping_structure = structure

    for symmop in symmetry_operations:
        cart_coords = coord_mapping(mapping_structure, symmop)
        new_structure = instantiate_structure(species=mapping_structure.species, coords=cart_coords)
        new_mapping = [x + 1 for x in list(mapping_list(new_structure, mapping_structure, atol))]
        key = tuple(new_mapping)
        if key not in seen:
            seen.add(key)
            mappings.append(new_mapping)

    return mappings


def space_group_symbol_from_structure(structure: Structure) -> str:
    """
    Returns the symbol for the space group defined by this structure.

    Args:
        structure: The input structure.

    Returns:
        The space group symbol.
    """
    symmetry_analyzer = SpacegroupAnalyzer(structure)
    symbol = symmetry_analyzer.get_space_group_symbol()
    return symbol


def space_group_from_structure(
    structure: Structure,
    subset: list[int] | None = None,
    atol: float = 1e-5
) -> SpaceGroup:
    """
    Generates a ``SpaceGroup`` object from a `pymatgen` ``Structure``.

    Args:
        structure: Structure to be used to define the :any:`SpaceGroup`.
        subset: List of atom indices to be used for generating the symmetry operations.
        atol: Tolerance factor for the ``pymatgen`` `coordinate mapping`_ under each 
            symmetry operation.

    Returns:
        A new :any:`SpaceGroup` instance

    .. _coordinate mapping:
        http://pymatgen.org/pymatgen.util.coord_utils.html#pymatgen.util.coord_utils.coord_list_mapping_pbc
    """
    mappings = unique_symmetry_operations_as_vectors_from_structure(structure, subset=subset, atol=atol)
    symmetry_operations = [SymmetryOperation.from_vector(m) for m in mappings]
    return SpaceGroup(symmetry_operations=symmetry_operations)


def point_group_from_molecule(
    molecule: Molecule,
    subset: list[int] | None = None,
    atol: float = 1e-5
) -> PointGroup:
    """
    Generates a ``PointGroup`` object from a `pymatgen` ``Molecule``.

    Args:
        molecule: Molecule to be used to define the :any:`PointGroup`.
        subset: List of atom indices to be used for generating the symmetry operations.
        atol: Tolerance factor for the ``pymatgen`` `coordinate mapping`_ under each 
            symmetry operation.

    Returns:
        A new :any:`PointGroup` instance

    .. _coordinate mapping:
        http://pymatgen.org/pymatgen.util.coord_utils.html#pymatgen.util.coord_utils.coord_list_mapping
    """
    molecule = Molecule(molecule.species, molecule.cart_coords - molecule.center_of_mass)
    mappings = unique_symmetry_operations_as_vectors_from_structure(molecule, subset=subset, atol=atol)
    symmetry_operations = [SymmetryOperation.from_vector(m) for m in mappings]
    return PointGroup(symmetry_operations=symmetry_operations)


def configuration_space_from_structure(
    structure: Structure,
    subset: list[int] | None = None,
    atol: float = 1e-5
) -> ConfigurationSpace:
    """
    Generate a ``ConfigurationSpace`` object from a `pymatgen` ``Structure``.

    Args:
        structure: Structure to be used to define the :any:`ConfigurationSpace`.
        subset: List of atom indices to be used for generating the configuration space.
        atol: Tolerance factor for the ``pymatgen`` `coordinate mapping`_ under each 
            symmetry operation.

    Returns:
        A new :any:`ConfigurationSpace` instance.

    .. _coordinate mapping:
        http://pymatgen.org/pymatgen.util.coord_utils.html#pymatgen.util.coord_utils.coord_list_mapping_pbc
    """
    space_group = space_group_from_structure(structure, subset=subset, atol=atol)
    if subset is None:
        subset = list(range(1, len(structure) + 1))
    config_space = ConfigurationSpace(objects=subset, symmetry_group=space_group)
    return config_space


def configuration_space_from_molecule(
    molecule: Molecule,
    subset: list[int] | None = None,
    atol: float = 1e-5
) -> ConfigurationSpace:
    """
    Generate a ``ConfigurationSpace`` object from a `pymatgen` ``Molecule``.

    Args:
        molecule: Molecule to be used to define the :any:`ConfigurationSpace`.
        subset: List of atom indices to be used for generating the configuration space.
        atol: Tolerance factor for the ``pymatgen`` `coordinate mapping`_ under each 
            symmetry operation.

    Returns:
        A new :any:`ConfigurationSpace` instance.

    .. _coordinate mapping:
        http://pymatgen.org/pymatgen.util.coord_utils.html#pymatgen.util.coord_utils.coord_list_mapping
    """
    point_group = point_group_from_molecule(molecule, subset=subset, atol=atol)
    if subset is None:
        subset = list(range(1, len(molecule) + 1))
    config_space = ConfigurationSpace(objects=subset, symmetry_group=point_group)
    return config_space


def unique_structure_substitutions(
    structure: Structure | Molecule,
    to_substitute: str,
    site_distribution: dict[str, int],
    verbose: bool = False,
    atol: float = 1e-5,
    show_progress: bool | str = False
) -> list[Any]:
    """
    Generate all symmetry-unique structures formed by substituting a set of sites in a 
    `pymatgen` structure.

    Args:
        structure: The parent structure.
        to_substitute: Atom label for the sites to be substituted.
        site_distribution: A dictionary that defines the number of each substituting element.
        verbose: Verbose output.
        atol: Tolerance factor for the ``pymatgen`` `coordinate mapping`_ under each 
            symmetry operation. Default=1e-5.
        show_progress: Show a progress bar. Setting to `True` gives a simple progress bar.
            Setting to `"notebook"` gives a Jupyter notebook compatible progress bar.

    Returns:
        A list of Structure objects for each unique substitution.

    Notes:
        The number of symmetry-equivalent configurations for each structure 
        is stored in the `number_of_equivalent_configurations` attribute.

        If the parent structure was previously generated using this function
        (as part of a sequence of substitutions) the full configuration
        degeneracy of each symmetry inequivalent configuration is stored in
        the `full_configuration_degeneracy` attribute. If the parent structure
        is a standard Pymatgen Structure object, `number_of_equivalent_configurations`
        and `full_configuration_degeneracy` will be equal.

    .. _coordinate mapping:
        http://pymatgen.org/pymatgen.util.coord_utils.html#pymatgen.util.coord_utils.coord_list_mapping_pbc
    """
    site_substitution_index = list(structure.indices_from_symbol(to_substitute))
    if len(site_substitution_index) != sum(site_distribution.values()):
        raise ValueError("Number of sites from index does not match number from site distribution")
    
    if isinstance(structure, Structure):
        config_space = configuration_space_from_structure(structure, subset=site_substitution_index, atol=atol)
    elif isinstance(structure, Molecule):
        structure = Molecule(structure.species, structure.cart_coords - structure.center_of_mass)
        config_space = configuration_space_from_molecule(structure, subset=site_substitution_index, atol=atol)
    else:
        raise ValueError("pymatgen Structure or Molecule object expected")
    
    numeric_site_distribution, numeric_site_mapping = parse_site_distribution(site_distribution)
    unique_configurations = config_space.unique_configurations(
        numeric_site_distribution, verbose=verbose, show_progress=show_progress
    )
    new_structures = [
        new_structure_from_substitution(
            structure, site_substitution_index, [numeric_site_mapping[k] for k in c.tolist()]
        )
        for c in unique_configurations
    ]
    
    if hasattr(structure, 'number_of_equivalent_configurations'):
        for s, c in zip(new_structures, unique_configurations):
            s.number_of_equivalent_configurations = c.count 
            s.full_configuration_degeneracy = c.count * structure.full_configuration_degeneracy  # type: ignore[union-attr]
    else:
        for s, c in zip(new_structures, unique_configurations):
            s.number_of_equivalent_configurations = c.count
            s.full_configuration_degeneracy = c.count
    
    return new_structures


def parse_site_distribution(
    site_distribution: dict[str, int]
) -> tuple[dict[int, int], dict[int, str]]:
    """
    Converts a site distribution using species labels into one using integer labels.

    Args:
        site_distribution: e.g. `{'Mg': 1, 'Li': 3}`

    Returns:
        Tuple of (numeric_site_distribution, numeric_site_mapping)
        e.g. ({0: 3, 1: 1}, {0: 'Li', 1: 'Mg'})
    """
    numeric_site_distribution: dict[int, int] = {}
    numeric_site_mapping: dict[int, str] = {}
    for i, k in enumerate(site_distribution.keys()):
        numeric_site_distribution[i] = site_distribution[k]
        numeric_site_mapping[i] = k
    return numeric_site_distribution, numeric_site_mapping


def new_structure_from_substitution(
    parent_structure: Structure | Molecule,
    site_substitution_index: list[int],
    new_species_list: list[str]
) -> Any:
    """
    Generate a new pymatgen ``Structure`` from site substitution parameters.

    Args:
        parent_structure: The parent pymatgen ``Structure`` or ``Molecule`` object.
        site_substitution_index: The list of site indices to be substituted.
        new_species_list: A list of the replacement atomic species.

    Returns:
        The new pymatgen ``Structure`` or ``Molecule``.

    Notes:
        pymatgen ``Structure`` and ``Molecule`` classes both subclass ``SiteCollection``.
        This function will also accept a parent ``Molecule`` object, and return a new
        ``Molecule``.
    """
    if len(site_substitution_index) != len(new_species_list):
        raise ValueError("site_substitution_index and new_species_list must have same length")
    if any(i >= len(parent_structure) for i in site_substitution_index):
        raise ValueError("site_substitution_index contains invalid site indices")
    
    s = parent_structure.copy()
    for i, spec in zip(site_substitution_index, new_species_list):
        s[i] = spec
    return s


def unique_structure_substitutions_by_composition(
    structure: Structure | Molecule,
    to_substitute: str,
    species_list: list[str],
    bounds: dict[str, tuple[int | None, int | None]] | None = None,
    verbose: bool = False,
    atol: float = 1e-5,
    show_progress: bool | str = False
) -> dict[tuple[int, ...], list[Any]]:
    """
    Generate all symmetry-unique structures for all compositions of substituting species.

    Args:
        structure: The parent structure.
        to_substitute: Atom label for the sites to be substituted.
        species_list: List of species to substitute, e.g., ['Li', 'Na', 'Mg'].
            Order determines composition tuple indices.
        bounds: Occupancy bounds for each species. Keys are species names, values are 
            (min, max) tuples. e.g., {'Li': (1, 3), 'Na': (0, 2)}
        verbose: Verbose output.
        atol: Tolerance factor for coordinate mapping. Default=1e-5.
        show_progress: Show progress bars. Can be True, False, or "notebook".

    Returns:
        Mapping from composition tuples to lists of Structure objects.
        Keys are tuples like (2, 1, 1) corresponding to species_list order.
        Each Structure has a `number_of_equivalent_configurations` attribute.

    Example:
        >>> results = unique_structure_substitutions_by_composition(
        ...     structure, 'X', ['Li', 'Na'])
        >>> # results[(2, 2)] gives structures with 2 Li and 2 Na
    """
    # Get sites to substitute
    site_substitution_index = list(structure.indices_from_symbol(to_substitute))
    n_sites = len(site_substitution_index)

    # Create configuration space
    if isinstance(structure, Structure):
        config_space = configuration_space_from_structure(structure, subset=site_substitution_index, atol=atol)
    elif isinstance(structure, Molecule):
        structure = Molecule(structure.species, structure.cart_coords - structure.center_of_mass)
        config_space = configuration_space_from_molecule(structure, subset=site_substitution_index, atol=atol)
    else:
        raise ValueError("pymatgen Structure or Molecule object expected")

    # Create species name to index mapping
    species_to_index = {species: i for i, species in enumerate(species_list)}

    # Convert bounds from species names to indices if provided
    bounds_numeric: dict[int, tuple[int | None, int | None]] | None = None
    if bounds is not None:
        bounds_numeric = {}
        for species, (min_val, max_val) in bounds.items():
            if species not in species_to_index:
                raise ValueError(f"Species '{species}' in bounds not found in species_list")
            index = species_to_index[species]
            bounds_numeric[index] = (min_val, max_val)

    # Get unique configurations by composition
    configs_by_composition = config_space.unique_configurations_by_composition(
        n_species=len(species_list),
        bounds=bounds_numeric,
        verbose=verbose,
        show_progress=show_progress
    )

    # Convert configurations to structures
    results: dict[tuple[int, ...], list[Any]] = {}
    for composition_tuple, configurations in configs_by_composition.items():
        structures: list[Any] = []
        for config in configurations:
            # Map configuration indices to species names
            species_for_sites = [species_list[species_idx] for species_idx in config.tolist()]

            # Create new structure
            new_structure = new_structure_from_substitution(
                structure,
                site_substitution_index,
                species_for_sites
            )

            # Add metadata
            new_structure.number_of_equivalent_configurations = config.count

            structures.append(new_structure)

        results[composition_tuple] = structures

    return results
    
def random_unique_structure_substitutions(
    structure,
    to_substitute,
    site_distribution,
    n,
    sampling='degeneracy_weighted',
    seed=None,
    atol=1e-5,
    exclude_file=None,
    output_file=None,
):
    """
    Generate n random symmetry-unique structures by substituting sites in a pymatgen structure.

    Args:
        structure (pymatgen.Structure): The parent structure.
        to_substitute (str): Atom label for the sites to be substituted.
        site_distribution (dict): Dictionary mapping species to counts, e.g. {'O': 8, 'F': 8}.
        n (int): Number of unique structures to generate.
        sampling (str): Sampling method. Either 'degeneracy_weighted' (default) or 'uniform'.
            'degeneracy_weighted' samples configurations with probability proportional
            to their degeneracy. 'uniform' samples uniformly over equivalence classes.
        seed (int, optional): Random seed for reproducibility.
        atol (float): Tolerance factor for coordinate mapping. Default=1e-5.
        exclude_file (str or list[str], optional): Path(s) to JSON file(s) of configurations to exclude.
        output_file (str, optional): Path to save generated configurations.

    Returns:
        list[Structure]: A list of n unique Structure objects. Each has a
            `number_of_equivalent_configurations` attribute.
    """
    site_substitution_index = list(structure.indices_from_symbol(to_substitute))
    n_sites = len(site_substitution_index)
    
    config_space = configuration_space_from_structure(
        structure,
        subset=site_substitution_index,
        atol=atol
    )
    
    numeric_site_distribution, numeric_site_mapping = parse_site_distribution(
        site_distribution
    )
    
    # Load excluded configurations
    exclude = None
    if exclude_file is not None:
        exclude = []
        if isinstance(exclude_file, str):
            exclude_file = [exclude_file]
        for filename in exclude_file:
            loaded = load_configurations(filename)
            for config in loaded:
                if len(config.tolist()) != n_sites:
                    raise ValueError(
                        f"Configuration length {len(config.tolist())} in '{filename}' "
                        f"does not match number of sites {n_sites}"
                    )
            exclude.extend(loaded)
    
    configurations = config_space.random_unique_configurations(
        site_distribution=numeric_site_distribution,
        n=n,
        sampling=sampling,
        seed=seed,
        exclude=exclude,
    )
    
    # Save configurations if requested
    if output_file is not None:
        save_configurations(configurations, output_file)
    
    unique_structures = []
    for config in configurations:
        species_for_sites = [numeric_site_mapping[i] for i in config.tolist()]
        new_structure = new_structure_from_substitution(
            structure,
            site_substitution_index,
            species_for_sites
        )
        new_structure.number_of_equivalent_configurations = config.count
        unique_structures.append(new_structure)
    
    return unique_structures