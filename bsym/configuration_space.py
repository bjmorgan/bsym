from bsym.permutations import flatten_list, unique_permutations, number_of_unique_permutations
from bsym.partitions import compute_mapping_vector
from bsym import Configuration, SymmetryGroup, SymmetryOperation
from bsym.partitions import generate_partitions, satisfies_bounds
import numpy as np
from itertools import combinations_with_replacement
from collections import Counter
from tqdm import tqdm, tqdm_notebook
from typing import Iterator


class ConfigurationSpace:

    def __init__(self,
        objects: list,
        symmetry_group: SymmetryGroup | None = None) -> None:
        """
        Create a :any:`ConfigurationSpace` object.
  
        Args:
            objects (list): The set of objects that define the vector space of this configuration space.
            symmetry_group (:any:`SymmetryGroup`): The set of symmetry operations describing the symmetries of this configuration space.

        Returns:
            None
        """
        # Check that all properties have compatible dimensions
        self.dim = len(objects)
        self.objects = objects
        if symmetry_group:
            for so in symmetry_group.symmetry_operations:
                if so.matrix.shape[0] != self.dim:
                    raise ValueError
            self.symmetry_group = symmetry_group
        else:
            self.symmetry_group = SymmetryGroup(
                symmetry_operations=[
                    SymmetryOperation(np.identity(self.dim, dtype=int), label='E')
                ]
            )

    def __repr__(self) -> str:
        to_return: str
        to_return = "ConfigurationSpace\n"
        to_return += self.objects.__repr__() + "\n"
        to_return += "\n".join(self.symmetry_group.__repr__().split("\n")[1:])
        return to_return

    def enumerate_configurations(self, generator, verbose=False):
        """
        Find all symmetry inequivalent configurations within the set produced by
        `generator`.
        
        Args:
            generator (:obj:`generator`): Generator object, that yields the configurations
                to search through.
            verbose (opt:default=False): Print verbose output.
        
        Returns:
            unique_configurations (list): A list of :any:`Configuration` objects, for each symmetry
                                        inequivalent configuration. 
        """
        seen = set()
        unique_configurations = []
        using_tqdm = hasattr(generator, 'postfix')
        
        for new_permutation in generator:
            perm_as_bytes = Configuration.tuple_to_bytes(new_permutation)
            if perm_as_bytes not in seen:
                config = Configuration.from_tuple(new_permutation)
                byte_equivalents = config.get_byte_equivalents(self.symmetry_group)
                config.count = len(byte_equivalents)
                seen.update(byte_equivalents)
                unique_configurations.append(config)
                if using_tqdm:
                    generator.set_postfix(found=len(unique_configurations))
        
        if verbose:
            print('unique configurations: {} / {}'.format(len(unique_configurations), len(seen)))
        
        return unique_configurations

    def unique_configurations(self,
        site_distribution,
        verbose=False,
        show_progress=False):
        """
        Find the symmetry inequivalent configurations for a given population of objects.

        Args:
            site_distribution (dict): A dictionary that defines the number of each object 
                                      to be arranged in this system.

                                      e.g. for a system with four sites, with two occupied (denoted `1`)
                                      and two unoccupied (denoted `0`)::

                                          { 1: 2, 0: 2 }
            verbose (opt:default=False): Print verbose output.
            show_progress (opt:default=False): Show a progress bar.
                                      Setting to `True` gives a simple progress bar.
                                      Setting to `"notebook"` gives a Jupyter notebook compatible progress bar.

        Returns:
            unique_configurations (list): A list of :any:`Configuration` objects, for each symmetry 
                                          inequivalent configuration. 
        """
        s = flatten_list([[key] * site_distribution[key] for key in site_distribution])
        total_permutations = number_of_unique_permutations(s)
        if verbose:
            print('total number of sites: ' + str( sum( site_distribution.values())))
            print('using {:d} symmetry operations.'.format( len( self.symmetry_group.symmetry_operations)))
            print('evaluating {:d} unique permutations.'.format( total_permutations))
        generator: Iterator[tuple[int, ...]] = unique_permutations(s)
        if show_progress:
            TqdmClass = tqdm_notebook if show_progress == 'notebook' else tqdm
            generator = TqdmClass( # type: ignore[assignment]
                generator,
                total=total_permutations,
                unit=' permutations',
                mininterval=0.1
            )  
        return self.enumerate_configurations(generator, verbose=verbose)

    def unique_colourings(self, colours, verbose=False):
        """
        Find the symmetry inequivalent colourings for a given number of 'colours'.

        Args:
            colours (list): A list of each object that may be arranged zero or more times in this system.
            verbose (opt:default=False): Print verbose output.

        Returns:
            unique_colours (list): A list of :any:`Configuration` objects, for each symmetry
                                   inequivalent colouring.
        """
        generator = colourings_generator( colours, self.dim )
        return self.enumerate_configurations( generator, verbose=verbose )
        
    def unique_configurations_by_composition(self,
        n_species: int,
        bounds: dict[int, tuple[int|None, int|None]] | None = None,
        verbose: bool = False,
        show_progress: bool | str = False
    ) -> dict[tuple[int, ...], list[Configuration]]:
        """[docstring unchanged]"""
        from bsym.partitions import compute_mapping_vector
        
        n_sites = self.dim
        all_partitions = generate_partitions(n_sites, n_species)
        
        # Initialize progress bar without pre-counting (avoids iterator exhaustion)
        if show_progress:
            TqdmClass = tqdm_notebook if show_progress == 'notebook' else tqdm
            progress_bar = TqdmClass(
                desc="Compositions",
                unit=" compositions",
                mininterval=0.1
            )
        
        results = {}
        partitions_analyzed = 0
        
        for partition in all_partitions:
            canonical = partition
            
            # Get all permutations for this partition
            all_perms = list(unique_permutations(partition))
            
            # Filter by bounds
            valid_perms = []
            for perm in all_perms:
                composition_dict = {i: count for i, count in enumerate(perm)}
                if bounds is None or satisfies_bounds(composition_dict, bounds):
                    valid_perms.append(perm)
            
            if not valid_perms:
                continue
            
            partitions_analyzed += 1
            
            if verbose and not show_progress:
                print(f"Processing partition {partition}...")
            
            # Build site_distribution for canonical
            site_distribution = {
                species: count
                for species, count in enumerate(canonical)
                if count > 0
            }
            
            # Analyze canonical - pass through show_progress
            canonical_configs = self.unique_configurations(
                site_distribution=site_distribution,
                verbose=False,
                show_progress=show_progress  # Changed: pass through
            )
            
            if verbose and not show_progress:
                print(f"  Found {len(canonical_configs)} unique configurations")
            
            # Add results for each valid permutation
            for perm in valid_perms:
                if perm == canonical:
                    results[canonical] = canonical_configs
                else:
                    mapping = compute_mapping_vector(canonical, perm)
                    relabeled = [apply_species_mapping(config, mapping) 
                               for config in canonical_configs]
                    results[perm] = relabeled
                
                if show_progress:
                    progress_bar.update(1)
        
        if show_progress:
            progress_bar.close()
        
        if verbose:
            print(f"\nSummary:")
            print(f"  Analyzed {partitions_analyzed} partitions")
            print(f"  Generated {len(results)} compositions")
            print(f"  Total unique configurations: {sum(len(configs) for configs in results.values())}")
        
        return results

def apply_species_mapping(config, mapping_vector):
    """
    Apply species permutation to a Configuration.
    
    Creates a new Configuration where each species index is remapped according
    to the mapping vector. This is used to generate configurations for non-canonical
    compositions by relabeling species from canonical composition results.
    
    Args:
        config (Configuration): Configuration object with species indices.
        mapping_vector (list[int]): 0-indexed list where mapping_vector[i] gives 
                                    the new species index for current species i.
                                    e.g., [1, 0] swaps species 0 and 1.
    
    Returns:
        Configuration: New Configuration with relabeled species. The count
                      (degeneracy) is preserved from the original configuration.
    
    Example:
        >>> config = Configuration([0, 0, 1])
        >>> config.count = 2
        >>> mapping = [1, 0]  # Swap species 0 and 1
        >>> result = apply_species_mapping(config, mapping)
        >>> list(result)
        [1, 1, 0]
        >>> result.count
        2
    """
    new_config_list = [mapping_vector[species_idx] for species_idx in config.vector]
    new_config = Configuration(new_config_list)
    new_config.count = config.count
    return new_config
        
def colourings_generator( colours, dim ):
    for s in combinations_with_replacement( colours, dim ):
        for new_permutation in unique_permutations( s ):
            yield new_permutation 
        
def permutation_as_config_number(p):
    """
    A numeric representation of a numeric list.

    Example:
        >>> permutation_as_config_number( [ 1, 1, 0, 0, 1 ] )
        11001
    """
    tot = 0
    for num in p:
        tot *= 10
        tot += num
    return tot
