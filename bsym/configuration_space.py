from bsym.permutations import flatten_list, unique_permutations, number_of_unique_permutations
from bsym import Configuration, SymmetryGroup, SymmetryOperation
from bsym.partitions import generate_partitions, satisfies_bounds
import numpy as np
from itertools import combinations_with_replacement
from collections import Counter
from tqdm import tqdm, tqdm_notebook

class ConfigurationSpace:

    def __init__( self, objects, symmetry_group=None ):
        """
        Create a :any:`ConfigurationSpace` object.
  
        Args:
            objects (list): The set of objects that define the vector space of this configuration space.
            symmetry_group (:any:`SymmetryGroup`): The set of symmetry operations describing the symmetries of this configuration space.

        Returns:
            None
        """
        # Check that all properties have compatible dimensions
        self.dim = len( objects )
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

    def __repr__( self ):
        to_return = "ConfigurationSpace\n"
        to_return += self.objects.__repr__() + "\n"
        to_return += "\n".join( self.symmetry_group.__repr__().split("\n")[1:] )
        return to_return

    def enumerate_configurations(
        self,
        generator,
        verbose=False):
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
        working = True
        seen = set()
        unique_configurations = []
        using_tqdm = hasattr( generator, 'postfix' )
        for new_permutation in generator:
            if permutation_as_config_number( new_permutation ) not in seen:
                config = Configuration.from_tuple( new_permutation )
                numeric_equivalents = set( config.numeric_equivalents( self.symmetry_group.symmetry_operations ) )
                config.count = len(numeric_equivalents)
                [seen.add( i ) for i in numeric_equivalents]
                unique_configurations.append( config )
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
        generator = unique_permutations(s)
        if show_progress:
            if show_progress=='notebook':
                generator = tqdm_notebook(generator, total=total_permutations, unit=' permutations')
            else:
                generator = tqdm(generator, total=total_permutations, unit=' permutations')
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
    ) -> dict[tuple[int, int], list[Configuration]]:
        """
        Find all symmetry inequivalent configurations for all compositions within a range.
        
        Args:
            n_species (int): Number of species to distribute across sites.
            bounds (dict, optional): Occupancy bounds for each species. 
                                    Keys are species indices (int), values are (min, max) tuples.
                                    Species not listed default to (0, n_sites).
            verbose (bool): Print per-composition and summary information.
            show_progress (bool): Show progress bars. Can be True, False, or "notebook".
        
        Returns:
            dict: Mapping from composition tuples to lists of Configuration objects.
                Keys are tuples of length n_species with counts for each species.
                
        Example:
            For a 4-site system with 2 species:
            
            >>> config_space.unique_configurations_by_composition(n_species=2)
            {
                (4, 0): [Configuration(...)],
                (3, 1): [Configuration(...), Configuration(...)],
                (2, 2): [Configuration(...)],
                (1, 3): [Configuration(...), Configuration(...)],
                (0, 4): [Configuration(...)]
            }
        """
        n_sites = self.dim
        
        # Generate all partitions
        all_partitions = generate_partitions(n_sites, n_species)
        
        # Generate all composition tuples
        all_compositions = []
        for partition in all_partitions:
            for composition_tuple in unique_permutations(partition):
                all_compositions.append(composition_tuple)
        
        # Filter by bounds if provided
        if bounds is not None:
            valid_compositions = []
            for composition_tuple in all_compositions:
                composition_dict = {i: count for i, count in enumerate(composition_tuple)}
                if satisfies_bounds(composition_dict, bounds):
                    valid_compositions.append(composition_tuple)
        else:
            valid_compositions = all_compositions
        
        # Setup progress bar
        if show_progress:
            if show_progress == 'notebook':
                progress_bar = tqdm_notebook(
                    total=len(valid_compositions),
                    desc="Compositions",
                    unit=" compositions"
                )
            else:
                progress_bar = tqdm(
                    total=len(valid_compositions),
                    desc="Compositions", 
                    unit=" compositions"
                )
        
        # Container for results
        results = {}
        
        # Process each valid composition
        for composition_tuple in valid_compositions:
            
            # Verbose output
            if verbose:
                print(f"Processing composition {composition_tuple}...")
            
            # Create site_distribution (exclude zero counts)
            site_distribution = {
                species: count
                for species, count in enumerate(composition_tuple)
                if count > 0
            }
            
            # Get unique configurations for this composition
            unique_configs = self.unique_configurations(
                site_distribution=site_distribution,
                verbose=False,
                show_progress=show_progress
            )
            
            # Store results
            results[composition_tuple] = unique_configs
            
            # Verbose output
            if verbose:
                print(f"  Found {len(unique_configs)} unique configurations")
            
            # Update progress bar
            if show_progress:
                progress_bar.update(1)
                progress_bar.set_postfix(
                    total_configs=sum(len(configs) for configs in results.values())
                )
        
        # Close progress bar
        if show_progress:
            progress_bar.close()
        
        # Summary output
        if verbose:
            print(f"\nSummary:")
            print(f"  Evaluated {len(valid_compositions)} compositions")
            if bounds:
                print(f"  (filtered from {len(all_compositions)} total compositions)")
            print(f"  Total unique configurations: {sum(len(configs) for configs in results.values())}")
        
        return results

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
