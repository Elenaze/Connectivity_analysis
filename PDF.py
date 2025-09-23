import numpy as np
import pandas as pd
from datetime import datetime   
from scipy import stats, signal
from pathlib import Path
from typing import Dict, List, Tuple, Set
from results import AnalysisResults
from matplotlib import pyplot as plt
import seaborn as sns
from ase.data import covalent_radii, chemical_symbols



class PDFAnalyzer:
    def __init__(self, xyz_files):
        # Handle both single file (string) and multiple files (list)
        if isinstance(xyz_files, (str, Path)):
            self.xyz_files = [Path(xyz_files)]
        else:
            self.xyz_files = [Path(f) for f in xyz_files]
            
        self.systems: Dict[int, pd.DataFrame] = {}
        self.unique_atoms: Set[str] = set()
        self.pair_distances: Dict[Tuple[str, str], List[float]] = {}
        self.kde_data: Dict[Tuple[str, str], pd.DataFrame] = {}
        self.discarded_pairs: List[Tuple[str, str]] = []
        self.X_max: Dict[Tuple[str, str], float] = {} 
        self.X_min: Dict[Tuple[str, str], float] = {}
        
        # Reference radii from ASE
        self.ASE_radii = {symbol: covalent_radii[Z] 
                  for Z, symbol in enumerate(chemical_symbols) if Z > 0}
        
    def load_systems(self) -> None:
        """Load all systems from XYZ files into DataFrames"""
        system_counter = 0
        
        for xyz_file in self.xyz_files:
            # Check if it's a multi-structure file by looking for multiple atom count lines
            if self._is_multi_structure_file(xyz_file):
                system_counter = self._load_multi_xyz(xyz_file, system_counter)
            else:
                # Handle single-structure XYZ files
                self._load_single_xyz(xyz_file, system_counter)
                system_counter += 1
    
    def _load_single_xyz(self, xyz_file: Path, system_id: int) -> None:
        """Load a single XYZ structure file"""
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
            
        natoms = int(lines[0].strip())
        data = []
        
        for line in lines[2:2+natoms]:  # Skip header lines
            parts = line.split()
            if len(parts) >= 4:
                atom, x, y, z = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
                data.append([atom, x, y, z])
        
        self.systems[system_id] = pd.DataFrame(
            data, 
            columns=['Atom type', 'x', 'y', 'z']
        )
    
    def _load_multi_xyz(self, xyz_file: Path, start_system_id: int) -> int:
        """Load multi-structure XYZ file (original behavior)"""
        natoms = []
        
        # First get number of atoms for each system
        with open(xyz_file, 'r') as f:
            for line in f:
                if line.strip().isdigit():
                    natoms.append(int(line.strip()))
        
        # Now extract coordinates
        start = 2
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
            system_id = start_system_id
            
            for atoms in natoms[::2]:  # Every other entry (skip duplicates)
                end = start + atoms
                data = []
                
                for line in lines[start:end]:
                    atom, x, y, z = line.split()
                    data.append([atom, float(x), float(y), float(z)])
                
                self.systems[system_id] = pd.DataFrame(
                    data, 
                    columns=['Atom type', 'x', 'y', 'z']
                )
                system_id += 1
                start = end + 2
                
        return system_id
    
    def _is_multi_structure_file(self, xyz_file: Path) -> bool:
        """Check if XYZ file contains multiple structures"""
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < 3:
            return False
            
        try:
            first_natoms = int(lines[0].strip())
            # Check if there's another atom count line after the expected position
            expected_next_structure_line = first_natoms + 2
            if expected_next_structure_line < len(lines):
                try:
                    int(lines[expected_next_structure_line].strip())
                    return True  # Found another atom count, it's multi-structure
                except ValueError:
                    return False  # Not a number, probably single structure
        except ValueError:
            return False
            
        return False
    
    def find_unique_atoms(self) -> Set[str]:
        """Get set of unique atoms across all systems"""
        for system in self.systems.values():
            self.unique_atoms.update(system['Atom type'].unique())
        return self.unique_atoms
    
    def calculate_all_pair_distances(self) -> Dict[Tuple[str, str], List[float]]:
        """Calculate distances for all possible atom type pairs"""
        atoms_list = list(self.unique_atoms)
        
        for i, atom1 in enumerate(atoms_list):
            for j in range(i, len(atoms_list)):  # Include self-pairs
                atom2 = atoms_list[j]
                distances = []
                
                # Loop through all systems
                for system in self.systems.values():
                    coords1 = system[system['Atom type'] == atom1][['x','y','z']].values
                    coords2 = system[system['Atom type'] == atom2][['x','y','z']].values
                    
                    # Calculate distances between all pairs
                    for c1 in coords1:
                        for c2 in coords2:
                            if not np.array_equal(c1, c2):  # Skip same atom
                                dist = np.linalg.norm(c1 - c2)
                                distances.append(dist)
                
                if distances:  # Store only if pairs were found
                    self.pair_distances[(atom1, atom2)] = distances
        
        return self.pair_distances
    
    def compute_kde(self, bw_method: float = 0.04) -> None:
        """Compute KDE for all atom pairs"""
        for key, distances in self.pair_distances.items():
            kde = stats.gaussian_kde(distances, bw_method=bw_method)
            x = np.linspace(min(distances), max(distances), 1000)
            y = kde(x)
            self.kde_data[key] = pd.DataFrame({'X': x, 'Y': y})

    
    def find_max_min(self, min_peak_height: float = 0.002) -> Tuple[Dict[Tuple[str, str], float], Dict[Tuple[str, str], float]]:
        """Find first significant maximum and following minimum in KDE for each atom pair"""
        for key in self.kde_data:
            # Find maxima
            max_index = signal.find_peaks_cwt(self.kde_data[key]['Y'], np.arange(1,10))
            if max_index.size != 0:
                for idx in max_index:
                    if self.kde_data[key].iloc[idx]['Y'] > min_peak_height:
                        x_max = self.kde_data[key].iloc[idx]['X']
                        self.X_max[key] = x_max
                        break
            
            # Find minima after maximum
            if key in self.X_max:
                min_index, _ = signal.find_peaks(-self.kde_data[key]['Y'])
                if min_index.size != 0:
                    for idx in min_index:
                        if self.kde_data[key].iloc[idx]['X'] > self.X_max[key]:
                            x_min = self.kde_data[key].iloc[idx]['X']
                            self.X_min[key] = x_min
                            break

        return self.X_max, self.X_min
    

    def determine_alpha(self, radii: Dict[str, float],
                       min_samples: int = 200,
                       max_alpha: float = 2.0) -> pd.DataFrame:
        """Calculate alpha values using stored maxima and minima"""
        alpha_values = {}

        for key in self.kde_data:
            if len(self.pair_distances[key]) < min_samples:
                print(f"Not enough statistics for {key}")
                continue

            if key not in self.X_min or key not in self.X_max:
                continue

            # Check for zero radii to avoid division by zero
            sum_radii = radii[key[0]] + radii[key[1]]
            if sum_radii == 0:
                print(f"Pair: {key}, Alpha: inf (zero radii sum)")
                continue
                
            alpha = self.X_min[key] / sum_radii
            
            if alpha <= max_alpha:
                if alpha <= (radii[key[0]] + radii[key[1]]):
                    alpha_values[key] = alpha
                else:
                    alpha_values[key] = np.nan
                    
            print(f"Pair: {key}, Alpha: {alpha:.3f}")

        # Create DataFrame and sort
        df = pd.DataFrame({'Pair': list(alpha_values.keys()), 
                          'Alpha': list(alpha_values.values())})
        df.sort_values(by='Alpha', inplace=True)
        
        return df.dropna()
    

    def determine_same_atom_radii(self, min_atom_count: int = 100, min_peak_height: float = 0.002, max_alpha: float = 1.2) -> Dict[str, List[float]]:
        """Determine atomic radii from same-atom pair distances"""

        radii = {}
        self.discarded_pairs = []
    
        for key in self.kde_data:
            if key[0] != key[1]:  # Skip different atom pairs
                continue
            
            atom_type = key[0]
            atom_count = sum(1 for sys in self.systems.values() 
                        for atom in sys['Atom type'] if atom == atom_type)
        
        # Check statistics
            if atom_count < min_atom_count:
                print(f"Not enough {atom_type} atoms ({atom_count})")
                radii[atom_type] = [0]
                self.discarded_pairs.append(key)
                continue
            
        # Find first maximum
            max_index, _ = signal.find_peaks(self.kde_data[key]['Y'])
            x_max = None
        
            if max_index.size != 0:
                for idx in max_index:
                    if self.kde_data[key].iloc[idx]['Y'] > min_peak_height:
                        x_max = self.kde_data[key].iloc[idx]['X']
                        break
                    
                if x_max is not None:
                    alpha = x_max / (2 * self.ASE_radii[atom_type])
                    if alpha <= max_alpha:
                        radii[atom_type] = [x_max / 2]
                        print(f"Radius for {atom_type}: {x_max/2:.3f} Å "
                              f"(ASE: {self.ASE_radii[atom_type]:.3f} Å)")
                    else:
                        print(f"Distance too large for {atom_type}")
                        radii[atom_type] = [0]
                        self.discarded_pairs.append(key)
                    
        return radii

    def determine_different_atom_radii(self, radii: Dict[str, List[float]],
                                 min_samples: int = 100,
                                 min_peak_height: float = 0.002,
                                 max_alpha: float = 1.1) -> Dict[str, List[float]]:
        """
    Determine remaining atomic radii using different-atom pairs
    
    Args:
        radii: Dictionary of already determined radii
        min_samples: Minimum number of distance samples
        min_peak_height: Minimum peak height for detection
        max_alpha: Maximum allowed ratio vs ASE radii
        
    Returns:
        Updated radii dictionary
        """
        for key in self.kde_data:
            if key[0] == key[1]:  # Skip same atom pairs
                continue
            
            atom1, atom2 = key
            
            # Ensure both atoms are in radii dict (add with 0 if missing)
            if atom1 not in radii:
                radii[atom1] = [0]
            if atom2 not in radii:
                radii[atom2] = [0]
                
            if ((radii[atom1][0] != 0 and radii[atom2][0] == 0) or
                (radii[atom1][0] == 0 and radii[atom2][0] != 0)):
            
            # Check statistics
                if len(self.pair_distances[key]) < min_samples:
                    print(f"Not enough samples for {key}")
                    self.discarded_pairs.append(key)
                    continue
                
            # Find first maximum
                max_index, _ = signal.find_peaks(self.kde_data[key]['Y'])
                if max_index.size == 0:
                    continue
                
                for idx in max_index:
                    if self.kde_data[key].iloc[idx]['Y'] > min_peak_height:
                        x_max = self.kde_data[key].iloc[idx]['X']
                        break
                else:
                    continue
                
            # Calculate new radius
                if radii[atom1][0] == 0:
                    alpha = x_max / (self.ASE_radii[atom1] + radii[atom2][0])
                    if alpha <= max_alpha:
                        new_radius = x_max - radii[atom2][0]
                        radii[atom1].append(new_radius)
                        print(f"Radius for {atom1}: {new_radius:.3f} Å")
                    else:
                        self.discarded_pairs.append(key)
                else:
                    alpha = x_max / (radii[atom1][0] + self.ASE_radii[atom2])
                    if alpha <= max_alpha:
                        new_radius = x_max - radii[atom1][0]
                        radii[atom2].append(new_radius)
                        print(f"Radius for {atom2}: {new_radius:.3f} Å")
                    else:
                        self.discarded_pairs.append(key)
                    
        return radii
    
    
    def save_results(self, output_dir: str):
        """Save analysis results"""
        results = AnalysisResults(
            timestamp=datetime.now().isoformat(),
            xyz_file=str(self.xyz_file),
            pair_distances=self.pair_distances,
            kde_data=self.kde_data,
            X_max=self.X_max,
            radii=self.radii,
            alpha_values=self.alpha_values
        )
        results.save(output_dir)

    def plot_pair_distributions(self, output_dir: str | None = None):
        """
        Plot pair distribution functions with maxima and minima markers
    
        Args:
            output_dir: Optional directory to save plots
        """
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
        plt.rcParams['font.size'] = 16
    
        bins = np.linspace(0, 10, 101)
    
        for key in self.kde_data.keys():
            if key in self.pair_distances:
            # Create figure
                fig, ax = plt.subplots()
            
            # Plot histogram and KDE
                sns.histplot(self.pair_distances[key], bins=bins, 
                        stat='density', alpha=0.4, element='step')
                ax.plot(self.kde_data[key]['X'], self.kde_data[key]['Y'])
            
            # Set plot style
                ax.set_title(f'{key[1]}-{key[0]}', fontsize=16)
                ax.set_xlabel(r'd / $\mathrm{\AA}$')
                ax.set_facecolor('None')
                ax.set_xlim(0, 10)
            
            # Add vertical lines for maxima and minima
                if key in self.X_max:
                    ax.axvline(x=self.X_max[key], linestyle='--', 
                          color='black', label='X max')
                
                if key in self.X_min:
                    ax.axvline(x=self.X_min[key], linestyle='--',
                          color='salmon', label='X min')
            
                ax.legend()
            
            # Save plot if output directory specified
                if output_dir:
                    out_path = Path(output_dir)
                    out_path.mkdir(parents=True, exist_ok=True)
                    plt.savefig(out_path / f'pdf_{key[0]}_{key[1]}.pdf',
                           bbox_inches='tight', dpi=300)
            
                plt.close()

    
    def analyze(self) -> Dict[Tuple[str, str], List[float]]:
        """Run complete analysis workflow"""
        print("Loading systems...")
        self.load_systems()
        
        print("Finding unique atoms...")
        self.find_unique_atoms()
        
        print("Calculating pair distances...")
        self.calculate_all_pair_distances()
        
        print("Computing KDE...")
        self.compute_kde()
        
        print("Finding maximum and minum...")
        self.find_max_min()
        
        return self.pair_distances