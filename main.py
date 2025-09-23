from PDF import PDFAnalyzer, AnalysisResults
from datetime import datetime
from pathlib import Path

def main():
    # Get the directory containing main.py
    script_dir = Path(__file__).parent.absolute()
    
    # Create results directory if it doesn't exist
    output_dir = script_dir / "results"
    output_dir.mkdir(exist_ok=True)
    
    # Find all XYZ files in the data directory
    data_dir = script_dir / "data"
    xyz_files = list(data_dir.glob("*.xyz"))
    
    if not xyz_files:
        raise FileNotFoundError(f"No .xyz files found in: {data_dir}")
    
    print(f"Found {len(xyz_files)} .xyz files in {data_dir}")
    for xyz_file in xyz_files:
        print(f"  - {xyz_file.name}")
    
    # Initialize analyzer with all XYZ files
    analyzer = PDFAnalyzer(xyz_files)
    
    # Run analysis
    pair_distances = analyzer.analyze()
    
    # Get all results
    radii = analyzer.determine_same_atom_radii()
    radii = analyzer.determine_different_atom_radii(radii)
    X_max, X_min = analyzer.find_max_min(min_peak_height=0.002)
    
    # Clean up radii - use ASE radii as fallback for zero values
    for atom, radius_list in radii.items():
        radius = next((r for r in radius_list if r != 0), 0)
        if radius == 0 and atom in analyzer.ASE_radii:
            radius = analyzer.ASE_radii[atom]
            print(f"Using ASE fallback radius for {atom}: {radius:.3f} Å")
        radii[atom] = radius
    
    # Get alpha values
    df_alpha = analyzer.determine_alpha(radii)
    
    # Create summary data
    pair_summary = {}
    for (atom1, atom2), distances in pair_distances.items():
        pair_summary[(atom1, atom2)] = {
            'count': len(distances),
            'min_distance': min(distances),
            'max_distance': max(distances)
        }
    
    # Create results object
    results = AnalysisResults(
        timestamp=datetime.now().isoformat(),
        xyz_file=xyz_files,  # Now a list of files
        pair_distances=pair_distances,
        kde_data=analyzer.kde_data,
        X_max=X_max,
        X_min=X_min,
        radii=radii,
        alpha_values=df_alpha,
        pair_summary=pair_summary
    )
    
    # Save results
    results.save(output_dir)
    
    # Print summary
    print("\nAnalysis Results Summary")
    print("=" * 40)
    print(f"\nPair Distance Statistics:")
    for (atom1, atom2), stats in pair_summary.items():
        print(f"\n{atom1}-{atom2} pairs:")
        print(f"  Count: {stats['count']}")
        print(f"  Range: {stats['min_distance']:.3f} - {stats['max_distance']:.3f} Å")
    
    print("\nAtomic Radii:")
    for atom, radius in radii.items():
        print(f"  {atom}: {radius:.3f} Å")
    
    print("\nAlpha Values:")
    print(df_alpha)
    
    print(f"\nResults saved to: {output_dir}")
     # Create plots subdirectory
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    analyzer.plot_pair_distributions(plots_dir)

if __name__ == "__main__":
    main()