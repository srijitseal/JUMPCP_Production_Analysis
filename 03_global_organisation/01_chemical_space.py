import marimo

__generated_with = "0.18.4"
app = marimo.App(width="full")


@app.cell
def _():
    import pandas as pd
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from scipy.spatial.distance import pdist
    from tqdm import tqdm
    return (
        Chem,
        Descriptors,
        PCA,
        StandardScaler,
        np,
        pd,
        pdist,
        plt,
        sns,
        tqdm,
    )


@app.cell
def _(pd):
    # Load compound SMILES data
    df = pd.read_csv("../data/compound_smilesr.csv.gz")
    print(f"Loaded {len(df)} compounds")
    df.head()
    return (df,)


@app.cell
def _(Chem, Descriptors, df, pd, tqdm):
    """Calculate physicochemical properties for all compounds"""

    def calculate_properties(smiles):
        """Calculate physicochemical descriptors from SMILES"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        return {
            'MW': Descriptors.MolWt(mol),
            'LogP': Descriptors.MolLogP(mol),
            'HBA': Descriptors.NumHAcceptors(mol),
            'HBD': Descriptors.NumHDonors(mol),
            'TPSA': Descriptors.TPSA(mol),
            'RotBonds': Descriptors.NumRotatableBonds(mol)
        }

    # Calculate properties
    print("Calculating physicochemical properties...")
    properties_list = []
    valid_indices = []

    for idx, smiles in tqdm(enumerate(df['smiles_r']), total=len(df)):
        props = calculate_properties(smiles)
        if props is not None:
            properties_list.append(props)
            valid_indices.append(idx)

    properties_df = pd.DataFrame(properties_list)
    df_valid = df.iloc[valid_indices].reset_index(drop=True)
    df_with_props = pd.concat([df_valid, properties_df], axis=1)

    print(f"Calculated properties for {len(properties_df)} valid compounds")
    return (properties_df,)


@app.cell
def _(PCA, StandardScaler, properties_df):
    """PCA on physicochemical properties"""

    print("Performing PCA on physicochemical properties...")

    # Standardize the properties
    scaler = StandardScaler()
    properties_scaled = scaler.fit_transform(properties_df)

    # Perform PCA
    pca = PCA(n_components=2, random_state=42)
    pca_coords = pca.fit_transform(properties_scaled)

    # Get variance explained
    var_explained = pca.explained_variance_ratio_ * 100

    print(f"PC1 variance explained: {var_explained[0]:.2f}%")
    print(f"PC2 variance explained: {var_explained[1]:.2f}%")
    print(f"Total variance explained: {var_explained.sum():.2f}%")
    return pca_coords, var_explained


@app.cell
def _(pca_coords, plt, var_explained):
    """PCA visualization"""

    fig, ax = plt.subplots(figsize=(10, 8))

    # Scatter plot
    ax.scatter(pca_coords[:, 0], pca_coords[:, 1], alpha=0.5, s=1)

    # Add labels with variance explained
    ax.set_xlabel(f'PC1 ({var_explained[0]:.1f}% variance)', fontsize=12)
    ax.set_ylabel(f'PC2 ({var_explained[1]:.1f}% variance)', fontsize=12)
    ax.set_title('Chemical Space Visualization (PCA on Physicochemical Properties)', fontsize=14)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig
    return


@app.cell
def _(np, pdist, properties_df):
    """Calculate pairwise Euclidean distances in property space"""

    print("Calculating pairwise distances in property space...")
    # Sample for computational efficiency
    n_sample = min(5000, len(properties_df))
    indices = np.random.choice(len(properties_df), n_sample, replace=False)
    props_sample = properties_df.iloc[indices].values

    # Calculate Euclidean distances
    distances = pdist(props_sample, metric='euclidean')

    print(f"Calculated {len(distances)} pairwise distances")
    return (distances,)


@app.cell
def _(distances, np, plt):
    """Distance distribution histogram"""

    fig_dist, ax_dist = plt.subplots(figsize=(10, 5))

    # Histogram
    ax_dist.hist(distances, bins=50, alpha=0.7, edgecolor='black')

    # Add statistics
    mean_dist = distances.mean()
    median_dist = np.median(distances)

    ax_dist.axvline(mean_dist, color='red', linestyle='--', linewidth=2,
                    label=f'Mean: {mean_dist:.2f}')
    ax_dist.axvline(median_dist, color='blue', linestyle='--', linewidth=2,
                    label=f'Median: {median_dist:.2f}')

    ax_dist.set_xlabel('Euclidean Distance', fontsize=12)
    ax_dist.set_ylabel('Count', fontsize=12)
    ax_dist.set_title('Pairwise Euclidean Distance Distribution', fontsize=14)
    ax_dist.legend()
    ax_dist.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_dist
    return


@app.cell
def _(plt, properties_df, sns):
    """Molecular property distributions - pairplot"""

    # Sample for visualization efficiency
    n_viz = min(10000, len(properties_df))
    viz_sample = properties_df.sample(n=n_viz, random_state=42)

    # Create pairplot
    fig_props = sns.pairplot(
        viz_sample,
        vars=['MW', 'LogP', 'HBA', 'HBD', 'TPSA', 'RotBonds'],
        diag_kind='hist',
        plot_kws={'alpha': 0.3, 's': 5},
        diag_kws={'bins': 30, 'edgecolor': 'black'}
    )

    fig_props.fig.suptitle('Molecular Property Distributions', y=1.02, fontsize=16)
    plt.tight_layout()
    fig_props
    return


@app.cell
def _(properties_df):
    """Summary statistics"""

    summary_stats = properties_df.describe()
    summary_stats
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
