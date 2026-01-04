import marimo

__generated_with = "0.18.4"
app = marimo.App()


@app.cell
def _():
    import pandas as pd
    return (pd,)


@app.cell
def _(pd):
    average_precision = pd.read_csv("../01_map/output/harmony_activity_map_results.csv")
    average_precision
    return (average_precision,)


@app.cell
def _(pd):
    cell_count_info = pd.read_parquet("../data/profiles_var_mad_int_featselect_harmony.parquet")[['Metadata_JCP2022', "X_623"]]
    cell_count_info = cell_count_info.groupby('Metadata_JCP2022').median().reset_index()
    cell_count_info
    return (cell_count_info,)


@app.cell
def _():
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


@app.cell
def _(average_precision, cell_count_info, pd):
    df = pd.merge(average_precision, cell_count_info, on=['Metadata_JCP2022'])
    df
    return (df,)


@app.cell
def _():
    #Analysis of phenotypic activity rates (per compound)
    return


@app.cell
def _():
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scipy.stats as stats
    return plt, sns, stats


@app.cell
def _(df, plt, sns):
    # Plot correlation
    plt.figure(figsize=(6, 4))
    sns.scatterplot(x=df['mean_average_precision'], y=df['X_623'], s=1, alpha=0.5)
    plt.xlabel('mean_average_precision')
    plt.ylabel('X_623')
    plt.show()
    return


@app.cell
def _(mo):
    mo.md(r"""
    Mean Average Precision appears to be higher when the cell count is lower
    """)
    return


@app.cell
def _(df, pd, plt, sns):
    # Define custom bins for mean_average_precision
    bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    labels = [f'{bins[i]}-{bins[i + 1]}' for i in range(len(bins) - 1)]
    df['map_bin'] = pd.cut(df['mean_average_precision'], bins=bins, labels=labels, include_lowest=True)
    plt.figure(figsize=(10, 5))
    # Create figure for transparent boxplot with scatter points behind
    sns.boxplot(x=df['map_bin'], y=df['X_623'], showfliers=False, boxprops={'alpha': 0.5})
    plt.xticks(rotation=90)
    # Scatter plot (background points)
    # sns.stripplot(x=df["map_bin"], y=df["X_623"], color="black", alpha=1, size=0.4, jitter=True)
    plt.xlabel('Mean Average Precision (Binned)')
    # Transparent boxplot
    plt.ylabel('X_623')
    plt.title('Transparent Boxplot with Scatter of X_623 by Mean Average Precision Bins')
    # Rotate x-axis labels for better readability
    # Labels and title
    plt.show()
    return


@app.cell
def _():
    #Analysis of phenotypic activity rates (per compound) and cell count
    return


@app.cell
def _(df, stats):
    (_correlation, _p_val) = stats.pearsonr(df['X_623'], df['mean_average_precision'])
    print(_correlation, _p_val)
    return


@app.cell
def _(df, plt, sns):
    import ptitprince as pt
    plt.figure(figsize=(6, 4))
    # Create figure
    ax = pt.RainCloud(x='below_corrected_p', y='X_623', data=df, palette='Set2', bw=0.2, width_viol=0.6, orient='v', alpha=0.65)
    sns.boxplot(x='below_corrected_p', y='X_623', data=df, showcaps=True, boxprops={'facecolor': 'None'}, showfliers=False, whiskerprops={'linewidth': 2})
    # Create raincloud plot
    plt.xlabel('Below Corrected p-value')
    plt.ylabel('X_623')
    plt.title('Raincloud Plot of X_623 by Below Corrected p-value')
    # Add statistical annotation
    # Add labels
    plt.show()
    return


@app.cell
def _(df, stats):
    (_correlation, _p_val) = stats.pearsonr(df['below_corrected_p'], df['X_623'])
    print(_correlation, _p_val)
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


@app.cell
def _(mo):
    mo.md(r"""
    Negative correlation between cell count and below_corrected_p
    """)
    return


@app.cell
def _(df):
    # Compare cell count between below_corrected_p groups
    grouped_counts = df.groupby('below_corrected_p')['X_623'].mean()
    grouped_counts
    return


@app.cell
def _():
    #Analysis of phenotypic activity rates at well level. Plot average mAP based on well position
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
