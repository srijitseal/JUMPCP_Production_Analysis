import marimo

__generated_with = "0.18.4"
app = marimo.App()


@app.cell
def _():
    import pandas as pd
    return (pd,)


@app.cell
def _(pd):
    cell_count_info = pd.read_parquet("../data/profiles_var_mad_int_featselect_harmony.parquet")
    cell_count_info = cell_count_info[['Metadata_Source', 'Metadata_JCP2022', 'Metadata_Plate', 'Metadata_Well', 'X_623']]
    cell_count_info
    return (cell_count_info,)


@app.cell
def _(pd):
    average_precision = pd.read_csv("../01_map/output/harmony_activity_ap_scores.csv")
    #average_precision = pd.read_csv("output/phenotypic-activity-compoundslevel-profiles_wellpos_cc_var_mad_outlier_featselect_sphering_harmony.csv.gz")
    average_precision
    return (average_precision,)


@app.cell
def _(average_precision, cell_count_info, pd):
    df_1 = pd.merge(average_precision, cell_count_info, on=['Metadata_JCP2022', 'Metadata_Well', 'Metadata_Plate'])
    return (df_1,)


@app.cell
def _(df_1):
    df_2 = df_1[~df_1['average_precision'].isna()]
    df_2
    return (df_2,)


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
def _(df_2, plt, sns):
    # Plot correlation
    plt.figure(figsize=(6, 4))
    sns.scatterplot(x=df_2['average_precision'], y=df_2['X_623'], s=1, alpha=0.5)
    plt.xlabel('average_precision')
    plt.ylabel('X_623')
    plt.show()
    return


@app.cell
def _(df_2, plt, sns):
    # Plot correlation
    plt.figure(figsize=(12, 4))
    sns.scatterplot(x=df_2['Metadata_Source_x'], y=df_2['X_623'], s=1, alpha=0.5)
    plt.xlabel('Metadata_Source_x')
    plt.ylabel('X_623')
    plt.show()
    return


@app.cell
def _(df_2, plt, sns):
    # Plot correlation
    plt.figure(figsize=(12, 4))
    sns.boxplot(x=df_2['Metadata_Source_x'], y=df_2['average_precision'])
    plt.xlabel('Metadata_Source_x')
    plt.ylabel('average_precision')
    plt.show()
    return


@app.cell
def _():
    return


@app.cell
def _(df_2, pd, plt, sns):
    # Define custom bins for mean_average_precision
    bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    labels = [f'{bins[i]}-{bins[i + 1]}' for i in range(len(bins) - 1)]
    df_2['map_bin'] = pd.cut(df_2['average_precision'], bins=bins, labels=labels, include_lowest=True)
    plt.figure(figsize=(10, 5))
    # Create figure for transparent boxplot with scatter points behind
    sns.boxplot(x=df_2['map_bin'], y=df_2['X_623'], showfliers=False, boxprops={'alpha': 0.5})
    plt.xticks(rotation=90)
    # Scatter plot (background points)
    # sns.stripplot(x=df["map_bin"], y=df["X_623"], color="black", alpha=1, size=0.4, jitter=True)
    plt.xlabel('Average Precision (Binned)')
    # Transparent boxplot
    plt.ylabel('X_623')
    plt.title('Transparent Boxplot with Scatter of X_623 by Mean Average Precision Bins')
    # Rotate x-axis labels for better readability
    # Labels and title
    plt.show()
    return


@app.cell
def _(mo):
    mo.md(r"""
    Higher AP is associated with lower cell count feature values
    """)
    return


@app.cell
def _():
    #Analysis of phenotypic activity rates (per compound) and cell count
    return


@app.cell
def _(df_2, stats):
    (_correlation, _p_val) = stats.pearsonr(df_2['X_623'], df_2['average_precision'])
    print(_correlation, _p_val)
    return


@app.cell
def _(df_2):
    df_2['high_ap'] = df_2['average_precision'] > 0.3
    df_2
    return


@app.cell
def _(df_2, plt, sns):
    import ptitprince as pt
    plt.figure(figsize=(6, 4))
    # Create figure
    ax = pt.RainCloud(x='high_ap', y='X_623', data=df_2, palette='Set2', bw=0.2, width_viol=0.6, orient='v', alpha=0.65)
    sns.boxplot(x='high_ap', y='X_623', data=df_2, showcaps=True, boxprops={'facecolor': 'None'}, showfliers=False, whiskerprops={'linewidth': 2})
    # Create raincloud plot
    plt.xlabel('high_ap')
    plt.ylabel('X_623')
    plt.title('Raincloud Plot of X_623 by high_ap')
    # Add statistical annotation
    # Add labels
    plt.show()
    return


@app.cell
def _(df_2, stats):
    (_correlation, _p_val) = stats.pearsonr(df_2['high_ap'], df_2['X_623'])
    print(_correlation, _p_val)
    return


@app.cell
def _(df_2):
    # Compare cell count between below_corrected_p groups
    grouped_counts = df_2.groupby('high_ap')['X_623'].mean()
    grouped_counts
    return


@app.cell
def _():
    #Analysis of phenotypic activity rates at well level. Plot average mAP based on well position
    return


@app.cell
def _(df_2):
    df_2["Metadata_Column"] = df_2["Metadata_Well"].str[:-2]
    df_2['Metadata_Row'] = df_2["Metadata_Well"].str[-2:].astype(int)
    df_3 = df_2[df_2["Metadata_Column"].isin(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M','N', 'O', 'P'])]
    df_3 = df_3[df_3["Metadata_Row"]<=16]
    return (df_3,)


@app.cell
def _(df_3):
    df_3.Metadata_Column.unique()
    return


@app.cell
def _(df_3, plt, sns):
    # Plot correlation
    plt.figure(figsize=(12, 4))
    sns.boxplot(x=df_3['Metadata_Column'], y=df_3['average_precision'])
    plt.xlabel('Metadata_Column')
    plt.ylabel('average_precision')
    plt.show()
    return


@app.cell
def _(df_3, plt, sns):
    # Plot correlation
    plt.figure(figsize=(30, 4))
    sns.boxplot(x=df_3['Metadata_Row'], y=df_3['average_precision'])
    plt.xlabel('Metadata_Row')
    plt.ylabel('average_precision')
    plt.show()
    return


@app.cell
def _(mo):
    mo.md(r"""
    Row 1 appears to be DMSO with highest MAP. No other trends in 364 well plates overall.
    """)
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
