import marimo

__generated_with = "0.18.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import pandas as pd
    return (pd,)


@app.cell
def _(pd):
    average_precision = pd.read_csv("../01_map/output/harmony_activity_ap_scores.csv")
    average_precision
    return (average_precision,)


@app.cell
def _(pd):
    cell_count_info = pd.read_parquet("../data/profiles_var_mad_int_featselect_harmony.parquet")
    cell_count_info
    return (cell_count_info,)


@app.cell
def _(average_precision, cell_count_info, pd):
    df = pd.merge(cell_count_info [['Metadata_Source', 'Metadata_Plate', 'Metadata_Well','Metadata_JCP2022', "X_623"]], average_precision)

    #Remove 1600 well plates
    df = df[~df["Metadata_Source"].isin(["source_1", "source_9"])].reset_index(drop=True)

    df = df[~df["average_precision"].isna()].reset_index(drop=True)
    df
    return (df,)


@app.cell
def _(df):
    df.Metadata_Source.value_counts()
    return


@app.cell
def _(df):
    df.Metadata_pert_type.value_counts()
    return


@app.cell
def _(df):
    # Define edge wells
    edge_rows = {"A",  "P"} #Maybe add "B", "O"?
    #edge_columns = {"1", "24"} #Remove these, they are controls! #####REMOVE
    edge_columns = {}

    # Function to check if a well is on the edge
    def is_edge_well(well):
        row, col = well[0], well[1:]
        return row in edge_rows or col in edge_columns

    # Apply function
    df["Edge_Well"] = df["Metadata_Well"].apply(is_edge_well)
    df
    return


@app.cell
def _():
    return


@app.cell
def _(df):
    df_active = df[df["average_precision"]>0.4]
    df_active
    return


@app.cell
def _(df):
    import ptitprince as pt
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Create figure
    plt.figure(figsize=(6, 4))

    # Create raincloud plot
    ax = pt.RainCloud(x="Edge_Well", y="average_precision", data=df, 
                      bw=0.2, width_viol=0.6, orient="v", alpha=0.65)

    # Add statistical annotation
    sns.boxplot(x="Edge_Well", y="average_precision", data=df, 
                showcaps=True, boxprops={'facecolor':'None'}, showfliers=False, 
                whiskerprops={'linewidth':2})

    # Add labels
    plt.xlabel("Edge_Well")
    plt.ylabel("average_precision")
    plt.title("Raincloud Plot of average_precision by Edge_Well")

    plt.show()
    return plt, pt, sns


@app.cell
def _(df):
    import scipy.stats as stats

    correlation, p_val = stats.pearsonr(df["Edge_Well"], df["average_precision"])
    print(correlation, p_val)
    return (stats,)


@app.cell
def _(df):
    # Compare cell count between below_corrected_p groups
    grouped_counts = df.groupby("Edge_Well")["average_precision"].mean()
    grouped_counts
    return


@app.cell
def _():
    #Analysis of phenotypic activity rates at well level. Plot average mAP based on well position, edge or not
    return


@app.cell
def _(df, plt, pt, sns):
    # Create figure
    plt.figure(figsize=(6, 4))

    # Create raincloud plot
    ax2 = pt.RainCloud(x="Edge_Well", y="X_623", data=df, 
                      bw=0.2, width_viol=0.6, orient="v", alpha=0.65)

    # Add statistical annotation
    sns.boxplot(x="Edge_Well", y="X_623", data=df, 
                showcaps=True, boxprops={'facecolor':'None'}, showfliers=False, 
                whiskerprops={'linewidth':2})

    # Add labels
    plt.xlabel("Edge_Well")
    plt.ylabel("X_623")
    plt.title("Raincloud Plot of X_623 by Edge_Well")

    plt.show()
    return


@app.cell
def _(df, stats):
    stats.pearsonr(df["Edge_Well"], df["X_623"])
    return


@app.cell
def _(df):
    # Compare cell count between below_corrected_p groups
    df.groupby("Edge_Well")["X_623"].mean()
    return


@app.cell
def _(mo):
    mo.md(r"""
    There no significant difference in cell count or AP based on whether the compound is in the edge well or not.
    """)
    return


@app.cell
def _():
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
