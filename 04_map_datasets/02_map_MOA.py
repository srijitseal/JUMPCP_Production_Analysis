import marimo

__generated_with = "0.18.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import pandas as pd
    import utils
    import numpy as np
    import warnings
    from tqdm import tqdm
    from copairs import map
    from copairs.map import average_precision, mean_average_precision
    return mean_average_precision, np, pd, utils


@app.cell
def _(np, pd):
    filename = "../data/moa/moa-rf_curated"

    moa_annotations_df = pd.read_csv(
        f"{filename}.csv"
    ).set_index("Metadata_JCP2022")

    moa_annotations_df = (
        (
            moa_annotations_df.stack()
            .reset_index()
            .rename(
                columns={0: "connection", "level_1": "gene", "JUMP_ID": "Metadata_JCP2022"}
            )
            .query("connection != 0")
            .drop(columns=["connection"])
        )
        .groupby("Metadata_JCP2022")
        .apply(lambda x: "|".join(np.unique(x["gene"])))
        .reset_index()
        .rename(columns={0: "Metadata_moa_list"})
    )

    moa_annotations_df
    return (moa_annotations_df,)


@app.cell
def _():
    return


@app.cell
def _(moa_df):
    moa_df
    return


@app.cell
def _(pd):
    compound_phenotypic_activity_df = pd.read_csv("../01_mAP/output/harmony_activity_map_results.csv")
    moa_df = pd.read_csv("../data/moa/moa-rf_curated.csv")

    print(f"Number of total compounds in this dataset: {len(compound_phenotypic_activity_df)}")

    print(f"Number of total compounds with a target in this dataset: {len(moa_df)}")

    print(f"Number of total compounds with phenotypic activity in this dataset: {len(compound_phenotypic_activity_df.query('below_corrected_p==True'))}")

    moa_df = moa_df.set_index("Metadata_JCP2022")

    print(moa_df.shape)
    return (moa_df,)


@app.cell
def _(moa_df):
    moa_df
    return


@app.cell
def _():
    return


@app.cell
def _():
    operations = "var_mad_int_featselect_harmony"
    batch_size = 20000
    null_size = 20000
    fdr = 0.05

    multi_label_col = "Metadata_target_list"

    pos_sameby = [f"{multi_label_col}"]
    pos_diffby = []
    neg_sameby = []
    neg_diffby = [f"{multi_label_col}"]
    return (
        batch_size,
        fdr,
        multi_label_col,
        neg_diffby,
        neg_sameby,
        null_size,
        pos_diffby,
        pos_sameby,
    )


@app.cell
def _(moa_annotations_df, multi_label_col, pd, utils):
    compound_df = pd.read_parquet("../data/profiles_var_mad_int_featselect_harmony.parquet")

    print(compound_df.shape)

    compound_df = compound_df[compound_df["Metadata_JCP2022"] != 'JCP_2022_033924'].reset_index(drop=True)

    print(compound_df.shape)

    compound_df = utils.remove_nan_features(compound_df)

    print(compound_df.shape)

    compound_df = compound_df.merge(
        moa_annotations_df, on="Metadata_JCP2022", how="inner"
    )

    ### Expand target list

    compound_df = (
        compound_df.assign(col = lambda x: x["Metadata_moa_list"].str.split("|"))
        .rename(columns={"col": f"{multi_label_col}"})
    )

    print(compound_df.shape)
    compound_df
    return (compound_df,)


@app.cell
def _(compound_df):
    compound_df.columns
    return


@app.cell
def _(compound_df, moa_annotations_df, multi_label_col):
    median_compound_df = compound_df.drop(['Metadata_Source',
           'Metadata_Plate', 'Metadata_Well', "Metadata_moa_list", "Metadata_target_list"], axis=1).groupby('Metadata_JCP2022').median().reset_index()
    median_compound_df

    median_compound_df = median_compound_df.merge(
        moa_annotations_df, on="Metadata_JCP2022", how="inner"
    )

    ### Expand target list
    median_compound_df = (
        median_compound_df.assign(col = lambda x: x["Metadata_moa_list"].str.split("|"))
        .rename(columns={"col": f"{multi_label_col}"})
    )
    median_compound_df
    return (median_compound_df,)


@app.cell
def _(compound_df):
    compound_df.columns[-15:]
    return


@app.cell
def _(pd):
    phenotypic_activity_df = pd.read_csv(
        f"../01_mAP/output/harmony_activity_map_results.csv",
        usecols=["Metadata_JCP2022", "below_corrected_p"],
    ).query("below_corrected_p==True")
    phenotypic_activity_df
    return (phenotypic_activity_df,)


@app.cell
def _(median_compound_df, phenotypic_activity_df):
    median_compound_df_phenotactive = median_compound_df.merge(
        phenotypic_activity_df, on="Metadata_JCP2022", how="inner"
    ).drop(columns=["below_corrected_p"])
    print(median_compound_df_phenotactive.shape)
    return (median_compound_df_phenotactive,)


@app.cell
def _(median_compound_df_phenotactive, utils):
    consensus_df = utils.consensus(median_compound_df_phenotactive, "Metadata_JCP2022")
    consensus_df
    return (consensus_df,)


@app.cell
def _(
    batch_size,
    consensus_df,
    fdr,
    mean_average_precision,
    multi_label_col,
    neg_diffby,
    neg_sameby,
    null_size,
    pos_diffby,
    pos_sameby,
    utils,
):
    from copairs.map import multilabel

    metadata_df = utils.get_metadata(consensus_df)
    feature_df = utils.get_featuredata(consensus_df)
    feature_values = feature_df.values

    result = multilabel.average_precision(
            metadata_df,
            feature_values,
            pos_sameby,
            pos_diffby,
            neg_sameby,
            neg_diffby,
            batch_size=batch_size,
            multilabel_col=multi_label_col,
        )

    agg_result = mean_average_precision(
            result, pos_sameby, null_size, threshold=fdr, seed=12527
        )

    agg_result.to_csv(f"output/phenotypic-consistency-target-retrieval.csv.gz", index=False)
    return (agg_result,)


@app.cell
def _(agg_result):
    df = agg_result
    df["plate_type"] = "COMPOUND"

    n_retrieved = (
        df.groupby("plate_type")
        .below_corrected_p.sum()
        .reset_index()
        .rename(columns={"below_corrected_p": "n_retrieved"})
    )

    compounds_n_retrieved = n_retrieved[
        n_retrieved.plate_type == "COMPOUND"
    ].n_retrieved.values[0]
    return df, n_retrieved


@app.cell
def _(n_retrieved):
    n_retrieved
    return


@app.cell
def _(df):
    df.sort_values("mean_average_precision", ascending=False)
    return


@app.cell
def _(df, np):
    import matplotlib.pyplot as plt

    # Thresholds
    effect_size_threshold = 0.1  # Threshold on mean average precision
    pval_threshold = -np.log10(0.05)  # Threshold on corrected p-value (corresponding to p = 0.05)

    # Get data
    effect_size = df['mean_average_precision']
    neg_log_pvalue = -np.log10(df['corrected_p_value'])  # Transform p-values

    # Create figure
    plt.figure(figsize=(8, 6))

    # Plot all points
    plt.scatter(effect_size, neg_log_pvalue, color='grey', alpha=0.6, edgecolors='w', s=40)

    # Identify significant points
    significant = (effect_size > effect_size_threshold) & (neg_log_pvalue > pval_threshold)

    # Plot significant points
    plt.scatter(effect_size[significant], neg_log_pvalue[significant], color='red', label='Significant', s=50)

    # Threshold lines
    plt.axhline(y=pval_threshold, color='black', linestyle='--', linewidth=1)
    plt.axvline(x=effect_size_threshold, color='black', linestyle='--', linewidth=1)

    # Labels and title
    plt.xlabel('Mean Average Precision (Effect Size)', fontsize=12)
    plt.ylabel('-Log$_{10}$ Corrected p-value', fontsize=12)
    plt.title('Volcano Plot (Positive Effect Sizes Only)', fontsize=14)
    plt.legend(frameon=True)
    plt.grid(True, linestyle='--', alpha=0.5)

    # Aesthetic improvements
    plt.tight_layout()
    plt.show()
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
