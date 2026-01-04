import marimo

__generated_with = "0.18.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import pandas as pd
    import utilscp
    import numpy as np
    import warnings
    from tqdm import tqdm
    from copairs import map
    from copairs.map import average_precision, mean_average_precision

    operations = "var_mad_int_featselect_harmony"
    batch_size = 20000
    null_size = 20000
    fdr = 0.1


    jump_df = pd.read_parquet(
        "../data/profiles_var_mad_int_featselect_harmony.parquet"
    )

    jump_df= jump_df[["Metadata_Source",
    "Metadata_Plate",
    "Metadata_Well",
    "Metadata_JCP2022",
    "X_623"
    ]]

    print(jump_df.shape)

    compoundinfo = pd.read_csv("../data/compound.csv.gz")
    plateinfo = pd.read_csv("../data/plate.csv.gz")

    jump_df = utilscp.remove_nan_features(jump_df)
    jump_df = jump_df[jump_df.Metadata_JCP2022.isin(compoundinfo.Metadata_JCP2022.to_list())]
    print(jump_df.shape)
    return (
        average_precision,
        batch_size,
        compoundinfo,
        fdr,
        jump_df,
        mean_average_precision,
        np,
        null_size,
        plateinfo,
        utilscp,
    )


@app.cell
def _(jump_df):
    jump_df
    return


@app.cell
def _():
    import altair as alt
    return


@app.cell
def _(compoundinfo):
    compoundinfo
    return


@app.cell
def _(plateinfo):
    plateinfo
    return


@app.cell
def _(
    average_precision,
    batch_size,
    fdr,
    jump_df,
    mean_average_precision,
    np,
    null_size,
    utilscp,
):


    # Adding a new column for negative control
    jump_df["Metadata_negcon"] = np.where(jump_df["Metadata_JCP2022"] == "JCP2022_033924", 1, 0) #DMSO
    pos_sameby = ["Metadata_JCP2022"]
    pos_diffby = []
    neg_sameby = ["Metadata_Plate"]
    neg_diffby = ["Metadata_negcon"]
    metadata_df = utilscp.get_metadata(jump_df)
    feature_df = utilscp.get_featuredata(jump_df)
    feature_values = feature_df.values

    result = average_precision(
        metadata_df, feature_values, pos_sameby, pos_diffby, neg_sameby, neg_diffby, batch_size=batch_size
    )
    # Remove negcon
    result = result.query('Metadata_pert_type!="negcon"').reset_index(drop=True)

    agg_result = (
        mean_average_precision(result, pos_sameby, null_size=null_size, threshold=fdr, seed=12527)
        .rename(columns={'average_precision': 'mean_average_precision'})
    )
    agg_result.to_csv(f"output/phenotypic-activity-profiles_var_mad_int_featselect_sphering_harmony.csv.gz", index=False)

    result.to_csv(f"output/phenotypic-activity-compoundslevel-profiles_var_mad_int_featselect_sphering_harmony.csv.gz", index=False)
    return


if __name__ == "__main__":
    app.run()
