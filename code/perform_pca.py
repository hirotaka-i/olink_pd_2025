import os
import pandas as pd
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

def perform_pca(
    df,
    id_name,
    prot_column_num,
    missingness_threshold=0.05,
    n_components=10,
    outlier_sd_threshold=4,
    output_scores="temp/olink_pcs.csv",
    plot_prefix="temp/olink_pca"
):
    # Ensure output directories exist
    os.makedirs(os.path.dirname(output_scores), exist_ok=True)
    os.makedirs(os.path.dirname(plot_prefix), exist_ok=True)

    # 1) Check ID uniqueness
    if not df[id_name].is_unique:
        raise ValueError(f"ID column '{id_name}' is not unique.")

    # 2) Subset to proteomic columns and set index
    Xp = df.iloc[:, prot_column_num:].set_index(df[id_name])
    print("Input data shape:", Xp.shape)

    # 3) Filter out high-missingness probes (columns)
    col_missing = Xp.isna().mean()
    Xp = Xp.loc[:, col_missing < missingness_threshold]
    print(f"After column filtering (missingness < {missingness_threshold}):", Xp.shape)

    # 4) Filter out high-missingness samples (rows)
    row_missing = Xp.isna().mean(axis=1)
    Xp = Xp.loc[row_missing <= missingness_threshold]
    print(f"After row filtering (missingness <= {missingness_threshold}):", Xp.shape)

    # 5) Iterative PCA with outlier removal
    iteration = 0
    while True:
        iteration += 1
        print(f"PCA iteration {iteration}")

        # 5a) Impute missing values
        imputer = SimpleImputer(strategy="mean")
        X_imp = imputer.fit_transform(Xp)

        # 5b) Standardize data
        scaler = StandardScaler()
        X_imp = scaler.fit_transform(X_imp)

        # 5c) Perform PCA
        max_comp = min(n_components, X_imp.shape[1])
        pca = PCA(n_components=max_comp)
        scores = pca.fit_transform(X_imp)

        # 5d) Create DataFrame of PC scores
        pcs = pd.DataFrame(
            scores,
            index=Xp.index,
            columns=[f"pPC{i+1}" for i in range(max_comp)]
        )

        # 5e) Identify outliers on PC1 and PC2
        outlier_mask = pd.Series(False, index=pcs.index)
        for pc in ["pPC1", "pPC2"]:
            mean = pcs[pc].mean()
            sd = pcs[pc].std()
            mask_pc = ~pcs[pc].between(mean - outlier_sd_threshold * sd,
                                        mean + outlier_sd_threshold * sd)
            outlier_mask = outlier_mask | mask_pc

        outlier_ids = pcs.index[outlier_mask]

        # 5f) Plot PC1 vs PC2 for this iteration, labeling outliers
        fig, ax = plt.subplots()
        ax.scatter(pcs['pPC1'], pcs['pPC2'], alpha=0.7)
        if len(outlier_ids) > 0:
            ax.scatter(pcs.loc[outlier_ids, 'pPC1'], pcs.loc[outlier_ids, 'pPC2'], s=50)
            for oid in outlier_ids:
                ax.annotate(str(oid), (pcs.loc[oid, 'pPC1'], pcs.loc[oid, 'pPC2']))
        ax.set(xlabel='PC1', ylabel='PC2', title=f'PCA Iteration {iteration}')
        fig.savefig(f"{plot_prefix}_iter{iteration}.png")
        plt.close(fig)

        # 5g) If no outliers, break; otherwise remove and repeat
        if len(outlier_ids) == 0:
            break
        print(f"Iteration {iteration}, removing outliers: {list(outlier_ids)}")
        Xp = Xp.drop(index=outlier_ids)

    # 6) Save final PCA scores
    pcs.to_csv(output_scores)
    print("Final PCA scores saved to", output_scores)

    # 7) Final PC1 vs PC2 plot (no iteration label)
    if 'pPC1' in pcs.columns and 'pPC2' in pcs.columns:
        fig, ax = plt.subplots()
        ax.scatter(pcs['pPC1'], pcs['pPC2'], alpha=0.7)
        ax.set(xlabel='PC1', ylabel='PC2', title='Proteomic PCA')
        fig.savefig(f"{plot_prefix}.png")
        plt.close(fig)

    # 8) Scree plot
    fig, ax = plt.subplots()
    ax.plot(np.arange(1, pcs.shape[1] + 1), pca.explained_variance_ratio_, marker='o')
    ax.set(xlabel='PC', ylabel='Explained Variance Ratio', title='Scree Plot')
    fig.savefig(f"{plot_prefix}_scree.png")
    plt.close(fig)

    return pcs, pca

if __name__ == "__main__":
    df = pd.read_excel("data/OLINK Disaease Blood Atlas.xlsx")
    pcs, pca = perform_pca(
        df,
        id_name='DAid',
        prot_column_num=9,
        missingness_threshold=0.05,
        n_components=10,
        outlier_sd_threshold=4
    )
    #  move file from temp to report
    os.makedirs("report", exist_ok=True)
    os.rename("temp/olink_pca.png", "report/olink_pca.png")
    os.rename("temp/olink_pca_scree.png", "report/olink_pca_scree.png")

    # make a dataframe for analysis
    df_var = df.iloc[:, 1:9].set_index(df['DAid']) 
    df_var = df_var.dropna(subset=['LEDD']) # Remove 4 LEDD missing individuals
    df_var['LEDD_CAT'] = ['DenovoPD' if (ledd==0) & (diagnosis=='PD') else 
                          'PD' if diagnosis=='PD' else 'Control'
                          for ledd, diagnosis in zip(df_var['LEDD'], df_var['diagnosis'])]
    df_var['LEDD_sq'] = df_var['LEDD']**2
    df_var = df_var[pd.notna(df_var['LEDD'])].copy() 
    df_prot = df.iloc[:, 9:].set_index(df['DAid'])
    pcs_df = pd.concat([df_var, pcs, df_prot], axis=1, join='inner')
    #  remove spaces in the data
    pcs_df.columns = pcs_df.columns.str.replace(' ', '_')
    #  remove spaces in the data
    str_cols = pcs_df.select_dtypes(include=['object']).columns
    for col in str_cols:
        pcs_df[col] = pcs_df[col].str.replace(' ', '_')
    pcs_df.to_csv("temp/ds.csv")
