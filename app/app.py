"""
Gene Expression Data Explorer
==============================
A Streamlit app for exploring CSV-based gene expression data — built for a
Python & Biopython workshop.

Run with:
    streamlit run app/app.py

Expected CSV shape:
    One row per gene, one column per sample, plus a gene identifier column.
    Example:
        GeneID, GeneSymbol, Control_1, Control_2, ..., Treatment_1, ...
"""

import io

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from scipy import stats

st.set_page_config(
    page_title="Gene Expression Data Explorer",
    page_icon="🧬",
    layout="wide",
)

SAMPLE_CSV_PATH = "data/gene_expression_sample.csv"


# --------------------------------------------------------------------------
# Data loading helpers
# --------------------------------------------------------------------------
@st.cache_data
def load_csv(file_or_path) -> pd.DataFrame:
    return pd.read_csv(file_or_path)


def guess_id_columns(df: pd.DataFrame) -> list[str]:
    """Guess which columns are gene identifiers (non-numeric columns)."""
    return [c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c])]


def compute_de_stats(
    df: pd.DataFrame, group_a_cols: list[str], group_b_cols: list[str]
) -> pd.DataFrame:
    """Compute mean expression, log2 fold change, and t-test p-value per gene."""
    a = df[group_a_cols].to_numpy(dtype=float)
    b = df[group_b_cols].to_numpy(dtype=float)

    mean_a = a.mean(axis=1)
    mean_b = b.mean(axis=1)
    log2fc = mean_b - mean_a  # data already assumed log2-scale

    t_stat, p_val = stats.ttest_ind(b, a, axis=1, equal_var=False)
    p_val = np.nan_to_num(p_val, nan=1.0)
    neg_log10_p = -np.log10(np.clip(p_val, 1e-300, None))

    out = pd.DataFrame(
        {
            "mean_groupA": mean_a,
            "mean_groupB": mean_b,
            "log2FC": log2fc,
            "p_value": p_val,
            "neg_log10_p": neg_log10_p,
        },
        index=df.index,
    )
    return out


# --------------------------------------------------------------------------
# Sidebar — data source & group selection
# --------------------------------------------------------------------------
st.sidebar.title("🧬 Data Setup")

data_source = st.sidebar.radio(
    "Data source", ["Use sample dataset", "Upload my own CSV"], index=0
)

if data_source == "Upload my own CSV":
    uploaded = st.sidebar.file_uploader("Upload gene expression CSV", type=["csv"])
    if uploaded is not None:
        df = load_csv(uploaded)
    else:
        st.sidebar.info("Upload a CSV to get started, or switch to the sample dataset.")
        st.stop()
else:
    df = load_csv(SAMPLE_CSV_PATH)

id_cols_guess = guess_id_columns(df)
numeric_cols = [c for c in df.columns if c not in id_cols_guess]

st.sidebar.markdown("---")
gene_label_col = st.sidebar.selectbox(
    "Gene label column",
    options=id_cols_guess if id_cols_guess else df.columns.tolist(),
    index=(id_cols_guess.index("GeneSymbol") if "GeneSymbol" in id_cols_guess else 0),
)

st.sidebar.markdown("---")
st.sidebar.subheader("Sample Groups")
default_a = [c for c in numeric_cols if "control" in c.lower()]
default_b = [c for c in numeric_cols if "treatment" in c.lower()]

group_a_cols = st.sidebar.multiselect(
    "Group A samples", options=numeric_cols, default=default_a or numeric_cols[: len(numeric_cols) // 2]
)
group_b_cols = st.sidebar.multiselect(
    "Group B samples",
    options=[c for c in numeric_cols if c not in group_a_cols],
    default=default_b or [c for c in numeric_cols if c not in (default_a or [])],
)

# --------------------------------------------------------------------------
# Header
# --------------------------------------------------------------------------
st.title("🧬 Gene Expression Data Explorer")
st.caption(
    "Explore CSV-based gene expression data: browse, search, visualize, and "
    "compare expression between sample groups."
)

col1, col2, col3, col4 = st.columns(4)
col1.metric("Genes", f"{df.shape[0]:,}")
col2.metric("Samples", f"{len(numeric_cols):,}")
col3.metric("Group A samples", len(group_a_cols))
col4.metric("Group B samples", len(group_b_cols))

tabs = st.tabs(
    ["📋 Data Table", "🔍 Gene Lookup", "🌋 Differential Expression", "🔥 Heatmap", "📊 PCA"]
)

# --------------------------------------------------------------------------
# Tab 1: Data table
# --------------------------------------------------------------------------
with tabs[0]:
    st.subheader("Browse & Filter")
    search = st.text_input("Search by gene label", "")
    view_df = df.copy()
    if search:
        view_df = view_df[
            view_df[gene_label_col].astype(str).str.contains(search, case=False, na=False)
        ]
    st.dataframe(view_df, use_container_width=True, height=450)
    st.download_button(
        "Download filtered data as CSV",
        data=view_df.to_csv(index=False).encode("utf-8"),
        file_name="filtered_gene_expression.csv",
        mime="text/csv",
    )

# --------------------------------------------------------------------------
# Tab 2: Single gene lookup
# --------------------------------------------------------------------------
with tabs[1]:
    st.subheader("Look up a single gene")
    gene_options = df[gene_label_col].astype(str).tolist()
    selected_gene = st.selectbox("Select a gene", options=gene_options)

    row = df[df[gene_label_col].astype(str) == selected_gene].iloc[0]
    plot_data = []
    for c in group_a_cols:
        plot_data.append({"Sample": c, "Expression": row[c], "Group": "Group A"})
    for c in group_b_cols:
        plot_data.append({"Sample": c, "Expression": row[c], "Group": "Group B"})
    plot_df = pd.DataFrame(plot_data)

    left, right = st.columns([2, 1])
    with left:
        fig = px.box(
            plot_df,
            x="Group",
            y="Expression",
            points="all",
            color="Group",
            title=f"Expression of {selected_gene} across groups",
        )
        st.plotly_chart(fig, use_container_width=True)
    with right:
        st.write("**Raw values**")
        st.dataframe(plot_df.set_index("Sample"), use_container_width=True)
        if group_a_cols and group_b_cols:
            mean_a = row[group_a_cols].mean()
            mean_b = row[group_b_cols].mean()
            st.metric("Mean Group A", f"{mean_a:.2f}")
            st.metric("Mean Group B", f"{mean_b:.2f}")
            st.metric("log2 Fold Change (B vs A)", f"{mean_b - mean_a:+.2f}")

# --------------------------------------------------------------------------
# Tab 3: Differential expression / volcano plot
# --------------------------------------------------------------------------
with tabs[2]:
    st.subheader("Differential Expression: Group B vs Group A")
    if len(group_a_cols) < 2 or len(group_b_cols) < 2:
        st.warning("Select at least 2 samples per group in the sidebar to compute statistics.")
    else:
        de = compute_de_stats(df, group_a_cols, group_b_cols)
        de[gene_label_col] = df[gene_label_col].values

        c1, c2 = st.columns(2)
        fc_threshold = c1.slider("log2FC threshold", 0.0, 5.0, 1.0, 0.1)
        p_threshold = c2.slider("p-value threshold", 0.001, 0.1, 0.05, 0.001)

        def classify(r):
            if r["p_value"] < p_threshold and r["log2FC"] >= fc_threshold:
                return "Up in Group B"
            elif r["p_value"] < p_threshold and r["log2FC"] <= -fc_threshold:
                return "Down in Group B"
            return "Not significant"

        de["Status"] = de.apply(classify, axis=1)

        fig = px.scatter(
            de,
            x="log2FC",
            y="neg_log10_p",
            color="Status",
            hover_name=gene_label_col,
            color_discrete_map={
                "Up in Group B": "#d62728",
                "Down in Group B": "#1f77b4",
                "Not significant": "#B0B0B0",
            },
            labels={"neg_log10_p": "-log10(p-value)"},
            title="Volcano Plot",
        )
        fig.add_vline(x=fc_threshold, line_dash="dash", line_color="gray")
        fig.add_vline(x=-fc_threshold, line_dash="dash", line_color="gray")
        fig.add_hline(y=-np.log10(p_threshold), line_dash="dash", line_color="gray")
        st.plotly_chart(fig, use_container_width=True)

        st.write(f"**{(de['Status'] != 'Not significant').sum()} genes** pass both thresholds.")
        st.dataframe(
            de[de["Status"] != "Not significant"]
            .sort_values("p_value")[[gene_label_col, "mean_groupA", "mean_groupB", "log2FC", "p_value", "Status"]],
            use_container_width=True,
        )

# --------------------------------------------------------------------------
# Tab 4: Heatmap of top variable genes
# --------------------------------------------------------------------------
with tabs[3]:
    st.subheader("Heatmap of Top Variable Genes")
    n_top = st.slider("Number of top variable genes", 5, 50, 20)
    all_sample_cols = group_a_cols + group_b_cols
    if len(all_sample_cols) < 2:
        st.warning("Select samples in the sidebar to build a heatmap.")
    else:
        variances = df[all_sample_cols].var(axis=1)
        top_idx = variances.sort_values(ascending=False).head(n_top).index
        heat_df = df.loc[top_idx, all_sample_cols]
        heat_df.index = df.loc[top_idx, gene_label_col]

        # z-score per gene (row) for visual comparability
        z = heat_df.sub(heat_df.mean(axis=1), axis=0).div(heat_df.std(axis=1) + 1e-9, axis=0)

        fig = go.Figure(
            data=go.Heatmap(
                z=z.values,
                x=z.columns,
                y=z.index.astype(str),
                colorscale="RdBu_r",
                zmid=0,
                colorbar=dict(title="z-score"),
            )
        )
        fig.update_layout(height=max(400, n_top * 20), title="Top variable genes (row z-score)")
        st.plotly_chart(fig, use_container_width=True)

# --------------------------------------------------------------------------
# Tab 5: PCA
# --------------------------------------------------------------------------
with tabs[4]:
    st.subheader("Principal Component Analysis (Samples)")
    all_sample_cols = group_a_cols + group_b_cols
    if len(all_sample_cols) < 3:
        st.warning("Select at least 3 samples in the sidebar to run PCA.")
    else:
        X = df[all_sample_cols].to_numpy(dtype=float).T  # samples x genes
        X_centered = X - X.mean(axis=0)

        # Simple PCA via SVD (no sklearn dependency needed)
        U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
        scores = U * S
        explained_var = (S**2) / np.sum(S**2)

        pca_df = pd.DataFrame(
            {
                "PC1": scores[:, 0],
                "PC2": scores[:, 1] if scores.shape[1] > 1 else np.zeros(len(scores)),
                "Sample": all_sample_cols,
                "Group": ["Group A"] * len(group_a_cols) + ["Group B"] * len(group_b_cols),
            }
        )

        fig = px.scatter(
            pca_df,
            x="PC1",
            y="PC2",
            color="Group",
            text="Sample",
            title="PCA of Samples",
            labels={
                "PC1": f"PC1 ({explained_var[0]*100:.1f}% variance)",
                "PC2": f"PC2 ({explained_var[1]*100:.1f}% variance)" if len(explained_var) > 1 else "PC2",
            },
        )
        fig.update_traces(textposition="top center", marker=dict(size=12))
        st.plotly_chart(fig, use_container_width=True)

st.sidebar.markdown("---")
st.sidebar.caption("Built for a Python & Biopython workshop 🧬")
