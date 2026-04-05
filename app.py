"""
GlycoAntibody Studio - Main Application
糖鎖抗原を標的とした抗体設計・シミュレーションプラットフォーム
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from modules.engine import GlycoEngine
from modules.utils import search_refined_pdb, validate_smiles, parse_candidates_csv
from modules.reports import generate_pareto_chart, RankingReport
from modules.visualization import render_3d_structure

# ── Page Configuration ──────────────────────────────────────────────
st.set_page_config(
    page_title="GlycoAntibody Studio",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Session State Initialization ────────────────────────────────────
_DEFAULT_STATE = {
    "results_df": None,
    "running": False,
    "selected_candidate": None,
    "pdb_structures": {},
    "error_message": None,
}
for key, default in _DEFAULT_STATE.items():
    if key not in st.session_state:
        st.session_state[key] = default


# ── Helper Functions ────────────────────────────────────────────────
def _run_batch_processing(
    target_smiles: str, candidates_df: pd.DataFrame
) -> pd.DataFrame:
    """
    候補CDR3配列のバッチ処理を実行する。

    Parameters
    ----------
    target_smiles : str
        ターゲット糖鎖のSMILES文字列
    candidates_df : pd.DataFrame
        CDR3候補配列を含むDataFrame

    Returns
    -------
    pd.DataFrame
        スコアリング結果を含むDataFrame
    """
    engine = GlycoEngine(target_smiles=target_smiles)
    results = []

    progress_bar = st.progress(0, text="Initializing...")
    total = len(candidates_df)

    for idx, row in candidates_df.iterrows():
        cdr3_seq = row["cdr3_seq"]
        progress_bar.progress(
            (idx + 1) / total,
            text=f"Processing {cdr3_seq} ({idx + 1}/{total})",
        )

        try:
            result = engine.run_pipeline(cdr3_seq)
            results.append(result)
        except Exception as e:
            st.warning(f"⚠️ {cdr3_seq} の処理中にエラー: {e}")
            results.append(
                {
                    "candidate_id": f"CAN_{idx + 1:03d}",
                    "cdr3_seq": cdr3_seq,
                    "calibrated_score": np.nan,
                    "hydrophobicity": np.nan,
                    "h_bonds": np.nan,
                    "status": "failed",
                }
            )

    progress_bar.empty()
    return pd.DataFrame(results)


def _identify_pareto_front(df: pd.DataFrame) -> pd.Series:
    """
    パレート最適解を特定する（結合親和性の最大化 × 疎水性の最小化）。

    Parameters
    ----------
    df : pd.DataFrame
        calibrated_score と hydrophobicity 列を含むDataFrame

    Returns
    -------
    pd.Series
        各行がパレートフロントに含まれるかどうかのbool Series
    """
    scores = df[["calibrated_score", "hydrophobicity"]].dropna().values
    is_pareto = np.ones(len(scores), dtype=bool)

    for i, (score_i, hydro_i) in enumerate(scores):
        if not is_pareto[i]:
            continue
        for j, (score_j, hydro_j) in enumerate(scores):
            if i == j:
                continue
            # score_j がより低い（より強い結合）かつ hydro_j も低い場合、iは支配される
            if score_j <= score_i and hydro_j <= hydro_i and (
                score_j < score_i or hydro_j < hydro_i
            ):
                is_pareto[i] = False
                break

    result = pd.Series(False, index=df.index)
    valid_idx = df[["calibrated_score", "hydrophobicity"]].dropna().index
    result.loc[valid_idx] = is_pareto
    return result


# ── Sidebar ─────────────────────────────────────────────────────────
with st.sidebar:
    st.header("🔬 Project Setup")

    target_smiles = st.text_input(
        "Target Glycan (SMILES)",
        placeholder="例: OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
        help="ターゲット糖鎖のSMILES表記を入力してください。",
    )

    uploaded_file = st.file_uploader(
        "Upload Candidates CSV",
        type=["csv"],
        help="CDR3配列を含むCSVファイル（cdr3_seq 列が必要）。",
    )

    st.divider()

    # Advanced Settings
    with st.expander("⚙️ Advanced Settings"):
        n_conformers = st.slider(
            "Glycan Conformers", min_value=5, max_value=50, value=20,
            help="糖鎖コンフォメーションアンサンブルの数",
        )
        use_conserved_waters = st.checkbox(
            "Include Conserved Waters", value=True,
            help="保存水を考慮したドッキングを実行",
        )
        scoring_weight_hbond = st.slider(
            "H-bond Weight", min_value=0.5, max_value=2.0, value=1.5, step=0.1,
            help="水素結合スコアの重み係数",
        )

    st.divider()

    # Start Button
    can_start = bool(target_smiles) and uploaded_file is not None
    if st.button(
        "🚀 Start Batch Processing",
        disabled=not can_start,
        use_container_width=True,
        type="primary",
    ):
        # 入力バリデーション
        if not validate_smiles(target_smiles):
            st.error("❌ 無効なSMILES文字列です。")
        else:
            candidates_df = parse_candidates_csv(uploaded_file)
            if candidates_df is None:
                st.error("❌ CSVに `cdr3_seq` 列が見つかりません。")
            elif len(candidates_df) == 0:
                st.error("❌ 有効な候補配列がありません。")
            else:
                st.session_state.running = True
                st.session_state.error_message = None

    if not can_start:
        st.caption("SMILES とCSVファイルを入力してください。")


# ── Main Content ────────────────────────────────────────────────────
st.title("🧪 GlycoAntibody Studio")
st.caption("糖鎖抗原ターゲット抗体の自動設計・解析プラットフォーム")

# Batch Processing Execution
if st.session_state.running:
    with st.status("⏳ Batch Processing...", expanded=True) as status:
        try:
            candidates_df = parse_candidates_csv(uploaded_file)
            results_df = _run_batch_processing(target_smiles, candidates_df)
            results_df["is_pareto"] = _identify_pareto_front(results_df)
            st.session_state.results_df = results_df
            status.update(label="✅ Analysis Complete!", state="complete")
        except Exception as e:
            st.session_state.error_message = str(e)
            status.update(label="❌ Analysis Failed", state="error")
        finally:
            st.session_state.running = False

if st.session_state.error_message:
    st.error(f"エラーが発生しました: {st.session_state.error_message}")

# ── Tabs ────────────────────────────────────────────────────────────
tab_analytics, tab_3d, tab_report = st.tabs(
    ["📊 Analytics", "🧬 3D Viewer", "📑 Report"]
)

df = st.session_state.results_df

# ── Tab 1: Analytics ────────────────────────────────────────────────
with tab_analytics:
    if df is None:
        st.info("👈 サイドバーからターゲット糖鎖と候補CSVを設定して解析を開始してください。")
    else:
        col_chart, col_table = st.columns([3, 2])

        with col_chart:
            st.subheader("Pareto Analysis")
            fig = px.scatter(
                df.dropna(subset=["calibrated_score", "hydrophobicity"]),
                x="calibrated_score",
                y="hydrophobicity",
                color="h_bonds",
                symbol="is_pareto",
                symbol_map={True: "star", False: "circle"},
                hover_data=["candidate_id", "cdr3_seq"],
                labels={
                    "calibrated_score": "Calibrated Binding Score (kcal/mol)",
                    "hydrophobicity": "Hydrophobicity Index",
                    "h_bonds": "H-bonds",
                    "is_pareto": "Pareto Optimal",
                },
                color_continuous_scale="Viridis",
            )
            fig.update_layout(
                xaxis_title="Calibrated Score (lower = stronger binding)",
                yaxis_title="Hydrophobicity (lower = better)",
                template="plotly_white",
            )
            st.plotly_chart(fig, use_container_width=True)

        with col_table:
            st.subheader("Ranking Table")
            display_df = df.sort_values("calibrated_score").reset_index(drop=True)
            display_df.index += 1  # 1-based ranking
            display_df.index.name = "Rank"

            st.dataframe(
                display_df[
                    ["candidate_id", "cdr3_seq", "calibrated_score", "hydrophobicity", "h_bonds", "is_pareto"]
                ].style.format(
                    {"calibrated_score": "{:.2f}", "hydrophobicity": "{:.1f}"}
                ).map(
                    lambda v: "background-color: #e6f3ff" if v else "",
                    subset=["is_pareto"],
                ),
                use_container_width=True,
            )

        # Summary metrics
        st.divider()
        valid_df = df.dropna(subset=["calibrated_score"])
        pareto_count = df["is_pareto"].sum() if "is_pareto" in df.columns else 0

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Total Candidates", len(df))
        m2.metric("Successful", len(valid_df))
        m3.metric("Pareto Optimal", int(pareto_count))
        m4.metric(
            "Best Score",
            f"{valid_df['calibrated_score'].min():.2f}" if len(valid_df) > 0 else "N/A",
        )


# ── Tab 2: 3D Viewer ───────────────────────────────────────────────
with tab_3d:
    if df is None:
        st.info("解析完了後に3D構造を表示できます。")
    else:
        selected_id = st.selectbox(
            "Select Candidate",
            options=df["candidate_id"].tolist(),
            format_func=lambda x: f"{x} — {df.loc[df['candidate_id'] == x, 'cdr3_seq'].values[0]}",
        )

        if selected_id:
            st.session_state.selected_candidate = selected_id
            candidate_row = df[df["candidate_id"] == selected_id].iloc[0]

            col_info, col_view = st.columns([1, 3])

            with col_info:
                st.markdown(f"**Candidate:** {candidate_row['candidate_id']}")
                st.markdown(f"**CDR3:** `{candidate_row['cdr3_seq']}`")
                st.markdown(f"**Score:** {candidate_row['calibrated_score']:.2f} kcal/mol")
                st.markdown(f"**Hydrophobicity:** {candidate_row['hydrophobicity']:.1f}")
                st.markdown(f"**H-bonds:** {candidate_row['h_bonds']}")
                if candidate_row.get("is_pareto", False):
                    st.success("⭐ Pareto Optimal")

            with col_view:
                try:
                    render_3d_structure(selected_id, st.session_state.pdb_structures)
                except Exception:
                    st.warning(
                        "3D構造の表示にはPDB構造データが必要です。"
                        "実際のパイプライン実行時に生成されます。"
                    )


# ── Tab 3: Report ───────────────────────────────────────────────────
with tab_report:
    if df is None:
        st.info("解析完了後にレポートを生成できます。")
    else:
        st.subheader("📄 Analysis Report")

        col_preview, col_download = st.columns([3, 1])

        with col_preview:
            st.markdown("### Summary")
            valid_df = df.dropna(subset=["calibrated_score"])
            pareto_df = df[df["is_pareto"]] if "is_pareto" in df.columns else pd.DataFrame()

            st.markdown(
                f"- **解析候補数:** {len(df)}\n"
                f"- **成功:** {len(valid_df)} / **失敗:** {len(df) - len(valid_df)}\n"
                f"- **パレート最適候補:** {len(pareto_df)}\n"
                f"- **最良スコア:** {valid_df['calibrated_score'].min():.2f} kcal/mol"
                if len(valid_df) > 0
                else "有効な結果がありません。"
            )

            if len(pareto_df) > 0:
                st.markdown("### Pareto Optimal Candidates")
                st.dataframe(
                    pareto_df[["candidate_id", "cdr3_seq", "calibrated_score", "hydrophobicity", "h_bonds"]],
                    use_container_width=True,
                )

        with col_download:
            st.markdown("### Export")

            # CSV Download
            csv_data = df.to_csv(index=False).encode("utf-8")
            st.download_button(
                "📥 Download CSV",
                data=csv_data,
                file_name="glycoantibody_results.csv",
                mime="text/csv",
                use_container_width=True,
            )

            # PDF Report
            try:
                report = RankingReport(df, target_smiles=target_smiles)
                pdf_bytes = report.generate_pdf()
                st.download_button(
                    "📥 Download PDF Report",
                    data=pdf_bytes,
                    file_name="glycoantibody_report.pdf",
                    mime="application/pdf",
                    use_container_width=True,
                )
            except Exception as e:
                st.warning(f"PDF生成エラー: {e}")
