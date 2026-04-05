"""
GlycoAntibody Studio - Report Generation
パレート解析チャートおよびPDFレポートの生成
"""

from __future__ import annotations

import io
import logging
from datetime import datetime
from typing import Optional

import matplotlib
matplotlib.use("Agg")  # 非GUIバックエンド

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

logger = logging.getLogger(__name__)


# ── Pareto Chart ────────────────────────────────────────────────────

def generate_pareto_chart(
    df: pd.DataFrame,
    x_col: str = "calibrated_score",
    y_col: str = "hydrophobicity",
    color_col: str = "h_bonds",
    pareto_col: str = "is_pareto",
) -> go.Figure:
    """
    パレート解析の散布図を生成する。

    Parameters
    ----------
    df : pd.DataFrame
        解析結果DataFrame
    x_col : str
        X軸カラム（結合スコア）
    y_col : str
        Y軸カラム（疎水性）
    color_col : str
        色分けカラム（水素結合数）
    pareto_col : str
        パレートフロントフラグのカラム

    Returns
    -------
    go.Figure
        Plotly Figureオブジェクト
    """
    plot_df = df.dropna(subset=[x_col, y_col]).copy()

    if plot_df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No valid data to display",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False, font_size=16,
        )
        return fig

    # Non-Pareto points
    non_pareto = plot_df[~plot_df.get(pareto_col, False)]
    pareto = plot_df[plot_df.get(pareto_col, False)]

    fig = go.Figure()

    # Non-Pareto散布点
    if len(non_pareto) > 0:
        fig.add_trace(
            go.Scatter(
                x=non_pareto[x_col],
                y=non_pareto[y_col],
                mode="markers",
                marker=dict(
                    size=10,
                    color=non_pareto[color_col],
                    colorscale="Viridis",
                    colorbar=dict(title="H-bonds"),
                    opacity=0.6,
                    line=dict(width=0.5, color="white"),
                ),
                text=non_pareto.apply(
                    lambda r: f"{r['candidate_id']}<br>{r['cdr3_seq']}", axis=1
                ),
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    f"{x_col}: %{{x:.2f}}<br>"
                    f"{y_col}: %{{y:.1f}}<br>"
                    "<extra></extra>"
                ),
                name="Candidates",
            )
        )

    # Pareto最適点
    if len(pareto) > 0:
        fig.add_trace(
            go.Scatter(
                x=pareto[x_col],
                y=pareto[y_col],
                mode="markers",
                marker=dict(
                    size=14,
                    color=pareto[color_col],
                    colorscale="Viridis",
                    symbol="star",
                    line=dict(width=1.5, color="gold"),
                ),
                text=pareto.apply(
                    lambda r: f"{r['candidate_id']}<br>{r['cdr3_seq']}", axis=1
                ),
                hovertemplate=(
                    "<b>⭐ %{text}</b><br>"
                    f"{x_col}: %{{x:.2f}}<br>"
                    f"{y_col}: %{{y:.1f}}<br>"
                    "<extra></extra>"
                ),
                name="Pareto Optimal",
            )
        )

        # パレートフロントの線
        pareto_sorted = pareto.sort_values(x_col)
        fig.add_trace(
            go.Scatter(
                x=pareto_sorted[x_col],
                y=pareto_sorted[y_col],
                mode="lines",
                line=dict(color="rgba(255, 215, 0, 0.5)", width=2, dash="dash"),
                showlegend=False,
                hoverinfo="skip",
            )
        )

    fig.update_layout(
        title="Pareto Analysis: Affinity vs. Hydrophobicity",
        xaxis_title="Calibrated Binding Score (kcal/mol, lower = stronger)",
        yaxis_title="Hydrophobicity Index (lower = less aggregation risk)",
        template="plotly_white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        height=500,
    )

    # 理想領域のハイライト
    x_range = plot_df[x_col]
    y_range = plot_df[y_col]
    fig.add_shape(
        type="rect",
        x0=x_range.min() - 0.5,
        y0=y_range.min() - 1,
        x1=x_range.quantile(0.25),
        y1=y_range.quantile(0.25),
        fillcolor="rgba(0, 255, 0, 0.05)",
        line=dict(color="rgba(0, 200, 0, 0.3)", dash="dot"),
    )
    fig.add_annotation(
        x=x_range.quantile(0.15),
        y=y_range.quantile(0.15),
        text="Ideal Region",
        showarrow=False,
        font=dict(size=10, color="green"),
        opacity=0.6,
    )

    return fig


def generate_pareto_chart_matplotlib(
    df: pd.DataFrame,
    x_col: str = "calibrated_score",
    y_col: str = "hydrophobicity",
) -> bytes:
    """
    Matplotlib版パレート図（PDF埋め込み用）。

    Returns
    -------
    bytes
        PNGイメージのバイト列
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    plot_df = df.dropna(subset=[x_col, y_col])
    is_pareto = plot_df.get("is_pareto", pd.Series(False, index=plot_df.index))

    # 通常点
    non_pareto = plot_df[~is_pareto]
    ax.scatter(
        non_pareto[x_col], non_pareto[y_col],
        c=non_pareto.get("h_bonds", "steelblue"),
        cmap="viridis", alpha=0.6, s=60, edgecolors="white", linewidth=0.5,
        label="Candidates",
    )

    # パレート点
    pareto = plot_df[is_pareto]
    if len(pareto) > 0:
        ax.scatter(
            pareto[x_col], pareto[y_col],
            c=pareto.get("h_bonds", "gold"),
            cmap="viridis", marker="*", s=200, edgecolors="gold", linewidth=1.5,
            label="Pareto Optimal",
        )

    ax.set_xlabel("Calibrated Score (kcal/mol)", fontsize=11)
    ax.set_ylabel("Hydrophobicity Index", fontsize=11)
    ax.set_title("Pareto Analysis", fontsize=13, fontweight="bold")
    ax.legend(frameon=True, fancybox=True)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf.read()


# ── PDF Report ──────────────────────────────────────────────────────

class RankingReport:
    """
    解析結果のPDFレポートを生成するクラス。

    Parameters
    ----------
    results_df : pd.DataFrame
        解析結果DataFrame
    target_smiles : str
        ターゲット糖鎖のSMILES
    """

    def __init__(self, results_df: pd.DataFrame, target_smiles: str = ""):
        self.df = results_df
        self.target_smiles = target_smiles
        self.timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def generate_pdf(self) -> bytes:
        """
        PDF形式のレポートを生成する。

        Returns
        -------
        bytes
            PDFファイルのバイト列
        """
        try:
            return self._generate_fpdf()
        except ImportError:
            logger.warning("fpdf not available; generating text-based report")
            return self._generate_text_report().encode("utf-8")

    def _generate_fpdf(self) -> bytes:
        """FPDFを使用したPDF生成。"""
        from fpdf import FPDF

        pdf = FPDF()
        pdf.set_auto_page_break(auto=True, margin=15)

        # ── Cover Page ──
        pdf.add_page()
        pdf.set_font("Helvetica", "B", 24)
        pdf.cell(0, 40, "", ln=True)
        pdf.cell(0, 15, "GlycoAntibody Studio", ln=True, align="C")

        pdf.set_font("Helvetica", "", 14)
        pdf.cell(0, 10, "Analysis Report", ln=True, align="C")
        pdf.cell(0, 8, self.timestamp, ln=True, align="C")

        pdf.set_font("Helvetica", "", 10)
        pdf.cell(0, 20, "", ln=True)
        pdf.cell(0, 6, f"Target SMILES: {self.target_smiles[:80]}", ln=True, align="C")

        # ── Summary Page ──
        pdf.add_page()
        pdf.set_font("Helvetica", "B", 16)
        pdf.cell(0, 12, "1. Analysis Summary", ln=True)
        pdf.ln(5)

        valid_df = self.df.dropna(subset=["calibrated_score"])
        pareto_df = self.df[self.df.get("is_pareto", False)] if "is_pareto" in self.df.columns else pd.DataFrame()

        summary_items = [
            ("Total Candidates", str(len(self.df))),
            ("Successful Analyses", str(len(valid_df))),
            ("Failed Analyses", str(len(self.df) - len(valid_df))),
            ("Pareto Optimal", str(len(pareto_df))),
        ]

        if len(valid_df) > 0:
            summary_items.extend([
                ("Best Score", f"{valid_df['calibrated_score'].min():.3f} kcal/mol"),
                ("Worst Score", f"{valid_df['calibrated_score'].max():.3f} kcal/mol"),
                ("Mean Score", f"{valid_df['calibrated_score'].mean():.3f} kcal/mol"),
                ("Mean H-bonds", f"{valid_df['h_bonds'].mean():.1f}"),
            ])

        pdf.set_font("Helvetica", "", 11)
        for label, value in summary_items:
            pdf.cell(80, 8, f"  {label}:", ln=False)
            pdf.set_font("Helvetica", "B", 11)
            pdf.cell(0, 8, value, ln=True)
            pdf.set_font("Helvetica", "", 11)

        # ── Pareto Chart (embedded image) ──
        pdf.ln(10)
        pdf.set_font("Helvetica", "B", 16)
        pdf.cell(0, 12, "2. Pareto Analysis", ln=True)
        pdf.ln(3)

        try:
            chart_bytes = generate_pareto_chart_matplotlib(self.df)
            chart_path = "/tmp/pareto_chart.png"
            with open(chart_path, "wb") as f:
                f.write(chart_bytes)
            pdf.image(chart_path, x=15, w=180)
        except Exception as e:
            pdf.set_font("Helvetica", "I", 10)
            pdf.cell(0, 8, f"(Chart generation failed: {e})", ln=True)

        # ── Results Table ──
        pdf.add_page()
        pdf.set_font("Helvetica", "B", 16)
        pdf.cell(0, 12, "3. Ranking Table", ln=True)
        pdf.ln(5)

        # Table header
        col_widths = [20, 35, 35, 30, 25, 25]
        headers = ["Rank", "ID", "CDR3", "Score", "Hydro.", "H-bonds"]
        pdf.set_font("Helvetica", "B", 9)
        pdf.set_fill_color(220, 235, 250)
        for w, h in zip(col_widths, headers):
            pdf.cell(w, 8, h, border=1, fill=True, align="C")
        pdf.ln()

        # Table body
        pdf.set_font("Helvetica", "", 8)
        sorted_df = valid_df.sort_values("calibrated_score").reset_index(drop=True)

        for rank, (_, row) in enumerate(sorted_df.iterrows(), 1):
            is_pareto = row.get("is_pareto", False)
            if is_pareto:
                pdf.set_fill_color(255, 250, 205)
            else:
                pdf.set_fill_color(255, 255, 255)

            values = [
                str(rank),
                str(row.get("candidate_id", "")),
                str(row.get("cdr3_seq", ""))[:12],
                f"{row['calibrated_score']:.2f}",
                f"{row.get('hydrophobicity', 0):.1f}",
                str(int(row.get("h_bonds", 0))),
            ]
            for w, v in zip(col_widths, values):
                pdf.cell(w, 7, v, border=1, fill=is_pareto, align="C")
            pdf.ln()

            if rank >= 50:
                pdf.set_font("Helvetica", "I", 8)
                pdf.cell(0, 7, f"  ... and {len(sorted_df) - 50} more candidates", ln=True)
                break

        # ── Pareto Optimal Details ──
        if len(pareto_df) > 0:
            pdf.ln(10)
            pdf.set_font("Helvetica", "B", 16)
            pdf.cell(0, 12, "4. Pareto Optimal Candidates", ln=True)
            pdf.ln(3)

            pdf.set_font("Helvetica", "", 10)
            for _, row in pareto_df.sort_values("calibrated_score").iterrows():
                pdf.set_font("Helvetica", "B", 10)
                pdf.cell(0, 7, f"{row['candidate_id']}  |  CDR3: {row['cdr3_seq']}", ln=True)
                pdf.set_font("Helvetica", "", 9)
                pdf.cell(
                    0, 6,
                    f"    Score: {row['calibrated_score']:.3f} kcal/mol  |  "
                    f"Hydrophobicity: {row.get('hydrophobicity', 0):.2f}  |  "
                    f"H-bonds: {int(row.get('h_bonds', 0))}",
                    ln=True,
                )
                pdf.ln(2)

        # ── Footer ──
        pdf.ln(10)
        pdf.set_font("Helvetica", "I", 8)
        pdf.cell(
            0, 5,
            "Generated by GlycoAntibody Studio | "
            "Scores are computational predictions and require experimental validation.",
            ln=True,
            align="C",
        )

        return pdf.output()

    def _generate_text_report(self) -> str:
        """テキスト形式のフォールバックレポート。"""
        lines = [
            "=" * 60,
            "GlycoAntibody Studio - Analysis Report",
            f"Generated: {self.timestamp}",
            f"Target: {self.target_smiles}",
            "=" * 60,
            "",
        ]

        valid_df = self.df.dropna(subset=["calibrated_score"])
        lines.append(f"Total Candidates: {len(self.df)}")
        lines.append(f"Successful: {len(valid_df)}")

        if len(valid_df) > 0:
            lines.append(f"Best Score: {valid_df['calibrated_score'].min():.3f} kcal/mol")
            lines.append("")
            lines.append("Top 10 Candidates:")
            lines.append("-" * 50)
            for rank, (_, row) in enumerate(
                valid_df.sort_values("calibrated_score").head(10).iterrows(), 1
            ):
                lines.append(
                    f"  {rank:2d}. {row['candidate_id']} | "
                    f"{row['cdr3_seq']} | "
                    f"Score: {row['calibrated_score']:.3f}"
                )

        return "\n".join(lines)
