"""
GlycoAntibody Studio - 3D Visualization
Py3Dmol / stmol を用いた抗体-糖鎖複合体の3D表示
"""

from __future__ import annotations

import logging
from typing import Optional

import streamlit as st

logger = logging.getLogger(__name__)


# ── 3D Structure Rendering ──────────────────────────────────────────

def render_3d_structure(
    candidate_id: str,
    pdb_structures: dict[str, str],
    height: int = 500,
    width: int = 700,
    style: str = "cartoon",
) -> None:
    """
    Streamlit上にPDB構造の3Dインタラクティブビューを表示する。

    stmol (Py3Dmol wrapper) が利用可能な場合はインタラクティブ表示、
    なければNGL Viewerのiframeフォールバックを使用する。

    Parameters
    ----------
    candidate_id : str
        表示する候補のID
    pdb_structures : dict[str, str]
        candidate_id → PDB文字列のマッピング
    height : int
        ビューアの高さ（px）
    width : int
        ビューアの幅（px）
    style : str
        表示スタイル（'cartoon', 'stick', 'sphere', 'surface'）
    """
    pdb_string = pdb_structures.get(candidate_id)

    if pdb_string is None:
        _render_placeholder(candidate_id)
        return

    try:
        _render_with_stmol(pdb_string, height, width, style)
    except ImportError:
        try:
            _render_with_py3dmol_html(pdb_string, height, width, style)
        except Exception:
            _render_pdb_text(pdb_string, candidate_id)


def _render_with_stmol(
    pdb_string: str,
    height: int,
    width: int,
    style: str,
) -> None:
    """stmol (streamlit-py3dmol) による3D表示。"""
    from stmol import showmol
    import py3Dmol

    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_string, "pdb")

    # スタイル適用
    style_map = _build_style_config(style)
    view.setStyle(style_map)

    # CDR領域のハイライト（IMGT番号付け基準）
    _highlight_cdr_regions(view)

    view.zoomTo()
    view.setBackgroundColor("white")

    showmol(view, height=height, width=width)


def _render_with_py3dmol_html(
    pdb_string: str,
    height: int,
    width: int,
    style: str,
) -> None:
    """py3Dmolを直接HTMLとして埋め込む方式。"""
    import py3Dmol

    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_string, "pdb")

    style_map = _build_style_config(style)
    view.setStyle(style_map)
    _highlight_cdr_regions(view)

    view.zoomTo()
    view.setBackgroundColor("white")

    # HTML文字列として取得し、components.htmlで表示
    html = view._make_html()
    st.components.v1.html(html, height=height + 20, width=width + 20)


def _render_pdb_text(pdb_string: str, candidate_id: str) -> None:
    """テキストベースのフォールバック表示。"""
    st.warning("3Dビューアライブラリが利用できません。PDB内容をテキストで表示します。")
    st.caption(f"Structure: {candidate_id}")

    # 基本統計を抽出
    lines = pdb_string.strip().split("\n")
    atom_lines = [l for l in lines if l.startswith("ATOM") or l.startswith("HETATM")]
    residues = set()
    chains = set()
    for line in atom_lines:
        if len(line) >= 22:
            chains.add(line[21])
        if len(line) >= 26:
            residues.add(line[17:20].strip() + line[22:26].strip())

    col1, col2, col3 = st.columns(3)
    col1.metric("Atoms", len(atom_lines))
    col2.metric("Residues", len(residues))
    col3.metric("Chains", len(chains))

    with st.expander("PDB Content (first 50 lines)"):
        st.code("\n".join(lines[:50]), language="text")


def _render_placeholder(candidate_id: str) -> None:
    """PDB構造がない場合のプレースホルダー表示。"""
    st.info(
        f"📦 {candidate_id} の3D構造はまだ生成されていません。\n\n"
        "実際のパイプライン実行（OpenMM + Vinaドッキング）により生成されます。\n"
        "現在はシミュレーションモードで動作しています。"
    )

    # テンプレート構造（トラスツズマブ）の表示を提供
    with st.expander("💡 テンプレート構造（Trastuzumab Fab: 1N8Z）を表示"):
        _render_template_structure()


def _render_template_structure() -> None:
    """RCSB PDBからテンプレート構造を取得して表示。"""
    try:
        import py3Dmol

        view = py3Dmol.view(width=600, height=400)
        view.addModel("", "pdb")
        # PDB IDから直接取得
        view.setQuery({"url": "https://files.rcsb.org/download/1N8Z.pdb"})
        view.setStyle({"cartoon": {"color": "spectrum"}})
        view.zoomTo()

        try:
            from stmol import showmol
            showmol(view, height=400, width=600)
        except ImportError:
            html = view._make_html()
            st.components.v1.html(html, height=420, width=620)
    except ImportError:
        st.caption(
            "py3Dmol がインストールされていません。\n"
            "`pip install py3Dmol stmol` で3Dビューアを有効化できます。"
        )


# ── Style Configuration ────────────────────────────────────────────

def _build_style_config(style: str) -> dict:
    """Py3Dmolスタイル設定を構築する。"""
    configs = {
        "cartoon": {"cartoon": {"color": "spectrum", "opacity": 0.85}},
        "stick": {"stick": {"colorscheme": "amino", "radius": 0.15}},
        "sphere": {"sphere": {"colorscheme": "amino", "scale": 0.3}},
        "surface": {"cartoon": {"color": "white", "opacity": 0.3}},
        "ball_and_stick": {
            "stick": {"radius": 0.1},
            "sphere": {"scale": 0.25, "colorscheme": "amino"},
        },
    }
    return configs.get(style, configs["cartoon"])


def _highlight_cdr_regions(view) -> None:
    """
    CDR領域をハイライト表示する。

    IMGT番号付けに基づくCDR位置:
      CDR1: 27-38, CDR2: 56-65, CDR3: 105-117
    """
    cdr_colors = {
        "CDR1": {"start": 27, "end": 38, "color": "#FF6B6B"},   # 赤系
        "CDR2": {"start": 56, "end": 65, "color": "#4ECDC4"},   # 青緑
        "CDR3": {"start": 105, "end": 117, "color": "#FFE66D"},  # 黄色
    }

    for cdr_name, info in cdr_colors.items():
        view.addStyle(
            {"resi": list(range(info["start"], info["end"] + 1))},
            {
                "cartoon": {"color": info["color"]},
                "stick": {"color": info["color"], "radius": 0.15},
            },
        )


# ── Style Selector Widget ──────────────────────────────────────────

def style_selector() -> str:
    """Streamlitウィジェットで3D表示スタイルを選択する。"""
    return st.selectbox(
        "Display Style",
        options=["cartoon", "stick", "sphere", "ball_and_stick", "surface"],
        format_func=lambda x: {
            "cartoon": "🎗️ Cartoon (Ribbon)",
            "stick": "🪵 Stick",
            "sphere": "⚪ Sphere (CPK)",
            "ball_and_stick": "🔵 Ball & Stick",
            "surface": "🫧 Surface",
        }.get(x, x),
    )


# ── Contact Map ─────────────────────────────────────────────────────

def render_contact_map(
    pdb_string: str,
    distance_cutoff: float = 4.5,
) -> None:
    """
    残基間コンタクトマップを表示する。

    Parameters
    ----------
    pdb_string : str
        PDBファイルの内容
    distance_cutoff : float
        コンタクト判定距離（Å）
    """
    try:
        import numpy as np
        import plotly.express as px

        # CA原子の座標を抽出
        coords = []
        residue_labels = []
        for line in pdb_string.split("\n"):
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                resname = line[17:20].strip()
                resnum = line[22:26].strip()
                coords.append([x, y, z])
                residue_labels.append(f"{resname}{resnum}")

        if len(coords) < 2:
            st.warning("コンタクトマップ生成に十分な残基が見つかりません。")
            return

        coords = np.array(coords)
        n = len(coords)

        # 距離行列の計算
        dist_matrix = np.sqrt(
            np.sum((coords[:, np.newaxis, :] - coords[np.newaxis, :, :]) ** 2, axis=2)
        )
        contact_matrix = (dist_matrix < distance_cutoff).astype(float)

        fig = px.imshow(
            contact_matrix,
            x=residue_labels,
            y=residue_labels,
            color_continuous_scale="Blues",
            labels=dict(color="Contact"),
            title=f"Contact Map (cutoff: {distance_cutoff} Å)",
        )
        fig.update_layout(height=500, width=500)
        st.plotly_chart(fig, use_container_width=True)

    except Exception as e:
        st.warning(f"コンタクトマップの生成に失敗しました: {e}")
