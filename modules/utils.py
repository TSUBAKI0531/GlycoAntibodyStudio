"""
GlycoAntibody Studio - Utility Functions
入力バリデーション、PDB検索、配列処理などの共通ユーティリティ
"""

from __future__ import annotations

import io
import logging
import re
from typing import Optional

import numpy as np
import pandas as pd
import requests
from Bio.SeqUtils.ProtParam import ProteinAnalysis

logger = logging.getLogger(__name__)

# ── SMILES Validation ───────────────────────────────────────────────

def validate_smiles(smiles: str) -> bool:
    """
    SMILES文字列の妥当性を検証する。

    RDKitが利用可能な場合はMolオブジェクトへの変換を試行し、
    利用不可の場合は基本的な文字列パターンチェックを行う。

    Parameters
    ----------
    smiles : str
        検証対象のSMILES文字列

    Returns
    -------
    bool
        有効な場合True
    """
    if not smiles or not smiles.strip():
        return False

    try:
        from rdkit import Chem

        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            logger.warning("RDKit could not parse SMILES: %s", smiles)
            return False
        return True
    except ImportError:
        logger.info("RDKit not available; using regex-based SMILES validation")
        return _validate_smiles_regex(smiles.strip())


def _validate_smiles_regex(smiles: str) -> bool:
    """RDKit非依存のSMILES簡易バリデーション。"""
    # 許容文字パターン（原子記号、結合、環番号、分岐、立体化学）
    allowed = re.compile(
        r'^[A-Za-z0-9@+\-\[\]\(\)\.\#\=\\/\:%]+$'
    )
    if not allowed.match(smiles):
        return False

    # 括弧の対応チェック
    depth_round = 0
    depth_square = 0
    for ch in smiles:
        if ch == '(':
            depth_round += 1
        elif ch == ')':
            depth_round -= 1
        elif ch == '[':
            depth_square += 1
        elif ch == ']':
            depth_square -= 1
        if depth_round < 0 or depth_square < 0:
            return False

    return depth_round == 0 and depth_square == 0


# ── CSV Parsing ─────────────────────────────────────────────────────

def parse_candidates_csv(
    uploaded_file,
    seq_column: str = "cdr3_seq",
) -> Optional[pd.DataFrame]:
    """
    アップロードされたCSVファイルをパースし、CDR3候補配列を抽出する。

    Parameters
    ----------
    uploaded_file : UploadedFile
        StreamlitのUploadedFileオブジェクト
    seq_column : str
        CDR3配列が格納されている列名（デフォルト: 'cdr3_seq'）

    Returns
    -------
    Optional[pd.DataFrame]
        パース結果のDataFrame。seq_column が見つからない場合はNone。
    """
    if uploaded_file is None:
        return None

    try:
        uploaded_file.seek(0)
        df = pd.read_csv(uploaded_file)
    except Exception as e:
        logger.error("CSV parsing failed: %s", e)
        return None

    # カラム名の正規化（空白除去、小文字化）
    df.columns = df.columns.str.strip().str.lower()
    target_col = seq_column.lower()

    if target_col not in df.columns:
        # 類似カラム名を探す（cdr3, sequence, seq など）
        candidates = [c for c in df.columns if "cdr3" in c or "seq" in c]
        if candidates:
            logger.info("Column '%s' not found; using '%s'", target_col, candidates[0])
            df = df.rename(columns={candidates[0]: target_col})
        else:
            logger.error("No sequence column found in CSV. Columns: %s", list(df.columns))
            return None

    # 配列バリデーション
    df[target_col] = df[target_col].astype(str).str.strip().str.upper()
    valid_mask = df[target_col].apply(_is_valid_amino_acid_seq)
    n_invalid = (~valid_mask).sum()
    if n_invalid > 0:
        logger.warning("Removed %d rows with invalid amino acid sequences", n_invalid)
    df = df[valid_mask].reset_index(drop=True)

    # candidate_id の付与（なければ生成）
    if "candidate_id" not in df.columns:
        df["candidate_id"] = [f"CAN_{i + 1:03d}" for i in range(len(df))]

    return df


def _is_valid_amino_acid_seq(seq: str) -> bool:
    """アミノ酸配列の妥当性を検証する。"""
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    return len(seq) >= 3 and all(c in valid_aa for c in seq)


# ── PDB Utilities ───────────────────────────────────────────────────

# トラスツズマブFab構造 (PDB: 1N8Z) の参照配列
TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)

TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSG"
    "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)


def search_refined_pdb(
    pdb_id: str = "1N8Z",
    timeout: int = 10,
) -> Optional[str]:
    """
    RCSB PDBからPDBファイルを取得する。

    Parameters
    ----------
    pdb_id : str
        PDB ID（デフォルト: 1N8Z = Trastuzumab Fab）
    timeout : int
        リクエストタイムアウト（秒）

    Returns
    -------
    Optional[str]
        PDBファイルの内容文字列。取得失敗時はNone。
    """
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        logger.info("Successfully fetched PDB: %s", pdb_id)
        return response.text
    except requests.RequestException as e:
        logger.error("Failed to fetch PDB %s: %s", pdb_id, e)
        return None


def fetch_pdb_structure(pdb_id: str) -> Optional[str]:
    """PDB構造をmmCIF形式で取得する。"""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        return response.text
    except requests.RequestException as e:
        logger.error("Failed to fetch mmCIF %s: %s", pdb_id, e)
        return None


# ── Sequence Analysis ───────────────────────────────────────────────

def compute_hydrophobicity(sequence: str) -> float:
    """
    Kyte-Doolittle法による平均疎水性（GRAVY）を算出する。

    Parameters
    ----------
    sequence : str
        アミノ酸配列

    Returns
    -------
    float
        GRAVY値（正 = 疎水的、負 = 親水的）
    """
    try:
        analysis = ProteinAnalysis(sequence)
        return analysis.gravy()
    except Exception:
        # 手動フォールバック
        kd_scale = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
        }
        values = [kd_scale.get(aa, 0.0) for aa in sequence.upper()]
        return np.mean(values) if values else 0.0


def compute_isoelectric_point(sequence: str) -> float:
    """等電点を算出する。"""
    try:
        analysis = ProteinAnalysis(sequence)
        return analysis.isoelectric_point()
    except Exception:
        return 0.0


def compute_molecular_weight(sequence: str) -> float:
    """分子量（Da）を算出する。"""
    try:
        analysis = ProteinAnalysis(sequence)
        return analysis.molecular_weight()
    except Exception:
        # 平均アミノ酸分子量による概算
        return len(sequence) * 110.0


# ── CDR Annotation ──────────────────────────────────────────────────

def annotate_cdr_regions(
    vh_sequence: str,
    scheme: str = "imgt",
) -> Optional[dict]:
    """
    ANARCIを用いてCDR領域をアノテーションする。

    Parameters
    ----------
    vh_sequence : str
        VH鎖のアミノ酸配列
    scheme : str
        ナンバリングスキーム（'imgt', 'chothia', 'kabat'）

    Returns
    -------
    Optional[dict]
        CDR1, CDR2, CDR3 の配列とIMGT位置情報を含む辞書。
        ANARCI利用不可時はNone。
    """
    try:
        from anarci import anarci as run_anarci

        results = run_anarci(
            [("query", vh_sequence)],
            scheme=scheme,
            output=False,
        )

        if results is None or results[0] is None:
            logger.warning("ANARCI returned no results for the input sequence")
            return None

        numbering = results[0][0][0]  # 最初の配列の最初のドメイン

        # IMGT定義によるCDR位置
        cdr_ranges = {
            "imgt": {
                "CDR1": (27, 38),
                "CDR2": (56, 65),
                "CDR3": (105, 117),
            },
            "chothia": {
                "CDR1": (26, 32),
                "CDR2": (52, 56),
                "CDR3": (95, 102),
            },
            "kabat": {
                "CDR1": (31, 35),
                "CDR2": (50, 65),
                "CDR3": (95, 102),
            },
        }

        ranges = cdr_ranges.get(scheme, cdr_ranges["imgt"])
        cdrs = {}
        for cdr_name, (start, end) in ranges.items():
            residues = [
                aa for (pos, insertion), aa in numbering
                if start <= pos <= end and aa != "-"
            ]
            cdrs[cdr_name] = "".join(residues)

        return cdrs

    except ImportError:
        logger.warning("ANARCI is not installed; CDR annotation unavailable")
        return None
    except Exception as e:
        logger.error("CDR annotation failed: %s", e)
        return None


# ── Glycan Utilities ────────────────────────────────────────────────

def generate_glycan_conformers(
    smiles: str,
    n_conformers: int = 20,
    random_seed: int = 42,
) -> list:
    """
    RDKitのETKDG法で糖鎖コンフォメーションアンサンブルを生成する。

    Parameters
    ----------
    smiles : str
        糖鎖のSMILES表記
    n_conformers : int
        生成するコンフォーマーの数
    random_seed : int
        乱数シード

    Returns
    -------
    list
        RDKit MolオブジェクトのコンフォーマーID一覧
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        mol = Chem.AddHs(mol)

        params = AllChem.ETKDGv3()
        params.randomSeed = random_seed
        params.numThreads = 0  # 自動検出
        params.pruneRmsThresh = 0.5  # 類似コンフォーマーの除去閾値

        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conformers, params=params)

        if len(conf_ids) == 0:
            logger.warning("No conformers generated; falling back to single conformer")
            AllChem.EmbedMolecule(mol, params)
            conf_ids = [0]

        # MMFF94力場で最適化
        for conf_id in conf_ids:
            try:
                AllChem.MMFFOptimizeMolecule(mol, confId=conf_id, maxIters=500)
            except Exception:
                pass

        logger.info("Generated %d conformers for glycan", len(conf_ids))
        return list(conf_ids)

    except ImportError:
        logger.error("RDKit is required for conformer generation")
        return []
