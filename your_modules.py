"""
GlycoAntibody Studio - Core Modules
3D構造生成（ComplexBuilder）および CDR移植エンジン（AntibodyGraftingEngine）
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ── Constants ───────────────────────────────────────────────────────

# 3文字 → 1文字 アミノ酸コード変換
AA_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
AA_1TO3 = {v: k for k, v in AA_3TO1.items()}

VALID_AA = set(AA_1TO3.keys())

# ヘリックス構造パラメータ
HELIX_RADIUS = 2.3          # Å（α-ヘリックス半径）
HELIX_PITCH = 1.5           # Å / residue（ライズ）
HELIX_THETA_STEP_DEG = 100  # °/ residue（回転角）

# バックボーン原子の相対座標オフセット（CA基準）
BACKBONE_ATOMS = (
    (" N  ", -0.5, +1.2, -1.0, "N"),
    (" CA ",  0.0,  0.0,  0.0, "C"),
    (" C  ", +0.5, -1.2, +1.0, "C"),
    (" O  ", +1.5, -1.5, +1.0, "O"),
)


# ── Data Classes ────────────────────────────────────────────────────

@dataclass
class CDRSet:
    """CDR配列セット（重鎖または軽鎖）。"""
    cdr1: str
    cdr2: str
    cdr3: str

    def __post_init__(self):
        for name, seq in [("CDR1", self.cdr1), ("CDR2", self.cdr2), ("CDR3", self.cdr3)]:
            if not seq or not all(c in VALID_AA for c in seq.upper()):
                raise ValueError(f"Invalid amino acid sequence in {name}: '{seq}'")
        self.cdr1 = self.cdr1.upper()
        self.cdr2 = self.cdr2.upper()
        self.cdr3 = self.cdr3.upper()

    def as_list(self) -> list[str]:
        return [self.cdr1, self.cdr2, self.cdr3]


@dataclass
class GraftResult:
    """CDR移植の結果。"""
    heavy_chain: str
    light_chain: str
    heavy_cdrs: CDRSet
    light_cdrs: CDRSet

    @property
    def heavy_length(self) -> int:
        return len(self.heavy_chain)

    @property
    def light_length(self) -> int:
        return len(self.light_chain)

    def summary(self) -> dict:
        return {
            "heavy_chain_length": self.heavy_length,
            "light_chain_length": self.light_length,
            "H-CDR1": self.heavy_cdrs.cdr1,
            "H-CDR2": self.heavy_cdrs.cdr2,
            "H-CDR3": self.heavy_cdrs.cdr3,
            "L-CDR1": self.light_cdrs.cdr1,
            "L-CDR2": self.light_cdrs.cdr2,
            "L-CDR3": self.light_cdrs.cdr3,
        }


# ── ComplexBuilder ──────────────────────────────────────────────────

class ComplexBuilder:
    """
    Tab 1: アプリ内 3D 表示用のヘリックス構造生成。

    抗体VH配列からα-ヘリックス近似の骨格座標を生成し、
    PDB形式の文字列として出力する。リガンド・糖鎖SMILESは
    REMARKレコードに記録する。

    Note
    ----
    本クラスはビジュアライゼーション用の簡易モデルであり、
    物理的に正確なフォールディングは行わない。
    実構造が必要な場合は OpenMM / ESMFold 等を使用すること。
    """

    def __init__(
        self,
        radius: float = HELIX_RADIUS,
        pitch: float = HELIX_PITCH,
        theta_step_deg: float = HELIX_THETA_STEP_DEG,
    ):
        self.radius = radius
        self.pitch = pitch
        self.theta_step = np.radians(theta_step_deg)

    def _generate_coords(
        self,
        sequence: str,
        chain_id: str = "A",
        offset_z: float = 0.0,
        start_atom_idx: int = 1,
    ) -> tuple[str, int]:
        """
        アミノ酸配列からα-ヘリックス近似のATOMレコードを生成する。

        Parameters
        ----------
        sequence : str
            アミノ酸配列（1文字コード）
        chain_id : str
            PDBチェインID
        offset_z : float
            Z軸方向のオフセット（Å）。マルチチェイン時に使用。
        start_atom_idx : int
            開始ATOM番号

        Returns
        -------
        tuple[str, int]
            (PDB ATOMレコード文字列, 次のATOM番号)
        """
        if not sequence:
            raise ValueError("Empty sequence provided")

        lines = []
        atom_idx = start_atom_idx

        for i, aa in enumerate(sequence.upper()):
            resname = AA_1TO3.get(aa, "UNK")
            theta = i * self.theta_step

            # Cα座標（ヘリックス軌道上）
            ca_x = self.radius * np.cos(theta)
            ca_y = self.radius * np.sin(theta)
            ca_z = offset_z + i * self.pitch

            for name, dx, dy, dz, elem in BACKBONE_ATOMS:
                x, y, z = ca_x + dx, ca_y + dy, ca_z + dz
                lines.append(
                    f"ATOM  {atom_idx:>5} {name} {resname:>3} {chain_id}{i + 1:>4}    "
                    f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           {elem}"
                )
                atom_idx += 1

        lines.append(f"TER   {atom_idx:>5}      {resname:>3} {chain_id}{i + 1:>4}")
        atom_idx += 1

        return "\n".join(lines), atom_idx

    def build_pdb(
        self,
        vh_sequence: str,
        ligand_smiles: str = "",
        glycan_smiles: str = "",
        vl_sequence: Optional[str] = None,
    ) -> str:
        """
        抗体-糖鎖複合体のPDB文字列を生成する。

        Parameters
        ----------
        vh_sequence : str
            重鎖可変領域の配列
        ligand_smiles : str
            リガンドSMILES
        glycan_smiles : str
            糖鎖SMILES
        vl_sequence : str, optional
            軽鎖可変領域の配列（指定時はB鎖として追加）

        Returns
        -------
        str
            PDB形式の文字列
        """
        # バリデーション
        _validate_sequence(vh_sequence, "VH")
        if vl_sequence:
            _validate_sequence(vl_sequence, "VL")

        # ヘッダー
        header_lines = [
            "REMARK   1 GlycoAntibody Studio - Model Structure",
            "REMARK   2 NOTE: Alpha-helix approximation for visualization only",
        ]
        if ligand_smiles or glycan_smiles:
            combined = ".".join(filter(None, [ligand_smiles, glycan_smiles]))
            header_lines.append(f"REMARK   3 SMILES: {combined}")

        header = "\n".join(header_lines)

        # 重鎖（Chain A）
        vh_records, next_idx = self._generate_coords(
            vh_sequence, chain_id="A", offset_z=0.0
        )

        # 軽鎖（Chain B）
        vl_records = ""
        if vl_sequence:
            z_offset = len(vh_sequence) * self.pitch + 10.0  # 十分な間隔
            vl_records_str, next_idx = self._generate_coords(
                vl_sequence, chain_id="B", offset_z=z_offset,
                start_atom_idx=next_idx,
            )
            vl_records = "\n" + vl_records_str

        return f"{header}\n{vh_records}{vl_records}\nEND"


# ── AntibodyGraftingEngine ──────────────────────────────────────────

class AntibodyGraftingEngine:
    """
    Tab 2: トラスツズマブのフレームワークへの CDR 移植。

    トラスツズマブ（Herceptin）の安定したフレームワーク領域（FR1-4）を
    スキャフォールドとして使用し、標的糖鎖に対して予測されたCDR配列を
    移植することで、完全な重鎖・軽鎖可変領域配列を生成する。

    References
    ----------
    - Carter P, et al. (1992) PNAS 89:4285-4289 (Trastuzumab humanization)
    - Lefranc MP, et al. (2003) Dev Comp Immunol 27:55-77 (IMGT numbering)
    """

    # トラスツズマブ FR 配列（Kabatナンバリング準拠）
    HEAVY_FR = {
        "FR1": "EVQLVESGGGLVQPGGSLRLSCAAS",
        "FR2": "WVRQAPGKGLEWVA",
        "FR3": "RFTISADTSKNTAYLQMNSLRAEDTAVYYC",
        "FR4": "WGQGTLVTVSS",
    }
    LIGHT_FR = {
        "FR1": "DIQMTQSPSSLSASVGDRVTITC",
        "FR2": "WYQQKPGKAPKLLIY",
        "FR3": "GVPSRFSGSGSGTDFTLTISSLQPEDFATYYC",
        "FR4": "FGQGTKVEIK",
    }

    # 糖鎖タイプ別のデフォルトCDRテンプレート（デモ用）
    _CDR_TEMPLATES = {
        "default": {
            "heavy": CDRSet(cdr1="GFTFSRYT", cdr2="ISSSGGST", cdr3="ARTVRYGMDV"),
            "light": CDRSet(cdr1="QSVSSY",   cdr2="DAS",       cdr3="QQRSSWPFT"),
        },
        "tn_antigen": {
            "heavy": CDRSet(cdr1="GYTFTSYW", cdr2="IYPGSGNT", cdr3="ARDYYGSSY"),
            "light": CDRSet(cdr1="RASQSV",   cdr2="AAS",       cdr3="QQYGSSPWT"),
        },
        "sialyl_tn": {
            "heavy": CDRSet(cdr1="GFTFSSYA", cdr2="ISGSGGST", cdr3="AKDYGSSWY"),
            "light": CDRSet(cdr1="QSVLYS",   cdr2="GAS",       cdr3="QQYGSSPIT"),
        },
    }

    def __init__(self):
        logger.info("AntibodyGraftingEngine initialized (scaffold: Trastuzumab)")

    def predict_cdrs(
        self,
        glycan_smiles: str,
        template: str = "default",
    ) -> tuple[CDRSet, CDRSet]:
        """
        標的糖鎖に対する CDR 配列を予測する。

        Parameters
        ----------
        glycan_smiles : str
            標的糖鎖のSMILES表記
        template : str
            使用するCDRテンプレート名
            ('default', 'tn_antigen', 'sialyl_tn')

        Returns
        -------
        tuple[CDRSet, CDRSet]
            (重鎖CDRセット, 軽鎖CDRセット)

        Note
        ----
        現在はテンプレートベースの予測。将来的には深層学習モデル
        （例: AbLang, IgFold等）との統合を予定。
        """
        if not glycan_smiles or not glycan_smiles.strip():
            raise ValueError("Glycan SMILES is required for CDR prediction")

        entry = self._CDR_TEMPLATES.get(template)
        if entry is None:
            available = ", ".join(self._CDR_TEMPLATES.keys())
            raise ValueError(
                f"Unknown template '{template}'. Available: {available}"
            )

        logger.info(
            "CDR prediction: template='%s', glycan='%s'",
            template, glycan_smiles[:40],
        )
        return entry["heavy"], entry["light"]

    def graft(
        self,
        heavy_cdrs: CDRSet,
        light_cdrs: CDRSet,
    ) -> GraftResult:
        """
        予測 CDR をトラスツズマブFRに移植して完全なFv配列を生成する。

        Parameters
        ----------
        heavy_cdrs : CDRSet
            移植する重鎖CDR配列セット
        light_cdrs : CDRSet
            移植する軽鎖CDR配列セット

        Returns
        -------
        GraftResult
            移植結果（重鎖・軽鎖の完全配列とCDR情報）
        """
        heavy = self._assemble_chain(self.HEAVY_FR, heavy_cdrs)
        light = self._assemble_chain(self.LIGHT_FR, light_cdrs)

        logger.info(
            "Grafting complete: VH=%d aa (CDR3=%s), VL=%d aa (CDR3=%s)",
            len(heavy), heavy_cdrs.cdr3,
            len(light), light_cdrs.cdr3,
        )

        return GraftResult(
            heavy_chain=heavy,
            light_chain=light,
            heavy_cdrs=heavy_cdrs,
            light_cdrs=light_cdrs,
        )

    @staticmethod
    def _assemble_chain(fr: dict[str, str], cdrs: CDRSet) -> str:
        """FR-CDR-FR-CDR-FR-CDR-FR の順に結合する。"""
        return (
            f"{fr['FR1']}{cdrs.cdr1}"
            f"{fr['FR2']}{cdrs.cdr2}"
            f"{fr['FR3']}{cdrs.cdr3}"
            f"{fr['FR4']}"
        )

    def graft_from_smiles(self, glycan_smiles: str, template: str = "default") -> GraftResult:
        """
        SMILES → CDR予測 → 移植 のワンステップ実行。

        Parameters
        ----------
        glycan_smiles : str
            標的糖鎖のSMILES
        template : str
            CDRテンプレート名

        Returns
        -------
        GraftResult
            移植結果
        """
        h_cdrs, l_cdrs = self.predict_cdrs(glycan_smiles, template)
        return self.graft(h_cdrs, l_cdrs)

    def get_available_templates(self) -> list[str]:
        """利用可能なCDRテンプレート名の一覧を返す。"""
        return list(self._CDR_TEMPLATES.keys())


# ── Utility Functions ───────────────────────────────────────────────

def _validate_sequence(sequence: str, label: str = "sequence") -> None:
    """アミノ酸配列のバリデーション。"""
    if not sequence:
        raise ValueError(f"{label}: empty sequence")
    invalid = set(sequence.upper()) - VALID_AA
    if invalid:
        raise ValueError(
            f"{label}: invalid amino acid characters: {', '.join(sorted(invalid))}"
        )
