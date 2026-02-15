🧪 GlycoAntibody Studio
GlycoAntibody Studio は、癌関連糖鎖抗原（TACA）などの糖鎖を標的とした抗体を自動設計・シミュレーションするための計算プラットフォームです。トラスツズマブ（Trastuzumab）の安定したフレームワークを基盤とし、標的糖鎖に対する相補性決定領域（CDR）の最適化、水和を考慮したドッキング、および多目的最適化による候補選択を統合しています。

🔬 科学的背景：糖鎖抗原と抗体の相互作用解析
糖鎖はタンパク質と異なり、高度な柔軟性（Conformational Flexibility）と複雑な水和シェルを持つため、抗体による特異的認識の難易度が高いターゲットです。本アプリでは、以下の科学的知見に基づいた解析手法を採用しています。

1. 糖鎖のコンフォメーション多様性（Ensemble Approach）
糖鎖は溶液中で複数の立体配座を取り得るため、単一の静的な構造（Native state）のみをターゲットとすることは不十分です。本ツールでは、RDKitを用いたETKDG法により糖鎖のアンサンブルを生成し、多様な配座に対してドッキングを行うことで、真の結合ポーズの捕捉率を向上させています。

2. 保存水（Conserved Waters）による架橋作用
糖鎖の表面は水酸基（-OH）が豊富であり、抗体との結合界面には水分子が介在することが多々あります。本アプリでは、高分解能PDB構造から同定された「保存水」を受容体の一部として保持したままドッキングを行うことで、水分子を介した水素結合ネットワークを考慮した高精度な予測を実現しています。

3. 脱水和コストと水素結合の再重み付け（Scoring Calibration）
糖鎖ドッキングの精度を向上させるため、AutoDock Vinaのスコアを補正しています。糖鎖結合において支配的な水素結合の寄与を高め、水溶液中からポケットへ移動する際の脱水和（Desolvation）コストを物理化学的パラメータに基づいて再評価しています。

🚀 主な機能
CDR Grafting & Indels Handling: IMGT番号付けに基づき、トラスツズマブのフレームワークへCDRを移植。ループ長の変化（挿入・欠損）にも対応。

Structural Optimization: OpenMMを用いたエネルギー最小化（AMBER14力場）により、置換に伴う立体衝突を解消。

Hydrated Ensemble Docking: 保存水と糖鎖アンサンブルを用いた高度なシミュレーション。

Pareto Analysis: 結合親和性（Affinity）と疎水性（Hydrophobicity/Aggregation risk）のトレードオフを可視化し、ドラッグライクなリード抗体を特定。

Automated Reporting: 解析結果を比較テーブルとパレート図を含むPDF形式で出力。

🛠 インストールと実行
依存ライブラリのインストール
Bash
conda install -c bioconda anarci
conda install -c conda-forge openmm pdbfixer
pip install -r requirements.txt
アプリの起動
Bash
streamlit run app.py
📈 ワークフロー
Input: ターゲット糖鎖のSMILESと、設計したいCDR3配列のリスト（CSV）をアップロード。

Analysis: バックエンドで「移植→最適化→ドッキング→スコアリング」をバッチ処理。

Selection: パレート図を用いて、親和性と物性のバランスが優れた候補を選択。

Export: 3D構造の確認と、解析レポートのダウンロード。

📅 今後の展望
深層学習を用いたCDR配列の自動生成エンジンの統合。

分子動力学（MD）シミュレーションによる結合自由エネルギーのより精密な見積もり（TI/FEP法など）。