import matplotlib.pyplot as plt
import pandas as pd
from fpdf import FPDF

class RankingReport(FPDF):
    def header(self):
        self.set_font('Arial', 'B', 16)
        self.cell(0, 10, 'Antibody Candidate Ranking Report', 0, 1, 'C')
        self.ln(5)

def generate_pareto_chart(df, output_img):
    plt.figure(figsize=(8, 6))
    plt.scatter(df['calibrated_score'], df['hydrophobicity'], alpha=0.6, c='blue', s=100)
    plt.title('Affinity vs Hydrophobicity')
    plt.xlabel('Binding Affinity (kcal/mol)')
    plt.ylabel('Hydrophobicity Index')
    plt.grid(True)
    plt.savefig(output_img)