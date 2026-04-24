import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__))))
from utils import read_config, savefigs

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--path_sims',  required=True)
parser.add_argument('-t', '--path_stats', required=True)
parser.add_argument('-o', '--path_out',   required=True)
args = parser.parse_args()

config = read_config()
METHOD_NAMES = config['method_names']
COLORS       = config['colors']['nets']
ORDERED_METHODS = [
    'celloracle', 'dictys', 'figr', 'granie', 'linger', 'linger_baseline', 'pando',
    'collectri', 'dorothea', 'random', 'scenic', 'scenicplus',
]

# load & clean
sims  = pd.read_csv(args.path_sims)
stats = pd.read_csv(args.path_stats)

for col in ('name_a', 'name_b'):
    sims[col] = sims[col].str.split('.').str[0].str.replace('o_', '', regex=False)
stats['name'] = stats['name'].str.split('.').str[0].str.replace('o_', '', regex=False)

present = [m for m in ORDERED_METHODS if m in stats['name'].values]
pretty  = [METHOD_NAMES.get(m, m) for m in present]


def make_sym_matrix(df, value_col):
    pivot = df.pivot_table(index='name_a', columns='name_b', values=value_col, fill_value=0)
    all_names = pd.unique(df[['name_a', 'name_b']].values.ravel())
    pivot = pivot.reindex(index=all_names, columns=all_names, fill_value=0)
    mat = pivot + pivot.T
    np.fill_diagonal(mat.values, 1.0)
    return mat.loc[present, present]


def heatmap_ax(ax, mat, title, show_yticklabels=True, show_cbar=False):
    sns.heatmap(
        mat, ax=ax,
        cmap='Purples', vmin=0, vmax=1,
        xticklabels=pretty,
        yticklabels=pretty if show_yticklabels else False,
        linewidths=0.3, linecolor='white',
        cbar=show_cbar,
        cbar_kws={'shrink': 0.6, 'label': 'Overlap\nCoefficient'} if show_cbar else {},
        square=True,
    )
    ax.set_title(title, fontsize=9)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=7)
    if show_yticklabels:
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=7)
    ax.set_xlabel('')
    ax.set_ylabel('')


def barplot_ax(ax, col, xlabel):
    sub = stats[stats['name'].isin(present)].copy()
    sub = sub.set_index('name').loc[present].reset_index()
    sub['pretty'] = sub['name'].map(METHOD_NAMES).fillna(sub['name'])
    sub['color']  = sub['name'].map(COLORS).fillna('gray')
    ax.barh(sub['pretty'], sub[col], color=sub['color'], edgecolor='white', height=0.7)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.tick_params(axis='both', labelsize=7)
    ax.invert_yaxis()


# --- panel a: bar charts (all methods, no separation) -----------------------
stat_cols = [
    ('n_tfs',    'Number TFs'),
    ('n_edges',  'Number Edges'),
    ('n_targets','Number Genes'),
    ('odegree',  'Regulon size'),
]
n = len(stat_cols)
h = max(3, len(present) * 0.35)
fig_a, axes = plt.subplots(1, n, figsize=(n * 2.2 + 2.5, h))
for i, (ax, (col, label)) in enumerate(zip(axes, stat_cols)):
    barplot_ax(ax, col, label)
    if i > 0:
        ax.tick_params(labelleft=False)
fig_a.tight_layout(w_pad=0.3)

# --- panel c: TFs / CREs / Genes together, then Edges alone ----------------
sz = max(3, len(present) * 0.35)

# TFs, CREs, Genes side by side — shared colorbar on the right
fig_c1, axes_c1 = plt.subplots(1, 3, figsize=(sz * 3 + 1.2, sz))
for i, (ax, col, title) in enumerate(zip(
    axes_c1,
    ['tf_oc', 'cre_oc', 'target_oc'],
    ['TFs', 'CREs', 'Genes'],
)):
    show_y = (i == 0)
    heatmap_ax(ax, make_sym_matrix(sims, col), title, show_yticklabels=show_y, show_cbar=False)

# shared colorbar
import matplotlib as mpl
cbar_ax = fig_c1.add_axes([0.92, 0.15, 0.015, 0.7])
norm = mpl.colors.Normalize(vmin=0, vmax=1)
sm   = mpl.cm.ScalarMappable(cmap='Purples', norm=norm)
fig_c1.colorbar(sm, cax=cbar_ax, label='Overlap\nCoefficient')
fig_c1.subplots_adjust(right=0.90, wspace=0.05)

# Edges alone
fig_c2, ax_e = plt.subplots(1, 1, figsize=(sz + 1.2, sz))
heatmap_ax(ax_e, make_sym_matrix(sims, 'edge_oc'), 'Edges', show_yticklabels=True, show_cbar=True)
fig_c2.tight_layout()

savefigs([fig_a, fig_c1, fig_c2], args.path_out)
