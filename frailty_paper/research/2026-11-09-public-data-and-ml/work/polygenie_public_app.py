# =========================
# Imports
# =========================
import os
import sys
import logging
from pathlib import Path

import dash
from dash import Dash, html, dcc, Input, Output, State, dash_table
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import seaborn as sns
import statsmodels.api as sm
import matplotlib.colors as mcolors

# Ensure project root is importable
root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(root_dir)

from sqlitedb.db_handler import DBHandler

# =========================
# Logging configuration
# =========================

log_level = os.getenv('LOG_LEVEL') or ('DEBUG' if os.getenv('VERBOSE') else 'INFO')
logging.basicConfig(level=getattr(logging, log_level))
log = logging.getLogger(__name__)

# =========================
# Constants
# =========================

ALPHA = 0.05
LOESS_FRAC = 0.3

PHEWAS_HEIGHT = 600
PREVALENCE_HEIGHT = 300

Y_CAP_PERCENTILE = 99

# =========================
# Initialize DB and App
# =========================

# DB
database_file = 'db/polygenie.db'
db = DBHandler(database_file)

# App
app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY])
server = app.server

app.title = "PolyGenie"

# =========================
# Palette and Controls
# =========================

palette = sns.color_palette('colorblind')
colorblind_palette_hex = [mcolors.to_hex(c) for c in palette]

# Controls
gwas_names = db.get_gwas_names()
disease_options = [{'label': d, 'value': d} for d in gwas_names]

reference_options = [
    {'label': 'Bottom', 'value': 'low'},
    {'label': 'Bottom and Middle', 'value': 'rest'},
]

# Default division options (will be updated based on PRS when possible)
division_options = [{'label': '4 (quartiles)', 'value': '4'}, {'label': '10 (deciles)', 'value': '10'}]

significance_options = [
    {'label': 'Bonferroni (domain-wise)', 'value': 'bonferroni_domain'},
    {'label': 'Bonferroni (global)', 'value': 'bonferroni_global'},
    {'label': 'FDR (global)', 'value': 'fdr_global'},
]

# =========================
# Targe classes for tabs
# =========================

# Build dynamic tabs from the DB's target classes
_tc_df = db.get_target_classes()

target_classes = (
    _tc_df['target_class'].astype(str).tolist()
    if not _tc_df.empty
    else ['Phecodes','ICD_codes','Metabolites','Other Variables']
)

tabs_children = [dcc.Tab(label=tc, value=tc) for tc in target_classes]
default_tab = target_classes[0] if target_classes else 'ICD_codes'

# =========================
# Helper functions
# =========================

def empty_figure(message: str):
    """Return a placeholder empty figure."""
    return px.scatter(title=message)


def add_bonferroni_lines(fig, element_count, alpha=ALPHA):
    """Add Bonferroni significance lines and annotations."""
    if not element_count:
        return None, None

    pos = -np.log10(alpha / element_count)
    neg = np.log10(alpha / element_count)

    for y in (pos, neg):
        fig.add_shape(
            type='line',
            x0=0, x1=1, xref='paper',
            y0=y, y1=y, yref='y',
            line=dict(color='red', width=2, dash='dash')
        )
        fig.add_annotation(
            xref='paper',
            x=1,
            y=y,
            text=f"{y:.3f}",
            showarrow=False,
            font=dict(color='red', size=12),
            align='right',
            yshift=10
        )

    return pos, neg


def add_loess(fig, df, x, y, label, visible):
    """Add a LOESS-smoothed line to a Plotly figure."""
    if df.empty:
        return

    loess = sm.nonparametric.lowess(df[y], df[x], frac=LOESS_FRAC)

    fig.add_trace(
        go.Scatter(
            x=loess[:, 0],
            y=loess[:, 1],
            mode='lines',
            name=f'LOESS {label}',
            line=dict(width=3),
            visible=visible,
            showlegend=True
        )
    )


def compute_significance_cutoff(pvalues, method='bonferroni', alpha=ALPHA, m=None):
    """Return a p-value cutoff for Bonferroni or FDR correction."""
    pvalues = np.asarray(pvalues, dtype=float)
    pvalues = pvalues[np.isfinite(pvalues)]
    if m is None:
        m = len(pvalues)

    if m == 0:
        return None

    if method == 'fdr':
        sorted_p = np.sort(pvalues)
        thresholds = np.arange(1, m + 1) * alpha / m
        valid = sorted_p <= thresholds
        if not np.any(valid):
            return alpha / m
        return float(sorted_p[valid].max())

    return alpha / m


def get_significance_cutoffs(filtered_data, selection):
    """Compute a p-value cutoff from filtered data based on the selected threshold strategy."""
    if filtered_data.empty or 'P' not in filtered_data.columns:
        return []

    if selection == 'bonferroni_domain':
        domain_count = filtered_data['domain'].nunique() if 'domain' in filtered_data.columns else len(filtered_data)
        cutoff = compute_significance_cutoff(filtered_data['P'].values, method='bonferroni', m=domain_count)
        return [('Bonferroni (domain-wise)', cutoff, domain_count)]

    if selection == 'bonferroni_global':
        cutoff = compute_significance_cutoff(filtered_data['P'].values, method='bonferroni', m=len(filtered_data))
        return [('Bonferroni (global)', cutoff, len(filtered_data))]

    if selection == 'fdr_global':
        cutoff = compute_significance_cutoff(filtered_data['P'].values, method='fdr', m=len(filtered_data))
        return [('FDR (global)', cutoff, len(filtered_data))]

    return []


def get_target_class_for_phewas(prs_code):
    """Choose a sensible default target class for a PRS."""
    df = db.get_target_classes_for_prs(prs_code)
    if df.empty:
        return 'Phecodes'

    classes = df['target_class'].tolist()
    for c in classes:
        if 'phecode' in c.lower():
            return c
    for c in classes:
        if 'icd' in c.lower():
            return c
    return classes[0]

# =========================
# App Layout
# =========================

# Layout: Only Phewas (Phecodes) and prevalence + table
app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.Div([
                html.H1([
    			html.Img(src="/assets/PolyGenie.png", height="80px", style={"verticalAlign": "middle"})
		]),
                html.P("Unearthing Genetic Links with Polygenic Scores", className="subheading"),
                html.Hr(className="short-sep-line"),
            ], className="prs-info-box"),
            html.P('Select PRS'),
            dcc.Dropdown(
                id='disease-dropdown',
                options=disease_options,
                value=(gwas_names[0] if gwas_names else None)
            ),

            html.P('Reference'),
            dcc.Dropdown(
                id='reference-dropdown',
                options=reference_options,
                value='low'
            ),

            html.P('Division (n groups)'),
            dcc.Dropdown(
                id='division-dropdown',
                options=division_options,
                value='10'
            ),

            html.P('Significance threshold'),
            dcc.Dropdown(
                id='significance-selection-dropdown',
                options=significance_options,
                value='bonferroni_domain'
            ),
        ], width=3),

        dbc.Col([
            dcc.Tabs(id='tabs', value=default_tab, children=tabs_children),

            dcc.Graph(
                id='phewas-graph',
                style={'height': f'{PHEWAS_HEIGHT}px'},
                config={'displayModeBar': False}
            ),

            html.Hr(),
            html.H4('Prevalence by percentile'),

            dcc.Graph(
                id='prevalence-graph',
                style={'height': f'{PREVALENCE_HEIGHT}px'},
                config={'displayModeBar': False}
            ),

            html.Hr(),
            html.H4('Top hits'),

            dash_table.DataTable(
                id='prs-table',
                page_size=10,
                style_table={'overflowX': 'auto'},
                style_cell={'textAlign': 'left', 'padding': '5px'},
                style_header={'fontWeight': 'bold'},
            ),
        ], width=9),
    ]),

    dcc.Store(id='clicked-data-store')
], fluid=True)

# =========================
# Callbacks
# =========================


@app.callback(
    Output('division-dropdown', 'options'), Output('division-dropdown', 'value'),
    Input('disease-dropdown', 'value')
)
def update_divisions(prs_display_name):
    if not prs_display_name:
        return division_options, '10'
    prs = db.get_gwas_code_from_name(prs_display_name)
    groups = db.get_prs_n_groups(prs)
    if not groups:
        return division_options, '10'
    opts = [{'label': f"{n} bins", 'value': str(n)} for n in sorted(groups)]
    return opts, str(sorted(groups)[-1])


@app.callback(
    Output('reference-dropdown','options'), Output('reference-dropdown','value'),
    Input('disease-dropdown','value')
)
def update_reference_options(prs_display_name):
    """Update Reference dropdown to include 'Bottom and Middle' when PRS runs include intermediates."""
    base = [{'label': 'Bottom', 'value': 'low'}]
    if not prs_display_name:
        return base, 'low'
    prs = db.get_gwas_code_from_name(prs_display_name)
    try:
        if db.get_prs_include_intermediates(prs):
            opts = base + [{'label': 'Bottom and Middle', 'value': 'rest'}]
            return opts, 'rest'
    except Exception:
        # If something goes wrong, fall back to base option
        pass
    return base, 'low'

@app.callback(
    Output('phewas-graph', 'figure'),
    [Input('disease-dropdown', 'value'), Input('reference-dropdown', 'value'), Input('division-dropdown', 'value'), Input('tabs', 'value'), Input('significance-selection-dropdown', 'value')]
)
def update_phewas(prs_display_name, reference_value, division_value, tab, significance_selection):
    """Use the plotting code from `old/app.py` but allow tab-based selection of target class."""
    if not prs_display_name:
        return px.scatter(title='Select PRS')
    prs = db.get_gwas_code_from_name(prs_display_name)
    # Determine target class preference: tab overrides auto-detection; tabs carry exact DB target_class names
    if tab:
        target_type = tab
    else:
        target_type = get_target_class_for_phewas(prs)
    log.debug("Using target_type '%s' for PRS %s (tab=%s)", target_type, prs, tab)
    # Fetch data (allow for target_type being a tuple)
    if isinstance(target_type, tuple):
        df = pd.concat([
            db.get_correlations(prs, reference_value, division_value, target_type[0]),
            db.get_correlations(prs, reference_value, division_value, target_type[1])
        ], ignore_index=True)
    else:
        df = db.get_correlations(prs, reference_value, division_value, target_type)
    if df.empty:
        return px.scatter(title='No data available')

    # same processing as full app
    df = df.dropna(subset=['P'])
    df = df.sort_values(by=['class', 'domain', 'logpxdir'], ascending=[True, True, True]).reset_index(drop=True)

    # map each row to an integer x position so we can place domain separators precisely
    df['xpos'] = df.index.astype(float)
    x = 'xpos'
    y = 'logpxdir'

    unique_categories = df['domain'].unique() if 'domain' in df.columns else []
    color_map = {category: colorblind_palette_hex[i % len(colorblind_palette_hex)] for i, category in enumerate(unique_categories)}

    # Add lowercase code for hover
    if 'Code' in df.columns and 'code' not in df.columns:
        df['code'] = df['Code']

    # Build scatter on numeric x so domain separators (vertical lines) can be placed at fractional positions
    fig = px.scatter(df, 
                     x=x, 
                     y=y, 
                     color='domain', 
                     color_discrete_map=color_map, 
                     template='plotly_white',
                     custom_data=['Code','description','domain'],
                     hover_data={'code': True, 
                                 'description': True,
                                 'class': True,
                                 'odds_ratio': ':.2f',
                                 'beta': ':.2f', 
                                 'logpxdir': False, 
                                 'P': ':.2e',
                                 'domain': True,
                                 'xpos': False,
                                 })


    # Threshold lines and annotations
    sig_cutoffs = get_significance_cutoffs(df, significance_selection)
    if sig_cutoffs:
        label, cutoff, count = sig_cutoffs[0]
        if cutoff is not None and cutoff > 0:
            positive_sig = -np.log10(cutoff)
            negative_sig = np.log10(cutoff)
            for y in (positive_sig, negative_sig):
                fig.add_shape(
                    type='line',
                    x0=0, x1=1, xref='paper',
                    y0=y, y1=y, yref='y',
                    line=dict(color='red', width=2, dash='dash')
                )
            fig.add_annotation(
                xref='paper',
                x=1,
                y=positive_sig,
                text=f"{label} ({count} tests): {cutoff:.2e}",
                showarrow=False,
                font=dict(color='red', size=12),
                align='right',
                yshift=10
            )

            outliers = df[(df['logpxdir'] < negative_sig) | (df['logpxdir'] > positive_sig)]
            for index, row in outliers.iterrows():
                xpos_val = row['xpos'] if 'xpos' in row else row[x]
                desc = row.get('description', '')
                fig.add_annotation(x=xpos_val, y=row['logpxdir'], text=f"{desc}", showarrow=False, font=dict(color='black', size=12), align='center', yshift=10)

    # Cap y-axis
    #yvals = df['logpxdir'].abs().dropna()
    #if len(yvals) > 0:
    #    ycap = max(1.0, float(np.nanpercentile(yvals, Y_CAP_PERCENTILE)))
    #    y_lim = ycap * 1.2
    #    fig.update_yaxes(range=[-y_lim, y_lim])

    fig.update_layout(
            xaxis_title='Diseases',
            yaxis_title='log10(p-value) × Effect size direction',
            xaxis_visible=False,
            xaxis_showticklabels=False,
            showlegend=False,
            autosize=False,
            height=PHEWAS_HEIGHT,
            margin=dict(l=40, r=40, t=40, b=80)
        )

    # baseline
    fig.add_shape(type='line', x0=0, x1=1, xref='paper', y0=0, y1=0, yref='y', line=dict(color='black', width=1))
    fig.update_xaxes(range=[-0.5, len(df) - 0.5])
    return fig


@app.callback(Output('clicked-data-store', 'data'), Input('phewas-graph', 'clickData'))
def store_clicked(clickData):
    if not clickData or 'points' not in clickData or len(clickData['points']) == 0:
        return None
    pt = clickData['points'][0]
    custom = pt.get('customdata') or []
    code = custom[0] if len(custom) > 0 else pt.get('x')
    desc = custom[1] if len(custom) > 1 else (pt.get('text') or code)
    try:
        # Ensure serializable strings
        code = str(code)
    except Exception:
        code = None
    try:
        desc = str(desc)
    except Exception:
        desc = None
    return {'code': code, 'description': desc}


@app.callback(
    Output('prevalence-graph', 'figure'),
    [
        Input('clicked-data-store', 'data'),
        Input('disease-dropdown', 'value')
    ]
)
def update_prevalence(clicked, disease_value):

    if not clicked or not clicked.get('code'):
        return px.scatter(title='Click a point on the main plot to show prevalence by percentile')

    prs_code = db.get_gwas_code_from_name(disease_value)
    df = db.get_prevalences(prs_code, clicked['code'])
    ttype = db.get_target_type(clicked['code'])

    if df.empty:
        return px.scatter(title='No prevalence data available for selected target')

    # Map PRS column → display label
    df['Sex'] = df['prs_column'].map({
        'PRS_agg': 'All',
        'PRS_male': 'Male',
        'PRS_female': 'Female'
    })

    df = df.dropna(subset=['Sex'])

    if ttype == 'continuous':
        ylabel = 'Mean value'
    else:
        ylabel = 'Prevalence'

    n_groups = df['percentile'].max() +1
    if n_groups == 100:
        xlabel = "Percentile"
    elif n_groups == 10:
        xlabel = "Decile"
    elif n_groups == 4:
        xlabel = "Quartile"
    else:
        xlabel = f"{n_groups} risk score groups"

    fig = px.scatter(
        df,
        x='percentile',
        y='prevalence',
        color='Sex',
        opacity=0.4,
        template='plotly_white'
    )

    # Hide Male/Female points by default
    fig.for_each_trace(
        lambda t: t.update(visible=True if t.name == 'All' else 'legendonly')
    )

    # Apply LOESS
    add_loess(fig, df[df['Sex'] == 'All'], 'percentile', 'prevalence', 'All', True)
    add_loess(fig, df[df['Sex'] == 'Female'], 'percentile', 'prevalence', 'Female', 'legendonly')
    add_loess(fig, df[df['Sex'] == 'Male'], 'percentile', 'prevalence', 'Male', 'legendonly')

    fig.update_layout(
        xaxis_title=f'{xlabel} ({disease_value})',
        yaxis_title=ylabel,
        height=PREVALENCE_HEIGHT,
        margin=dict(l=40, r=40, t=20, b=40)
    )

    return fig


@app.callback(Output('prs-table', 'data'), [Input('disease-dropdown', 'value'), Input('reference-dropdown', 'value'), Input('division-dropdown', 'value'), Input('tabs', 'value')])
def update_table(disease_value, reference_value, division_value, tab):
    if not disease_value:
        return []
    prs = db.get_gwas_code_from_name(disease_value)
    # choose target_type directly from tab if supplied
    if tab:
        target_type = tab
    else:
        target_type = get_target_class_for_phewas(prs)
    log.debug("Table query using target_type '%s' for PRS %s (tab=%s)", target_type, prs, tab)
    df = db.get_correlations(prs, reference_value, division_value, target_type)
    if df.empty:
        return []
    df = df.dropna(subset=['P'])
    numeric_columns = df.select_dtypes(include=['number']).columns
    df[numeric_columns] = df[numeric_columns].round(6)
    df = df.replace({pd.NA: None})
    table_df = df[['Code', 'description', 'class', 'beta', 'odds_ratio', 'P']]
    table_df.sort_values(by='P', inplace=True)
    return table_df.to_dict('records')


if __name__ == '__main__':
    app.run(debug=True)
