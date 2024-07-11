#!/usr/bin/env python3

#Author: Jared Adolf-Bryfogle

import plotly.figure_factory as ff
import plotly.graph_objects as go
import plotly.express as px
import dash
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output
import jade2.rosetta_jade.ScoreFiles as scores
import jade2.rosetta_jade.score_util as score_util
from jade2.pymol_jade import PyMolScriptWriter
import jade2.basic.path as jp
from jade2.basic.sequence import calculate_mw
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from typing import *
import numpy as np
import pandas
import os,sys, pathlib, re
from collections import defaultdict

def get_directories(inpath):
    """
    Get a list of directories recursively in a path.  Skips hidden directories.
    :param inpath: str
    :rtype: list
    """

    all_dirs = []
    for root, dirs, files in os.walk(inpath):
        all_dirs.extend([root+"/"+d for d in dirs if d[0] != '.'])
    return all_dirs

def add_annotations(fig, x_pos, y_pos, angle, text):
    fig.update_layout(
        # keep the original annotations and add a list of new annotations:
        annotations = list(fig.layout.annotations) +
                      [go.layout.Annotation(
                          x=x_pos,
                          y=y_pos,
                          font=dict(
                              size=14
                          ),
                          showarrow=False,
                          text=text,
                          textangle=angle,
                          xref="paper",
                          yref="paper"
                      )
                      ]
    )
    return fig

def fix_facetting(fig, x, y, facet_row, facet_col):
    #Fix facetting: https://stackoverflow.com/questions/58167028/single-axis-caption-in-plotly-express-facet-plot

    if not facet_row and not facet_col:
        return fig

    o = -.07 #Origin 1D
    m = .5   #Middle

    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    if facet_row and facet_col:
        for axis in fig.layout:
            if type(fig.layout[axis]) == go.layout.YAxis or type(fig.layout[axis]) == go.layout.XAxis :
                fig.layout[axis].title.text = ''

        fig = add_annotations(fig, o, m, -90, y)
        fig =  add_annotations(fig, m, o, 0, x)
        return fig

    elif facet_row and y:
        for axis in fig.layout:
            if type(fig.layout[axis]) == go.layout.YAxis:
                fig.layout[axis].title.text = ''
        return add_annotations(fig, o, m, -90, y)

    elif facet_col and x:
        for axis in fig.layout:
            if type(fig.layout[axis]) == go.layout.XAxis:
                fig.layout[axis].title.text = ''
        return add_annotations(fig, m, o, 0, x)
    else:
        return fig

    #Add annotation:



def change_font_dec(func):
    def font_up(*args, **kwargs):
        fig = func(*args, **kwargs)
        fig.update_layout(
            font=dict(
                size=15,
            )
        )
        return fig
    return font_up

dfs = 15

def change_font(fig, font_size):
    if not font_size:
        fs = dfs
    else:
        fs = font_size

    fig.update_layout(
        font=dict(
            size=fs,
        )
    )

def order_categories(data, categories):
    categories = [x for x in categories if x]
    orders = {c: sorted([str(s) for s in data[c].unique()]) for c in categories if c}
    return orders

class DashPlot:
    dropbox_dict = {
        "scatter": ["x", "y", "color", "facet_col", "facet_row"],
        "kde": ['x', 'groups'],
        "3d_scatter": ["x", "y", "z", "color"],
        "histogram": ["x", "y", "color"],
        "violin":["x", "y", "color", "facet_col", "facet_row"],
        "box":["x", "y", "color", "facet_col", "facet_row"]
    }

    extra_columns = {
        "kde": {
            "show_hist":['True', 'False'],
            "show_curve":['True','False'],
            "show_rug":['True', 'False']
        },
        "histogram": {
            "histfunc":['count', 'sum', 'avg', 'min','max'],
            "histnorm":['percent', 'probability', 'density','probability density'],
            "marginal": ['rug', 'box', 'violin', 'histogram']
        },
        "scatter": {
            "marginal_x": ['rug', 'box', 'violin', 'histogram'],
            "marginal_y": ['rug', 'box', 'violin', 'histogram'],
            "size": [i for i in np.arange(3.0, 6.25, .25)],
            "trendline ": ['Off','ols','lowess']
        },
        "3d_scatter":{
            "size": [i for i in np.arange(0, 3.25, .25)]
        },
        "violin":{
            "mode":['group', 'overlay'],
            "points":['outliers', 'suspectedoutliers','all','False'],
            "box" :['True', 'False']
        },
        "box": {
            "mode":['group', 'overlay'],
            "points":['outliers', 'suspectedoutliers','all','False'],
            "notched": ['True', 'False']
        }

    }

    style_options ={
        "style" : ["plotly", "plotly_white", "plotly_dark", "ggplot2", "seaborn", "simple_white", "none"],
        "colorscale" : px.colors.named_colorscales(),

        "font_size": [12, 13, 14, 15, 16, 17, 18, 19, 20]
    }

    outdir = "analysis_cache"

    def __init__(self, data, hover_name="decoy_path", outdir="analysis_cache", pymol_cutoff=100, skip_box_copy = False, combined_copy=False):
        self.data = data
        self.hover_name = hover_name
        self.outdir = outdir
        self.box_selection_pymol_cutoff = pymol_cutoff
        self.skip_box_copy = skip_box_copy

        self.figure_dict = {
            "scatter": self.make_scatter,
            "histogram": self.make_histogram,
            "kde": self.make_kde,
            "3d_scatter": self.make_3d_scatter,
            "violin": self.make_violin,
            "box": self.make_box
        }

        self.combined_copy = combined_copy

    def get_all_plot_inputs(self, type):
        return self.dropbox_dict[type] + [d for d in self.extra_columns[type]]

    def get_all_plot_inputs2(self):
        inputs = []
        for type in self.dropbox_dict:
            inputs.extend(self.dropbox_dict[type])
            if type in self.extra_columns:
                inputs.extend(self.extra_columns[type])
        return inputs


    # New Figures go here
    def make_scatter(self,x, y, color, facet_col, facet_row, marginal_x, marginal_y, size, trendline, style, colorscale, fs):
        if not trendline or trendline == "Off":
            trendline = None

        #Try/except for the statsmodels package that doesn't seem to be installed automatically.

        categories = [color, facet_row, facet_col]

        #print(sorted(self.data['yield'].unique()))

        try:
            fig =  px.scatter(
                self.data,
                x=x,
                y=y,
                color=color,
                facet_col=facet_col,
                facet_row=facet_row,
                height=700,
                hover_name=self.hover_name,
                marginal_x=marginal_x,
                marginal_y=marginal_y,
                trendline=trendline,
                color_continuous_scale=colorscale,
                color_discrete_sequence=px.colors.sequential.Jet,
                category_orders = order_categories(self.data, categories)
            )
        except ImportError:
            print("Import error: {0}".format(err))

        if not size:
            size = 5.5
        fig.update_traces(marker=dict(size=size))

        if style:
            fig.update_layout(template=style)

        change_font(fig, fs)
        if facet_col or facet_row:
            #return fig
            return fix_facetting(fig, x, y, facet_row, facet_col)
        else:
            return fig



    def make_histogram(self,x, y, color,histfunc,histnorm,marginal, style, colorscale, fs):
        fig =  px.histogram(
            self.data,
            x=x,
            y=y,
            color=color,
            height=700,
            hover_name=self.hover_name,
            histfunc=histfunc,
            histnorm=histnorm,
            marginal=marginal,
            category_orders = order_categories(self.data, [color])
        )
        change_font(fig, fs)
        if style:
            fig.update_layout(template=style)
        return fig


    def make_3d_scatter(self, x, y, z, color, size, style, colorscale, fs):

        #A bit easier to see as a nice default
        if x and y and z and not color:
            color = z

        fig =  px.scatter_3d(
            self.data,
            x = x,
            y= y,
            z = z,
            color = color,
            height=700,
            hover_name=self.hover_name,
            color_continuous_scale=colorscale,
            category_orders = order_categories(self.data, [color])
        )
        if not size:
            size = 1.75
        fig.update_traces(marker=dict(size=size))
        change_font(fig, fs)
        if style:
            fig.update_layout(template=style)
        return fig

    def make_violin(self, x, y, color, facet_row, facet_col, mode, points, box, style, colorscale, fs):
        if box:
            box = eval(box)
        else:
            box = False

        fig =  px.violin(
            self.data,
            x = x,
            y = y,
            height=700,
            color = color,
            facet_row = facet_row,
            facet_col = facet_col,
            hover_name = self.hover_name,
            violinmode=mode,
            points=points,
            box=box,
            category_orders = order_categories(self.data, [color, facet_row, facet_col])
        )
        fix_facetting(fig, x, y, facet_row, facet_col)
        change_font(fig, fs)
        if style:
            fig.update_layout(template=style)
        return fig

    def make_box(self, x, y, color, facet_row, facet_col, mode, points, notched, style, colorscale, fs):
        if notched:
            notched = eval(notched)

        fig = px.box(
            self.data,
            x = x,
            y = y,
            height = 700,
            color = color,
            facet_row=facet_row,
            facet_col=facet_col,
            hover_name=self.hover_name,
            boxmode=mode,
            points=points,
            notched=notched,
            category_orders = order_categories(self.data, [color, facet_row, facet_col])
        )

        change_font(fig, fs)
        fix_facetting(fig, x, y, facet_row, facet_col)
        if style:
            fig.update_layout(template=style)
        return fig

    def make_kde(self, x, group, hist, curve, rug, style, colorscale, fs):

        if group:
            data = []
            labels = []
            for name, d in self.data.groupby(group):
                labels.append(name)
                data.append(d[x])
        else:
            data = [self.data[x]]
            labels = [x]

        if hist:
            hist = eval(hist)
        else:
            hist = False

        if curve:
            curve = eval(curve)
        else:
            curve = True
        if rug:
            rug = eval(rug)
        else:
            rug = True

        fig = ff.create_distplot(data, labels, show_hist=hist, show_curve=curve, show_rug=rug, bin_size=.25)
        fig.update_layout(height=700)
        change_font(fig, fs)
        if group:
            fig.update_layout(xaxis_title=x)

        if style:
            fig.update_layout(template=style)

        return fig

    ### Extra functions that require class data ###
    def copy_decoy_scatter(self, x_scatter, y_scatter, selection):
        return self.copy_decoy(x_scatter, y_scatter, selection)

    def copy_decoy_3d_scatter(self, x_3d_scatter, y_3d_scatter, z_3d_scatter, selection):
        return self.copy_decoy(x_3d_scatter, y_3d_scatter, selection, z_3d_scatter)

    def copy_decoys_scatter(self, x_scatter, y_scatter, selection):
        return self.copy_decoys(x_scatter, y_scatter, selection)

    def copy_decoys_3d_scatter(self, x_3d_scatter, y_3d_scatter, z_3d_scatter, selection):
        return self.copy_decoys(x_3d_scatter, y_3d_scatter, selection, z_3d_scatter)

    def copy_decoy(self, x, y, selection, z=None):

        #We are using x,y as callbacks here because I do not yet know how to properly deal with the store stuff.
        # Once I figure that out, this will be cleaner.
        if selection is None:
            return {}

        ctx = dash.callback_context

        if not ctx.triggered:
            button_id = 'No clicks yet'
        else:
            button_id = ctx.triggered[0]['prop_id'].split('.')[0]
            if(button_id != "scatter_graph"):
                return {}

        outdir = self.outdir+'/clicked'
        if not os.path.exists(outdir):
            os.mkdir(outdir)


        print(selection)

        outlog = self.outdir+"/click_log.txt"
        pdb_path = selection['points'][0]['hovertext']
        pdb_path_new = pdb_path
        pdb_parent = pathlib.PurePath(pdb_path).parent.name
        if pdb_parent:
            #pdb_path_new=pdb_parent+"_"+os.path.basename(pdb_path)
            pdb_path_new = jp.get_decoy_name(pdb_path)+"_"+pdb_parent+jp.get_decoy_extension(pdb_path)

        x_val = float(selection['points'][0]['x'])
        y_val = float(selection['points'][0]['y'])
        if 'z' in selection['points'][0]:
            z_val = selection['points'][0]['z']
        else:
            z_val = 0

        #line = pdb_path+'\t'+x+"=%s"%.3x_val
        line = "{},{},{:.3f},{},{:.3f},{},{:.3f}\n".format(pdb_path_new,x,x_val,y,y_val,z,z_val)

        if not os.path.exists(outlog):
            LOG = open(outlog, 'w')
            LOG.write("#pdb,x,x_val,y,y_val,z,z_val\n")
        else:
            LOG = open(outlog, 'a')

        LOG.write(line)
        LOG.close()

        print("Copying ", pdb_path, " into "+outdir)
        os.system('cp '+selection['points'][0]['hovertext']+' '+outdir+"/"+pdb_path_new)
        sele = ['hide sticks, polymer.protein']

        colors = ['magenta', 'green', 'blue', 'red', 'orange']
        selections = 0
        for line in open(outdir+"/"+pdb_path_new):
            lineSP = line.strip().split()
            if len(lineSP) > 4 and lineSP[1] == "select":
                sele_name = lineSP[2].strip(',')
                sele.append(" ".join(lineSP[1:]))
                sele.append('center '+sele_name)
                sele.append(f'color {colors[selections]}, '+sele_name)
                sele.append('show sticks, '+sele_name)
                sele.append('zoom '+sele_name+',10')
                sele.append('color atomic, (not elem C)')
                selections+=1

        cmd = ";".join(sele)
        cmd+=";deselect"
        #os.system('pymol '+outdir+"/"+pdb_path_new+f" {native}f"+f" -d '{cmd}' &")
        if native:
            cmd = 'pymol ' + outdir + "/" + pdb_path_new + f" {native}" + f" -d '{cmd}' &"
            print(cmd)
            os.system(cmd)
        else:
            cmd = 'pymol ' + outdir + "/" + pdb_path_new + f" -d '{cmd}' &"
            print(cmd)
            os.system(cmd)

    def copy_decoys(self, x, y, selection, z = None):
        #print("args:",args)
        if selection is None:
            return {}

        ctx = dash.callback_context

        if not ctx.triggered:
            button_id = 'No clicks yet'
        else:
            button_id = ctx.triggered[0]['prop_id'].split('.')[0]
            if button_id != "scatter_graph":
                return {}
        points = selection['points']
        pdbs = []
        outdir = self.outdir+'/grouped'
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        #Get the number of directory
        dirs = get_directories(outdir)
        starting_dir = 0
        for d in dirs:
            if 'group' not in d: continue
            dSP = os.path.basename(d).split("_")
            #print(d)
            print(repr(dSP))
            if int(dSP[-1]) > starting_dir:
                starting_dir = int(dSP[-1])

        outdir = outdir+'/group_'+str(starting_dir+1)
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        xs = []; ys = []; zs = []

        #Get each PDB and numerical values of core xy(z) data.
        for coord in points:
            pdb = coord['hovertext']
            pdbs.append(pdb)
            x_val = float(coord['x'])
            y_val = float(coord['y'])

            if 'z' in coord:
                z_val = coord['z']
            else:
                z_val = 0

            xs.append(x_val)
            ys.append(y_val)
            zs.append(z_val)

        x_m = np.array(xs).mean()
        x_sd = np.array(ys).std()
        y_m = np.array(ys).mean()
        y_sd = np.array(ys).std()
        z_m = np.array(zs).mean()
        z_sd = np.array(zs).std()

        print("Updating group log. ")
        outlog = self.outdir+'/group_log.txt'

        line = "{},{},{:.3f},{:.3f},{},{:.3f},{:.3f},{},{:.3f},{:.3f}\n".format('group_'+str(starting_dir+1),x,x_m,x_sd,y,y_m,y_sd,z,z_m,z_sd)

        if not os.path.exists(outlog):
            LOG = open(outlog, 'w')
            LOG.write("#group,x,x_mean,x_sd,y,y_mean,y_sd,z,z_mean,z_sd\n")
        else:
            LOG = open(outlog, 'a')

        LOG.write(line)
        LOG.close()

        print("Copying ",len(pdbs),"to",outdir)
        model_log = outdir+'/log.txt'

        if not os.path.exists(model_log):
            LOG = open(model_log, 'w')
            LOG.write("#pdb,x,x_val,y,y_val,z,z_val\n")
        else:
            LOG = open(model_log, 'a')

        new_names = {}
        if not self.skip_box_copy:
            print("Copying box/lasso selection into new directory")
        for i, pdb_path in enumerate(pdbs):
            pdb_path_new = pdb_path
            pdb_parent = pathlib.PurePath(pdb_path).parent.name
            if pdb_parent:
                pdb_path_new=jp.get_decoy_name(pdb_path)+"_"+pdb_parent+jp.get_decoy_extension(pdb_path)

            new_names[pdb_path] = pdb_path_new
            x_val = xs[i]
            y_val = ys[i]
            z_val = zs[i]

            line = "{},{},{:.3f},{},{:.3f},{},{:.3f}\n".format(pdb_path_new,x,x_val,y,y_val,z,z_val)
            LOG.write(line)

            if not self.skip_box_copy:
                if not self.combined_copy:
                    os.system('cp '+pdb_path+' '+outdir+"/"+pdb_path_new)
                else:
                    os.system('cp ' + pdb_path + ' ' + outdir)
        LOG.close()

        #NOTE: 7/2024 JAB: With new pymol 3 it seems to call cd for some reason to the grouped directory.
        #So we will explicitly CD back out to the directory we should have been in the first place.
        if len(pdbs) <= self.box_selection_pymol_cutoff:
            #These could also have the same name from different groups.
            #So do 1->N
            scripter = PyMolScriptWriter(outdir)
            scripter.add_line(f'cd {os.getcwd()}') #Fix for PyMol3 wonkiness

            first = True
            first_name = new_names[pdbs[0]]

            if native and os.path.exists(native):
                scripter.add_load_pdb(native, "native")
                
            for i, pdb_path in enumerate(pdbs):
                new_name = jp.get_decoy_name(new_names[pdb_path])
                scripter.add_load_pdb(pdb_path, new_name)
                if not first:
                    scripter.add_align_to(new_name, first_name)
                    scripter.add_line("center "+new_name)
                first = False



            #scripter.add_show("sticks")
            #scripter.add_hide("(hydro)")
            scripter.add_line('import os')
            s = 'cat '+ model_log
            scripter.add_line("os.system(\'"+s+"\' )")

            s = 'head -n 1 ' +outlog
            scripter.add_line("os.system(\'"+s+"\' )")

            s = 'tail -n 1 ' +outlog
            scripter.add_line("os.system(\'"+s+"\' )")
            scripter.save_script()

            print(f'cwd {os.getcwd()}')
            print("Running"+f'pymol {os.getcwd()}/{outdir}/pml_script.pml &')
            os.system(f'pymol {os.getcwd()}/{outdir}/pml_script.pml &')
            os.system('cat '+model_log)

def create_figure_tab(col_options: List[Dict[Any,str]], label: str, pt: str) -> dbc.Tab:
    """
    Create a tab for a particular figure type. Use the label.
    pt is lowercase figure used in DashPlot.
    col_options are
    col_options: data labels for the main data-plotting options.
    """
    t = dbc.Tab(label=label, children=
    [
        html.Div([
            html.Div(
                [
                    html.P([dcc.Dropdown(id=pt+'_'+d, options=col_options, placeholder=d)])
                    for d in DashPlot.dropbox_dict[pt]
                ],
                style={"width": "25%", "float": "left"},
            ),
            dcc.Graph(id=pt+"_graph", style={"width": "75%", "display": "inline-block", "float": "right"}, config={'toImageButtonOptions':{'scale':4}}),
            html.Div(
                [
                    html.P([dcc.Dropdown(id=pt+'_'+d, options=[dict(label=x, value=x) for x in DashPlot.extra_columns[pt][d]], placeholder=d)])
                    for d in DashPlot.extra_columns[pt]
                ],
                style={"width": "25%", "float": "left"},
            ),
            dcc.Store(id=pt+'_store', storage_type='memory'),
            dcc.Store(id=pt+'_lasso_store', storage_type='memory'),
        ])
    ])
    return t


#https://pierpaolo28.github.io/blog/blog21/
def apply_layout(app, df_columns: List[str]):
    #print("df_columns", df_columns)
    pt = "scatter"
    col_options = [dict(label=x, value=x) for x in df_columns]
    #print(col_options)
    app.layout = html.Div([
        dbc.Row(
            [
                html.H1("Jade Integrated Analysis Environment")
            ], justify="center", align="center", className="h-50"
        ),
        dbc.Row(
            [
                dbc.Col(html.P([dcc.Dropdown(id=d, options=[dict(label=x, value=x) for x in DashPlot.style_options[d]], placeholder=d)]))
                for d in DashPlot.style_options
            ],
            #Had to remove no_gutters=True,
            style={'padding-left':'70%'},
        ),

        html.Div([
            dbc.Tabs([
                create_figure_tab(col_options, "Scatter", 'scatter'),
                create_figure_tab(col_options, "KDE",'kde'),
                create_figure_tab(col_options, "Histogram", 'histogram'),
                create_figure_tab(col_options, "Violin", "violin"),
                create_figure_tab(col_options, "Box", "box"),
                create_figure_tab(col_options, "3D Scatter", '3d_scatter')
            ])
        ])
    ])

def get_options():
    parser = ArgumentParser("This app is to analyze scorefiles as JSON. ", formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument("--scorefile", "-s",
                        help="The scorefile(s).  If multiple scorefiles, it will concat and add a subset field"
                             " The subset field will use the last directory in each path as the name OR can "
                             " be specified with scorefile_path=subset_name in order to name each subset. ",
                        default="test",
                        nargs="*")

    parser.add_argument("--highlights", "-e",
                        help="A TSV/CSV file with a PDB name and experimental metrics columns.  They can be numerical or otherwise. "
                             "Must have header. Used for subsetting or coloring based on results. Does not need every PDB. decoy should be the column with the decoy name, including 0001 or whatever.",
                        )

    parser.add_argument("--debug", '-d',
                        help = "Run in debug mode",
                        default = False,
                        action="store_true")

    parser.add_argument("--outdir", '-o',
                        help = "Output directory for decoys.",
                        default="analysis_cache")

    parser.add_argument("--box_selection_pymol_cutoff", '-c',
                        help = "The cutoff of the number of decoys for box selection (lasso/box) to not open a pymol window. Set to 0 if you do not wish to open one.",
                        default = 150)

    parser.add_argument("--skip_box_copy", '-b',
                        help = "Skip box copy. Useful if just browsing. ",
                        default = False,
                        action = "store_true")

    parser.add_argument("--clear_output_dir_on_load", '-l',
                        help = "Clear the output directory on load? Useful if you are just looking at structures. ",
                        default = False,
                        action="store_true")

    parser.add_argument("--hover_name",
                        default = "decoy_path",
                        help = "Hover name to use.  Change if not analyzing scorefiles.  Otherwise, will be None if not present.")

    parser.add_argument("--sequence_column", '-q',
                        default= "sequence",
                        help = "Sequence column used for calculating MW if present. ")

    parser.add_argument("--native",
                        help = "Load the native pose along with the designs")

    parser.add_argument("--split_copy",
                        help = "When doing box within directories, split by subdirectories in pdb name",
                        action = "store_true",
                        default = False)

    options = parser.parse_args()
    return options



def read_and_drop_unique_columns(scorefiles) ->pandas.DataFrame:
    """
    Here, we get the union of columns, drop any unique columns, and then concat the result
    I couldn't find a nice way to do this off the shelf.
    We also add either the parent name of the file or the named subset for analysis.
    :param scorefiles:
    :return:
    """
    new_dfs = []
    sets = []
    dfs = []
    for s in scorefiles:
        print("Reading",s)
        score_path = ""
        subset_name = ""
        if '=' not in s:
            score_path = s
            subset_name = pathlib.PurePath(s).parent.name
        else:
            s2 = s.split('=')
            score_path = s2[0]
            subset_name = s2[1]

        if score_path.split('.')[-1] == "csv":
            df = score_util.get_dataframe_from_csv(score_path)
        else:
            df = scores.get_dataframe(score_path, set_index=False)

        df['subset'] = subset_name
        dfs.append(df)
        sets.append(set(df.columns))

    cs_union = sets[0].intersection(*sets[1:])
    for df in dfs:
        to_drop = []
        #print("start: ",len(df.columns))
        for c in df.columns:
            if c not in cs_union:
                to_drop.append(c)
        new_df = df.drop(columns=to_drop)
        new_dfs.append(new_df)

    out_df = pandas.concat(new_dfs)
    #print("Final  :", len(out_df.columns))
    #print("Uniques:", out_df['subset'].unique())

    return out_df

if __name__ == "__main__":
    #app = dash.Dash(
    #    __name__, external_stylesheets=["https://codepen.io/chriddyp/pen/bWLwgP.css"]
    #)

    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
    app.title = "Jade IAE"

    options = get_options()

    global native
    native = options.native
    print("NATIVE: ",native)


    hover_name="decoy_path"
    if options.hover_name:
        hover_name = options.hover_name
    if options.scorefile == "test":
        df = px.data.tips()
        hover_name= "day"
    elif len(options.scorefile) > 1:
        df = read_and_drop_unique_columns(options.scorefile)
    else:
        s = options.scorefile[0]
        if '=' in s:
            s2 = s.split('=')
            score_path = s2[0]
            subset_name = s2[1]

        else:
            score_path = s
            subset_name = pathlib.PurePath(s).parent.name

        print("Reading Single Scorefile:", score_path)
        if score_path.split('.')[-1] == "csv":
            df = score_util.get_dataframe_from_csv(score_path)
        else:
            df = scores.get_dataframe(score_path, set_index=False)
        df['subset'] = subset_name

    #Sort columns
    df = df.reindex(sorted(df.columns), axis=1)

    #Used for classifiers:
    if 'accuracy' in df.columns:
        hover_name = 'accuracy'

    if hover_name not in df.columns:
        hover_name = None

    if options.clear_output_dir_on_load:
        if os.path.exists(options.outdir):
            os.system('rm -r '+options.outdir)


    #Drop 2D and 3D columns. Drop pymol selection.
    columns_3d = [x for x in df.columns if x.split("_")[-1] == "3D"]
    df = df.drop(columns_3d, axis=1)

    columns_2d = [x for x in df.columns if x.split("_")[-1] == "2D"]
    df = df.drop(columns_2d, axis=1)

    #Rename 1D columns if present
    new_names = defaultdict()
    for column in df.columns:
        if column.split("_")[-1] == "1D":
            new_names[column] = column.strip("_1D")
    df = df.rename(new_names, axis=1)

    if options.highlights:
        f = options.highlights
        if not os.path.exists(f):
            sys.exit("Hightlight path does not exist")

        print("Reading experimental results: ", f)
        if f.endswith('.tsv'):
            df2 = pandas.read_csv(f, sep='\t')
        else:
            df2 = score_util.get_dataframe_from_json_or_csv(f)
        if not 'decoy' in df2.columns:
            sys.exit("Hightlight file must have a decoy column that should be the same names as the decoys without .pdb")

        df = df.merge(df2, how="outer", on="decoy").fillna(0)

    if options.sequence_column in df.columns:
        df["MW (kDa)"] = df[options.sequence_column].apply(calculate_mw)

        print("Added molecular weight data as MW (kDa)")

    #print(df['decoy'])
    apply_layout(app, df)

    #DO ALL YOUR DF MANIPULATIONS ABOVE THIS LINE!
    p = DashPlot(df, hover_name, options.outdir, int(options.box_selection_pymol_cutoff), options.skip_box_copy, not options.split_copy)
    if not os.path.exists(p.outdir):
        os.mkdir(p.outdir)



    df.to_csv(p.outdir+"/"+"current_data.csv")
    col_options = [dict(label=x, value=x) for x in df.columns]

    #Decorate the given figure type

    for pt in DashPlot.dropbox_dict:
        app.callback(Output(pt+"_graph", "figure"), [Input(pt+"_"+d, "value") for d in p.get_all_plot_inputs(pt)], [Input(d, "value") for d in p.style_options])(p.figure_dict[pt])
        if pt == "scatter":
            app.callback(Output(pt+'_store', 'data'), [Input(pt+"_"+d, "value") for d in ["x", "y"]],Input(pt+"_graph", 'clickData'))(p.copy_decoy_scatter)
            app.callback(Output(pt+'_lasso_store', 'data'), [Input(pt+"_"+d, "value") for d in ["x", "y"]],Input(pt+"_graph", 'selectedData'))(p.copy_decoys_scatter)

        if pt == "3d_scatter":
            app.callback(Output(pt+'_store', 'data'), [Input(pt+"_"+d, "value") for d in ["x", "y","z"]], Input(pt+"_graph", 'clickData'))(p.copy_decoy_3d_scatter)
            app.callback(Output(pt+'_lasso_store', 'data'), [Input(pt+"_"+d, "value") for d in ["x", "y", "z"]],Input(pt+"_graph", 'selectedData'))(p.copy_decoys_3d_scatter)


    try:
        if options.debug:
            app.server.run()
            app.run_server(debug=True, port=8000, host='127.0.0.1')
        else:
            #sys.stdout = open(os.devnull, "w")
            sys.stderr = open(os.devnull, "w")
            app.run_server(debug=False, port=8000, host='127.0.0.1')
    except OSError:
        print("Jade IAE cannot run while previous pymol windows are open or another instance of the app is running.  Please close and restart the app")

    except ImportError as err:
        print("Import error: {0}".format(err))



#Notes:
 #CSS: https://www.w3schools.com/css/css_float_clear.asp
 #Color scales: https://plotly.com/python/builtin-colorscales/ (Maybe upper right)
 #Themes: https://plotly.com/python/templates/ (Maybe upper right)
 ## https://medium.com/plotly/introducing-plotly-py-theming-b644109ac9c7

 ## .iae directory or .iaerc file for personalization - especially color scales, pymol stuff, etc.
 ## Defaults dictionary with options perhaps
 ### colorscale, trendline, pymol stuff, dot width
 #Add click callback to violin/box plots to see outliers, etc.

 #Multiple hosts: https://community.plotly.com/t/open-dash-app-on-specified-address/5596/3
 #Changing the download location - can't find much to help with that, would need to be a new button.
 #Configuring download image: https://plotly.com/python/configuration-options/
 ## https://github.com/plotly/documentation/pull/1215/files
 #Configuring download paths: https://docs.faculty.ai/user-guide/apps/examples/dash_file_upload_download.html




