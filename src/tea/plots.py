# Venn diagrams
from random import sample
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted, venn3
from matplotlib.lines import Line2D
from pathlib import Path
from tea.utils import sort_for_var

# plotly
import plotly.express as px
import plotly.graph_objects as go
config = {
  'toImageButtonOptions': {
    'format': 'png', # one of png, svg, jpeg, webp
    'filename': 'MyImage',
    'scale': 3 # Multiply title/legend/axis/canvas sizes by this factor
  }
}

def plot_venn(
    working_dir, 
    sets, 
    labels, 
    sample_name = '',
    colors = ['orange','deepskyblue','red'], 
    plot_config = {'font.size': 10, 'font.family': 'Arial', 'font.weight': 'normal'}):
    '''
    Plots a venn diagram with the specified sets and labels
    note: only 2 or 3 sets are supported

    inputs:
    - working_dir: directory to put the "Venn" folder in and save the plot
    - sets: list of sets to be plotted
    - labels: list of labels for each set
    - sample_name: sample name to be used as title (default is none)
    - colors: list of colors for each set (default: ['orange','deepskyblue','red'])
    - plot_config: matplotlib plot configuration
    
    '''
    if not len(sets) == len(labels):
        print('length of sets and labels must be the same')

    plt.rcParams.update(plot_config)
    plt.figure(figsize=(5,5)) 
    ax = plt.gca() 
    
    if len(sets) == 2:
        #print(colors[:2])
        v = venn2(subsets = sets, 
            set_labels = ['',''], 
            ax = ax,
            set_colors = colors[:2],
            alpha=0.7)  
        
    elif len(sets) == 3:
        v = venn3(subsets = sets, 
            set_labels = ['','',''], 
            ax = ax,
            set_colors = colors,
            alpha=0.7)
    else:
        print('only 2 or 3 sets are supported')

    legend_elements = [Line2D([0], [0], marker='o', color='w', label=labels[idx],markerfacecolor=colors[idx], markersize=15) for idx in range(len(labels))]

    ax.legend(handles= legend_elements,
            loc='center', bbox_to_anchor=(0.5, -0.05)) 
      
    if not isinstance(working_dir, Path):
        working_dir = Path(working_dir)
    (working_dir).mkdir(exist_ok=True, parents=True)
    plt.savefig(
        working_dir / f"{sample_name}_{'_'.join(labels)}_Venn_diagram.png",
        dpi=300
        )
    plt.show()

def plot_var_sc_mut_prev_histogram(sample_obj, sample_name, **variants_to_highlight):

    ## STILL UNDER DEVELOPMENT

    total_cell_num = sample_obj.dna.shape[0]

    all_vars_mut_prev_array = sample_obj.dna.get_attribute(
        'mut', constraint='row'
    ).sum(axis=0)
    


    # (a) add known bulk somatic variants
    somatic_vars_mut_prev_array = all_vars_mut_prev_array[[i for i in known_bulk_somatic_vars if i in all_vars_mut_prev_array.index]]

    # (b) add known germline variants
    germline_vars_mut_prev_array = all_vars_mut_prev_array[[i for i in known_bulk_normal_vars_list if i in all_vars_mut_prev_array.index]]

    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x = all_vars_mut_prev_array,
        marker_color = '#6F6F6F',
        opacity=0.75,
        name = 'all variants',
        xbins=dict( # bins used for histogram
            start=0,
            end=total_cell_num*1.2,
            size=5
        ),        
        ))
    fig.add_trace(go.Histogram(
        x = somatic_vars_mut_prev_array,
        marker_color = '#FF6666',
        opacity=0.75,
        name = 'somatic variants identified by bulk',
        xbins=dict( # bins used for histogram
            start=0,
            end=total_cell_num*1.2,
            size=5
        ),   
        ))
    fig.add_trace(go.Histogram(
        x = germline_vars_mut_prev_array,
        marker_color = '#66CC00',
        opacity=0.75,
        name = 'germline variants identified by bulk',
        xbins=dict( # bins used for histogram
            start=0,
            end=total_cell_num*1.2,
            size=5
        ),           
        ))

    fig.update_yaxes(type="log")
    
    fig.update_layout(
        title = f'{sample_i} variant single-cell prevalence histogram',
        barmode='overlay', # <------- overlay histograms
    )

    # label total cell number
    fig.add_vline(
        total_cell_num, 
        line_color='green', 
        line_width=2, row=1, col=1, opacity=0.2,
        annotation_text=f"total cell # = {total_cell_num}",
    #     annotation_font_size=10,
    #     annotation_font_color="red",
    )
    count = 0
    y_val = 0.5
    ay_val = -30
    for var_i in known_bulk_somatic_vars:
        if var_i in all_vars_mut_prev_array.index:
            y_val += count * 0.15
            count += 1
            
            fig.add_annotation(x = all_vars_mut_prev_array[var_i], y = y_val,
                    text = str(known_bulk_somatic_vars[var_i]['Hugo_Symbol']) + ' ' + str(known_bulk_somatic_vars[var_i]['HGVSp_Short']),
                    xref = 'x', yref = 'y',
                    showarrow=True,
                    font=dict(
                        family="Courier New, monospace",
                        size=10,
                        color="#ffffff"
                        ),
                    align="center",
                    arrowhead=2,
                    arrowsize=1,
                    arrowwidth=0.5,
                    arrowcolor="#000000",
                    ax=0,
                    ay=ay_val,
                    bordercolor="#c7c7c7",
                    borderwidth=2,
                    borderpad=4,
                    bgcolor="#ff7f0e",
                    opacity=0.8,
                    captureevents = True
                    )
    mut_prev_histograms[sample_i] = fig

def plot_var_sc_mutational_prev(sample_obj, min_mut_prev_of_interest = 0.01,vars_to_highlight=None, sample_name = None):
    '''
    Plot the single-cell level mutation prevalence of the sample

    inputs:
    - sample_obj: missionbio.sample object

    - min_mut_prev_of_interest: minimum mutation prevalence to be plotted
    - vars_to_highlight: list of variable names to highlight
    
    '''
    if sample_name is None:
        # attempt to get sample_name from h5 metadata
        sample_name = sample_obj.dna.metadata['sample_name']

    mut_prev_array = sample_obj.dna.get_attribute(
        'mut', constraint='row'
    ).sum(axis=0)
    total_cell_num = sample_obj.dna.shape[0]

    fig = px.histogram(
        mut_prev_array,
        log_y = True,
        nbins = total_cell_num
    )

    fig.update_layout(width = 1000, height=400,\
                    title_text = f'{sample_name} variant mutational prevalence )single-cells) histogram<br><sup>total: {total_cell_num} cells</sup>',\
                    title_x = 0.5, title_xanchor = 'center', title_yanchor = 'top')
                    
    fig.update_yaxes(
        title = 'log count of number of variants',
        titlefont=dict(family='Arial', color='crimson', size=20)
    )
    fig.update_xaxes(
        title = 'number of cells mutated',
        range=(0, total_cell_num*1.2)
    )

    # label total cell number
    fig.add_vline(
        total_cell_num, 
        line_color='green', 
        line_width=2, row=1, col=1, opacity=0.2,
        annotation_text=f"total cell # = {total_cell_num}",
    #     annotation_font_size=10,
    #     annotation_font_color="red",
    )

    # label minimum threshold for prevalence
    fig.add_vline(
        total_cell_num * min_mut_prev_of_interest, 
        line_color='red', 
        line_width=2, row=1, col=1, opacity=0.2,
        annotation_text=f"mutation prevalence = {min_mut_prev_of_interest * 100}%",
        annotation_font_size=10,
        annotation_font_color="red",
    )

    # label known bulk somatic mutations
    if vars_to_highlight is not None:
        count = 0
        y_val = 0.5
        ay_val = -30 
        for var_i in vars_to_highlight:
            if var_i in mut_prev_array.index:
                y_val += count * 0.15
                count += 1
                
                fig.add_annotation(x = mut_prev_array[var_i], y = y_val,
                        text = str(vars_to_highlight[var_i]['Hugo_Symbol']) + ' ' + str(vars_to_highlight[var_i]['HGVSp_Short']),
                        xref = 'x', yref = 'y',
                        showarrow=True,
                        font=dict(
                            family="Courier New, monospace",
                            size=10,
                            color="#ffffff"
                            ),
                        align="center",
                        arrowhead=2,
                        arrowsize=1,
                        arrowwidth=0.5,
                        arrowcolor="#000000",
                        ax=0,
                        ay=ay_val,
                        bordercolor="#c7c7c7",
                        borderwidth=2,
                        borderpad=4,
                        bgcolor="#ff7f0e",
                        opacity=0.8,
                        captureevents = True
                        )

    return fig

def plot_sc_mutational_burden(sample_obj, sample_name = None):
    '''
    Plot the single-cell level mutation burden of the sample
    '''
    if sample_name is None:
        # attempt to get sample_name from h5 metadata
        sample_name = sample_obj.dna.metadata['sample_name']

    mut_prev_array = sample_obj.dna.get_attribute(
        'mut', constraint='row'
    ).sum(axis=1)
    #total_cell_num = sample_obj.dna.shape[0]

    fig = px.histogram(
        mut_prev_array,
        nbins = int(mut_prev_array.max() * 1.1 / 10)
    )

    fig.update_layout(width = 1000, height=400,\
                title_text = f'{sample_name} variant mutational prevalence (single-cells) histogram<br><sup>total: {total_cell_num} cells</sup>',\
                title_x = 0.5, title_xanchor = 'center', title_yanchor = 'top')

    fig.update_yaxes(
        title = 'number of cells',
        titlefont=dict(family='Arial', color='crimson', size=20)
    )
    fig.update_xaxes(
        title = 'number of short genomic variants carried',
        range=(0, mut_prev_array.max() * 1.1)
    )
    return fig

def plot_snv_clone(
    sample_obj,
    sample_name, 
    story_topic,
    voi,
    attribute,
    vars_to_sort_by = None,
    barcode_sort_method='hier',
    ann_map= None,
    ):
    '''
    Plot sc-heatmap for given SNV for given attribute (e.g. AF / NGT)

    inputs:
    - sample_obj: missionbio.sample object
    - sample_name: sample name
    - story_topic: story topic
    - voi: list of variants of interest
    - attribute: attribute to plot
    - vars_to_sort_by: variant to sort by
    - barcode_sort_method: 'hier'(default) or 'single_var'; method to sort barcodes
    - ann_map: annotation map

    '''
    if vars_to_sort_by is None:
        vars_to_sort_by = voi[0]

    sorted_bars = sort_for_var(
        dna = sample_obj.dna,
        vars = vars_to_sort_by,
        attribute = attribute,
        method = barcode_sort_method
        )
    #print(sorted_bars)
    fig = sample_obj.dna.heatmap(
        attribute,
        features = voi,
        bars_order = sorted_bars,
    )
    # alter colorscale for sc-AF heatmap
    if attribute == "AF_MISSING":
        fig.layout.coloraxis.colorscale = [
                [0, 'rgb(255,255,51)'], 
                [1/3, 'rgb(204,229,255)'], 
                [2/3, 'rgb(112,112,255)'],
                [1, 'rgb(255,0,0)']
            ] # @HZ: this colorscale is the best tried so far

    if ann_map is not None:
        new_x_axis_vals = []
        for var_i in fig.layout.xaxis.ticktext.tolist():
            ann = ann_map[var_i]['variant annotation']
            new_x_axis_vals.append(ann)
        fig.update_xaxes(
            ticktext= new_x_axis_vals,
            tickangle= 90, 
            tickfont=dict(
                family='Arial', 
                color='crimson', 
                size=15
                )
            )

    fig.update_layout(
        title= {
            'font': {'family': 'Arial','size': 15}, 
            'text': f'{sample_name} - {story_topic}',
            'x': 0.5,
            'y': 0.95
            }
    )

    return fig
