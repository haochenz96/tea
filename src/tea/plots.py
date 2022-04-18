# Venn diagrams
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted, venn3
from matplotlib.lines import Line2D
from pathlib import Path

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

def plot_var_sc_mutational_prev(sample_obj, known_bulk_somatic_vars, min_mut_prev_of_interest, sample_name = None):
    '''
    Plot the single-cell level mutation prevalence of the sample
    
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
    count = 0
    y_val = 0.5
    ay_val = -30 
    for var_i in known_bulk_somatic_vars:
        if var_i in mut_prev_array.index:
            y_val += count * 0.15
            count += 1
            
            fig.add_annotation(x = mut_prev_array[var_i], y = y_val,
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


    fig.update_layout(width = 1000, height=400,\
                    title_text = f'{sample_name} variant mutational prevalence )single-cells) histogram<br><sup>total: {total_cell_num} cells</sup>',\
                    title_x = 0.5, title_xanchor = 'center', title_yanchor = 'top')
    return fig

# def plot_sc_mutational_burden(sample_obj, sample_name = None):
#     '''
#     Plot the single-cell level mutation burden of the sample
#     '''
#     if sample_name is None:
#         # attempt to get sample_name from h5 metadata
#         sample_name = sample_obj.dna.metadata['sample_name']
