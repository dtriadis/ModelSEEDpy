# from django import shortcuts
from numpy import nan, isnan, unique
import jinja2
import os, re, math, json

package_dir = os.path.abspath(f"{os.path.dirname(__file__)}/..")

# Helper function
def round_float_str(s, precision=6):
    """
    Given a string containing floats delimited by whitespace
    round each float in the string and return the rounded string.

    Parameters
    ----------
    s : str
        String containing floats
    precision : int
        Number of decimals to round each float to
    """
    round_str = ''
    for token in s.split():
        try:
            f = round(float(token), precision)
            round_str += str(f)
        except:
            round_str += token
        round_str += ' '
    return round_str

# Helper functions for formating
missing_format = lambda x: x if x else '-'
nan_format = lambda x: '-' if isnan(x) else x
round_format = lambda x: nan_format(round(x, 6))
equation_format = lambda rxn: round_float_str(rxn.build_reaction_string(use_metabolite_names=True))

# Helper function to determine reaction class
def class_formater(rxnID, fva_sol):
    min_zero = math.isclose(fva_sol.minimum[rxnID], 0, abs_tol=1e-07)
    max_zero = math.isclose(fva_sol.maximum[rxnID], 0, abs_tol=1e-07)

    min_ = 0 if min_zero else fva_sol.minimum[rxnID]
    max_ = 0 if max_zero else fva_sol.maximum[rxnID]

    if min_zero and max_zero:
        return 'blocked'
    if min_ > 0 or max_ < 0:
        return 'essential'
    #if min_ < 0 or max_ > 0:
    return 'functional'

# Helper function to format reaction/ex-reaction display data
def reaction_formater(model, fba_sol, fva_sol, ex):
    """ex specifies exchange reaction"""
    if fva_sol is None:
        return json.dumps([])

    # Select either exchange or normal reactions
    if ex:
        rxnIDs = fva_sol.loc[fva_sol.index.str[:2] == 'EX'].index
    else:
        rxnIDs = fva_sol.loc[fva_sol.index.str[:2] != 'EX'].index

    # Get reaction objects from ids
    rxns = map(model.reactions.get_by_id, rxnIDs)

    return json.dumps([{'id': rxnID,
                        'flux': round_format(fba_sol.fluxes[rxnID]),
                        'min_flux': round_format(fva_sol.minimum[rxnID]),
                        'max_flux': round_format(fva_sol.maximum[rxnID]),
                        'class': class_formater(rxnID, fva_sol),
                        'equation': equation_format(rxn),
                        'name': missing_format(rxn.name)}
                       for rxnID, rxn in zip(rxnIDs, rxns)])

# Helper function for formatting model summary
def model_summary(model):
    def rxn_name(rxnID):
        if rxnID is nan:
            return '-'

        try:
            name = model.reactions.get_by_id('EX_' + rxnID).name
            return name + f'\n({rxnID})'
        except:
            try:
                name = model.reactions.get_by_id(rxnID).name
                return name + f'\n({rxnID})'
            except:
                try:
                    name = model.metabolites.get_by_id(rxnID).name
                    return name + f'\n({rxnID})'
                except:
                    return '-'

    df =  model.summary().to_frame().applymap(rxn_name)
    return [row[1:] for row in df.itertuples()]

# Helper function for formatting ATP summary
def atp_summary_formatter(model):
    """Returns list of ATP summary values if metabolites are found
       or a message stating they could not be found. Also return a
       bool specicifying whether or not summary exists."""

    # Select ATP metabolite
    if 'atp_c' in model.metabolites:
        df = model.metabolites.atp_c.summary().to_frame()
    elif 'ATP_c' in model.metabolites:
        df = model.metabolites.ATP_c.summary().to_frame()
    elif 'cpd00002_c0' in model.metabolites:
        df = model.metabolites.cpd00002_c0.summary().to_frame()
    else:
        # Empty ATP summary
        msg = 'None of the < atp_c, ATP_c, and cpd00002_c0 > metabolites are defined. ' \
              'Add one of these metabolites in the model to display an ATP summary.'
        return msg, False

    rxns = [(model.reactions.get_by_id(rxnID), rxnID)
            for rxnID in df.index.get_level_values(1)]

    df['FLUX'] = df['FLUX'].apply(round_format)
    df['PERCENT'] = df['PERCENT'].apply(lambda x: f'{round(x, 2)}%')
    df['SIDE']= df.index.get_level_values(0)
    df['NAME_ID'] = [rxn.name + f'\n({rxnID})' for rxn, rxnID in rxns]
    df['REACTION_STRING'] = [equation_format(rxn) for rxn, _ in rxns]

    atp_summary = [list(row) for row in zip(
        df['SIDE'], df['NAME_ID'], df['PERCENT'], df['FLUX'], df['REACTION_STRING'])]
    return atp_summary, True

# Helper function for formating essential genes
def essential_genes_formatter(model, essential_genes):
    if not essential_genes:
        return json.dumps([])
    return json.dumps([{'name': missing_format(gene.id),
                        'essential': 'Yes' if gene in essential_genes else 'No'}
                      for gene in model.genes])


######################## REPORT FUNCTIONS ########################

def fba_report(model, fba_sol, fva_sol, essential_genes, model_id,
               media_id, export_html_path="steadycom_report.html"):
    """Build output report and return string of html."""

    atp_summary, is_atp_summary = atp_summary_formatter(model)
    # these key-value pairs are the variables-values that are substituted into the HTML report
    content = {'summary':     model_summary(model),
               'atp_summary': {'is_atp_summary': is_atp_summary, 'summary': atp_summary},
               'overview':    [{'name': 'Model',                    'value': model_id},
                               {'name': 'Media',                    'value': media_id},
                               {'name': 'Optimization status',      'value': fba_sol.status},
                               {'name': 'Objective',                'value': model.objective},
                               {'name': 'Number of reactions',      'value': len(model.reactions)},
                               {'name': 'Number of compounds',      'value': len(model.metabolites)}],
               'reaction_tab': {
                   'is_reactions': fva_sol is not None,
                   'reactions': reaction_formater(model, fba_sol, fva_sol, ex=False),
                   'help': 'Select FVA setting and rerun to produce results.'
                },
               'ex_reaction_tab': {
                   'is_reactions': fva_sol is not None,
                   'reactions': reaction_formater(model, fba_sol, fva_sol, ex=True),
                   'help': 'Select FVA setting and rerun to produce results.'
                },
               'essential_genes_tab': {
                    'is_essential_genes': len(essential_genes) > 0,
                    'essential_genes': essential_genes_formatter(model, essential_genes),
                    'help': 'Select simulate all single KO to produce results.'
                }
           }
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(package_dir),
                             autoescape=jinja2.select_autoescape(['html', 'xml']))
    html_report = env.get_template(os.path.join("core", "fba_template.html")).render(content)
    with open(export_html_path, "w") as out:
        out.writelines(html_report)
    return html_report

def steadycom_report(flux_df, exMets_df, export_html_path="steadycom_report.html"):
    total_table = flux_df.iloc[-len(flux_df.columns):]
    total_table.index = [i.replace("zz_", "") for i in total_table.index]
    flux_df = flux_df.iloc[:-len(flux_df.columns)]
    content = {'flux_table': flux_df.style.background_gradient().to_html(),
               "total_table": total_table.style.background_gradient().to_html(),
               "exchanged_mets": exMets_df.to_html()}
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(package_dir),
                             autoescape=jinja2.select_autoescape(['html', 'xml']))
    html_report = env.get_template(os.path.join("community", "steadycom_template.html")).render(content)
    with open(export_html_path, "w") as out:
        out.writelines(html_report)
    return html_report

def smetana_report(df, mets, export_html_path="smetana_report.html"):
    # refine the DataFrame into a heatmap
    def quantify_MRO(element):
        print(element)
        return float(re.sub("(\s\(.+\))", "", str(element)))

    # construct the Heatmap
    heatmap_df = df.copy(deep=True) # takes some time
    heatmap_df_index = zip(heatmap_df["model1"].to_numpy(), heatmap_df["model2"].to_numpy())
    heatmap_df.index = [" ++ ".join(index) for index in heatmap_df_index]
    heatmap_df.index.name = "model1 ++ model2"
    heatmap_df = heatmap_df.loc[~heatmap_df.index.duplicated(), :]
    heatmap_df = heatmap_df.drop(["model1", "model2"], axis=1)
    heatmap_df["mro_model1"] = heatmap_df["mro_model1"].apply(quantify_MRO)
    heatmap_df["mro_model2"] = heatmap_df["mro_model2"].apply(quantify_MRO)
    heatmap_df = heatmap_df.astype(float)
    heatmap_df["mip"] = heatmap_df["mip"].astype(int)
    # populate the HTML template with the assembled simulation data from the DataFrame -> HTML conversion
    content = {'table': df.to_html(),
               "heatmap": heatmap_df.style.background_gradient().to_html(),
               "mets": mets}
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(package_dir),
                             autoescape=jinja2.select_autoescape(['html', 'xml']))
    html_report = env.get_template(os.path.join("community", "smetana_template.html")).render(content)
    with open(export_html_path, "w") as out:
        out.writelines(html_report)
    return html_report

def msdb_justification_report(proposed_changes, export_html_path="msdb_correction_report.html"):
    # construct the report tables
    from pandas import DataFrame
    proposed_changes_df = DataFrame(proposed_changes).T
    proposed_changes_df["normalized_change"] = abs(
        proposed_changes_df["original"] - proposed_changes_df["proposed"]
    ) / min([proposed_changes_df["original"], proposed_changes_df["proposed"]])
    heatmap_changes = proposed_changes_df.copy(deep=True)
    # construct and export the HTML files
    content = {"table": proposed_changes_df, "heatmap": heatmap_changes.style.background_gradient().to_html()}
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(package_dir),
                             autoescape=jinja2.select_autoescape(['html', 'xml']))
    html_report = env.get_template(os.path.join("biochem", "msdb_template.html")).render(content)
    with open(export_html_path, "w") as out:
        out.writelines(html_report)
    return html_report
