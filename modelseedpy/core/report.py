import os
import json
import math
import jinja2
import numpy as np
from django import shortcuts

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
nan_format = lambda x: '-' if np.isnan(x) else x
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
        if rxnID is np.nan:
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

# Call this function to build the report
# TODO The arguments that are specific for executing FBA models can be removed, to leave only rendering-specific arguments and logic
def build_report(model, fba_sol, fva_sol,
                 essential_genes, model_id, media_id):
    """Build output report and return string of html."""

    atp_summary, is_atp_summary = atp_summary_formatter(model)
    # these key-value pairs are the variables-values that are substituted into the HTML report
    context = {'summary':     model_summary(model),
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

    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(os.path.dirname(os.path.realpath(__file__))),
        autoescape=jinja2.select_autoescape(['html', 'xml']))
    # Return string of html
    return env.get_template('template.html').render(context)

def steadycom_report(flux_df, exMets_df):
    total_table = flux_df.iloc[-len(flux_df.columns):]
    total_table.index = [i.replace("zz_", "") for i in total_table.index]
    content = {'flux_table': flux_df.to_html(),
               "total_table": total_table.to_html(),
               "exchanged_mets": exMets_df.to_html()}
    package_dir = "/".join(os.path.dirname(os.path.realpath(__file__)).split("\\")[:-1])
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(package_dir),
        autoescape=jinja2.select_autoescape(['html', 'xml']))
    # print(os.path.dirname(__file__))
    # Return string of html
    # print(os.path.exists(os.path.join(os.path.dirname(__file__), '..', "community", "steadycom_output.html")))
    return env.get_template("/".join(["community", "steadycom_output.html"])).render(content)
