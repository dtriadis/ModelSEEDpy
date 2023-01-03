from icecream import ic

from modelseedpy import FBAHelper
from modelseedpy.core.exceptions import ObjectAlreadyDefinedError, ParameterError, NoFluxError
from modelseedpy.community.commhelper import CommHelper
from itertools import combinations
from optlang.symbolics import Zero
from pandas import DataFrame
from matplotlib import pyplot
from numpy import array
import networkx
import sigfig
import os, re


def add_collection_item(met_name, normalized_flux, flux_threshold, ignore_mets,
                        species_collection, first, second):
    if flux_threshold and normalized_flux <= flux_threshold:
        return species_collection
    if not any([re.search(x, met_name, flags=re.IGNORECASE) for x in ignore_mets]):
        species_collection[first][second].append(re.sub(r"(_\w\d$)", "", met_name))
    return species_collection


class MSSteadyCom:
    
    @staticmethod
    def compute(
            mscommodel,                          # The MSCommunity object of the model (mandatory to prevent circular imports)
            solution = None,                     # the COBRA simulation solution that will be parsed and visualized
            media=None,                          # The media in which the community model will be simulated
            # names=None, abundances=None,         # names and abundances of the community species
            flux_threshold: int = 1,             # The threshold of normalized flux below which a reaction is not plotted
            visualize: bool = True,              # specifies whether the net flux will be depicted in a network diagram
            filename: str = 'cross_feeding_diagram.svg',  # Cross-feeding figure export name
            export_directory: str = None,        # specifies the directory to which the network diagram and associated datatable will be exported, where None does not export the content
            node_metabolites: bool = True,       # specifies whether the metabolites of each node will be printed
            x_offset: float = 0.15,              # specifies the x-axis buffer between each species node and its metabolite list in the network diagram
            show_figure: bool = True,            # specifies whether the figure will be printed to the console
            ignore_mets=None                     # cross-fed exchanges that will not be displayed in the graphs
            ):
        #Check for solution
        solution = solution or mscommodel.run(media)
        if not solution:
            raise ParameterError("A solution must be provided. Interactions are computed from a solution.")
        if all(array(list(solution.fluxes.values)) == 0):
            print(list(solution.fluxes.values))
            raise NoFluxError("The simulation lacks any flux.")

        #Initialize data
        metabolite_data, species_data, species_collection = {}, {"Environment":{}}, {"Environment":{}}
        data = {"IDs":[],"Metabolites/Donor":[], "Environment":[]}
        species_list = {}

        # track extracellularly exchanged metabolites
        exchange_mets_list = mscommodel.util.exchange_mets_list()
        for met in exchange_mets_list:
            data["IDs"].append(met.id)
            data["Metabolites/Donor"].append(re.sub(r"(_\w\d$)", "", met.name))
            metabolite_data[met.id] = {"Environment": 0}
            metabolite_data[met.id].update({individual.id: 0 for individual in mscommodel.species})

        # computing net metabolite flux from each reaction
        for individual in mscommodel.species:
            species_data[individual.id], species_collection[individual.id] = {}, {}
            species_list[individual.index] = individual
            data[individual.id] = []
            for other in mscommodel.species:
                species_data[individual.id][other.id] = 0
                species_collection[individual.id][other.id] = []
            species_data["Environment"][individual.id] = species_data[individual.id]["Environment"] = 0
            species_collection["Environment"][individual.id] = []
            species_collection[individual.id]["Environment"] = []

        for rxn in mscommodel.util.model.reactions:
            if rxn.id[0:3] == "EX_":
                cpd = list(rxn.metabolites.keys())[0]
                # the Environment takes the opposite perspective to the members
                metabolite_data[cpd.id]["Environment"] += -solution.fluxes[rxn.id]
            rxn_index = int(FBAHelper.rxn_compartment(rxn)[1:])
            if not any([met not in exchange_mets_list for met in rxn.metabolites]) or rxn_index not in species_list:
                continue
            for met in rxn.metabolites:
                if met.id not in metabolite_data:
                    continue
                metabolite_data[met.id][species_list[rxn_index].id] += solution.fluxes[rxn.id]*rxn.metabolites[met]

        # translating net metabolite flux into species interaction flux
        ignore_mets = ignore_mets if ignore_mets is not None else ["h2o_e0", "co2_e0"]
        for met in exchange_mets_list:
            #Iterating through the metabolite producers
            # !!! Why are fluxes normalized?
            total = sum([max([metabolite_data[met.id][individual.id], 0]) for individual in mscommodel.species
                         ]) + max([metabolite_data[met.id]["Environment"], 0])
            for individual in mscommodel.species:
                ## calculate metabolic consumption of a species from the environment
                if metabolite_data[met.id][individual.id] < Zero:
                    if metabolite_data[met.id]["Environment"] <= Zero:
                        continue
                    normalized_flux = abs(metabolite_data[met.id][individual.id]
                                          * metabolite_data[met.id]["Environment"]) / total
                    species_data["Environment"][individual.id] += normalized_flux
                    species_collection = add_collection_item(met.name, normalized_flux, flux_threshold, ignore_mets,
                                                             species_collection, "Environment", individual.id)
                ## calculate and track metabolic donations between a member and another or the environment
                elif metabolite_data[met.id][individual.id] > Zero:
                    for other in mscommodel.species:
                        ### filter against organisms that do not consume
                        if metabolite_data[met.id][other.id] >= Zero:
                            continue
                        normalized_flux = abs(metabolite_data[met.id][individual.id]
                                              * metabolite_data[met.id][other.id])/total
                        species_data[individual.id][other.id] += normalized_flux
                        # print(met.name, normalized_flux > flux_threshold, not any(
                        #         [re.search(x, met.name, flags=re.IGNORECASE) for x in ignore_mets]))
                        species_collection = add_collection_item(met.name, normalized_flux, flux_threshold, ignore_mets,
                                                                 species_collection, individual.id, other.id)
                    ## calculate donations to the environment
                    if metabolite_data[met.id]["Environment"] >= Zero:
                        continue
                    normalized_flux = abs(metabolite_data[met.id][individual.id]
                                          * metabolite_data[met.id]["Environment"])/total
                    ic(normalized_flux)
                    species_data[individual.id]["Environment"] += normalized_flux
                    species_collection = add_collection_item(met.name, normalized_flux, flux_threshold, ignore_mets,
                                                             species_collection, individual.id, "Environment")

        # construct the dataframes
        for metID in metabolite_data:
            for individual in mscommodel.species:
                data[individual.id].append(metabolite_data[metID][individual.id])
            data["Environment"].append(metabolite_data[metID]["Environment"])

        ## process the fluxes dataframe
        data["IDs"].append("zz_Environment")
        data["Metabolites/Donor"].append(0)
        for individual in mscommodel.species:
            data[individual.id].append(species_data["Environment"][individual.id])
        data["Environment"].append(0)
        for individual in mscommodel.species:
            for other in mscommodel.species:
                data[individual.id].append(species_data[individual.id][other.id])
            data["Environment"].append(species_data[individual.id]["Environment"])
            data["IDs"].append(f"zz_Species{individual.index}")
            data["Metabolites/Donor"].append(individual.id)

        # if len(set(list(map(len, list(data.values()))))) != 1:
        #     print([(col, len(content)) for col, content in data.items()])
        cross_feeding_df = DataFrame(data)
        cross_feeding_df.index = list(map(str, cross_feeding_df["IDs"])) # sorted()
        cross_feeding_df.index.name = "Metabolite/Donor ID"
        cross_feeding_df.drop(['IDs', "Metabolites/Donor"], axis=1, inplace=True)
        cross_feeding_df = cross_feeding_df.loc[(cross_feeding_df != 0).any(axis=1)]
        cross_feeding_df.sort_index(inplace=True)

        ## process the identities dataframe
        exchanged_mets = {"Environment": [" "], "Donor ID": ["Environment"]}
        exchanged_mets.update({ind.id: [] for ind in mscommodel.species})
        for individual in mscommodel.species:
            ### environment exchanges
            exchanged_mets[individual.id].append("; ".join(species_collection["Environment"][individual.id]))
            exchanged_mets["Environment"].append("; ".join(species_collection[individual.id]["Environment"]))
            ### member exchanges
            exchanged_mets["Donor ID"].append(individual.id)
            for other in mscommodel.species:
                exchanged_mets[individual.id].append("; ".join(species_collection[individual.id][other.id]))

        # if len(set(list(map(len, list(exchanged_mets.values()))))) != 1:
        #     print([(col, len(content)) for col, content in exchanged_mets.items()])
        exMets_df = DataFrame(exchanged_mets)
        exMets_df.index = list(map(str, exMets_df["Donor ID"]))
        exMets_df.index.name = "Donor ID"
        exMets_df.drop(["Donor ID"], axis=1, inplace=True)
        exMets_df.sort_index(inplace=True)
        exMets_df.fillna(" ")
        # logger.info(cross_feeding_df)

        # graph the network diagram
        if visualize:
            MSSteadyCom._visualize_cross_feeding(
                cross_feeding_df, exMets_df, filename, export_directory, mscommodel.species,
                node_metabolites, x_offset, show_figure)

        return cross_feeding_df, exMets_df

    @staticmethod
    def _visualize_cross_feeding(cross_feeding_df, exMets_df, filename, export_directory,
            species, node_metabolites = True, x_offset = 0.15, show_figure = True):
        # define species and the metabolite fluxes
        graph = networkx.Graph()
        species_nums = {}
        for species in species:
            species_nums[species.index] = set()
            graph.add_node(species.index)
            for col in exMets_df.columns:
                if 'Species' in col and col != f"Species{species.index}":
                    species_nums[species.index].update(exMets_df.at[f"Species{species.index}", col].split('; '))

        # define the net fluxes for each combination of two species
        ID_prefix = "zz_Species" if any(["zz_Species" in x for x in list(cross_feeding_df.index)]) else "Species"
        for index1, index2 in combinations(list(species_nums.keys()), 2):
            species_1 = min([index1, index2]) ; species_2 = max([index1, index2])
            if not all([x in cross_feeding_df.index for x in [f'{ID_prefix}{species_1}', f'{ID_prefix}{species_2}']]):
                continue
            species_1_to_2 = cross_feeding_df.at[f'{ID_prefix}{species_1}', f'Species{species_2}']
            species_2_to_1 = cross_feeding_df.at[f'{ID_prefix}{species_2}', f'Species{species_1}']
            interaction_net_flux = sigfig.round(species_1_to_2 - species_2_to_1, 3)
            graph.add_edge(species_1, species_2, flux=interaction_net_flux)

        # compose the network diagram of net fluxes
        pos = networkx.circular_layout(graph)
        if node_metabolites:
            for species in pos:
                x, y = pos[species]
                metabolites = '\n'.join(species_nums[species])
                pyplot.text(x+x_offset, y, metabolites)
        networkx.draw_networkx(graph, pos)
        labels = networkx.get_edge_attributes(graph, 'flux')
        networkx.draw_networkx_edge_labels(graph, pos, edge_labels=labels)

        # export and view the figure
        filename = filename or 'cross_feeding_diagram.svg'
        export_directory = export_directory or os.getcwd()
        pyplot.savefig(os.path.join(export_directory, filename), bbox_inches="tight", transparent=True)
        csv_filename = re.sub(r"(\.\w+)", ".csv", filename)
        csv_filename = csv_filename.replace("_diagram", "")
        cross_feeding_df.to_csv(os.path.join(export_directory, csv_filename))
        if show_figure:
            pyplot.show()
