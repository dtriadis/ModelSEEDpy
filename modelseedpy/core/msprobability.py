from cobrakbase.core.kbasefba.fbamodel_from_cobra import CobraModelConverter
from cobra.io import write_sbml_model, read_sbml_model
from cobra import Model, Reaction
from numpy import load as npload
from os import path, environ
from json import load, dump
from glob import glob
from os import path
import cobrakbase
import re


with open(f"reaction_counts/all_models.json", "r") as jsonIn:  reaction_counts = load(jsonIn)
with open("model_gcf_mapping.json", "r") as jsonIn:   model_gcf_mapping = load(jsonIn)
with open("unique_asv_mappings.json", "r") as jsonIn:  unique_asv_mappings = load(jsonIn)
# with open("broken_model_gcfs.npy", "r") as aryIn:  broken_model_gcfs = npload(aryIn)


class MSProbability:
        
    @staticmethod
    def megaModel(kbase_api, clades_paths, threshold = .5, cobrakbase_path=""):
        # compute the reaction frequency of the models in a given clade
        if not kbase_api:
            with open(cobrakbase_path) as token_file:   kbase_api = cobrakbase.KBaseAPI(token_file.readline())
        broken_models = []
        # models_paths = glob(f"{models_path}/*.xml")
        for clade, paths in clades_paths.items():
            for model_path in paths:
                reaction_counts = {}
                try:
                    model = read_sbml_model(model_path)
                    for rxn in model.reactions:
                        if rxn.id in reaction_counts:  reaction_counts[rxn.id] += 1
                        else: reaction_counts[rxn.id] = 1
                        # TODO storing a list of the rxn objects will save computational effort in the subsequent step
                    print(f"{model.id}\t\t\t\t\t\t\t\t\t\t\t", end="\r")
                except Exception as e:
                    print("broken", model_path)
                    broken_models.append(model_path)
            reaction_counts = {"numMembers": len(paths)-len(broken_models)}
            reaction_counts.update({rxnID:(count/reaction_counts["numMembers"]) for rxnID,count in reaction_counts.items()})
            with open(f"reaction_counts/{clade}.json", "w") as jsonOut:   dump(reaction_counts, jsonOut, indent=3)

            # constructing the probabilistic clade model
            megaModel = CobraModelConverter(Model(clade, f"MegaModel for {clade} from {reaction_counts['numMembers']} members")).build()
            remaining_rxnIDs = set(list(reaction_counts.keys()))
            captured_reactions, captured_rxnIDs = [], set()
            print("\n", clade)
            for model_path in paths:
                try:
                    model = read_sbml_model(model_path)
                    print(model.id)
                except Exception as e:
                    print(f"Broken: {model_path}") ; continue
                captured_reactions.extend([rxn for rxn in model.reactions if rxn.id not in captured_rxnIDs])
                captured_rxnIDs.update([rxn.id for rxn in model.reactions])
                remaining_rxnIDs -= captured_rxnIDs
            if captured_reactions == []:
                print(f"No models for {clade} are defined.")
                continue
            megaModel.add_reactions(list(captured_reactions))
            for rxn in megaModel.reactions:
                rxn.probability = reaction_counts[rxn.id]
            print(len(megaModel.reactions), len(reaction_counts), megaModel.reactions[0].probability)
            write_sbml_model(megaModel, path.join(path.dirname(model_path), clade))