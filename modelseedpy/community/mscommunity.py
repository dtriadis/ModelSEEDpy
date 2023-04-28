# -*- coding: utf-8 -*-
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.community.mscompatibility import MSCompatibility
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.community.mssteadycom import MSSteadyCom
from modelseedpy.community.commhelper import build_from_species_models
from modelseedpy.core.exceptions import ObjectAlreadyDefinedError, FeasibilityError
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.fbahelper import FBAHelper
#from modelseedpy.fbapkg.gapfillingpkg import default_blacklist
from modelseedpy.core.msatpcorrection import MSATPCorrection
from cobra import Reaction
from cobra.core.dictlist import DictList
from cobra.io import save_matlab_model
from optlang.symbolics import Zero
from pandas import DataFrame
from pprint import pprint
import logging

# import itertools
import cobra
import re, os

logger = logging.getLogger(__name__)


class CommunityMembers:
    def __init__(self, community, biomass_cpd, name=None, index=None):
        self.community, self.biomass_cpd = community, biomass_cpd
        self.index = index or int(self.biomass_cpd.compartment[1:])
        self.abundance = 0
        if self.biomass_cpd in self.community.primary_biomass.metabolites:
            self.abundance = abs(self.community.primary_biomass.metabolites[self.biomass_cpd])
        if name:  self.id = name
        elif "species_name" in self.biomass_cpd.annotation:
            self.id = self.biomass_cpd.annotation["species_name"]
        else:  self.id = "Species"+str(self.index)

        logger.info("Making atp hydrolysis reaction for species: "+self.id)
        atp_rxn = self.community.util.add_atp_hydrolysis("c"+str(self.index))
        self.atp_hydrolysis = atp_rxn["reaction"]
        self.biomass_drain = None
        self.biomasses, self.reactions = [], []
        for rxn in self.community.util.model.reactions:
            rxnComp = FBAHelper.rxn_compartment(rxn)
            if not rxnComp:
                print(f"The reaction {rxn.id} strangely lacks a compartment.")
            elif int(rxnComp[1:]) == self.index and 'bio' not in rxn.name:
                self.reactions.append(rxn)
            if self.biomass_cpd in rxn.metabolites:
                if rxn.metabolites[self.biomass_cpd] == 1 and len(rxn.metabolites) > 1:
                    self.biomasses.append(rxn)
                elif len(rxn.metabolites) == 1 and rxn.metabolites[self.biomass_cpd] < 0:
                    self.biomass_drain = rxn

        if len(self.biomasses) == 0:
            logger.critical("No biomass reaction found for species " + self.id)
        if not self.biomass_drain:
            logger.info("Making biomass drain reaction for species: "+self.id)
            self.biomass_drain = Reaction(
                id="DM_"+self.biomass_cpd.id, name="DM_" + self.biomass_cpd.name, lower_bound=0, upper_bound=100)
            self.community.util.model.add_reactions([self.biomass_drain])
            self.biomass_drain.add_metabolites({self.biomass_cpd: -1})
            self.biomass_drain.annotation["sbo"] = 'SBO:0000627'

    def disable_species(self):
        for reaction in self.community.model.reactions:
            reaction_index = FBAHelper.rxn_compartment(reaction)[1:]
            if int(reaction_index) == self.index:
                reaction.upper_bound = reaction.lower_bound = 0

    def compute_max_biomass(self):
        if len(self.biomasses) == 0:
            logger.critical("No biomass reaction found for species "+self.id)
        self.community.util.add_objective(Zero, coef={self.biomasses[0].forward_variable:1})
        if self.community.lp_filename != None:
            self.community.print_lp(self.community.lp_filename+"_"+self.id+"_Biomass")
        return self.community.model.optimize()

    def compute_max_atp(self):
        if not self.atp_hydrolysis:
            logger.critical("No ATP hydrolysis found for species:" + self.id)
        self.community.util.add_objective(Zero, coef={self.atp_hydrolysis.forward_variable: 1})
        if self.community.lp_filename:
            self.community.print_lp(self.community.lp_filename + "_" + self.id + "_ATP")
        return self.community.model.optimize()


class MSCommunity:
    def __init__(self, model=None,            # the model that will be defined
                 models:list=None,            # the list of models that will be assembled into a community
                 names=None, abundances=None, # names and abundances of the community species
                 pfba = True,                 # specify whether parsimonious FBA will be simulated
                 lp_filename = None           # specify a filename to create an lp file
                 ):
        self.lp_filename, self.pfba = lp_filename, pfba
        self.gapfillings = {}

        #Define Data attributes as None
        self.solution = self.biomass_cpd = self.primary_biomass = self.biomass_drain = None
        self.msgapfill = self.element_uptake_limit = self.kinetic_coeff = self.msdb_path = None
        self.members = DictList()

        # defining the models
        model = model if not models else build_from_species_models(
            models, names=names, abundances=abundances, cobra_model=True)
        self.id = model.id
        self.util = MSModelUtil(model, True)
        self.pkgmgr = MSPackageManager.get_pkg_mgr(self.util.model)
        msid_cobraid_hash = self.util.msid_hash()
        if "cpd11416" not in msid_cobraid_hash:
            raise KeyError("Could not find biomass compound for the model.")
        other_biomass_cpds = []
        for self.biomass_cpd in msid_cobraid_hash["cpd11416"]:
            if self.biomass_cpd.compartment == "c0":
                for reaction in self.util.model.reactions:
                    if self.biomass_cpd not in reaction.metabolites:  continue
                    print(self.biomass_cpd, reaction)
                    if reaction.metabolites[self.biomass_cpd] == 1 and len(reaction.metabolites) > 1:
                        if self.primary_biomass:
                            raise ObjectAlreadyDefinedError(
                                f"The primary biomass {self.primary_biomass} is already defined,"
                                f"hence, the {reaction} cannot be defined as the model primary biomass.")
                        print('primary biomass defined', reaction)
                        self.primary_biomass = reaction
                    elif reaction.metabolites[self.biomass_cpd] < 0 and len(reaction.metabolites) == 1:
                        self.biomass_drain = reaction
            elif 'c' in self.biomass_cpd.compartment:  # else does not seem to capture built model members
                other_biomass_cpds.append(self.biomass_cpd)
        for memIndex, biomass_cpd in enumerate(other_biomass_cpds):
            name = names[memIndex]
            print(name, biomass_cpd)
            self.members.append(CommunityMembers(community=self, biomass_cpd=biomass_cpd, name=name))
        self.abundances_set = False
        if abundances:  self.set_abundance(abundances)
        self.set_objective()

    #Manipulation functions
    def set_abundance(self, abundances):
        #calculate the normalized biomass
        total_abundance = sum(list(abundances.values()))
        # map abundances to all species
        for species, abundance in abundances.items():
            if species in self.members:  self.members.get_by_id(species).abundance = abundance/total_abundance
        #remake the primary biomass reaction based on abundances
        if self.primary_biomass is None:  logger.critical("Primary biomass reaction not found in community model")
        all_metabolites = {self.primary_biomass.products[0]: 1}
        all_metabolites.update({mem.biomass_cpd: -abundances[mem.id]/total_abundance for mem in self.members})
        self.primary_biomass.add_metabolites(all_metabolites, combine=False)
        self.abundances_set = True

    def set_objective(self, target=None, targets=None, minimize=False):
        targets = targets or [self.util.model.reactions.get_by_id(target or self.primary_biomass.id).flux_expression]
        self.util.model.objective = self.util.model.problem.Objective(
            sum(targets), direction="max" if not minimize else "min")
        # TODO explore a multi-level objective

    def constrain(self, element_uptake_limit=None, kinetic_coeff=None,
                  thermo_params=None, msdb_path=None):
        if element_uptake_limit:
            self.element_uptake_limit = element_uptake_limit
            self.pkgmgr.getpkg("ElementUptakePkg").build_package(element_uptake_limit)
        if kinetic_coeff:
            self.kinetic_coeff = kinetic_coeff
            self.pkgmgr.getpkg("CommKineticPkg").build_package(kinetic_coeff, self)
        if thermo_params:
            if msdb_path:
                self.msdb_path = msdb_path
                thermo_params.update({'modelseed_db_path':msdb_path})
                self.pkgmgr.getpkg("FullThermoPkg").build_package(thermo_params)
            else:  self.pkgmgr.getpkg("SimpleThermoPkg").build_package(thermo_params)

    def interactions(self, solution=None, media=None, filename=None, export_directory=None,
                     node_metabolites=True, flux_threshold=1, visualize=True, ignore_mets=None):
        return MSSteadyCom.interactions(self, solution, media, flux_threshold, visualize, filename, export_directory,
                                        node_metabolites, show_figure=True, ignore_mets=ignore_mets)

    #Utility functions
    def print_lp(self, filename=None):
        filename = filename or self.lp_filename
        with open(filename+".lp", 'w') as out:
            out.write(str(self.util.model.solver))
            out.close()

    #Analysis functions
    def gapfill(self, media = None, target = None, minimize = False, default_gapfill_templates=None, default_gapfill_models=None,
                test_conditions=None, reaction_scores=None, blacklist=None, suffix = None, solver:str="glpk"):
        default_gapfill_templates = default_gapfill_templates or []
        default_gapfill_models = default_gapfill_models or []
        test_conditions = test_conditions or []
        reaction_scores = reaction_scores or {}
        blacklist = blacklist or []
        if not target:  target = self.primary_biomass.id
        self.set_objective(target, minimize)
        gfname = FBAHelper.medianame(media) + "-" + target
        if suffix:  gfname += f"-{suffix}"
        self.gapfillings[gfname] = MSGapfill(self.util.model, default_gapfill_templates, default_gapfill_models,
                                             test_conditions, reaction_scores, blacklist, solver)
        gfresults = self.gapfillings[gfname].run_gapfilling(media, target)
        if not gfresults:
            logger.critical("Gapfilling failed with the specified model, media, and target reaction.")
            return None
        return self.gapfillings[gfname].integrate_gapfill_solution(gfresults)

    def test_individual_species(self, media=None, allow_cross_feeding=True, run_atp=True, run_biomass=True):
        self.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        # Iterating over species and running tests
        data = {"Species": [], "Biomass": [], "ATP": []}
        for individual in self.members:
            data["Species"].append(individual.id)
            with self.util.model:
                #If no interaction allowed, iterate over all other species and disable them
                if not allow_cross_feeding:
                    for indtwo in self.members:
                        if indtwo != individual:  indtwo.disable_species()
                # If testing biomass, setting objective to individual species biomass and optimizing
                if run_biomass:   data["Biomass"].append(individual.compute_max_biomass())
                # If testing atp, setting objective to individual species atp and optimizing
                if run_atp:  data["ATP"].append(individual.compute_max_atp())
        df = DataFrame(data)
        logger.info(df)
        return df

    def atp_correction(self,core_template, atp_medias, atp_objective="bio2", max_gapfilling=None, gapfilling_delta=0):
        self.atpcorrect = MSATPCorrection(self.util.model,core_template, atp_medias,
                                          atp_objective="bio2", max_gapfilling=None, gapfilling_delta=0)

    # TODO evaluate the comparison of this method with MICOM
    def predict_abundances(self,media=None,pfba=True,kinetic_coeff = None):
        with self.util.model:   # WITH, here, discards changes after each simulation
            kinetic_coeff = kinetic_coeff or self.kinetic_coeff or 2000
            self.pkgmgr.getpkg("CommKineticPkg").build_package(kinetic_coeff,self)
            self.util.model.objective = self.util.model.problem.Objective(Zero,direction="max")
            self.util.model.objective .set_linear_coefficients(
                {species.biomasses[0].forward_variable: 1 for species in self.members})
            self.run_fba(media, pfba)
            return self._compute_relative_abundance_from_solution()

    def run_fba(self, media=None, pfba=False, fva_reactions=None):
        if media is not None:  self.util.add_medium(media)
        return self._set_solution(self.util.run_fba(None, pfba, fva_reactions))

    def _compute_relative_abundance_from_solution(self,solution = None):
        if not solution and not self.solution:
            logger.warning("The simulation lacks any flux.")
            return None
        data = {"Species": [], "Abundance": []}
        totalgrowth = sum([self.solution.fluxes[member.biomasses[0].id] for member in self.members])
        if totalgrowth == 0:
            logger.warning("The community did not grow!")
            return None
        for member in self.members:
            data["Member"].append(member.id)
            data["Abundance"].append(self.solution.fluxes[member.biomasses[0].id] / totalgrowth)
        df = DataFrame(data)
        logger.info(df)
        return df

    def _set_solution(self, solution):
        if solution.status != "optimal":
            FeasibilityError(f'The solution is sub-optimal, with a(n) {solution} status.')
            self.solution = None
            self.print_lp()
            save_matlab_model(self.util.model, self.util.model.name + ".mat")
        self.solution = solution
        logger.info(self.util.model.summary())
        return self.solution
