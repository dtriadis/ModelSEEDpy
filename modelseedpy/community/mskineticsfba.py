# -*- coding: utf-8 -*-

from scipy.constants import milli, hour, minute, day, femto
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy import MSModelUtil
from optlang import Constraint
from modelseedpy.core.fbahelper import FBAHelper
from collections import OrderedDict
from optlang.symbolics import Zero
from numpy import log10, nan, mean
from warnings import warn
from matplotlib import pyplot
from pprint import pprint
from datetime import date
from math import inf
import pandas
import json, re, os


class MetaboliteError(Exception):
    pass

def _x_axis_determination(total_time):
    time = total_time * minute
    if time <= 600:  return minute, "s"
    if time > 600:  return 1, "min"
    if time > 7200:  return 1/hour, "hr"
    return 1/day, "days"

class dFBAPkg(BaseFBAPkg):
    # TODO code conventionally defined dFBA (only exchange kinetics) or make this is a special case of the same package as the Kinetics model
    pass



class MSKineticsFBA(BaseFBAPkg):
    def __init__(self, model, msdb_path:str, warnings: bool = True, verbose: bool = False,
                 printing: bool = False, jupyter: bool = False):
        # define the parameter and variable dictionaries
        BaseFBAPkg.__init__(self, model, "BasedFBA", {"met": "metabolite"}, {"conc": "metabolite"})
        # self.pkgmgr.addpkgs(["FullThermoPkg"])
        # self.parameters["modelseed_api"] = FBAHelper.get_modelseed_db_api(self.parameters["msdb_path"])
        self.warnings, self.verbose, self.printing, self.jupyter = warnings, verbose, printing, jupyter
        self.model_util = MSModelUtil(model)
        if msdb_path:
            self.pkgmgr.addpkgs(["FullThermoPkg"])
            self.pkgmgr.getpkg("FullThermoPkg").build_package({"modelseed_db_path": msdb_path})
        self.met_ids = OrderedDict({met.id: met.name for met in self.model_util.model.metabolites})
        self.met_names = list(self.met_ids.values())

    # user functions
    def print_lp(self, filename=None):
        BaseFBAPkg.print_lp(self, filename.removesuffix(".lp"))

    def _check_datum(self, datum):
        if "substituted_rate_law" not in datum:
            print(f"RateLawError: The {datum} datum lacks a rate law.")
            return False
        remainder = re.sub("([0-9A-Za-z/()e\-\+\.\*])", "", datum["substituted_rate_law"])
        if remainder != "":
            print(f'RateLawError: The {datum["substituted_rate_law"]}'
                  f' rate law contains unknown characters: {remainder}')
            return False

    def simulate(self, kinetics_path: str = None, kinetics_data: dict = None, initial_concentrations_M: dict = None,  # a dictionary of the initial metabolic concentrations, which supplants concentrations from the defined kinetics data
                 total_time_min: float = 200, timestep_min: float = 20, export_name = None, export_directory = None,
                 chemostat_L: float = None, feed_profile: dict = None, chemo_exchange_rate_L_hr: float = None,
                 temperature: float = 25, p_h: float = 7, cellular_dry_mass_fg: float = 222, cellular_fL: float = 1,
                 conc_figure_title = "Metabolic perturbation", included_metabolites: list = None,labeled_plots = True,
                 visualize = True, export = True):
        # define the dataframe for the time series content
        feed_profile = feed_profile or {}
        included_metabolites, self.solutions = included_metabolites or [], []
        self.cellular_dry_mass_fg, self.cellular_fL = cellular_dry_mass_fg*femto, cellular_fL*femto
        self.parameters["timesteps"] = int(total_time_min/timestep_min)
        self.timestep_min, self.total_time_min = timestep_min, total_time_min
        self.constrained = OrderedDict()
        self.minimum = inf

        # define experimental conditions
        self.parameters["pH"], self.parameters["temperature"], self.variables["elapsed_time"] = p_h, temperature, 0
        # define initial concentrations
        self._initial_concentrations(initial_concentrations_M, kinetics_path, kinetics_data)
        chemostat_requirements = [chemostat_L is not None, feed_profile != {}, chemo_exchange_rate_L_hr is not None]
        # determine the reactions for which kinetics are predefined
        self.defined_reactions = {rxn.name:rxn for rxn in self.model_util.model.reactions if rxn.name in kinetics_data}
        # execute FBA for each timestep, then calculate custom fluxes, constrain the model, and update concentrations
        for self.timestep in range(1, self.parameters["timesteps"] + 1):
            colNum = self.timestep_min * self.timestep_min
            self.col, self.previous_col = f"{colNum} min", f"{colNum-1} min"
            self.concentrations[self.col] = [float(0)] * len(self.concentrations.index)
            self.fluxes[self.col] = [nan] * len(self.fluxes.index)
            # create a metabolite variable that prevents negative concentrations
            timestep_hr = self.timestep_min * (minute / hour)
            for met in self.model_util.model.metabolites:
                if met.id not in self.initial_concentrations_M:  continue
                coef = {}
                for rxn in met.reactions:
                    ## The product of the reaction stoichiometry and the timestep
                    stoich = abs(timestep_hr * rxn.metabolites[met])
                    coef[rxn.forward_variable], coef[rxn.reverse_variable] = stoich, -stoich
                ## build the metabolite constraint
                concentration = -self.concentrations.at[met.name, self.previous_col]
                if met.id in self.constraints["conc"]:
                    self.constraints["conc"][met.id].lb = concentration
                    continue
                self.constraints["conc"][met.id] = Constraint(Zero, lb=concentration, ub=None, name=f"{met.id}_conc")
                self.model_util.create_constraint(self.constraints["conc"][met.id], coef)
                # var = BaseFBAPkg.build_variable(self,"met",0,None,"continuous",met)
                # BaseFBAPkg.build_constraint(self,"conc",0,None,{met:1},met)
            # calculate the flux
            for reaction_name in self.kinetics_data:
                fluxes = []
                for source in self.kinetics_data[reaction_name]:
                    datum = self.kinetics_data[reaction_name][source]
                    if not self._check_datum(datum):  continue
                    # define each variable concentration
                    conc_dict = {var: self.concentrations.at[self.met_ids[datum["met_id"][var]], self.previous_col]*milli
                                 for var in datum["met_id"] if len(var) == 1}
                    # define rate law variables; calculate flux; average or overwrite the flux based on data criteria
                    locals().update(conc_dict)
                    flux = eval(datum["substituted_rate_law"])
                    if ("metadata" not in self.kinetics_data[reaction_name][source]
                        or self.__find_data_match(reaction_name, source) == 'a'):  fluxes.append(flux)
                    else:  fluxes = [flux]

                flux = mean(fluxes)
                if reaction_name in self.defined_reactions:
                    self.__set_constraints(reaction_name, flux)
                    self.fluxes.at[reaction_name, self.col] = flux
                    if self.printing:   print("\n")
                elif self.warnings:  warn(f"ReactionError: The {reaction_name} reaction, with a"
                                          f" flux of {flux}, is not described by the model.")
            # execute the COBRA model
            solution = self.model_util.model.optimize()
            self.solutions.append(solution)
            ## add previously undefined concentrations
            self.fluxes[self.col] = [solution.fluxes[rxn.id] for rxn in self.fluxes.index
                                     if not FBAHelper.isnumber(self.fluxes.at[rxn.name, self.col])]
            for met in self.model_util.model.metabolites:
                self.concentrations.at[self.met_ids[met.id], self.col] = 0
                for rxn in met.reactions:  # flux units: mmol/(g_(dry weight)*hour)
                    stoich = rxn.metabolites[met]
                    flux = self.fluxes.at[rxn.name, self.col]
                    delta_conc = stoich * flux * self.timestep_min * (minute / hour) * (
                                self.cellular_dry_mass_fg / self.cellular_fL)
                    self.concentrations.at[self.met_ids[met.id], self.col] += delta_conc
            if all(chemostat_requirements):
                self.chemical_moles[self.col] = (self.concentrations[self.col] * milli * chemostat_L)
                self._chemostat(feed_profile, chemo_exchange_rate_L_hr, chemostat_L)
            elif any(chemostat_requirements): warn("The < chemostat_L > , < feed_profile >, and < chemo_exchange_rate_L_hr >"
                                                   " parameters must all be defined to simulate a chemostat.")
            self.variables["elapsed_time"] += self.timestep_min
            if self.printing:
                print(f"\nobjective value for timestep {self.timestep_min}: ", self.solutions[-1].objective_value)

        # identify the chemicals that dynamically changed in concentrations
        self.changed, self.unchanged = set(), set()
        for met_name in self.met_names:
            first = self.concentrations.at[met_name, "0 min"]
            final = self.concentrations.at[met_name, self.col]
            if first != final:  self.changed.add(met_name)
            else:  self.unchanged.add(met_name)

        # visualize concentration changes over time
        if visualize:   self._visualize(conc_figure_title, included_metabolites, labeled_plots)
        if export:      self._export(export_name, export_directory)
        if self.verbose:  print("\n\nChanged metabolite  concentrations\n",
                                "=" * 2 * len("changed metabolites"), f"\n{self.changed}",
                                "\nConstrained reactions:", self.constrained.keys())
        elif self.printing:
            if self.jupyter:  pandas.set_option("max_rows", None)  ; display(self.concentrations, self.fluxes)
            if self.unchanged == set():  print("\nAll of the metabolites changed concentration over the simulation")
            else:  print("\n\nUnchanged metabolite concentrations\n",
                         "=" * 2 * len("unchanged metabolites"), f"\n{self.unchanged}")
        return self.concentrations, self.fluxes

    # utility functions
    def _initial_concentrations(self, initial_concentrations_M, kinetics_path: str = None, kinetics_data: dict = None):
        # define kinetics of the system
        self.kinetics_data = {}
        if kinetics_path:
            with open(kinetics_path) as data:  self.kinetics_data = json.load(data)
        elif kinetics_data:         self.kinetics_data = kinetics_data.copy()
        if not self.kinetics_data:  raise ValueError("Kinetics data must be defined.")

        # define the DataFrames
        self.col = "0 min"
        self.concentrations = pandas.DataFrame(data=[float(0)]*len(self.met_names),
                                               index=set(self.met_names), columns=[self.col])
        self.concentrations.index.name = "metabolite (\u0394mM)"
        self.chemical_moles = self.concentrations.copy()
        self.fluxes = pandas.DataFrame(index=set(rxn.name for rxn in self.model_util.model.reactions), columns=[self.col])
        self.fluxes.index.name = "reactions (mmol/(hr*g_(dw)))"

        # parse the kinetics data
        initial_concentrations = {}
        for reaction_name, content in self.kinetics_data.items():
            for condition, datum in content.items():
                for var, conc in datum["initial_concentrations_M"].items():
                    met_id = datum["met_id"][var]
                    name = self.met_ids[met_id]
                    if name in self.met_names:
                        self.concentrations.at[name, self.col] += conc/milli
                        initial_concentrations[met_id] = self.concentrations.at[name, self.col]
                    elif self.warnings:  warn(f"KineticsError: The {name} reagent ({var}) in the"
                                              f" {datum['substituted_rate_law']} rate law is not defined by the model.")

        # incorporate custom initial concentrations
        self.initial_concentrations_M = initial_concentrations_M
        if self.initial_concentrations_M:
            for met_id in self.initial_concentrations_M:
                met_name = self.met_ids[met_id]
                if met_name not in self.concentrations.index:
                    if self.warnings:  warn(f"InitialConcError: The {met_id} ({met_name})"
                                            f" metabolite is not defined by the model.")
                else:
                    self.concentrations.at[met_name, self.col] = (self.initial_concentrations_M[met_id] * milli)
                    initial_concentrations[met_id] = self.concentrations.at[name, self.col]
        self.initial_concentrations_M = initial_concentrations

    def _visualize(self, conc_fig_title, included_metabolites, labeled_plots):
        # define the figure
        pyplot.rcParams['figure.figsize'] = (11, 7)
        pyplot.rcParams['figure.dpi'] = 150
        self.figure, ax = pyplot.subplots()
        ax.set_title(conc_fig_title)
        ax.set_ylabel("Concentrations (mM)")

        x_axis_scalar, unit = _x_axis_determination(self.total_time_min)
        ax.set_xlabel("Time " + unit)
        legend_list = []
        times = [t * self.timestep_min * x_axis_scalar for t in range(self.parameters["timesteps"] + 1)]

        # determine the plotted metabolites and the scale of the figure axis
        bbox = (1, 1)
        if not included_metabolites:
            bbox = (1.7, 1)
            # 1e-2 is an arbitrary concentration threshold for plotting on the figure
            included_metabolites = [chem for chem in self.changed
                                    if max(self.concentrations.loc[[chem]].values[0].tolist()) > 1e-2]

        log_axis = False
        minimum, maximum = inf, -inf
        printed_concentrations = {}
        for chem in self.changed:
            if chem not in included_metabolites:    continue
            concentrations = self.concentrations.loc[[chem]].values[0].tolist()
            maximum = max(maximum, max([x if x > 1e-9 else 0 for x in concentrations]))
            minimum = min(minimum, min([x if x > 1e-9 else 0 for x in concentrations]))
            # plot chemicals with perturbed concentrations
            ax.plot(times, concentrations)
            if len(chem) > 25:      chem = list(self.met_ids.keys())[self.met_names.index(chem)]
            if not concentrations[0] < 1e-9:    legend_list.append(chem)
            else:       legend_list.append(f"(rel) {chem}")

            # design the proper location of the overlaid labels in the figure
            if not labeled_plots:   continue
            for i, conc in enumerate(concentrations):
                if conc <= 1e-9:    continue
                x_value = i * self.timestep_min
                vertical_adjustment = 0
                if x_value in printed_concentrations:
                    vertical_adjustment = (maximum - minimum) * 0.05
                    if log_axis:    vertical_adjustment = log10(maximum - minimum) / 3
                ax.text(x_value, conc + vertical_adjustment, f"{chem} - {round(conc, 4)}", ha="left")
                printed_concentrations[x_value] = conc
                break

        # finalize figure details
        if maximum > 10 * minimum:  ax.set_yscale("log")
        ax.set_xticks(times)
        ax.grid(True)
        ax.legend(legend_list, title="Changed chemicals", loc="upper right",
                  bbox_to_anchor=bbox, title_fontsize="x-large", fontsize="large")

    def _export(self, export_name: str, export_directory: str):
        # define a unique simulation name
        directory = os.getcwd() if export_directory is None else os.path.dirname(export_directory)
        if not export_name:     export_name = "-".join([re.sub(" ", "_", str(x)) for x in [
            date.today(), "dFBA", self.model_util.model.name, f"{self.total_time_min} min"] ])
        simulation_number = -1
        while os.path.exists(os.path.join(directory, export_name)):
            simulation_number += 1
            export_name = re.sub("(\-\d+$)", "", export_name)
            export_name = "-".join([export_name, str(simulation_number)])
        self.parameters["simulation_path"] = self.simulation_path = os.path.join(directory, export_name)
        os.mkdir(self.simulation_path)

        # export simulation content
        self.fluxes.to_csv(os.path.join(self.simulation_path, "fluxes.csv"))
        self.concentrations.to_csv(os.path.join(self.simulation_path, "concentrations.csv"))
        times = self.fluxes.columns
        with open(os.path.join(self.simulation_path, "objective_values.csv"), "w") as obj_val:
            obj_val.write("min,objective_value")
            for sol in self.solutions:
                time = re.sub("(\smin)", "", times[self.solutions.index(sol)])
                obj_val.write(f"\n{time},{sol.objective_value}")

        # export the parameters
        parameters = {"parameter": list(self.parameters.key()), "value": list(self.parameters.values())}
        parameters_table = pandas.DataFrame(parameters)
        parameters_table.to_csv(os.path.join(self.simulation_path, "parameters.csv"))

        # export the figure
        self.figure.savefig(os.path.join(self.simulation_path, "changed_concentrations.svg"))
        if self.verbose and not self.jupyter:   self.figure.show()

    def _chemostat(self, feed_profile:dict, chemo_exchange_rate_L_hr, chemostat_L):
        L_changed = chemo_exchange_rate_L_hr * self.timestep_min
        # chemostat addition
        for chem_id, conc in feed_profile.items():
            chem_name = self.met_ids[chem_id]
            self.chemical_moles.at[chem_name, self.col] += conc * L_changed
            self.concentrations.at[chem_name, self.col] = (
                self.chemical_moles.at[chem_name, self.col] / milli / chemostat_L)  # normalize to the chemostat volume

        # chemostat subtraction
        for met in self.model_util.model.metabolites:
            if met.compartment[0] != "e":   continue
            ## update the chemical moles
            self.chemical_moles.at[met.name, self.col] -= (self.concentrations.at[met.name, self.col] * L_changed)
            ## define the chemical concentration
            self.concentrations.at[met.name, self.col] = (
                    self.chemical_moles.at[met.name, self.col] / milli / chemostat_L)

    # nested functions
    def __find_data_match(self, reaction_name: str,
                          source: str,  # specifies which enzymatic datum among many will be used
                          ):
        # identifies the datum whose experimental conditions most closely matches the simulation conditions
        temperature_deviation = ph_deviation = 0
        if FBAHelper.isnumber(self.kinetics_data[reaction_name][source]["metadata"]["Temperature"]):
            temp = float(self.kinetics_data[reaction_name][source]["metadata"]["Temperature"])
            temperature_deviation = (abs(self.parameters["temperature"] - temp) / self.parameters["temperature"])
        if FBAHelper.isnumber(self.kinetics_data[reaction_name][source]["metadata"]["pH"]):
            pH = float(self.kinetics_data[reaction_name][source]["metadata"]["pH"])
            ph_deviation = (abs(self.parameters["pH"] - pH) / self.parameters["pH"])

        # equally weight between temperature and pH deviation from the simulation conditions
        old_minimum = self.minimum
        deviation = mean(temperature_deviation, ph_deviation)
        self.minimum = min(deviation, self.minimum)
        return "a" if old_minimum == self.minimum else "w" # append or write a list of data

    def __set_constraints(self, reaction_name: str, flux: float):
        rxn = self.defined_reactions[reaction_name]
        rxn_name = re.sub(" ", "_", rxn.name)
        if rxn_name in self.constrained:
            self.model_util.model.remove_cons_vars(self.constrained[rxn_name])
            self.model_util.model.solver.update()
        self.constrained[rxn_name] = self.model_util.model.problem.Constraint(
            rxn.flux_expression, lb=flux, ub=flux, name=f"{rxn_name}_kinetics")
        self.model_util.create_constraint(self.constrained[rxn_name])
        if self.verbose:    print(self.model_util.model.constraints[f"{rxn_name}_kinetics"])
