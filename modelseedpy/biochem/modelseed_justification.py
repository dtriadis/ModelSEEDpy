from modelseedpy.core.optlanghelper import OptlangHelper, Bounds, tupVariable, tupConstraint, tupObjective
from modelseedpy.core.exceptions import FeasibilityError, ObjectAlreadyDefinedError
from modelseedpy.biochem import from_local
from collections import Counter
from time import process_time
from json import dump
import re

def combine_elements(*ele_dics):
    ele_dic = {}
    for dic in ele_dics:
        # print(dic, ele_dics)
        for ele, count in dic.items():
            if ele in ele_dic:  ele_dic[ele] += count
            else:  ele_dic[ele] = count
    return ele_dic

def _check_names(name, names):
    if name in names:  raise ObjectAlreadyDefinedError(f"The {name} is already defined in the model.")
    names.append(name)
    return name

def justifyDB(msdb_path:str, changes_path:str=None):
    # db_element_counts = combine_elements([met.elements for rxn in msdb.reactions for met in rxn.metabolites])
    msdb = from_local(msdb_path)

    mass_variables, charge_variables, mass_constraints, charge_constraints = {}, {}, {}, {}
    variables, constraints, names = [], [], []
    objective = tupObjective("database_correction", [], "min")
    time1 = process_time()
    print("Defining variables and constraints")
    mets_frequency = Counter([met.id for rxn in msdb.reactions for met in rxn.metabolites])
    for rxn in msdb.reactions:
        rxn_element_counts = combine_elements(*[met.elements for met in rxn.metabolites])
        mass_variables[rxn.id], mass_constraints[rxn.id] = {}, {}
        charge_variables[rxn.id], charge_constraints[rxn.id] = {}, {}
        # print(f"{rxn.id}\t{rxn.reaction}\t\t\t\t\t\t\t\t\t\t", end="\r")
        for met in rxn.metabolites:
            # mass balance
            mass_variables[rxn.id][met.id] = {ele: tupVariable(
                _check_names(f"{rxn.id}_{met.id}~{ele}", names),
                tuple(sorted((met.elements[ele]*0.5, met.elements[ele]*1.5)))
            ) for ele in met.elements}
            # sum_{m,r,e}^{M,R,E} ( StoichVar_{m,r,e} * freq_{m} - stoich_{m,e} * freq_{m} )
            objective.expr.extend([{
                "elements": [
                    {"elements": [mass_variables[rxn.id][met.id][ele].name, mets_frequency[met.id]], "operation": "Mul"},
                    -met.elements[ele]*mets_frequency[met.id]],
                "operation": "Add"} for ele in mass_variables[rxn.id][met.id]])
            # charge balance
            charge_variables[rxn.id][met.id] = tupVariable(
                f"{rxn.id}_{met.id}~charge", tuple(sorted((met.charge*0.5, met.charge*1.5))))
            # sum_{m,r}^{M,R} ( ChargeVar_{m,r} * freq_{m} - charge_{m,r} * freq_{m} )
            objective.expr.extend([{
                "elements": [
                    {"elements": [charge_variables[rxn.id][met.id].name, mets_frequency[met.id]], "operation": "Mul"},
                    -met.charge*mets_frequency[met.id]],
                "operation": "Add"}])
            # add the variables
            variables.extend([*mass_variables[rxn.id][met.id].values(), charge_variables[rxn.id][met.id]])
        # add the constraints
        # TODO possibly constrain all compartment versions of a compound to be equivalent
        for ele in rxn_element_counts:
            # sum_m^M( Ele_(met,rxn) * n_(met,rxn) ) = 0 , for a given reaction rxn
            mass_constraints[rxn.id][ele] = tupConstraint(f"{rxn.id}_{ele}", expr={
                "elements": [{"elements": [mass_variables[rxn.id][met.id][ele].name, met.elements[ele]],
                              "operation": "Mul"} for met in rxn.metabolites if ele in met.elements],
                "operation": "Add"})
        # sum_m^M( Charge_(met,rxn) * n_(met,rxn) ) = 0 , for a given reaction rxn
        charge_constraints[rxn.id] = tupConstraint(f"{rxn.id}_charge", expr={
            "elements": [{"elements": [charge_variables[rxn.id][met.id].name, met.charge],
                          "operation": "Mul"} for met in rxn.metabolites],
            "operation": "Add"})
        constraints.extend(list(mass_constraints[rxn.id].values())+[charge_constraints[rxn.id]])

    # construct the model
    optlang_model = OptlangHelper.define_model("Correct_MSDB", variables, constraints, objective, True)
    with open("Correct_MSDB.lp", 'w') as lp:  lp.write(optlang_model.to_lp())
    with open("Correct_MSDB.json", 'w') as jsonOut:  dump(optlang_model.to_json(), jsonOut, indent=3)
    print(f"The Model is constructed: {(process_time()-time1)/60} minutes")

    # acquire the optimized minimum from the model
    print("Starting optimization.")
    before = process_time()
    solution = optlang_model.optimize()
    after = process_time()
    print(f"The optimization concluded:\t{(after-before)/60} minutes")
    if solution != "optimal":  FeasibilityError(f"The optimization is {solution}.")
    # export the changes
    if changes_path is not None:
        # evaluate proposed changes from the optimization
        proposed_changes = {}
        for varName, val in optlang_model.primal_values.items():
            # print(var, val)
            if "charge" not in varName:
                rxnID, metID_ele = varName.split("_", 1)
                metID, ele = metID_ele.split("~")
                # print(rxnID, metID, ele)
                met = msdb.compounds.get_by_id("_".join(metID.split("_")[:-1]))
                ele_stoich = met.elements[ele]
                if val != ele_stoich:
                    if metID not in proposed_changes:
                        proposed_changes[metID] = {}
                    proposed_changes[metID][ele] = {"original": ele_stoich, "proposed": val}
                continue
            # print(varName)
            rxnID, metID = varName.split("_", 1)
            metID = metID.split("~")[0]
            # print(rxnID, metID, ele)
            met = msdb.compounds.get_by_id("_".join(metID.split("_")[:-1]))
            if val != met.charge:
                if metID not in proposed_changes:
                    proposed_changes[metID] = {}
                proposed_changes[metID]["charge"] = {"original": met.charge, "proposed": val}
        # export the proposed changes
        with open(changes_path, "w") as out:
            dump(proposed_changes, out, indent=3)
        return optlang_model, proposed_changes
    return optlang_model

if __name__ == "__main__":
    # accept the command-line arguments
    from argparse import ArgumentParser
    from pathlib import Path

    parser = ArgumentParser()
    parser.add_argument("--msdb_path", type=Path, default=".")
    parser.add_argument("--changes_path", type=Path, default="MSDB_corrections.json")
    args = parser.parse_args()
    justifyDB(args.msdb_path, args.changes_path)
