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

def justifyDB(msdb_path:str, changes_path:str="MSDB_corrections.json"):
    msdb = from_local(msdb_path)

    db_element_counts = combine_elements(*[met.elements for rxn in msdb.reactions for met in rxn.metabolites])
    # print(len(msdb.compounds), db_element_counts)
    stoich_vars, charge_vars, stoich_diff_pos, stoich_diff_neg, charge_diff_pos, charge_diff_neg = {}, {}, {}, {}, {}, {}
    stoich_constraints, charge_constraints = {}, {}
    variables, constraints, names = [], [], []
    objective = tupObjective("database_correction", [], "min")
    time1 = process_time()
    print("Defining variables and objective", end="\t")
    mets_frequency = Counter([met.id for rxn in msdb.reactions for met in rxn.metabolites])
    for met in msdb.compounds:
        # mass balance
        # stoich_constraints[met.id], stoich_vars[met.id] = {}, {}
        stoich_diff_pos[met.id], stoich_diff_neg[met.id] = {}, {}
        met_elements = met.elements if met.elements != {} else list(db_element_counts.keys())
        for ele in met_elements:
            # stoich_vars[met.id][ele] = tupVariable(_check_names(f"{met.id}~{ele}", names),
            #                                        tuple(sorted((met.elements[ele] * 0.5, met.elements[ele] * 1.5))))
            stoich_diff_pos[met.id][ele] = tupVariable(f"{met.id}~{ele}_diffpos")
            stoich_diff_neg[met.id][ele] = tupVariable(f"{met.id}~{ele}_diffneg")
            ## StoichVar_{m,e} + (diffpos_{m,e} - diffneg_{m,e}) = stoich_{m,e}
            # stoich_constraints[met.id][ele] = tupConstraint(f"{met.id}_{ele}_diff", expr={
            #     "elements": [stoich_diff_pos[met.id][ele].name, stoich_vars[met.id][ele].name, -met.elements[ele],
            #                  {"elements": [stoich_diff_neg[met.id][ele].name, -1], "operation": "Mul"}],
            #     "operation": "Add"})
            objective.expr.extend([{"elements": [
                    {"elements": [stoich_diff_pos[met.id][ele].name, mets_frequency[met.id]], "operation": "Mul"},
                    {"elements": [stoich_diff_neg[met.id][ele].name, mets_frequency[met.id]], "operation": "Mul"}],
                "operation": "Add"}])
        # charge balance
        # TODO remove the bounded limitations
        # charge_vars[met.id] = tupVariable(f"{met.id}~charge", tuple(sorted((met.charge * 0.5, met.charge * 1.5))))
        charge_diff_pos[met.id] = tupVariable(f"{met.id}~charge_diffpos")
        charge_diff_neg[met.id] = tupVariable(f"{met.id}~charge_diffneg")
        ## ChargeVar_{m} + (diffpos_{m} - diffneg_{m}) = charge_{m}
        # charge_constraints[met.id] = tupConstraint(f"{met.id}_charge_diff", expr={
        #     "elements": [charge_diff_pos[met.id].name, charge_vars[met.id].name, -met.charge,
        #                  {"elements": [charge_diff_neg[met.id].name, -1], "operation": "Mul"}],
        #     "operation": "Add"})
        # update the objective expression and store the constructed variables and constraints
        objective.expr.extend([{"elements": [
                {"elements": [charge_diff_pos[met.id].name, mets_frequency[met.id]], "operation": "Mul"},
                {"elements": [charge_diff_neg[met.id].name, mets_frequency[met.id]], "operation": "Mul"}],
            "operation": "Add"}])
        variables.extend([#*stoich_vars[met.id].values(),
                          *stoich_diff_pos[met.id].values(), *stoich_diff_neg[met.id].values(),
                          charge_diff_pos[met.id], charge_diff_neg[met.id],
                          # charge_vars[met.id]
                          ])
        # constraints.extend([*stoich_constraints[met.id].values(), *charge_constraints.values()])
    time2 = process_time()
    print(f"Done after {(time2-time1)/60} minutes")
    print("Defining constraints", end="\t")
    # print(len(constraints))
    empty_reactions = []
    for rxn in msdb.reactions:
        rxn_element_counts = combine_elements(*[met.elements for met in rxn.metabolites])
        if rxn_element_counts == {}:  empty_reactions.append(rxn.id)
        charge_constraints[rxn.id] = {}
        for met in rxn.metabolites:
            # sum_m^M( (-diffpos_{m} + diffneg_{m}) * n_{met,rxn} ) = charge_{m} * n_{met,rxn} , per reaction rxn
            charge_constraints[rxn.id][met.id] = tupConstraint(f"{rxn.id}_{met.id}_charge", expr={
                "elements": [{"elements": [charge_diff_pos[re.sub(r"(_\d+)", "", met.id)].name,
                                           -rxn.metabolites[met]],
                              "operation": "Mul"},
                             {"elements": [charge_diff_neg[re.sub(r"(_\d+)", "", met.id)].name,
                                           rxn.metabolites[met]],
                              "operation": "Mul"},
                             met.charge * rxn.metabolites[met]],
                "operation": "Add"})
        # sum_m^M( (-diffpos_{m,e} + diffneg_{m,e}) * n_{met,rxn} ) = -stoich_{m,e} * n_{met,rxn} , per reaction rxn
        stoich_constraints[rxn.id] = {}
        for ele, count in rxn_element_counts.items():
            stoich_constraints[rxn.id][ele] = tupConstraint(
                f"{rxn.id}_{ele}", expr={"elements": [count], "operation":"Add"})
            for met in rxn.metabolites:
                if ele not in met.elements:  continue
                stoich_constraints[rxn.id][ele].expr["elements"].extend([
                    {"elements": [stoich_diff_pos[re.sub(r"(_\d+)", "", met.id)][ele].name, -rxn.metabolites[met]],
                     "operation": "Mul"},
                    {"elements": [stoich_diff_neg[re.sub(r"(_\d+)", "", met.id)][ele].name, rxn.metabolites[met]],
                     "operation": "Mul"}])
        constraints.extend([*stoich_constraints[rxn.id].values(), *charge_constraints[rxn.id].values()])
    if empty_reactions:
        print(f"The {empty_reactions} reactions lack any metabolites with "
              f"defined formula, and thus are not constrained.", end="\t")
    time3 = process_time()
    print(f"Done after {(time3-time2)/60} minutes")

    # construct the model
    print("Constructing the model", end="\t")
    print(list(map(len, [variables, constraints, stoich_constraints, charge_constraints, objective.expr])), end="\t")
    optlang_model = OptlangHelper.define_model("Correct_MSDB", variables, constraints, objective, True)
    with open("Correct_MSDB.lp", 'w') as lp:  lp.write(optlang_model.to_lp())
    with open("Correct_MSDB.json", 'w') as jsonOut:  dump(optlang_model.to_json(), jsonOut, indent=3)
    print(f"Done after {(process_time()-time3)/60} minutes")

    # acquire the optimized minimum from the model
    print("Starting optimization.", end="\t")
    before = process_time()
    solution = optlang_model.optimize()
    after = process_time()
    print(f"Done after \t{(after-before)/60} minutes")
    if solution != "optimal":  FeasibilityError(f"The optimization is {solution}.")
    return optlang_model
    print("Exporting primals", end="\t")
    # export the changes
    if changes_path is not None:
        # evaluate proposed changes from the optimization
        proposed_changes = {}
        for varName, val in optlang_model.primal_values.items():
            # print(var, val)
            content, diffName = varName.split("_")
            metID, context = content.split("~")
            # print(rxnID, metID, ele)
            met = msdb.compounds.get_by_id(metID)
            missing_formula = False
            if context == "charge":
                original_val = met.charge
            elif context in met.elements:
                original_val = met.elements[context]
            elif context not in met.elements:
                missing_formula = True ; original_val = None
                print(f"The {metID} metabolite is mis-categorized with the {context} element,"
                      f" or possibly the formula and thus elements attribute is not defined.")
            if val != 0 or missing_formula:
                val = -val if "neg" in diffName else val
                if metID not in proposed_changes:  proposed_changes[metID] = {}
                if context in proposed_changes[metID]:  proposed_changes[metID][context]["proposed"] += val
                else:  proposed_changes[metID][context] = {
                    "original": original_val, "proposed": val+original_val if original_val is not None else val}
        # export the proposed changes
        with open(changes_path, "w") as out:  dump(proposed_changes, out, indent=3)
        print(f"Done after {(process_time()-after)/60} minutes")  ;  return optlang_model, proposed_changes
    print(f"Done after {(process_time()-after)/60} minutes")  ;  return optlang_model

# if __name__ == "__main__":
#     # accept the command-line arguments
#     from argparse import ArgumentParser
#     from pathlib import Path
#
#     parser = ArgumentParser()
#     parser.add_argument("--msdb_path", type=Path, default=".")
#     parser.add_argument("--changes_path", type=Path, default="MSDB_corrections.json")
#     args = parser.parse_args()
#     justifyDB(args.msdb_path, args.changes_path)
