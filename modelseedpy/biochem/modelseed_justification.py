from modelseedpy.core.optlanghelper import OptlangHelper, Bounds, tupVariable, tupConstraint, tupObjective
from modelseedpy.core.exceptions import FeasibilityError
from modelseedpy.biochem import from_local
from time import process_time
import re

def combine_elements(*ele_dics):
    ele_dic = {}
    for dic in ele_dics:
        # print(dic, ele_dics)
        for ele, count in dic.items():
            if ele in ele_dic:  ele_dic[ele] += count
            else:  ele_dic[ele] = count
    return ele_dic

def parse_primals(primal_values):
    # TODO implement some refining logic for the primal values
    values = {var: val for var, val in primal_values.items()}
    return values

def _check_names(name, names):
    names.append(name)
    return name

def justifyDB(msdb_path:str, primals_path:str=None):
    # db_element_counts = combine_elements([met.elements for rxn in msdb.reactions for met in rxn.metabolites])
    msdb = from_local(msdb_path)

    mass_variables, charge_variables, mass_constraints, charge_constraints = {}, {}, {}, {}
    variables, constraints, names = [], [], []
    objective = tupObjective("database_correction", [], "min")
    for rxn in msdb.reactions:
        rxn_element_counts = combine_elements(*[met.elements for met in rxn.metabolites])
        mass_variables[rxn.id], mass_constraints[rxn.id] = {}, {}
        charge_variables[rxn.id], charge_constraints[rxn.id] = {}, {}
        print(rxn.reaction, end="\r")
        for met in rxn.metabolites:
            metID = re.sub("(\_\w?\d+)", "", met.id)
            # mass balance
            mass_variables[rxn.id][met.id] = {ele: tupVariable(_check_names(f"{rxn.id}_{metID}_{ele}", names))
                                              for ele in met.elements if f"{rxn.id}_{metID}_{ele}" not in names}
            objective.expr.extend([{"elements": [mass_variables[rxn.id][met.id][ele].name, -met.elements[ele]],
                                    "operation": "Add"} for ele in met.elements])
            # else:  print(rxn.reaction, end="\r")
            # charge balance
            charge_variables[rxn.id][met.id] = tupVariable(f"{rxn.id}_{metID}_charge")
            objective.expr.extend([{"elements": [charge_variables[rxn.id][met.id].name, -met.charge],
                                    "operation": "Add"}])
            # add the variables
            variables.extend([*mass_variables[rxn.id][met.id].values(), charge_variables[rxn.id][met.id]])
        # add the constraints
        for ele in rxn_element_counts:
            mass_constraints[rxn.id][ele] = tupConstraint(f"{rxn.id}_{ele}", expr={
                "elements": [{"elements": [mass_variables[rxn.id][met.id][ele].name, met.elements[ele]],
                              "operation": "Mul"} for met in rxn.metabolites if ele in met.elements],
                "operation": "Add"})
        charge_constraints[rxn.id] = tupConstraint(f"{rxn.id}_charge", expr={
            "elements": [{"elements": [charge_variables[rxn.id][met.id].name, met.charge],
                          "operation": "Mul"} for met in rxn.metabolites],
            "operation": "Add"})
        constraints.extend(list(mass_constraints[rxn.id].values())+[charge_constraints[rxn.id]])

    # construct the model
    optlang_model = OptlangHelper.define_model("Correct_MSDB", variables, constraints, objective, True)
    with open("Correct_MSDB.lp", 'w') as lp:  lp.write(optlang_model.to_lp())

    # acquire the optimized minimum from the model
    print("Starting optimization.")
    before = process_time()
    solution = optlang_model.optimize()
    after = process_time()
    print(f"The optimization concluded:\t{(after-before)/60} minutes")
    if solution != "optimal":
        FeasibilityError(f"The optimization is {solution}.")
    if primals_path is not None:
        from json import dump
        with open(primals_path, "w") as out:
            dump(optlang_model.primal_values, out, indent=3)
    return optlang_model

if __name__ == "__main__":
    # accept the command-line arguments
    from argparse import ArgumentParser
    from pathlib import Path

    parser = ArgumentParser()
    parser.add_argument("--msdb_path", type=Path, required=True)
    parser.add_argument("--primals_path", type=Path, required=False)
    args = parser.parse_args()

    # execute the optimization
    justifyDB(args.msdb_path, args.primals_path)
