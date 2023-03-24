from modelseedpy.core.optlanghelper import OptlangHelper, Bounds, tupVariable, tupConstraint, tupObjective, isIterable, define_term
from modelseedpy.core.exceptions import FeasibilityError, ParameterError, ObjectAlreadyDefinedError, NoFluxError
from modelseedpy.biochem import from_local
from pprint import pprint


def combine_elements(*ele_dics):
    ele_dic = {}
    for dic in ele_dics:
        for ele in dic:
            if ele in ele_dic:  ele_dic[ele] += dic[ele]
            else:  ele_dic[ele] = dic[ele]
    return ele_dic


def justifyDB(msdb_path:str):
    msdb = from_local(msdb_path)
    db_element_counts = combine_elements([met.elements for rxn in msdb.reactions for met in rxn.metabolites])

    mass_variables, charge_variables, mass_constraints, charge_constraints = {}, {}, {}, {}
    variables, constraints = [], []
    objective = tupObjective("database_correction", [], "min")
    for rxn in msdb.reactions:
        rxn_element_counts = combine_elements([met.elements for met in rxn.metabolites])
        mass_variables[rxn.id], mass_constraints[rxn.id] = {}, {}
        charge_variables[rxn.id], charge_constraints[rxn.id] = {}, {}
        for met in rxn.metabolites:
            # mass balance
            mass_variables[rxn.id][met.id] = {ele: tupVariable(f"{rxn.id}_{met.id}_{ele}") for ele in met.elements}
            objective.expr.extend([{"elements": [mass_variables[rxn.id][met.id][ele], -met.elements[ele]],
                                    "operation": "Add"}
                                   for ele in met.elements])
            mass_constraints[rxn.id][met.id] = tupConstraint(f"{rxn.id}_{met.id}_elements", expr={
                "elements": [{"elements": [mass_variables[rxn.id][met.id][ele], met.elements[ele]], "operation": "Mul"}
                             for ele in met.elements],
                "operation": "Add"})
            # charge balance
            charge_variables[rxn.id][met.id] = tupVariable(f"{rxn.id}_{met.id}_charge")
            objective.expr.extend([{"elements": [charge_variables[rxn.id][met.id], -met.charge], "operation": "Add"}])
            # store the constructed objects
            variables.extend([mass_variables[rxn.id][met.id], charge_variables[rxn.id][met.id]])
            constraints.extend([mass_constraints[rxn.id][met.id]])
        # charge balance constraint
        charge_constraints[rxn.id] = tupConstraint(f"{rxn.id}_charge", expr={
            "elements": [{"elements": [charge_variables[rxn.id][met.id], met.charge], "operation": "Mul"}
                         for met in rxn.metabolites],
            "operation": "Add"})
        constraints.extend([charge_constraints[rxn.id]])

    # construct the model
    optlang_model = OptlangHelper.define_model("Correct_MSDB", variables, constraints, objective, True)
    with open("Correct_MSDB.lp", 'w') as lp:
        lp.write(optlang_model.to_lp())

    # acquire the optimized minimum from the model
    return optlang_model.optimize()
