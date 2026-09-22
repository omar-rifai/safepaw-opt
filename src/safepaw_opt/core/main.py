import pulp
from src.safepaw_opt.core.optimization import declare_constraints, set_obj_fn
from src.safepaw_opt.core.utils.data_utils import package_results, get_summary_results
import typer
import os
from pathlib import Path
import time 
from typing import Any
from pydantic import BaseModel

 
class OptVars(BaseModel):
    P : Any
    P_gk: Any
    P_gkr: Any
    Q: Any
    Delta_plus: Any
    Delta_moins:Any
    z_hl_plus: Any
    z_hl_moins: Any
    s_hl: Any


def run_driver(params_system):
    LP = pulp.LpProblem('regional_case_mix', pulp.LpMaximize)

    P = {g: {k: {r: {a: {h: pulp.LpVariable(f"P_{g}_{k}_{r}_{a}_{h}", lowBound=0, upBound=0 if h not in params_system["O_gk"][g][k] else None)
                         for h in params_system["H"]}
                         for a in params_system["A_idx"][g][k]}
                         for r in params_system["R"]}
                         for k in params_system["K_idx"][g]}
                         for g in params_system["G"]}

    P_gkr = {g: {k: {r: pulp.LpVariable(f"P_gkr_{g}_{k}_{r}", lowBound=0)
                     for r in params_system["R"]} 
                     for k in params_system["K_idx"][g]}
                     for g in params_system["G"]}

    P_gk = {g: {k: pulp.LpVariable(f"P_gk_{g}_{k}", lowBound=0)
                for k in params_system["K_idx"][g]}
                for g in params_system["G"]}

    Q = {g: {k: {r: {a: {h: pulp.LpVariable(f"Q_{g}_{k}_{r}_{a}_{h}", lowBound=0)
                         for h in params_system["H"]}
                         for a in params_system["A_idx"][g][k]}
                         for r in params_system["R"]}
                         for k in params_system["K_idx"][g]}
                         for g in params_system["G"]}

    Delta_plus = {h: {l: pulp.LpVariable(f"Delta_plus_{h}_{l}", lowBound=0)
                      for l in params_system["L"]}
                      for h in params_system["H"]}

    Delta_moins = {h: {l: pulp.LpVariable(f"Delta_moins_{h}_{l}", lowBound=0)
                       for l in params_system["L"]}
                       for h in params_system["H"]}

    z_hl_plus = {h: {l: pulp.LpVariable(f"z_hl_plus_{h}_{l}", lowBound=0)
                     for l in params_system["L"]}
                     for h in params_system["H"]}

    z_hl_moins = {h: {l: pulp.LpVariable(f"z_hl_moins_{h}_{l}", lowBound=0) 
                      for l in params_system["L"]}
                      for h in params_system["H"]}

    s_hl = {h: {l: pulp.LpVariable(f"s_hl{h}_{l}") 
                      for l in params_system["L"]}
                      for h in params_system["H"]}

    
    vars_system = OptVars(P=P,P_gkr=P_gkr,P_gk=P_gk,Q=Q,Delta_plus=Delta_plus, Delta_moins=Delta_moins, z_hl_plus=z_hl_plus, z_hl_moins=z_hl_moins, s_hl=s_hl)
    start = time.time()
    set_obj_fn(LP, P_gk, P, Delta_plus, Delta_moins, s_hl, params_system)
    end = time.time()
    print(f"setting obj fn time {end-start}")
    print("Declaring Constraints...")
    start = time.time()
    declare_constraints(LP, vars_system, params_system)
    end = time.time()
    print(f"declaring constraints {end-start}")
    print("Running solver...")
    if "GRB_LICENSE_FILE" in os.environ:
        LP.solve(pulp.GUROBI(msg=True, timeLimit=600))
    else:
        print("Using HiGHs solver")
        LP.solve(pulp.HiGHS(msg=1,gapRel=1e-4, timeLimit=600))
    status =  pulp.LpStatus[LP.status]

    if LP.sol_status!=1:
        status = "Infeasible"
    objective = pulp.value(LP.objective)

    return status, objective, vars_system

    


def main(params_file : Path = typer.Argument(..., help="Path to JSON parameters file")):
    import json

    if not params_file.exists():
        typer.echo(f"Error: file {params_file} does not exist")
        raise typer.Exit(code=1)
    
    out_dir = Path("./results/")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"results_{params_file.name}"

    with open(params_file) as fp:
        params_system = json.load(fp)

    status, objective, vars_system = run_driver(params_system)
    dict_results = package_results(vars_system, params_system, out_file)

    return dict_results


if __name__ == "__main__":
    typer.run(main)



    
