
from pulp import *

# Define the objective function for LP
    
def set_obj_fn(LP, P_gk, P, Delta_plus, Delta_minus, s_hl, params_system):
    """
    Define a weighted-sum scalarized multi-objective problem
    """
    w_rh = params_system["w_rh"]

    indices_a = ((g, k, r, a, h)
    for g in params_system["G"]
    for k in params_system["K_idx"][g]
    for r in params_system["R"]
    for a in params_system["A_idx"][g][k]
    for h in params_system["H"])

    coeff_rh = {(r,h): params_system["D"] * w_rh[r][h] for r in params_system["R"] for h in params_system["H"]}

    match params_system["mode"]:
        case "slack":
            LP += lpSum(-s_hl[h][l] for h in params_system["H"] for l in params_system["L"])
        case "maternities":
          
            LP +=  lpSum( P[g][k][r][a][h] * coeff_rh[(r,h)] for (g,k,r,a,h) in indices_a)\
                                                    -  1e-8*lpSum(Delta_plus[h][l] for h in params_system["H"] for l in params_system["L"])\
                                                    -  1e-8*lpSum(Delta_minus[h][l] for h in params_system["H"] for l in params_system["L"])
        case _:
            LP += (1 - params_system["alpha"]) * lpSum(params_system["c_gk"][g][k] * params_system["D"] * P_gk[g][k] 
                                                        for g in params_system["G"]
                                                        for k in params_system["K_idx"][g]) \
                + params_system["alpha"] * lpSum(P[g][k][r][a][h] * coeff_rh[(r,h)]
                                                 for (g,k,r,a,h) in indices_a)\
                -  1e-6*lpSum(Delta_plus[h][l] for h in params_system["H"] for l in params_system["L"])\
                -  1e-6*lpSum(Delta_minus[h][l] for h in params_system["H"] for l in params_system["L"])





# Define the constraints on the model

# \sum_h P_{g,k,r,a,h} = P_{g,k,r}  for every sub-type (g,k,r,a)
def const_P_gkr(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """
    Define the number of patients assigned to (group, region, and pathway) as the sum of patients across health facilities for every activity.
    """
    for g in params_system["G"]:
        for k in params_system["K_idx"][g]:
            for r in params_system["R"]:
                for a in params_system["A_idx"][g][k]:
                    fast_add(LP, lpSum(vars_system.P[g][k][r][a][h] for h in params_system["H"]) == vars_system.P_gkr[g][k][r])



#\sum_r P_{g,k,r} = P_{g,k} For every g \in G, k \in K_g
def const_P_gk(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None :
    """
    Define the number of patients for (group, pathway) as a sum across all regions
    """
    for g in params_system["G"]:
        for k in params_system["K_idx"][g]: 
            LP += lpSum(vars_system.P_gkr[g][k][r] for r in params_system["R"]) == vars_system.P_gk[g][k]


#\sum_k P_{g,k,r} >= d_{g,r} for every g \in G, r \in R
def def_const_d_gr(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """
    Define the casemix lowerbound constraint
    """
    for g in params_system["G"]:
        for r in params_system["R"]:
            LP +=  lpSum(vars_system.P_gkr[g][k][r] for k in params_system["K_idx"][g]) >= params_system["d_gr"][g][r]

     


# Used for maternity instance
#\sum_k P_{g,k,r} == d_{g,r} for for every g \in G, r \in R
def const_d_gr_maternity(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """
    Define the strict equality demande satisfaction constraint 
    """
    for g in params_system["G"]:
        for r in params_system["R"]:
            LP +=  lpSum(vars_system.P_gkr[g][k][r] for k in params_system["K_idx"][g]) ==  params_system["d_gr"][g][r]



#\sum_k P_{g,k} \geq \underline{q}_g \cdot (\sum_{g2} \sum_{k2} P_{g2,k2}) for every g \in G
def const_lowerbound_qg(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """Define a lower bound on the number of patients in group g
    """
    for g in params_system["G"]:
        LP +=  lpSum(vars_system.P_gk[g][k] for k in params_system["K_idx"][g]) >= \
            params_system["Under_q_g"][g] * lpSum(vars_system.P_gk[g][k]  \
                                              for g in params_system["G"] for k in  params_system["K_idx"][g])



#\sum_k P_{g,k} \leq \overline{q}_g \cdot (\sum_{g'} \sum_{k'} P_{g',k'}) for every g
def const_upperbound_qg(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """Define an upper bound over_q_g on the number of patients of group g
    """
    for g in params_system["G"]:
        LP += lpSum(vars_system.P_gk[g][k] for k in params_system["K_idx"][g]) <= \
              params_system["Over_q_g"][g] * lpSum(vars_system.P_gk[g][k] for g in params_system["G"]for k in params_system["K_idx"][g])




#\sum_{k' \in M_{g,u} P_{g,k'} \geq \underline{q}_{g,u} \cdot (\sum_{k} P_{g,k})
def const_lowerbound_qgu(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """
    Define the lowerbound on the number of patients per group / quality level
    """
    for g in params_system["G"]:
        for u in params_system["U_idx"][g]:
            LP += lpSum(vars_system.P_gk[g][k] for k in params_system["I_gu"][g][u]) >= \
                params_system["Under_q_gu"][g][u] * lpSum(vars_system.P_gk[g][k] \
                                                            for k in params_system["K_idx"][g])



#\sum_{k' \in I_{g,u} P_{g,k'} \leq \overline{q}_{g,u} \cdot (\sum_{k} P_{g,k}) for every g \in G, u \in U
def const_upperbound_qgu(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """
    Define the upperbound on the number of patients per group / quality level
    """
    for g in params_system["G"]:
        for u in params_system["U_idx"][g]:

            LP += lpSum(vars_system.P_gk[g][k] for k in params_system["I_gu"][g][u]) <=\
                  params_system["Over_q_gu"][g][u] * lpSum(vars_system.P_gk[g][k]\
                                                            for k in params_system["K_idx"][g])


#P_{g,k,r,a,h} \leq \sum_{h' \in J_h} P_{g,k,r,a',h'} for every g \in G, k \in K_g, r \in R, a, a' \in A_{g,k}, a \neq a', h \in H
def def_const_J_h(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """ Define constraint on allowed transfers between facilities 
    """
    for g in params_system["G"]:
        for k in params_system["K_idx"][g]:
            for r in params_system["R"]:
                for a1 in params_system["A_idx"][g][k]:
                    if a1 in params_system["N_gka_1"][g][k]:
                        a2 = params_system["N_gka_2"][g][k][a1]
                        for h in params_system["H"]:
                            J2_h = params_system["J_h"][h]
                            fast_add(LP,vars_system.P[g][k][r][a1][h] <= lpSum([vars_system.P[g][k][r][a2][h2] for h2 in J2_h]))


                        
# Q_{g,k,r,a,h} \geq P_{g,k,r,a,h} - P_{g,k,r,a+1,h}
def const_Q_gkrah(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """ Define flow constraint on the patients transfered from one activity to the next
    """
    for g in params_system["G"]:
        for k in params_system["K_idx"][g]:
            N1 = params_system["N_gka_1"][g][k]
            N2 = params_system["N_gka_2"][g][k]
            for r in params_system["R"]:
                P_gkr = vars_system.P[g][k][r]
                Q_gkr = vars_system.Q[g][k][r]
                for a in N1:
                    ap = N2[a]
                    for h in params_system["H"]:
                    # Q >= P_a - P_ap  →  Q - P_a + P_ap >= 0
                        e = LpAffineExpression(
                            [(Q_gkr[a][h], 1.0), (P_gkr[a][h], -1.0), (P_gkr[ap][h], 1.0)],
                            constant=0
                        )
                        fast_add(LP, LpConstraint(e, LpConstraintGE, rhs=0))



#\sum_g \sum_k \sum_r \sum_a \sum_h Q_{g,k,r,a,h} \leq f \cdot (\sum_g \sum_k n_g \cdot P_{g,k})
def const_f(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """
    Define an upper bound on the number of patient transfers allowed
    """

    def lpsum_f(Q, P_gk, n_g, params_system):
        return lpSum(Q[g][k][r][a][h] for g in params_system["G"] for k in params_system["K_idx"][g] \
                    for r in params_system["R"] for a in params_system["N_gka_1"][g][k] for h in params_system["H"]) \
                        <= params_system["p_transf"] * lpSum(P_gk[g][k] * n_g[g] for g in params_system["G"] for k in params_system["K_idx"][g])

    n_g = {k: max(len(x) for x in v.values()) for k, v in params_system["N_gka_1"].items()}
    fast_add(LP,lpsum_f(vars_system.Q, vars_system.P_gk, n_g, params_system))
    

    
def const_m_hl(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """Define the resources constraint as consumption x volume <= capacity m_hl
    """ 
    for h in params_system["H"]:
        for l in params_system["L"]:
            usage = lpSum(
                params_system["D"] * vars_system.P[g][k][r][a][h] * params_system["t_gkal"][g][k][a][l]
                for g in params_system["G"]
                for k in params_system["K_idx"][g]
                for r in params_system["R"]
                for a in params_system["A_idx"][g][k])
            LP += ( usage - vars_system.Delta_plus[h][l] + vars_system.Delta_moins[h][l]
                <= params_system["m_hl"][h][l])


def const_s_hl(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """Define slack constraints to evaluate missing capacity 
    """
    for h in params_system["H"]:
        for l in params_system["L"]:
            LP += lpSum([params_system["D"] * vars_system.P[g][k][r][a][h] * params_system["t_gkal"][g][k][a][l] 
                  for g in params_system["G"] for k in params_system["K_idx"][g] for r in params_system["R"] for a in params_system["A_idx"][g][k]]) \
                    <= params_system["m_hl"][h][l] + vars_system.s_hl[h][l] + vars_system.Delta_plus[h][l] - vars_system.Delta_moins[h][l]



#\sum_h \Delta_{h,l}^{+} = \sum_h \Delta_{h,l}^{-} for every l \in L
def const_delta_zero(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """ Flow constraints on patients transfers
    """
    for l in params_system["L"]:
        LP += lpSum([vars_system.Delta_plus[h][l] for h in params_system["H"]]) - lpSum([vars_system.Delta_moins[h][l] for h in params_system["H"]]) == 0




#\Delta_{h,l}^{+} = P * z_{h,l}^{+} for every h \in H, l \in L
def const_delta_plus(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """
    Define Delta as a function of P
    """
    for h in params_system["H"]:
        for l in params_system["L"]:
            LP += vars_system.Delta_plus[h][l] == params_system["delta_l"][l] *  vars_system.z_hl_plus[h][l]

# \Delta_{h,l}^{-} = P * z_{h,l}^{-} for every h \in H, l \in L
def const_delta_moins(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """
    Define Delta as a function of P
    """
    for h in params_system["H"]:
        for l in params_system["L"]:
            LP += vars_system.Delta_moins[h][l] == params_system["delta_l"][l] * vars_system.z_hl_moins[h][l]

    
#\Delta_{h,l}^{+} \leq b_in{h,l} * m_{h,l} for every h \in H, l \in L
def const_upperbound_resources_transfers_in(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """
    Define upper bound on resource transfer variables (in)
    """
    for h in params_system["H"]:
        for l in params_system["L"]:
            LP += vars_system.Delta_plus[h][l] <= params_system["b_hl_in"][h][l] * params_system["m_hl"][h][l]

#\Delta_{h,l}^{+} \leq b_out{h,l} * m_{h,l} for every h \in H, l \in L
def const_upperbound_resources_transfers_out(LP: pulp.LpProblem, vars_system: dict, params_system: dict) -> None:
    """
    Define upper bound on resource transfer variables (out)
    """
    for h in params_system["H"]:
        for l in params_system["L"]:
            LP += vars_system.Delta_moins[h][l] <= params_system["b_hl_out"][h][l] * params_system["m_hl"][h][l]



# Impose exactly +11.5% of patients (c.f published article)

def def_const_demand(LP, vars_system, params_system):
            LP += lpSum(vars_system.P_gk[g][k]
                        for g in params_system["G"]
                        for k in params_system["K_idx"][g]) * params_system["D"] ==  params_system["D"]*1.115




#################################################
###    CHOOSE SET OF CONSTRAINTS TO INCLUDE   ###
#################################################



def declare_constraints(LP, vars_system, params_system):
    """
    Note: the constraints on the O_gk and positivity of Q are defined on variable initialization
    """

    CONSTRAINTS_DEFAULT=[const_P_gkr, const_P_gk, def_const_d_gr, const_lowerbound_qg, const_upperbound_qg, const_lowerbound_qgu, const_upperbound_qgu,
                           def_const_J_h, const_Q_gkrah, const_f, const_m_hl, const_delta_zero,
                           const_delta_plus, const_delta_moins,  const_upperbound_resources_transfers_in, const_upperbound_resources_transfers_out]
    
    CONSTRAINTS_SLACK=[const_P_gkr, const_P_gk, def_const_d_gr, const_lowerbound_qg, const_upperbound_qg, const_lowerbound_qgu, const_upperbound_qgu,
                           def_const_J_h, const_Q_gkrah, const_f,  const_s_hl, const_delta_zero,
                           const_delta_plus, const_delta_moins,  const_upperbound_resources_transfers_in, const_upperbound_resources_transfers_out, def_const_demand]

    CONSTRAINTS_MATERNITY=[const_P_gkr, const_P_gk, const_d_gr_maternity, const_lowerbound_qg, const_upperbound_qg, const_lowerbound_qgu,
                           def_const_J_h, const_m_hl, const_delta_zero, const_delta_plus, const_delta_moins,
                            const_upperbound_resources_transfers_in, const_upperbound_resources_transfers_out]

    match params_system["mode"]:
        case "slack":
            print("Defining set of constraints with slack capacity.")
            for fn in CONSTRAINTS_SLACK:
                fn(LP, vars_system, params_system)
        case "maternities":
            print("Defining constraints for maternity instance.")
            for fn in CONSTRAINTS_MATERNITY:
                fn(LP, vars_system, params_system)
        case _:
            print("Defining default set of constraints.")
            for fn in CONSTRAINTS_DEFAULT:
                fn(LP, vars_system, params_system)



def preregister_variables(LP, vars_system, params_system):
    """Register all variables into LP once, upfront"""
    all_vars = []

    for g in params_system["G"]:
        for k in params_system["K_idx"][g]:
            all_vars.append(vars_system.P_gk[g][k])
            for r in params_system["R"]:
                all_vars.append(vars_system.P_gkr[g][k][r])
                for a in params_system["A_idx"][g][k]:
                    for h in params_system["H"]:
                        all_vars.append(vars_system.P[g][k][r][a][h])
                        all_vars.append(vars_system.Q[g][k][r][a][h])
    for h in params_system["H"]:
        for l in params_system["L"]:
            all_vars.append(vars_system.Delta_plus[h][l])
            all_vars.append(vars_system.Delta_moins[h][l])
            all_vars.append(vars_system.z_hl_plus[h][l])
            all_vars.append(vars_system.z_hl_moins[h][l])

    LP.addVariables(all_vars)

def fast_add(LP, constraint):
    """Bypass addConstraint's addVariables call — variables pre-registered."""
    name = LP.unusedConstraintName()
    LP.constraints[name] = constraint