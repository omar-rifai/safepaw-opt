


def check_resource_saturation(params_system:dict, P:dict)-> dict:
    """Retruns a dict with the least saturated ward level for group g, pathway k and resource l"""
    import numpy as np
    l_map = {"cap_OT":[x for x in params_system["H"] if x.split("_")[1] == "OT"],
             "cap_ICU":[x for x in params_system["H"] if x.split("_")[1] == "ICU"],
             "cap_Ward":[x for x in params_system["H"] if x.split("_")[1] not in  ["ICU", "OT"]]}
    
    saturation = {l: {h: np.nan for h in l_map[l]} for l in params_system["L"]}
    for l in params_system["L"]:
        for h in params_system["H"]:
            if h in l_map[l]:
                consumption = 0
                capacity = 0
                for g in params_system["G"]:
                    for k in params_system["K_idx"][g]:
                        for r in params_system["R"]:
                            for a in params_system["A_idx"][g][k]:
                                if h in params_system["O_gk"][g][k]:
                                    consumption += params_system["t_gkal"][g][k][a][l] * P[(g, k, r, a, h)] * params_system["D"]        
                capacity = params_system["m_hl"][h][l]
                if capacity > 0:
                    saturation[l][h] = round(consumption / capacity,5)

    min_lgk = {l: {g :{k: 0 for k in params_system["K_idx"][g]} for g in params_system["G"]} for l in params_system["L"]}
    min_lg = {l: {g : 0 for g in params_system["G"]} for l in params_system["L"]}
    
    for l in params_system["L"]:
        for g in params_system["G"]:
            for k in params_system["K_idx"][g]:
                curr_min = 1
                for h in params_system["O_gk"][g][k]:
                    if h in l_map[l]:
                        if saturation[l][h] < curr_min:
                            curr_min = saturation[l][h]
                min_lgk [l][g][k] = curr_min
            min_lg[l][g] = min(min_lgk[l][g].values()) 
    return min_lgk


def get_u_from_k(g, k, I_gu):
    """Retrienve the quality level `u` which is related to the pathway `k`"""
    for u in I_gu[g]:
        if k in I_gu[g][u]:
            return u

def get_patients_blocking(params,P):
    """Get percentage of patients in saturated pathways"""
    saturation= check_resource_saturation(params,P)
    list_gk_saturated = {g:[] for g in params["G"]}
    for l in params["L"]:
        for g in params["G"]:
            for k in params["K_idx"][g]:
                if saturation[l][g][k] == 1:
                    list_gk_saturated[g].append(k)

    p_saturated = 0
    for g in list_gk_saturated:
        for k in list_gk_saturated[g]:
            u = get_u_from_k(g,k,params["I_gu"])
            p_saturated += params["Under_q_gu"][g][u]*params["Under_q_g"][g]

    return p_saturated
