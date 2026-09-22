import pulp
from pathlib import Path

def read_inputs(file_params_system):
    import json
    with open(file_params_system, "r") as f:
        params_system = json.load(f)

    return params_system

def get_var(curr_var, row, list_dims):
    v = curr_var[row["group"]][row["pathway"]]
    if "region" in list_dims:
        v = v[row["region"]]
    if "activity" in list_dims:
        v = v[row["activity"]]
    if "facility" in list_dims:
        v = v[row["facility"]]
    if "resource" in list_dims:
        v = v[row["resource"]]
    return pulp.value(v)


def hl_to_records(var_hl, params_system):
    """For variables only indexed by h (facilities) and l (resources) (e.g Delta, z, s)"""
    return [
        {"facility": h, "resource": l, "value": pulp.value(var_hl[h][l])}
        for h in params_system["H"]
        for l in params_system["L"]
    ]

def vars_to_records(curr_var, list_dims, params_system):
    records = []
    for g in params_system["G"]:
        for k in params_system["K_idx"][g]:
            for r in (params_system["R"] if "region" in list_dims else [None]):
                for a in (params_system["A_idx"][g][k] if "activity" in list_dims else [None]):
                    for h in (params_system["H"] if "facility" in list_dims else [None]):
                        for l in (params_system["L"] if "resource" in list_dims else [None]):
                                row = { "group": g, "pathway": k}
                                if "region" in list_dims: row["region"] = r
                                if "activity" in list_dims: row["activity"] = a
                                if "facility" in list_dims: row["facility"] = h
                                if "resource" in list_dims: row["resource"]= l
                                row["value"] = get_var(curr_var, row, list_dims)
                                records.append(row)
    return records


def package_results(vars_system, params_system):
    dict_results = {
        "P_gkrah": vars_to_records(vars_system.P, ["group","pathway","region","activity","facility"], params_system),
        "Q_gkrah": vars_to_records(vars_system.Q, ["group","pathway","region","activity","facility"], params_system),
        "P_gkr": vars_to_records(vars_system.P_gkr, ["group","pathway","region"], params_system),
        "P_gk": vars_to_records(vars_system.P_gk, ["group","pathway"], params_system),
        "Delta_plus" : hl_to_records(vars_system.Delta_plus, params_system),
        "Delta_moins": hl_to_records(vars_system.Delta_moins, params_system),
        "z_hl_plus": hl_to_records(vars_system.z_hl_plus, params_system),
        "z_hl_moins": hl_to_records(vars_system.z_hl_moins, params_system),
        "s_hl": hl_to_records(vars_system.s_hl, params_system)
    }
    return dict_results



def get_summary_results(dict_results, params_system, objective_value, output_path=None):
    """Returns summary file of results"""
    from src.core.utils.data_utils_Burdett import get_patients_blocking
    import pandas as pd
    import json
    P = pd.DataFrame(dict_results["P_gkrah"]).groupby(["group", "pathway", "region","activity", "facility"])["value"].sum().to_dict()
    
    Delta_plus = pd.DataFrame(dict_results["Delta_plus"]).groupby(["facility", "resource"])["value"].sum().to_dict()
    Delta_moins = pd.DataFrame(dict_results["Delta_moins"]).groupby(["facility", "resource"])["value"].sum().to_dict()
    P_gk = pd.DataFrame(dict_results["P_gk"]).groupby(["group", "pathway"])["value"].sum().to_dict()
    n_patients = get_n_patients(params_system, P_gk)
    
    l_usage = get_resources_usage(params_system, P, Delta_plus=Delta_plus, Delta_moins=Delta_moins) 
    k_usage = get_pathways_usage(params_system, P_gk)
    try:
        spi = get_SPI(params_system, P)
        qci = get_QCI(params_system, P_gk)  
    except Exception as e:
        spi = None
        qci = None
    try:
        # Only relevant for the Burdett dataset
        p_blocking = get_patients_blocking(params_system,P)
    except Exception as e:
        p_blocking = None

    results = {"obj":objective_value, "n_patients": n_patients, "spi": spi, "qci": qci, "resources_usage":l_usage, "pathway_distribution": k_usage, "saturation": p_blocking}
    
    if output_path:
        with open(output_path, "w") as fp:
            json.dump(results, fp)
    return results



def get_n_patients(params_system: dict, P_gk: dict):
    """Returns the total number of patients treated in the system"""
    n_patients = sum(params_system["D"] * P_gk[(g, k)] for g in params_system["G"]
                                                       for k in params_system["K_idx"][g])
    return n_patients


def get_SPI(params_system: dict, P: dict) -> float:
    """Returns a Spatial Proximity Index metric"""
    nominator = sum(params_system["D"] * params_system["w_rh"][r][h] * P[(g,k,r,"ANES",h)]
                    for g in params_system["G"]
                    for k in params_system["K_idx"][g]
                    for r in params_system["R"]
                    for h in params_system["H"] if h not in ["DOM", "ORTHf"])
    denominator = sum(params_system["D"] * P[(g,k,r,"ANES",h)]
                      for g in params_system["G"]
                      for k in params_system["K_idx"][g]
                      for r in params_system["R"]
                      for h in params_system["H"] if h not in ["DOM", "ORTHf"])
    return  nominator / denominator if denominator != 0 else 0

def get_QCI(params_system: dict, P_gk: dict) -> float:
    """Return a Quality of Care Index metric"""
    nominator = sum(params_system["D"] * params_system["c_gk"][g][k] * P_gk[(g,k)] 
                    for g in params_system["G"]
                    for k in params_system["K_idx"][g])
    denominator = sum(params_system["D"] * P_gk[(g,k)] 
                    for g in params_system["G"]
                    for k in params_system["K_idx"][g])

    return  nominator / denominator if denominator != 0 else 0
 


def get_resources_usage(params_system: dict, P: dict, Delta_plus: dict, Delta_moins: dict) -> dict:
    """Returns a dict with the usage percentage of each resource type"""
    utilisation = {l: 0 for l in params_system["L"]}
    for l in params_system["L"]:
        nominator = 0
        denominator = 0
        for g in params_system["G"]:
            for k in params_system["K_idx"][g]:
                for r in params_system["R"]:
                    for a in params_system["A_idx"][g][k]:
                        for h in params_system["H"]:
                            nominator += params_system["t_gkal"][g][k][a][l] * P[(g, k, r, a, h)] * params_system["D"]
        for h in params_system["H"]:
          
            denominator += params_system["m_hl"][h][l] + Delta_plus[h,l] - Delta_moins[h,l]
        utilisation[l] = nominator / denominator if denominator != 0 else 0
    return utilisation


def get_pathways_usage(params_system: dict, P_gk: dict) -> dict:
    """Returns the percentage of patients for each pathway"""
    n_total = sum(params_system["D"] * P_gk[(g, k)] for g in params_system["G"]
                                                       for k in params_system["K_idx"][g])
    all_pathways = list({k for g in params_system["G"] for k in params_system["K_idx"][g]})
    utilisation = {k: 0 for k in all_pathways}
    for k in all_pathways:
        n_k = 0
        for g in params_system["G"]:
            if k in params_system["K_idx"][g]:
                n_k += params_system["D"] * P_gk[(g, k)] 
        utilisation[k] = round(n_k / n_total,3) if n_total != 0 else 0
    return utilisation



def read_metadata(inputfile : str | Path):
    """
    Reads metadata from a JSON file.
    """
    import json

    with open(inputfile, "r") as f:
        metadata = json.load(f)

    return metadata


def read_configs(config_category, config_path="backend/config.yaml"):
    import yaml

    with open(config_path, "r") as f:
        configs = yaml.safe_load(f)

    return configs.get(config_category)



def read_geojson_projected(filename: str | Path):
    import geopandas as gpd
    gdf = gpd.read_file(Path(filename))
    if gdf.empty:
        raise ValueError(f"GeoJSON file {filename} contains no features.")
    gdf = gdf.to_crs(epsg = 2154)
    return gdf

    

def get_distance_to_dep(dep_name : str, coords: list) -> int:
    """ Returns the distance in kms between a french departments and a point [lon, lat]"""
    from shapely.geometry import Point
    import geopandas as gpd

    geo_deps = read_geojson_projected("/data/departements.geojson")
    try:
        dep_geo = geo_deps.loc[geo_deps["nom"] == dep_name, "geometry"].iloc[0]
    
    except IndexError:
        raise ValueError(f"Departments {dep_name} not found")

    gdf_point = gpd.GeoSeries(Point(coords), crs="EPSG:4326").to_crs(geo_deps.crs)
        
    distance_m = gdf_point.iloc[0].distance(dep_geo)
    distance_km = distance_m / 1000
    
    return distance_km

def get_department_coords(dep_name: str, dep_centroids):
    return dep_centroids[dep_name]

