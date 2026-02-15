import requests
import json
import os

def search_refined_pdb(resolution=2.5, chemical_id=None):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    nodes = [
        {"type": "terminal", "service": "text", "parameters": {"attribute": "struct_keywords.pdbx_keywords", "operator": "contains_phrase", "value": "ANTIBODY"}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.resolution_combined", "operator": "less_or_equal", "value": resolution}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.nonpolymer_bound_components", "operator": "contains_phrase", "value": "saccharide"}}
    ]
    if chemical_id:
        nodes.append({"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entity_nonpolymer_feature_summary.comp_id", "operator": "exact_match", "value": chemical_id}})
    
    query = {"query": {"type": "group", "logical_operator": "and", "nodes": nodes}, "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": 100}}}
    response = requests.post(url, json=query)
    return [item['identifier'] for item in response.json().get('result_set', [])] if response.status_code == 200 else []

def get_smiles_from_pdb(het_id):
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{het_id}"
    res = requests.get(url)
    if res.status_code == 200:
        return res.json().get("rcsb_chem_comp_descriptor", {}).get("smilesiso", "Unknown")
    return "Unknown"