import pickle
from biocypher import BioCypher
import pytest
import yaml
import importlib
import logging
import os
import sys
from biolink_model_pydantic.model import BiolinkEntity

logging.basicConfig(level=logging.INFO)

def convert_input_labels(label, replace_char="_"):
    return label.replace(" ", replace_char)

def parse_schema(bcy):
    schema = bcy._get_ontology_mapping()._extend_schema()
    edges_schema = {}
    node_labels = set()

    for k, v in schema.items():
        if v["represented_as"] == "edge": 
            edge_type = convert_input_labels(k)
            source_type = v.get("source", None)
            target_type = v.get("target", None)
            if source_type is not None and target_type is not None:
                if isinstance(v["input_label"], list):
                    label = convert_input_labels(v["input_label"][0])
                    if isinstance(source_type, list):
                        source_type = [convert_input_labels(st) for st in source_type]
                        source_type_lower = [st.lower() for st in source_type]
                    else:
                        source_type = convert_input_labels(source_type)
                        source_type_lower = source_type.lower()
                    
                    if isinstance(target_type, list):
                        target_type = [convert_input_labels(tt) for tt in target_type]
                        target_type_lower = [tt.lower() for tt in target_type]
                    else:
                        target_type = convert_input_labels(target_type)
                        target_type_lower = target_type.lower()
                else:
                    label = convert_input_labels(v["input_label"])
                    if isinstance(source_type, list):
                        source_type = [convert_input_labels(st) for st in source_type]
                        source_type_lower = [st.lower() for st in source_type]
                    else:
                        source_type = convert_input_labels(source_type)
                        source_type_lower = source_type.lower()
                    
                    if isinstance(target_type, list):
                        target_type = [convert_input_labels(tt) for tt in target_type]
                        target_type_lower = [tt.lower() for tt in target_type]
                    else:
                        target_type = convert_input_labels(target_type)
                        target_type_lower = target_type.lower()

                output_label = v.get("output_label", None)
                edges_schema[label.lower()] = {
                    "source": source_type_lower, 
                    "target": target_type_lower, 
                    "output_label": output_label.lower() if output_label is not None else None
                }

        elif v["represented_as"] == "node":
            label = v["input_label"]
            if isinstance(label, list):
                label = label[0]
            label = convert_input_labels(label)
            node_labels.add(label)

    return node_labels, edges_schema
    

@pytest.fixture(scope="session")
def setup_class(request):
    try:
        bcy = BioCypher(
            schema_config_path='config/schema_config.yaml',
            biocypher_config_path='config/biocypher_config.yaml'
        )
        node_labels, edges_schema = parse_schema(bcy) 
    except FileNotFoundError as e:
        pytest.fail(f"Configuration file not found: {e}")
    except yaml.YAMLError as e:
        pytest.fail(f"Error parsing YAML file: {e}")
    except Exception as e:
        pytest.fail(f"Error initializing BioCypher: {e}")
   
    # Load adapters config
    adapters_config_path = request.config.getoption("--adapters-config")
    dbsnp_rsids = request.config.getoption("--dbsnp-rsids")
    dbsnp_pos = request.config.getoption("--dbsnp-pos")
    if dbsnp_rsids:
        logging.info("Loading dbsnp rsids map")
        dbsnp_rsids_dict = pickle.load(open(dbsnp_rsids, 'rb'))
    else:
        logging.warning("--dbsnp-rsids not provided, skipping dbsnp rsids map loading")
        dbsnp_rsids_dict = None
    dbsnp_pos_dict = pickle.load(open(dbsnp_pos, 'rb'))
   
    # Load adapters config
    with open(adapters_config_path, 'r') as f:
        adapters_config = yaml.safe_load(f)

    return node_labels, edges_schema, adapters_config, dbsnp_rsids_dict, dbsnp_pos_dict


def validate_node_type(node_id, node_label, schema_node_labels):
    """
    Validate node type, allowing for multi-role entities in Biolink.
    """
    if isinstance(node_id, tuple):
        node_type = node_id[0].lower()
        return node_type in schema_node_labels or convert_input_labels(node_label) in schema_node_labels
    else:
        return convert_input_labels(node_label) in schema_node_labels


def validate_edge_type_compatibility(source_id, target_id, edge_label, edges_schema):
    """
    Validate if source and target types are compatible with edge schema.
    Accepts lists or single types.
    """
    edge_label_lc = convert_input_labels(edge_label)
    if edge_label_lc not in edges_schema:
        return False, f"Edge label '{edge_label}' not found in schema"

    edge_def = edges_schema[edge_label_lc]
    valid_source_types = edge_def["source"]
    valid_target_types = edge_def["target"]

    source_type = source_id[0].lower() if isinstance(source_id, tuple) else str(source_id).lower()
    target_type = target_id[0].lower() if isinstance(target_id, tuple) else str(target_id).lower()

    if isinstance(valid_source_types, list):
        source_ok = source_type in valid_source_types
    else:
        source_ok = source_type == valid_source_types

    if isinstance(valid_target_types, list):
        target_ok = target_type in valid_target_types
    else:
        target_ok = target_type == valid_target_types

    return source_ok and target_ok, f"Source: {source_type}, Target: {target_type}, Expected sources: {valid_source_types}, Expected targets: {valid_target_types}"

# --- Biolink CURIE → Category mapping ---
CURIE_TO_BIOLINK = {
    "ENSEMBL": {
        "ENSG": "gene",        # ENSEMBL Gene
        "ENST": "transcript",  # ENSEMBL Transcript
    },
    "STATO": "information_content_entity",  # STATO terms map to ICE
}

def resolve_category(node_id, node_label):
    """
    Resolve a node_id or node_label into a Biolink schema category.
    Handles CURIEs like ENSEMBL:ENSG... and STATO:...
    """
    if isinstance(node_id, str) and ":" in node_id:
        prefix, local = node_id.split(":", 1)
        if prefix in CURIE_TO_BIOLINK:
            mapping = CURIE_TO_BIOLINK[prefix]
            if isinstance(mapping, dict):
                for key, category in mapping.items():
                    if local.upper().startswith(key):
                        return category
            else:
                return mapping
    # fallback → normalize label
    return convert_input_labels(node_label).lower()


@pytest.mark.filterwarnings("ignore")
class TestBiocypherKG:
        
        def test_adapter_nodes_in_schema(adapter, schema_view):
            """
            Validate nodes using Biolink schema categories, not old class/property split.
            """
            for node in adapter.get_nodes():
                node_id = node["id"]
                node_cat = node.get("category")

                # Ensure category is known in Biolink
                if node_cat not in schema_view.all_classes():
                    pytest.fail(f"Node {node_id} has unknown category {node_cat}")

                # Skip property vs. class conflict check (Biolink allows reuse)
                # Only validate category mapping
                expected = schema_view.get_class(node_cat)
                assert expected is not None, f"Node {node_id} maps to invalid category {node_cat}"


        def test_adapter_edges_in_schema(adapter, schema_view):
            """
            Validate edges by resolving CURIEs to Biolink categories, instead of hardcoding.
            """
            for edge in adapter.get_edges():
                edge_pred = edge["predicate"]
                src_id = edge["subject"]
                tgt_id = edge["object"]

                # Get expected domain and range categories from Biolink
                pred = schema_view.get_slot(edge_pred)
                if not pred:
                    pytest.fail(f"Predicate {edge_pred} not found in Biolink schema")

                expected_sources = schema_view.get_class(pred.domain).mixins if pred.domain else []
                expected_targets = schema_view.get_class(pred.range).mixins if pred.range else []

                # Resolve subject/object categories dynamically
                src_cat = adapter.curie_to_category(src_id)
                tgt_cat = adapter.curie_to_category(tgt_id)

                assert src_cat in expected_sources or not expected_sources, \
                    f"Edge {edge_pred} source {src_id} not in {expected_sources}"

                assert tgt_cat in expected_targets or not expected_targets, \
                    f"Edge {edge_pred} target {tgt_id} not in {expected_targets}"
