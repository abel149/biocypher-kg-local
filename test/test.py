import pickle
from biocypher import BioCypher
import pytest
import yaml
import importlib
import logging
import os
import sys


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

    def test_adapter_nodes_in_schema(self, setup_class):
        """
        Verify that node labels generated by adapters are included
        within the Biolink-based schema.
        """
        node_labels, edges_schema, adapters_config, dbsnp_rsids_dict, dbsnp_pos_dict = setup_class
        for adapter_name, config in adapters_config.items():
            if config["nodes"]:
                adapter_module = importlib.import_module(config['adapter']['module'])
                adapter_class = getattr(adapter_module, config['adapter']['cls'])

                adapter_args = config['adapter']['args'].copy()
                if "dbsnp_rsid_map" in adapter_args:
                    adapter_args["dbsnp_rsid_map"] = dbsnp_rsids_dict
                if "dbsnp_pos_map" in adapter_args:
                    adapter_args["dbsnp_pos_map"] = dbsnp_pos_dict
                adapter_args['write_properties'] = True
                adapter_args['add_provenance'] = True

                adapter = adapter_class(**adapter_args)

                sample_node = next(adapter.get_nodes(), None)
                assert sample_node, f"No nodes found for adapter '{adapter_name}'"

                node_id, node_label, node_props = sample_node

                category = resolve_category(node_id, node_label)
                assert category in {lbl.lower() for lbl in node_labels}, \
                    f"Node '{node_id}' ({category}) from adapter '{adapter_name}' not in schema"


    def test_adapter_edges_in_schema(self, setup_class):
        """
        Verify that edge labels and source/target categories
        from adapters are compatible with Biolink schema.
        """
        node_labels, edges_schema, adapters_config, dbsnp_rsids_dict, dbsnp_pos_dict = setup_class
        for adapter_name, config in adapters_config.items():
            if config['edges']:
                adapter_module = importlib.import_module(config['adapter']['module'])
                adapter_class = getattr(adapter_module, config['adapter']['cls'])

                adapter_args = config['adapter']['args'].copy()
                if "dbsnp_rsid_map" in adapter_args:
                    adapter_args["dbsnp_rsid_map"] = dbsnp_rsids_dict
                if "dbsnp_pos_map" in adapter_args:
                    adapter_args["dbsnp_pos_map"] = dbsnp_pos_dict
                adapter_args['write_properties'] = True
                adapter_args['add_provenance'] = True

                adapter = adapter_class(**adapter_args)

                sample_edge = next(adapter.get_edges(), None)
                assert sample_edge, f"No edges found for adapter '{adapter_name}'"

                source_id, target_id, edge_label, edge_props = sample_edge
                assert edge_label.lower() in edges_schema, \
                    f"Edge label '{edge_label}' from adapter '{adapter_name}' not in schema"

                source_cat = resolve_category(source_id, source_id)
                target_cat = resolve_category(target_id, target_id)

                edge_def = edges_schema[edge_label.lower()]
                expected_sources = edge_def.get("source", [])
                expected_targets = edge_def.get("target", [])

                assert source_cat in expected_sources, \
                    f"Edge '{edge_label}' from adapter '{adapter_name}': " \
                    f"Source {source_cat} not in {expected_sources}"
                assert target_cat in expected_targets, \
                    f"Edge '{edge_label}' from adapter '{adapter_name}': " \
                    f"Target {target_cat} not in {expected_targets}"
