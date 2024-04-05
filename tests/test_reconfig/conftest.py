import os 
import sys
import pytest
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from dev_tools import orthogroups_compare


@pytest.fixture(scope="module")
def orthogroup_stats():
    orthostats = orthogroups_compare.orthogroups_analysis("ExampleData")
    return orthostats

@pytest.fixture(scope="module")
def unassigned_gens_stats():
    orthostats = orthogroups_compare.unassigned_gens_analysis("ExampleData")
    return orthostats

@pytest.fixture(scope="module")
def orthologues_stats():
    orthostats = orthogroups_compare.orthologues_analysis("ExampleData")
    return orthostats