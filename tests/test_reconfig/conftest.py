import os 
import sys
import pytest
import shutil
import pathlib
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from tests.test_reconfig import benchmark_compare

output_dir = pathlib.Path.cwd() / "ExampleData" / "OrthoFinder" 
# if os.path.exists(output_dir):
#     shutil.rmtree(output_dir)  

for entry in os.scandir(path):
    print(entry.stat().st_mtime)


"""
Before running this test, please make sure the test data is obtained with the `-fn`.
"""

@pytest.fixture(scope="module")
def orthogroup_stats():
    orthostats = benchmark_compare.orthogroups_analysis("ExampleData")
    return orthostats

@pytest.fixture(scope="module")
def unassigned_gens_stats():
    orthostats = benchmark_compare.unassigned_gens_analysis("ExampleData")
    return orthostats

@pytest.fixture(scope="module")
def orthologues_stats():
    orthostats = benchmark_compare.orthologues_analysis("ExampleData")
    return orthostats