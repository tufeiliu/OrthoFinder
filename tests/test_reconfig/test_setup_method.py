import pytest
import os 
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from dev_tools import orthogroups_compare


@pytest.fixture(scope="module")
def orthogroup_stats():
    orthostats = orthogroups_compare.orthogroups_analysis("update_vs_orthofinderv2")
    return orthostats


def test_orthogroups_left_diff(orthogroup_stats):
    left_diff_num = orthogroup_stats["left_diff"]
    left_label, rihgt_label = orthogroup_stats["labels"]
    assert left_diff_num == 0, f"The difference between {left_label} and {rihgt_label} is not empty!"


def test_orthogroups_right_diff(orthogroup_stats):
    right_diff_num = orthogroup_stats["right_diff"]
    left_label, rihgt_label = orthogroup_stats["labels"]
    assert right_diff_num == 0, f"The difference between {right_label} and {left_label} is not empty!"

def test_orthogroups_left_intersection(orthogroup_stats):
    intersection_num = orthogroup_stats["intersection"]
    left_diff_num = orthogroup_stats["left_diff"]
    left_num = intersection_num + left_diff_num
    left_label, _ = orthogroup_stats["labels"]
    assert intersection_num == left_num, f"The intersection is not equal to {left_label}!"


def test_orthogroups_right_intersection(orthogroup_stats):
    intersection_num = orthogroup_stats["intersection"]
    right_diff_num = orthogroup_stats["right_diff"]
    right_num = intersection_num + right_diff_num
    _, rihgt_label = orthogroup_stats["labels"]
    assert intersection_num == right_num, f"The intersection is not equal to {right_label}!"

def test_orthogroups_similarity(orthogroup_stats):
    similarity = orthogroup_stats["similarity"]
    
    assert (similarity > 99.0 and similarity <= 100.0), f"The similarity between {left_label} and {right_label} is outside of the range (99.0, 100.0)!"
