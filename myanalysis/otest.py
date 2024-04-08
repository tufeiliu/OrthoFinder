import pathlib
import os 
import pandas as pd 
import numpy as np 
from collections import deque, defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import functools
import itertools


# def compare_diff(left_sets, right_sets):

#     distinct_orthogroups = set()

#     if left_sets == right_sets:
#         return left_sets, 0, 0

#     left_diffs = left_sets - right_sets 
#     right_diffs = right_sets - left_sets
    
#     occurred = set()
#     gen_intersection_dict = defaultdict(lambda: defaultdict(dict))

#     for left_diff in left_diffs:

#         max_intersection = 0
#         left_gens = len(", ".join(left_diff).split(", "))
 
#         for right_diff in right_diffs:
#             right_gens = len(", ".join(right_diff).split(", "))
            
#             for left, right in itertools.product(left_diff, right_diff):

#                 left_com = left.count(", ")
#                 right_com = right.count(", ")

#                 if left_com > 0:
#                     left_set = set(left.split(": ")[1].split(", "))
                    
#                 else:
#                     left_set = {left.split(": ")[1]}

#                 if right_com > 0:
#                     right_set = set(right.split(": ")[1].split(", "))
#                 else:
#                     right_set = {right.split(": ")[1]}
                
#                 intersection = left_set & right_set 
                
#                 if len(intersection) > 0:
#                     print(left_diff, right_diff)
#                     occurred.add(right_diff)
        
#                     max_intersection = np.amax((max_intersection, len(intersection)))
#                     max_intersection_union = max_intersection / (left_gens + right_gens - len(intersection))
#                     max_intersection_left = max_intersection / left_gens
#                     max_intersection_right = max_intersection / right_gens

#         if max_intersection != 0:
#             gen_intersection_dict[left_diff]["max_intersection"] = max_intersection
#             gen_intersection_dict[left_diff]["max_intersection_divide_left"] = max_intersection_left
#             gen_intersection_dict[left_diff]["max_intersection_divide_right"] = max_intersection_right 
#             gen_intersection_dict[left_diff]["max_intersection_divide_union"] = max_intersection_union

#     distinct_orthogroups = left_sets | (right_diffs - occurred)       
#     distinct_right, orthogroup_overlap = len(right_diffs - occurred), len(occurred)

#     return distinct_orthogroups, distinct_right, orthogroup_overlap, gen_intersection_dict

def compare_diff(left_sets, right_sets):

    distinct_orthogroups = set()

    if left_sets == right_sets:
        return left_sets, 0, 0

    left_diffs = left_sets - right_sets 
    right_diffs = right_sets - left_sets
    
    occurred = set()
    gen_intersection_dict = defaultdict(lambda: defaultdict(dict))
    prev = set()

    for left, right in itertools.product(left_diffs, right_diffs):

        left_gens = len(left)
        right_gens = len(right)

        intersection = left & right
        print(f"intersection {intersection}")

        if len(intersection) > 0:
            occurred.add(right)

        if left not in prev:
            max_intersection = 0
            max_intersection_union = 0
            max_intersection_left = 0
            max_intersection_right = 0
            prev.add(left)

        max_intersection = np.amax((max_intersection, len(intersection)))
        max_intersection_union = np.amax((max_intersection_union, max_intersection / (left_gens + right_gens - len(intersection))))
        max_intersection_left = np.amax((max_intersection_left, max_intersection / left_gens))
        max_intersection_right = np.amax((max_intersection_right, max_intersection / right_gens))


        gen_intersection_dict[left]["max_intersection"] = max_intersection
        gen_intersection_dict[left]["max_intersection_divide_union"] = max_intersection_union
        gen_intersection_dict[left]["max_intersection_divide_left"]  = max_intersection_left
        gen_intersection_dict[left]["max_intersection_divide_right"]  = max_intersection_right
    
    distinct_orthogroups = left_sets | (right_diffs - occurred)       
    distinct_right, orthogroup_overlap = right_diffs - occurred, occurred

    return distinct_orthogroups, distinct_right, orthogroup_overlap, gen_intersection_dict




        # for left, right in itertools.product(left_diff, right_diff):
        #     print(left, right)

    #         left_com = left.count(", ")
    #         right_com = right.count(", ")

    #         if left_com > 0:
    #             left_set = set(left.split(": ")[1].split(", "))
                
    #         else:
    #             left_set = {left.split(": ")[1]}

    #         if right_com > 0:
    #             right_set = set(right.split(": ")[1].split(", "))
    #         else:
    #             right_set = {right.split(": ")[1]}
            
    #         intersection = left_set & right_set 
            
    #         if len(intersection) > 0:
    #             print(left_diff, right_diff)
    #             occurred.add(right_diff)
    
    #             max_intersection = np.amax((max_intersection, len(intersection)))
    #             max_intersection_union = max_intersection / (left_gens + right_gens - len(intersection))
    #             max_intersection_left = max_intersection / left_gens
    #             max_intersection_right = max_intersection / right_gens

    # if max_intersection != 0:
    #     gen_intersection_dict[left_diff]["max_intersection"] = max_intersection
    #     gen_intersection_dict[left_diff]["max_intersection_divide_left"] = max_intersection_left
    #     gen_intersection_dict[left_diff]["max_intersection_divide_right"] = max_intersection_right 
    #     gen_intersection_dict[left_diff]["max_intersection_divide_union"] = max_intersection_union

    # distinct_orthogroups = left_sets | (right_diffs - occurred)       
    # distinct_right, orthogroup_overlap = len(right_diffs - occurred), len(occurred)

    # return distinct_orthogroups, distinct_right, orthogroup_overlap, gen_intersection_dict

# left_sets = {frozenset({"A: a1, a2, a4", "B: b1, b2", "C: c1"}), 
#                frozenset({"A: a5, a6", "B: b3", "C: c2, c3"}),
#                frozenset({"A: a7, a8", "B: b4, b9", "C: c6, c7"})}


# right_sets = {frozenset({"A: a3", "B: b3, b5"}),
#                 frozenset({"A: a7", "B: b6", "C: c4, c5"})}


# distinct_orthogroups, distinct_right, orthogroup_overlap, gen_intersection_dict = compare_diff(left_sets, right_sets)

# print(distinct_right, orthogroup_overlap)

left_sets = {frozenset({("A", "a1"), ("A", "a2"), ("A", "a4"), ("B", "b1"), ("B", "b2"), ("C", "c1")}), 
               frozenset({("A", "a5"), ("A", "a6"), ("B", "b3"), ("C", "c2"), ("C", "c3")}),
               frozenset({("A", "a7"), ("A", "a8"), ("B", "b4"), ("B", "b9"), ("C", "c6"), ("C", "c7")})}


right_sets = {frozenset({("A", "a3"), ("B", "b3"), ("B", "b5")}),
                frozenset({("A", "a7"), ("B", "b6"), ("C", "c4"), ("C", "c5")}),
               frozenset({("A", "a10"), ("B", "b10"), ("C", "c10")}) }

compare_diff(left_sets, right_sets)

# print(len(gen_intersection_dict))
# for key, val in gen_intersection_dict.items():
#     print(key)
#     print(len(val))
#     print(val)
# print(gen_intersection_dict)
# right_diff = {frozenset({"A: a1, a2, a4", "B: b1, b2", "C: c1"}), 
#                frozenset({"A: a5, a6", "B: b3", "C: c2, c3"})}
# print(distinct_orthogroups)
# print(all_orthogroups)


A = 'a7, a8, b4, b9, c6, c7'
B = 'a7, b6, c4, c5'

C = 'a5, a6, b3, c2, c3'
D = 'a3, b3, b5'