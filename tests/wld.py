from weighted_levenshtein import LD, WLD

import numpy as np



# list of all string chains
test_data = ["abc", "ab", "Ec", "b*c", "def", "a9"]

# string union
test_strings = ["abcdef"]

# similarity matrix
test_sim_matrix = np.array(
    [[6, 1, 5, 1, 4, 0],
     [1, 6, 2, 5, 2, 0],
     [5, 2, 5, 2, 3, 0],
     [1, 5, 2, 7, 2, 1],
     [4, 2, 3, 2, 7, 0],
     [0, 0, 0, 1, 0, 10]], dtype=np.uint8)

# distance matrix
# test_sim_matrix.max()-test_sim_matrix
test_dist_matrix = np.array(
    [[4, 9, 5, 9, 6, 10],
     [9, 4, 8, 5, 8, 10],
     [5, 8, 5, 8, 7, 10],
     [9, 5, 8, 3, 8, 9],
     [6, 8, 7, 8, 3, 10],
     [10, 10, 10, 9, 10, 0]], dtype=np.uint8)


# test_functions = [levenshtein, weighted_levenshtein]

# score of insertion and deletion
# test_indel_score = [median, maximum]


# Levenshtein distance tests
# def test_mono_levenshtein():
#    string1 = "abc"
#    list_string2 = ["", "aB", "Ec", "b*c", "def"]
#    list_output = [3, 1, 2, 2, 3]
#    for string2, output in zip(list_string2, list_output):
#        assert monoLD(string1, string2) == output

# def test_pairwise_levenshtein():
#    a = LD(["fB", "abc", "def"])
#    list_output = [['fB', 'abc', 2],
#                   ['fB', 'def', 3],
#                   ['abc', 'def', 3]]
#    assert a.distance == list_output


def test_init_ok():
    LD(test_data)


def test_init_input_is_list():
    try:
        LD("abc")
        assert False
    except TypeError:
        assert True


def test_init_input_greater_than_1():
    try:
        LD(["abc"])
        assert False
    except ValueError:
        assert True


def test_init_input_all_str():
    try:
        LD(["abc", 0])
        assert False
    except TypeError:
        assert True


def test_init_input_no_empty():
    try:
        LD(["abc", ""])
        assert False
    except ValueError:
        assert True


def test_init_input_no_whitespace():
    try:
        LD(["abc", " ", "   "])
        assert False
    except ValueError:
        assert True


def test_sort_init_input():
    ld = LD(['1', 'Zq', 'Zq'])
    assert ld.words == ['1', 'Zq']


def test_get_unique_values_of_all_words():
    ld = LD(['fs', 'sdf', '1', '+', 'Zq'])
    assert ld.unique() == ['+', '1', 'd', 'f', 'q', 's', 'Z']


def test_get_unique_pairwises_of_words():
    ld = LD(['abc', '1', 'Zq', 'abc', 'Zq'])
    assert ld.pairwise() == [('1', 'abc'), ('1', 'Zq'), ('abc', 'Zq')]

