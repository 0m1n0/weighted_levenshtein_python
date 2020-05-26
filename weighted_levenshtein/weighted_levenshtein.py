import itertools
from typing import List
from tqdm import tqdm
import numpy as np


class WLD:
    """Weighted Levenshtein distance
    """

    class WeightType:
        DEFAULT = 0
        MIN = 1
        MAX = 2
        MEAN = 3
        MEDIAN = 4
        CUSTOM = 5

    weight_substitution_default = 1

    def __init__(self,
                 words: List[str],
                 string_list: List[str] = None,
                 weight_substitution_matrix: np.ndarray = None,
                 weight_insertion_type: WeightType = WeightType.DEFAULT,
                 weight_deletion_type: WeightType = WeightType.DEFAULT,
                 custom_weight_insertion: int = 1,
                 custom_weight_deletion: int = 1):
        """

        Args:
            words (list): list of words in string.
        """
        if type(words) != list:
            raise TypeError("Must be a list.")
        if len(words) <= 1:
            raise ValueError("Must contain at least 2 words.")
        for word in words:
            if type(word) != str:
                raise TypeError("All words must be string.")
            elif not word:
                raise ValueError("Empty string.")
            elif word.isspace():
                raise ValueError("Only whitespace string.")

        # remove " " and sort by unique values
        no_empty = [x for x in words if not x.isspace()]
        self.words = sorted(set(no_empty), key=str.lower)

        # substitution
        if weight_substitution_matrix is None:
            self.strings = self.unique()
            # substitution matrix filled by default value
            self.sub_matrix = np.full((len(self.strings), len(self.strings)),
                                      self.weight_substitution_default, dtype=np.uint8)
            # diagonal = 0
            np.fill_diagonal(self.sub_matrix, 0)
        else:
            self.strings = string_list
            self.sub_matrix = weight_substitution_matrix

        print(self.strings)
        print(self.sub_matrix)

        # insertion
        if weight_insertion_type == WLD.WeightType.CUSTOM:
            self.w_ins = custom_weight_insertion
        else:
            self.w_ins = self.sub_matrix.max()

        # deletion
        if weight_deletion_type == WLD.WeightType.CUSTOM:
            self.w_del = custom_weight_deletion
        else:
            self.w_del = self.sub_matrix.max()

    def unique(self):
        """

        Returns:
            list: unique values of all words, one by one.

        """
        merged_words = ''.join(self.words)
        list_unique = sorted(list(set(merged_words)), key=str.lower)
        return list_unique

    def pairwise(self):
        """

        Returns:
            list: all unique pairwise combinations of word list.

        """
        list_pairwise = list(itertools.combinations(self.words, 2))
        return list_pairwise

    def mono_levenshtein(self, str1, str2):
        print("\n", str1, str2)

        # add empty row and column
        n_row = len(str1) + 1
        n_col = len(str2) + 1

        # generate zero matrix
        matrix = np.zeros((n_row, n_col), dtype=np.uint8)

        # first row and column, start from zero and indexed increasingly by insertion/deletion value
        for i in range(n_row):
            matrix[i, 0] = i * self.w_ins
        for j in range(n_col):
            matrix[0, j] = j * self.w_del

        # two loops: row-wise and column-wise
        for i in range(1, n_row):
            for j in range(1, n_col):
                # fill new value at position matrix[i, j] = x11
                # x11 is calculated from the three surrounding values (x00, x01, x10) previously completed
                # |-----------|
                # | x00 | x01 |
                # | x10 | x11 |
                # |-----------|
                x00 = matrix[i - 1, j - 1]  # substitution
                x01 = matrix[i, j - 1]  # insertion
                x10 = matrix[i - 1, j]  # deletion

                # get substitution weight for str1 and str2
                w_sub = self.sub_matrix[self.strings.index(str1[i - 1]), self.strings.index(str2[j - 1])]
                # # if two letters are identical, the substitution costs 0 otherwise a positive value
                # w_sub = (str1[i - 1] != str2[j - 1]) * w_sub

                matrix[i, j] = min(x00 + w_sub, x01 + self.w_ins, x10 + self.w_del)

        print(matrix)
        distance = matrix[n_row - 1, n_col - 1]
        return distance

    def levenshtein(self):
        all_list_distance = []
        for str1, str2 in tqdm(self.pairwise()):
            # list_distance = [str1, str2, self.mono_levenshtein(str1, str2)]
            # all_list_distance.append(list_distance)
            yield str1, str2, self.mono_levenshtein(str1, str2)


dist_matrix = np.array(
    [[4, 9, 5, 9, 6, 10],
     [9, 4, 8, 5, 8, 10],
     [5, 8, 5, 8, 7, 10],
     [9, 5, 8, 3, 8, 9],
     [6, 8, 7, 8, 3, 10],
     [10, 10, 10, 9, 10, 0]], dtype=np.uint8)

# a = WLD(['abc', 'ab', 'ec'],
#         string_list=['a', 'b', 'c', 'd', 'e'],
#         weight_substitution_matrix=dist_matrix)
a = WLD(['abc', 'def'],
        string_list=['a', 'b', 'c', 'd', 'e', 'f'],
        weight_substitution_matrix=dist_matrix)
print("deletion", a.w_del)
print("insertion", a.w_ins)

print(list(a.levenshtein()))

# WLD(words, weight_insertion_type=WLD.WeightType.CUSTOM, custom_weight_insertion=2)
