import itertools
from typing import List
from tqdm import tqdm
import numpy as np


class LD:
    """Levenshtein distance
    """

    def __init__(self, words: List[str]):
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

        # first row and column, start from zero and indexed increasingly
        for i in range(n_row):
            matrix[i, 0] = i
        for j in range(n_col):
            matrix[0, j] = j

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

                # if two letters are identical, the substitution costs 0 otherwise 1
                w_sub = 0 if str1[i - 1] == str2[j - 1] else 1
                # weight of insertion and deletion
                w_ins, w_del = 1, 1

                matrix[i, j] = min(x00 + w_sub, x01 + w_ins, x10 + w_del)

        print(matrix)
        distance = matrix[n_row - 1, n_col - 1]
        return distance

    def levenshtein(self):
        all_list_distance = []
        for str1, str2 in tqdm(self.pairwise()):
            list_distance = [str1, str2, self.mono_levenshtein(str1, str2)]
            all_list_distance.append(list_distance)
        print(all_list_distance)


a = LD(['fs', 'sdf', 'sdf', '1'])
a.levenshtein()


class WLD:
    pass


# wld = WL(Matrix, str)
# wld.process(aa)
#
# WL(Matrix, str, aa)
