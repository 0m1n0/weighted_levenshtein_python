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
        # print(str1, str2)

        # add empty row and column
        n_rows = len(str1) + 1
        n_cols = len(str2) + 1

        # generate zero matrix
        # matrix = np.zeros((size_x, size_y))

        # dist = [[0 for x in range(n_cols)] for x in range(n_rows)]
        # print(dist)

        return "mono"

    def levenshtein(self):
        all_list_distance = []
        for str1, str2 in tqdm(self.pairwise()):
            list_distance = [str1, str2, self.mono_levenshtein(str1, str2)]
            all_list_distance.append(list_distance)
        print(all_list_distance)

a = LD(['fs', 'sdf', 'sdf', '1'])

print(a.pairwise())
a.levenshtein()

class WLD:
    pass


# wld = WL(Matrix, str)
# wld.process(aa)
#
# WL(Matrix, str, aa)
