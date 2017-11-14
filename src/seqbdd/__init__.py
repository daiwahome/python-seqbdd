"""The sequence binary decision diagram (SeqBDD) module \
implemented in C++."""

from typing import Callable, List, Tuple
import collections as _collections

import _seqbdd

AlignmentResult = _collections.namedtuple('AlignmentResult',
                                          ('sequence1',
                                           'sequence2',
                                           'aligned_sequence1',
                                           'aligned_sequence2',
                                           'span_sequence1',
                                           'span_sequence2',
                                           'score'))
AlignmentResult.__doc__ = 'The result of a pairwise sequence alignment.'
AlignmentResult.sequence1.__doc__ = 'A sequence.'
AlignmentResult.sequence2.__doc__ = 'A sequence.'
AlignmentResult.aligned_sequence1.__doc__ = 'A aligned sequence.'
AlignmentResult.aligned_sequence2.__doc__ = 'A aligned sequence.'
AlignmentResult.span_sequence1.__doc__ = \
    'The index pair of ``sequence1`` which is corresponding to a aligned sequence.'
AlignmentResult.span_sequence2.__doc__ = \
    'The index pair of ``sequence2`` which is corresponding to a aligned sequence.'
AlignmentResult.score.__doc__ = 'An alignment score.'

class SeqBDD:
    """The sequence binary decision diagram (SeqBDD) class."""
    def __init__(self, sequences: List[str]) -> None:
        self._initialize(_seqbdd.from_sequences, sequences)


    def _initialize(self, f: Callable, *args) -> 'SeqBDD':
        self._root: _seqbdd.Node = f(*args)
        self._failure: _seqbdd.Failure = None
        self._nodes: _seqbdd.NodeVector = None
        self._best_score: _seqbdd.BestScore = None


    @classmethod
    def empty(cls) -> 'SeqBDD':
        """Get the empty set.
        :return:
            The empty set.
        """
        sdd = cls.__new__(cls)
        sdd._initialize(_seqbdd.zero_terminal)

        return sdd


    @classmethod
    def null(cls) -> 'SeqBDD':
        """Get the a set only has empty sequence.
        :return:
            A set only has empty sequence.
        """
        sdd = cls.__new__(cls)
        sdd._initialize(_seqbdd.one_terminal)

        return sdd


    def values(self) -> List[str]:
        """Get sequences from a SeqBDD.
        :return:
            Sequences are obteined from a SeqBDD.
        """
        return _seqbdd.values(self._root)


    def __contains__(self, sequence: str) -> bool:
        """Whether a SeqBDD has a sequence.
        :param sequence:
            A sequence.
        :return:
            Whether a sequence is contained in a SeqBDD.
        """
        return _seqbdd.has_sequence(self._root, sequence)


    def __eq__(self, other: 'SeqBDD') -> bool:
        """Whether a SeqBDD is identical to an other.
        :param other:
            An other SeqBDD.
        :return:
            Wether two SeqBDDs are equal.
        """
        return _seqbdd.equal_nodes(self._root, other._root)


    def __ne__(self, other: 'SeqBDD') -> bool:
        """Whether a SeqBDD is not identical to an other.
        :param other:
            An other SeqBDD.
        :return:
            Wether two SeqBDDs are not equal.
        """
        return _seqbdd.not_equal_nodes(self._root, other._root)


    def __or__(self, right: 'SeqBDD') -> 'SeqBDD':
        """Get a union set.
        :param right:
            A SeqBDD.
        :return:
            A SeqBDD represents a union set.
        """
        sdd = self.__new__(type(self))
        sdd._initialize(_seqbdd.union_, self._root, right._root)

        return sdd


    def __and__(self, right: 'SeqBDD') -> 'SeqBDD':
        """Get an intersection set.
        :param right:
            A SeqBDD.
        :return:
            A SeqBDD represents an intersection set.
        """
        sdd = self.__new__(type(self))
        sdd._initialize(_seqbdd.intersection, self._root, right._root)

        return sdd


    def __sub__(self, right: 'SeqBDD') -> 'SeqBDD':
        """Get a difference set.
        :param right:
            A SeqBDD.
        :return:
            A SeqBDD represents a difference set.
        """
        sdd = self.__new__(type(self))
        sdd._initialize(_seqbdd.difference, self._root, right._root)

        return sdd


    def __xor__(self, right: 'SeqBDD') -> 'SeqBDD':
        """Get a symmetric difference set.
        :param right:
            A SeqBDD.
        :return:
            A SeqBDD represents a symmetric difference set.
        """
        sdd = self.__new__(type(self))
        sdd._initialize(_seqbdd.symmetric_difference, self._root, right._root)

        return sdd


    def onset(self, label: str) -> 'SeqBDD':
        """Get a onset.
        :param label:
            A character.
        :return:
            A SeqBDD.
        """
        sdd = self.__new__(type(self))
        sdd._initialize(_seqbdd.onset, self._root, label)

        return sdd


    def push(self, label: str) -> 'SeqBDD':
        """Get a SeqBDD is pushed a character.
        :param label:
            A character to push.
        :return:
            A SeqBDD.
        """
        sdd = self.__new__(type(self))
        sdd._initialize(_seqbdd.push, self._root, label)

        return sdd


    def top(self) -> str:
        """Get the character of root node.
        :return:
            A character.
        """
        return _seqbdd.top(self._root)


    def count(self) -> int:
        """Count sequences are contained in a SeqBDD.
        :return:
            The number of sequences.
        """
        return _seqbdd.count(self._root)


    def count_node(self) -> int:
        """Count nodes in a SeqBDD.
        :return:
            The number of nodes.
        """
        return _seqbdd.count_node(self._root)


    @classmethod
    def suffixdd(cls, sequence: str) -> 'SeqBDD':
        """Generate a SuffixDD from a sequence.
        :param sequence:
            A sequence to build all substrings index.
        :return:
            A SuffixDD.
        """
        sdd = cls.__new__(cls)
        sdd._initialize(_seqbdd.suffixdd, sequence)

        return sdd


    def search(self, sequence: str) -> List[Tuple[int, str]]:
        """Search a sequence by a SeqBDD.
        :param sequence:
            A sequence to search.
        :return:
            A SeqBDD.
        """
        if self._failure is None:
            self._failure = _seqbdd.make_failure(self._root)

        return _seqbdd.search(self._root, self._failure, sequence)


    def glocal_alignment(self, sequence: str, gapopen: int, gapext: int,
                         p: int) -> List['AlignmentResult']:
        """Align a sequence with a SeqBDD glocally.
        :param sequence:
            A sequence to align locally.
        :return:
            A SeqBDD to align globally.
        """
        if self._nodes is None:
            self._nodes = _seqbdd.get_nodes(self._root)
        if self._best_score is None:
            self._best_score = _seqbdd.get_best_score(self._root, gapopen, gapext)

        results = _seqbdd.glocal_alignment(sequence, self._root, self._nodes,
                                           self._best_score, gapopen, gapext, p)

        return [AlignmentResult(*result) for result in results]
