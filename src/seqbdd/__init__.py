"""The sequence binary decision diagram (SeqBDD) module \
implemented in C++."""

from typing import Callable, List, Tuple

import _seqbdd

class SeqBDD:
    """The sequence binary decision diagram (SeqBDD) class."""
    def __init__(self, sequences: List[str]) -> None:
        self._initialize(_seqbdd.from_sequences, sequences)


    def _initialize(self, f: Callable, *args) -> 'SeqBDD':
        self.root = f(*args)
        self.failure = None


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
        return _seqbdd.values(self.root)


    def __contains__(self, sequence: str) -> bool:
        """Whether a SeqBDD has a sequence.
        :param sequence:
            A sequence.
        :return:
            Whether a sequence is contained in a SeqBDD.
        """
        return _seqbdd.has_sequence(self.root, sequence)


    def __eq__(self, other: 'SeqBDD') -> bool:
        """Whether a SeqBDD is identical to an other.
        :param other:
            An other SeqBDD.
        :return:
            Wether two SeqBDDs are equal.
        """
        return _seqbdd.equal_nodes(self.root, other.root)


    def __ne__(self, other: 'SeqBDD') -> bool:
        """Whether a SeqBDD is not identical to an other.
        :param other:
            An other SeqBDD.
        :return:
            Wether two SeqBDDs are not equal.
        """
        return _seqbdd.not_equal_nodes(self.root, other.root)


    def __or__(self, right: 'SeqBDD') -> 'SeqBDD':
        """Get a union set.
        :param right:
            A SeqBDD.
        :return:
            A SeqBDD represents a union set.
        """
        sdd = self.__new__(type(self))
        sdd._initialize(_seqbdd.union_, self.root, right.root)

        return sdd


    def __and__(self, right: 'SeqBDD') -> 'SeqBDD':
        """Get an intersection set.
        :param right:
            A SeqBDD.
        :return:
            A SeqBDD represents an intersection set.
        """
        sdd = self.__new__(type(self))
        sdd._initialize(_seqbdd.intersection, self.root, right.root)

        return sdd


    def __sub__(self, right: 'SeqBDD') -> 'SeqBDD':
        """Get a difference set.
        :param right:
            A SeqBDD.
        :return:
            A SeqBDD represents a difference set.
        """
        sdd = self.__new__(type(self))
        sdd._initialize(_seqbdd.difference, self.root, right.root)

        return sdd


    def __xor__(self, right: 'SeqBDD') -> 'SeqBDD':
        """Get a symmetric difference set.
        :param right:
            A SeqBDD.
        :return:
            A SeqBDD represents a symmetric difference set.
        """
        sdd = self.__new__(type(self))
        sdd._initialize(_seqbdd.symmetric_difference, self.root, right.root)

        return sdd


    def onset(self, label: str) -> 'SeqBDD':
        """Get a onset.
        :param label:
            A character.
        :return:
            A SeqBDD.
        """
        sdd = self.__new__(type(self))
        sdd._initialize(_seqbdd.onset, self.root, label)

        return sdd


    def push(self, label: str) -> 'SeqBDD':
        """Get a SeqBDD is pushed a character.
        :param label:
            A character to push.
        :return:
            A SeqBDD.
        """
        sdd = self.__new__(type(self))
        sdd._initialize(_seqbdd.push, self.root, label)

        return sdd


    def top(self) -> str:
        """Get the character of root node.
        :return:
            A character.
        """
        return _seqbdd.top(self.root)


    def count(self) -> int:
        """Count sequences are contained in a SeqBDD.
        :return:
            The number of sequences.
        """
        return _seqbdd.count(self.root)


    def count_node(self) -> int:
        """Count nodes in a SeqBDD.
        :return:
            The number of nodes.
        """
        return _seqbdd.count_node(self.root)


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
        if self.failure is None:
            self.failure = _seqbdd.make_failure(self.root)

        return _seqbdd.search(self.root, self.failure, sequence)
