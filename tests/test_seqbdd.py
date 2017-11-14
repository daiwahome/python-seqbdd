"""The unit test for seqbdd module."""
from unittest import TestCase

class SeqBDDTestCase(TestCase):
    """The unit test for seqbdd.SeqBDD."""
    def setUp(self):
        from seqbdd import SeqBDD

        self.seq1 = {'aa', 'ab'}
        self.seq2 = {'aa', 'ac'}
        self.sdd1 = SeqBDD(list(self.seq1))
        self.sdd2 = SeqBDD(list(self.seq2))


    def test_instance(self):
        """Generate an instance of the SeqBDD"""
        from seqbdd import SeqBDD

        self.assertEqual(SeqBDD(['a']).values(), ['a'])


    def test_empty(self):
        """Generate the empty set"""
        from seqbdd import SeqBDD

        self.assertFalse(bool(SeqBDD.empty().values()))


    def test_null(self):
        """Generate a set only has empty sequence"""
        from seqbdd import SeqBDD

        self.assertEqual(SeqBDD.null().values(), [''])


    def test_contains(self):
        """Work `in` statements"""
        from seqbdd import SeqBDD

        sdd = SeqBDD(['aa', 'ab'])
        self.assertTrue('aa' in sdd)
        self.assertFalse('ac' in sdd)


    def test_or(self):
        """Get a union set"""
        union = self.sdd1 | self.sdd2
        self.assertEqual(union.values(), sorted(list(self.seq1|self.seq2)))


    def test_and(self):
        """Get an intersection set"""
        inter = self.sdd1 & self.sdd2
        self.assertEqual(inter.values(), sorted(list(self.seq1&self.seq2)))


    def test_sub(self):
        """Get an difference set"""
        diff = self.sdd1 - self.sdd2
        self.assertEqual(diff.values(), sorted(list(self.seq1-self.seq2)))


    def test_xor(self):
        """Get an symmetric difference set"""
        diff = self.sdd1 ^ self.sdd2
        self.assertEqual(diff.values(), sorted(list(self.seq1^self.seq2)))


    def test_onset(self):
        """Get an onset"""
        from seqbdd import SeqBDD

        self.assertEqual(SeqBDD(['bc', 'bd']).onset('b').values(), ['c', 'd'])
        self.assertEqual(SeqBDD(['bc', 'bd']).onset('a').values(), [])


    def test_push(self):
        """Get a SeqBDD is pushed a character"""
        from seqbdd import SeqBDD

        self.assertEqual(SeqBDD(['bc', 'bd']).push('a').values(),
                         ['abc', 'abd'])


    def test_top(self):
        """Get the character of root node"""
        from seqbdd import SeqBDD

        self.assertEqual(SeqBDD.empty().top(), '0')
        self.assertEqual(SeqBDD(['ab', 'cd']).top(), 'c')


    def test_count(self):
        """Count sequences are contained in a SeqBDD"""
        from seqbdd import SeqBDD

        self.assertEqual(SeqBDD.empty().count(), 0)
        self.assertEqual(SeqBDD.null().count(), 1)
        self.assertEqual(SeqBDD(['ab', 'cd']).count(), 2)


    def test_count_node(self):
        """Count nodes in a SeqBDD"""
        from seqbdd import SeqBDD

        self.assertEqual(SeqBDD.empty().count_node(), 0)
        self.assertEqual(SeqBDD.null().count_node(), 0)
        self.assertEqual(SeqBDD(['aa', 'ab']).count_node(), 3)


    def test_suffixdd(self):
        """Get a SuffixDD"""
        from seqbdd import SeqBDD

        self.assertEqual(SeqBDD.suffixdd('abc').values(),
                         ['', 'a', 'ab', 'abc', 'b', 'bc', 'c'])


    def test_search(self):
        """Search a sequence by a SeqBDD"""
        from seqbdd import SeqBDD

        sdd = SeqBDD(['abc', 'def'])
        self.assertEqual(sdd.search('abcdef'), [(0, 'abc'), (3, 'def')])
        self.assertEqual(sdd.search('abcdefabcdef'),
                         [(0, 'abc'), (3, 'def'), (6, 'abc'), (9, 'def')])

        sdd = SeqBDD(['ab', 'bc', 'bab', 'd', 'abcde'])
        self.assertEqual(sdd.search('xbabcdex'),
                         [(1, 'bab'), (2, 'ab'), (2, 'abcde'),
                          (3, 'bc'), (5, 'd')])
        self.assertEqual(sdd.search('abc'), [(0, 'ab'), (1, 'bc')])
        self.assertEqual(sdd.search('ab'), [(0, 'ab')])

        sdd = SeqBDD(['a', 'ab', 'bc', 'bca', 'c', 'caa'])
        self.assertEqual(sdd.search('abccab'),
                         [(0, 'a'), (0, 'ab'), (1, 'bc'), (2, 'c'),
                          (3, 'c'), (4, 'a'), (4, 'ab')])
