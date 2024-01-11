import unittest

import main


class TestFasttreeMethods(unittest.TestCase):

    def testCreateProfile(self):
        targetProfile = [{'A': 1, 'C': 0, 'G': 0, 'T': 0}, {'A': 1, 'C': 0, 'G': 0, 'T': 0}, {'A': 0, 'C': 0, 'G': 0, 'T': 1}]
        self.assertEqual(targetProfile, main.initializeProfile('AAT', 3, 'ACGT'))




if __name__ == '__main__':
    unittest.main()
