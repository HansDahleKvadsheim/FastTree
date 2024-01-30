import unittest
import main


class TestFasttreeMethods(unittest.TestCase):

    def testInitializeProfile(self):
        targetProfile = [
            {'A': 0, 'C': 0, 'G': 0, 'T': 1},
            {'A': 1, 'C': 0, 'G': 0, 'T': 0},
            {'A': 0, 'C': 0, 'G': 1, 'T': 0},
            {'A': 1, 'C': 0, 'G': 0, 'T': 0}
        ]
        self.assertCountEqual(targetProfile, main.initializeProfile('TAGA', 4, 'ACGT'))

        targetProfile = [
            {'A': 0, 'C': 0, 'G': 0, 'T': 0},
            {'A': 0, 'C': 0, 'G': 0, 'T': 0},
            {'A': 0, 'C': 0, 'G': 0, 'T': 0},
            {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        ]
        self.assertEqual(targetProfile, main.initializeProfile('', 4, 'ACGT'))

    def testComputeProfileDistance(self):
        profile1 = main.initializeProfile('ACGTTT', 6, 'ACGT')
        profile2 = main.initializeProfile('AGTCAT', 6, 'ACGT')
        self.assertEqual(4 / 6, main.profileDistance(profile1, profile2))

    def testComputeUpDistance(self):
        nodes = {
            0: main.Node(0, -1, profile=None),
            1: main.Node(1, 0, profile=main.initializeProfile('ACGT', 4, 'ACGT')),
            2: main.Node(2, 0, profile=main.initializeProfile('ATAT', 4, 'ACGT'))
        }
        nodes[0].children = [1, 2]

        self.assertEqual(0, main.calculateUpDistance(nodes[1], nodes))
        self.assertEqual(0.25, main.calculateUpDistance(nodes[0], nodes))

    def testMergeProfiles(self):
        profile1 = main.initializeProfile('AAGC', 4, 'ACGT')
        profile2 = main.initializeProfile('GACC', 4, 'ACGT')
        targetProfile = [
            {'A': 0.5, 'C': 0, 'G': 0.5, 'T': 0},
            {'A': 1, 'C': 0, 'G': 0, 'T': 0},
            {'A': 0, 'C': 0.5, 'G': 0.5, 'T': 0},
            {'A': 0, 'C': 1, 'G': 0, 'T': 0}
        ]
        self.assertEqual(targetProfile, main.mergeProfiles(profile1, profile2))

    def testNodeDistance(self):
        nodes = {
            # Leaf nodes
            1: main.Node(1, 5, profile=main.initializeProfile('ACGT', 4, 'ACGT')),
            2: main.Node(2, 5, profile=main.initializeProfile('ACGA', 4, 'ACGT')),
            3: main.Node(3, 6, profile=main.initializeProfile('AGCG', 4, 'ACGT')),
            4: main.Node(4, 6, profile=main.initializeProfile('AGCC', 4, 'ACGT'))
        }
        nodes[5] = main.Node(5, -1, profile=main.mergeProfiles(nodes[1].profile, nodes[2].profile), children=[1, 2])
        nodes[6] = main.Node(6, -1, profile=main.mergeProfiles(nodes[3].profile, nodes[4].profile), children=[3, 4])
        nodes[5].upDistance = main.calculateUpDistance(nodes[5], nodes)
        nodes[6].upDistance = main.calculateUpDistance(nodes[6], nodes)

        # Test for leaf nodes
        self.assertEqual(0.25, main.nodeDistance(nodes[1], nodes[2]))

        # Test for leaf node and inner node
        self.assertEqual(0.625, main.nodeDistance(nodes[5], nodes[3]))

        # Test for inner nodes
        self.assertEqual(0.5, main.nodeDistance(nodes[5], nodes[6]))

    def testSeedTopHitsLength(self):
        nodes = {
            1: main.Node(1, -1, profile=main.initializeProfile('ACGT', 4, 'ACGT')),
            2: main.Node(2, -1, profile=main.initializeProfile('ACGA', 4, 'ACGT')),
            3: main.Node(3, -1, profile=main.initializeProfile('ACCA', 4, 'ACGT')),
            4: main.Node(4, -1, profile=main.initializeProfile('AAAC', 4, 'ACGT'))
        }

        m = 2
        activeNodes = [1, 2, 3, 4]
        top_hits = nodes[1].initialize_top_hits(nodes, activeNodes, m)
        self.assertEqual(4, len(top_hits))

    def testSeedTopHits(self):
        nodes = {
            1: main.Node(1, -1, profile=main.initializeProfile('ACGT', 4, 'ACGT')),
            2: main.Node(2, -1, profile=main.initializeProfile('ACGA', 4, 'ACGT')),
            3: main.Node(3, -1, profile=main.initializeProfile('ACCA', 4, 'ACGT')),
            4: main.Node(4, -1, profile=main.initializeProfile('AAAC', 4, 'ACGT'))
        }

        m = 2
        active_nodes = [1, 2, 3, 4]
        top_hits = nodes[1].initialize_top_hits(nodes, active_nodes, m)
        self.assertEqual(1, top_hits[0][0])
        self.assertEqual(2, top_hits[1][0])
        self.assertEqual(3, top_hits[2][0])

        self.assertNotEqual(1, nodes[1].topHits[0][0])
        self.assertEqual(2, nodes[1].topHits[0][0])

    def testTopHitsApproximation(self):
        nodes = {
            1: main.Node(1, -1, profile=main.initializeProfile('ACGT', 4, 'ACGT')),
            2: main.Node(2, -1, profile=main.initializeProfile('ACGA', 4, 'ACGT')),
            3: main.Node(3, -1, profile=main.initializeProfile('AGGT', 4, 'ACGT')),
            4: main.Node(4, -1, profile=main.initializeProfile('CGGT', 4, 'ACGT')),
            5: main.Node(5, -1, profile=main.initializeProfile('AAAT', 4, 'ACGT'))
        }

        m = 2
        active_nodes = [1, 2, 3, 4, 5]
        top_hits = nodes[1].initialize_top_hits(nodes, active_nodes, m)
        nodes[2].approximate_top_hits(top_hits, m, nodes)

        self.assertEqual(1, nodes[2].topHits[0][0])
        self.assertEqual(3, nodes[2].topHits[1][0])

    def testInitializeTopHits(self):
        nodes = {
            0: main.Node(0, -1, main.initializeProfile('', 4, 'ACGT')),
            1: main.Node(1, 0, profile=main.initializeProfile('AAAA', 4, 'ACGT')),
            2: main.Node(2, 0, profile=main.initializeProfile('AATC', 4, 'ACGT')),
            3: main.Node(3, 0, profile=main.initializeProfile('AACG', 4, 'ACGT')),
            4: main.Node(4, 0, profile=main.initializeProfile('AAGT', 4, 'ACGT')),
            5: main.Node(5, 0, profile=main.initializeProfile('ACTA', 4, 'ACGT')),
        }
        nodes[0].children = [1, 2, 3, 4, 5]
        m = 2
        activeNodes = [1, 2, 3, 4, 5]

        main.initialize_top_hits(m, nodes, activeNodes)
        for leaf in activeNodes:
            self.assertEqual(2, len(nodes[leaf].topHits))

    def testComputeTotalProfile(self):
        innerProfile = [
            {'A': 0.5, 'C': 0, 'G': 0, 'T': 0.5},
            {'A': 1, 'C': 0, 'G': 0, 'T': 0},
            {'A': 0, 'C': 0.5, 'G': 0.5, 'T': 0},
            {'A': 0, 'C': 0, 'G': 1, 'T': 0}
        ]
        nodes = {
            0: main.Node(nodeId=0, parent=-1, profile=main.initializeProfile('', 4, 'ACGT')),
            1: main.Node(nodeId=1, parent=0, profile=innerProfile),
            2: main.Node(nodeId=2, parent=1, profile=main.initializeProfile('AACG', 4, 'ACGT')),
            3: main.Node(nodeId=3, parent=1, profile=main.initializeProfile('TAGG', 4, 'ACGT')),
            4: main.Node(nodeId=4, parent=0, profile=main.initializeProfile('TAAC', 4, 'ACGT')),
            5: main.Node(nodeId=5, parent=0, profile=main.initializeProfile('CACA', 4, 'ACGT'))
        }
        nodes[1].children = [2, 3]
        nodes[0].children = [1, 4, 5]

        expectedProfile = [
            {'A': 0.5 / 3, 'C': 1 / 3, 'G': 0, 'T': 0.5},
            {'A': 1, 'C': 0, 'G': 0, 'T': 0},
            {'A': 1 / 3, 'C': 0.5, 'G': 0.5 / 3, 'T': 0},
            {'A': 1 / 3, 'C': 1 / 3, 'G': 1 / 3, 'T': 0}
        ]
        self.assertEqual(expectedProfile, main.computeTotalProfile(nodes))

    def testUpdateTotalProfileWithoutMerge(self):
        totalProfile = [
            {'A': 0.5, 'C': 0, 'G': 0, 'T': 0.5},
            {'A': 1, 'C': 0, 'G': 0, 'T': 0},
            {'A': 0, 'C': 0.5, 'G': 0, 'T': 0.5},
            {'A': 0, 'C': 0, 'G': 0, 'T': 1}
        ]
        newProfile = main.initializeProfile('ACGT', 4, 'ACGT')
        expectedProfile = [
            {'A': 2 / 3, 'C': 0.0, 'G': 0.0, 'T': 1 / 3},
            {'A': 2 / 3, 'C': 1 / 3, 'G': 0.0, 'T': 0.0},
            {'A': 0.0, 'C': 1 / 3, 'G': 1 / 3, 'T': 1 / 3},
            {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 1.0}
        ]
        totalProfile = main.updateTotalProfile(3, newProfile, totalProfile)
        for i in range(4):
            for c in 'ACGT':
                self.assertAlmostEqual(expectedProfile[i][c], totalProfile[i][c])

    def testUpdateTotalProfileWithMerge(self):
        totalProfile = [
            {'A': 1/3, 'C': 0, 'G': 0, 'T': 2/3},
            {'A': 2/3, 'C': 1/3, 'G': 0, 'T': 0},
            {'A': 1/3, 'C': 1/3, 'G': 0, 'T': 1/3},
            {'A': 0, 'C': 1/3, 'G': 0, 'T': 2/3}
        ]
        oldProfile1 = main.initializeProfile('AACT', 4, 'ACGT')
        oldProfile2 = main.initializeProfile('TATT', 4, 'ACGT')
        mergedProfile = [
            {'A': 0.5, 'C': 0, 'G': 0, 'T': 0.5},
            {'A': 1, 'C': 0, 'G': 0, 'T': 0},
            {'A': 0, 'C': 0.5, 'G': 0, 'T': 0.5},
            {'A': 0, 'C': 0, 'G': 0, 'T': 1}
        ]
        expectedProfile = [
            {'A': 0.25, 'C': 0, 'G': 0, 'T': 0.75},
            {'A': 0.5, 'C': 0.5, 'G': 0, 'T': 0},
            {'A': 0.5, 'C': 0.25, 'G': 0, 'T': 0.25},
            {'A': 0, 'C': 0.5, 'G': 0, 'T': 0.5}
        ]
        totalProfile = main.updateTotalProfile(3, mergedProfile, totalProfile, oldProfile1, oldProfile2)
        for i in range(4):
            for c in 'ACGT':
                self.assertAlmostEqual(expectedProfile[i][c], totalProfile[i][c])


if __name__ == '__main__':
    unittest.main()
