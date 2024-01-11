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
        self.assertEqual(4/6, main.profileDistance(profile1, profile2))

    def testComputeUpDistanceForLeafNode(self):
        innerNode = main.Node(nodeId=0, parent=None, profile=None)
        leafNode = main.Node(nodeId=1, parent=innerNode, profile=main.initializeProfile('ACGT', 4, 'ACGT'))
        leafNode2 = main.Node(nodeId=2, parent=innerNode, profile=main.initializeProfile('ATAT', 4, 'ACGT'))
        leafNode.children = []
        leafNode2.children = []
        innerNode.children = [leafNode, leafNode2]
        self.assertEqual(0, main.upDistance(leafNode))
        self.assertEqual(0.25, main.upDistance(innerNode))

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
        # Test for leaf nodes
        leaf1 = main.Node(nodeId=0, parent=None, profile=main.initializeProfile('ACGT', 4, 'ACGT'))
        leaf2 = main.Node(nodeId=1, parent=None, profile=main.initializeProfile('ACGA', 4, 'ACGT'))
        self.assertEqual(0.25, main.nodeDistance(leaf1, leaf2))

        # Test for leaf node and inner node
        innerNode1 = main.Node(nodeId=2, parent=None, profile=main.mergeProfiles(leaf1.profile, leaf2.profile))
        innerNode1.children = [leaf1, leaf2]
        leaf1.parent = innerNode1
        leaf2.parent = innerNode1
        leaf3 = main.Node(nodeId=3, parent=None, profile=main.initializeProfile('AGCG', 4, 'ACGT'))
        self.assertEqual(0.625, main.nodeDistance(innerNode1, leaf3))

        # Test for inner nodes
        leaf4 = main.Node(nodeId=4, parent=None, profile=main.initializeProfile('AGCC', 4, 'ACGT'))
        innerNode2 = main.Node(nodeId=5, parent=None, profile=main.mergeProfiles(leaf3.profile, leaf4.profile))
        innerNode2.children = [leaf3, leaf4]
        leaf3.parent = innerNode2
        leaf4.parent = innerNode2
        self.assertEqual(0.5, main.nodeDistance(innerNode1, innerNode2))

    def testSeedTopHitsLength(self):
        seed = main.Node(nodeId=0, parent=None, profile=main.initializeProfile('ACGT', 4, 'ACGT'))
        leaf1 = main.Node(nodeId=0, parent=None, profile=main.initializeProfile('ACGA', 4, 'ACGT'))
        leaf2 = main.Node(nodeId=1, parent=None, profile=main.initializeProfile('ACCA', 4, 'ACGT'))
        leaf3 = main.Node(nodeId=3, parent=None, profile=main.initializeProfile('AAAC', 4, 'ACGT'))
        m = 2
        nodes = [seed, leaf1, leaf2, leaf3]
        top_hits = seed.initialize_top_hits(nodes, m)
        self.assertEqual(3, len(top_hits))
        

    def testSeedTopHits(self):
        seed = main.Node(nodeId=0, parent=None, profile=main.initializeProfile('ACGT', 4, 'ACGT'))
        leaf1 = main.Node(nodeId=0, parent=None, profile=main.initializeProfile('ACGA', 4, 'ACGT'))
        leaf2 = main.Node(nodeId=1, parent=None, profile=main.initializeProfile('ACCA', 4, 'ACGT'))
        leaf3 = main.Node(nodeId=3, parent=None, profile=main.initializeProfile('AAAC', 4, 'ACGT'))
        m = 2
        nodes = [seed, leaf1, leaf2, leaf3]
        top_hits = seed.initialize_top_hits(nodes, m)
        self.assertEqual(leaf1, top_hits[0][0])
        self.assertEqual(leaf2, top_hits[1][0])
        self.assertEqual(leaf3, top_hits[2][0])

    def testTopHitsApproximation(self):
        seed = main.Node(nodeId=0, parent=None, profile=main.initializeProfile('ACGT', 4, 'ACGT'))
        leaf1 = main.Node(nodeId=0, parent=None, profile=main.initializeProfile('ACGA', 4, 'ACGT'))
        leaf2 = main.Node(nodeId=1, parent=None, profile=main.initializeProfile('AGGT', 4, 'ACGT'))
        leaf3 = main.Node(nodeId=3, parent=None, profile=main.initializeProfile('CGGT', 4, 'ACGT'))
        leaf4 = main.Node(nodeId=3, parent=None, profile=main.initializeProfile('AAAT', 4, 'ACGT'))
        leaf5 = main.Node(nodeId=3, parent=None, profile=main.initializeProfile('AAAG', 4, 'ACGT'))
        leaf6 = main.Node(nodeId=3, parent=None, profile=main.initializeProfile('AAAC', 4, 'ACGT'))
        leaf7 = main.Node(nodeId=3, parent=None, profile=main.initializeProfile('AACC', 4, 'ACGT'))
        leaf8 = main.Node(nodeId=3, parent=None, profile=main.initializeProfile('GAAA', 4, 'ACGT'))
        m = 3
        nodes = [seed, leaf1, leaf2, leaf3, leaf4, leaf5, leaf6, leaf7, leaf8]
        top_hits = seed.initialize_top_hits(nodes, m)
        leaf1.approximate_top_hits(top_hits, m)

        self.assertEqual(seed, leaf1.topHits[0][0])
        

if __name__ == '__main__':
    unittest.main()
