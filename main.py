# Neighbour joining with profiles
# Nearest neighbour interchange
# Top hits heuristic

# Typehints
profile = list[dict[str, float]]

# Globals
ALPHABET = 'ACGT'
DATA_FILE = 'test-small.alm'
L = -1                   # Will be set dynamically on loading the data

class Node:
    def __init__(self, nodeId: int):
        self.nodeId = nodeId
        self.children = []
        self.parent = None

def readFile(fileName: str) -> list[str]:
    with open(fileName, 'r') as file:
        rawData = file.readlines()

    return [line for line in rawData if line[0] != '>']

def init():
    data = readFile('')

def initializeProfile(length: int) -> profile:
    '''
    Create a profile matrix of specified length initialized with all zeros.
    This function only covers the 4 base pairs ACGT
    '''
    return [{base: 0 for base in ALPHABET} for _ in range(length)]

def profileDistance(profile1: profile, profile2: profile):
    '''
    The profile distance of the new Node formed by joined profile 1 and profile 2 (A and B)
    '''

    pass


def internalNodeDistance(node1, node2):
    '''
    The distance between two internal nodes
    '''

    pass

def upDistance(node):
    '''
    Average distance of the node from it's children. 0 for leaf nodes
    '''
    pass

def neighBourJoiningCriterion(i, j):
    pass

def averageOutDistance(i):
    pass
