# Neighbour joining with profiles
# Nearest neighbour interchange
# Top hits heuristic
import copy
import random
import math

# Typehints
profile = list[dict[str, float]]

class Node:
    def __init__(self, nodeId: int, parent: 'Node' = None, profile: profile = None):
        self.nodeId = nodeId
        self.children = []
        self.parent = parent
        self.profile = profile
        self.active = True
        self.bestHit = None
        self.topHits = []
        self.bestHitDistance = float('inf')

    def findActiveAncestor(self) -> 'Node':
        currentNode = self
        while not currentNode.active:
            currentNode = currentNode.parent
        return currentNode

    def initialize_top_hits(self, active_nodes, m):
        top_hits = []
        for seq in active_nodes:
            if seq is not self:
                score = profileDistance(self.profile, seq.profile)
                top_hits.append((seq, score))

        # Sort based on score
        top_hits.sort(key=lambda x: x[1])
        self.topHits = top_hits[:m].copy()
        return [(self, 0)] + top_hits
    
    def approximate_top_hits(self, seed_top_hits, m):
        top_hits = []
        for hit, _ in seed_top_hits[:2*m]:
            if hit is not self:
                score = profileDistance(self.profile, hit.profile)
                top_hits.append((hit, score))
        top_hits.sort(key=lambda x: x[1])
        self.topHits = top_hits[:m]


# =========================== Globals =======================================
ALPHABET = 'ACGT'
DATA_FILE = 'test-small.aln'
GENOME_LENGTH = 0
MAX_ID = 0
m = 0
TOTAL_PROFILE = None
ROOT_NODE = None

# ======================= Util functions ====================================
def readFile(fileName: str) -> list[str]:
    """
    Reads the file with the input data and parses it
    :param fileName: The name of the file to where the input data resides
    :return: the data as a list of all sequences
    """
    with open(fileName, 'r') as file:
        rawData = file.readlines()

    return [line.strip() for line in rawData if line[0] != '>']

def createSeed(alphabet: str = ALPHABET, length: int = GENOME_LENGTH) -> str:
    return ''.join(random.choice(alphabet) for _ in range(length))

def createNewick(node: Node) -> str:
    output = '('
    for child in node.children:
        if not child.children:
            output += str(child.nodeId) + ','
        else:
            output += createNewick(child)
    output = output[:-1] + ');'
    return output

# ======================= Algorithm functions ================================
def initializeProfile(genome: str, length: int, alphabet: str = ALPHABET) -> profile:
    if genome == '':
        return [{base: 0 for base in alphabet} for _ in range(length)]

    return [{base: float(genome[i] == base) for base in alphabet} for i in range(len(genome))]

def computeTotalProfile(nodes: list[Node]) -> profile:
    genomeLength = len(nodes[0].profile)
    alphabet = nodes[0].profile[0].keys()
    totalProfile = initializeProfile('', genomeLength, alphabet)

    for child in nodes:
        for i in range(genomeLength):
            for key in child.profile[i].keys():
                totalProfile[i][key] += child.profile[i][key]

    for i in range(genomeLength):
        for key in totalProfile[i]:
            totalProfile[i][key] = totalProfile[i][key] / len(nodes)

    return totalProfile

def updateTotalProfile(amountOfTerms: int, newProfile: profile, totalProfile: profile = TOTAL_PROFILE) -> profile:
    genomeLength = len(newProfile)

    for i in range(genomeLength):
        for key in totalProfile[i]:
            totalProfile[i][key] = totalProfile[i][key] + (newProfile[i][key] - totalProfile[i][key]) / amountOfTerms
    return totalProfile

def mergeProfiles(profile1: profile, profile2: profile):
    genomeLength = len(profile1)
    alphabet = profile1[0].keys()
    return [{base: (profile1[i][base] + profile2[i][base]) / 2 for base in alphabet} for i in range(genomeLength)]

def profileDistance(i: profile, j: profile) -> float:
    total = 0
    genomeLength = len(i)
    # Go over every position
    for l in range(genomeLength):
        # Calculate the distance between all possible combinations
        for key in i[l].keys():
            for otherKey in j[l].keys():
                total += i[l][key] * j[l][otherKey] * int(key != otherKey)

    return total / genomeLength

def nodeDistance(node1: Node, node2: Node) -> float:
    return profileDistance(node1.profile, node2.profile) - upDistance(node1) - upDistance(node2)

def upDistance(node: Node) -> float:
    # Leafs are by definition 0
    if not node.children:
        return 0

    return profileDistance(node.children[0].profile, node.children[1].profile) / 2

# TODO, implement this correctly
def mergeNodes(node1: Node, node2: Node):
    global MAX_ID
    # Merge nodes and add them to root
    newNode = Node(MAX_ID, ROOT_NODE, mergeProfiles(node1.profile, node2.profile))
    ROOT_NODE.children.remove(node1)
    ROOT_NODE.children.remove(node2)
    ROOT_NODE.children.append(newNode)
    newNode.children = [node1, node2]
    node1.parent = newNode
    node2.parent = newNode
    node1.active = False
    node2.active = False
    MAX_ID += 1

def outDistance(node: Node) -> float:
    if node.children:
        return len(ROOT_NODE.children) * profileDistance(node.profile, TOTAL_PROFILE) - profileDistance(node.children[0].profile, node.children[1].profile)

    return len(ROOT_NODE.children) * profileDistance(node.profile, TOTAL_PROFILE)

def initialize_top_hits(m: int, root: Node):
    seedSequences = (root.children).copy()
    while seedSequences:
        seed = random.choice(seedSequences)
        seedSequences.remove(seed)
        top_hits = seed.initialize_top_hits(root.children, m)
        for seq, _ in top_hits[:m]:
            #if topHits is empty and neighbours are close enough.
            if not seq.topHits and profileDistance(seed.profile, seq.profile)/profileDistance(seed.profile, top_hits[2*m-1][0].profile) < 0.75:
                #below code never hits. Should be a valid criteria. Might be because of test size. 
                seq.approximate_top_hits(top_hits, m)
                seedSequences.remove(seq)

# =============================== Algorithm =======================================


if __name__ == '__main__':
    data = readFile(DATA_FILE)
    GENOME_LENGTH = len(data[0])

    # Create initial star topology
    ROOT_NODE = Node(-1, None, initializeProfile('', GENOME_LENGTH, ALPHABET))
    for genome in data:
        ROOT_NODE.children.append(Node(MAX_ID, ROOT_NODE, initializeProfile(genome, GENOME_LENGTH, ALPHABET)))
        MAX_ID += 1

    # Create total profile
    TOTAL_PROFILE = computeTotalProfile()

    #create initial top Hits
    m = math.ceil(len(data) ** 0.5) 
    initialize_top_hits(m, ROOT_NODE)

    for seq in ROOT_NODE.children:
        print(seq.topHits)
