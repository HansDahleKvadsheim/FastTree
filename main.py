# Neighbour joining with profiles
# Nearest neighbour interchange
# Top hits heuristic
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
        self.topHits = []

    def findActiveAncestor(self) -> 'Node':
        currentNode = self
        while not currentNode.active:
            currentNode = currentNode.parent
        return currentNode
    
    def initialize_top_hits(self, active_nodes, m):
        # Assuming 'sequences' is a list of sequences and self.profile is the profile of the current node
        #active_nodes cannot contain seed
        top_hits = []
        for seq in active_nodes:
            score = profileDistance(self.profile, seq.profile)
            top_hits.append((seq, score))

        # Sort based on score
        top_hits.sort(key=lambda x: x[1], reverse=True)
        self.topHits = top_hits
        self.topHits = top_hits[:m]
    
    def approximate_top_hits_for_neighbors(self, neighbors, seed_top_hits, m):
        # Assuming 'neighbors' is a list of neighbor nodes (which we get from the top_hits seed)
        # and 'seed_top_hits' is the top 2m hits of the seed
        for neighbor in neighbors:
            neighbor_hits = []
            for hit, _ in seed_top_hits[:2*m]:
                score = profileDistance(neighbor.profile, hit.profile)
                neighbor_hits.append((hit, score))

            neighbor_hits.sort(key=lambda x: x[1], reverse=True)
            neighbor.topHits = neighbor_hits[:m]

    def refresh_top_hits(self, active_nodes, m):
        # Re-compute the top-hit list
        new_top_hits = []
        for node in active_nodes:
            if node is not self:
                score = profileDistance(self.profile, node.profile)
                new_top_hits.append((node, score))

        new_top_hits.sort(key=lambda x: x[1], reverse=True)
        self.topHits = new_top_hits[:m]

        # Update the close neighbors' top-hit lists
        for neighbor, _ in self.topHits[:m]:
            neighbor.approximate_top_hits_for_neighbors([self], self.topHits, m)
    


def update_top_hits_on_merge(sequence, other_sequence, m):
    #updates the top hits of seq A and B to the top m hits combined. 
    merged_hits = sorted(sequence.topHits + other_sequence.topHits, key=lambda x: x[1], reverse=True)

    combined_hits = []
    seen_seqs = set()
    for seq, score in merged_hits:
        if seq not in seen_seqs and seq not in {sequence, other_sequence}:
            combined_hits.append((seq, score))
            seen_seqs.add(seq)

        # Break if combined_hits reaches desired length
        if len(combined_hits) == m:
            break
    return combined_hits
    



# Globals
ALPHABET = 'ACGT'
DATA_FILE = 'test-small.aln'
JOIN_ITERATIONS = 200
JOINS_DONE = 0
GENOME_LENGTH = 0
MAX_ID = 0
m = 0
ACTIVE_NODES = 0

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


# ======================= Algorithm functions ================================
def initializeProfile(genome: str = '') -> profile:
    if genome == '':
        return [{base: 0 for base in ALPHABET} for _ in range(GENOME_LENGTH)]

    return [{base: float(genome[i] == base) for base in ALPHABET} for i in range(GENOME_LENGTH)]

def computeTotalProfile() -> profile:
    nodesToCheck = [ROOT_NODE]
    newProfile = initializeProfile()
    global ACTIVE_NODES
    ACTIVE_NODES = -1

    while nodesToCheck:
        currentNode = nodesToCheck.pop(0)
        for i in range(GENOME_LENGTH):
            for key in currentNode.profile[i]:
                newProfile[i][key] += currentNode.profile[i][key]

        for child in currentNode.children:
            if child.active:
                nodesToCheck.append(child)
        ACTIVE_NODES += 1

    for i in range(GENOME_LENGTH):
        for key in newProfile[i]:
            newProfile[i][key] = newProfile[i][key] / ACTIVE_NODES

    return newProfile

def mergeProfiles(profile1: profile, profile2: profile):
    return [{base: (profile1[i][base] + profile2[i][base]) / 2 for base in ALPHABET} for i in range(GENOME_LENGTH)]

def profileDistance(i: profile, j: profile) -> float:
    total = 0
    # Go over every position
    for l in range(GENOME_LENGTH):
        # Calculate the distance between all possible combinations
        for key in i[l].keys():
            for otherKey in j[l].keys():
                total += i[l][key] * j[l][otherKey] * int(key != otherKey)

    return total / GENOME_LENGTH

def nodeDistance(node1: Node, node2: Node) -> float:
    return profileDistance(node1.profile, node2.profile) - upDistance(node1) - upDistance(node2)

def upDistance(node: Node) -> float:
    # Leafs are by definition 0
    if not node.children:
        return 0

    return profileDistance(node.children[0].profile, node.children[1].profile) / 2

def mergeNodes(node1: Node, node2: Node):
    global MAX_ID
    MAX_ID += 1
    newNode = Node(MAX_ID, ROOT_NODE, mergeProfiles(node1.profile, node2.profile))
    ROOT_NODE.children.remove(node1)
    ROOT_NODE.children.remove(node2)
    ROOT_NODE.children.append(newNode)
    newNode.children = [node1, node2]
    node1.active = False
    node2.active = False
    global JOINS_DONE
    JOINS_DONE += 1

    if JOINS_DONE >= JOIN_ITERATIONS:
        global TOTAL_PROFILE
        TOTAL_PROFILE = computeTotalProfile()

def outDistance(node: Node) -> float:
    if node.children:
        return ACTIVE_NODES * profileDistance(node.profile, TOTAL_PROFILE) - profileDistance(node.children[0].profile, node.children[1].profile)

    return ACTIVE_NODES * profileDistance(node.profile, TOTAL_PROFILE)

def topHits(node: Node):
    pass


# =============================== Algorithm =======================================
data = readFile(DATA_FILE)
GENOME_LENGTH = len(data[0])


# Create initial topology
ROOT_NODE = Node(-1, None, initializeProfile(''))

for genome in data:
    ROOT_NODE.children.append(Node(MAX_ID, ROOT_NODE, initializeProfile(genome)))
    MAX_ID += 1

# Create total profile
TOTAL_PROFILE = computeTotalProfile()

# Find the top hits for each sequence
m = len(data) ** 0.5
seed = createSeed(ALPHABET, GENOME_LENGTH)
results = []

for child in ROOT_NODE.children:
    results.append(outDistance(child))

# TODO work it out