# Neighbour joining with profiles
# Nearest neighbour interchange
# Top hits heuristic
import random

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
        self.bestHitDistance = float('inf')

    def findActiveAncestor(self) -> 'Node':
        currentNode = self
        while not currentNode.active:
            currentNode = currentNode.parent
        return currentNode


# =========================== Globals =======================================
ALPHABET = 'ACGT'
DATA_FILE = 'test-small.aln'
JOIN_ITERATIONS = 200
JOINS_DONE = 0
GENOME_LENGTH = 0
MAX_ID = 0
m = 0

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
def initializeProfile(genome: str, length: int, alphabet: str) -> profile:
    if genome == '':
        return [{base: 0 for base in alphabet} for _ in range(length)]

    return [{base: float(genome[i] == base) for base in alphabet} for i in range(len(genome))]

def computeTotalProfile() -> profile:
    nodesToCheck = [ROOT_NODE]
    newProfile = initializeProfile('', GENOME_LENGTH, ALPHABET)
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

    global TOTAL_PROFILE
    for i in range(GENOME_LENGTH):
        for letter in ALPHABET:
            TOTAL_PROFILE[i][letter] = ((MAX_ID - 1) * TOTAL_PROFILE[i][letter] + newNode.profile[i][letter]) / MAX_ID

    # Update the best hits
    for child in ROOT_NODE.children:
        if child.nodeId == newNode.nodeId:
            continue
        distance = nodeDistance(child, newNode)
        if distance < newNode.bestHitDistance:
            newNode.bestHit = child
        if distance < child.bestHitDistance:
            child.bestHit = newNode
            child.bestHitDistance = distance
        if not child.bestHit.active:
            child.bestHit = child.findActiveAncestor()
            child.bestHitDistance = distance

    # Recompute total profile after a certain number of iterations
    global JOINS_DONE
    JOINS_DONE += 1

    if JOINS_DONE >= JOIN_ITERATIONS:
        TOTAL_PROFILE = computeTotalProfile()

def outDistance(node: Node) -> float:
    if node.children:
        return len(ROOT_NODE.children) * profileDistance(node.profile, TOTAL_PROFILE) - profileDistance(node.children[0].profile, node.children[1].profile)

    return len(ROOT_NODE.children) * profileDistance(node.profile, TOTAL_PROFILE)

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

    # Create top hits list for each node
    while len(ROOT_NODE.children) > 2:
        bestMatch = (None, None)
        bestValue = float("inf")
        for child in ROOT_NODE.children:
            for otherChild in ROOT_NODE.children:
                if child.nodeId != otherChild.nodeId:
                    continue
                dist = nodeDistance(child, otherChild)
                if dist < bestValue:
                    bestMatch = (child, otherChild)
                    bestValue = dist

        mergeNodes(bestMatch[0], bestMatch[1])


    print(createNewick(ROOT_NODE))
