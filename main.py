# Neighbour joining with profiles
# Nearest neighbour interchange
# Top hits heuristic
import copy
import queue
import random
import math

# Typehints
profile = list[dict[str, float]]
topHitsList = list[tuple['Node', float]]

class Node:
    def __init__(self, nodeId: int, parent: 'Node' = None, profile: profile = None):
        self.nodeId = nodeId
        self.children = []
        self.parent = parent
        self.profile = profile
        self.upDistance = 0
        self.topHits = []
        self.active = True
        self.age = 0

    def __lt__(self, other):
        return self.nodeId > other.nodeId

    def __repr__(self):
        return ("{id: " + str(self.nodeId) +
                " parent: " + str(self.parent.nodeId) +
                "\nchildren: " + ' '.join([str(n.nodeId) for n in self.children]) +
                "\ntophits: " + ' '.join([str(n[0].nodeId) for n in self.topHits]) + '}')

    def findActiveAncestor(self) -> 'Node':
        currentNode = self
        while not currentNode.active:
            currentNode = currentNode.parent
        return currentNode

    def initialize_top_hits(self, active_nodes: list['Node'], m: int) -> topHitsList:
        # Keep a list of top hits for this sequence
        top_hits = []


        for seq in active_nodes:
            if seq is not self:
                score = profileDistance(self.profile, seq.profile)
                top_hits.append((seq, score))

        # Sort based on score
        top_hits.sort(key=lambda x: x[1])
        self.topHits = top_hits[:m].copy()
        return [(self, 0)] + top_hits
    
    def approximate_top_hits(self, seed_top_hits: list['Node'], m: int) -> None:
        top_hits = []
        for hit, _ in seed_top_hits[:2*m]:
            if hit is not self:
                score = profileDistance(self.profile, hit.profile)
                top_hits.append((hit, score))
        top_hits.sort(key=lambda x: x[1])
        self.topHits = top_hits[:m]


# =========================== Globals =======================================
# Constants
ALPHABET = 'ACGT'
DATA_FILE = 'test-small.aln'
ROOT_ID = -1
# Non-constants
GENOME_LENGTH = 0
MAX_ID = 0
TOTAL_PROFILE = None

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
    """
    Creates a new profile given a certain genome
    :param genome: The genome to convert to a profile, an empty genome creates a profile with all zeros
    :param length: The length of the profile, necessary when making an empty profile
    :param alphabet: The characters to use for the profile
    :return: The initialized profile
    """
    if genome == '':
        return [{base: 0 for base in alphabet} for _ in range(length)]

    return [{base: float(genome[i] == base) for base in alphabet} for i in range(len(genome))]

def computeTotalProfile(nodes: list[Node]) -> profile:
    """
    Compute the total profile using the profiles from certain nodes
    :param nodes: The nodes to use to compute the total profile
    :return: The total profile
    """
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
    """
    Updates the total profile with the new profile
    :param amountOfTerms: The amount of profiles which have been used to compute the total so far
    :param newProfile: The new profile to update the total profile with
    :param totalProfile: The current total profile
    :return: The updated total profile
    """
    genomeLength = len(newProfile)

    for i in range(genomeLength):
        for key in totalProfile[i]:
            totalProfile[i][key] = totalProfile[i][key] + (newProfile[i][key] - totalProfile[i][key]) / amountOfTerms
    return totalProfile

def mergeProfiles(profile1: profile, profile2: profile) -> profile:
    """
    Calculates the average of two profiles
    :param profile1: The first profile to merge
    :param profile2: The second profile to merge
    :return: The merged profile as the average of the two
    """
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
    """
    Calculates the updistance for a given node
    :param node: The node to calculate the updistance for
    :return: The updistance for the node
    """

    # Leafs are by definition 0
    if not node.children:
        return 0

    return profileDistance(node.children[0].profile, node.children[1].profile) / 2

def mergeNodes(node1: Node, node2: Node, m: int) -> Node:
    """
    Takes two nodes and combines them according to the fast tree algorithm
    :param node1: First node to merge
    :param node2: Second node to merge
    :return: The newly made node having both param nodes as children
    """
    # Create new node
    newNode = Node(0, None, mergeProfiles(node1.profile, node2.profile))

    # Add old nodes as children
    newNode.children = [node1, node2]
    node1.parent = newNode
    node2.parent = newNode

    # Calculate the updistance
    newNode.upDistance = upDistance(newNode)

    # Calculate the tophits list
    node1.topHits.remove(node2)
    node2.topHits.remove(node1)
    combinedList = node1.topHits
    combinedList = node1.topHits + node2.topHits

    newNode.initialize_top_hits()

    # Set old nodes to inactive
    node1.active = False
    node2.active = False

    return newNode

def outDistance(node: Node, totalProfile: profile = TOTAL_PROFILE) -> float:
    if not node.children:
        return len(root_node.children) * profileDistance(node.profile, totalProfile)

    return len(root_node.children) * profileDistance(node.profile, totalProfile) - profileDistance(node.children[0].profile, node.children[1].profile)


def initialize_top_hits(m: int, root: Node):
    seedSequences = copy.copy(root.children)
    while seedSequences:
        # Take an arbitrary seed sequence
        seed = random.choice(seedSequences)
        seedSequences.remove(seed)

        # Generate a top hits list for that sequence by comparing it with the neighbour joining criterion
        top_hits = seed.initialize_top_hits(root.children, m)

        # For each neighbour in the top m hits (the closest m neighbours)
        for neighbour, _ in top_hits[:m]:

            # If the top hits of the neighbour is still empty
            # Possible add an addition check which also must be true:
            # profileDistance(seed.profile, neighbour.profile)/profileDistance(seed.profile, top_hits[2*m-1][0].profile) < 0.75
            if not neighbour.topHits:
                neighbour.approximate_top_hits(top_hits, m)
                seedSequences.remove(neighbour)

def findBestJoin(topHits: topHitsList):
    bestCandidate = None
    bestCriterion = 1.0
    for hit in topHits:
        node = hit[1].findActiveAncestor()
        node.age += 1
        distance = hit[1]
        if distance < bestCriterion:
            bestCandidate = node
            bestCriterion = distance

    return bestCandidate



# =============================== Algorithm =======================================


if __name__ == '__main__':
    data = readFile(DATA_FILE)
    GENOME_LENGTH = len(data[0])
    amountOfGenomes = len(data)

    # Create initial star topology
    root_node = Node(ROOT_ID, None, initializeProfile('', GENOME_LENGTH, ALPHABET))
    for genome in data:
        root_node.children.append(Node(MAX_ID, root_node, initializeProfile(genome, GENOME_LENGTH, ALPHABET)))
        MAX_ID += 1

    # Create total profile
    TOTAL_PROFILE = computeTotalProfile(root_node.children)

    # create initial top Hits
    m = math.ceil(amountOfGenomes ** 0.5)
    initialize_top_hits(m, root_node)

    for seq in root_node.children:
        print(str(seq.nodeId) + ' ', end='')
        for topHit in seq.topHits:
            print('(' + str(topHit[0].nodeId) + " " + str(topHit[1]) + ')', end="")
        print()

    # Get the top m joins
    q = queue.PriorityQueue()
    for seq in root_node.children:
        for topHit in seq.topHits:
            q.put((topHit[1], seq, topHit[0]))  # (Distance, Seed node, neighbour node)

    bestJoins = [q.get() for _ in range(m)]

    for bestJoin in bestJoins:
        pass
        # print(bestJoin)

    while len(root_node.children) >= 3:
        bestJoin = findBestJoin(bestJoins)
        break
