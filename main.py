# Neighbour joining with profiles
# Nearest neighbour interchange
# Top hits heuristic
import queue
import random
import math

# Typehints
profile = list[dict[str, float]]
topHitsList = list[tuple[int, float]]
nodeList = dict[int, 'Node']

class Node:
    def __init__(self, nodeId: int, parent: int, profile: profile = None, children=None):
        if children is None:
            children = []
        self.nodeId = nodeId
        self.children = children
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
                " parent: " + str(self.parent) +
                "\nchildren: " + ' '.join([str(n) for n in self.children]) +
                "\ntophits: " + ' '.join([str(n[0]) for n in self.topHits]) + '}')

    def findActiveAncestor(self, nodes: nodeList) -> int:
        """
        Finds the active ancestor of this node
        :return: The first active ancestor encountered or this node if this node is active
        """
        currentNode = self
        while not currentNode.active:
            currentNode = nodes[currentNode.parent]
        return currentNode.nodeId

    def initialize_top_hits(self, nodes: nodeList, active_nodes: list[int], m: int) -> topHitsList:
        """
        Initializes the top hits list for this node
        :param nodes:   The list of all nodes
        :param active_nodes:    The list of current active nodes
        :param m:   m
        :return: The list of top hits for this node
        """
        # Keep a list of top hits for this sequence
        top_hits = []

        for nodeId in active_nodes:
            if nodeId == self.nodeId:
                continue
            seq = nodes[nodeId]
            score = profileDistance(self.profile, seq.profile)
            top_hits.append((nodeId, score))

        # Sort based on score
        top_hits.sort(key=lambda x: x[1])
        self.topHits = top_hits[:m]
        return [(self.nodeId, 0)] + top_hits
    
    def approximate_top_hits(self, seed_top_hits: topHitsList, m: int, nodes: nodeList) -> None:
        top_hits = []
        for hit, _ in seed_top_hits[:JOIN_SAFETY_FACTOR * m]:
            if hit == self.nodeId:
                continue
            hitNode = nodes[hit]
            score = profileDistance(self.profile, hitNode.profile)
            top_hits.append((hit, score))

        top_hits.sort(key=lambda x: x[1])
        self.topHits = top_hits[:m]


# =========================== Globals =======================================
# Constants
ALPHABET = 'ACGT'
DATA_FILE = 'test-small.aln'
JOIN_SAFETY_FACTOR = 2
ROOT_NODE_ID = 0

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

def computeTotalProfile(nodes: dict[int: Node]) -> profile:
    """
    Compute the total profile using the profiles from certain nodes
    :param nodes: The nodes to use to compute the total profile
    :return: The total profile
    """
    genomeLength = len(nodes[ROOT_NODE_ID].profile)
    alphabet = nodes[ROOT_NODE_ID].profile[0].keys()
    totalProfile = initializeProfile('', genomeLength, alphabet)
    activeNodes = nodes[ROOT_NODE_ID].children

    for nodeId in activeNodes:
        child = nodes[nodeId]
        for i in range(genomeLength):
            for key in child.profile[i].keys():
                totalProfile[i][key] += child.profile[i][key]

    for i in range(genomeLength):
        for key in totalProfile[i]:
            totalProfile[i][key] = totalProfile[i][key] / len(activeNodes)

    return totalProfile

def updateTotalProfile(amountOfTerms: int, newProfile: profile, totalProfile: profile) -> profile:
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
    """
    Calculates the distance between two nodes
    :param node1: The first node
    :param node2: The second node
    :return: The distance
    """
    # d_u(i, j) = P(i, j) - u_1 - u_2
    return profileDistance(node1.profile, node2.profile) - node1.upDistance - node2.upDistance

def calculateUpDistance(node: Node, nodes: nodeList) -> float:
    """
    Calculates the updistance for a given node
    :param node: The node to calculate the updistance for
    :param nodes: The list of all nodes
    :return: The updistance for the node
    """
    # Leafs are by definition 0
    if not node.children:
        return 0

    # u(ij) = P(i, j) / 2
    node1 = nodes[node.children[0]]
    node2 = nodes[node.children[1]]
    return profileDistance(node1.profile, node2.profile) / 2

def mergeNodes(node1: Node, node2: Node, m: int, nodes: nodeList) -> Node:
    """
    Takes two nodes and combines them according to the fast tree algorithm
    :param node1: First node to merge
    :param node2: Second node to merge
    :return: The newly made node having both param nodes as children
    """
    # Create new node
    newNode = Node(len(nodes), 0, mergeProfiles(node1.profile, node2.profile))

    # Add old nodes as children
    newNode.children = [node1.nodeId, node2.nodeId]
    node1.parent = newNode.nodeId
    node2.parent = newNode.nodeId

    # Calculate the updistance
    newNode.upDistance = calculateUpDistance(newNode, nodes)

    #Update the age of the newNode
    newNode.age = 1 + max(node1.age, node2.age)

    # Remove children nodes from top Hits list
    #node1.topHits.remove(node2.nodeId)
    #node2.topHits.remove(node1.nodeId)

    #Combine top Hits list and remove duplicates
    combinedList = list(dict.fromkeys(node1.topHits + node2.topHits))

    combinedTophits = []
    for child in combinedList:
        node = nodes[child[0]]
        score = profileDistance(node.profile, newNode.profile)
        combinedTophits.append((child[0], score))

    combinedTophits.sort(key=lambda x: x[1])

    #In case tophits node1 == node2
    if len(combinedTophits) == len(node1.topHits) and len(combinedTophits) == len(node2.topHits):
        newNode.topHits = combinedTophits[:m-1]
    else:
        newNode.topHits = combinedTophits[:m]

    # Set old nodes to inactive
    node1.active = False
    node2.active = False

    return newNode

def initialize_top_hits(m: int, nodes: nodeList, activeNodes: list[int]):
    seedSequences = nodes[ROOT_NODE_ID].children

    while seedSequences != []:
        # Take an arbitrary seed sequence
        seed = random.choice(seedSequences)
        seedSequences.remove(seed)
        seedNode = nodes[seed]

        # Generate a top hits list for that sequence by comparing it with the neighbour joining criterion
        top_hits = seedNode.initialize_top_hits(nodes, activeNodes, m)

        # For each neighbour in the top m hits (the closest m neighbours)
        for neighbour, _ in top_hits[:m]:
            neighbourNode = nodes[neighbour]
            # If the top hits of the neighbour is still empty
            # Possible add an addition check which also must be true:
            #
            if (not neighbourNode.topHits and
                    profileDistance(seedNode.profile, neighbourNode.profile)/profileDistance(seedNode.profile, nodes[top_hits[2*m-1][0]].profile) < 0.50):
                # The neighbour is a 'close neighbour', so we estimate the top hits for it
                neighbourNode.approximate_top_hits(top_hits, m, nodes)
                seedSequences.remove(neighbour)

def findBestJoin(topHits: topHitsList):
    bestCandidate = None
    bestCriterion = float("inf")
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
    genomeLength = len(data[0])
    amountOfGenomes = len(data)

    # Create initial star topology
    nodes = {0: Node(0, -1, initializeProfile('', genomeLength, ALPHABET))}
    activesNodes = []
    for genome in data:
        newNode = Node(len(nodes), 0, initializeProfile(genome, genomeLength, ALPHABET))
        nodes[0].children.append(newNode.nodeId)
        nodes[newNode.nodeId] = newNode
        activesNodes.append(newNode.nodeId)

    # Create total profile
    totalProfile = computeTotalProfile(nodes)

    # create initial top Hits
    m = math.ceil(amountOfGenomes ** 0.5)
    initialize_top_hits(m, nodes, activesNodes)

    # Do N - 3 joins
    while len(activesNodes) > 3:
        q = queue.PriorityQueue()
        for i in activesNodes:
            for topHit in nodes[i].topHits:
                nodeId1 = nodes[i].findActiveAncestor(nodes)
                nodeId2 = nodes[topHit[0]].findActiveAncestor(nodes)
                if nodeId1 != nodeId2:
                    q.put((topHit[1], nodeId1, nodeId2))    #  (Distance, node from, node to)

        # Get the best m hits over all nodes
        toppesthits = []
        for _ in range(m):
            toppesthits.append(q.get())

        # Compute the score for all best hits and take the best one
        bestHit = None
        bestValue = float('inf')
        for _, nodeFrom, nodeTo in toppesthits:
            node1 = nodes[nodeFrom]
            node2 = nodes[nodeTo]
            score = nodeDistance(node1, node2)
            if score < bestValue:
                bestHit = (node1, node2)
                bestValue = score

        # Merge the nodes
        mergedNode = mergeNodes(bestHit[0], bestHit[1], m, nodes)
        activesNodes.remove(bestHit[0].nodeId)
        activesNodes.remove(bestHit[1].nodeId)
        activesNodes.append(mergedNode.nodeId)
        nodes[mergedNode.nodeId] = mergedNode
