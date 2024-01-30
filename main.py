import queue
import random
import math

# ========================= Typehints ===============================
profile = list[dict[str, float]]
topHitsList = list[tuple[int, float]]
nodeList = dict[int, 'Node']

# ========================= Classes =============================
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
        # Used for the priority queue when determining the top hits order in case there are top hits with the same score
        return self.nodeId > other.nodeId

    def __repr__(self):
        """
        Used for easy string representation of the node should we desire to print the contents.
        Mainly meant for debugging purposes
        """
        return ("id: " + str(self.nodeId) +
                ", parent: " + str(self.parent) +
                ", children: " + ' '.join([str(n) for n in self.children] if self.children else ['none']) +
                ", tophits: " + ' '.join([str(n[0]) for n in self.topHits]))

    def findActiveAncestor(self, nodes: nodeList) -> int:
        """
        Finds the active ancestor of this node
        :return: The first active ancestor encountered or this node if this node is active
        """
        currentNode = self
        while not currentNode.active:
            currentNode = nodes[currentNode.parent]
        return currentNode.nodeId

    def initialize_top_hits(self, nodes: nodeList, active_nodes: list[int], m: int, totalProfile: profile) -> topHitsList:
        """
        Initializes the top hits list for this node
        :param nodes:   The list of all nodes
        :param active_nodes:    The list of current active nodes
        :param m:   m
        :param totalProfile: The average profile of all active nodes
        :return: The list of top hits for this node
        """
        # Keep a list of top hits for this sequence
        top_hits = []

        for nodeId in active_nodes:
            if nodeId == self.nodeId:
                continue
            seq = nodes[nodeId]
            # we are only calculating the distance now, not the criterion that should be minimized
            # criterion = d_u(i,j) - r(i) - r(j)
            score = nodeDistance(self, seq) - profileDistance(self.profile, totalProfile) - profileDistance(seq.profile, totalProfile)
            top_hits.append((nodeId, score))

        # Sort based on score
        top_hits.sort(key=lambda x: x[1])
        self.topHits = top_hits[:m]
        return [(self.nodeId, 0)] + top_hits
    
    def approximate_top_hits(self, seed_top_hits: topHitsList, m: int, nodes: nodeList, totalProfile: profile):
        """
        Approximates the top hits list for this node based on the top hits list from the seed
        :param seed_top_hits: The top hits from the seeds
        :param m: m
        :param nodes: The collection of all nodes in the tree
        :param totalProfile: The current total profile over all active nodes
        """
        top_hits = []
        for hit, _ in seed_top_hits[:JOIN_SAFETY_FACTOR * m]:
            if hit == self.nodeId:
                continue
            hitNode = nodes[hit]
            # Take the neighbour criterion for score: d_u(i, j) - r(i) - r(j),  r(x) = p(x, T)
            score = nodeDistance(self, hitNode) - profileDistance(self.profile, totalProfile) - profileDistance(hitNode.profile, totalProfile)
            top_hits.append((hit, score))

        top_hits.sort(key=lambda x: x[1])
        self.topHits = top_hits[:m]



# =========================== Globals =======================================
# Constants
ALPHABET = 'ACGT'
DATA_FILE = 'fasttree-input.aln'
JOIN_SAFETY_FACTOR = 2
TOP_HITS_CLOSENESS = 0.5
ROOT_NODE_ID = 0
VERBOSE = True

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

def createNewick(nodes, currentNode: int = ROOT_NODE_ID) -> str:
    """
    Recursively builds a newick tree from the given data
    :param nodes: the list of all nodes in the tree
    :param currentNode: the root node of the tree, default 0
    :return: the newick tree in string format
    """
    output = '('
    currentNode = nodes[currentNode]
    for child in currentNode.children:
        childNode = nodes[child]
        if not childNode.children:
            # subtract 1 because our tree uses 0 for the root node and the data uses 0 for the first genome
            output += str(child-1) + ','
        else:
            output += createNewick(nodes, childNode.nodeId) + ','
    output = output[:-1] + ')'
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

def computeTotalProfile(nodes: nodeList) -> profile:
    """
    Compute the total profile using the profiles from certain nodes
    :param nodes: The nodes to use to compute the total profile
    :return: The total profile
    """
    genomeLength = len(nodes[ROOT_NODE_ID].profile)
    alphabet = nodes[ROOT_NODE_ID].profile[0].keys()
    totalProfile = initializeProfile('', genomeLength, alphabet)
    activeNodes = nodes[ROOT_NODE_ID].children
    # Add all frequencies to a single profile
    for nodeId in activeNodes:
        child = nodes[nodeId]
        for i in range(genomeLength):
            for key in child.profile[i].keys():
                totalProfile[i][key] += child.profile[i][key]

    # Divide every frequency by the number of active nodes
    for i in range(genomeLength):
        for key in totalProfile[i]:
            totalProfile[i][key] = totalProfile[i][key] / len(activeNodes)

    return totalProfile

def updateTotalProfile(amountOfTerms: int, newProfile, totalProfile, oldProfile1, oldProfile2) -> profile:
    """
    Updates the total profile with the new profile in O(La) time
    :param amountOfTerms: The amount of profiles which have been used to compute the total so far
    :param newProfile: The new profile to update the total profile with
    :param totalProfile: The current total profile
    :param oldProfile1: Old profile to remove (optional)
    :param oldProfile2: Old profile to remove (optional)
    :return: The updated total profile
    """
    genomeLength = len(newProfile)
    for i in range(genomeLength):
        for key in totalProfile[i]:
            # Remove oldProfile1 from the total profile
            totalProfile[i][key] = ((totalProfile[i][key] * amountOfTerms) - oldProfile1[i][key]) / (amountOfTerms - 1)
            # Remove oldProfile2 from the total profile
            totalProfile[i][key] = ((totalProfile[i][key] * (amountOfTerms-1)) - oldProfile2[i][key]) / (amountOfTerms - 2)
            # Add the merged profile to the total profile
            totalProfile[i][key] = totalProfile[i][key] + (newProfile[i][key] - totalProfile[i][key]) / (amountOfTerms - 1)

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
    """
    Computes the distance between two profiles using a %-difference notion
    :param i: the first profile
    :param j: the second profile
    :return: the distance between profile i and profile j
    """
    total = 0
    genomeLength = len(i)
    for l in range(genomeLength):
        for key in i[l].keys():
            for otherKey in j[l].keys():
                # Only add the product of frequencies if the bases are different
                total += i[l][key] * j[l][otherKey] * int(key != otherKey)

    return total / genomeLength

def nodeDistance(node1: Node, node2: Node) -> float:
    """
    Calculates the distance between two nodes using the formula
    d_u(i, j) = P(i, j) - u_1 - u_2
    :param node1: The first node
    :param node2: The second node
    :return: The distance between the 2 nodes
    """
    return profileDistance(node1.profile, node2.profile) - node1.upDistance - node2.upDistance

def calculateUpDistance(node: Node, nodes: nodeList) -> float:
    """
    Calculates the updistance for a given node, with 0 for leaves and
    u(ij) = P(i, j) / 2 for inner nodes
    :param node: The node to calculate the updistance for
    :param nodes: The list of all nodes
    :return: The updistance for the node
    """
    if not node.children:
        return 0

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
    newNode = Node(len(nodes), ROOT_NODE_ID, mergeProfiles(node1.profile, node2.profile))

    # Remove previous nodes as children form their previous parents
    nodes[ROOT_NODE_ID].children.append(newNode.nodeId)
    nodes[node1.parent].children.remove(node1.nodeId)
    nodes[node2.parent].children.remove(node2.nodeId)

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

def initialize_top_hits(m: int, nodes: nodeList, activeNodes: list[int], totalProfile: profile):
    seedSequences = nodes[ROOT_NODE_ID].children.copy()

    while seedSequences != []:
        # Take an arbitrary seed sequence
        seed = random.choice(seedSequences)
        seedSequences.remove(seed)
        seedNode = nodes[seed]

        # Generate a top hits list for that sequence by comparing it with the neighbour joining criterion
        top_hits = seedNode.initialize_top_hits(nodes, activeNodes, m, totalProfile)

        # For each neighbour in the top m hits (the closest m neighbours)
        for neighbour, _ in top_hits[:m]:
            neighbourNode = nodes[neighbour]
            # If the top hits of the neighbour is still empty
            # Possible add an addition check which also must be true:
            #
            if (not neighbourNode.topHits and
                    profileDistance(seedNode.profile, neighbourNode.profile)/profileDistance(seedNode.profile, nodes[top_hits[2*m-1][0]].profile) < TOP_HITS_CLOSENESS):
                # The neighbour is a 'close neighbour', so we estimate the top hits for it
                neighbourNode.approximate_top_hits(top_hits, m, nodes, totalProfile)
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

def log_corrected_distance(profile1, profile2):
    # Implement the log-corrected distance calculation between two profiles
    # This is a simplified placeholder; you'll need to implement the actual calculation
    return math.log(1 + profileDistance(profile1, profile2))


def perform_nni(node1: Node, nodes: nodeList) -> None:
    """
    Performs the nearest neighbor interchange on a given node (F1) and its parent (F2).
    :param node: a child node (F2) to perform NNI on. will not work on root
    :param nodes: The list of all nodes.
    """
    parentId = node1.parent
    if parentId < 0:
        print("cannot do NNI on root node")
        return
    
    node2 = nodes[parentId]

    if len(node1.children) < 2 or len(node2.children) < 2:
        print("cannot preform NNI on leaf nodes")
        return
    
    # if node1.nodeId in node2.children:
    A, B = nodes[node1.children[0]], nodes[node1.children[1]]

    siblings = []
    for sibling in node2.children:
        if sibling != node1.nodeId:
            siblings.append(sibling)

    C = nodes[siblings[0]]
    if parentId == 0:
        D = nodes[siblings[1]]
    else:
        D = nodes[node2.parent]
    
    # Calculate distances for current and alternate topologies
    current_distance = log_corrected_distance(A.profile, B.profile) + log_corrected_distance(C.profile, D.profile)
    alt_distance_1 = log_corrected_distance(A.profile, C.profile) + log_corrected_distance(B.profile, D.profile)
    alt_distance_2 = log_corrected_distance(B.profile, C.profile) + log_corrected_distance(A.profile, D.profile)

    # Determine if an alternate topology has a lower distance
    if min(alt_distance_1, alt_distance_2) < current_distance:
        if alt_distance_1 < alt_distance_2:
            # Perform swap for the first alternate topology
            node1.children = [A.nodeId, C.nodeId]
            if parentId == 0:
                node2.children = [node1.nodeId, B.nodeId, D.nodeId]
            else:
                node2.children = [node1.nodeId, B.nodeId]
        else:
            # Perform swap for the second alternate topology
            node1.children = [B.nodeId, C.nodeId]
            if parentId == 0:
                node2.children = [node1.nodeId, A.nodeId, D.nodeId]
            else:
                node2.children = [node1.nodeId, A.nodeId]
        # Recompute profiles for affected nodes
        for n in range(len(node1.children)):
            nodes[node1.children[n]].parent = node1.nodeId
        for n in range(len(node2.children)):
            nodes[node2.children[n]].parent = node2.nodeId
        
        node1.profile = mergeProfiles(nodes[node1.children[0]].profile, nodes[node1.children[1]].profile)
        if parentId != 0:
            node2.profile = mergeProfiles(nodes[node2.children[0]].profile, nodes[node2.children[1]].profile)  


def perform_nni_rounds(nodes: nodeList, rounds: int) -> None:
    """
    Applies NNI to all applicable nodes in the tree.
    :param nodes: The list of all nodes.
    :param activeNodes: The list of active nodes.
    """
    for _ in range(rounds):
        for node_id, node in nodes.items():
            if len(node.children) == 2 and node_id != 0:  # Ensure it's an internal node with two children (not root)
                perform_nni(node, nodes)


# =============================== Algorithm =======================================


if __name__ == '__main__':
    data = readFile(DATA_FILE)
    genomeLength = len(data[0])
    amountOfGenomes = len(data)

    if VERBOSE:
        print('found {} genomes of length {}'.format(amountOfGenomes, genomeLength))
        print('initializing star topology...')

    # Create initial star topology
    nodes = {ROOT_NODE_ID: Node(ROOT_NODE_ID, -1, initializeProfile('', genomeLength, ALPHABET))}
    activesNodes = []
    for genome in data:
        newNode = Node(len(nodes), ROOT_NODE_ID, initializeProfile(genome, genomeLength, ALPHABET))
        nodes[ROOT_NODE_ID].children.append(newNode.nodeId)
        nodes[newNode.nodeId] = newNode
        activesNodes.append(newNode.nodeId)

    if VERBOSE:
        print('computing total profile...')
    # Create total profile
    totalProfile = computeTotalProfile(nodes)

    # create initial top Hits
    m = math.ceil(amountOfGenomes ** 0.5)
    if VERBOSE:
        print('initializing top hits...')
    initialize_top_hits(m, nodes, activesNodes, totalProfile)

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
        toppestHits = []
        for _ in range(m):
            toppestHits.append(q.get())

        # Compute the score for all best hits and take the best one
        bestHit = None
        bestValue = float('inf')
        for _, nodeFrom, nodeTo in toppestHits:
            node1 = nodes[nodeFrom]
            node2 = nodes[nodeTo]
            score = nodeDistance(node1, node2) - profileDistance(node1.profile, totalProfile) - (
                profileDistance(node2.profile, totalProfile))
            if score < bestValue:
                bestHit = (node1, node2)
                bestValue = score

        if VERBOSE:
            print('>> joining nodes {} and {}...'.format(bestHit[0].nodeId, bestHit[1].nodeId))

        # Merge the nodes
        mergedNode = mergeNodes(bestHit[0], bestHit[1], m, nodes)
        activesNodes.remove(bestHit[0].nodeId)
        activesNodes.remove(bestHit[1].nodeId)
        activesNodes.append(mergedNode.nodeId)
        nodes[mergedNode.nodeId] = mergedNode
        nodes[mergedNode.nodeId] = mergedNode
        #update the total profile
        totalProfile = updateTotalProfile((len(activesNodes)), mergedNode.profile, totalProfile, bestHit[0].profile, bestHit[1].profile)

    print(createNewick(nodes))
