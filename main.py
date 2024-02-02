import queue
import random
import math
import sys

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
        self.label = str(nodeId)
        self.distanceToParent = 1

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
            score = nodeDistance(self, seq) - calculateOutDistance(self, active_nodes, nodes, totalProfile) - calculateOutDistance(seq, active_nodes, nodes, totalProfile)
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
        activeNodes = nodes.keys()
        for hit, _ in seed_top_hits[:min(JOIN_SAFETY_FACTOR * m, len(seed_top_hits))]:
            if hit == self.nodeId:
                continue
            hitNode = nodes[hit]
            # Take the neighbour criterion for score: d_u(i, j) - r(i) - r(j),  r(x) = p(x, T)
            score = nodeDistance(self, hitNode) - calculateOutDistance(self, activeNodes,nodes, totalProfile) - calculateOutDistance(hitNode, activeNodes, nodes,totalProfile)
            top_hits.append((hit, score))

        top_hits.sort(key=lambda x: x[1])
        self.topHits = top_hits[:m]
    

    def merge_top_hits(self, seed_top_hits: topHitsList, m: int, nodes: nodeList, activesNodes: list[int], totalProfile: profile) -> None:
        """
        Approximates a new top hits list for this node based on the top hits list from a seed and an allready existing topHits
        :param seed_top_hits: The top hits from the seeds
        :param m: (minimum of sqrt(N) or active enodes left)
        :param nodes: The collection of all nodes in the tree
        :param totalProfile: The current total profile over all active nodes
        """
        new_top_hits = []
        #Get the score of entries seed top hits. 
        for hit, _ in seed_top_hits:
            if hit == self.nodeId:
                continue
            hitNode = nodes[hit]
            score = nodeDistance(self, hitNode) - calculateOutDistance(self, activesNodes, nodes, totalProfile) - calculateOutDistance(hitNode, activesNodes, nodes, totalProfile)
            new_top_hits.append((hit, score))

        new_top_hits.sort(key=lambda x: x[1])

        combinedList = list(dict.fromkeys(self.topHits[:m] + new_top_hits[:m]))

        #combine new top hit listwith existing one. 
        combinedTophits = []
        for child in combinedList:
            node = nodes[child[0]]
            score = nodeDistance(node, self) - calculateOutDistance(self, activesNodes, nodes, totalProfile) - calculateOutDistance(node, activesNodes, nodes, totalProfile)
            combinedTophits.append((child[0], score))

        combinedTophits.sort(key=lambda x: x[1])

        self.topHits = combinedTophits[:min(m, len(combinedTophits))]



# =========================== Globals =======================================
# Constants
ALPHABET = 'ACGT'
DATA_FILE = 'fasttree-input.aln'
JOIN_SAFETY_FACTOR = 2
TOP_HITS_CLOSENESS = 0.5
ROOT_NODE_ID = 0
VERBOSE = True
REFRESH_FACTOR = 0.8

# ======================= Util functions ====================================
def readFile(fileName: str) -> list[tuple[str, str]]:
    """
    Reads the file with the input data and parses it. If the data is invalid we terminate the script early.
    :param fileName: The name of the file to where the input data resides
    :return: the data as a list of a combination of sequences and their labels
    """
    try:
        with open(fileName, 'r') as file:
            rawData = file.readlines()
        data = []
        for i in range(0, len(rawData), 2):
            data.append((rawData[i].strip()[1:], rawData[i+1].strip()))

        return data
    except:
        print('please provide a valid data file')
        sys.exit(1)

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
            output += childNode.label + ':' + str(childNode.distanceToParent) + ','
        else:
            output += createNewick(nodes, childNode.nodeId) + ':' + str(childNode.distanceToParent) + ','
    output = output[:-1] + ')'
    return output

# ======================= Algorithm functions ================================

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


def calculateOutDistance(node: Node, activeNodes: list[int], nodes: nodeList, totalProfile) -> float:
    sum = 0
    n = len(activeNodes)
    for j in activeNodes:
        sum += nodes[j].upDistance

    outDistance = n*profileDistance(node.profile, totalProfile) - profileDistance(node.profile, node.profile) \
    - (n-1)*node.upDistance + node.upDistance - sum
    return outDistance/(n-2)




def mergeNodes(node1: Node, node2: Node, m: int, nodes: nodeList, totalProfile: profile) -> Node:
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

    #Remove children nodes from top Hits list
    node1_topHits = list(filter(lambda x: x[0] != node2.nodeId, node1.topHits))
    node2_topHits = list(filter(lambda x: x[0] != node1.nodeId, node2.topHits))

    #Combine top Hits list and remove duplicates
    combinedList = list(dict.fromkeys(node1_topHits + node2_topHits))
    nodes[newNode.nodeId] = newNode
    combinedTophits = []
    active_nodes = nodes[ROOT_NODE_ID].children
    for child in combinedList:
        node = nodes[child[0]]
        score = nodeDistance(node, newNode) -calculateOutDistance(node, active_nodes, nodes, totalProfile) - calculateOutDistance(newNode, active_nodes, nodes,  totalProfile)
        combinedTophits.append((child[0], score))

    combinedTophits.sort(key=lambda x: x[1])

    if len(combinedTophits) < m:
        newNode.topHits = combinedTophits
    else:
        newNode.topHits = combinedTophits[:m]

    # Set old nodes to inactive
    node1.active = False
    node2.active = False
    return newNode

def initialize_top_hits(m: int, nodes: nodeList, activeNodes: list[int], totalProfile: profile) -> None:
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

def refresh_top_hits(node: Node, m: int, nodes: nodeList, activeNodes: list[int], totalProfile: profile) -> None:
    """
    Creates a new top_hits list for a node, and refreshed the top hits of all m (or less if len(active nodes) < m) hits.
    :param node: the Node to refresh
    :param m: how many nodes to keep in the tophits
    :param nodes: The list of all nodes.
    :param activeNodes: The list of active nodes.
    :param totalProfile: The current total profile over all active nodes
    """
    top_hits = node.initialize_top_hits(nodes, activeNodes, m, totalProfile)
    for neighbour, _ in top_hits[:m]:
        neighbourNode = nodes[neighbour]
        neighbourNode.merge_top_hits(top_hits, m, nodes, activeNodes, totalProfile)
    node.age = 0

def log_corrected_distance(node1: Node, node2: Node) -> float:
    # Implement the log-corrected distance calculation between two profiles
    profileDist = profileDistance(node1.profile, node2.profile)
    val = -1.3 * math.log(1 - (profileDist if profileDist < 1 else 0.999))
    return val


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
    current_distance = log_corrected_distance(A, B) + log_corrected_distance(C, D)
    alt_distance_1 = log_corrected_distance(A, C) + log_corrected_distance(B, D)
    alt_distance_2 = log_corrected_distance(B, C) + log_corrected_distance(A, D)

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
    :param rounds: how many times NNI is to be applied.
    """
    for _ in range(rounds):
        for node_id, node in nodes.items():
            if len(node.children) == 2 and node_id != 0:  # Ensure it's an internal node with two children (not root)
                perform_nni(node, nodes)

def branchLength(parent: Node, child: Node, nodes):
    #In case the child node is a leaf. (Parent node will never be a leaf)
    if len(child.children) == 0:
        nodeA = child
        if parent.parent < 0:
            return 0
        nodeR = nodes[parent.parent]
        nodeB = nodes[[c for c in parent.children if c != child.nodeId][0]]
        branchLength = (log_corrected_distance(nodeA, nodeR) + log_corrected_distance(nodeA, nodeB) -
                        log_corrected_distance(nodeB, nodeR)) / 2
    else:
        nodeA = nodes[child.children[0]]
        nodeB = nodes[child.children[1]]
        nodeR = parent
        nodeC = nodes[[c for c in nodeR.children if c != child.nodeId][0]]
        branchLength = (log_corrected_distance(nodeA, nodeR) + log_corrected_distance(nodeA, nodeC) +
                        log_corrected_distance(nodeB, nodeR) + log_corrected_distance(nodeB, nodeC)) / 4 - (
                        log_corrected_distance(nodeA, nodeB) + log_corrected_distance(nodeR, nodeC)) / 2

    return branchLength if branchLength >= 0 else 0

# =============================== Algorithm =======================================


if __name__ == '__main__':
    if len(sys.argv) > 1:
        user_data = sys.argv[1]
        if user_data != '':
            DATA_FILE = user_data
    data = readFile(DATA_FILE)
    genomeLength = len(data[0][1])
    amountOfGenomes = len(data)

    if VERBOSE:
        print('found {} genomes of length {}'.format(amountOfGenomes, genomeLength))
        print('initializing star topology...')

    # Create initial star topology
    nodes = {ROOT_NODE_ID: Node(ROOT_NODE_ID, -1, initializeProfile('', genomeLength, ALPHABET))}
    activesNodes = []
    for label, genome in data:
        newNode = Node(len(nodes), ROOT_NODE_ID, initializeProfile(genome, genomeLength, ALPHABET))
        newNode.label = label
        nodes[ROOT_NODE_ID].children.append(newNode.nodeId)
        nodes[newNode.nodeId] = newNode
        activesNodes.append(newNode.nodeId)

    if VERBOSE:
        print('computing total profile...')
    # Create total profile
    totalProfile = computeTotalProfile(nodes)

    # create initial top Hits
    m = math.ceil(amountOfGenomes ** 0.5)
    max_age = int(1 + math.log2(m))
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
        for _ in range(min(m, q.qsize())):
            toppestHits.append(q.get())

        # Compute the score for all best hits and take the best one
        bestHit = None
        bestValue = float('inf')
        for _, nodeFrom, nodeTo in toppestHits:
            node1 = nodes[nodeFrom]
            node2 = nodes[nodeTo]
            score = nodeDistance(node1, node2) - calculateOutDistance(node1, activesNodes,nodes, totalProfile) - calculateOutDistance(node2, activesNodes, nodes,totalProfile)
            if score < bestValue:
                bestHit = (node1, node2)
                bestValue = score

        if VERBOSE:
            print('>> joining nodes {} and {}...'.format(bestHit[0].nodeId, bestHit[1].nodeId))

        # Merge the nodes
        node = mergeNodes(bestHit[0], bestHit[1], m, nodes, totalProfile)
        activesNodes.remove(bestHit[0].nodeId)
        activesNodes.remove(bestHit[1].nodeId)
        activesNodes.append(node.nodeId)
        nodes[node.nodeId] = node
        #update the total profile
        totalProfile = updateTotalProfile((len(activesNodes)), node.profile, totalProfile, bestHit[0].profile, bestHit[1].profile)
        if ((len(node.topHits) <= REFRESH_FACTOR*m or node.age >= max_age)) :
            if VERBOSE:
                print('Refreshing node {}...'.format(node.nodeId))
                refresh_top_hits(node, min(m, len(activesNodes)-1), nodes, activesNodes, totalProfile)

    if VERBOSE:
        print('performing NNI...')
    perform_nni_rounds(nodes, math.ceil(math.log2(len(nodes)) + 1))


    if VERBOSE:
        print('computing branch lengths')

    branchLengths = []
    for node in nodes:
        if nodes[node].parent <0:
            branchLengths.append((node, node, 0))
            nodes[node].distanceToParent = 0
        else:
            parent = nodes[nodes[node].parent]
            branchLengths.append((node, parent.nodeId, branchLength(parent, nodes[node], nodes)))
            nodes[node].distanceToParent = branchLength(parent, nodes[node], nodes)

    print('newick output:')
    print(createNewick(nodes))
