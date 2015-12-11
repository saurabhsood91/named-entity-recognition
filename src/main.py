#! /usr/bin/python
import pprint
pp = pprint.PrettyPrinter(indent=4)

def ComputeTransitionProbabilities():
    transition_counts = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    with open("gene.train-1.txt") as train_file:
        lines = train_file.readlines()
    prev = 0
    curr = 1
    while curr < len(lines):
        line1 = lines[prev]
        line2 = lines[curr]
        if line2 == "\n":
            curr += 2
            prev += 2
            continue
        prev_tag = line1.split()[1]
        curr_tag = line2.split()[1]
        # print prev_tag, " ", curr_tag
        prev += 1
        curr += 1

        x = 0
        y = 0
        if prev_tag == "O":
            x = 1
        elif prev_tag == "B":
            x = 2
        if curr_tag == "O":
            y = 1
        elif curr_tag == "B":
            y = 2
        transition_counts[x][y] += 1
    i_counts = sum(transition_counts[0])
    o_counts = sum(transition_counts[1])
    b_counts = sum(transition_counts[2])

    transition_counts[0] = [float(x) / i_counts for x in transition_counts[0]]
    transition_counts[1] = [float(x) / o_counts for x in transition_counts[1]]
    transition_counts[2] = [float(x) / b_counts for x in transition_counts[2]]

    # print transition_counts
    return transition_counts
    # print transition_counts

def GetObservationCounts():
    observation_counts = {}
    tag_counts = {"I": 0, "O": 0, "B": 0}
    # word_counts = {}
    with open("gene.train.processed.txt") as train_file:
        lines = train_file.readlines()
    for line in lines:
        word, tag = line.split()
        if (word, tag) in observation_counts:
            observation_counts[(word, tag)] += 1
        else:
            observation_counts[(word, tag)] = 1
        # if word in word_counts:
        #     word_counts[word] += 1
        # else:
        #     word_counts[word] = 1
        tag_counts[tag] += 1
    # print tag_counts

    observation_counts_unk = {('UNK', 'I'): 0, ('UNK', 'O'): 0, ('UNK', 'B'): 0}
    for (word, tag) in observation_counts:
        count = observation_counts[(word, tag)]
        if count <= 1:
            observation_counts_unk[('UNK', tag)] += count
        else:
            observation_counts_unk[(word, tag)] = count

    # print observation_counts_unk[('UNK', 'I')]
    # print observation_counts_unk[('UNK', 'O')]
    # print observation_counts_unk[('UNK', 'B')]

    for (word, tag) in observation_counts_unk:
        observation_counts_unk[(word, tag)] = float(observation_counts_unk[(word, tag)]) / tag_counts[tag]
    # print observation_counts_unk[('UNK', 'I')]
    # print observation_counts_unk[('UNK', 'O')]
    # print observation_counts_unk[('UNK', 'B')]
    return observation_counts_unk

def viterbi(observation_list, transition_probabilities, observation_likelihoods):
    observation = " ".join(observation_list)
    # print observation
    obs_words = observation.split()
    obs_length = len(obs_words)
    viterbi_table = [[0, 0, 0] for i in xrange(0, obs_length)]
    backpointers = [[0, 0, 0] for i in xrange(0, obs_length)]

    # Initialization
    if (obs_words[0], "I") in observation_likelihoods:
        viterbi_table[0][0] = observation_likelihoods[(obs_words[0], "I")]
    else:
        viterbi_table[0][0] = observation_likelihoods[("UNK", "I")]
    if (obs_words[0], "O") in observation_likelihoods:
        viterbi_table[0][1] = observation_likelihoods[(obs_words[0], "O")]
    else:
        viterbi_table[0][1] = observation_likelihoods[("UNK", "O")]
    if (obs_words[0], "B") in observation_likelihoods:
        viterbi_table[0][2] = observation_likelihoods[(obs_words[0], "B")]
    else:
        viterbi_table[0][2] = observation_likelihoods[("UNK", "B")]

    # Continuation step
    for t in xrange(1, len(obs_words)):
        for i in [0, 1, 2]:
            if i == 0:
                s = "I"
            elif i == 1:
                s = "O"
            else:
                s = "B"
            if (obs_words[t], s) in observation_likelihoods:
                prob = observation_likelihoods[(obs_words[t], s)]
            else:
                prob = observation_likelihoods[('UNK', s)]
            # for j in xrange(0, 3):
            #     print viterbi_table[t - 1][j], " ", transition_probabilities[j][i]
            # print ""
            # print viterbi_table[t - 1][0], transition_probabilities[0][i]
            base_0 = viterbi_table[t - 1][0] * transition_probabilities[0][i]
            base_1 = viterbi_table[t - 1][1] * transition_probabilities[1][i]
            base_2 = viterbi_table[t - 1][2] * transition_probabilities[2][i]
            incoming_0 = base_0 * prob
            incoming_1 = base_1 * prob
            incoming_2 = base_2 * prob

            viterbi_table[t][i] = max(incoming_0, incoming_1, incoming_2)

            # compute backpointers
            if base_0 > base_1 and base_0 > base_2:
                backpointers[t][i] = 0
            if base_1 > base_0 and base_1 > base_2:
                backpointers[t][i] = 1
            if base_2 > base_0 and base_2 > base_1:
                backpointers[t][i] = 2

    # pp.pprint(zip(*backpointers))
    # pp.pprint(zip(*viterbi_table))
    tags = []
    arr = viterbi_table[-1]
    if arr[0] >= arr[1] and arr[0] >= arr[2]:
        tags.append("I")
        state = 0
    elif arr[1] > arr[0] and arr[1] > arr[2]:
        tags.append("O")
        state = 1
    elif arr[2] > arr[0] and arr[2] > arr[1]:
        tags.append("B")
        state = 2

    for i in xrange(len(backpointers) - 1, 0, -1):
        if backpointers[i][state] == 0:
            state = 0
            tags.append("I")
        elif backpointers[i][state] == 1:
            state = 1
            tags.append("O")
        elif backpointers[i][state] == 2:
            state = 2
            tags.append("B")
    tags.reverse()
    # print observation
    # print tags
    # print len(observation_list)
    # print len(tags)
    for i in xrange(0, len(observation_list)):
        print observation_list[i], "\t", tags[i]
    print
    # print "\n"
    # exit(1)

def ReadTestFile():
    with open("test.txt") as test_file:
        observations = []
        words = []
        for line in test_file:
            if line == "\n":
                observations.append(words)
                words = []
            else:
                word = line.split()[0]
                words.append(word)
                # print line.split()
        # print observations[0]
        return observations

if __name__ == "__main__":
    transition_probabilities = ComputeTransitionProbabilities()
    # ComputeTransitionProbabilities()
    observation_likelihoods = GetObservationCounts()
    observations = ReadTestFile()
    for element in observations:
        viterbi(element, transition_probabilities, observation_likelihoods)
    # viterbi("Comparison with alkaline phosphatases", transition_probabilities, observation_likelihoods)
