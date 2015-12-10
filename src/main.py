#! /usr/bin/python

def ComputeTransitionProbabilities():
    transition_counts = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    with open("gene.train.processed.txt") as train_file:
        lines = train_file.readlines()
    prev = 0
    curr = 1
    while curr < len(lines):
        line1 = lines[prev]
        line2 = lines[curr]
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

    for (word, tag) in observation_counts:
        observation_counts[(word, tag)] = float(observation_counts[(word, tag)]) / tag_counts[tag]
    return observation_counts

def viterbi(observation, transition_probabilities, observation_likelihoods):
    obs_words = observation.split()
    obs_length = len(obs_words)
    viterbi_table = [[0, 0, 0] for i in xrange(0, obs_length)]
    backpointers = [[0, 0, 0] for i in xrange(0, obs_length)]

    # Initialization
    if (obs_words[0], "I") in observation_likelihoods:
        viterbi_table[0][0] = observation_likelihoods[(obs_words[0], "I")]
    else:
        viterbi_table[0][0] = 0
    if (obs_words[0], "O") in observation_likelihoods:
        viterbi_table[0][1] = observation_likelihoods[(obs_words[0], "O")]
    else:
        viterbi_table[0][1] = 0
    if (obs_words[0], "I") in observation_likelihoods:
        viterbi_table[0][2] = observation_likelihoods[(obs_words[0], "B")]
    else:
        viterbi_table[0][2] = 0

    # Continuation step
    for t in xrange(1, len(obs_words)):
        for i in [0, 1, 2]:
            if i == 0:
                s = "I"
            elif i == 1:
                s = "O"
            else:
                s == "B"
            if (obs_words[t], s) in observation_likelihoods:
                prob = observation_likelihoods[(obs_words[t], s)]
            else:
                prob = 0
            base_0 = viterbi_table[t - 1][0] * transition_probabilities[0][i]
            base_1 = viterbi_table[t - 1][1] * transition_probabilities[1][i]
            base_2 = viterbi_table[t - 1][2] * transition_probabilities[2][i]
            incoming_0 = base_0 * prob
            incoming_1 = base_1 * prob
            incoming_2 = base_2 * prob
            # viterbi_table[t][i] = max(incoming_0, incoming_1, incoming_2)

            if incoming_0 > incoming_1 and incoming_0 > incoming_2:
                viterbi_table[t][i] = incoming_0
            elif incoming_1 > incoming_0 and incoming_1 > incoming_2:
                viterbi_table[t][i] = incoming_1
            elif incoming_2 > incoming_0 and incoming_2 > incoming_1:
                viterbi_table[t][i] = incoming_2

            # compute backpointers
            if base_0 > base_1 and base_0 > base_2:
                backpointers[t][i] = 0
            if base_1 > base_0 and base_1 > base_2:
                backpointers[t][i] = 1
            if base_2 > base_0 and base_2 > base_1:
                backpointers[t][i] = 2

    print backpointers
    print viterbi_table

if __name__ == "__main__":
    transition_probabilities = ComputeTransitionProbabilities()
    # ComputeTransitionProbabilities()
    observation_likelihoods = GetObservationCounts()
    viterbi("Comparison with alkaline phosphatases and 5-nucleotidase", transition_probabilities, observation_likelihoods)
