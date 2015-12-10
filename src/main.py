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
    with open("gene.train.processed.txt") as train_file:
        lines = train_file.readlines()
    for line in lines:
        word, tag = line.split()
        if (word, tag) in observation_counts:
            observation_counts[(word, tag)] += 1
        else:
            observation_counts[(word, tag)] = 1
    # print observation_counts

if __name__ == "__main__":
    transition_counts = ComputeTransitionProbabilities()
    # ComputeTransitionProbabilities()
    GetObservationCounts()
