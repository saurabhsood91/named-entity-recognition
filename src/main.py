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

    return transition_counts[x][y] += 1
    # print transition_counts


if __name__ == "__main__":
    transition_counts = ComputeTransitionProbabilities()