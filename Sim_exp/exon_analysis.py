# text file contains the start and stop positions of the exons in two columns
# find the sensitivity and specificity of the exon prediction

# function to read the positions of the exons from the text file
def read_exon_positions(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    exon_positions = []
    for line in lines[1:]:  # Skip the header line
        start, end = map(int, line.strip().split('\t'))
        exon_positions.append((start, end))
    return exon_positions

# function to calculate the sensitivity and specificity of the exon prediction
def calculate_sensitivity_specificity(predicted_exons, true_exons):
    true_positive = 0
    false_positive = 0
    false_negative = 0

    # Calculate true positives and false negatives
    # the exons are a match if they overlap and the difference is less than 15% of the length of the exon
    for start, end in predicted_exons:
        for true_start, true_end in true_exons:
            # define overlap if abs(true_start - start) < 0.15 * (end - start) and true_end - end < 0.15 * (end - start)
            if (abs(true_start - start) < 0.15 * (end - start) and abs(true_end - end) < 0.15 * (end - start)):
                true_positive += 1
                break
        else:
            false_positive += 1
    for true_start, true_end in true_exons:
        for start, end in predicted_exons:
            if (abs(true_start - start) < 0.15 * (end - start) and abs(true_end - end) < 0.15 * (end - start)):
                break
        else:
            false_negative += 1
    # Calculate sensitivity and specificity

    sensitivity = true_positive / (true_positive + false_negative) if (true_positive + false_negative) > 0 else 0
    specificity = true_positive / (true_positive + false_positive) if (true_positive + false_positive) > 0 else 0

    return sensitivity, specificity

# run the code
if __name__ == "__main__":
    # Read the predicted exon positions from the text file
    predicted_exon_file = "exons.txt"
    predicted_exons = read_exon_positions(predicted_exon_file)

    # Read the true exon positions from the text file
    true_exon_file = "exon_positions.txt"
    true_exons = read_exon_positions(true_exon_file)

    # Calculate sensitivity and specificity
    sensitivity, specificity = calculate_sensitivity_specificity(predicted_exons, true_exons)

    print(f"Sensitivity: {sensitivity:.2f}")
    print(f"Specificity: {specificity:.2f}")
    