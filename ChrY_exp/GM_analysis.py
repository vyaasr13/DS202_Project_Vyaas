# read a gtf file and a gff3 file and store all the exon coordinates
# store them as a tuple
# calculate the sensitivity and specificity of the gtf file against the gff3 file

import pandas as pd
import numpy as np
import gffutils as gff
import gtfparse as gtf
import tempfile
import os

# function to read the gtf file and store the exon coordinates
def extract_exons_from_gff(file_path):
    # Create a temporary database in memory
    db_path = tempfile.NamedTemporaryFile(delete=False).name
    db = gff.create_db(file_path, dbfn=db_path, disable_infer_transcripts=True, disable_infer_genes=True, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)

    exons = []
    for exon in db.features_of_type('exon'):
        start = exon.start
        end = exon.end
        exons.append((start, end))
    
    db.conn.close()  # Close the database connection
    os.remove(db_path)  # Clean up
    return exons

# function to calculate the sensitivity and specificity of the gtf file against the gff3 file
def calculate_sensitivity_specificity(gtf_exons, gff_exons):
    true_positive = 0
    false_positive = 0
    false_negative = 0

    matched_gff = set()

    for i, (start, end) in enumerate(gtf_exons):
        matched = False
        for j, (true_start, true_end) in enumerate(gff_exons):
            if j in matched_gff:
                continue
            if abs(true_start - start) < 0.15 * (end - start) and abs(true_end - end) < 0.15 * (end - start):
                true_positive += 1
                matched_gff.add(j)
                matched = True
                break
        if not matched:
            false_positive += 1

    for j, (true_start, true_end) in enumerate(gff_exons):
        if j not in matched_gff:
            false_negative += 1

    sensitivity = true_positive / (true_positive + false_negative) if (true_positive + false_negative) > 0 else 0
    specificity = true_positive / (true_positive + false_positive) if (true_positive + false_positive) > 0 else 0

    return sensitivity, specificity

# run the code
if __name__ == "__main__":
    # Read the exon positions from the gtf file
    gtf_file = "GM_annot.gtf"
    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file {gtf_file} not found.")
    gtf_exons = extract_exons_from_gff(gtf_file)

    # Read the exon positions from the gff3 file
    gff_file = "Y_annot.gff3"
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"GFF file {gff_file} not found.")
    gff_exons = extract_exons_from_gff(gff_file)

    # Calculate sensitivity and specificity
    sensitivity, specificity = calculate_sensitivity_specificity(gtf_exons, gff_exons)

    print(f"Sensitivity: {sensitivity:.2f}")
    print(f"Specificity: {specificity:.2f}")