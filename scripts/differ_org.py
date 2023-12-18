import difflib

# Function to compare two files and generate differences
def compare_files(file1, file2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

    # Compute the difference between the two files
    differ = difflib.Differ()
    diff = list(differ.compare(lines1, lines2))

    # Extract differing lines from each file
    diff_file1 = [line[2:] for line in diff if line.startswith('- ')]
    diff_file2 = [line[2:] for line in diff if line.startswith('+ ')]

    return diff_file1, diff_file2

# Input file paths
file1_path = 'org_olp.txt'
file2_path = 'org_ilp.txt'

# Compare the files
differ_1, differ_2 = compare_files(file1_path, file2_path)

# Write the differing content to new files
with open('differ_org_olp.txt', 'w') as df1, open('differ_org_ilp.txt', 'w') as df2:
    df1.writelines(differ_1)
    df2.writelines(differ_2)

print("Differences saved in differ_1.txt and differ_2.txt")
