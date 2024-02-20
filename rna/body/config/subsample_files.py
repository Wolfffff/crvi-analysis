import pandas as pd
import random

# File paths
input_file1 = 'samples_all.tsv'
input_file2 = 'units_all.tsv'
output_file1 = 'samples.tsv'
output_file2 = 'units.tsv'

# Read the files
df1 = pd.read_csv(input_file1, sep='\t')
df2 = pd.read_csv(input_file2, sep='\t')

# Sample both where tissue in samples is body
selector = df1['tissue'] == 'body'
df1 = df1[selector]
df2 = df2[selector]

# Ensure there are enough rows
min_size = min(len(df1), len(df2))
if min_size < 25:
    raise ValueError("Files do not contain enough lines. Need at least 25 lines each.")

# Generate 24 unique random indices, starting from the second row
random_indices = random.sample(range(1, min_size), 24)


# Select the rows
selected_df1 = df1.iloc[random_indices]
selected_df2 = df2.iloc[random_indices]

# Write to new files
selected_df1.to_csv(output_file1, sep='\t', index=False)
selected_df2.to_csv(output_file2, sep='\t', index=False)
