import pandas as pd

# Read the CSV file
df = pd.read_csv('/Users/satvikverma/Workspace/CSC-746/cp4/calc-speedup.csv')

# Strip leading and trailing spaces from column names
df.columns = df.columns.str.strip()

# Print the column names to verify
print("Columns in the DataFrame:", df.columns)

# Check if 'serial' and 'OMP' columns exist
if 'serial' in df.columns and 'OMP' in df.columns:
    df['speedup'] = df['serial'] / df['OMP']
    print(df)
else:
    print("Required columns 'serial' and/or 'OMP' are missing in the CSV file.")