filepath = ("proteinGroups.txt")
import pandas as pd
df = pd.read_csv(filepath, delimiter= "\t")
print(df)