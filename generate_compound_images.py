import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import base64
import io
import json
import os

# Load CSV exported from R
df = pd.read_csv("/Users/zisan/Documents/GitHub/CGIP2025/data/lipinski.csv")

compound_images = {}

for _, row in df.iterrows():
    compound = row["Compounds"]
    smiles = row["SMILES"]

    if pd.notna(smiles) and len(smiles) > 0:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(400, 400))
                buffered = io.BytesIO()
                img.save(buffered, format="PNG")
                img_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")
                compound_images[compound] = f"data:image/png;base64,{img_base64}"
        except Exception as e:
            print(f"Error processing {compound}: {e}")

# Save to JSON
output_path = "/Users/zisan/Documents/GitHub/CGIP2025/data/compound_images.json"
with open(output_path, "w") as f:
    json.dump(compound_images, f)

print(f"Saved base64 images to: {output_path}")

