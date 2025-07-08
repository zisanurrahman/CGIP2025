import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

# Load SMILES dataset
smiles_df = pd.read_csv('https://zenodo.org/records/13117525/files/unique_compounds_smiles.csv')

# Function to calculate Lipinski descriptors
def calculate_lipinski_params(mol):
    if mol is None:
        return {}  # Return an empty dictionary for invalid molecules
    
    params = {
        'Molecular Weight (Da)': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'Fraction SP3': Lipinski.FractionCSP3(mol),
        'Heavy Atom Count': Lipinski.HeavyAtomCount(mol),
        'N/O Count': Lipinski.NOCount(mol),
        'NH/OH Count': Lipinski.NHOHCount(mol),
        'Aliphatic Carbocycles': Lipinski.NumAliphaticCarbocycles(mol),
        'Aliphatic Heterocycles': Lipinski.NumAliphaticHeterocycles(mol),
        'Aliphatic Rings': Lipinski.NumAliphaticRings(mol),
        'Aromatic Carbocycles': Lipinski.NumAromaticCarbocycles(mol),
        'Aromatic Heterocycles': Lipinski.NumAromaticHeterocycles(mol),
        'Aromatic Rings': Lipinski.NumAromaticRings(mol),
        'Hydrogen Bond Acceptors': Lipinski.NumHAcceptors(mol),
        'Hydrogen Bond Donors': Lipinski.NumHDonors(mol),
        'Heteroatoms': Lipinski.NumHeteroatoms(mol),
        'Rotatable Bonds': Lipinski.NumRotatableBonds(mol),
        'Saturated Carbocycles': Lipinski.NumSaturatedCarbocycles(mol),
        'Saturated Heterocycles': Lipinski.NumSaturatedHeterocycles(mol),
        'Saturated Rings': Lipinski.NumSaturatedRings(mol),
        'Total Rings': Lipinski.RingCount(mol)
    }
    return params

# Convert SMILES to RDKit Mol objects and handle invalid SMILES
smiles_df['Mol'] = smiles_df['SMILES'].apply(lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None)

# Compute Lipinski properties for each molecule
lipinski_data = smiles_df['Mol'].apply(lambda mol: calculate_lipinski_params(mol) if mol else {})

# Convert the Lipinski results into a DataFrame (ignore None values)
lipinski_df = pd.DataFrame(list(lipinski_data))

# Merge the Lipinski parameters back into the original dataset
smiles_df = pd.concat([smiles_df, lipinski_df], axis=1)

# Drop the temporary 'Mol' column
smiles_df.drop(columns=['Mol'], inplace=True)

# Save the resulting DataFrame as a CSV file
output_file = 'lipinski_parameters_250308.csv'
smiles_df.to_csv(output_file, index=False)

# Print the path where the file was saved
print(f"File saved as: {output_file}")
