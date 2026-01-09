import xml.etree.ElementTree as ET
import csv
import sys

def extract_peptide_data(idxml_file, output_file=None):
    """
    Extract PeptideHit data from idXML file.
    
    Args:
        idxml_file: Path to the idXML file
        output_file: Optional output CSV file path
    """
    
    # Parse the XML file
    tree = ET.parse(idxml_file)
    root = tree.getroot()
    
    # Dictionary to store extracted data
    protein_data = {}
    
    # Find all ProteinHit elements
    for protein_hit in root.iter('ProteinHit'):
        # Extract basic attributes
        id = protein_hit.get('id')
        accession = protein_hit.get('accession')
        protein_data[id] = accession

    # List to store extracted data
    peptide_data = []
    
    # Find all PeptideHit elements
    for peptide_hit in root.iter('PeptideHit'):
        # Extract basic attributes
        score = peptide_hit.get('score', '')
        sequence = peptide_hit.get('sequence', '')
        protein_refs = peptide_hit.get('protein_refs', '')
        proteins = protein_refs.split(' ')
        protein_accessions = [protein_data[ref] for ref in proteins]
        protein_accession = ', '.join(protein_accessions)
        
        # Initialize data dictionary
        hit_data = {
            'score': score,
            'sequence': sequence,
            'protein': protein_accession,
            'target_decoy': '',
            'fragment_annotation': '',
            'NuXL_isXL': ''
        }
        
        # Extract UserParam values
        for user_param in peptide_hit.findall('UserParam'):
            param_name = user_param.get('name', '')
            param_value = user_param.get('value', '')
            
            if param_name == 'target_decoy':
                hit_data['target_decoy'] = param_value
            elif param_name == 'fragment_annotation':
                hit_data['fragment_annotation'] = param_value
            elif param_name == 'NuXL:isXL':
                hit_data['NuXL_isXL'] = param_value
            elif param_name == 'NuXL:score':
                hit_data['score'] = param_value
        
        peptide_data.append(hit_data)
    
    # Output results
    if output_file:
        # Write to CSV file
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['score', 'sequence', 'protein', 'target_decoy', 
                         'NuXL_isXL', 'fragment_annotation']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
            
            writer.writeheader()
            for data in peptide_data:
                writer.writerow(data)
        
        print(f"Data extracted to {output_file}")
    else:
        # Print to console
        print("Extracted PeptideHit data:")
        print("-" * 80)
        for i, data in enumerate(peptide_data, 1):
            print(f"Hit {i}:")
            for key, value in data.items():
                print(f"  {key}: {value}")
            print("-" * 40)
    
    return peptide_data

def main():
    if len(sys.argv) < 2:
        print("Usage: python extract_peptide_data.py <idxml_file> [output_tsv]")
        print("Example: python extract_peptide_data.py data.idXML output.tsv")
        return
    
    idxml_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    try:
        extract_peptide_data(idxml_file, output_file)
    except FileNotFoundError:
        print(f"Error: File '{idxml_file}' not found.")
    except ET.ParseError as e:
        print(f"Error parsing XML file: {e}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()