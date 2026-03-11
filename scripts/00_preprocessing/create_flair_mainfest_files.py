import json
import os
import re

# --- Configuration Variables ---
# NOTE: Ensure config.json is accessible at this path when you run the script.
CONFIG_FILE = "config.json" 
OUTPUT_DIR = "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/metadata/" 

MANIFEST1_NAME = "flair_combine_manifest.tsv"
MANIFEST2_NAME = "flair_quantify_manifest.tsv"

# List of keys to exclude when parsing the JSON, as they appear to be comments
# and cause the standard Python json parser to fail if they are not the last 
# item in an object (due to trailing commas) or if they contain invisible characters.
COMMENT_KEYS = [
    "Commment_Input_Output_Directories", # Note the misspelling 'Commment'
    "Comment_Reference",
    "Comment_Sample_Info"
]

def load_json_with_comments(file_path):
    """
    Loads a JSON file, filtering out lines containing known comment keys
    and fixing non-standard syntax (like single quotes) for compatibility 
    with Python's strict json parser.
    """
    try:
        # 1. Read the file content
        # Using 'utf-8-sig' helps strip the Byte Order Mark (BOM) if present.
        with open(file_path, 'r', encoding='utf-8-sig') as f:
            lines = f.readlines()
            
        filtered_lines = []
        for line in lines:
            # Check if the line contains any of the comment keys
            is_comment_line = any(f'"{key}"' in line for key in COMMENT_KEYS)
            
            if not is_comment_line:
                filtered_lines.append(line)
        
        # 2. Join the filtered lines and clean up residual issues
        content = "".join(filtered_lines)
        
        # FIX 1: Remove double/triple commas caused by deleting a line that had a comma.
        content = re.sub(r',\s*,', ',', content)

        # FIX 2: Remove the trailing comma before a closing bracket/brace.
        content = re.sub(r',\s*([\]}])', r'\1', content, flags=re.DOTALL)

        # FIX 3 (New): Replace single quotes used for array elements with double quotes.
        # This targets non-standard JSON like: ['item1', 'item2']
        # This is the most likely cause of the "Expecting value" error.
        content = re.sub(r"'([^']*)'", r'"\1"', content) 
        
        # 3. Parse the cleaned content
        config = json.loads(content)
        return config

    except FileNotFoundError:
        print(f"ERROR: Config file not found at {file_path}")
        return None
    except json.JSONDecodeError as e:
        print(f"ERROR: Still failed to decode JSON after fixing quotes. Original error: {e}")
        print("Please check config.json manually for non-standard characters or misplaced commas.")
        return None


def generate_manifest_files(config_file, output_dir):
    """
    Parses the config data (from the filtered JSON) to create two manifest files.
    """
    # Load configuration using the custom function
    config = load_json_with_comments(config_file)
    if config is None:
        return

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Get base directories from config
    collapsed_dir = config.get("collapsed", "")
    fastq_dir = config.get("fastq", "")
    
    # Initialize content for both manifest files
    manifest1_lines = []
    manifest2_lines = []

    # Get the list of samples to iterate over
    sample_names = config.get("sample_names", [])

    # 1. Iterate through samples and build manifest lines
    for full_sample_name in sample_names:
        sample_info = config.get(full_sample_name)
        if not sample_info:
            print(f"Warning: Information for sample {full_sample_name} not found in config.")
            continue
        
        # Core info extracted from the per-sample dictionaries
        sample_id = sample_info.get("ID")  
        type = "isoform" # type (isoform or fusionisoform)
        group = sample_info.get("group")       
        batch = sample_info.get("PU")
        
        # --- Manifest 1: flair combine (ID, BED path, FASTA path) ---
        # Used to merge individual transcriptomes (after collapse)
        bed_filename = f"{full_sample_name}.isoforms.bed"
        fa_filename = f"{full_sample_name}.isoforms.fa"
        map_filename = f"{full_sample_name}.combined.isoform.read.map.txt "

        bed_path = os.path.join(collapsed_dir, bed_filename)
        fa_path = os.path.join(collapsed_dir, fa_filename)
        map_path = os.path.join(collapsed_dir, map_filename)
   
        # Output format: ID \t BED_path \t FASTA_path
        manifest1_lines.append(f"{sample_id}\t{type}\t{bed_path}\t{fa_path}\t{map_path}")
        
        # --- Manifest 2: flair quantify (ID, Condition, Batch, FASTQ path) ---
        # Used for expression quantification
        fastq_filename = f"{full_sample_name}.fastq.gz"
        fastq_path = os.path.join(fastq_dir, fastq_filename)

        # Output format: ID \t Condition \t Batch \t FASTQ_path
        manifest2_lines.append(f"{sample_id}\t{group}\t{batch}\t{fastq_path}")

    # 2. Write Manifest 1 (flair combine input)
    manifest1_path = os.path.join(output_dir, MANIFEST1_NAME)
    with open(manifest1_path, 'w') as f:
        f.write('\n'.join(manifest1_lines) + '\n') 
    print(f"✅ Success! Wrote Manifest 1 (flair_combine_manifest.tsv) to: {manifest1_path}")

    # 3. Write Manifest 2 (flair quantify input)
    manifest2_path = os.path.join(output_dir, MANIFEST2_NAME)
    with open(manifest2_path, 'w') as f:
        f.write('\n'.join(manifest2_lines) + '\n') 
    print(f"✅ Success! Wrote Manifest 2 (flair_quantify_manifest.tsv) to: {manifest2_path}")

# --- Execute the main function ---
generate_manifest_files(CONFIG_FILE, OUTPUT_DIR)

