#!/usr/bin/env python
"""
Script to fill in a Template.csv based on input Excel/CSV files
Standalone version with own configuration
"""

import os
import sys
import json
import platform
import contextlib
import pandas as pd
import logging
import re
import datetime

# Configure logging - only errors and warnings by default
logging.basicConfig(
    level=logging.WARNING,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("TemplateCreator")

# Enable debug mode with --debug flag
if len(sys.argv) > 1 and '--debug' in sys.argv:
    logger.setLevel(logging.DEBUG)

# ============ Configuration ============

CONFIG_DIR = os.path.join(os.path.expanduser("~"), ".template_creator")
CONFIG_FILE = os.path.join(CONFIG_DIR, "config.json")

@contextlib.contextmanager
def _silence_stderr():
    """Temporarily redirect stderr to /dev/null to suppress wxPython warnings."""
    if platform.system() == "Darwin":
        old_fd = os.dup(2)
        try:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, 2)
            os.close(devnull)
            yield
        finally:
            os.dup2(old_fd, 2)
            os.close(old_fd)
    else:
        yield

def get_config():
    """Load configuration from file."""
    config = {}
    try:
        if os.path.exists(CONFIG_FILE):
            with open(CONFIG_FILE, 'r') as f:
                config = json.load(f)
    except Exception as e:
        logger.debug(f"Error loading config: {str(e)}")
    return config

def save_config(config):
    """Save configuration to file."""
    try:
        if not os.path.exists(CONFIG_DIR):
            os.makedirs(CONFIG_DIR, mode=0o755, exist_ok=True)
        
        with open(CONFIG_FILE, 'w') as f:
            json.dump(config, f)
            f.flush()
            os.fsync(f.fileno())
        
    except Exception as e:
        logger.debug(f"Error saving config: {str(e)}")

def find_default_directory():
    """Find a sensible default directory based on the OS."""
    potential_paths = []
    home_dir = os.path.expanduser("~")
    
    # Add OS-specific common locations
    if sys.platform == 'win32':  # Windows
        potential_paths.extend([
            os.path.join(home_dir, "Downloads"),
            os.path.join(home_dir, "Documents"),
            os.path.join(home_dir, "Desktop"),
            "C:\\Data"
        ])
    elif sys.platform == 'darwin':  # macOS
        potential_paths.extend([
            os.path.join(home_dir, "Downloads"),
            os.path.join(home_dir, "Documents"),
            os.path.join(home_dir, "Desktop"),
            os.path.join(home_dir, "Library", "Mobile Documents"),
            "/Volumes"
        ])
    else:  # Linux/Unix
        potential_paths.extend([
            os.path.join(home_dir, "Downloads"),
            os.path.join(home_dir, "Documents"),
            os.path.join(home_dir, "Desktop"),
            "/mnt",
            "/media"
        ])
    
    potential_paths.insert(0, os.getcwd())
    
    for path in potential_paths:
        if os.path.exists(path) and os.path.isdir(path):
            return path
    
    return home_dir

def select_file(default_path=None, title="Select file", wildcard="CSV files (*.csv)|*.csv"):
    """Display a file selection dialog and return the selected path."""
    config = get_config()
    last_dir = config.get('last_directory')
    
    if default_path is None:
        if last_dir and os.path.isdir(last_dir):
            parent_dir = os.path.dirname(last_dir)
            if parent_dir and os.path.isdir(parent_dir):
                default_path = parent_dir
            else:
                default_path = last_dir
        else:
            default_path = find_default_directory()
    
    # Try using wxPython dialog
    try:
        import wx
        
        _wx_app = wx.App(False)
        
        with _silence_stderr():
            style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
            
            dlg = wx.FileDialog(
                None,
                message=title,
                defaultDir=default_path if default_path else "",
                defaultFile="",
                wildcard=wildcard,
                style=style
            )
            
            file_path = None
            if dlg.ShowModal() == wx.ID_OK:
                file_path = dlg.GetPath()
                
                config['last_directory'] = os.path.dirname(file_path)
                save_config(config)
            
            dlg.Destroy()
            
            if platform.system() == "Darwin":
                try:
                    import AppKit
                    NSApplication = AppKit.NSApplication.sharedApplication()
                    NSApplication.setActivationPolicy_(AppKit.NSApplicationActivationPolicyAccessory)
                except:
                    pass
        
        _wx_app = None
        
        if file_path:
            return file_path
            
    except ImportError:
        logger.info("wxPython not available, falling back to console input")
    except Exception as e:
        logger.error(f"Error in GUI dialog: {str(e)}")
        logger.info("Falling back to console input")
    
    # Fallback to manual input
    print(f"\nEnter full path to file ({wildcard.split('|')[0]}):")
    if default_path:
        print(f"[Previous directory: {default_path}]")
    file_path = input("> ").strip()
    
    if not os.path.isfile(file_path):
        logger.error(f"Invalid file: {file_path}")
        print(f"Error: '{file_path}' is not a valid file.")
        return None
    
    config['last_directory'] = os.path.dirname(file_path)
    save_config(config)
    
    return file_path

# ============ Template Functions ============

def create_template_header():
    """Create the header lines for the template file."""
    now = datetime.datetime.now()
    header_lines = [
        [
            "ddplate - DO NOT MODIFY THIS LINE",
            "Version=1",
            "ApplicationName=QX Manager Standard Edition",
            "ApplicationVersion=2.3.0.32",
            "ApplicationEdition=ResearchEmbedded",
            "User=\\QX User",
            f"CreatedDate={now.strftime('%m/%d/%Y %H:%M:%S')}",
            ""
        ],
        [""],
        ["PlateSize=GCR96"],
        ["PlateNotes="],
        [
            "Well", "Perform Droplet Reading", "ExperimentType", "Sample description 1",
            "Sample description 2", "Sample description 3", "Sample description 4",
            "SampleType", "SupermixName", "AssayType", "TargetName", "TargetType",
            "Signal Ch1", "Signal Ch2", "Reference Copies", "Well Notes", "Plot?",
            "RdqConversionFactor"
        ]
    ]
    return header_lines

def get_template_row_template():
    """Get a default template row structure."""
    return [
        "",  # Well - to be filled
        "Yes",  # Perform Droplet Reading
        "Direct Quantification (DQ)",  # ExperimentType
        "",  # Sample description 1 - to be filled
        "",  # Sample description 2 - to be filled
        "",  # Sample description 3 - to be filled
        "",  # Sample description 4 - to be filled
        "Unknown",  # SampleType - to be filled
        "ddPCR Supermix for Probes (No dUTP)",  # SupermixName
        "Probe Mix Triplex",  # AssayType
        "1",  # TargetName (Concentration)
        "Unknown",  # TargetType
        "FAM",  # Signal Ch1
        "HEX",  # Signal Ch2
        "",  # Reference Copies
        "",  # Well Notes
        "False",  # Plot?
        ""  # RdqConversionFactor
    ]

def detect_control_type(sample_name):
    """
    Detect if a sample is a control type based on its name.
    
    Args:
        sample_name (str): The sample description
        
    Returns:
        str: "NegCtrl", "PosCtrl", or "Unknown"
    """
    if not sample_name:
        return "Unknown"
    
    name_lower = str(sample_name).lower()
    
    neg_patterns = [
        r'negative\s*control',
        r'neg\s*ctrl',
        r'nc',
        r'blank',
        r'water'
    ]
    
    pos_patterns = [
        r'positive\s*control',
        r'pos\s*ctrl',
        r'pc'
    ]
    
    for pattern in neg_patterns:
        if re.search(pattern, name_lower):
            return "NegCtrl"
    
    for pattern in pos_patterns:
        if re.search(pattern, name_lower):
            return "PosCtrl"
    
    return "Unknown"

def calculate_wells_per_row(total_rows):
    """
    Calculate how many wells to fill per row based on total number of samples.
    
    Args:
        total_rows (int): Total number of samples
        
    Returns:
        int: Number of wells per row
    """
    return total_rows // 8

def read_input_file(file_path):
    """
    Read CSV or Excel file and extract sample information.
    
    Args:
        file_path (str): Path to input file
        
    Returns:
        pd.DataFrame: DataFrame with sample information
    """
    file_extension = os.path.splitext(file_path)[1].lower()
    
    try:
        if file_extension == '.csv':
            # Read without headers first to check if file has them
            df_check = pd.read_csv(file_path, header=None, nrows=1)
            has_header = any(df_check.iloc[0].astype(str).str.contains('Sample description', case=False))
            
            if has_header:
                df = pd.read_csv(file_path)
            else:
                # Read without headers and assign column names
                df = pd.read_csv(file_path, header=None)
                column_names = [f"Sample description {i+1}" for i in range(len(df.columns))]
                df.columns = column_names
        
        elif file_extension in ['.xlsx', '.xls']:
            # Read Excel file first line to check for headers
            df_check = pd.read_excel(file_path, header=None, nrows=1)
            has_header = any(df_check.iloc[0].astype(str).str.contains('Sample description', case=False))
            
            if has_header:
                df = pd.read_excel(file_path)
            else:
                # Read without headers and assign column names
                df = pd.read_excel(file_path, header=None)
                column_names = [f"Sample description {i+1}" for i in range(len(df.columns))]
                df.columns = column_names
        else:
            raise ValueError(f"Unsupported file format: {file_extension}")
        
        return df
    except Exception as e:
        logger.error(f"Error reading file: {str(e)}")
        raise

def process_template(input_file_path):
    """
    Process the input file and generate the filled template.
    
    Args:
        input_file_path (str): Path to input file
    """
    # Read input file
    df = read_input_file(input_file_path)
    
    total_rows = len(df)
    wells_per_row = calculate_wells_per_row(total_rows)
    
    print(f"Processing {total_rows} samples...")
    
    # Check if number of rows is a multiple of 8
    if total_rows % 8 != 0:
        print("--------------")
        print(f"WARNING: Number of samples ({total_rows}) is not a multiple of 8!")
        print(f"Some wells could not be filled.")
        print("--------------")
    
    print(f"Wells per row: {wells_per_row}")
    
    # Generate output filename - same as input but with .csv extension
    input_dir = os.path.dirname(input_file_path)
    input_filename = os.path.basename(input_file_path)
    name_without_ext = os.path.splitext(input_filename)[0]
    output_file_path = os.path.join(input_dir, name_without_ext + ".csv")
    
    # Create header
    output_lines = create_template_header()
    
    # Fill all wells up to 12 columns
    for row_letter in 'ABCDEFGH':
        for col_num in range(1, 13):
            well_id = f"{row_letter}{col_num:02d}"
            
            # Calculate the data index using column-major order
            row_idx = ord(row_letter) - ord('A')
            col_idx = col_num - 1
            
            # Column-major order: data_idx = column * 8 + row
            data_idx = col_idx * 8 + row_idx
            
            # Check if this well has data
            if data_idx < total_rows:
                # This well has data
                row = df.iloc[data_idx]
                
                # Extract sample descriptions from the known columns
                sample_desc = ""
                additional_desc = ["", "", ""]
                
                # Get values from columns
                if "Sample description 1" in df.columns:
                    sample_desc = str(row.get("Sample description 1", ""))
                    
                    for i in range(2, 5):
                        col_name = f"Sample description {i}"
                        if col_name in df.columns and pd.notna(row.get(col_name)):
                            additional_desc[i-2] = str(row.get(col_name))
                
                # Fallback if no descriptive columns found
                if not sample_desc and len(df.columns) > 0:
                    sample_desc = str(row.iloc[0]) if pd.notna(row.iloc[0]) else ""
                    for i in range(1, min(4, len(df.columns))):
                        if i < len(row) and pd.notna(row.iloc[i]):
                            additional_desc[i-1] = str(row.iloc[i])
                
                # Determine control type
                control_type = detect_control_type(sample_desc)
                
                # Create row based on template with data
                output_row = get_template_row_template()
                
                # Fill in the data
                output_row[0] = well_id  # Well
                output_row[3] = sample_desc  # Sample description 1
                output_row[4] = additional_desc[0]  # Sample description 2
                output_row[5] = additional_desc[1]  # Sample description 3
                output_row[6] = additional_desc[2]  # Sample description 4
                output_row[7] = control_type  # SampleType
            else:
                # Create empty row
                output_row = [
                    well_id,  # Well
                    "No",     # Perform Droplet Reading
                    "",       # ExperimentType
                    "",       # Sample description 1
                    "",       # Sample description 2
                    "",       # Sample description 3
                    "",       # Sample description 4
                    "",       # SampleType
                    "",       # SupermixName
                    "",       # AssayType
                    "",       # TargetName
                    "",       # TargetType
                    "",       # Signal Ch1
                    "",       # Signal Ch2
                    "",       # Reference Copies
                    "",       # Well Notes
                    "",       # Plot?
                    ""        # RdqConversionFactor
                ]
            
            output_lines.append(output_row)
    
    # Create DataFrame and save
    import csv
    
    with open(output_file_path, 'w', newline='') as f:
        writer = csv.writer(f)
        for line in output_lines:
            writer.writerow(line)
    
    print(f"Template saved to: {output_file_path}")
    print("")
    
    # Delete the input file only if the total rows is a multiple of 8
    if total_rows % 8 == 0:
        try:
            os.remove(input_file_path)
        except Exception as e:
            logger.error(f"Could not delete input file: {str(e)}")
            print(f"Warning: Could not delete input file: {str(e)}")
    else:
        print("Input file kept!")
        print("")
    
    return output_file_path

def main():
    """Main function to run the template creator."""
    print("\n=== Template Creator ===")
    
    config = get_config()
    default_path = config.get('last_directory')
    if not default_path:
        default_path = find_default_directory()
    
    input_file = select_file(
        default_path=default_path,
        title="Select input file (CSV or Excel)",
        wildcard="Supported files (*.csv;*.xlsx)|*.csv;*.xlsx|CSV files (*.csv)|*.csv|Excel files (*.xlsx)|*.xlsx"
    )
    
    if not input_file:
        print("No file selected. Exiting.")
        return
    
    try:
        output_file = process_template(input_file)
    
    except Exception as e:
        logger.error(f"Error processing file: {str(e)}")
        print(f"\nError processing file: {str(e)}")
        print("Please check the file format and try again.")

if __name__ == "__main__":
    main()