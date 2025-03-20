import streamlit as st
import subprocess
import sys
import pandas as pd
import os
import base64
import re
from openpyxl import Workbook
from openpyxl.drawing.image import Image as XLImage
from openpyxl.utils import get_column_letter

# Check if rdkit is installed
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
except ImportError:
    st.error("""
    RDKit is not installed. Please follow these steps:
    
    1. Run in terminal:
    ```bash
    conda install -c conda-forge rdkit
    ```
    
    2. Restart the application after installation.
    """)
    st.stop()

from PIL import Image
import io
import zipfile

def is_valid_smiles(smiles):
    """Check if string is a valid SMILES"""
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        return mol is not None
    except:
        return False

def find_smiles_column(df):
    """Automatically identify column containing SMILES"""
    for column in df.columns:
        # Check if column name contains 'smiles' (case insensitive)
        if 'smiles' in column.lower():
            # Verify if column contains valid SMILES
            valid_count = sum(1 for x in df[column] if is_valid_smiles(x))
            if valid_count > len(df) * 0.5:  # If more than 50% values are valid SMILES
                return column
    return None

def save_to_excel(df, smiles_column, output_path):
    """Save data to Excel file with structure images"""
    wb = Workbook()
    ws = wb.active
    
    # Write column headers
    for col, header in enumerate(df.columns, 1):
        ws.cell(row=1, column=col, value=header)
    
    # Create temporary directory for images
    if not os.path.exists('temp_images'):
        os.makedirs('temp_images')
    
    # Process each row
    for row_idx, row in df.iterrows():
        excel_row = row_idx + 2  # Excel rows start at 1, header takes row 1
        
        # Write data
        for col_idx, value in enumerate(row, 1):
            ws.cell(row=excel_row, column=col_idx, value=value)
        
        # Process SMILES and add image
        smiles = str(row[smiles_column])
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                # Generate molecular image
                img = Draw.MolToImage(mol, size=(400, 400))
                img_path = f'temp_images/molecule_{row_idx}.png'
                img.save(img_path)
                
                # Add image to Excel
                img = XLImage(img_path)
                # Adjust image size
                img.width = 200
                img.height = 200
                # Add image after SMILES column
                ws.add_image(img, f'{get_column_letter(len(df.columns) + 1)}{excel_row}')
        except Exception as e:
            st.warning(f"Error processing SMILES in row {row_idx + 1}: {str(e)}")
    
    # Adjust column widths
    for col in range(1, len(df.columns) + 2):  # +2 because of image column
        ws.column_dimensions[get_column_letter(col)].width = 15
    
    # Save Excel file
    wb.save(output_path)
    
    # Clean up temporary files
    for file in os.listdir('temp_images'):
        os.remove(os.path.join('temp_images', file))
    os.rmdir('temp_images')

st.set_page_config(page_title="SMILES Structure Viewer", layout="wide")

st.title("SMILES Structure Viewer")
st.write("Enter SMILES string or upload CSV file to view chemical structures")

# Create two tabs
tab1, tab2 = st.tabs(["Single Conversion", "Batch Conversion"])

with tab1:
    # Create two-column layout
    col1, col2 = st.columns(2)

    with col1:
        smiles_input = st.text_input("Enter SMILES string", "C1=CC=CC=C1")
        
        if smiles_input:
            try:
                # Create RDKit molecule object
                mol = Chem.MolFromSmiles(smiles_input)
                if mol is None:
                    st.error("Invalid SMILES string")
                else:
                    # Generate molecular image
                    img = Draw.MolToImage(mol, size=(400, 400))
                    
                    # Convert image to bytes
                    img_byte_arr = io.BytesIO()
                    img.save(img_byte_arr, format='PNG')
                    img_byte_arr = img_byte_arr.getvalue()
                    
                    # Display image
                    st.image(img_byte_arr, caption="Chemical Structure", use_container_width=True)
                    
                    # Add download button
                    st.download_button(
                        label="Download Structure Image",
                        data=img_byte_arr,
                        file_name="molecule.png",
                        mime="image/png"
                    )
            except Exception as e:
                st.error(f"Error processing: {str(e)}")

    with col2:
        st.markdown("""
        ### Instructions
        1. Enter SMILES string in the input box on the left
        2. System will automatically display the chemical structure
        3. Click "Download Structure Image" to save the structure
        
        ### Example SMILES
        - Benzene: `C1=CC=CC=C1`
        - Ethanol: `CCO`
        - Aspirin: `CC(=O)OC1=CC=CC=C1C(=O)O`
        """)

with tab2:
    st.markdown("""
    ### Batch Conversion Instructions
    1. Upload a CSV file containing SMILES column
    2. System will automatically identify the SMILES column
    3. System will generate and display all structure images
    4. You can download individual structure images or export to Excel
    """)
    
    uploaded_file = st.file_uploader("Upload CSV file", type=['csv'])
    
    if uploaded_file is not None:
        try:
            # Read CSV file
            df = pd.read_csv(uploaded_file)
            
            # Display data preview
            st.subheader("Data Preview")
            st.dataframe(df.head())
            
            # Automatically identify SMILES column
            smiles_column = find_smiles_column(df)
            
            if smiles_column:
                st.success(f"Automatically identified SMILES column: {smiles_column}")
            else:
                smiles_column = st.selectbox("SMILES column not automatically identified, please select manually", df.columns)
            
            if smiles_column:
                # Create table data
                table_data = []
                
                # Process each SMILES
                success_count = 0
                error_count = 0
                error_smiles = []
                
                for idx, row in df.iterrows():
                    smiles = str(row[smiles_column])
                    try:
                        mol = Chem.MolFromSmiles(smiles)
                        if mol is not None:
                            # Generate molecular image
                            img = Draw.MolToImage(mol, size=(400, 400))
                            
                            # Convert image to bytes
                            img_byte_arr = io.BytesIO()
                            img.save(img_byte_arr, format='PNG')
                            img_byte_arr = img_byte_arr.getvalue()
                            
                            # Create table row data
                            row_data = {}
                            for col in df.columns:
                                if col == smiles_column:
                                    row_data['Structure'] = img_byte_arr
                                else:
                                    row_data[col] = row[col]
                            
                            table_data.append(row_data)
                            success_count += 1
                        else:
                            error_count += 1
                            error_smiles.append(smiles)
                    except Exception as e:
                        error_count += 1
                        error_smiles.append(smiles)
                
                # Display statistics
                st.success(f"Processing complete! Success: {success_count}, Failed: {error_count}")
                
                # Display failed SMILES if any
                if error_smiles:
                    with st.expander("View Failed SMILES"):
                        st.write("The following SMILES failed to process:")
                        for smiles in error_smiles:
                            st.code(smiles)
                
                # Display data in table format
                if table_data:
                    st.subheader("Structure Preview")
                    for row in table_data:
                        # Create column layout
                        cols = st.columns([1] + [2] * (len(row) - 1))
                        
                        # Display structure image
                        with cols[0]:
                            st.image(row['Structure'], use_container_width=True)
                        
                        # Display other data
                        for i, (key, value) in enumerate(row.items()):
                            if key != 'Structure':
                                with cols[i]:
                                    st.write(f"**{key}:**")
                                    st.write(value)
                        
                        # Add separator
                        st.markdown("---")
                
                # Add Excel export functionality
                if st.button("Export to Excel"):
                    try:
                        # Create temporary Excel file
                        excel_path = "molecules_with_structures.xlsx"
                        save_to_excel(df, smiles_column, excel_path)
                        
                        # Read Excel file and create download button
                        with open(excel_path, 'rb') as f:
                            st.download_button(
                                label="Download Excel File",
                                data=f,
                                file_name="molecules_with_structures.xlsx",
                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                            )
                        
                        # Delete temporary file
                        os.remove(excel_path)
                    except Exception as e:
                        st.error(f"Error exporting to Excel: {str(e)}")
                
        except Exception as e:
            st.error(f"Error processing CSV file: {str(e)}") 