import os
import sys
import re
import numpy as np

'''
import_temp_data.py

Imports the climate forcing temperature provided by Tamsin.

Parameters:
filename        Name of the file to import

Return:
data_dict       Dictionary of data using the header line as keys. The key for the data 
                is 'data'.

'''


def import_temp_data(filename):
    # Initialize the data dictionary structure
    data_dict = {}
    
    # Load the emulated data file
    with open(filename, 'r') as f:
        # Get the header info
        header_line = f.readline().rstrip()
        header = np.array(header_line.split(","))
        
        # Find the index of the first year dynamically
        first_year_idx = np.flatnonzero([h.isdigit() for h in header])[0]
        
        # Define the data indices
        data_idx = np.arange(first_year_idx, len(header))
        
        # Define the meta data indices
        meta_idx = np.arange(first_year_idx)
        
        # Extract the years from the header line
        data_dict['years'] = np.array([int(x) for x in header[data_idx]])
        
        # Initialize the data dictionary
        for i in meta_idx:
            data_dict[header[i]] = []
        data_dict['data'] = []
        
        # Loop through the lines in the file
        for line in f:
            # Get the line pieces
            line = line.rstrip()
            line_pieces = line.split(",")
            
            # Replace "NA" and empty strings with numpy nan
            numeric_line = []
            for i, value in enumerate(line_pieces):
                if value == "NA" or value == "":
                    numeric_line.append(np.nan)
                else:
                    try:
                        numeric_line.append(float(value))
                    except ValueError:
                        numeric_line.append(value)  # Keep as string if not convertible
            
            line_pieces = np.array(numeric_line, dtype=object)  # Use dtype=object to handle mixed types
            
            # Extract the data from the line
            data_dict['data'].append([line_pieces[i] for i in data_idx])
            
            # Extract the meta data
            for i in meta_idx:
                data_dict[header[i]].append(line_pieces[i])
    
    # Convert everything into numpy arrays
    for this_key in data_dict.keys():
        try:
            data_dict[this_key] = np.array(data_dict[this_key], dtype=float)
        except ValueError:
            data_dict[this_key] = np.array(data_dict[this_key], dtype=object)
    
    # Return the data dictionary
    return data_dict


if __name__ == "__main__":
    filename = "CLIMATE_FORCING_1850.csv"
    
    x = import_temp_data(filename)
    
    print(x['data'].shape)
    
    exit()
