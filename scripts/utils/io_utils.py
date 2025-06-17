import pandas as pd
#Load, Save and preprocess data 

def load_data(file_path):
    """
    Load data from a CSV file into a pandas DataFrame.
    
    Parameters:
    file_path (str): Path to the CSV file containing the data.
    
    Returns:
    pd.DataFrame: DataFrame containing the loaded data.
    """
    try:
        df = pd.read_csv(file_path)
        print(f"Data loaded successfully from {file_path}")
        return df
    
    except FileExistsError as e:
        print(f"Error loading data: {e}")
        return None
    
    except FileExistsError as e:
        print(f"File not found: {e}")
        return None

def insertTime(data):
    """
    Insert a new column 'Time' into the DataFrame, converting 'Time (ms)' to seconds.

    Parameters:
    data (Data frame): DataFrame containing the 'Time (ms)' Time points as raws and ROI (cell) column.

    Returns:
    pd.DataFrame: DataFrame with the new 'Time' column inserted as the first column.
    """
    if 'Time (ms)' not in data.columns:
        print("Column 'Time (ms)' not found in DataFrame.")
        return data
    if 'Time' in data.columns:
        print("Column 'Time' already exists in DataFrame. Skipping insertion.")
        return data
    
    # Convert 'Time (ms)' to seconds and insert as a new column
    data['Time'] = data['Time (ms)'] / 1000.0
    # Move 'Time' column to the second position

    if 'Time' in data.columns:
        data = data[['Time'] + [col for col in data.columns if col != 'Time']]
    else:        
        print("Failed to insert 'Time' column.")
        return data
    print("Column 'Time' inserted successfully, retun Data Frame object")   

    return data

def cleanVarNames(df):
    """ Clean variable names in the DataFrame by removing unwanted characters and renaming columns.

    Parameters:
    df (pd.DataFrame): DataFrame containing the raw data.  

    Returns:
    pd.DataFrame: DataFrame with cleaned variable names.
    """
    # Remove unwanted characters from column names
    df.drop([var for var in df.columns if var[0] == 'U'],axis = 1, inplace = True)
    df.drop(['Frame','Time (ms)'], axis=1, inplace=True)
    return df
