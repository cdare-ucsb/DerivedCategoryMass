import pandas as pd
import re
import os

def parse_file(filename):
    with open(filename, 'r') as file:
        data = file.read()

    # Split data into blocks (each entry starts with "Num")
    blocks = data.split('Num    : ')[1:]

    parsed_data = []
    
    for block in blocks:
        entry = {}
        lines = block.strip().split('\n')

        # Extract simple fields
        entry['Num'] = int(re.search(r'\d+', lines[0]).group())
        entry['NumPs'] = int(re.search(r'\d+', lines[1]).group())
        entry['NumPol'] = int(re.search(r'\d+', lines[2]).group())
        entry['Eta'] = int(re.search(r'-?\d+', lines[3]).group())
        entry['H11'] = int(re.search(r'\d+', lines[4]).group())
        entry['H21'] = int(re.search(r'\d+', lines[5]).group())

        # Extract C2 (list of integers)
        entry['C2'] = [int(x) for x in re.findall(r'\d+', lines[6])]

        # Extract Redun (matrix)
        redun_start = 7
        redun_matrix = []
        while redun_start < len(lines) and '{' in lines[redun_start]:
            redun_matrix.append([int(x) for x in re.findall(r'-?\d+', lines[redun_start])])
            redun_start += 1

        entry['Redun'] = redun_matrix
        
        parsed_data.append(entry)

    # Convert the list of dictionaries to a Pandas DataFrame
    df = pd.DataFrame(parsed_data)
    return df



if __name__ == "__main__":
    filename = 'data/cicylist.txt'
    df = parse_file(filename)
    print(df.head(6))  # Display the first 6 entries
