import subprocess

def parse_file(filepath):
    data = {}
    with open(filepath, 'r') as file:
        for line in file:
            parts = line.split()
            key = parts[0]
            # Check if there are multiple values for this key
            if len(parts) > 2:
                # Convert all to float and store as a tuple or list
                values = tuple(float(value) for value in parts[1:])
            else:
                # Try to convert the single value to float, otherwise keep as string
                try:
                    value = float(parts[1])
                    # Handle integer values
                    if value.is_integer():
                        value = int(value)
                except ValueError:
                    value = parts[1]
                values = value
            data[key] = values
    return data

filepath = 'option.txt'
parsed_data = parse_file(filepath)
print(parsed_data)

# Change "rodRadius" to add it 0.00002 while it is less than 0.0051
while parsed_data['rodRadius'] < 0.0045:
    parsed_data['rodRadius'] = parsed_data['rodRadius'] + 0.00002
    print(parsed_data)
    # Write the updated data to a new file
    with open('option.txt', 'w') as file:
        for key, value in parsed_data.items():
            # Check if the value is a tuple or list
            if isinstance(value, (tuple, list)):
                # Convert all values to strings and join them with a space
                value_str = ' '.join(str(v) for v in value)
            else:
                value_str = str(value)
            file.write(f'{key} {value_str}\n')

    # Run the simulation with the new rodRadius
    # Specify the path to your text file
    file_path = 'Commands.txt'

    print('Running the commands from the file:', file_path)
    # Read the content of the file
    with open(file_path, 'r') as file:
        # Read each line from the file
        for line in file:
            # Remove leading and trailing whitespaces
            command = line.strip()
            # Execute the command using subprocess
            subprocess.run(['bash', '-c', command])

    # Copy the lines 16951 to 17289 from the file "configData" and use them to update the entire file tied-R1
    with open('datafiles/configData_numFlagella_1_EA_135400_dt_0.001_totalTime_0.511_imc_1_friction_0.txt', 'r') as file:
        configData = file.readlines()
    with open('knot_configurations/tied-R1', 'r') as file:
        tiedR1 = file.readlines()
    with open('knot_configurations/tied-R1', 'w') as file:
        for i, line in enumerate(tiedR1):
            index_in_file = i + 16950
            if 1 <= index_in_file <= len(configData):
            # remove the first "0 " from the first line and the last line of tiedR1
                if i == 0:
                    file.write(configData[index_in_file][2:])
                elif i == len(tiedR1) - 1:
                    file.write(configData[index_in_file][2:])
                else:
                    file.write(configData[index_in_file])
