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

# Change "rodRadius" to add it 0.000025 while it is less than 0.006
while parsed_data['rodRadius'] < 0.006:
    parsed_data['rodRadius'] = parsed_data['rodRadius'] + 0.000025
    print(parsed_data)
    
    # Run the simulation with the new rodRadius
    # Specify the path to your text file
    file_path = 'Commands.txt'

    # Read the content of the file
    with open(file_path, 'r') as file:
        # Read each line from the file
        for line in file:
            # Remove leading and trailing whitespaces
            command = line.strip()

            # Execute the command using subprocess
            subprocess.run(['bash', '-c', command])
    
    # Copy the lines 16951 to 17289 from the file "configData" and use them to update the entire file tied-R1
    with open('configData_numFlagella_1_EA_135400_dt_0.001_totalTime_0.511_imc_1_friction_0.txt', 'r') as file:
        configData = file.readlines()
    with open('tied-R1', 'r') as file:
        tiedR1 = file.readlines()
    with open('tied-R1', 'w') as file:
        for i, line in enumerate(tiedR1):
            if 16951 <= i <= 17289:
                file.write(configData[i])
                # remove the first "0 " from the first line and the last line
                if i == 16951 or i == 17289:
                    file.write(line[2:])
            else:
                file.write(line)



