# Selenium Config
# Module for functions to deal with the simulation configuration files

def readConfigFile(filename):
	"""
	Reads the config.txt file containing relevant simulation information
	Inputs:
		filename - string to the path/name of the config file
	Outputs:
		settings - dictionary containing information from the file
	"""

	settings = {}

	with open(filename, 'r') as file:
		for line in file:
			if line[0] == '#' or line[0] == '\n':
				pass
			else:
				name, value = line.split('=')
				value = castType(value)
				settings[name] = value

	return settings

def castType(stringData):
	"""
	Attempts to cast data read as a string to the appropriate data type
	"""

	try:
		castData = int(stringData)
	except ValueError:
		try:
			castData = float(stringData)
		except ValueError:
			castData = stringData.split('\n')[0]

	return castData

if __name__ == '__main__':
	settings = readConfigFile('./config.txt')
	print(settings)
