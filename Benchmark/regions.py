import sys

end = 0 
for line in open(sys.argv[1], 'r'):
	tokens = line.split()
	end = int(tokens[0])

out = open(sys.argv[1] + '.color', 'w')
for line in open(sys.argv[1], 'r'):
	tokens = line.split()
	pos, hapMap, IBDMap = float(tokens[0]), float(tokens[1]), float(tokens[2])
	color = 'black'
	if pos <= 2e6 or pos + 2e6 >= end: color = 'red'
	elif pos < 5e6 or pos + 5e6 >= end: color = 'blue'
	out.write(str(pos) + ' ' + str(hapMap) + ' ' + str(IBDMap) + ' ' + color + '\n')	

