import sys

out = open(sys.argv[1] + '.recomb', 'w')
prevPos, prevHapMap, prevIBDMap = -1, -1, -1
for line in open(sys.argv[1], 'r'):
	tokens = line.split()
	pos, hapMap, IBDMap = float(tokens[0]), float(tokens[1]), float(tokens[2])
	if prevPos != -1:
		hapMapRecomb = (hapMap - prevHapMap) / (pos - prevPos)
		IBDMapRecomb = (IBDMap - prevIBDMap) / (pos - prevPos)
		out.write(str(int(prevPos)) + ' ' + str(hapMapRecomb) + ' ' + str(IBDMapRecomb) + '\n')

	prevPos, prevHapMap, prevIBDMap = pos, hapMap, IBDMap
out.write(str(int(prevPos)) + ' 0 0\n')
