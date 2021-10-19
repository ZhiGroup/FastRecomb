import sys

hapMap = open('sim1M_chr20.gmap', 'r')
hapMap.readline() # remove header line

mapping = open(sys.argv[1] + '.comp', 'w')
prevPos, prevMap = 0, 0
curPos, curMap = 0, 0
for line in open(sys.argv[1], 'r'):
	tokens = line.split()
	if len(tokens) == 0: continue
	pos, gmap = float(tokens[0]), float(tokens[1])	
	
	while curPos < pos:
		tokens2 = hapMap.readline().split()
		if len(tokens2) == 0: break # need to extrapolate
		prevPos, prevMap = curPos, curMap
		curPos, curMap = float(tokens2[1]), float(tokens2[3])
	if curPos == pos:
		mapping.write(str(int(pos)) + ' ' + str(curMap) + ' ' + str(gmap) + '\n')
	else:
		interMap = prevMap + (curMap - prevMap) / (curPos - prevPos) * (pos - prevPos)	
		mapping.write(str(int(pos)) + ' ' + str(interMap) + ' ' + str(gmap) + '\n')
	
