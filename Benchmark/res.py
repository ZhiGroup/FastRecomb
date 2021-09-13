import sys

res = int(sys.argv[2])

end = 0
first = True
for line in open(sys.argv[1], 'r'):
	if first:
		first = False
		continue
	tokens = line.split()
	end = int(tokens[0])
end = end // res

fastRecomb = open(sys.argv[1], 'r')
fastRecomb.readline()
first = True
mapping = open('result.' + str(res) + '.map', 'w')
prevPos, prevMap = 0, 0
curPos, curMap = 0, 0
for i in range(end + 1):
	pos = i * res
	while curPos < pos:
		prevPos, prevMap = curPos, curMap
		tokens = fastRecomb.readline().split()
		curPos, curMap = float(tokens[0]), float(tokens[2])
	if first:
		first = False
		continue

	if curPos == pos:
		mapping.write(str(pos) + ' ' + str(curMap) + '\n')
	else:
		interMap = prevMap + (curMap - prevMap) / (curPos - prevPos) * (pos - prevPos)
		mapping.write(str(pos) + ' ' + str(interMap) + '\n')
	
