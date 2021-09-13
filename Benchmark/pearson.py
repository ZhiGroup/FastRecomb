import sys
from scipy import stats

red_x, red_y, blue_x, blue_y, black_x, black_y, overall_x, overall_y = [], [], [], [], [], [], [], []
for line in open(sys.argv[1], 'r'):
	tokens = line.split()
	pos, hapMap, IBDMap, color = float(tokens[0]), float(tokens[1]), float(tokens[2]), tokens[3]
	overall_x.append(hapMap); overall_y.append(IBDMap)
	if color == 'red':
		red_x.append(hapMap)
		red_y.append(IBDMap)
	if color == 'blue':
		blue_x.append(hapMap)
		blue_y.append(IBDMap)
	if color == 'black':
		black_x.append(hapMap)
		black_y.append(IBDMap)

r, p = stats.pearsonr(red_x, red_y)
print(str(r) + ' ', end = '')
r, p = stats.pearsonr(blue_x, blue_y)
print(str(r) + ' ', end = '')
r, p = stats.pearsonr(black_x, black_y)
print(str(r) + ' ', end = '')
r, p = stats.pearsonr(overall_x, overall_y)
print(str(r) + ' ')
