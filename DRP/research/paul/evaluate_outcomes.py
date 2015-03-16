import sys

if __name__ == "__main__":
	if len(sys.argv) != 2:
		raise Exception("Not the right number of arguments")
	with open(sys.argv[1]) as outcome:
		for i in range(5): # waste the crap at the top of the output
			outcome.readline()
		max_conf = 0.0
		cnt = 0
		indices = []
		for line in outcome:
			if "+" not in line and "4" in line:
				conf = float(line.rstrip().split()[-1])
				if conf > max_conf:
					max_conf = conf
					indices = [cnt] 
				elif conf == max_conf:
					indices.append(cnt)
			cnt += 1

		if len(indices) == 0:
			print str(max_conf) + " EMPTY"
		else:
			print str(max_conf) + " " + ",".join([str(x) for x in indices])
