import sys

if __name__ == "__main__":
	if len(sys.argv) != 2:
		raise Exception("Not the right number of arguments")
	with open(sys.argv[1]) as outcome:
		for i in range(5): # waste the crap at the top of the output
			outcome.readline()
		total = 0
		cnt = 0
		for line in outcome:
			if "+" not in line and "4" in line:
				cnt += 1
			total += 1


		print cnt/total
