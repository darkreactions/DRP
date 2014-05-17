import itertools, sys, os

def command_help():
	print "Command syntax: python model_performance.py <prefix>"
	
if __name__ == "__main__":
	if len(sys.argv) != 2:
		command_help()
	else: 
	
		#Assume that script is in scripts folder.
		CHEMML_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
		prefix = "{0}/{1}/{1}".format(CHEMML_DIR, sys.argv[1])
	
		total = 0
		out_good = 0
		pur_good = 0
		both_good = 0
		with open("%s_out.out" % prefix) as out, open("%s_pur.out" % prefix) as pur:
			for out_line, pur_line in itertools.izip(out,pur):
				out_notwrong = not "+" in out_line
				pur_notwrong = not "+" in pur_line
				total += 1
				if out_notwrong:
					out_good += 1
				if pur_notwrong:
					pur_good += 1
				if out_notwrong and pur_notwrong:
					both_good += 1
		
		perc_out_good = out_good/float(total)
		perc_pur_good = pur_good/float(total)
		perc_both_good = both_good/float(total)
		
		print "-- \033[4m Model performance: \033[0m" #Underline for style-points...
		print "---- Overall:	{}".format(total)
		print "---- Out Good:	{0}	{1:.4%}".format(out_good, perc_out_good)
		print "---- Pur Good:	{0}	{1:.4%}".format(pur_good, perc_pur_good)
		print "---- Both Good:	{0}	{1:.4%}".format(both_good, perc_both_good)
