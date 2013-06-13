#quick_validation should be asserted to check for validity
#False returns indicate invalid data; True returns are valid data.
def quick_validation(field, dirty_datum):
	#Import the ranges that will be used to determine data legality. 
	try:
		from data_ranges import *
	except:
		print "ERROR LOADING DATA_RANGES"
		
	try:
		if field in {"quantity_1","quantity_2", "quantity_3", "quantity_4", "quantity_5"}:
			dirty_datum = float(dirty_datum)
			assert(len(str(dirty_datum))<10) #Prevent super long decimals.
			assert type(dirty_datum) == type(1.0)
			assert quantity_range[0] < dirty_datum < quantity_range[1]
		if field in {"unit_1","unit_2", "unit_3", "unit_4", "unit_5"}:
			dirty_datum = str(dirty_datum).lower()
			assert dirty_datum in {"ml", "g","d"} ###AUTO-Gen?
		elif field=="temp":
			dirty_datum = float(dirty_datum)
			assert(temp_range[0] < dirty_datum <= temp_range[1]) #No absolute zero.
		elif field=="pH":
			dirty_datum = float(dirty_datum)
			assert(pH_range[0] <= dirty_datum <= pH_range[1])
		elif field=="outcome":
			dirty_datum = int(dirty_datum)
			assert(outcome_range[0] <= dirty_datum <= outcome_range[1])
		elif field=="purity":
			dirty_datum = int(dirty_datum)
			assert(purity_range[0] <= dirty_datum <= purity_range[1])
		elif field=="time":
			dirty_datum = int(dirty_datum)
			assert(time_range[0] <= dirty_datum <= time_range[1])
		elif field=="slow_cool" or field=="leak":
			dirty_datum = str(dirty_datum).lower()
			assert(dirty_datum in {"yes","no","?"})###AUTO-Gen?
		elif field=="ref":
			dirty_datum = str(dirty_datum)
			assert(ref_range[0] < len(dirty_datum) <= ref_range[1])
		return True #Dirty datum passes the quick verification test. 
	except:
		print "FAILED SERVER VALIDATION"###
		return False
	
