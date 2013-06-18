function quick_validation(field, dirty_datum) {
	#Import the ranges that will be used to determine data legality. 
	try {
		$.getJSON('{{ STATIC_URL }}data_ranges.json', function(data) {
		  alert(data);
		});
	}
	catch (error) {
		alert ("ERROR LOADING DATA_RANGES");
	}
		
	//try {
		//if (["mass_1","mass_2", "mass_3", "mass_4", "mass_5"].indexOf(field) != -1) {
			//dirty_datum = float(dirty_datum)
			//assert(len(str(dirty_datum))<10) #Prevent super long decimals.
			//assert type(dirty_datum) == type(1.0)
			//assert mass_range[0] < dirty_datum < mass_range[1]
		//} else if (field=="temp" {
			//dirty_datum = float(dirty_datum)
			//assert(temp_range[0] < dirty_datum < temp_range[1])
		//} else if (field=="pH" {
			//dirty_datum = float(dirty_datum)
			//assert(pH_range[0] < dirty_datum < pH_range[1])
		//} else if (field=="outcome" {
			//dirty_datum = int(dirty_datum)
			//assert(outcome_range[0] <= dirty_datum <= outcome_range[1])
		//} else if (field=="purity" {
			//dirty_datum = int(dirty_datum)
			//assert(purity_range[0] <= dirty_datum <= purity_range[1])
		//} else if (field=="time" {
			//dirty_datum = int(dirty_datum)
			//assert(time_range[0] <= dirty_datum <= time_range[1])
		//} else if (field=="slow_cool" or field=="leak" {
			//dirty_datum = str(dirty_datum).lower()
			//if "yes" == dirty_datum or "true" == dirty_datum { ###convert "y" to "yes" in jQuery
				//return True
			//} else if ("no" == dirty_datum or "false" == dirty_datum {
				//return True
			//else {
				//raise Exception
		//} else if (field=="ref" {
			//dirty_datum = str(dirty_datum)
			//assert(ref_range[0] < len(dirty_datum) <= ref_range[1])
		//return True #Dirty datum passes the quick verification test. 
	//except {
		//return False
