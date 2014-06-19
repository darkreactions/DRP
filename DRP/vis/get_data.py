mport csv
import time

"""
This library contains some functions you may find useful when working with the
given data.  Feel free to modify this to be more general and/or to point to 
other data you're looking at.  Or feel free not to use these at all if you don't
 find them useful.
"""

# reads in all the fields of the file as values (in order), 
# with their keys being the header of their respective columns
def get_data_list_of_dicts(data_file):
    list = []
    with open(data_file) as f:
        f_csv = csv.DictReader(f)
	for row in f_csv:
	    list.append(row)
	    if len(list) > 100000:
	    	break	
    return list

def get_headers(data_file):
    with open(data_file) as f:
        f_csv = csv.reader(f)
		# iterates over first row of file, appening each cell (column head) into list of headers
        headers = next(f_csv)
    return headers

def write_data_dicts(filename, headers, rows_list_of_dicts):
    with open(filename,'w') as f:
        f_csv = csv.DictWriter(f, headers)
        f_csv.writeheader()
        f_csv.writerows(rows_list_of_dicts)

def write_data(filename, headers, rows_list_of_lists):
    with open(filename,'w') as f:
        f_csv = csv.writer(f)
        f_csv.writerow(headers)
        f_csv.writerows(rows_list_of_lists)

def get_data_slice(column_name, list_of_dicts):
    list = []
    for dict in list_of_dicts:
        list.append(dict[column_name])
    return list

def date_convert(date_string):
    date = time.strptime(date_string, "%d/%m/%Y %H:%M:%S") 
    return time.mktime(date)

if __name__ == "__main__": _test()	
