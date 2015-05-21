


def get_matrix(filename):
  import csv

  with open(filename) as f:
    reader = csv.reader(f)
    data = [row for row in reader]
    headers = data.pop(0)
    return headers, data


def write_csv(new_file, data):
  import csv

  with open(new_file, "w") as f:
    writer = csv.writer(f)

    for row in data:
      writer.writerow(row)

  print "Written to {}".format(new_file)


def remove_duplicates(data):
  seen = set()
  ref_index = 1
  unique = []

  for row in data:
    if row[ref_index] not in seen:
      seen.add( row[ref_index] )

      unique.append(row)

  return unique


def main():
  import sys

  # Variable Setup
  sys.argv.pop(0) # Remove this script from the arguments.
  end_csv = sys.argv.pop()
  total_data = []
  total_headers = []

  for filename in sys.argv:
    print "Merging {}...".format(filename)

    headers, data = get_matrix(filename)

    if not total_headers:
      total_headers.extend(headers)

    total_data.extend(data)

  unique_data = remove_duplicates(total_data)

  end_data = [total_headers] + unique_data
  write_csv(end_csv, end_data)



if __name__=="__main__":
  main()
