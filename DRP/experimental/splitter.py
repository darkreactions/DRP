def split_data(data, headers, response, split=0.5):
  from sklearn.cross_validation import train_test_split

  header_index = headers.index(response)
  X = [[elem for i, elem in enumerate(row) if i!=header_index] for row in data]
  y = [row[header_index] for row in data]


  # TODO: Want to split data so equal representations in train+test
  return splits

