def default_preprocessor(data):
  def is_num(elem):
    try:
      float(elem)
      return True
    except:
      return False

  def norm(elems):
    import operator
    elems = map(float,elems)
    avg = float(reduce(operator.add, elems))/len(elems)

    return [elem/avg for elem in elems]

  # Variable Setup
  normalize = False

  headers = data.pop(0)

  # Set outcomes of 0 to be 1.
  outcome_index = headers.index("outcome")
  outcomes = [int(row[outcome_index]) for row in data]
  outcomes = [elem if elem>0 else 1 for elem in outcomes]
  for i, outcome in enumerate(outcomes):
    data[i][outcome_index] = str(outcome)

  # Normalize everything except the categorized columns.
  categories = [
    "outcome", "purity"
  ]

  cat_indexes = {headers.index(cat) for cat in categories if cat in headers}

  if normalize:
    cols = [[unicode(row[i]) for row in data] for i, h in enumerate(headers)]
    cols = [norm(col) if (all(map(is_num, col)) and i not in cat_indexes) else col
                    for i, col in enumerate(cols)]
    data = [[col[i] for col in cols] for i in xrange(len(data))]

  return [headers]+data


# TODO: Un-tested and perhaps unnecessary...
def purging_preprocessor(data):
  def remover(row):
    to_remove = {
      "water"
    }
    return not any(unicode(elem).lower() in to_remove for elem in row)

  headers = data.pop(0)

  data = filter(remover, data)

  return [headers] + data




