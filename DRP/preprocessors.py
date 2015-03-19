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
  data = [d.get_calculations_list() for d in data]

  # Set outcomes of 0 to be 1.
  outcome_index = headers.index("outcome")
  outcomes = [int(row[outcome_index]) for row in data]
  outcomes = [elem if elem>0 else 1 for elem in outcomes]
  for i, outcome in enumerate(outcomes):
    data[i][outcome_index] = str(outcome)


  if normalize:

    # Normalize everything except the categorized columns.
    categories = [
        "outcome", "purity"
    ]
    cat_indexes = {headers.index(cat) for cat in categories if cat in headers}

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
  data = [d.get_calculations_list() for d in data]

  data = filter(remover, data)

  return [headers] + data


def category_preprocessor(data):
  def categorize(data):
    return "Se" in data.atoms and "V" in data.atoms

  new_data = default_preprocessor(data)

  headers = new_data.pop(0)
  new_data = [row+[categorize(datum)] for row, datum in zip(new_data, data)]

  return [headers] + new_data

