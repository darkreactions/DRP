
def get_avg_confusion_dict(model_stats, table="test"):
  """
  Returns an average confusion dict from a list of model_stats
  """

  conf_dicts = [m.load_confusion_dict(table=table) for m in model_stats]
  return _avg_confusion_dicts(conf_dicts)

def _avg_confusion_dicts(conf_dicts):
  """
  Returns an average confusion dict from a list of confusion dicts.
  """

  avg_dict = {}

  for conf_dict in conf_dicts:
    for guess, actuals in conf_dict.items():

      if not guess in avg_dict:
        avg_dict[guess] = {}

      for actual, occurrences in actuals.items():
        if not actual in avg_dict[guess]:
          avg_dict[guess][actual] = 0

        avg_dict[guess][actual] += occurrences

  num_models = float(len(conf_dicts))
  trunc_divide = lambda x,y: float("{:.3f}".format(x/y))

  avg_dict ={guess:{actual:trunc_divide(count,num_models)
                      for actual,count in actuals.items()}
                        for guess, actuals in avg_dict.items()}

  return avg_dict


def make_confusion_dict(guesses, actuals):
  if len(guesses)==0 or len(actuals)==0:
    raise Exception("Either `guesses` or `actuals` is empty!")

  if not (len(guesses)==len(actuals)):
    raise Exception("`guesses` and `actuals` are of different sizes!")

  possible_vals = sorted(list(set(guesses+actuals)))
  cm = {r1:{r2:0 for r2 in possible_vals} for r1 in possible_vals}

  for guess, actual in zip(guesses, actuals):
    cm[actual][guess] += 1

  return cm

