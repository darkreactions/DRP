

# Set the Python path so that it has access to the Django settings.
import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


from DRP.models import ModelStats
from DRP.model_building.confusion_table import make_confusion_dict, _avg_confusion_dicts


def test1():
  m = ModelStats()

  m.set_correct_vals(["1"])
  guesses = [1,1,1,1,1,0,0,0,0,0]
  actual = [1,1,1,1,1,0,0,0,0,0]
  cm =  make_confusion_dict(guesses, actual)
  m.set_confusion_table(cm)

  result = (m.accuracy()==1 and m.precision()==1 and
           m.total("test")==len(guesses))

  m.save()
  m.delete()

  return result

def test2():
  m = ModelStats()

  m.set_correct_vals(["1"])
  guesses = [1,1,1,1,1,0,0,0,0,0]
  actual = [0,0,0,0,0,0,0,0,0,0]
  cm =  make_confusion_dict(guesses, actual)
  m.set_confusion_table(cm)

  result = (m.accuracy()==0.5 and m.precision()==0 and
           m.false_positives(normalize=True)==0.5 and
           m.total("test")==len(guesses))

  m.save()
  m.delete()

  return result

def test3():
  m = ModelStats()

  m.set_correct_vals(["1"])
  guesses = [1,1,1,1,1,0,0,0,0,0]
  actual = [1,1,1,1,1,1,1,1,1,1]
  cm =  make_confusion_dict(guesses, actual)
  m.set_confusion_table(cm)

  result = (m.accuracy()==0.5 and m.precision()==1 and
           m.recall()==0.5 and
           m.false_positives(normalize=True)==0 and
           m.true_positives(normalize=True)==0.5 and
           m.false_negatives(normalize=True)==0.5 and
           m.total("test")==len(guesses))

  m.save()
  m.delete()

  return result

def test4():
  m = ModelStats()

  m.set_correct_vals(["1", "2"])
  guesses = [2,2,2,2,2,0,0,0,0,0]
  actual = [1,1,1,1,1,1,1,1,1,1]
  cm =  make_confusion_dict(guesses, actual)
  m.set_confusion_table(cm)

  result = (m.accuracy()==0.5 and m.precision()==1 and
           m.recall()==0.5 and
           m.false_positives(normalize=True)==0 and
           m.true_positives(normalize=True)==0.5 and
           m.true_positives(ranges=False, normalize=True)==0 and
           m.accuracy(ranges=False)==0 and
           m.false_negatives(normalize=True)==0.5 and
           m.total("test")==len(guesses))

  m.save()
  m.delete()

  return result

def test5():
  guesses = [1,1,1,1,1]
  actual = [1,1,1,1,1]
  cm1 =  make_confusion_dict(guesses, actual)

  guesses = [0,0,0,0,0]
  actual = [1,1,1,1,1]
  cm2 =  make_confusion_dict(guesses, actual)

  avg = _avg_confusion_dicts([cm1, cm2])

  m = ModelStats()
  m.set_correct_vals(["1"])
  m.set_confusion_table(avg)

  print "avg", avg
  print "TP", m.true_positives()
  print "FP", m.false_positives()
  print "acc", m.accuracy()
  print "prec", m.precision()
  result = m.accuracy()==0.5 and m.recall()==0.5 and m.precision()==1

  m.save()
  m.delete()

  return result




def main():
  tests = [test1, test2, test3, test4, test5]
  successes = 0

  for test in tests:
    if test():
      successes += 1
    else:
      print "-- '{}' failed!".format(test.__name__)

  print "{}/{} Passed!".format(successes, len(tests))


if __name__=="__main__":
  main()

