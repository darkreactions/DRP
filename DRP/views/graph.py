import os, sys

full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


def get_graph(lines, base, xLabel="Percentage"):
  import matplotlib.pyplot as plt

  for header, y in lines.items():
    plt.plot(base, y, "o", label=header)
    plt.plot(base, y, label=header)
    plt.xlabel(xLabel)
    plt.ylabel(header)

  return plt.figure()


def graph(request):
  from django.http import HttpResponse
  from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

  from DRP.retrievalFunctions import get_usable_models
  from DRP.filters import apply_filters

  models = get_usable_models()
  models = apply_filters(request, models, model="ModelStats")

  classes = "2"

  stats = []
  totals = []
  for model in models:
    try:
      stats.append( model.stats()[classes] )
      totals.append( model.total() )
    except:
      pass

  fields = stats[0].keys()

  lines = {field:[stat[field] for stat in stats] for field in fields}

  fig = get_graph(lines, totals)

  response= HttpResponse(mimetype='image/png')

  # Be wary, the 'runserver' uses threading by default, which `print_png`
  #  seems to dislike. I'd recommend using the "--nothreading" option
  #  in the runserver.
  canvas = FigureCanvas(fig)
  canvas.print_png(response)

  return response

