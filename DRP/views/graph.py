import os, sys
from django.contrib.auth.decorators import login_required

full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


def get_graph(lines, base, xLabel="Percentage"):
  def frange(step):
    bottom = 0
    top = 1

    while bottom<top:
      yield bottom
      bottom += step

    yield top

  import matplotlib.pyplot as plt
  plt.ioff()

  figure = plt.figure(figsize=(20,10))
  figure.patch.set_alpha(0)

  ax = figure.add_subplot(1,1,1)

  for header, y in lines.items():
    plt.plot(base, y, label=header)

    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
            fancybox=True, shadow=True, ncol=10)

  plt.xlabel(xLabel.capitalize())

  ax.set_yticks(list(frange(0.1)))
  ax.set_yticks(list(frange(0.025)), minor=True)
  ax.grid(which='minor', alpha=0.2)
  ax.grid(which='major', alpha=0.5)

  return figure


@login_required
def graph(request, base="test"):
  import matplotlib
  matplotlib.use("Agg")
  matplotlib.rc("font", size=14)

  from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

  from django.http import HttpResponse, HttpResponseNotFound
  from DRP.retrievalFunctions import get_usable_models
  from DRP.filters import apply_filters

  models = get_usable_models()
  models = apply_filters(request, models, model="ModelStats")

  if models.count()==0:
    return HttpResponseNotFound()

  classes = "2"
  to_ignore = {
    "% TP",
    "% TN",
    "% FP",
    "% FN",
    "User Satisfaction",
  }

  stats = []
  baseline = []
  for model in models:
    try:
      this_stats = model.stats()[classes]

      if base=="train":
        x_val = model.train_size
      elif base=="test":
        x_val = model.total()

      elif base=="time":
        import datetime
        epoch = datetime.datetime.utcfromtimestamp(0)
        x_val = (model.end_time - epoch).total_seconds()
      else:
        return HttpResponseNotFound()

      del this_stats["Test Size"]
      del this_stats["Train Size"]

      stats.append( this_stats )
      baseline.append( x_val )
    except:
      pass

  if base=="time":
    sorted_baseline = sorted(baseline)
    baseline = [sorted_baseline.index(val) for val in baseline]


  fields = stats[0].keys()

  lines = {field:[stat[field] for stat in stats] for field in fields}
  lines = {field:vals for field,vals in lines.items() if field not in to_ignore}

  fig = get_graph(lines, baseline, xLabel=base)

  response= HttpResponse(mimetype='image/png')

  # Be wary, the 'runserver' uses threading by default, which `print_png`
  #  seems to dislike. I'd recommend using the "--nothreading" option
  #  in the runserver.
  canvas = FigureCanvas(fig)
  canvas.print_png(response)

  return response

