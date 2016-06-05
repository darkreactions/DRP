import os
import sys
from django.contrib.auth.decorators import login_required

full_path = os.path.dirname(os.path.realpath(__file__)) + "/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
    sys.path = [django_path] + sys.path
    os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


@login_required
def graph(request, base="test"):
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rc("font", size=14)

    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

    from django.http import HttpResponse, HttpResponseNotFound
    from DRP.retrievalFunctions import get_usable_models
    from DRP.filters import apply_filters
    from DRP.graph import get_graph

    models = get_usable_models()
    models = apply_filters(request, models, model="ModelStats")

    if models.count() == 0:
        return HttpResponseNotFound()

    category = "2-test"
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
            this_stats = model.stats()[category]

            if base == "train" or base == "test":
                x_val = model.total(table=base)

            elif base == "time":
                import datetime
                epoch = datetime.datetime.utcfromtimestamp(0)
                x_val = (model.end_time - epoch).total_seconds()
            else:
                return HttpResponseNotFound()

            del this_stats["Test Size"]
            del this_stats["Train Size"]

            stats.append(this_stats)
            baseline.append(x_val)
        except:
            pass

    if base == "time":
        sorted_baseline = sorted(baseline)
        baseline = [sorted_baseline.index(val) for val in baseline]

    fields = stats[0].keys()

    lines = {field: [stat[field] for stat in stats] for field in fields}
    lines = {field: vals for field, vals in lines.items() if field not in to_ignore}

    fig = get_graph(lines, baseline, xLabel=base)

    response = HttpResponse(mimetype='image/png')

    # Be wary, the 'runserver' uses threading by default, which `print_png`
    #  seems to dislike. I'd recommend using the "--nothreading" option
    #  in the runserver.
    canvas = FigureCanvas(fig)
    canvas.print_png(response)

    return response
