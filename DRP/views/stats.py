

from django.http import HttpResponse
from django.shortcuts import render
from django.contrib.auth.decorators import login_required


@login_required
def model_stats(request, model_id):

  def un_unicode(raw):
    if type(raw)==dict:
      return {un_unicode(key):un_unicode(val) for key, val in raw.items()}
    elif type(raw)==list:
      return [un_unicode(elem) for elem in raw]
    else:
      return str(raw)

  from DRP.models import ModelStats

  try:
    model = ModelStats.objects.get(id=model_id)

    train_cm = un_unicode(model.load_confusion_table(table="train"))
    test_cm = un_unicode(model.load_confusion_table(table="test"))

    return render(request, "model_stats.html", {
                              "model":model,
                              "train_cm":train_cm,
                              "test_cm":test_cm,
                                                })

  except Exception as e:
    print e
    return HttpResponse("Failed to summarize model...")


