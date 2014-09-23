
def build_previous_model(model_name, model_description, date, data=None):
  """
  Constructs a model from the data available on a given date.
  """

  from DRP.model_building import model_methods
  from DRP.retrievalFunctions import filter_by_date
  from DRP.models import Data

  if not data:
    data = Data.objects.all()

  filtered = filter_by_date(data, date, "previous")
  model_methods.gen_model(model_name, model_description, filtered)


