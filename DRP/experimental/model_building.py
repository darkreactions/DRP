
def build_previous_model(model_name, model_description, date):
  """
  Constructs a model from the data available on a given date.
  """

  from DRP.model_building import model_methods
  from DRP.retrievalFunctions import filter_by_date
  from DRP.models import Data

  filtered = filter_by_date(Data.objects.all(), date, "previous")
  model_methods.gen_model(model_name, model_description, filtered)

