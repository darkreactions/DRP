import importlib

class ModelFactory():
  """
  The ModelFactory is responsible for constructing and preparing
  new machine-learning models. It does not actually test the model;
  instead, it prepares the "recipe" for the model and saves it.

  modelLibrary === The library in which the modelType will be found.
  modelType === The type of Machine-Learning model that will be created.

  """

  def build(self, reactions, predictors, responses, modelLibrary="weka", modelType="svm",
                  splitterType="MutualInfoSplitter", debug=True):

    model = self._getModelVisitor(modelLibrary, modelType)

    if debug:
      model.enableDebug()

    if debug:
      print "\nStarting model generation in debug mode..."

    splitter = self._getSplitter(splitterType)
    model.setSplitter(splitter)

    model.setPredictors(predictors)
    model.setResponses(responses)

    training, testing = splitter.split(reactions)
    model.setTrainingData(training)
    model.setTestingData(testing)

    if debug:
      print "\nTraining the model on {} entries...".format(training.count())
    model._train()

    if debug:
      print "\nTesting the model on {} entries...".format(testing.count())
    model._test()

    return model

  def _getModelVisitor(self, library, tool):
    """ Returns a ModelVisitor object associated with this library and tool."""
    name = "{}_{}".format(library, tool)
    try:
      mod = importlib.import_module("DRP.ml_models.model_visitors.{}".format(name))
      return mod.ModelVisitor()
    except ImportError:
      error = "Model Visitor \"{}\" is not supported by the ModelFactory.".format(name)
      raise NotImplementedError(error)


  def _getSplitter(self, name):
    """ Returns the splitter object associated with this splitter namename."""
    try:
      mod = importlib.import_module("DRP.ml_models.splitters.{}".format(name))
      return mod.Splitter()
    except ImportError:
      error = "Splitter \"{}\" is not supported by the ModelFactory.".format(name)
      raise NotImplementedError(error)
