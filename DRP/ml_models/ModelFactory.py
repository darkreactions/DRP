import importlib

# TODO: Mebbe change dis to ModelContainor
class ModelFactory():
  """
  The ModelFactory is responsible for constructing and preparing
  new machine-learning models. It does not actually test the model;
  instead, it prepares the "recipe" for the model and saves it.

  modelLibrary === The library in which the modelType will be found.
  modelType === The type of Machine-Learning model that will be created.

  """

  # TODO: `build` should accept training_reactions_reactions and test_reactions and *not* a splitter.
  def build(self, training_reactions, test_reactions, predictors, responses, modelLibrary="weka", modelType="svm", debug=True):
    """Constructs, trains, and then tests a ML-model using a ModelVisitor
       of the given modelLibrary, modelType, and splitterType."""

    model = self._getModelVisitor(modelLibrary, modelType)

    if debug:
      model.enableDebug() #TODO:  Turn this into a list of strings/exceptions and attach to self.

    if debug:
      print "\nStarting model generation in debug mode..."

    model.setPredictors(predictors) #TODO:Assume predictors/resposnes are QSs. :)
    model.setResponses(responses)

    model.settraining_reactionsData(training_reactions)
    model.settest_reactionsData(test_reactions)

    if debug:
      print "\ntraining_reactions the model on {} entries...".format(training_reactions.count())
    model._train()

    if debug:
      print "\ntest_reactions the model on {} entries...".format(test_reactions.count())
    model._test()

    return model

  def _getModelVisitor(self, library, tool):
    """Returns a ModelVisitor object associated with this library and tool."""
    name = "{}_{}".format(library, tool)
    try:
      mod = importlib.import_module("DRP.ml_models.model_visitors.{}".format(name))
      return mod.ModelVisitor()
    except ImportError:
      error = "Model Visitor \"{}\" is not supported by the ModelFactory.".format(name)
      raise NotImplementedError(error)
