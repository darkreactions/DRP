import importlib

class ModelFactory():
  """
  The ModelFactory is responsible for constructing and preparing
  new machine-learning models. It does not actually test the model;
  instead, it prepares the "recipe" for the model and saves it.

  modelLibrary === The library in which the modelType will be found.
  modelType === The type of Machine-Learning model that will be created.

  """

  def build(self, reactions, headers, modelLibrary="weka", modelType="svm",
                  splitterType="MutualInfoSplitter", debug=True):

    model = self._getModelVisitor(modelLibrary, modelType)

    if debug:
      model.enableDebug()

    splitter = self._getSplitter(splitterType)
    model.setSplitter(splitter)

    training, testing = splitter.split(reactions)
    model.setTrainingData(training)
    model.setTestingData(testing)

    if debug:
      print "Training the model on {} entries...".format(training.count())
    model._train()

    if debug:
      print "Testing the model on {} entries...".format(testing.count())
    model._test()

    return model

  def _getModelVisitor(self, library, tool):
    name = "{}_{}".format(library, tool)
    try:
      mod = importlib.import_module("DRP.ml_models.model_visitors.{}".format(name))
      return mod.ModelVisitor()
    except Exception as e:
      print e
      error = "Model Visitor \"{}\" is not supported by the ModelFactory.".format(name)
      raise NotImplementedError(error)


  def _getSplitter(self, name):
    try:
      mod = importlib.import_module("DRP.ml_models.splitters.{}".format(name))
      return mod.Splitter()
    except:
      error = "Splitter \"{}\" is not supported by the ModelFactory.".format(name)
      raise NotImplementedError(error)
