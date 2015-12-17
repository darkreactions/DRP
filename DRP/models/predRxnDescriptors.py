from django.db import models
from rxnDescriptors import BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor
from ModelContainer import ModelContainer
from StatsModel import StatsModel

class PredictedDescriptor(models.Model):
    modelContainer = models.ForeignKey(ModelContainer)
    statsModel = models.ForeignKey(StatsModel, null=True)

    class Meta:
        app_label="DRP"
        abstract = True


class PredBoolRxnDescriptor(BoolRxnDescriptor, PredictedDescriptor):
    predictionOf = models.ForeignKey(BoolRxnDescriptor, related_name="prediction_of")

    class Meta:
        app_label='DRP'
        verbose_name = 'Predicted Boolean Rxn Descriptor'


class PredOrdRxnDescriptor(OrdRxnDescriptor, PredictedDescriptor):
    predictionOf = models.ForeignKey(OrdRxnDescriptor, related_name="predition_of")

    class Meta:
        app_label='DRP'
        verbose_name = 'Predicted Ordinal Rxn Descriptor'

    def summarize(self, model):
        return "Accuracy: {}".format(self.accuracy(model))

    def accuracy(self, model):
      conf = self.getConfusionMatrix(model)

      correct = 0.0
      total = 0.0
      for true, guesses in conf.items():
        for guess, count in guesses.items():
          if true == guess: correct += count
          total += count
      return correct/total

    def getConfusionMatrix(self, model):
        """Returns a dicionary of dictionaries of dictionaries, where the outer keys
           are the outcomeDescriptors, the middle keys are the "correct" or "true"
           values, the innermost keys are the "guessed" values that occurred, and
           the value is the integer number of occurrences of that guess when the
           true descriptor was the middle key.

           IE: {true: {guess:#, guess':#},
                true': {guess:#, guess':#}}
           Eg: {"1": {"1": 10
                      "2": 10
                      "3": 13
                      "4": 0
                     }
               , ...
               }
              } """
        matrix = {}
        for true, guess in self.getPredictionTuples(model):
            if true not in matrix:
                matrix[true] = {}

            matrix[true][guess] = matrix[true][guess]+1 if guess in matrix[true] else 1

        return matrix


class PredNumRxnDescriptor(NumRxnDescriptor, PredictedDescriptor):
    predictionOf = models.ForeignKey(NumRxnDescriptor, related_name="prediction_of")

    class Meta:
        app_label='DRP'
        verbose_name = 'Predicted Numeric Rxn Descriptor'


class PredCatRxnDescriptor(CatRxnDescriptor, PredictedDescriptor):
    predictionOf = models.ForeignKey(CatRxnDescriptor, related_name="prediction_of")

    class Meta:
        app_label='DRP'
        verbose_name = 'Predicted Categorical Rxn Descriptor'

