from django.db import models
from rxnDescriptors import BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor
from rxnDescriptorValues import BoolRxnDescriptorValue, OrdRxnDescriptorValue, NumRxnDescriptorValue, CatRxnDescriptorValue
from StatsModel import StatsModel

class PredictedDescriptor(models.Model):
    stats_model = models.ForeignKey(StatsModel)

    class Meta:
        app_label="DRP"
        abstract = True

    def getPredictionTuples(self, model):
        """Returns a dictionary of lists of outcome tuples, where the keys are the
           outcomeDescriptors and the outcomes are in the format (correct, guess).

           IE: [(true,guess),(true',guess'),(true'',guess'')]
           EG: [(1,2),(1,1),(2,2),(3,2),(4,3)]"""


        if isinstance(self, BoolRxnDescriptor):
            valueType = BoolRxnDescriptorValue
        elif isinstance(self, OrdRxnDescriptor):
            valueType = OrdRxnDescriptorValue
        elif isinstance(self, CatRxnDescriptor):
            valueType = CatRxnDescriptorValue
        elif isinstance(self, NumRxnDescriptor):
            valueType = NumRxnDescriptorValue
        else:
            error = "Unknown descriptorValue type for '{}'".format(self)
            raise NotImplementedError(error)

        predictions = []
        for prediction in valueType.objects.filter(model=model, descriptor=self):
            true = valueType.objects.get(reaction=prediction.reaction,
                                         descriptor=self.prediction_of).value
            guess = prediction.value
            predictions.append( (true, guess) )

        return predictions

class PredBoolRxnDescriptor(BoolRxnDescriptor, PredictedDescriptor):
    prediction_of = models.ForeignKey(BoolRxnDescriptor, related_name="prediction_of")

    class Meta:
        app_label='DRP'
        verbose_name = 'Predicted Boolean Reaction Descriptor'


class PredOrdRxnDescriptor(OrdRxnDescriptor, PredictedDescriptor):
    prediction_of = models.ForeignKey(OrdRxnDescriptor, related_name="predition_of")

    class Meta:
        app_label='DRP'
        verbose_name = 'Predicted Ordinal Reaction Descriptor'

    def summarize(self, model):
        return "Four-Class Accuracy: {}".format(self.fourClassAccuracy(model))

    def fourClassAccuracy(self, model):
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
    prediction_of = models.ForeignKey(NumRxnDescriptor, related_name="prediction_of")

    class Meta:
        app_label='DRP'
        verbose_name = 'Predicted Numeric Reaction Descriptor'


class PredCatRxnDescriptor(CatRxnDescriptor, PredictedDescriptor):
    prediction_of = models.ForeignKey(CatRxnDescriptor, related_name="prediction_of")

    class Meta:
        app_label='DRP'
        verbose_name = 'Predicted Categorical Reaction Descriptor'

