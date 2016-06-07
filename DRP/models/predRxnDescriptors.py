"""Module for information about predicted reaction descriptors."""
from django.db import models
from rxnDescriptors import BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor
from rxnDescriptorValues import BoolRxnDescriptorValue
from descriptors import DescriptorManager
from ModelContainer import ModelContainer
from StatsModel import StatsModel
from PerformedReaction import PerformedReaction
from django.db.models import F


class PredictedDescriptor(models.Model):

    """The general case of a predicted descriptor."""
    modelContainer = models.ForeignKey(ModelContainer)
    statsModel = models.ForeignKey(StatsModel, null=True)

    class Meta:
        app_label = "DRP"
        abstract = True

    objects = DescriptorManager()


class PredBoolRxnDescriptor(BoolRxnDescriptor, PredictedDescriptor):

    """Reaction boolean descriptor which has been predicted by a model."""

    predictionOf = models.ForeignKey(BoolRxnDescriptor, related_name="prediction_of")

    class Meta:
        app_label = 'DRP'
        verbose_name = 'Predicted Boolean Rxn Descriptor'

    objects = DescriptorManager()

    def summarize(self):
        """Return the accuracy of the boolean predictions."""
        return "Accuracy: {}".format(self.accuracy())

    def accuracy(self):
        """Calculate the accuracy of the boolean predictions."""
        conf = self.getConfusionMatrix()

        correct = 0.0
        total = 0.0
        for true, guesses in conf.items():
            for guess, count in guesses.items():
                if true == guess:
                    correct += count
                total += count
        return correct / total

    def oldGetConfusionMatrix(self):
        """
        The old method of returning confusion matrices.

        Returns a dicionary of dictionaries, where the outer keys are the "correct" or "true"
       values, the inner keys are the "guessed" values that occurred, and
       the value is the integer number of occurrences of that guess when the
       true descriptor was the second key.

       Eg: {true: {guess:#, guess':#},
            true': {guess:#, guess':#}}
       Eg: {"1": {"1": 10
                  "2": 10
                  "3": 13
                  "4": 0
                 }
           , ...
           }
          }
        """
        matrix = {
            True: {True: 0, False: 0},
            False: {True: 0, False: 0}
        }
        for true, guess in self.getPredictionTuples():
            matrix[true][guess] += 1

        return matrix

    def getConfusionMatrix(self, reactions=None):
        """
        Return the confusion matrices.

        Return a dicionary of dictionaries, where the outer keys are the "correct" or "true"
        values, the inner keys are the "guessed" values that occurred, and
        the value is the integer number of occurrences of that guess when the
        true descriptor was the second key.

        Eg: {true: {guess:#, guess':#},
            true': {guess:#, guess':#}}
        Eg: {"1": {"1": 10
                  "2": 10
                  "3": 13
                  "4": 0
                 }
           , ...
           }
          }
        """
        if reactions is None:
            reactions = PerformedReaction.objects.filter(boolrxndescriptorvalue__descriptor=self).distinct()

        permittedValues = [True, False]

        reactions = reactions.filter(boolrxndescriptorvalue__descriptor=self).annotate(predicted_val=F('boolrxndescriptorvalue__value'))
        reactions = reactions.filter(boolrxndescriptorvalue__descriptor=self.predictionOf).annotate(actual_val=F('boolrxndescriptorvalue__value'))

        matrix = {true: {guess: reactions.filter(predicted_val=guess, actual_val=true).count() for guess in permittedValues} for true in permittedValues}

        return matrix

    def getPredictionTuples(self):
        """
        Return the prediction tuples.       

        Return a list of tuples where the first value is the actual value for
        a descriptor of a reaction and the second value is the predicted value
        of that descriptor in the same reaction.
        EG: [(True,True), (False,True), (False,True), (True,True), (True,True)] for a model that always
            predicts "True"
        """

        actualDescValues = self.predictionOf.boolrxndescriptorvalue_set.all()
        predictedDescValues = self.boolrxndescriptorvalue_set.all()

        predictionTuples = []
        for prediction in predictedDescValues:
            actual = actualDescValues.get(reaction=prediction.reaction)
            predictionTuples.append((actual.value, prediction.value))

        return predictionTuples


class PredOrdRxnDescriptor(OrdRxnDescriptor, PredictedDescriptor):

    """Ordinal case of the predicted reaction descriptor."""

    predictionOf = models.ForeignKey(OrdRxnDescriptor, related_name="predition_of")

    class Meta:
        app_label = 'DRP'
        verbose_name = 'Predicted Ordinal Rxn Descriptor'

    objects = DescriptorManager()

    def summarize(self):
        """Return the accuracy of the predicted reaction descriptor."""
        return "Accuracy: {}".format(self.accuracy())

    def accuracy(self):
        """Calculate the accurate of the predicted reaction descriptor."""
        conf = self.getConfusionMatrix()

        correct = 0.0
        total = 0.0
        for true, guesses in conf.items():
            for guess, count in guesses.items():
                if true == guess:
                    correct += count
                total += count
        return correct / total

    def getConfusionMatrix(self):
        """
        Return a confusion matrix.

        Return a dicionary of dictionaries of dictionaries, where the outer keys
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

        matrix = {
            true: {guess: 0 for guess in xrange(self.minimum, self.maximum + 1)}
            for true in xrange(self.minimum, self.maximum + 1)
        }
        for true, guess in self.getPredictionTuples():
            matrix[true][guess] += 1
        return matrix

    def getPredictionTuples(self):
        """
        Return the prediction tuples.
 
        Return a list of tuples where the first value is the actual value for
        a descriptor of a reaction and the second value is the predicted value
        of that descriptor in the same reaction.
        EG: [(1,1), (2,1), (4,1), (3,1), (1,1)] for a model that always
            predicts "1" if there are 4 different values for a descriptor.
        """

        actualDescValues = self.predictionOf.ordrxndescriptorvalue_set.all()
        predictedDescValues = self.ordrxndescriptorvalue_set.all()

        predictionTuples = []
        for prediction in predictedDescValues:
            actual = actualDescValues.get(reaction=prediction.reaction)
            predictionTuples.append((actual.value, prediction.value))

        return predictionTuples


class PredNumRxnDescriptor(NumRxnDescriptor, PredictedDescriptor):

    """Numeric Predicted Reaction Descriptor."""

    predictionOf = models.ForeignKey(NumRxnDescriptor, related_name="prediction_of")

    class Meta:
        app_label = 'DRP'
        verbose_name = 'Predicted Numeric Rxn Descriptor'

    objects = DescriptorManager()

    def getPredictionTuples(self):
        """
        Return the prediction tuples.
 
        Return a list of tuples where the first value is the actual value for
        a descriptor of a reaction and the second value is the predicted value
        of that descriptor in the same reaction.
        EG: [(1,1), (2,1), (4,1), (3,1), (1,1)] for a model that always
            predicts "1" if there are 4 different values for a descriptor.
        """

        actualDescValues = self.predictionOf.ordrxndescriptorvalue_set.all()
        predictedDescValues = self.ordrxndescriptorvalue_set.all()

        predictionTuples = []
        for prediction in predictedDescValues:
            actual = actualDescValues.get(reaction=prediction.reaction)
            predictionTuples.append((actual.value, prediction.value))

        return predictionTuples


class PredCatRxnDescriptor(CatRxnDescriptor, PredictedDescriptor):

    """The categorical case of the predicted reaction descriptor."""

    predictionOf = models.ForeignKey(CatRxnDescriptor, related_name="prediction_of")

    class Meta:
        app_label = 'DRP'
        verbose_name = 'Predicted Categorical Rxn Descriptor'

    objects = DescriptorManager()
