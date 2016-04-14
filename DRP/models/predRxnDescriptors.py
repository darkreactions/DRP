from django.db import models
from rxnDescriptors import BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor
from rxnDescriptorValues import BoolRxnDescriptorValue
from descriptors import DescriptorManager
from ModelContainer import ModelContainer
from StatsModel import StatsModel
from PerformedReaction import PerformedReaction
from django.db.models import F

class PredictedDescriptor(models.Model):
    modelContainer = models.ForeignKey(ModelContainer)
    statsModel = models.ForeignKey(StatsModel, null=True)

    class Meta:
        app_label="DRP"
        abstract = True

    objects = DescriptorManager()


class PredBoolRxnDescriptor(BoolRxnDescriptor, PredictedDescriptor):
    predictionOf = models.ForeignKey(BoolRxnDescriptor, related_name="prediction_of")

    class Meta:
        app_label='DRP'
        verbose_name = 'Predicted Boolean Rxn Descriptor'

    objects = DescriptorManager()

 
    def summarize(self):
        return "Accuracy: {}".format(self.accuracy())

    def accuracy(self):
      conf = self.getConfusionMatrix()

      correct = 0.0
      total = 0.0
      for true, guesses in conf.items():
        for guess, count in guesses.items():
          if true == guess: correct += count
          total += count
      return correct/total

    def oldGetConfusionMatrix(self):
        """
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
        if reactions is None:
            reactions = PerformedReaction.objects.filter(boolrxndescriptorvalue__descriptor=self).distinct()
            
            
        matrix = {
                    True: {True: 0, False: 0},
                    False: {True: 0, False: 0}
                    }

        reactions = reactions.filter(boolrxndescriptorvalue__descriptor=self).annotate(predicted_val=F('boolrxndescriptorvalue__value'))
        reactions = reactions.filter(boolrxndescriptorvalue__descriptor=self.predictionOf).annotate(actual_val=F('boolrxndescriptorvalue__value'))

        for rxn in reactions:
            actual_val = rxn.actual_val
            predicted_val = rxn.predicted_val
            #if not rxn.predicted_val:
                #warnings.warn('Reaction {} does not have a predicted value for this descriptor'.format(rxn))
            #if not rxn.actual_val:
                #warnings.warn('Reaction {} does not have an actual value for this descriptor'.format(rxn))
            #if len(rxn.predicted_val) > 1:
                #raise RuntimeError('More than one predicted value for this reaction. This should be impossible. Your code is wrong.')
            #if len(rxn.actual_val) > 1:
                #raise RuntimeError('More than one actual value for this reaction. This should be impossible. Your code is wrong.')
            #actual_val = rxn.actual_val[0].value
            #predicted_val = rxn.predicted_val[0].value
            if predicted_val is not None and actual_val is not None:
                matrix[actual_val][predicted_val] += 1

        return matrix
        
    def getPredictionTuples(self):
        """"
        Returns a list of tuples where the first value is the actual value for
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
            predictionTuples.append( (actual.value, prediction.value) )

        return predictionTuples



class PredOrdRxnDescriptor(OrdRxnDescriptor, PredictedDescriptor):
    predictionOf = models.ForeignKey(OrdRxnDescriptor, related_name="predition_of")

    class Meta:
        app_label='DRP'
        verbose_name = 'Predicted Ordinal Rxn Descriptor'

    objects = DescriptorManager()
 
    def summarize(self):
        return "Accuracy: {}".format(self.accuracy())

    def accuracy(self):
      conf = self.getConfusionMatrix()

      correct = 0.0
      total = 0.0
      for true, guesses in conf.items():
        for guess, count in guesses.items():
          if true == guess: correct += count
          total += count
      return correct/total

    def getConfusionMatrix(self):
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
        
        matrix = {
                    true: {guess: 0 for guess in xrange(self.minimum, self.maximum+1)}
                    for true in xrange(self.minimum, self.maximum+1)
                    }
        for true, guess in self.getPredictionTuples():
            #if true not in matrix:
                #matrix[true] = {}
            try:
                matrix[true][guess] += 1 #matrix[true][guess]+1 if guess in matrix[true] else 1
            except KeyError:
                print true, guess, type(true), type(guess)
                exit()

        return matrix

    def getPredictionTuples(self):
        """"
        Returns a list of tuples where the first value is the actual value for
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
            predictionTuples.append( (actual.value, prediction.value) )

        return predictionTuples


class PredNumRxnDescriptor(NumRxnDescriptor, PredictedDescriptor):
    predictionOf = models.ForeignKey(NumRxnDescriptor, related_name="prediction_of")

    class Meta:
        app_label='DRP'
        verbose_name = 'Predicted Numeric Rxn Descriptor'

    objects = DescriptorManager()
    

    def getPredictionTuples(self):
        """"
        Returns a list of tuples where the first value is the actual value for
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
            predictionTuples.append( (actual.value, prediction.value) )

        return predictionTuples

class PredCatRxnDescriptor(CatRxnDescriptor, PredictedDescriptor):
    predictionOf = models.ForeignKey(CatRxnDescriptor, related_name="prediction_of")

    class Meta:
        app_label='DRP'
        verbose_name = 'Predicted Categorical Rxn Descriptor'

    objects = DescriptorManager()

