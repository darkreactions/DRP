"""Weka's model visitor library."""
from DRP.ml_models.model_visitors.weka.AbstractWekaModelVisitor import AbstractWekaModelVisitor


class J48(AbstractWekaModelVisitor):

    """J48 implementation of C4.5 decision tree."""

    wekaCommand = "weka.classifiers.trees.J48"


class KNN(AbstractWekaModelVisitor):

    """K nearest neighbours classifier."""

    wekaCommand = "weka.classifiers.lazy.IBk"


class LinearRegression(AbstractWekaModelVisitor):

    """Weka LinearRegression for predicting numeric responses."""

    wekaCommand = "weka.classifiers.functions.LinearRegression"


class BayesianLogisticRegression(AbstractWekaModelVisitor):

    """Weka bayesian regression for binary outcomes."""

    wekaCommand = "weka.classifiers.bayes.BayesianLogisticRegression"


class LogisticRegression(AbstractWekaModelVisitor):

    """Standard regression for binary outcomes."""

    wekaCommand = "weka.classifiers.functions.Logistic"


class M5P(AbstractWekaModelVisitor):

    """Weka M5P for predicting numeric responses."""

    wekaCommand = "weka.classifiers.trees.M5P"


class NaiveBayes(AbstractWekaModelVisitor):

    """Naive Bayes predictor."""

    wekaCommand = "weka.classifiers.bayes.NaiveBayes"


class RandomForest(AbstractWekaModelVisitor):

    """Weka random forest classifier."""

    wekaCommand = "weka.classifiers.trees.RandomForest"


class SVM_PUK(AbstractWekaModelVisitor):

    """PUK Kernel Support Vector Machine from Weka."""

    wekaCommand = "weka.classifiers.functions.SMO"

    def __init__(self, puk_omega=1, puk_sigma=1, *args, **kwargs):
        """Additional setup specific to support vector machines."""
        super(SVM_PUK, self).__init__(*args, **kwargs)
        self.puk_omega = puk_omega
        self.puk_sigma = puk_sigma

    @property
    def wekaTrainOptions(self):
        """Attach the additional parameters to the command line function."""
        kernel = '"weka.classifiers.functions.supportVector.Puk -O {} -S {}"'.format(
            self.puk_omega, self.puk_sigma)
        return "-K {}".format(kernel)
