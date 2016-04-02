from DRP.ml_models.model_visitors.weka.AbstractWekaModelVisitor import AbstractWekaModelVisitor
    
class J48(AbstractWekaModelVisitor):
    wekaCommand = "weka.classifiers.trees.J48"

class KNN(AbstractWekaModelVisitor):
    wekaCommand = "weka.classifiers.lazy.IBk"

class LinearRegression(AbstractWekaModelVisitor):
    """Weka LinearRegression for predicting numeric responses"""
    wekaCommand = "weka.classifiers.functions.LinearRegression"

class BayesianLogisticRegression(AbstractWekaModelVisitor):
    wekaCommand = "weka.classifiers.bayes.BayesianLogisticRegression"

class LogisticRegression(AbstractWekaModelVisitor):
    wekaCommand = "weka.classifiers.functions.Logistic"

class M5P(AbstractWekaModelVisitor):
    """Weka M5P for predicting numeric responses"""
    wekaCommand = "weka.classifiers.trees.M5P"

class NaiveBayes(AbstractWekaModelVisitor):
    wekaCommand = "weka.classifiers.bayes.NaiveBayes"

class RandomForest(AbstractWekaModelVisitor):
    wekaCommand = "weka.classifiers.trees.RandomForest"

class SVM_PUK(AbstractWekaModelVisitor):
    wekaCommand = "weka.classifiers.functions.SMO"
    
    def __init__(self, *args, **kwargs):
        super(SVM_PUK, self).__init__(*args, **kwargs)

        self.PUK_OMEGA = 1 #0.5
        self.PUK_SIGMA = 1 #7.0

    def wekaTrainOptions(self):
        kernel = '"weka.classifiers.functions.supportVector.Puk -O {} -S {}"'.format(self.PUK_OMEGA, self.PUK_SIGMA)
        return "-K {}".format(kernel)

