from AbstractWekaFeatureVisitor import AbstractWekaFeatureVisitor

class CFS(AbstractWekaFeatureVisitor):
    wekaCommand = "weka.attributeSelection.CfsSubsetEval -s weka.attributeSelection.BestFirst"

class InfoGain(AbstractWekaFeatureVisitor):
    wekaCommand = "weka.attributeSelection.InfoGainAttributeEval -s weka.attributeSelection.Ranker"
    
class ChiSquared(AbstractWekaFeatureVisitor):
    wekaCommand = "weka.attributeSelection.ChiSquaredAttributeEval -s weka.attributeSelection.Ranker"
