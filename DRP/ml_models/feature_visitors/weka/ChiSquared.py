from AbstractWekaFeatureVisitor import AbstractWekaFeatureVisitor

class ChiSquared(AbstractWekaFeatureVisitor):
    wekaCommand = "weka.attributeSelection.ChiSquaredAttributeEval -s weka.attributeSelection.Ranker"

    def wekaTrainCommand(self, arff_file, response_index):
        command = "java weka.attributeSelection.ChiSquaredAttributeEval -s weka.attributeSelection.Ranker -i {} -c {}".format(arff_file, response_index)
        return command
