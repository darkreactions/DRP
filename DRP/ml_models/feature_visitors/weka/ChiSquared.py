from AbstractWekaFeatureVisitor import AbstractWekaFeatureVisitor

class ChiSquared(AbstractWekaFeatureVisitor):

    maxResponseCount = 1

    def wekaTrainCommand(self, arff_file, response_index):
        command = "java weka.attributeSelection.ChiSquaredAttributeEval -s weka.attributeSelection.Ranker -i {} -c {}".format(arff_file, response_index)
        return command
