from AbstractWekaFeatureVisitor import AbstractWekaFeatureVisitor

class InfoGain(AbstractWekaFeatureVisitor):

    maxResponseCount = 1

    def wekaTrainCommand(self, arff_file, response_index):
        command = "java weka.attributeSelection.InfoGainAttributeEval -s weka.attributeSelection.Ranker -i {} -c {}".format(arff_file, response_index)
        return command


