from AbstractWekaFeatureVisitor import AbstractWekaFeatureVisitor

class InfoGain(AbstractWekaFeatureVisitor):

    def wekaTrainCommand(self, arff_file, response_index):
        command = "java weka.attributeSelection.InfoGainAttributeEval -s weka.attributeSelection.Ranker -i {} -c {}".format(arff_file, response_index)
        return command


