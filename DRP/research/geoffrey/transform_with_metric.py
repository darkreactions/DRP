import django
django.setup()
import build_metric
from sys import argv
from DRP.models import MetricContainer

if __name__ == '__main__':
    pk = int(argv[1])
    testSetName = argv[2]
    outfile = argv[3]
    container = MetricContainer.objects.get(pk=pk)
    print "Calling transform"
    build_metric.transform_rxns(
        container, testSetName=testSetName, outfile=outfile, verbose=True)
