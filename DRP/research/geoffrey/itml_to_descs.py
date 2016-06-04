from DRP.models import MetricContainer

containers = MetricContainer.objects.filter(built=True)

for i, container in enumerate(containers):
    print i, container.description


def print_descs(filename, container):
    with open(filename, 'wb') as f:
        for t in container.transformedRxnDescriptors.all():
            f.write(t.heading)
            f.write('\n')
