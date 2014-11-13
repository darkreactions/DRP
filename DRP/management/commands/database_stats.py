from django.core.management.base import BaseCommand, CommandError
from DRP.models import *
from DRP.retrievalFunctions import get_latest_ModelStats, get_active_recommendations

class Command(BaseCommand):
  def handle(self, *args, **kwargs):
    def test_queries(objects, queries):
      for query in queries:
        self.stdout.write("\t{}={}: {}".format(query.keys()[0],
                                               query.values()[0],
                                               objects.filter(**query).count()
                                               )
                                              )

    import datetime

    self.stdout.write("_________________"*3)
    dt = datetime.datetime.now()
    self.stdout.write("___ Database Stats ({}) ___".format(dt))


    objects = Data.objects.all()
    queries = [{"is_valid":True}]
    self.stdout.write("Data: {}".format(objects.count()))
    test_queries(objects, queries)

    objects = CompoundEntry.objects.all()
    queries = [{"calculations_failed":True}]
    self.stdout.write("CompoundEntry: {}".format(objects.count()))
    test_queries(objects, queries)

    objects = CG_calculations.objects.all()
    self.stdout.write("CG_calculations: {}".format(objects.count()))

    objects = Recommendation.objects.all()
    queries = [{"nonsense":True}, {"hidden":True}, {"saved":True}]
    self.stdout.write("Recommendation: {}".format(objects.count()))
    test_queries(objects, queries)

    objects = get_active_recommendations()
    queries = [{"nonsense":True}, {"hidden":True}, {"saved":True}]
    self.stdout.write("Active Recommendation: {}".format(objects.count()))
    test_queries(objects, queries)

    objects = ModelStats.objects.all()
    self.stdout.write("Models: {}".format(objects.count()))

    # Detail the statistics of the latest model.
    self.stdout.write("Current Active Model...")
    get_latest_ModelStats().summary()



    self.stdout.write("\n")

