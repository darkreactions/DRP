from django.core.management.base import BaseCommand
from DRP.models import Reaction
from multiprocessing import Process
from django import db


def reSave(startIndex, endIndex, stdout):
    for r in Reaction.objects.all()[startIndex: endIndex]:
        r.save()
        stdout.write(str(r.id))


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('threads', type=int, help='Number of threads to spawn.')

    def handle(self, *args, **kwargs):
        threads = kwargs['threads']
        examined = 0
        while examined * 20 < Reaction.objects.all().count():
            processes = []
            for i in range(0, threads):
                db.close_old_connections()
                processes.append(Process(target=reSave, args=((examined * 20) + (i * 20), (examined * 20) + ((i + 1) * 20), self.stdout)))
            for i in range(0, threads):
                processes[i].start()
            for i in range(0, threads):
                processes[i].join()
            examined += threads
