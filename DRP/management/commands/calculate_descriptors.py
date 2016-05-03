from django.core.management.base import BaseCommand
from DRP.models import Reaction, Compound
from django import db


class Command(BaseCommand):

    def handle(self, *args, **kwargs):
        Compound.objects.all().calculate_descriptors(verbose=True)
        Reaction.objects.all().calculate_descriptors(verbose=True)
