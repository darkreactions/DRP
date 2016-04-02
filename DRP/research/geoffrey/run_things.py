import django
django.setup()
from DRP.models import ModelContainer
import build_model
from sys import argv


m = ModelContainer.objects.order_by('-pk')[0]

m.create_duplicate()
