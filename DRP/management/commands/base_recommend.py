import os

import django

import DRP.recommender as rec

import os

os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

django.setup()


from DRP.models.performedReaction import PerformedReaction


from DRP.plugins.rxndescriptors.rxnhash import calculate_many



def do_calc():
    calculate_many(PerformedReaction.objects.all())






if __name__ == '__main__':
    print('hey')
    pass
    do_calc()
