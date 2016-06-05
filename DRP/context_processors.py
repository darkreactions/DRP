from django.conf import settings


def testing(request):
    if settings.TESTING:
        return {'testing': True}
    else:
        return {'testing': False}
