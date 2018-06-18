"""A module containing the views and associated functions for the dashboard"""
from django.shortcuts import render, redirect
from django.template import RequestContext, Context
from DRP.models import PerformedReaction


def dashboard(request):
    # get the number of experiments
    num_experiments = len(PerformedReaction.objects.all())
    num_experiments_public = len(PerformedReaction.objects.filter(public=True))
    num_experiments_private = num_experiments  - num_experiments_public

    # get the dates the experiments were ENTERED
    with open('static/dateEntered.csv', 'w') as f:
        f.write("date, performedBy\n")
        for reaction in PerformedReaction.objects.all():
            user = reaction.user
            labGroups = []
            if user:
                for labGroup in user.labgroup_set.all():
                    labGroups.append(labGroup)
                if len(labGroups) > 1:
                    print(labGroups)
            labGroups = list(map(str, labGroups))
            f.write(str(reaction.insertedDateTime)[:10] + ","
                + str(labGroups) + "\n")

    # get the date the experiments were PERFORMED
    with open('static/datePerformed.csv', 'w') as f:
        f.write("date\n")
        for reaction in PerformedReaction.objects.all():
            f.write(str(reaction.performedDateTime)[:10]+"\n")

    return render(request, template_name='dashboard.html',
                           context={'num_experiments' : num_experiments,
                                    'num_experiments_public' : num_experiments_public,
                                    'num_experiments_private' : num_experiments_private})
