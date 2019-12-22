"""A module containing the views and associated functions for the dashboard."""
from django.shortcuts import render, redirect
from django.template import RequestContext, Context
from DRP.models import PerformedReaction, LabGroup, ModelContainer, BoolRxnDescriptorValue, Reaction
from DRP.management.commands.build_model import display_model_results
from DRP.views.dashboardView_plot import *
from DRP.views.dashboardView_csv import *


import datetime
import requests
import numpy as np


# Set to False for faster, non-dynamic page loads (good for testing
# JS/Frontend)
generate_csvs = True


def dashboard(request):
    """The view for dashboard. Gathers and preprocesses some data for visualization."""
    # get the number of experiments
    num_experiments = len(PerformedReaction.objects.all())
    num_experiments_public = len(PerformedReaction.objects.filter(public=True))
    num_experiments_private = num_experiments - num_experiments_public

    if generate_csvs:
        today = datetime.datetime.today()
        timeframe = datetime.timedelta(days=365)
        # to show past year, the call looks like this
        # make_dates_csv_weekly('dateEntered.csv', 
        # 'inserted', dateRange=(today - timeframe, today))

        # get the dates the experiments were INSERTED into the datebase
        make_dates_csv_weekly('dateEntered.csv', 'inserted')

        # get the date the experiments were PERFORMED
        make_dates_csv_weekly('datePerformed.csv',
                              'performed', count_no_data=True)

        # get the date the experiments were PERFORMED, count the num
        # experiments cumulatively
        make_dates_csv_weekly(
            'datePerformedCumulativeByLab.csv', 'performed', cumulative=True)

        # get the date the experiments were PERFORMED, count the num 
        # experiments cumulatively, not by lab though
        # this makes the tiny line graph
        today = datetime.datetime.today()
        timeframe = datetime.timedelta(days=30)
        make_dates_csv_no_lab(
            'datePerformedCumulative.csv',
            'performed',
            cumulative=True,
            add_inbetween_times=True,
            dateRange=(
                today - timeframe,
                today))

        make_valid_reaction_csv("validReactions.csv")

    print("Next will print all the confusion matrices objects.")
    all_models = ModelContainer.objects.all()
    built_models = [model for model in all_models if model.built]

    stats_arr = np.array([
        display_model_results(model, verbose=False) for model in built_models
    ])
    avg_acc, avg_bcr, avg_matthews = np.average(stats_arr, axis=0)

    context = {
        'num_experiments': num_experiments,
        'num_experiments_public': num_experiments_public,
        'num_experiments_private': num_experiments_private,
        'dateEntered': make_stacked_bar("dateEntered.csv"),
        'datePerformed': make_stacked_bar("dateEntered.csv"),
        'noPerformedDate': make_stacked_bar("noPerformedDate.csv"),
        'datePerformedCumulativeByLab': make_area_chart("datePerformedCumulativeByLab.csv"),
        'validReactions': make_valid_reaction_bar_char("validReactions.csv"),
        'avg_acc': avg_acc,
        'avg_bcr': avg_bcr,
        'avg_matthews': avg_matthews,
    }

    return render(request, template_name='dashboard.html', context=context)
