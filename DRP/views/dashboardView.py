"""A module containing the views and associated functions for the dashboard."""
from django.shortcuts import render, redirect
from django.template import RequestContext, Context
from DRP.models import PerformedReaction, LabGroup, ModelContainer
from operator import add

js_testing = False


def dashboard(request):
    """The view for dashboard. Gathers and preprocesses some data for visualization."""
    # get the number of experiments
    num_experiments = len(PerformedReaction.objects.all())
    num_experiments_public = len(PerformedReaction.objects.filter(public=True))
    num_experiments_private = num_experiments - num_experiments_public
    if not js_testing:
        # get the dates the experiments were INSERTED into the datebase
        make_dates_csv('dateEntered.csv', 'inserted')

        # get the date the experiments were PERFORMED
        make_dates_csv('datePerformed.csv', 'performed')

        # get the date the experiments were PERFORMED, count the num experiments cumulatively
        make_dates_csv('datePerformedCumulative.csv', 'performed', cumulative=True)

    # Another thing to add would be the confusion matrices.
    # Unfortunately, there is no data for these in the test dataset I have, so this just prints an empty list
    print(ModelContainer.objects.all())
    for model in ModelContainer.objects.all():
        print(model)
        print(model.getOverallConfusionMatrices())

    return render(request, template_name='dashboard.html',
                  context={'num_experiments': num_experiments,
                           'num_experiments_public': num_experiments_public,
                           'num_experiments_private': num_experiments_private})


def make_dates_csv(csv_name, inserted_or_performed, cumulative=False):
    """Gather data from a csv.

    End result is a csv in static/csv_name with headers:
    date, Norquist Group,...,allthelabgroups
    """
    with open('static/' + csv_name, 'w') as f:
        # get all the potential labgroups
        all_labGroups = list(map(str, [lab_group for lab_group in LabGroup.objects.all()]))

        # give them an index; this is how we will differentiate them in a list later
        lab_group_index_dict = {}
        index = 0
        for lab_group in all_labGroups:
            lab_group_index_dict[lab_group] = index
            index += 1

        # generate a dictionary of dates (year month day format) that has the counts of each of the labgroups as values
        date_dictionary = {}
        for reaction in PerformedReaction.objects.all():
            if inserted_or_performed == "inserted":
                date = str(reaction.insertedDateTime)[:10]
            elif inserted_or_performed == "performed":
                date = str(reaction.performedDateTime)[:10]
            else:
                raise ValueError("Arguement 'inserted_or_performed' must be 'inserted' or 'performed'")

            # many dates are 'None', I am just skipping those
            if date == 'None':
                continue
            user = reaction.user
            labGroups = []
            if user:
                for labGroup in user.labgroup_set.all():
                    labGroups.append(labGroup)
            labGroups = list(map(str, labGroups))
            if date not in date_dictionary:
                date_dictionary[date] = len(all_labGroups) * [0]
            for lab_group in labGroups:
                        lab_group_index = lab_group_index_dict[lab_group]
                        date_dictionary[date][lab_group_index] += 1

        # write the headers
        f.write("date," + ",".join(all_labGroups) + "\n")
        if cumulative:
            previous = len(all_labGroups) * [0]
            ordered_dates = []
            # write the dates and corresponding entry counts
            for date in date_dictionary:
                ordered_dates.append((date, date_dictionary[date]))
            ordered_dates.sort()
            for date in ordered_dates:
                cumulative_counts = list(map(add, date[1], previous))
                previous = cumulative_counts
                lab_group_counts = map(str, cumulative_counts)
                f.write(date[0] + "," + ",".join(lab_group_counts) + "\n")
        else:
            # write the dates and corresponding entry counts
            for date in date_dictionary:
                lab_group_counts = map(str, date_dictionary[date])
                f.write(date + "," + ",".join(lab_group_counts) + "\n")
