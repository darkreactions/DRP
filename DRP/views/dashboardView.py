"""A module containing the views and associated functions for the dashboard."""
from django.shortcuts import render, redirect
from django.template import RequestContext, Context
from DRP.models import PerformedReaction, LabGroup, ModelContainer
from operator import add
import datetime

# Set to False for faster, non-dynamic page loads (good for testing JS/Frontend)
generate_csvs = True

# Labs included in this list will be totally ignored in the end visualization
DISCLUDED_LABS = ['default_amines']

def dashboard(request):
    """The view for dashboard. Gathers and preprocesses some data for visualization."""
    # get the number of experiments
    num_experiments = len(PerformedReaction.objects.all())
    num_experiments_public = len(PerformedReaction.objects.filter(public=True))
    num_experiments_private = num_experiments - num_experiments_public
    if generate_csvs:
        # get the dates the experiments were INSERTED into the datebase
        make_dates_csv('dateEntered.csv', 'inserted')

        # get the date the experiments were PERFORMED
        make_dates_csv('datePerformed.csv', 'performed', count_no_data=True)

        # get the date the experiments were PERFORMED, count the num experiments cumulatively
        make_dates_csv('datePerformedCumulativeByLab.csv', 'performed', cumulative=True)

        # get the date the experiments were PERFORMED, count the num experiments cumulatively, not by lab though
        today = datetime.datetime.today()
        timeframe = datetime.timedelta(days=30)
        make_dates_csv_no_lab('datePerformedCumulative.csv', 'performed', cumulative=True, add_inbetween_times=True, dateRange=(today - timeframe, today))

    # Another thing to add would be the confusion matrices.
    # Unfortunately, there is no data for these in the test dataset I have, so this just prints an empty list
    print("Next will print all the ModelContainer objects, which might be an empty list.")
    print(ModelContainer.objects.all())
    for model in ModelContainer.objects.all():
        print(model)
        print(model.getOverallConfusionMatrices())

    return render(request, template_name='dashboard.html',
                  context={'num_experiments': num_experiments,
                           'num_experiments_public': num_experiments_public,
                           'num_experiments_private': num_experiments_private})


def make_dates_csv(csv_name, inserted_or_performed, cumulative=False, count_no_data=False):
    """Gather data from a csv.

    End result is a csv in static/csv_name with headers:
    date, Norquist Group,...,allthelabgroups
    """
    no_datePerformed = {}

    # get all the potential labgroups ignoring group in DISCLUDED_LABS
    all_labGroups = []
    for lab_group in LabGroup.objects.all():
        lab_group = str(lab_group)
        if lab_group in DISCLUDED_LABS:
            continue
        all_labGroups.append(lab_group)

    # give them an index; this is how we will differentiate them in a list later
    lab_group_index_dict = {}
    index = 0
    for lab_group in all_labGroups:
        lab_group_index_dict[lab_group] = index
        index += 1



    with open('static/' + csv_name, 'w') as f:
        # generate a dictionary of dates (year month day format) that has the counts of each of the labgroups as values
        date_dictionary = {}
        for reaction in PerformedReaction.objects.all():
            if inserted_or_performed == "inserted":
                date = str(reaction.insertedDateTime)[:10]
            elif inserted_or_performed == "performed":
                date = str(reaction.performedDateTime)[:10]
            else:
                raise ValueError("Argument 'inserted_or_performed' must be 'inserted' or 'performed'")

            # many dates are 'None', I am just skipping those
            if date == 'None' and count_no_data:
                date = str(reaction.insertedDateTime)[:10]
                if date not in no_datePerformed:
                    no_datePerformed[date] = len(all_labGroups) * [0]
                lab_group = str(reaction.labGroup)
                if lab_group in DISCLUDED_LABS:
                    continue
                lab_group_index = lab_group_index_dict[lab_group]
                no_datePerformed[date][lab_group_index] += 1
                continue
            elif date == 'None':
                continue

            if date not in date_dictionary:
                date_dictionary[date] = len(all_labGroups) * [0]

            lab_group = str(reaction.labGroup)
            if lab_group in DISCLUDED_LABS:
                continue
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
    if count_no_data:
        with open("static/noPerformedDate.csv", 'w') as f:
            # write the headers
            f.write("date," + ",".join(all_labGroups) + "\n")
            # write the dates and corresponding entry counts
            for date in no_datePerformed:
                lab_group_counts = map(str, no_datePerformed[date])
                f.write(date + "," + ",".join(lab_group_counts) + "\n")

def make_dates_csv_no_lab(csv_name, inserted_or_performed, cumulative=False, add_inbetween_times=False, dateRange=None):
    """Generate a dates csv without accounting for the different labs.

    Used the create the line graph next to the number of experiments.
    Not tested throughly, but tested.
    """
    with open('static/' + csv_name, 'w') as f:
        date_dictionary = {}
        if add_inbetween_times:
            if dateRange:
                d1 = dateRange[0]
                d2 = dateRange[1]
            else:
                d1 = datetime.date(2015, 7, 8)
                d2 = datetime.date(2018, 1, 21)
            delta = d2 - d1

            for i in range(delta.days + 1):
                date_dictionary[str(d1 + datetime.timedelta(i))[:10]] = 0
        for reaction in PerformedReaction.objects.all():
            if inserted_or_performed == "inserted":
                date = reaction.insertedDateTime
            elif inserted_or_performed == "performed":
                date = reaction.performedDateTime
            else:
                raise ValueError("Argument 'inserted_or_performed' must be 'inserted' or 'performed'")

            if date is None or (dateRange and (date < dateRange[0] or date > dateRange[1])):
                continue

            date = str(date)[:10]
            if date not in date_dictionary:
                date_dictionary[date] = 0
            date_dictionary[date] += 1

        # write the headers
        f.write("date,count\n")
        if cumulative:
            previous = 0
            ordered_dates = []
            # write the dates and corresponding entry counts
            for date in date_dictionary:
                ordered_dates.append((date, date_dictionary[date]))
            ordered_dates.sort()
            for date in ordered_dates:
                cumulative_count = date[1] + previous
                previous = cumulative_count
                f.write(date[0] + "," + str(cumulative_count) + "\n")

        else:
            for date in date_dictionary:
                f.write(date + "," + str(date_dictionary[date]) + "\n")
