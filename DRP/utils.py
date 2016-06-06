"""Miscellaneous utility functions for use in DRP"""
from math import sqrt


def average_normalized_conf(confs):
    """
    Turn a list of confusion matrices into a single normalized confusion matrix.
    First normalize all matrices so their entries sum to 1, then average them.
    """
    possible_vals = set(confs[0].keys())
    sum_conf = {true: {guess: 0.0 for guess in possible_vals} for true in possible_vals}
    for conf in confs:
        total = 0.0
        for true, guesses in conf.items():
            for guess, count in guesses.items():
                total += count
        if total == 0:
            raise RuntimeError('Confusion matrix values sum to 0')
        for true, guesses in conf.items():
            for guess, count in guesses.items():
                sum_conf[true][guess] += count / total

    average_conf = {true: {guess: count / len(confs) for guess, count in guesses.items()} for true, guesses in sum_conf.items()}
    return average_conf


def accuracy(conf):
    """
    Compute the accuracy given a confusion matrix.
    """
    correct = 0.0
    total = 0.0
    for true, guesses in conf.items():
        for guess, count in guesses.items():
            if true == guess:
                correct += count
            total += count
    return (correct / total if total != 0 else 0.0)


def BCR(conf):
    """
    Compute the balanced classification rate given a confusion matrix.
    """
    class_accuracy_sum = 0.0
    num_classes = 0.0
    for true, guesses in conf.items():
        class_correct = 0.0
        class_total = 0.0
        for guess, count in guesses.items():
            if true == guess:
                class_correct += count
            class_total += count
        if class_total != 0:
            class_accuracy_sum += class_correct / class_total
            num_classes += 1

    return (class_accuracy_sum / num_classes if num_classes else 0.0)


def Matthews(conf):
    """
    Compute the Matthews coefficient given a confusion matrix.
    Only works for two-class confusion matrices.
    """
    class_accuracy_sum = 0.0
    num_classes = 0.0
    if len(conf) != 2:
        raise NotImplementedError("Matthews Correlation Coefficient is only defined for two-class confusion matrices. For a generalization to multi-class problems, investigate markedness and informedness.")
    # Note it doesn't matter which of these is actually considered true and which false as Matthews Correlation Coefficient is symmetrical
    true, true_guesses = conf.items()[0]
    false, false_guesses = conf.items()[1]

    if len(true_guesses) != 2 or len(false_guesses) != 2:
        raise NotImplementedError("Matthews Correlation Coefficient is only defined for two-class confusion matrices. For a generalization to multi-class problems, investigate markedness and informedness.")

    TP = true_guesses[true]
    FN = true_guesses[false]
    TN = false_guesses[false]
    FP = false_guesses[true]

    PP = TP + FP
    AP = TP + FN
    AN = TN + FP
    PN = TN + FN

    if 0 in [PP, AP, AN, PN]:
        # the correct value when any part of the denominator is zero
        return 0.0

    return (TP * TN - FP * FN) / sqrt(PP * AP * AN * PN)


def confusionMatrixString(confusionMatrix, headers=True):
    """
    Returns a string that will display a confusionMatrix.
    If headers=True, includes the headers as the first row and first column.
    """
    table = confusionMatrixTable(confusionMatrix, headers)
    return ('\n'.join([''.join(['{:^6}'.format(item) for item in row]) for row in table]))


def confusionMatrixTable(confusionMatrix, headers=True):
    """
    Converts a confusion matrix dictionary to a list of lists.
    Primarily for display purposes.
    Each list corresponds to a single true value and contains the
    counts for each predicted value.
    If headers=True, includes the headers as the first row and first column.
    """
    values = confusionMatrix.keys()
    table = [[0 for predicted in values] for true in values]

    for i, true in enumerate(values):
        for j, predicted in enumerate(values):
            table[j][i] = confusionMatrix[true][predicted]

    if headers:
        for j, predicted in enumerate(values):
            table[j].insert(0, str(predicted))
        table.insert(0, [""] + map(str, values))

    return table
