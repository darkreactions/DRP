def accuracy(conf):
    correct = 0.0
    total = 0.0
    for true, guesses in conf.items():
        for guess, count in guesses.items():
            if true == guess:
                correct += count
            total += count
    return (correct/total if total != 0 else 0.0)

    
def BCR(conf):
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
            class_accuracy_sum += class_correct/class_total
            num_classes += 1
    
    return (class_accuracy_sum/num_classes if num_classes else 0.0)

def confusionMatrixString(confusionMatrix, headers=True):
    """
    Returns a string that will display a confusionMatrix
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
