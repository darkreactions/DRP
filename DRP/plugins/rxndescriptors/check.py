a = []

with open('/home/h205c/DRP/DRP/descdescdescs.csv', 'r') as f:
    text = f.read()

not_there = []


with open('/home/h205c/Downloads/orthogonal_plus_noCA_info.dsc', 'r') as f:
    for line in f.readlines():
        if line[:-1] in text:
            not_there.append(line[:-1])

print(not_there)