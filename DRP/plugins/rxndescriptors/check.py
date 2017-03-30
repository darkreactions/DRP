a = []

with open('/home/h205c/DRP/DRP/descdescdescs.csv', 'r') as f:
    for line in f.readlines():
        count = 0
        found = False
        this_string = ''
        while not found:
                char = line[count]
                if char == '\t' or char == '{':
                    not_found = True
                else:
                    this_string += char
                count += 1

        a.append(this_string)


print(a)