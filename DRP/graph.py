
def graph(lines, base, xLabel="Percentage"):
  import matplotlib.pyplot as plt

  for header, y in lines.items():
    plt.plot(base, y, "o", label=header)
    plt.plot(base, y, label=header)
    plt.xlabel(xLabel)
    plt.ylabel(header)
    plt.show()

if __name__=="__main__":
    lines = {}
    base_line = []

    graph(lines, base_line, )

