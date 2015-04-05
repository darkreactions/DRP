def get_graph(lines, base,
              xLabel="Percentage",
              yLabel=None,
              tick_range=(0,1),
              major_tick=0.1, minor_tick=0.025,
              show_minor=True,
              show_legend=True,
              show_mean=False,
              ):

  def frange(step):
    bottom, top = tick_range

    while bottom<top:
      yield bottom
      bottom += step

    yield top

  import matplotlib.pyplot as plt
  plt.ioff()

  figure = plt.figure(figsize=(20,10))
  figure.patch.set_alpha(0)

  ax = figure.add_subplot(1,1,1)

  for header, y in lines.items():
    plt.plot(base, y, label=header)

    if show_mean:
      y_mean = [sum(y)/float(len(y))]*len(base)
      offset_index = int(len(base)/10)

      print "{}: {}".format(header, y_mean)

      plt.plot(base, y_mean, linestyle='--')

      plt.annotate("Mean={}".format(y_mean[offset_index]),
                   (base[offset_index], y_mean[offset_index]),
                   textcoords='offset points')

    if show_legend:
      plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                 fancybox=True, shadow=True, ncol=10)

  if xLabel: plt.xlabel(xLabel)

  if yLabel: plt.ylabel(yLabel)

  ax.get_yaxis().get_major_formatter().set_useOffset(False)

  ax.set_yticks(list(frange(major_tick)))
  ax.grid(which='major', alpha=0.5)

  if show_minor:
    ax.set_yticks(list(frange(minor_tick)), minor=True)
    ax.grid(which='minor', alpha=0.3)

  return figure


