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

  colors = {
            "ms115.6":"#1f78b4",
            "Average ms115.6 Spawn":"#1f78b4",

            "jho213.20":"#33a02c",
            "Average jho213.20 Spawn":"#33a02c",

            "jho252.5":"#e31a1c",
            "Average jho252.5 Spawn":"#e31a1c",

            "jho148.2":"#ff7f00",
            "Average jho148.2 Spawn":"#ff7f00",

            "jho148.2":"#6a3d9a",
            "Average jho148.2 Spawn":"#6a3d9a",

            "Average VSeOx (Intuition)":"#b15928",
            "Average VSeOx (Model)":"#b15928",
            "Average VSeOx Overall":"#b15928",

            "Average VTeOx (Intuition)":"#a6cee3",
            "Average VTeOx (Model)":"#a6cee3",
            "Average VTeOx Overall":"#a6cee3",

            "Average Overall":"#000000",
           }


  for header, y in lines.items():

    color = colors[header] if header in colors else None


    if "overall" in header.lower():
      linestyle = "-."
      weight = 2.0
    elif "intuition" in header.lower():
      linestyle = ":"
      weight = 2.0
    elif "average" in header.lower() or "model" in header.lower():
      linestyle = "--"
      weight = 1.0
    else:
      linestyle = "-"
      weight = 1.0

    if color:
      plt.plot(base, y, label=header, linestyle=linestyle, linewidth=weight, c=color)
    else:
      plt.plot(base, y, label=header, linestyle=linestyle, linewidth=weight)

    if show_mean:
      y_mean = [sum(y)/float(len(y))]*len(base)
      offset_index = int(len(base)/10)

      plt.plot(base, y_mean, linestyle='--')

      plt.annotate("Mean={}".format(y_mean[offset_index]),
                   (base[offset_index], y_mean[offset_index]),
                   textcoords='offset points')

    if show_legend:
      plt.legend(ncol=1, bbox_to_anchor=(1,1),
                 loc="upper right", borderaxespad=0.0,
                 )


  # Taken mostly from: http://stackoverflow.com/questions/22263807/how-is-order-of-items-in-matplotlib-legend-determined
  handles, labels = ax.get_legend_handles_labels()
  labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0], reverse=True))
  ax.legend(handles, labels)



  if xLabel: plt.xlabel(xLabel)

  if yLabel: plt.ylabel(yLabel)

  ax.get_yaxis().get_major_formatter().set_useOffset(False)

  ax.set_yticks(list(frange(major_tick)))
  ax.grid(which='major', alpha=0.5)

  if show_minor:
    ax.set_yticks(list(frange(minor_tick)), minor=True)
    ax.grid(which='minor', alpha=0.3)

  plt.tight_layout()

  return figure


