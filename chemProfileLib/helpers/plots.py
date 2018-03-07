import matplotlib.pyplot as plt

def plot_common_xaxis(plots, title = None, filename = None, type="plot"):
    num_of_plots = len(plots)

    fig, axes = plt.subplots(nrows=num_of_plots, ncols=1, sharex=True)

    for i in range(num_of_plots):
        item = plots[i]
        if isinstance(item, dict):
            if type == "plot":
                axes[i].plot(item["data"], linewidth=0.5)
            elif type == "bar":
                axes[i].bar(range(len(item["data"])), item["data"])
            axes[i].set_ylabel(item["label"])
            axes[i].get_yaxis().set_label_coords(1.05, 0.5)
        else:
            axes[i].plot(item)

    if title is not None:
        axes[0].set_title(title)

    if filename is not None:
        plt.savefig(filename, dpi=600)
    else:
        plt.show()

    plt.close(fig)


class StackedPlots():
    title = None
    traces = []
    type = None

    def __init__(self, title=None, type="bar"):
        self.title = title
        self.traces = []
        self.type = type

    def addTrace(self, y, label=None):
        self.traces.append({
            "data": y,
            "label": label
        })

    def plot(self, filename = None):
        plot_common_xaxis(self.traces, self.title, filename, type=self.type)