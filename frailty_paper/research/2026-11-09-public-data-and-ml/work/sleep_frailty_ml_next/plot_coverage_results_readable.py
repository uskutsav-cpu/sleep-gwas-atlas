"""Post-evaluation label layout correction; frozen analysis code stays unchanged."""
import matplotlib
matplotlib.use('Agg')
from matplotlib.axes import Axes
import plot_coverage_results

original_set_xticks=Axes.set_xticks

def readable_labels(self,ticks,labels=None,**kwargs):
    if labels==['Independent publication\nBlood K562','Held library in same paper\nOL33 three cell lines']:
        labels=['Independent blood\nK562','Same paper library\nOL33']
    return original_set_xticks(self,ticks,labels,**kwargs)

if __name__=='__main__':
    Axes.set_xticks=readable_labels
    plot_coverage_results.main()
