LABEL_SIZE  = 20
TITLE_SIZE  = 22
LEGEND_SIZE = 14
TICK_SIZE   = 16
FIG_SIZE    = (6, 6)


def apply_style(ax, xlabel, ylabel, title,
                xlim=None, ylim=None,
                legend_loc='upper right',
                log_y=False):
    ax.set_xlabel(xlabel, fontsize=LABEL_SIZE)
    ax.set_ylabel(ylabel, fontsize=LABEL_SIZE)
    ax.set_title(title, fontsize=TITLE_SIZE)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    if log_y:
        ax.set_yscale('log')
    ax.tick_params(axis='both', labelsize=TICK_SIZE)
    ax.grid(True, linestyle=':', alpha=0.4)
    ax.legend(loc=legend_loc, frameon=False, fontsize=LEGEND_SIZE)
