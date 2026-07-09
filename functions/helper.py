import matplotlib.pyplot as plt
import numpy as np
import logging
from tqdm import tqdm
import inspect
from datetime import datetime
import os

def save_plot(fig, file_path: str, save_plots: bool, file_types=["png"], dpi: float = 200, transparent=False, metadata=False):
    """Saves the current plot to specified file formats.

    Parameters
    ----------
    file_path : str
        The base file path for the saved plot.
    save_plots : bool
        Whether to save the plot.
    file_types : list, optional
        A list of file extensions to save (e.g., ["png", "pdf"]). Default is ["png"].
    dpi : float, optional
        The dots per inch for image resolution. Defaults to 200.

    Raises
    ------
    ValueError
        If no plot exists to save.
    """

    if not plt.gcf().get_axes():
        raise ValueError("No plot to save.")
    
    if metadata:
        caller_frame = inspect.stack()[1]
        script_name = os.path.basename(caller_frame.filename)
        now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        # First line: script + date/time
        fig.text(
            0.5, 0.01,
            f"Script: {script_name} | Generated: {now_str}",
            ha="center", va="bottom", fontsize=8
        )

    if save_plots:
        if "png" in file_types:
            plt.savefig(f"{file_path}.png", dpi=dpi, bbox_inches="tight", transparent=transparent)
            print(f"plot was saved as: {file_path}.png")
        if "pdf" in file_types:
            plt.savefig(f"{file_path}.pdf", bbox_inches="tight", transparent=transparent)
            print(f"plot was saved as: {file_path}.pdf")


def add_plot_footer(fig=None, data="unknown source", displacement=-0.1):
    if fig is None:
        fig = plt.gcf()

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    fig.text(
        0.5, displacement,
        f"Created: {timestamp} | Source: {data}",
        fontsize=8,
        alpha=0.6,
        ha="center",
        va="bottom"
    )

def despine(ax):
    """
    Remove the top and right spines from an ax or array of axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes or array-like
        A single axes object or an array of axes (e.g. from plt.subplots).
    """
    axes = np.asarray(ax).flatten()
    for a in axes:
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)

class TqdmLoggingHandler(logging.StreamHandler):
    def emit(self, record):
        msg = self.format(record)
        tqdm.write(msg)

def get_logger(log_level: str, name: str = None) -> logging.Logger:
    """Create and return a logger with the specified level and optional name."""
    level = getattr(logging, log_level.upper(), None)
    if level is None:
        raise ValueError(f"Invalid log level: '{log_level}'")
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    logger.handlers.clear()#if not logger.handlers:
    logger.addHandler(TqdmLoggingHandler())
    
    return logger