import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

def set_nature_style():
    mplstyle_path = Path(__file__).parent.resolve() / "nature.mplstyle"
    plt.style.use(mplstyle_path)
    
    
def plot_isotopologue_dist(dist: np.array,  prediction: np.array, uptake: float, p: float, name: str = "sample_fa", plots_path=Path("./")):
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    ax.bar(x=range(len(dist)), height=dist, label="Raw data", color="sandybrown")
    ax.set_title(name)
    ax.set_xlabel("Isotopologue number")
    ax.set_ylabel("Normalized intensity")
    
    ax.plot(np.arange(0, len(dist), 2)[1:], prediction[1:], label=f"Fit: uptake={uptake:1.2f}, p={p:1.2f}", c="cornflowerblue")
    ax.scatter(x=[0], y=prediction[0], alpha=1, color="cornflowerblue", marker="_", s=50)
    
    xtick_pos = range(len(dist))
    ax.set_xticks(xtick_pos)
    xtick_labels_m = ["M+%d"%d for d in xtick_pos]
    ax.set_xticklabels(xtick_labels_m, fontsize=6)
    ax.tick_params(axis="x", labelrotation=45)
    
    ax.legend()

    plt.savefig(plots_path / f"{name}.png")
    plt.savefig(plots_path / f"{name}.svg")