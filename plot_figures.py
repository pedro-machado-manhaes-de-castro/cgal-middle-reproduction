# Figures 1-3: Delaunay scaling and paired effects from the recorded measurements.
from pathlib import Path
import argparse
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter


NAMES = {
    "uniform": "Uniforme", "clustered": "Aglomerados",
    "concentrated": "Concentração extrema", "grid": "Grade",
    "collinear": "Colineares", "duplicates": "Duplicatas",
    "exponential": "Escalas exponenciais", "graded": "Escalas múltiplas (2D)",
    "roads_ny": "Nova York", "roads_cal": "Califórnia / Nevada",
    "roads_w": "Oeste dos EUA", "shorelines_europe": "GSHHG: recorte europeu",
}
POLICIES = ["middle-original", "middle-corrected", "median"]
LABELS = {"middle-original": "Middle original", "middle-corrected": "Middle corrigido", "median": "Median"}
COLORS = {"middle-original": "#BD4B38", "middle-corrected": "#167D9A", "median": "#4B8253"}


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def save(figure, output, name):
    for extension in ("pdf", "svg", "png"):
        figure.savefig(output / (name + "." + extension))
    plt.close(figure)


def scaling(rows, datasets, title, output, name, real=False):
    figure, axes = plt.subplots(2, 2, figsize=(11.8, 8.0), layout="constrained")
    for letter, axis, dataset in zip("abcd", axes.flat, datasets):
        maximum = max(int(row["n"]) for row in rows if row["dataset"] == dataset)
        for variant in ("median", "middle-corrected", "middle-original"):
            selected = sorted((row for row in rows if row["dataset"] == dataset and row["variant"] == variant),
                              key=lambda row: int(row["n"]))
            axis.plot([int(row["n"]) for row in selected], [float(row["total_seconds"]) for row in selected],
                      label=LABELS[variant], color=COLORS[variant],
                      linestyle={"middle-original": (0, (4, 2)), "middle-corrected": "-", "median": "-."}[variant],
                      linewidth=1.8, marker={"middle-original": "o", "middle-corrected": "s", "median": "^"}[variant],
                      markersize=7 if variant == "middle-original" else 4.5,
                      markerfacecolor="none" if variant == "middle-original" else COLORS[variant],
                      markeredgewidth=1.4)
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        axis.set_title(f"({letter}) {NAMES[dataset]}", loc="left", pad=10)
        if real:
            axis.text(0.02, 0.97, "Conjunto completo: " + f"{maximum:,}".replace(",", ".") + " pontos",
                      transform=axis.transAxes, va="top", fontsize=9, color="#414141")
            axis.set_ylim(top=axis.get_ylim()[1] * 1.3)
        axis.set_xlabel("Número de pontos")
        axis.set_ylabel("Ordenação + inserção (s)")
        axis.grid(True, alpha=0.17, which="major")
        axis.grid(True, alpha=0.08, which="minor", axis="y")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    lookup = dict(zip(labels, handles))
    figure.legend([lookup[LABELS[policy]] for policy in POLICIES], [LABELS[policy] for policy in POLICIES],
                  loc="outside lower center", ncol=3, frameon=False)
    figure.suptitle(title + "\nMedianas das sete sementes; ambos os eixos em escala logarítmica", fontsize=13)
    save(figure, output, name)


def effects(rows, output):
    maximum = {name: max(int(row["n"]) for row in rows if row["dataset"] == name) for name in NAMES}
    figure, axes = plt.subplots(1, 3, figsize=(15.0, 8.4), sharey=True, layout="constrained")
    for axis, metric, title, left in zip(axes, ("insertion_path", "insert_seconds", "total_seconds"),
                                         ("(a) Percurso", "(b) Inserção", "(c) Tempo total"), (0.013, 0.003, 0.0037)):
        for index, dataset in enumerate(NAMES):
            row = next(row for row in rows if row["dataset"] == dataset and int(row["n"]) == maximum[dataset]
                       and row["metric"] == metric and row["numerator"] == "middle-corrected"
                       and row["denominator"] == "middle-original")
            ratio, low, high = (float(row[field]) for field in ("ratio", "ci95_low", "ci95_high"))
            color, marker = ("#167D9A", "D") if row["kind"] == "real" else ("#855D9A", "o")
            axis.errorbar(ratio, index, xerr=[[ratio - low], [high - ratio]], fmt=marker,
                          color=color, capsize=3, markersize=5.5, linewidth=1.5, zorder=3)
            label = f"{ratio:.4f}" if ratio < 0.1 else f"{ratio:.3f}"
            if dataset == "concentrated" and metric == "insertion_path":
                label = "≈ 1"
            axis.text(1.03, index, label, transform=axis.get_yaxis_transform(), ha="left", va="center",
                      fontsize=9.5, color="#242424")
        axis.axvline(1, color="#444444", linestyle="--", linewidth=2.0, zorder=1)
        axis.set_xscale("log")
        axis.set_xlim(left, 1.5)
        axis.xaxis.set_major_locator(FixedLocator([0.01, 0.1, 1]))
        axis.xaxis.set_major_formatter(FuncFormatter(lambda value, position: f"{value:g}"))
        axis.xaxis.set_minor_formatter(NullFormatter())
        axis.set_title(title, loc="left", pad=10)
        axis.set_xlabel("Razão corrigido / original\nMenor que 1 favorece o corrigido")
        axis.grid(axis="x", which="major", alpha=0.2)
        axis.axhline(7.5, color="#999999", linewidth=0.7, alpha=0.6)
    axes[0].set_yticks(range(len(NAMES)), [f"{NAMES[name]}\nn = {maximum[name]:,}".replace(",", ".") for name in NAMES])
    axes[0].invert_yaxis()
    figure.suptitle("Maiores casos: percurso e custo das inserções\nMédia geométrica de sete razões pareadas e IC individual de 95%", fontsize=13)
    legend = [Line2D([], [], color="#855D9A", marker="o", linestyle="none", label="Dados sintéticos"),
              Line2D([], [], color="#167D9A", marker="D", linestyle="none", label="Dados reais"),
              Line2D([], [], color="#444444", linestyle="--", linewidth=2, label="Razão 1: mesmo valor")]
    figure.legend(handles=legend, loc="outside lower center", ncol=3, frameon=False)
    save(figure, output, "figure3_paired_effects")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.labelsize": 10,
                         "axes.titlesize": 11, "figure.dpi": 150, "savefig.dpi": 220,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "axes.spines.top": False,
                         "axes.spines.right": False})
    rows = read_csv(args.analysis / "timings_and_paths.csv")
    scaling(rows, ["uniform", "clustered", "concentrated", "graded"], "Tempo total em dados sintéticos",
            args.output, "figure1_synthetic_scaling")
    scaling(rows, ["roads_ny", "roads_cal", "roads_w", "shorelines_europe"], "Tempo total em dados reais",
            args.output, "figure2_real_scaling", real=True)
    effects(read_csv(args.analysis / "paired_ratios.csv"), args.output)


if __name__ == "__main__":
    main()
