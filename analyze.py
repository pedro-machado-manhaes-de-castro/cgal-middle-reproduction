# Tables 2-3 and Figures 1-3: paired Delaunay policy comparisons.
from collections import defaultdict
from pathlib import Path
import csv
import hashlib
import json
import math
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUTPUT = ROOT / "analysis"
OUTPUT.mkdir(exist_ok=True)
POLICIES = ["middle-original", "middle-corrected", "median"]
LABELS = {"middle-original": "Middle original", "middle-corrected": "Middle corrigido", "median": "Median"}
COLORS = {"middle-original": "#BD4B38", "middle-corrected": "#167D9A", "median": "#4B8253"}
NAMES = {"uniform": "Uniforme", "clustered": "Aglomerados", "concentrated": "Concentração extrema",
         "grid": "Grade", "collinear": "Colineares", "duplicates": "Duplicatas", "exponential": "Escalas exponenciais",
         "graded": "Escalas múltiplas (2D)",
         "roads_ny": "Nova York", "roads_cal": "Califórnia / Nevada", "roads_w": "Oeste dos EUA",
         "shorelines_europe": "GSHHG: recorte europeu"}
METRICS = ["sort_seconds", "insert_seconds", "total_seconds", "insertion_path",
           "sort_cpu_seconds", "insert_cpu_seconds", "total_cpu_seconds"]


def write_csv(name, rows):
    with (OUTPUT / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def interval(values, label):
    values = np.asarray(values, dtype=float)
    logs = np.log(values)
    seed = int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "little")
    random = np.random.default_rng(seed)
    draws = random.integers(0, len(values), size=(20000, len(values)))
    means = np.exp(logs[draws].mean(axis=1))
    low, high = np.quantile(means, [0.025, 0.975])
    return float(np.exp(logs.mean())), float(low), float(high)


def main():
    progress = json.loads((RESULTS / "progress.json").read_text())
    schedule = json.loads((RESULTS / "schedule.json").read_text())
    rows = [json.loads(line) for line in (RESULTS / "samples.jsonl").read_text().splitlines()]
    expected = {(x["dataset"], x["n"], x["seed"], x["variant"]): x["repeats"] for x in schedule["schedule"]}
    groups = defaultdict(list)
    for row in rows:
        key = (row["dataset"], row["n"], row["seed"], row["variant"])
        groups[key].append(row)
        if row["input_hash"] != row["output_hash"] or not row["valid"] or not row["validity_checked"]:
            raise RuntimeError("Failed preservation or validity check")
        if row["vertices"] != row["expected_vertices"]:
            raise RuntimeError("Distinct-point coverage mismatch")
        if not math.isclose(row["sort_seconds"] + row["insert_seconds"], row["total_seconds"], rel_tol=1e-12):
            raise RuntimeError("Phase times do not sum to total")
    incomplete = [key for key, count in expected.items() if len(groups[key]) != count]
    if progress["failed"] or incomplete:
        (OUTPUT / "INCOMPLETE.json").write_text(json.dumps({"failed": progress["failed"], "incomplete": incomplete}, indent=2))
        raise RuntimeError("Campaign contains failures or incomplete pairs; preserve censored cases")
    seeds = []
    for key, trials in sorted(groups.items()):
        for field in ["input_hash", "input_order_hash", "insertion_order_hash", "vertices", "dimension", "insertion_path"]:
            if len({row[field] for row in trials}) != 1:
                raise RuntimeError("Technical trials changed geometry or insertion order")
        dataset, n, seed, variant = key
        seeds.append({"dataset": dataset, "kind": trials[0]["kind"], "n": n, "seed": seed, "variant": variant,
                      "trials": len(trials), **{metric: statistics.median(x[metric] for x in trials) for metric in METRICS},
                      "input_path": trials[0]["input_path"], "vertices": trials[0]["vertices"],
                      "input_hash": trials[0]["input_hash"], "input_order_hash": trials[0]["input_order_hash"],
                      "insertion_order_hash": trials[0]["insertion_order_hash"],
                      "peak_rss_kib": max(x["peak_rss_kib"] for x in trials)})
    write_csv("seed_measurements.csv", seeds)
    lookup = {(x["dataset"], x["n"], x["seed"], x["variant"]): x for x in seeds}
    for case in schedule["cases"]:
        for seed in schedule["seeds"]:
            triplet = [lookup[(case["dataset"], case["n"], seed, variant)] for variant in POLICIES]
            for field in ["input_hash", "input_order_hash", "vertices"]:
                if len({x[field] for x in triplet}) != 1:
                    raise RuntimeError("Policies did not use identical inputs and vertex sets")
    summaries = []
    ratios = []
    for case in schedule["cases"]:
        dataset, n = case["dataset"], case["n"]
        for variant in POLICIES:
            subset = [lookup[(dataset, n, seed, variant)] for seed in schedule["seeds"]]
            summaries.append({"dataset": dataset, "kind": case["kind"], "n": n, "variant": variant,
                              "seeds": len(subset), "vertices_min": min(x["vertices"] for x in subset),
                              "vertices_max": max(x["vertices"] for x in subset),
                              **{metric: statistics.median(x[metric] for x in subset) for metric in METRICS},
                              "peak_rss_kib": max(x["peak_rss_kib"] for x in subset)})
        for numerator, denominator in [("middle-corrected", "middle-original"), ("median", "middle-original"),
                                        ("middle-corrected", "median")]:
            for metric in METRICS:
                values = [lookup[(dataset, n, seed, numerator)][metric] / lookup[(dataset, n, seed, denominator)][metric]
                          for seed in schedule["seeds"]]
                estimate, low, high = interval(values, f"{dataset}-{n}-{numerator}-{denominator}-{metric}")
                ratios.append({"dataset": dataset, "kind": case["kind"], "n": n, "numerator": numerator,
                               "denominator": denominator, "metric": metric, "ratio": estimate,
                               "ci95_low": low, "ci95_high": high, "pairs": len(values),
                               "min_pair_ratio": min(values), "max_pair_ratio": max(values),
                               "lower_in_pairs": sum(x < 1 - 1e-12 for x in values)})
    write_csv("timings_and_paths.csv", summaries)
    write_csv("paired_ratios.csv", ratios)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.labelsize": 10,
                         "axes.titlesize": 11, "figure.dpi": 150, "savefig.dpi": 220,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "axes.spines.top": False,
                         "axes.spines.right": False})
    plot_scaling(summaries, ["uniform", "clustered", "concentrated", "graded"],
                 "figure1_synthetic_scaling", "Tempo total de construção em dados sintéticos")
    plot_scaling(summaries, ["roads_ny", "roads_cal", "roads_w", "shorelines_europe"],
                 "figure2_real_scaling", "Tempo total de construção em dados reais")
    plot_effects(ratios)
    max_by_dataset = {x["dataset"]: max(y["n"] for y in summaries if y["dataset"] == x["dataset"]) for x in summaries}
    largest = [x for x in ratios if x["n"] == max_by_dataset[x["dataset"]]
               and x["numerator"] == "middle-corrected" and x["denominator"] == "middle-original"
               and x["metric"] in ["sort_seconds", "insert_seconds", "total_seconds", "insertion_path"]]
    write_csv("largest_case_ratios.csv", largest)
    report = ["Comparação de políticas na construção de Delaunay 2D", "",
              "A campanha mede se a ordenação corrigida reduz o percurso usado pelas inserções e se essa redução compensa seu custo. Cada comparação usa os mesmos pontos e a mesma permutação inicial para middle original, middle corrigido e median.", "",
              "O percurso registrado é a soma das distâncias euclidianas consecutivas na sequência realmente inserida, sem fechar um ciclo. Os relógios cobrem ordenação, inserção e a união contígua dessas fases. Leitura, projeção, amostragem, cópia de entrada, aquecimento, medição de percurso, verificação e destruição ficam fora desses intervalos.", "",
              "Cada caso usa sete sementes pareadas. Há cinco repetições técnicas até 65.536 pontos, três até 262.144 e uma acima desse tamanho. Os controles inteiramente colineares usam uma repetição por semente e tamanhos 4.096 e 16.384: a localização em dimensão 1 ignora a dica e percorre as arestas. A unidade de análise é a semente, usando a mediana das repetições técnicas. As razões são médias geométricas das razões pareadas; os intervalos percentis de bootstrap usam 20.000 reamostragens e confiança individual de 95%, sem correção de multiplicidade. Eles descrevem a variação nas sementes e ordens medidas, com limites para generalização a outras máquinas e populações.", "",
              f"Cobertura: {len(schedule['cases'])} casos de distribuição e tamanho, {len(schedule['schedule'])} execuções de política e {len(rows)} construções cronometradas. Todas passaram por is_valid(), contagem de vértices distintos, preservação dos pontos e soma dos tempos. Pico de memória das execuções: {max(x['peak_rss_kib'] for x in rows) / 1024:.1f} MiB.", "",
              "Tabela 1. As quatro fontes reais, as contagens, licenças, projeções e hashes estão em data/datasets.json. As redes DIMACS fornecem apenas coordenadas dos vértices. As distâncias são medidas na projeção plana indicada, em metros. Dados sintéticos usam as unidades de seus geradores.", "",
              "Tabela 2. Maiores casos; medianas dos tempos em segundos e do percurso nas unidades do conjunto.",
              "Conjunto | n | Política | Ordenação | Inserção | Total | Percurso"]
    for dataset, n in max_by_dataset.items():
        for x in summaries:
            if x["dataset"] == dataset and x["n"] == n:
                report.append(f"{NAMES[dataset]} | {n} | {LABELS[x['variant']]} | {x['sort_seconds']:.6f} | {x['insert_seconds']:.6f} | {x['total_seconds']:.6f} | {x['insertion_path']:.9g}")
    report += ["", "Tabela 3. Razões corrigido/original nos maiores casos. Valores abaixo de 1 favorecem a correção.",
               "Conjunto | Medida | Razão | IC individual de 95% | Pares com redução"]
    for x in largest:
        report.append(f"{NAMES[x['dataset']]} | {x['metric']} | {x['ratio']:.6f} | [{x['ci95_low']:.6f}, {x['ci95_high']:.6f}] | {x['lower_in_pairs']}/{x['pairs']}")
    report += ["", "As medidas separadas permitem observar associação entre percurso e tempo. Elas não identificam isoladamente o custo de localização, das mudanças topológicas e da memória. Os testes incluem casos colineares e duplicatas, nos quais a dimensão ou o número de vértices distintos pode ser pequeno.", "",
               "A comparação usa spatial_sort com seus parâmetros padrão e inserção ponto a ponto com a face do último ponto como dica. A inserção padrão por intervalo de Delaunay_triangulation_2 utiliza median. A mudança avaliada afeta middle em 2D; as conclusões quantitativas desta campanha têm esse alcance.", "",
               "Reprodução: configure CGAL_BASE com o checkout 6.1.x identificado em experiment.json e execute o script python reproduce.py em Linux com CMake, Ninja, GCC, Python, NumPy e Matplotlib. As fontes brutas e os pares projetados permanecem em data; o manifesto registra suas identidades. Os comandos, cabeçalhos efetivamente incluídos, validações e resultados por execução estão em results."]
    (OUTPUT / "RELATORIO.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
    audit = {"status": "COMPLETE_VALIDATED", "cases": len(schedule["cases"]), "policy_processes": len(schedule["schedule"]),
             "timed_constructions": len(rows), "paired_input_checks": len(schedule["cases"]) * len(schedule["seeds"]),
             "maximum_n": max(x["n"] for x in rows), "peak_rss_kib": max(x["peak_rss_kib"] for x in rows),
             "incomplete": incomplete, "instrumentation_validation": json.loads((RESULTS / "instrumentation_validation.json").read_text())["status"],
             "inference": "Individual descriptive bootstrap intervals; no multiplicity correction; no causal isolation of point location.",
             "largest_case_corrected_over_original": largest}
    (OUTPUT / "VALIDATION.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({key: value for key, value in audit.items() if key != "largest_case_corrected_over_original"}, indent=2))


def plot_scaling(summaries, datasets, filename, title):
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.3), layout="constrained")
    for letter, axis, dataset in zip("abcd", axes.flat, datasets):
        for variant in POLICIES:
            rows = sorted([x for x in summaries if x["dataset"] == dataset and x["variant"] == variant], key=lambda x: x["n"])
            axis.plot([x["n"] for x in rows], [x["total_seconds"] for x in rows],
                      marker={"middle-original": "o", "middle-corrected": "s", "median": "^"}[variant],
                      lw=1.8, markersize=4.8, color=COLORS[variant], label=LABELS[variant])
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        axis.set_title(f"({letter}) {NAMES[dataset]}", loc="left")
        axis.set_xlabel("Número de pontos")
        axis.set_ylabel("Ordenação + inserção (s)")
        axis.grid(True, alpha=0.18, which="both")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3, frameon=False)
    fig.suptitle(title, fontsize=14)
    for extension in ["pdf", "svg", "png"]:
        fig.savefig(OUTPUT / f"{filename}.{extension}")
    plt.close(fig)


def plot_effects(ratios):
    maximum = {name: max(x["n"] for x in ratios if x["dataset"] == name) for name in NAMES}
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 7.2), sharey=True, layout="constrained")
    for axis, metric, title in zip(axes, ["insertion_path", "insert_seconds", "total_seconds"],
                                   ["(a) Percurso", "(b) Inserção", "(c) Tempo total"]):
        for i, dataset in enumerate(NAMES):
            row = next(x for x in ratios if x["dataset"] == dataset and x["n"] == maximum[dataset]
                       and x["metric"] == metric and x["numerator"] == "middle-corrected" and x["denominator"] == "middle-original")
            axis.errorbar(row["ratio"], i, xerr=[[row["ratio"] - row["ci95_low"]], [row["ci95_high"] - row["ratio"]]],
                          color="#167D9A" if row["kind"] == "real" else "#855D9A", fmt="o", capsize=3, markersize=5)
        axis.axvline(1, color="#444444", linestyle="--", linewidth=1.5)
        axis.set_xscale("log")
        axis.set_title(title, loc="left")
        axis.set_xlabel("Razão corrigido / original")
        axis.grid(axis="x", which="both", alpha=0.18)
    axes[0].set_yticks(range(len(NAMES)), [f"{NAMES[name]}\nn = {maximum[name]:,}".replace(",", ".") for name in NAMES])
    axes[0].invert_yaxis()
    fig.suptitle("Maiores casos: percurso e custo das inserções\nRazões pareadas e intervalos individuais de 95%", fontsize=13)
    for extension in ["pdf", "svg", "png"]:
        fig.savefig(OUTPUT / f"figure3_paired_effects.{extension}")
    plt.close(fig)


if __name__ == "__main__":
    main()
