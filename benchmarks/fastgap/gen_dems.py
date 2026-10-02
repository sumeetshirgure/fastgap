"""Writes the benchmark DEMs (rotated surface-code memory, uniform circuit noise) to dems/.

Usage: python benchmarks/fastgap/gen_dems.py
"""

import json
import pathlib

import stim

HERE = pathlib.Path(__file__).parent
CONFIG = json.loads((HERE / "index_config.json").read_text())


def main():
    out = HERE / "dems"
    out.mkdir(exist_ok=True)
    for cfg in CONFIG["dems"]:
        p = cfg["p"]
        circuit = stim.Circuit.generated(
            cfg["task"],
            distance=cfg["d"],
            rounds=cfg.get("rounds", cfg["d"]),
            after_clifford_depolarization=p,
            before_round_data_depolarization=p,
            before_measure_flip_probability=p,
            after_reset_flip_probability=p,
        )
        dem = circuit.detector_error_model(decompose_errors=True)
        path = out / f"{cfg['name']}.dem"
        dem.to_file(path)
        print(path, dem.num_detectors, "detectors")


if __name__ == "__main__":
    main()
