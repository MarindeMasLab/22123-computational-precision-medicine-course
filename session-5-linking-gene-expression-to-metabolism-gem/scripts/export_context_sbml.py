#!/usr/bin/env python3
"""Export a scenario's effective exchange bounds into an SBML model."""

from __future__ import annotations

import argparse
from pathlib import Path

import cobra
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--bounds", required=True, type=Path)
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    model = cobra.io.read_sbml_model(str(args.source))
    bounds = pd.read_csv(args.bounds)
    required = {"context", "reaction_id", "lower_bound", "upper_bound"}
    missing = required.difference(bounds.columns)
    if missing:
        raise ValueError(f"Bounds CSV is missing columns: {sorted(missing)}")
    scenario = bounds[bounds["context"] == args.scenario]
    if scenario.empty:
        raise ValueError(f"No bounds found for scenario: {args.scenario}")

    for row in scenario.itertuples(index=False):
        reaction_id = row.reaction_id
        cobra_reaction_id = reaction_id[2:] if reaction_id.startswith("R_") else reaction_id
        if cobra_reaction_id not in model.reactions:
            raise ValueError(f"Unknown reaction in bounds CSV: {reaction_id}")
        reaction = model.reactions.get_by_id(cobra_reaction_id)
        reaction.lower_bound = float(row.lower_bound)
        reaction.upper_bound = float(row.upper_bound)

    model.id = args.scenario
    model.name = f"Recon1 context-specific model: {args.scenario}"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cobra.io.write_sbml_model(model, str(args.output))


if __name__ == "__main__":
    main()
