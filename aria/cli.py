"""
ARIA Command Line Interface

Main entry point for running ARIA analyses.
"""

import click
import logging
from pathlib import Path

from .core.orchestrator import WorkflowOrchestrator
from .core.reasoning import ReasoningEngine


@click.group()
@click.version_option(version="0.1.0")
def main():
    """ARIA: Adaptive Reasoning for Integrated Analysis"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )


@main.command()
@click.option("--input", "-i", required=True, help="Samplesheet CSV or counts matrix")
@click.option("--output", "-o", default="output", help="Output directory")
@click.option("--genome", "-g", default="GRCm39", help="Reference genome")
@click.option("--llm", default="claude", help="LLM backend (claude/gpt4)")
@click.option("--species", default="Mus musculus", help="Species for gene sets")
@click.option("--control", default="WT", help="Control group prefix")
@click.option("--treatment", default="KO", help="Treatment group prefix")
def analyze(input, output, genome, llm, species, control, treatment):
    """Run full adaptive RNA-seq analysis pipeline."""
    click.echo(f"ARIA v0.1.0 — Adaptive Reasoning for Integrated Analysis")
    click.echo(f"Input: {input}")
    click.echo(f"Output: {output}")
    click.echo(f"Genome: {genome}")
    click.echo(f"LLM backend: {llm}")

    config = {
        "input": input,
        "output_dir": output,
        "genome": genome,
        "species": species,
        "control_prefix": control,
        "treatment_label": treatment,
    }

    engine = ReasoningEngine(llm_backend=llm)
    orchestrator = WorkflowOrchestrator(config, engine)
    orchestrator.run(input)

    click.echo("Analysis complete.")


@main.command()
@click.option("--dataset", "-d", required=True,
              type=click.Choice(["seqc", "bottomly", "airway", "fmr1"]),
              help="Benchmark dataset to run")
@click.option("--output", "-o", default="benchmarks/results", help="Output directory")
def benchmark(dataset, output):
    """Run ARIA on a benchmark dataset."""
    click.echo(f"Running benchmark: {dataset}")
    click.echo(f"Output: {output}")
    # TODO: Implement benchmark runner
    click.echo("Benchmark runner not yet implemented.")


if __name__ == "__main__":
    main()
