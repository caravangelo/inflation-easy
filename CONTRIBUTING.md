# Contributing to InflationEasy

## Scope

Contributions are welcome for bug fixes, performance improvements, documentation, and new examples.

## Development Setup

1. Clone the repository.
2. Build the code:

```bash
make
```

3. Run a short simulation:

```bash
./inflation_easy
```

## Pull Request Guidelines

- Keep pull requests focused on one logical change.
- Preserve physical conventions and output formats unless a change is explicitly intended.
- Update documentation when behavior, parameters, or outputs change.
- Include a short validation note describing what was run and what changed.

## Style

- C++17.
- Favor clear, deterministic behavior over micro-optimizations.
- Keep comments concise and specific to non-obvious logic.

## Reporting Issues

When opening an issue, include:

- Platform and compiler version.
- Active `parameters.h` settings (or a minimal subset).
- Exact command run and relevant logs (`results/output.txt`, `results/info.dat`).
