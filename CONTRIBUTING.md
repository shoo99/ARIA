# Contributing to ARIA

Thank you for your interest in contributing to ARIA! This document provides guidelines for contributing to the project.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/YOUR-USERNAME/ARIA.git
   cd ARIA
   ```
3. Install in development mode:
   ```bash
   pip install -e ".[dev]"
   ```
4. Create a branch for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## Development Setup

### Requirements

- Python 3.10+
- Docker (for running R-based analysis modules)
- An Anthropic API key (for LLM-dependent features)

### Running Tests

```bash
pytest tests/ -v
```

### Code Style

We use [Ruff](https://github.com/astral-sh/ruff) for linting and formatting:

```bash
ruff check .
ruff format .
```

## Types of Contributions

### Bug Reports

Open an issue with:
- A clear description of the bug
- Steps to reproduce
- Expected vs actual behavior
- Python version and OS

### Feature Requests

Open an issue describing:
- The problem you're trying to solve
- Your proposed solution
- Alternatives you've considered

### Code Contributions

1. Ensure tests pass before submitting
2. Add tests for new functionality
3. Update documentation if needed
4. Submit a pull request with a clear description

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By participating, you agree to uphold this code.

## Questions?

Open an issue or contact the maintainer at shoo99@gmail.com.
