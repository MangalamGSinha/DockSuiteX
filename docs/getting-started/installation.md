# Installation

## Requirements

- **Python**: 3.8 or higher
- **Operating System**: Windows (currently Windows-only)
- **Dependencies**: Automatically installed via pip

## Prerequisites

**Java (JDK/JRE)** is required for **P2Rank** pocket prediction.
Download the latest LTS release from [Adoptium Temurin](https://adoptium.net/temurin/releases).
After installation, make sure Java is added to your **System PATH** so it can be accessed from the command line.

**Verify Java installation:**

```bash
java -version
```

## Install from PyPI

```bash
pip install docksuitex
```

## Install from Source

```bash
git clone https://github.com/MangalamGSinha/DockSuiteX.git
cd DockSuiteX
pip install .
```

## Binary Dependencies

DockSuiteX automatically downloads required binary executables on first import:

- MGLTools (AutoDockTools)
- AutoDock Vina
- AutoDock4 & AutoGrid
- P2Rank
- Open Babel

The binaries are downloaded from the GitHub repository and extracted to the package's `bin/` directory.

If binary download fails, you can manually download from:
```
https://github.com/MangalamGSinha/DockSuiteX_Binaries
```

Extract to: `<package_location>/docksuitex/bin/`