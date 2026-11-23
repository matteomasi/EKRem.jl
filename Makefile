# EKRem Development Makefile
# Provides common development commands for the project.
#
# All commands call Julia, so the Makefile works in PowerShell, cmd and POSIX shells alike.
# The task logic is in scripts/dev.jl. The development tools (Aqua, Coverage,
# JuliaFormatter) are installed by `make install-dev` in the shared Julia environment
# `ekrem-dev`; they are not dependencies of EKRem.
#
# The PhreeqcRM tests require the PhreeqcRM library, set once with
# `set_phreeqcrm_library(path)` (saved in LocalPreferences.toml). It can be overridden on
# the command line:
#     make test EKREM_PHREEQCRM=C:/path/to/PhreeqcRM.dll

# Variables
JULIA := julia
DEV := $(JULIA) --project=. scripts/dev.jl

.DEFAULT_GOAL := help

# Help target
.PHONY: help
help:  ## Show this help message
	@$(DEV) help

# Setup and Installation
.PHONY: install-dev
install-dev:  ## Instantiate the package environment and install the development tools
	$(DEV) install

# Testing
.PHONY: test
test:  ## Run the test suite
	$(JULIA) --project=. -e "using Pkg; Pkg.test()"

.PHONY: coverage
coverage:  ## Run the tests with coverage tracking; report src/ coverage and write lcov.info
	$(DEV) coverage

# Code Quality
.PHONY: lint
lint:  ## Check package quality with Aqua.jl
	$(DEV) lint

.PHONY: format
format:  ## Format code with JuliaFormatter
	$(DEV) format

.PHONY: format-check
format-check:  ## Check formatting without modifying files
	$(DEV) format-check

# Release
.PHONY: bump-patch
bump-patch:  ## Bump patch version in Project.toml
	$(DEV) bump patch

.PHONY: bump-minor
bump-minor:  ## Bump minor version in Project.toml
	$(DEV) bump minor

.PHONY: bump-major
bump-major:  ## Bump major version in Project.toml
	$(DEV) bump major

# Utility targets
.PHONY: clean
clean:  ## Remove coverage and allocation files
	$(DEV) clean

.PHONY: all
all: install-dev lint test  ## Run common development tasks (install, lint, test)
