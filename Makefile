# SimEK Development Makefile
# Provides common development commands for the project.
#
# All commands call Julia, so the Makefile works in PowerShell, cmd and POSIX shells alike.
# The task logic is in scripts/dev.jl. The development tools (Aqua, Coverage,
# JuliaFormatter) are installed by `make install-dev` in the shared Julia environment
# `simek-dev`; they are not dependencies of SimEK.
#
# The PhreeqcRM tests require the PhreeqcRM library, set once with
# `set_phreeqcrm_library(path)` (saved in LocalPreferences.toml). It can be overridden on
# the command line:
#     make test SIMEK_PHREEQCRM=C:/path/to/PhreeqcRM.dll

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
install-dev:  ## Instantiate the package, examples and docs environments and install the development tools
	$(DEV) install

# Testing
.PHONY: test
test:  ## Run the unit tests directly in the package environment (fast)
	$(JULIA) --project=. test/runtests.jl

.PHONY: test-pkg
test-pkg:  ## Run the unit tests with Pkg.test, in an isolated test environment
	$(JULIA) --project=. -e "using Pkg; Pkg.test()"

.PHONY: examples
examples:  ## Run all examples (figures in the output/ folder of each group)
	$(DEV) examples

# Documentation
.PHONY: docs
docs:  ## Build the documentation (site in docs/build)
	$(JULIA) --project=docs docs/make.jl

.PHONY: docs-rerun
docs-rerun:  ## Build the documentation, running all examples again
	$(JULIA) --project=docs docs/make.jl rerun

.PHONY: docs-serve
docs-serve:  ## Serve the documentation locally (LiveServer.jl), rebuilt when a page or an example changes
	$(JULIA) --project=docs docs/make.jl serve

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
