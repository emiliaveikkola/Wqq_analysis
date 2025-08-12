

# Wqq_analysis

A ROOT-based analysis project. C++ headers live in `interface/`, sources in `src/`, and ROOT macros in `make/`. A Makefile builds a shared library from `src/` and then runs macros that drive the analysis.

## Repository layout
```
interface/     # C++ headers (.h / .hh)
src/           # C++ sources (.cc)
lib/           # Built shared library (auto-created)
obj/           # Object files (auto-created)
make/          # ROOT macros to compile+run analysis
```

## Build the C++ library
Compiles `src/*.cc` with headers from `interface/` and links to:
```
lib/libAnalysis.so
```
Run:
```bash
make build
```

## Available jobs & how they run
Jobs run in the order: **build → compile-macro → run-macro**. The current jobs are:

- **tagandprobe**
  - compile macro: `make/mk_tpcompile.C`
  - run macro:     `make/mk_TagandProbe.C`
- **wqq**
  - compile macro: `make/mk_Wqqcompile.C`
  - run macro:     `make/mk_Wqq.C`
- **strangejet**
  - compile macro: `make/mk_compile.C`
  - run macro:     `make/mk_StrangeJet.C`

### Run everything
```bash
make
```

### Run a single job
```bash
make tagandprobe
make wqq
make strangejet
```

### Only compile or only run a job
```bash
make compile-tagandprobe
make run-tagandprobe
```

### List job → macro mapping
```bash
make list
```

## Adding a new job
1. Put your macros in `make/`, e.g. `make/mk_MyJobcompile.C` and `make/mk_MyJob.C`.
2. Edit the `Makefile` and append:
   - `JOBS += myjob`
   - `COMPILE_myjob := make/mk_MyJobcompile.C`
   - `RUN_myjob     := make/mk_MyJob.C`
3. Run it:
```bash
make myjob
```

## Cleaning build artifacts
Removes `obj/`, `lib/`, and typical ROOT build files:
```bash
make clean
```
---
**Quick start:**
```bash
make build
make tagandprobe   # or: make wqq, make strangejet
```