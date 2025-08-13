# Wqq Analysis

## Repository Layout
- `make/` — compile and run macros
- `src/` — C++ source files
- `interface/` — header files
- other folders — scripts, data, results

## Build & Run
Compile and run a job (e.g. tagandprobe):
```bash
make compile-tagandprobe
make run-tagandprobe
```

Other jobs:
- `strangejet`
- `wqq`

List all jobs:
```bash
make list
```

Clean build products:
```bash
make clean
```