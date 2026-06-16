# WoP: Walk on Probes

This repository contains the official implementation of our **SIGGRAPH 2026** paper
"[Probe-based Walk on Spheres for Efficient Path Reusing](https://t7imal.github.io/projects/2026wop/)".

## Build

Our code was developed and tested on Windows 11 with Clang 20.1.0.
To build the project, initialize the submodules and run CMake from the
repository's root directory:

```bash
git submodule update --init --recursive
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## Usage

We provide demo scenes under `data/exp_scenes/` and solver configurations under
`configs/`. For example, to run the algorithm described in the paper, you can execute:

```bash
./build/bin/Release/demo.exe --config configs/wopc.json
```

You can also override the parameters defined in the configuration file. For example:

```bash
./build/bin/Release/demo.exe --config configs/wopc.json -s data/exp_scenes/2d_1_engine
```

## Contact

If you find any compilation issues or bugs, please feel free to open an issue.
