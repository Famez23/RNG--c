# RNG-cc: Parallel Construction of Relative Neighborhood Graphs

A hybrid CPU-GPU implementation for constructing Relative Neighborhood Graphs (RNG) using CUDA and OpenMP, developed as part of a Master's thesis in High Performance Computing at Universidad de A Coruña.

## Overview

This project implements a parallel pipeline for RNG construction over large 2D point sets. The approach uses Delaunay triangulation as a supergraph to reduce the candidate edge set, then validates RNG edges in parallel on the GPU.

## Components

| File | Description                                       |
|------|---------------------------------------------------|
| `pt.cpp` | Point generation — produces the input dataset |
| `cgal.cpp` |RNG using Delaunay triangulation using CGAL library  (REFRENCE) |
| `dc.cpp` | RNG based sequential processing    ( BASELINE  )      |
| `testing.cu` |RNG using CUDA and OpenMP  |

## Requirements
-**this project was developed in the FINITERRAE III supercomputer of CESGA using 32 cores , a TESLA T4 and 128GB of memory**
- **C++17** compiler with OpenMP support
- **CUDA Toolkit** (for `testing.cu`)
- **CGAL** library with GMP and MPFR
- **jemalloc** (memory allocator for the CUDA executable)

## Build & Run

The Makefile runs the pipeline in order: point generation → triangulation → processing → GPU validation.
```bash
# Build and run the full pipeline
make

make run-pt        # Generate points
make run-cgal      # RNG BASED CGAL
make run-dc        # RNG BASED Personal implementation
make run-testing   # RNG BASED OpenMP-GPU 

# Clean build artifacts and generated CSVs
make clean

# Full rebuild
make rebuild
```

## Pipeline
```
pt.cpp 
  → points.csv
    → cgal.cpp 
      → dc.cpp 
        → testing.cu 
```

## Performance

Tested on datasets up to 30.5M points, achieving up to **52× GPU acceleration** for RNG validation and **12.8× end-to-end speedup** over sequential implementations.