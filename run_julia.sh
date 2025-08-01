#!/bin/bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=20
julia src/run_static_benchmarks.jl