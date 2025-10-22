# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build System

This project uses autotools (autoconf/automake) as the build system:

```bash
./configure [options]
make
make install
```

**Common Configure Options:**
- `--with-sst-core=<path>` - Enable SST simulation framework integration
- `--with-sst-elements=<path>` - Required if using SST core
- `--with-throttler-headers=<path>` - Enable throttler support (not tested in a decade)
- Various MPI, PAPI, and EPA tools flags are automatically detected

**Key Make Targets:**
- `make clean` - Remove build artifacts
- `make install` - Install libraries to configured directory

## Architecture Overview

This is the **epa-inst-libs** project - instrumentation libraries for PEBIL (Performance Engineering Binary Instrumentation Library) and EPAX. The architecture centers around dynamic binary instrumentation for performance analysis.

### Core Components

**InstrumentationCommon** (`InstrumentationCommon.hpp/.cpp`)
- Foundation for all instrumentation libraries
- Defines tool lifecycle hooks: `tool_dynamic_init`, `tool_image_init`, `tool_thread_init`, etc.
- Handles MPI/SHMEM initialization wrapping
- Provides task ID management and output utilities
- Central error handling and logging macros

**AddressStreamDriver** (`AddressStreamDriver.hpp/.cpp`)
- Main orchestrator for memory address analysis tools
- Manages multiple analysis tools simultaneously (cache simulation, reuse distance, spatial locality, etc.)
- Handles sampling methods and dynamic instrumentation points
- Supports both code-centric and data-centric analysis modes
- Coordinates thread-safe buffer processing

**Analysis Tools Architecture:**
- **AddressStreamBase** - Base class for memory stream analysis tools
- **CacheSimulation** - Multi-level cache hierarchy modeling with various replacement policies
- **ReuseDistance** - Memory reuse pattern analysis (separate subproject in `ReuseDistance/`)
- **AddressRange** - Memory address range tracking and statistics
- **SpatialLocality** - Spatial memory access pattern analysis
- **ScatterGatherLength** - Vector operation analysis
- **CounterFunctions** - Performance counter integration
- **TimerFunctions** - High-precision timing instrumentation
- **PAPI Integration** - Performance API wrapper functions

### Modular Design

The system uses a plugin-like architecture where analysis tools inherit from base classes and register with the AddressStreamDriver. Tools can be enabled/disabled at configure time or runtime via environment variables.

**Conditional Compilation Features:**
- `HAS_EPA_TOOLS` - Extended EPA analytics tools
- `HAS_DATA_STRUCTURE_MODULE` - Data structure analysis capabilities  
- `HAS_ARIEL_FRONTEND` - SST simulation integration
- `HAVE_MPI` - MPI support for parallel applications
- `HAVE_SHMEM` - SHMEM support

### Threading and Synchronization

The system is designed for multi-threaded applications with thread-safe data collection:
- Per-thread data structures managed by DataManager templates
- FastData template for high-performance buffering
- Thread lifecycle management through tool_thread_init/fini hooks
- Dynamic instrumentation point management for live code modification

### ReuseDistance Subproject

The `ReuseDistance/` directory contains a separate autotools-configured library for memory reuse distance analysis. It has its own build system and can be used independently.

## Environment Variables

- `PEBIL_OUTPUT_PREFIX` - Controls output file naming
- `PEBIL_ROOT` - Root directory for PEBIL installation
- Various tool-specific environment variables control analysis behavior

## Library Output

The build produces several shared libraries installed to the configured library directory:
- `libpebilruntime.so` - Core runtime support
- `libaddressstream.so` / `libaddressstream_omp.so` - Memory analysis tools
- `libtimer.so`, `libcounter.so` - Performance measurement
- `libpapifunc.so` - PAPI integration (if enabled)
