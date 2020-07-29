A Progressive Approach to Scalar Field Topology
===============================================

This repository contains the proposed implementation described in manuscript:

"A Progressive Approach to Scalar Field Topology"

Authors: Jules Vidal (Sorbonne Université), Pierre Guillou (Sorbonne
Université), Julien Tierny (CNRS, Sorbonne Université)

Link to the companion video: https://youtu.be/mDE3738UfWM

## Examine the code

The implementation code is based on the Topology ToolKit (TTK):
https://topology-tool-kit.github.io/ . Most of the added code is
available in the archive under the src/core/base/ directory hierarchy
and consists in the following folders:
* dynamicTree
* multiresTriangulation
* persistenceDiagram
* scalarFieldCriticalPoint

Please note that persistenceDiagram and scalarFieldCriticalPoint also
contain TTK's implementations.

## Build & Install

The present build instructions have been tested onto an up-to-date
Ubuntu 18.04 Linux distribution. Other distributions may need specific
instructions.

The first step is to install the software dependencies. Copy the
following shell statement into a terminal (omit the $ character):

```bash
$ sudo apt install -y g++ cmake libvtk7-dev python3
```

This might take some time, depending on the number of already
installed packages on the target system.

Now, move the compressed archive to your working directory and
decompress it.

```bash
$ unzip implementation.zip
```

To build the proposed implementation and install it in a local
directory, please use the following commands:

```bash
$ cd implementation/src
$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX=../../install ..
$ make -j$(nproc) install
```

## Fetching the data

Download the compressed archive containing the data, available at this link:

https://nuage.lip6.fr/s/qZpXrbnTbxSnYzq/download

Decompress the archive and move the resulting `data` folder into the
`implementation` folder.

```bash
$ unzip data.zip
$ mv data/ implementation/
```




## Reproducing the paper results

The *data* folder contains all the datasets used in this submission.
The *scripts* folder contains Shell and Python scripts used to generate the quantitative
results proposed (tables 2 and 3).

To reproduce those results, please move to the *scripts* folder and run the
benchmarks (it takes a couple hours, as results are averaged on twelve runs
each): 

```bash
$ cd scripts
$ bash bench_pd.sh
$ bash benchs_cp.sh
```

Once the computation is done, you may run the Python scripts to generate the results, directly as LATEX tables:

```bash
$ python generateTables2and3.py
```

The TEX files *table2.tex* and *table3.tex* will be generated.

## Using the standalone executables

The above installation has provided you with two executables files:
**ttkPersistenceDiagramCmd** and **ttkScalarFieldCriticalPointsCmd**,
available in `install/bin/`.

They can be used directly to try out our implementation of our
method. Simply entering one of those commands should list you the
needed parameters. The outputs are VTK unstructured grid files using
the '.vtu' extension.

For the progressive persistence diagram computation:

```
[CommandLine] Missing mandatory argument:
[CommandLine]   -i <{Input data-sets (*.vti, *vtu, *vtp)}>
[CommandLine]
[CommandLine] Usage:
[CommandLine]   ./ttkScalarFieldCriticalPointsCmd
[CommandLine] Argument(s):
[CommandLine]   [-d <Global debug level (default: 3)>]
[CommandLine]   [-t <Global thread number (default: 8)>]
[CommandLine]   [-P <Use progressive algo (default: 0)>]
[CommandLine]   [-F <Input scalar field identifier (default: 0)>]
[CommandLine]   [-O <Input vertex offset field identifier (default: -1)>]
[CommandLine]   [-S <start decimation (default: 0)>]
[CommandLine]   [-E <end decimation (default: 0)>]
[CommandLine]   [-I <Use Integration Tracking (default: 0)>]
[CommandLine]   -i <{Input data-sets (*.vti, *vtu, *vtp)}>
[CommandLine]   [-o <Output file name base (no extension) (default: `output')>]
[CommandLine] Option(s):
```

For the progressive critical points computation:

```
[CommandLine] Missing mandatory argument:
[CommandLine]   -i <{Input data-sets (*.vti, *vtu, *vtp)}>
[CommandLine]
[CommandLine] Usage:
[CommandLine]   ./ttkPersistenceDiagramCmd
[CommandLine] Argument(s):
[CommandLine]   [-d <Global debug level (default: 3)>]
[CommandLine]   [-t <Global thread number (default: 8)>]
[CommandLine]   [-F <Input scalar field identifier (default: 0)>]
[CommandLine]   [-S <starting decimation level (default: 0)>]
[CommandLine]   [-E <stopping decimation level (default: 0)>]
[CommandLine]   [-P <Use progressive algorithm (default: 0)>]
[CommandLine]   [-O <Input vertex offset field identifier (default: -1)>]
[CommandLine]   -i <{Input data-sets (*.vti, *vtu, *vtp)}>
[CommandLine]   [-o <Output file name base (no extension) (default: `output')>]
[CommandLine] Option(s):
```

A standard exemple for the progressive persistence diagram computation
with minimal verbosity is :

```bash
$ ttkPersistenceDiagramCmd -i input.vti -P 1 -S 10 -E 0 -t 1 -d 1
```

Look at the provided scripts for more examples on how to use those
commands.
