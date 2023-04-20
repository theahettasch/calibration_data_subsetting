# druglogics-synergy

<!-- badges: start -->
[![Java CI with Maven](https://github.com/druglogics/druglogics-synergy/workflows/Java%20CI%20with%20Maven/badge.svg)](https://github.com/druglogics/druglogics-synergy/actions)
[![License](https://img.shields.io/github/license/druglogics/druglogics-synergy)](https://github.com/druglogics/druglogics-synergy/blob/master/LICENSE)
<!-- badges: end -->

This package is a Java launcher for sequentially executing the [gitsbe](https://github.com/druglogics/gitsbe) and [drabme](https://github.com/druglogics/drabme) packages.
It gets all necessary input from a directory (prior knowledge and molecular observations) and produces synergy predictions for the specified perturbations.

## Tutorial

For examples on how to run this package and visualize the performance of the output synergy predictions see: https://druglogics.github.io/synergy-tutorial/

## Install

Prerequisites: `maven 3.6.0` and `Java 8`.

Follow the installation guide for each respective package (and in that order):

- [gitsbe](https://github.com/druglogics/gitsbe)
- [drabme](https://github.com/druglogics/drabme)

Then:
```
git clone https://github.com/druglogics/druglogics-synergy
cd druglogics-synergy
mvn clean install
```

The above command creates a package `<name>-jar-with-dependencies.jar` file with all dependencies installed, in the `target` directory.

Alternatively, you could just use directly one of the [released packages](https://github.com/druglogics/druglogics-synergy/packages/).

## Example

The recommended way to run this package is to use its `Launcher`.
From the root directory of the repo run (remember to change the `{version}` to the appropriate one, e.g. `1.2.1`):

```
java -cp ./target/synergy-{version}-jar-with-dependencies.jar eu.druglogics.synergy.Launcher --project=test --inputDir=example_run_ags
```

or run the mvn profile directly (same input as the command above through the `pom.xml`):

```
mvn compile -P runExampleAGS
```

You can use the script [run_druglogics_synergy.sh](https://github.com/druglogics/druglogics-synergy/blob/master/run_druglogics_synergy.sh) to run the above java command with input files from the `ags_cascade_1.0` and `ags_cascade_2.0` directories as well.
This script also offers the possibility to test various input configurations, namely changing the number of simulations, the attractor tool, the choice of training data and the method used in `drabme` for the synergy calculations.

## Docker image

Download:
```
docker pull bblodfon/druglogics-synergy:1.0
```

Run example:
```
git clone https://github.com/druglogics/druglogics-synergy
cd druglogics-synergy

docker run -v ${PWD}:/data bblodfon/druglogics-synergy:1.0 --inputDir=/data/example_run_ags
```

So you just need to mount a directory (like `druglogics-synergy` above) to a path in the docker image (`/data`) and use any sub-directory with input files (`/data/example_run_ags`) to run the `druglogics-synergy` software in one go.

In case you need to access the image directly use this command:

```
docker run -v ${PWD}:/data -it --entrypoint /bin/bash bblodfon/druglogics-synergy:1.0
```

## Input

Running the `druglogics-synergy` Launcher with no parameters, generates a usage message with the available options. 

The only required parameter is the `inputDir`, which is the directory that should contain all relevant input files for the Gitsbe and Drabme packages to execute. 
The input files are identified by a case-insensitive *substring* in their name (and thus must be unique).
Below is a list of identifiers for the input files in the given directory:

Input filename identifier | Example filename
------------------------- | ----------------
`training` | toy_ags_training_data.tab
`config` | toy_ags_config.tab
`network` | toy_ags_network.sif
`modeloutputs` | toy_ags_modeloutputs.tab
`drugpanel` | toy_ags_drugpanel.tab
`perturbations` | toy_ags_perturbations.tab

The only non-required parameters is the `project` (project name) that is used as the name of the directory where the output files will be stored.
