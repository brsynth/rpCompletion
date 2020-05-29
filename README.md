# rpCompletion

Completes mono-component reactions output by RetroPath2.0 with the appropriate cofactors. Creates sub-paths when multiple reaction rules are associated with a single reaction. Input is a single pathways file produced by RP2Paths. It stands on rpCache which store pre-computed data.

## Input

Required:
* **-rp2paths_pathways**: (string) Path to the rp2paths pathways file
* **-rp2paths_compounds**: (string) Path to the rp2paths compounds file
* **-rp2_pathways**: (string) Path to the RetroPath2.0 pathways file

Advanced options:
* **-upper_flux_bound**: (integer, default=9999) Upper flux bound value
* **-lower_flux_bound**: (integer, default=0) Lower flux bound value
* **-maxSubPaths_filter**: (integer, default=10) Number of subpaths per path
* **-pathway_id**: (string, default=rp_pathway) ID of the heterologous pathway
* **-compartment_id**: (string, default=MNXC3 (i.e. cytoplasm)) Heterologous pathway compartment ID
* **-species_group_id**: (string, default=central_species) ID of the central species, i.e. not cofactors, in the heterologous reactions
* **-sm**: (optional, string, default: file) Store mode. If 'file', rpCache is supposed to be stored in files. Else, the rpCache is supposed to be stored in a Redis database which the name is the value of this input field.

## Output

* **-outputTar**: (string) Path to the output tar.xz file


## Standalone Mode

### Prerequisites

* Python 3 (with `requests` module)

### Quick start
The main code is `src/rpCompletion.py`. Once a scope has been produced by RP2Paths, a typical command line for completing the pathways from the results is:
```
python3 ../src/rpCompletion.py \
  -rp2_pathways rp2_pathways.csv \
  -rp2paths_pathways rp2paths_pathways.csv \
  -rp2paths_compounds rp2paths_compounds.csv \
  -output out \
  -maxSubPaths_filter 10 \
  -sm db
```
where:
- `-rp2_pathways` is the metabolic space outputted by the RetroPath2.0 workflow.
- `-rp2paths_pathways` is the set of pathways outputted by the RP2Paths tool.
- `-rp2paths_compounds` is the set of compounds  outputted by the RP2Paths tool.
- `--output` specify the directory in which all files will be outputted.

### Test
Some tests can be runned. To do so, please follow insructions below:
```
cd test
./run.sh [small | normal | big] <max_subpaths>
```

The first parameter indicates the dataset to use, the second indicates the `maxSubPaths_filter` option.


## Docker Mode

rpCompletion can be run into a docker container.

### Prerequisites

* Docker - [Install](https://docs.docker.com/install/)
* rpCache: `brsynth/rpCache <https://hub.docker.com/r/brsynth/rpcache>`_

### Build image
Before running the container, the image has to be built with:
```
cd docker
docker-compose build
```

### Run
Then, the tool is runnable by:
```
cd docker
./rpCompletion.sh <absolute_indata_folder>
```

Inside the container, rpCompletion can be run following the Standalone Mode.


## Test
All modes can be tested with:
```
cd test
./run[-in-docker].sh
```

## Authors

* **Melchior du Lac**
* **Joan Hérisson**

## Acknowledgments

* Thomas Duigou




## How to cite rpCompletion?

## Licence
rpCompletion is released under the MIT licence. See the LICENCE.txt file for details.


* Docker image: [brsynth/rpcompletion](https://hub.docker.com/r/brsynth/rpcompletion)
