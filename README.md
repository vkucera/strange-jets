# Strange-jet analysis

## Get input files

* Use the `download_from_grid.sh` script to download files from the Grid.
  Execute:
  ```bash
  bash download_from_grid.sh
  ```
  and follow the instructions.
* Create a text file with paths to the input files.

## Prepare the steering script

* Customise the `runLocal-11-real.sh` script.

## Execute the steering script to run the analysis

* Load the AliPhysics environment:
  ```bash
  alienv enter AliPhysics/latest
  ```
* Execute the steering script:
  ```bash
  bash runLocal-(...).sh <number of files>
  ```
* See the output in `FinalOutputs.root`.
  ```bash
  rootbrowse FinalOutputs.root
  ```
