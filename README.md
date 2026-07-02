# ImageJ Scripts Repository

This repository contains ImageJ/Fiji macros I've written for various scientific image processing tasks for both my own research and others. Broadly, these scripts automate workflows for confocal fluorescence microscopy image analysis, including: biosensor analysis, neurite tracing, axonal degeneration quantification and batch processing of image files.

Most scripts have some attempt at documentation in the header comments...

## Directory Structure

### biosensors/

Scripts for biosensor analyses.

- **camkar/**
  - [`camkar_nn_analysis.ijm`](biosensors/camkar/camkar_nn_analysis.ijm) - Automated analysis of [CaMKAR (CaMKII Activity Reporter)](https://doi.org/10.1126/scitranslmed.abq7839) timelapse data; segments cells, tracks them via nearest-neighbour, and exports per-cell intensity measurements.
  - [`camkar_nn_analysis_batch.ijm`](biosensors/camkar/camkar_nn_analysis_batch.ijm) - Batch-processing version of `camkar_nn_analysis.ijm` for processing multiple images in a directory.

### neurodegeneration/

Scripts for neurite analysis and axonal degeneration quantification.

- [`motor_neuron_puncta.ijm`](neurodegeneration/motor_neuron_puncta.ijm) - Automated quantification of ChAT-associated puncta in motor neuron somata from spinal cord sections.
- [`tuj1_axonal_degeneration.ijm`](neurodegeneration/tuj1_axonal_degeneration.ijm) - Quantifies neurite fragmentation from TUJ1-stained images using a Degeneration Index based on fragment area.

### utils/

Utility scripts.

- [`batch_transform.ijm`](utils/batch_transform.ijm) - Applies geometric transformations (rotate/flip) to all .tif images in a folder.
- [`czi_to_tif.ijm`](utils/czi_to_tif.ijm) - Batch converts .czi files to .tif format using Bio-Formats.
- [`event_stamper.ijm`](utils/event_stamper.ijm) - Adds text overlay labels to a range of frames in a hyperstack (e.g. to mark drug additions).
- [`plot_intensity_time.ijm`](utils/plot_intensity_time.ijm) - Generates an intensity-over-time plot from an ROI or entire image.
- [`rotate_align.ijm`](utils/rotate_align.ijm) - Rotates an image to align with a drawn straight line selection.

## Attribution & Licensing

If you find anything in this repository useful to your own work, acknowledgement/citation is greatly appreciated.
Please see the Author line in script headers or [citation.cff](citation.cff) for more information.

All files within this repository are subject to the [BSD 3-Clause License](LICENSE).

## Disclaimer

Everything in this repo is provided as-is and is to be used at your own risk.
If you run into any bugs please feel free to submit an issue to this repository or contact me directly.
