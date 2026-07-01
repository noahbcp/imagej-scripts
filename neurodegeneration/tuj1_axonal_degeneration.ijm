/* tuj1_axonal_degeneration.ijm
 * Author: Tayla Gibson, Noah B.C. Piper
 * SPDX-License-Identifier: BSD-3-Clause
 * ==============================================================================
 * Implementation of a TUJ1-stained neurite degeneration quantification workflow
 * inspired by the Degeneration Index (DI) methodology described by Clements et al.
 * (2022), eNeuro. https://doi.org/10.1523/ENEURO.0327-21.2022
 *
 * This macro quantifies neurite fragmentation from fluorescence images by:
 *   1. Thresholding TUJ1-labelled neurites.
 *   2. Measuring total TUJ1-positive area.
 *	 3. Identifying neurite fragments using ImageJ particle analysis.
 *   4. Calculating DI based on fragment area relative to the total TUJ1-positive area.
 *
 * Higher DI values indicate a greater proportion of fragmentation.
 *
 * This implementation does not perform soma removal and is intended for datasets
 * where cell bodies contribute minimally to measurements such as images of
 * distal neurites.
 *
 * Thresholding may be performed:
 *   - Automatically on an image-by-image basis using the Phansalkar method.
 *   - Manually, by selecting a threshold from a representative control image and
 *     applying it to all subsequent images in the batch.
 *
 * Fragment identification is performed using two circularity cut-offs
 * (>=0.20 and >=0.30). Fragment size limits follow recommendations from Clements et al.
 * Very small objects (<=4 pixels) are excluded as debris or thresholding artefacts.
 *
 * IMPORTANT:
 *   Particle size thresholds are defined in pixel units. Parameters were tuned
 *   for 1024 × 1024 pixel confocal images (~320 × 320 µm FoV).
 *   Fragment size limits will require adjustment for different image dimensions,
 *   microscope objectives, zoom factors, or pixel sizes.
 * ==============================================================================
 * Input Requirements
 * ==============================================================================
 * Directory containing:
 *   - .tif or .czi images readable through Bio-Formats.
 *
 * Images should:
 *   - Contain fluorescence-labelled neurites (TUJ1).
 *   - Be single-plane.
 * ==============================================================================
 * Output
 * ==============================================================================
 * Tabulated results are written to the input directory as `TUJ1_circularity.csv`.
 *   filename
 *     Source image filename.
 *
 *   total_tuj1_area
 *     Total thresholded TUJ1-positive area after exclusion of particles
 *     <=4 pixels.
 *
 *   frag_X_area
 *     Combined area of particles with circularity X–1.00 and size
 *     5–10000 pixels.
 *
 *   frag_X_count
 *     Number of particles meeting the circularity >=X criteria.
 *
 *   frag_X_DI
 *     Degeneration Index calculated as: frag_X_area / total_tuj1_area
 * ==============================================================================
 */

macro "TUJ1 Fragment Quant." {
	// Init params.
	inputDir = getDirectory("Choose a Directory");
	fileList = getFileList(inputDir);
	outputFile = inputDir + "TUJ1_circularity.csv";
	tuj1Channel = getNumber("Which channel is TUJ1?", 1);
	fileHeader = "filename,well,total_tuj1_area,frag_0.2_area,frag_0.2_count,frag_0.2_DI,frag_0.3_area,frag_0.3_count,frag_0.3_DI";
	manualThresh = getBoolean("Normalise to control threshold?"); // Option to threshold a control and use that threshold for all future images parsed.

	// Write header only once (if file doesn't exist).
	if (!File.exists(outputFile)) {
		File.append(fileHeader, outputFile);
	}

	if (manualThresh) {
		waitForUser("Open and set threshold for a control image. Then hit Ok!");
		getThreshold(lower, upper);
		run("Close All");
	}

	// Enable headless mode
	setBatchMode("hide");

	// Start looping through files in target directory
	for (i = 0; i < fileList.length; i++) {
		if (endsWith(fileList[i], ".czi") || endsWith(fileList[i], ".tif") || endsWith(fileList[i], ".tiff")) {
			roiManager("reset");
			run("Clear Results");
			run("Bio-Formats Importer", "open=[" + inputDir + fileList[i] + "] autoscale color_mode=Colorized view=Hyperstack stack_order=XYCZT");

			Stack.getDimensions(width, height, channels, slices, frames);
			if (slices > 1) {
				continue;
			}
			if (channels > 1) {
				Stack.setChannel(tuj1Channel);
			}

			// Only convert to 8-bit if running auto local threshold.
			if (!manualThresh) {
				run("8-bit");
			}
			stackID = getImageID();

			// Threshold
			if (manualThresh) {
				setThreshold(lower, upper);
			} else {
				run("Auto Local Threshold", "method=Phansalkar radius=10 parameter_1=0 parameter_2=0 white");
			}
			thresholdID = getImageID();

			// Measure total area, subtracting small (area <= 4px) particles.
			// Calculate area of particles to be removed
			run("Analyze Particles...", "size=0-4 pixel circularity=0.0-1.00 show=Masks");
			run("Invert");
			run("Create Selection");
			run("Measure");
			areaToSubtract = getResult("Area", 0);
			run("Clear Results");
			// Now calculate total area
			selectImage(thresholdID);
			run("Create Selection");
			run("Measure");
			totalArea = getResult("Area", 0) - areaToSubtract;
			run("Clear Results");

			// Measure fragmenets (0.2 circularity)
			selectImage(thresholdID);
			run("Analyze Particles...", "size=5-10000 pixel circularity=0.2-1.00 show=Masks display");
			fragCountOne = nResults;
			run("Clear Results");
			run("Invert");
			run("Create Selection");
			run("Measure");
			fragAreaOne = getResult("Area", 0);
			run("Clear Results");
			// DI is total particle area / total area
			degenIndexOne = fragAreaOne / totalArea;

			// Measure fragmenets (0.3 circularity)
			selectImage(thresholdID);
			run("Analyze Particles...", "size=5-10000 pixel circularity=0.3-1.00 show=Masks display");
			fragCountTwo = nResults;
			run("Clear Results");
			run("Invert");
			run("Create Selection");
			run("Measure");
			fragAreaTwo = getResult("Area", 0);
			run("Clear Results");
			degenIndexTwo = fragAreaTwo / totalArea;

			File.append(fileName + "," + totalArea + "," + fragAreaOne + "," + fragCountOne + "," + degenIndexOne + "," + fragAreaTwo + "," + fragCountTwo + "," + degenIndexTwo, outputFile);

			run("Close All");
			run("Clear Results");
		}
	}
	run("Collect Garbage");
}
