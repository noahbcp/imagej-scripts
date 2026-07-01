/* camkar_nn_analysis_batch.ijm
 * Author: Noah B.C. Piper
 * SPDX-License-Identifier: BSD-3-Clause
 * ==============================================================================
 * Batch-processing version of `camkar_nn_analysis.ijm`.
 *
 * Automated analysis of CaMKAR (CaMKII Activity Reporter) timelapse microscopy data.
 * See Reyes Gaido et al., 2023., SciTranslMed (https://doi.org/10.1126/scitranslmed.abq7839)
 *
 * This workflow is designed for experiments using the CaMKAR reporter in
 * live-cell imaging. 488-ex and 405-ex emission intensity from segmented cells
 * is measured across time. Each cell is annoted and tracked in successive 
 * frames via a nearest-neighbour algorithm.
 *
 * The macro automatically:
 *   1. Segments cells using an Otsu-thresholded mask derived from a user-
 *      specified reference channel (e.g. 488-ex emission).
 *   2. Measures fluorescence intensity from two excitation channels for each
 *      detected cell in every frame.
 *   3. Tracks cells through time using centroid-based nearest-neighbour
 *      assignment.
 *   4. Preserves cell identities between frames where centroid displacement
 *      remains below a user-defined threshold.
 *   5. Exports per-cell intensity measurements and spatial coordinates as a
 *      tabulated ImageJ results table.
 *
 * ------------------------------------------------------------------------------
 * Detection and Tracking Strategy
 * ------------------------------------------------------------------------------
 * Cell Segmentation
 *   - Gaussian smoothing (user-defined sigma)
 *   - Otsu automatic thresholding
 *   - Particle analysis to identify cellular ROIs
 *   - Minimum particle size filtering
 *
 * Cell Tracking
 *   - Centroids are calculated for each ROI.
 *   - Objects are linked between adjacent frames using Euclidean distance.
 *   - Existing identities are propagated when the nearest object falls within
 *     the specified maximum tracking distance.
 *   - Objects appearing after frame 1 are discarded to minimise tracking errors
 *     caused by segmentation artefacts or transient detections.
 *
 * ==============================================================================
 * Input Requirements
 * ==============================================================================
 * Directory of images (.czi, .tif) with the following:
 *   - Time dimension (required)
 *   - Two fluorescence channels corresponding to the CaMKAR imaging workflow
 *   - Cells that remain sufficiently stationary between consecutive frames for
 *     nearest-neighbour tracking
 *
 * ==============================================================================
 * Output
 * ==============================================================================
 * A results table `CaMKAR.csv` is exported to the root directory.
 *   - cell_ID: Unique identifier assigned to each tracked cell.
 *   - frame_n: Timelapse frame number.
 *   - intensity_488: Mean intensity from the 488-excited channel.
 *   - intensity_405: Mean intensity from the 405-excited channel.
 *   - position_x: ROI centroid X coordinate (pixels) <for QC>
 *	 - position_y: ROI centroid Y coordinate (pixels) <for QC>
 *
 * These measurements can be used to calculate excitation ratios, normalised
 * reporter responses, and CaMKII activity dynamics over time.
 *
 * ==============================================================================
 * Assumptions and Limitations
 * ==============================================================================
 * Tracking is based solely on centroid proximity and assumes:
 *
 *   - Limited cell movement between frames.
 *   - Minimal cell overlap.
 *   - Minimal appearance or disappearance of cells during acquisition.
 *
 * Particle size thresholds and maximum tracking distance were empirically chosen
 * and should be validated for each experimental dataset.
 *
 * Fluorescence measurements are reported as raw mean intensities and no
 * background subtraction, bleaching correction, or ratio calculations are
 * performed automatically.
 * ==============================================================================
 */

macro "CaMKAR Analysis" {
	
	// Environment
	PARTICLE_SIZE = 150 // Microns. Tunable. Lower limit for particle analysis.
	MAX_DIST = 5  // Pixels. Tunable. Max-distance to consider for nearest-neighbour analysis.
	MAX_VALUE = 2^16 // Integer. Static. Initial length of storage arrays and `bestDist`.
	
	// Initialisation
	
	// BATCH folder select
	inputDir = getDirectory("Choose a Directory");
	list = getFileList(inputDir);
	outputFile = inputDir + "CaMKAR.csv";
	
	// BATCH write output header
	// Write header only once (if file doesn't exist)
	if (!File.exists(outputFile)) {
		File.append("filename,cell_ID,frame_n,intensity_488,intensity_405,position_x,position_y", outputFile);
	}
				
	// BATCH loop thru files in folder
	for (openFile = 0; openFile < list.length; openFile++) {
	    if (endsWith(fileList[i], ".czi") || endsWith(fileList[i], ".tif") || endsWith(fileList[i], ".tiff")) {
	    		run("Bio-Formats Importer", "open=[" + inputDir + list[openFile] + "] autoscale color_mode=Colorized view=Hyperstack stack_order=XYCZT");
	
				imgFilename = getInfo("image.filename");
				Stack.getDimensions(width, height, channels, slices, frames);
				
				// User input, only ask if first-loop of batch.
				if (openFile == 0) {
					channelOne = getNumber("488-excited channel:", 1);
					channelTwo = getNumber("405-excited channel:", 2);
					sigma = getNumber("Sigma radius (Gaussian):", 1);
				}
				// Prepare workspace
				run("ROI Manager...");
				roiManager("reset");
				run("Set Measurements...", "mean centroid redirect=None decimal=3");
				run("Select None");
				
				// Start batch mode
				setBatchMode("hide");
				
				
				// Prepare image
				initID = getImageID();
				run("Duplicate...", "duplicate");
				stackID = getImageID();
				run("Gaussian Blur...", "sigma=" + sigma + " stack");
				run("Duplicate...", "title=488ex duplicate channels=" + channelOne);
				channelOneID = getImageID();
				
				// Generate masks for each frame using channelOne & autothreshold (Otsu)
				selectImage(channelOneID);
				run("Duplicate...", "duplicate title=threshold");
				run("Auto Threshold", "method=Otsu white stack");
				thresholdID = getImageID();
				
				// Storage arrays
				// Current frame
				cellID = newArray(MAX_VALUE); // Stores cell IDs. Each cell has a unique ID.
				cellX = newArray(MAX_VALUE); // X-axis position.
				cellY = newArray(MAX_VALUE); // Y-axis position.
				cellFrame = newArray(MAX_VALUE); // Timelapse frame.
				cellIntensityOne = newArray(MAX_VALUE);
				cellIntensityTwo = newArray(MAX_VALUE);
				// Previous frame
				prevX = newArray(MAX_VALUE);
				prevY = newArray(MAX_VALUE);
				prevID = newArray(MAX_VALUE);
				prevCount = 0;
				
				nextID = 1;
				nCells = 0;
				
				
				for (f = 1; f <= frames; f++) {
					roiManager("reset");
					
					// Segment cells
					selectImage(thresholdID);
					setSlice(f);
					run("Analyze Particles...", "size=" + PARTICLE_SIZE + "-Infinity add clear slice");
					n = roiManager("count"); // Number of cells in current frame
					
					// Current frame arrays
					currX = newArray(n);
					currY = newArray(n);
					currIntensityOne = newArray(n);
					currIntensityTwo = newArray(n);
					currID = newArray(n);
					
					// Select image to measure
					selectImage(stackID);
					
					// Measure all ROIs
					for (i = 0; i < n; i++) {
						run("Clear Results");
						selectImage(stackID);
						roiManager("select", i);
						Stack.setChannel(channelOne);
						run("Measure");
						Stack.setChannel(channelTwo);
						run("Measure");
						
						// Fetch data from results table
						x = getResult("X", 0);
						y = getResult("Y", 0);
						intensityOne = getResult("Mean", 0);
						intensityTwo = getResult("Mean", 1);
						
						// Store in frame arrays
						currX[i] = x;
						currY[i] = y;
						currIntensityOne[i] = intensityOne;
						currIntensityTwo[i] = intensityTwo;
					}
					
					// Tracking via nearest-neighbour
					for (i = 0; i < n; i++) {
						bestDist = MAX_VALUE;
						bestIdx = -1; // The index of the cell in the previous frame that best matches the current cell in the current frame.
						
						// Compare to previous frame
						for (j = 0; j < prevCount; j++) {
							// The Euclidean distance (d) describes the shortest distance between two points on a 2D-plane.
							// d(x,y) = sqrt((x2 - x1)^2 + (y2 - y1)^2)
							dX = currX[i] - prevX[j];
							dY = currY[i] - prevY[j];
							delta = sqrt(dX * dX + dY * dY);
							
							if (delta < bestDist) {
								bestDist = delta;
								bestIdx = j;
							}
						}
						
						if (f == 1) {
						        // Cells present in frame 1 are assigned IDs.
						        currID[i] = nextID;
						        nextID++;
						    } 
						    else if (bestDist < MAX_DIST && bestIdx != -1) {
						        // Only reuse existing IDs
						        // IDs == 0 have been previously skipped and are discarded.
						        if (prevID[bestIdx] != 0) {
						        		currID[i] = prevID[bestIdx];
						        	} else {
						        		continue
				        			}
						    } 
						    else {
						        // If new 'cells' are identified by analyse particles, they are discarded.
						        continue;
						    }
			
						// Store to global storage arrays
				        cellX[nCells] = currX[i];
				        cellY[nCells] = currY[i];
				        cellFrame[nCells] = f;
				        cellIntensityOne[nCells] = currIntensityOne[i];
				        cellIntensityTwo[nCells] = currIntensityTwo[i];
				        cellID[nCells] = currID[i];
				        nCells++;
					}
					
				    // update previous frame
				    prevCount = n;
				    for (i = 0; i < n; i++) {
				        prevX[i] = currX[i];
				        prevY[i] = currY[i];
				        prevID[i] = currID[i];
				    }
				}
				
				// Trim arrays from init length;
				cell_ID = Array.slice(cellID, 0, nCells - 1);
				frame_n = Array.slice(cellFrame, 0, nCells - 1);
				intensity_488 = Array.slice(cellIntensityOne, 0, nCells - 1);
				intensity_405 = Array.slice(cellIntensityTwo, 0, nCells - 1);
				position_x = Array.slice(cellX, 0, nCells - 1);
				position_y = Array.slice(cellY, 0, nCells - 1);
			
				selectImage(initID);
				run("Collect Garbage");
				close("Results");
				outputName = split(imgFilename, ".");
				arrayLength = frame_n.length;
				
				for (idx = 0; idx < arrayLength; idx++) {
					File.append(imgFilename + "," + cell_ID[idx] + "," + frame_n[idx] + "," + intensity_488[idx] + "," + intensity_405[idx] + "," + position_x[idx] + "," + position_y[idx], outputFile);
				}
	    }
	}
}
