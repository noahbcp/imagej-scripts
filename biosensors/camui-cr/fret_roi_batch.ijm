// Macro for calculating FRET (sensitised emission) from 2/3-channel stacks.
// Inputs are donor-em (donor-ex), fret-em (donor-ex) & acceptor-em (acceptor-ex) if available.
// Returns FRET, donor and acceptor emission intensity at drawn ROIs or across an entire (automatically thresholded) field of view.
// In auto-ROI mode, measurements can alternatively be returned per individual ROI: each frame is segmented with
// 'Analyze Particles' and ROIs are matched across the timelapse by centroid nearest-neighbour.
// If provided, accepts alpha & beta factors to correct spectral bleed through.
// DA Ratio is calculated pixel-wise: DA Ratio = Donor_intensity / FRET_intensity
// Author: Noah B.C. Piper

macro "FRET ROI Measurement (BATCH)" {
	

	// Init params.
	inputDir = getDirectory("Choose a Directory");
	list = getFileList(inputDir);
	outputFile = inputDir + "FRET_results.csv";
	donor_channel = getNumber("Donor fluorescence channel:", 1);
	fret_channel = getNumber("FRET fluorescence channel:", 2);
	acceptor_channel = getNumber("Acceptor fluorescence channel:", 3);
	sbt_correction = getBoolean("Correct for spectral bleed through? Requires α and β.");
	if (sbt_correction) {
		alpha = getNumber("Donor BT (α): ", 0);
		beta = getNumber("Acceptor BT (β): ", 0);
	}
	auto_roi = getBoolean("Use which ROI Method?", "Auto (Otsu)", "Manual");
	individual_roi = false;
	if (auto_roi) {
		min_roi_area = getNumber("Minimum ROI area? (microns)", 10);
		// Toggle: segment and measure every individual ROI in each frame (tracked across the timelapse by nearest-neighbour)
		individual_roi = getBoolean("Measure individual ROIs? (nearest-neighbour tracked)", "Individual ROIs", "Whole field of view");
		if (individual_roi) {
			max_dist = getNumber("Max ROI tracking distance? (pixels)", 5);
		}
	}
	gaussian_sigma = getNumber("Sigma radius (Gaussian blur):", 0.25);
	exportMode = getBoolean("Export FRET images?");
	run("Set Measurements...", "area mean min redirect=None decimal=3"); // Validate environment.
	
	// Write header only once (if file doesn't exist). All files are measured
	// frame-by-frame (single-frame files simply produce one row), so the
	// frame/time_s columns are always present. Individual-ROI mode appends the
	// ROI centroid coordinates (position_x/position_y) for QC.
	if (!File.exists(outputFile)) {
		if (individual_roi) {
			File.append("filename,donor_em,acceptor_em,fret_em,da_ratio,roi,frame,time_s,position_x,position_y", outputFile);
		} else {
			File.append("filename,donor_em,acceptor_em,fret_em,da_ratio,roi,frame,time_s", outputFile);
		}
	}

for (openFile = 0; openFile < list.length; openFile++) {
    if (endsWith(list[openFile], ".czi") | endsWith(list[openFile], ".tif")) {
    		run("Bio-Formats Importer", "open=[" + inputDir + list[openFile] + "] autoscale color_mode=Colorized view=Hyperstack stack_order=XYCZT");
			// Init.
			img_name = getInfo("image.filename");
			Stack.getDimensions(width, height, n_channels, n_slices, n_frames);
			Stack.getPosition(channel, zstack, frame);
			
			// Measure all frames (multi-frame files = timecourse; single-frame
			// files just produce one row). No frame range parameter.
			timecourse_mode = n_frames > 1;
			frame_range = "1-" + n_frames;
			frame_start = 1;
			frame_end = n_frames;
			frames_to_process = frame_end - frame_start + 1;  // inclusive
	
			Stack.setFrame(frame_start);
			run("ROI Manager...");
			
			// ROI input
			if (auto_roi == false) {
				run("Z Project...", "projection=[Average Intensity]");
				setTool("freehand");
				roiManager("reset");
		  		roiManager("show all with labels");
				waitForUser("Add ROIs to the image and the ROI Manager (t) then click OK [" + img_name + "].");
				close();
				if (roiManager("count") == 0) {
					exit("No ROIs in the ROI Manager for " + img_name + ".");
				}
				}
		
			// Disable GUI for faster processing.
			setBatchMode("hide");
			
			// Workspace prep
			run("Set Measurements...", "mean");
			roiManager("deselect");
			run("Select None");
			stack_title = getTitle();
			stack_id = getImageID();
			selectImage(stack_id);
			run("Gaussian Blur...", "sigma=" + gaussian_sigma + " stack");
			run("Select None");
			run("Duplicate...", "title=donor_img duplicate channels=" + donor_channel + " slices=" + zstack + " frames=" + frame_range);
			donor_id = getImageID();
			// SBT correction using alpha value
			if (sbt_correction) {
				run("Duplicate...", "duplicate title=donor_alpha");
				run("Multiply...", "value=" + alpha + " stack");
			}
			selectImage(stack_id);
			run("Duplicate...", "title=acceptor_img duplicate channels=" + acceptor_channel + " slices=" + zstack + " frames=" + frame_range);
			acceptor_id = getImageID();
			selectImage(stack_id);
			run("Duplicate...", "title=fret_img duplicate channels=" + fret_channel + " slices=" + zstack + " frames=" + frame_range);
			if (sbt_correction) {
				// SBT correction follows: FRET.corrected = (1 - beta)*D_A - (alpha * D_D)
				// Completed pixel-wise.
				imageCalculator("Subtract stack", "fret_img", "donor_alpha");
				run("Multiply...", "value=" + (1 - beta) + " stack");
			}
			fret_id = getImageID();
			
		
			// Get image frame count from donor image (assumed representative)
			selectImage(donor_id);
			getDimensions(width, height, channels, slices, frames);
			
			// Automatically generate masks using the donor (presumably the 'brightest' ch) if auto_roi mode is selected
			if (auto_roi)  {
				roiManager("reset");
				selectImage(donor_id);
				run("Duplicate...", "duplicate title=donor_threshold");
				threshold_id = getImageID();
				selectImage(threshold_id);
				run("Auto Threshold", "method=Otsu white stack use_stack_histogram");
				if (individual_roi) {
					// Individual-ROI mode: add every particle of every frame to
					// the ROI Manager and record the Manager index range that
					// belongs to each frame.
					frame_first = newArray(frame_end + 1);
					frame_last = newArray(frame_end + 1);
					for (f = frame_start; f <= frame_end; f++) {
						frame_first[f] = roiManager("count");
						setSlice(f);
						run("Analyze Particles...", "size=" + min_roi_area + "-Infinity add slice");
						frame_last[f] = roiManager("count") - 1;
					}
					if (roiManager("count") == 0) {
						exit("No ROIs were detected in " + img_name + ".\nTry lowering the minimum ROI area, or use manual ROIs.");
					}

					// Nearest-neighbour tracking
					MAX_VALUE = 65536; // Static. Initial length of `bestDist` and tracking arrays.
					total_rois = roiManager("count");
					measManIdx = newArray(total_rois); // ROI Manager index of each tracked measurement.
					measID = newArray(total_rois);     // Tracked ROI ID of each tracked measurement.
					measFrame = newArray(total_rois);  // Frame of each tracked measurement.
					measX = newArray(total_rois);
					measY = newArray(total_rois);
					nextID = 1;
					nMeas = 0;
					prevCount = 0;
					prevX = newArray(MAX_VALUE); // Previous-frame centroids.
					prevY = newArray(MAX_VALUE);
					prevID = newArray(MAX_VALUE); // Previous-frame IDs (0 = discarded/untracked).

					for (f = frame_start; f <= frame_end; f++) {
						n = frame_last[f] - frame_first[f] + 1; // ROIs segmented in this frame.
						currX = newArray(n);
						currY = newArray(n);
						currID = newArray(n);

						for (i = 0; i < n; i++) {
							roiManager("select", frame_first[f] + i);
							getSelectionBounds(bx, by, bw, bh);
							currX[i] = bx + bw / 2; // Centroid (bounding-box centre) of this ROI.
							currY[i] = by + bh / 2;

							// Compare to previous frame; find nearest ROI.
							// The Euclidean distance (d) describes the shortest distance
							// between two points on a 2D-plane: d(x,y) = sqrt((x2 - x1)^2 + (y2 - y1)^2)
							bestDist = MAX_VALUE;
							bestIdx = -1; // Index of the previous-frame ROI that best matches this ROI.
							for (j = 0; j < prevCount; j++) {
								dX = currX[i] - prevX[j];
								dY = currY[i] - prevY[j];
								delta = sqrt(dX * dX + dY * dY);
								if (delta < bestDist) {
									bestDist = delta;
									bestIdx = j;
								}
							}

							if (f == frame_start) {
								// ROIs present in the first tracked frame are assigned IDs.
								currID[i] = nextID;
								nextID++;
							} else if (bestDist < max_dist && bestIdx != -1 && prevID[bestIdx] > 0) {
								// Only reuse existing IDs; IDs == 0 were previously
								// discarded and are never propagated.
								currID[i] = prevID[bestIdx];
							} else {
								// If 'new' ROIs appear after the first tracked frame, they are discarded.
								continue;
							}

							// Store tracking record.
							measManIdx[nMeas] = frame_first[f] + i;
							measID[nMeas] = currID[i];
							measFrame[nMeas] = f;
							measX[nMeas] = currX[i];
							measY[nMeas] = currY[i];
							nMeas++;
						}

						// Update previous-frame arrays.
						prevCount = n;
						for (i = 0; i < n; i++) {
							prevX[i] = currX[i];
							prevY[i] = currY[i];
							prevID[i] = currID[i];
						}
					}
					if (nMeas == 0) {
						exit("No ROIs could be tracked in " + img_name + ".\nTry increasing the max tracking distance, or use manual ROIs.");
					}
					// Capture the frame interval now; the source stacks are closed later
					// and its calibration may not survive the image calculations.
					selectImage(stack_id);
					dt = Stack.getFrameInterval(); // seconds per frame
				} else {
					// Whole-field mode: one merged ROI per frame (existing behaviour).
					run("Analyze Particles...", "size=" + min_roi_area + "-Infinity show=Masks stack");
					run("Invert", "stack");
					frame_roi_idx = newArray(frame_end + 1);
					for (f = frame_start; f <= frame_end; f++) {
						frame_roi_idx[f] = -1;
						setSlice(f);
//						run("Median", "radius=10"); // Helps smooth out masks produced by poor auto-thresholding.
						run("Create Selection");
						if (selectionType() == -1) {
							run("Select None");
							continue; // no usable ROI for this frame
						}
						roiManager("add");
						run("Select None");
						frame_roi_idx[f] = roiManager("count") - 1;
					}
					if (roiManager("count") == 0) {
						exit("No ROIs were detected in " + img_name + ".\nTry lowering the minimum ROI area, or use manual ROIs.");
					}
				}
			}
		
		
			// Array prep
			if (individual_roi) {
				array_length = nMeas; // one row per tracked ROI per frame
			} else if (auto_roi) {
				roi_count = 1;
				array_length = roi_count * frames_to_process;
			} else {
				roi_count = roiManager("count");
				array_length = roi_count * frames_to_process;
			}
			
			name_df		= newArray(array_length);
			donor_df = newArray(array_length);
			acceptor_df = newArray(array_length);
			fret_df = newArray(array_length);
			frame_df = newArray(array_length);
			roi_df = newArray(array_length);
			time_s = newArray(array_length);
		
			
			// Processing
			// Fill arrays in "frame-major" order:
			// For each frame (f), all ROI values are stored consecutively:
			// indices [ (f - frame_start)*roi_count  ...  (f - frame_start)*roi_count + (roi_count-1) ]
			
			// Donor
			selectImage(donor_id);
			if (individual_roi) {
				// Measure each tracked ROI, one row per measurement.
				for (m = 0; m < nMeas; m++) {
					donor_df[m] = measureROI(donor_id, measManIdx[m], measFrame[m]);
				}
			}
			if (individual_roi == false) {
			for (f = frame_start; f <= frame_end; f++) {
				// Manager index of this frame's ROI (auto mode). Manual mode has
				// one ROI per ROI Manager entry, so the 0..roi_count-1 loop below
				// selects directly.
				if (auto_roi) {
					roi_idx = frame_roi_idx[f];
					} else {
						roi_idx = 0;
					}
	
			    // Copy to flattened arrays using `base` as an offset (i.e. `base` becomes the starting loop iteration instead of 0).
			    base = (f - frame_start) * roi_count;
			    for (i = 0; i < roi_count; i++) {
			        idx = base + i;
			        if (auto_roi && roi_idx < 0) {
			            donor_df[idx] = NaN; // no ROI for this frame
			            } else if (auto_roi) {
			                donor_df[idx] = measureROI(donor_id, roi_idx, f);
			            } else {
			                donor_df[idx] = measureROI(donor_id, i, f);
			            }
			        frame_df[idx] = f;
			        roi_df[idx] = i + 1;
			        name_df[idx] = img_name;
			
			        dt = Stack.getFrameInterval(); // seconds per frame
			        time_s[idx] = (f - 1) * dt; // returns the time (s) that a frame was imaged at. For frame 1, this time is always 0.
			    }
			    run("Clear Results");
			}
			}
			
			// Acceptor
			selectImage(acceptor_id);
			if (individual_roi) {
				for (m = 0; m < nMeas; m++) {
					acceptor_df[m] = measureROI(acceptor_id, measManIdx[m], measFrame[m]);
				}
			}
			if (individual_roi == false) {
			for (f = frame_start; f <= frame_end; f++) {
				if (auto_roi) {
					roi_idx = frame_roi_idx[f];
					} else {
						roi_idx = 0;
					}
			    
			    base = (f - frame_start) * roi_count;
			    for (i = 0; i < roi_count; i++) {
			        idx = base + i;
			        if (auto_roi && roi_idx < 0) {
			            acceptor_df[idx] = NaN;
			            } else if (auto_roi) {
			                acceptor_df[idx] = measureROI(acceptor_id, roi_idx, f);
			            } else {
			                acceptor_df[idx] = measureROI(acceptor_id, i, f);
			            }
			    }
				run("Clear Results");
			}
			}
			
			// FRET
			selectImage(fret_id);
			if (individual_roi) {
				for (m = 0; m < nMeas; m++) {
					fret_df[m] = measureROI(fret_id, measManIdx[m], measFrame[m]);
				}
			}
			if (individual_roi == false) {
			for (f = frame_start; f <= frame_end; f++) {
				if (auto_roi) {
					roi_idx = frame_roi_idx[f];
					} else {
						roi_idx = 0;
					}
			    
			    base = (f - frame_start) * roi_count;
			    for (i = 0; i < roi_count; i++) {
			        idx = base + i;
			        if (auto_roi && roi_idx < 0) {
			            fret_df[idx] = NaN;
			            } else if (auto_roi) {
			                fret_df[idx] = measureROI(fret_id, roi_idx, f);
			            } else {
			                fret_df[idx] = measureROI(fret_id, i, f);
			            }
			    }
			    run("Clear Results");
			}
			}
			
			imageCalculator("Divide create 32-bit stack", "donor_img", "fret_img");
			run("Rainbow RGB");
			close("\\Others");
			fret_calc = getImageID();
			getDimensions(width, height, channels, slices, frames);
		
			// Crops movie to the ROI only if automated whole-field ROIs are generated.
			// (Individual-ROI mode keeps the full field of view in the FRET export.)
			if (auto_roi && individual_roi == false) {
				for (i = 0; i < frames; i++) {
					roi_idx = frame_roi_idx[frame_start + i];
					if (roi_idx < 0) {
						continue; // keep this frame uncropped (no ROI)
					}
					selectImage(fret_calc);
					roiManager("select", roi_idx); // ROIs are 0-based
					setSlice(i + 1); // slices are 1-based
					run("Clear Outside", "slice"); // keep only inside the ROI
					run("Select None");
				}
				fret_calc = getImageID();
			}
			
			selectImage(fret_calc);
			if (exportMode) {
				run("Duplicate...", "duplicate title=" + img_name + "_FRET_Results");
				save(inputDir + img_name + "_FRET.tif");
			}
			close("Results");
		
			// Calculate D/A emission ratio
			da_ratio_df = newArray(array_length);
			for (i = 0; i < array_length; i++) {
				da_ratio_df[i] = donor_df[i]/fret_df[i];
				}
			
			// Individual-ROI mode: fill the bookkeeping arrays from the tracking records.
			if (individual_roi) {
				for (m = 0; m < array_length; m++) {
					frame_df[m] = measFrame[m];
					roi_df[m] = measID[m];
					time_s[m] = (measFrame[m] - 1) * dt; // frame 1 is always t = 0
					name_df[m] = img_name;
				}
			}
			
			// Array.show uses variable names as col. titles. Rename variables to things more sensible.
			donor_em = donor_df;
			acceptor_em = acceptor_df;
			fret_em = fret_df;
			da_ratio = da_ratio_df;
			roi = roi_df;
			frame_arr = frame_df;
			filename = name_df;
			
			if (individual_roi) {
				pos_x = measX;
				pos_y = measY;
				for (idx = 0; idx < array_length; idx++) {
					File.append(filename[idx] + "," + donor_em[idx] + "," + acceptor_em[idx] + "," + fret_em[idx] + "," + da_ratio[idx] + "," + roi[idx] + "," + frame_arr[idx] + "," + time_s[idx] + "," + pos_x[idx] + "," + pos_y[idx], outputFile);
				}
			} else {
				for (idx = 0; idx < array_length; idx++) {
					File.append(filename[idx] + "," + donor_em[idx] + "," + acceptor_em[idx] + "," + fret_em[idx] + "," + da_ratio[idx] + "," + roi[idx] + "," + frame_arr[idx] + "," + time_s[idx], outputFile);
				}
			}
		}
	}

	// Re-enable GUI; shows the FRET result of the last processed file.
	setBatchMode("show");
}

// Measures the mean intensity of the ROI at ROI Manager index `roi_idx` on the
// given frame of the image `image_id`, and returns the value.
// `roiManager("select")` restores the ROI together with its own stored
// hyperstack position, so the frame must be (re)set AFTER selecting. The plain
// "Measure" command is used instead of `roiManager("Measure")`, which otherwise
// measures at the ROI's stored position regardless of the current frame. Each
// measurement appends one row to the Results table; the freshly appended row is
// therefore read via the last-row index.
function measureROI(image_id, roi_idx, frame) {
	selectImage(image_id);
	roiManager("select", roi_idx);
	Stack.setFrame(frame);
	run("Measure");
	return getResult("Mean", nResults - 1);
}
