// Macro for calculating FRET (sensitised emission) from 2/3-channel stacks.
// Inputs are donor-em (donor-ex), fret-em (donor-ex) & acceptor-em (acceptor-ex) if available.
// Returns FRET, donor and acceptor emission intensity at drawn ROIs or across an entire (automatically thresholded) field of view.
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
	if (auto_roi) {
		min_roi_area = getNumber("Minimum ROI area? (microns)", 10);
	}
	gaussian_sigma = getNumber("Sigma radius (Gaussian blur):", 0.25);
	exportMode = getBoolean("Export FRET images?");
	
	// Write header only once (if file doesn't exist). All files are measured
	// frame-by-frame (single-frame files simply produce one row), so the
	// frame/time_s columns are always present.
	if (!File.exists(outputFile)) {
		File.append("filename,donor_em,acceptor_em,fret_em,da_ratio,roi,frame,time_s", outputFile);
	}

for (openFile = 0; openFile < list.length; openFile++) {
    if (endsWith(list[openFile], ".czi")) {
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
				run("Analyze Particles...", "size=" + min_roi_area + "-Infinity show=Masks stack");
				run("Invert", "stack");
				frame_roi_idx = newArray(frame_end + 1);
				for (f = frame_start; f <= frame_end; f++) {
					frame_roi_idx[f] = -1;
					setSlice(f);
					run("Median", "radius=10"); // Helps smooth out masks produced by poor auto-thresholding.
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
		
		
			// Array prep
			if (auto_roi) {
				roi_count = 1;
				} else {
					roi_count = roiManager("count");
					}
			array_length = roi_count * frames_to_process;
			
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
			for (f = frame_start; f <= frame_end; f++) {
				// Manager index of this frame's ROI (auto mode). Manual mode has
				// one ROI per ROI Manager entry, so the 0..roi_count-1 loop below
				// selects directly.
				if (auto_roi) {
					roi_idx = frame_roi_idx[f];
					} else {
						roi_idx = 0;
					}
				if (auto_roi == false) {
				    // Measure all ROIs, appending to Results table.
				    for (i = 0; i < roi_count; i++) {
				        roiManager("Select", i);
				        Stack.setFrame(f);
				        roiManager("Update");
				        roiManager("Measure");
			    	}
				} else if (roi_idx >= 0) {
					roiManager("Select", roi_idx); //1st ROI in manager is index = 0 so need to correct this.
					Stack.setFrame(f);
					roiManager("Update");
					roiManager("Measure");
					}
		
			    // Copy to flattened arrays using `base` as an offset (i.e. `base` becomes the starting loop iteration instead of 0).
			    base = (f - frame_start) * roi_count;
			    for (i = 0; i < roi_count; i++) {
			        idx = base + i;
			        if (auto_roi && roi_idx < 0) {
			            donor_df[idx] = NaN; // no ROI for this frame
			            } else {
			                donor_df[idx] = getResult("Mean", i);
			            }
			        frame_df[idx] = f;
			        roi_df[idx] = i + 1;
			        name_df[idx] = img_name;
			
			        dt = Stack.getFrameInterval(); // seconds per frame
			        time_s[idx] = (f - 1) * dt; // returns the time (s) that a frame was imaged at. For frame 1, this time is always 0.
			    }
			    run("Clear Results");
			}
			
			// Acceptor
			selectImage(acceptor_id);
			for (f = frame_start; f <= frame_end; f++) {
				if (auto_roi) {
					roi_idx = frame_roi_idx[f];
					} else {
						roi_idx = 0;
					}
				if (auto_roi == false) {
					for (i = 0; i < roi_count; i++) {
				        roiManager("Select", i);
				        Stack.setFrame(f);
				        roiManager("Update");
				        roiManager("Measure");
			    	}
				} else if (roi_idx >= 0) {
					roiManager("Select", roi_idx);
					Stack.setFrame(f);
					roiManager("Update");
					roiManager("Measure");
					}
			    
			    base = (f - frame_start) * roi_count;
			    for (i = 0; i < roi_count; i++) {
			        idx = base + i;
			        if (auto_roi && roi_idx < 0) {
			            acceptor_df[idx] = NaN;
			            } else {
			                acceptor_df[idx] = getResult("Mean", i);
			            }
			    }
				run("Clear Results");
			}
			
			// FRET
			selectImage(fret_id);
			for (f = frame_start; f <= frame_end; f++) {
				if (auto_roi) {
					roi_idx = frame_roi_idx[f];
					} else {
						roi_idx = 0;
					}
				if (auto_roi == false) {
					for (i = 0; i < roi_count; i++) {
				        roiManager("Select", i);
				        Stack.setFrame(f);
				        roiManager("Update");
				        roiManager("Measure");
			    	}
				} else if (roi_idx >= 0) {
					roiManager("Select", roi_idx);
					Stack.setFrame(f);
					roiManager("Update");
					roiManager("Measure");
				}
			    
			    base = (f - frame_start) * roi_count;
			    for (i = 0; i < roi_count; i++) {
			        idx = base + i;
			        if (auto_roi && roi_idx < 0) {
			            fret_df[idx] = NaN;
			            } else {
			                fret_df[idx] = getResult("Mean", i);
			            }
			    }
			    run("Clear Results");
			}
			
			imageCalculator("Divide create 32-bit stack", "donor_img", "fret_img");
			run("Rainbow RGB");
			close("\\Others");
			fret_calc = getImageID();
			getDimensions(width, height, channels, slices, frames);
		
			// Crops movie to ROIs only if automated ROIs are generated.
			if (auto_roi) {
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
			
			// Array.show uses variable names as col. titles. Rename variables to things more sensible.
			donor_em = donor_df;
			acceptor_em = acceptor_df;
			fret_em = fret_df;
			da_ratio = da_ratio_df;
			roi = roi_df;
			frame_arr = frame_df;
			filename = name_df;
			
			for (idx = 0; idx < array_length; idx++) {
				File.append(filename[idx] + "," + donor_em[idx] + "," + acceptor_em[idx] + "," + fret_em[idx] + "," + da_ratio[idx] + "," + roi[idx] + "," + frame_arr[idx] + "," + time_s[idx], outputFile);
			}
		}
	}

	// Re-enable GUI; shows the FRET result of the last processed file.
	setBatchMode("show");
}
