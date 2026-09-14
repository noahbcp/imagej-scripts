// Macro for calculating FRET (sensitised emission) from 2/3-channel stacks.
// Inputs are donor-em (donor-ex), fret-em (donor-ex) & acceptor-em (acceptor-ex) if available.
// Returns FRET, donor and acceptor emission intensity at drawn ROIs or across an entire (automatically thresholded) field of view.
// If provided, accepts alpha & beta factors to correct spectral bleed through.
// DA Ratio is calculated pixel-wise: DA Ratio = Donor_intensity / FRET_intensity
// Author: Noah B.C. Piper

macro "FRET ROI Measurement" {
	
	// Init.
	if (nImages == 0) {
		exit("Please open an image before running intensity_roi");
		}
	
	if (nImages > 1) {
		// Warn that all other image windows will be closed by the macro.
		if (getBoolean("All other images except the active window will be closed. Proceed?") == false) {
			exit;
		} else {
			close("\\Others");
			}
	}
	img_name = getInfo("image.filename");
	Stack.getDimensions(width, height, n_channels, n_slices, n_frames);
	Stack.getPosition(channel, zstack, frame);
	
	// Stack param input
	donor_channel = getNumber("Donor fluorescence channel:", 1);
	fret_channel = getNumber("FRET fluorescence channel:", 2);
	acceptor_channel = getNumber("Acceptor fluorescence channel:", 3);
	
	// SBT input
	sbt_correction = getBoolean("Correct for spectral bleed through? Requires α and β.");
	if (sbt_correction) {
		alpha = getNumber("Donor BT (α): ", 0);
		beta = getNumber("Acceptor BT (β): ", 0);
	}

	if (n_frames > 1) {
		timecourse_mode = getBoolean("Measure at multiple timepoints?");
		if (timecourse_mode) {
			frame_range = getString("Which frames to measure (range; X-Y)?", frame + "-" + n_frames);
			str_parts = split(frame_range, "-");
			frame_start = parseInt(str_parts[0]);
			frame_end = parseInt(str_parts[1]);
			} else {
				frame_range = "" + frame + "-" + frame;
				frame_start = frame;
				timecourse_mode = false;
				}
		} else {
			frame_range = "" + frame + "-" + frame;
			frame_start = frame;
			timecourse_mode = false;
		}
	
	Stack.setFrame(frame_start);
	run("ROI Manager...");
	
	// Check if pre-populated ROIs should be used.
	if (roiManager("count") > 0) {
		reuse_roi = getBoolean("Use ROIs already in ROI Manager?");
		auto_roi = false;
		} else {
			reuse_roi = false;
			}
			
	// ROI input
	if (reuse_roi == false) {
		auto_roi = getBoolean("Use which ROI Method?", "Auto (Otsu)", "Manual");
		
		if (auto_roi == false) {
			setTool("freehand");
			roiManager("reset");
	  		roiManager("show all with labels");
			waitForUser("Add ROIs to the image and the ROI Manager (t) then click OK.");
			}
		}
	
	gaussian_sigma = getNumber("Sigma radius (Gaussian blur):", 0.25);
	
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
	
	// If timecourse, compute an inclusive frame range
	if (timecourse_mode) {
	    frames_to_process = frame_end - frame_start + 1;  // inclusive
	} else {
	    // Non-timecourse: treat as one frame
	    frame_start = frame;
	    frame_end   = frame;
	    frames_to_process = 1;
	}
	
	// Automatically generate masks using the donor (presumably the 'brightest' ch) if auto_roi mode is selected
	if (auto_roi)  {
		roiManager("reset");
		selectImage(donor_id);
		run("Duplicate...", "duplicate title=donor_threshold");
		treshold_id = getImageID();
		selectImage(treshold_id);
		run("Auto Threshold", "method=Otsu white stack");
		for (f = frame_start; f <= frame_end; f++) {
			setSlice(f);
			run("Median", "radius=10"); // Helps smooth out masks produced by poor auto-thresholding.
			run("Create Selection");
			roiManager("add");
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
	donor_df    = newArray(array_length);
	acceptor_df = newArray(array_length);
	fret_df     = newArray(array_length);
	frame_df    = newArray(array_length);
	roi_df 		= newArray(array_length);
	if (timecourse_mode) {
	    time_s = newArray(array_length);
	}

	
	// Processing
	// Fill arrays in "frame-major" order:
	// For each frame (f), all ROI values are stored consecutively:
	// indices [ (f - frame_start)*roi_count  ...  (f - frame_start)*roi_count + (roi_count-1) ]
	
	// Donor
	selectImage(donor_id);
	for (f = frame_start; f <= frame_end; f++) {
		if (auto_roi == false) {
		    // Measure all ROIs, appending to Results table.
		    for (i = 0; i < roi_count; i++) {
		        roiManager("Select", i);
		        Stack.setFrame(f);
		        roiManager("Update");
		        roiManager("Measure");
	    	}
		} else {
			roiManager("Select", f - 1); //1st ROI in manager is index = 0 so need to correct this.
			Stack.setFrame(f);
			roiManager("Update");
			roiManager("Measure");
			}

	    // Copy to flattened arrays using `base` as an offset (i.e. `base` becomes the starting loop iteration instead of 0).
	    base = (f - frame_start) * roi_count;
	    for (i = 0; i < roi_count; i++) {
	        idx = base + i;
	        donor_df[idx] = getResult("Mean", i);
	        frame_df[idx] = f;
	        roi_df[idx] = i + 1;
	        name_df[idx] = img_name;
	
	        if (timecourse_mode) {
	            dt = Stack.getFrameInterval(); // seconds per frame
	            time_s[idx] = (f - 1) * dt; // returns the time (s) that a frame was imaged at. For frame 1, this time is always 0.
	        }
	    }
	    run("Clear Results");
	}
	
	// Acceptor
	selectImage(acceptor_id);
	for (f = frame_start; f <= frame_end; f++) {
		if (auto_roi == false) {
			for (i = 0; i < roi_count; i++) {
		        roiManager("Select", i);
		        Stack.setFrame(f);
		        roiManager("Update");
		        roiManager("Measure");
	    	}
		} else {
			roiManager("Select", f - 1);
			Stack.setFrame(f);
			roiManager("Update");
			roiManager("Measure");
			}
	    
	    base = (f - frame_start) * roi_count;
	    for (i = 0; i < roi_count; i++) {
	        idx = base + i;
	        acceptor_df[idx] = getResult("Mean", i);
	    }
		run("Clear Results");
	}
	
	// FRET
	selectImage(fret_id);
	for (f = frame_start; f <= frame_end; f++) {
		if (auto_roi == false) {
			for (i = 0; i < roi_count; i++) {
		        roiManager("Select", i);
		        Stack.setFrame(f);
		        roiManager("Update");
		        roiManager("Measure");
	    	}
		} else {
			roiManager("Select", f - 1);
			Stack.setFrame(f);
			roiManager("Update");
			roiManager("Measure");
		}
	    
	    base = (f - frame_start) * roi_count;
	    for (i = 0; i < roi_count; i++) {
	        idx = base + i;
	        fret_df[idx] = getResult("Mean", i);
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
			selectImage(fret_calc);
		    roiManager("select", i);        // ROIs are 0-based
		    setSlice(i + 1);                  // slices are 1-based
		    run("Clear Outside", "slice");           // keep only inside the ROI
		    run("Select None");
		}
		fret_calc = getImageID();
	}
	
	
	selectImage(fret_calc); // setBatchMode closes all but the last selected image when toggled back to 'show'.
	run("Duplicate...", "duplicate title=" + img_name + "_FRET_Results");
	setBatchMode("show");
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
	frame = frame_df;
	filename = name_df;
	
	if (timecourse_mode) {
	    Array.show(img_name + "_fl_raw", filename, donor_em, acceptor_em, fret_em, da_ratio, roi, frame, time_s);
	} else {
	    Array.show(img_name + "_fl_raw", filename, donor_em, acceptor_em, fret_em, da_ratio, roi, frame);
	}

  	
}