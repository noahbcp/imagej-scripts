/* spectral_plot.ijm
 * Author: Noah B.C. Piper
 * SPDX-License-Identifier: BSD-3-Clause
 * ==============================================================================
 * Generates a 'spectra' from Zeiss Lambda (or similar) spectral scans.
 * ==============================================================================
 */

macro "Plot Spectra" {
	// Close if no images are open.
	if (nImages == 0) {
		exit;
	}
	// User input.
	Dialog.create("Plot Spectra");
	Dialog.addNumber("Spectral width (nm)", 1);
	Dialog.addNumber("Range starting wavelength (nm)", 450);
	Dialog.addCheckbox("Display intensity as % maximum", true);
	Dialog.show();
	specWidth = Dialog.getNumber();
	specStart = Dialog.getNumber();
	normaliseMode = Dialog.getCheckbox();

	// Assume that `channels` contains spectral scan data.
	Stack.getDimensions(width, height, channels, slices, frames);

	// Initialise storage variables
	intensityArray = newArray(channels);
	wlArray = newArray(channels);
	maxIntensity = 0;

	run("Clear Results");
	for (i = 0; i < channels; i++) {
		Stack.setChannel(i + 1); // Channels are base-1.
		run("Measure");
		meanIntensity = getResult("Mean", i);
		intensityArray[i] = meanIntensity;
		if (meanIntensity > maxIntensity) {
			maxIntensity = meanIntensity;
		}
		wlArray[i] = specStart;
		specStart = specStart + specWidth;
	}
	if (normaliseMode) {
		for (i = 0; i < intensityArray.length; i++) {
			intensityArray[i] = intensityArray[i] / maxIntensity;	
		}	
	}
	Plot.create(getTitle() + " Spectra", "Wavelength (nm)", "Intensity", wlArray, intensityArray);
	run("Clear Results");
}
