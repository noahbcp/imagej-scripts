/* plot_intensity_time
 * Author: Noah B.C. Piper
 * SPDX-License-Identifier: BSD-3-Clause
 * ==============================================================================
 * Generates an intensity-time plot from an ROI or the entire image.
 * ==============================================================================
 */

macro "Plot Intensity Over Time" {
    getDimensions(width, height, channels, slices, frames);
    Stack.getPosition(channel, slice, frame);
    if (frames < 2) {
        exit("Image has no time-dimension.");
    }
    if (nResults > 0) {
        showMessageWithCancel("Plot Intensity Over Time", "Results will be cleared!");
    }
    run("Clear Results");
    intensityArray = newArray(frames);
    timeArray = newArray(frames);
    if (Stack.getFrameInterval() > 0) {
        dt = Stack.getFrameInterval();
        xTitle = "Time (s)";
    } else {
        dt = 1; // Fallback if instrument hasn't written frameInterval to file.
        xTitle = "Frame";
    }
    for (i = 0; i < frames; i++) {
        Stack.setFrame(i + 1); // Frames are base-1
        Stack.setChannel(channel);
        run("Measure");
        intensityArray[i] = getResult("Mean", i);
        timeArray[i] = i * dt;
    }
    Plot.create(getTitle() + " Intensity Plot", xTitle, "Intensity", timeArray, intensityArray);
    run("Clear Results");
}
