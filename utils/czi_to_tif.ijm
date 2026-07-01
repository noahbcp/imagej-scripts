/* czi_to_tif.ijm
 * Author: Noah B.C. Piper
 * SPDX-License-Identifier: BSD-3-Clause
 * ==============================================================================
 * Adds a text label to a range of frames in a hyperstack. Useful for showing
 * when a condition begins and ends (e.g. a drug addition).
 * Does not 'burn in' the overlay.
 * ==============================================================================
 * Input Requirements: Directory with .czi files.
 * Outputs: .tif files written to the input directory as `filename.tif`.
 * ==============================================================================
 */

macro "CZI to TIF" {
    dir = getDirectory("Choose folder with .czi files");
    list = getFileList(dir);

    setBatchMode(true);

    for (i = 0; i < list.length; i++) {

        if (endsWith(list[i], ".czi")) {

            inputPath = dir + list[i];
            outputPath = dir + replace(list[i], ".czi", ".tif");

            print("Converting: " + inputPath);

            // Import using Bio-Formats
            run("Bio-Formats Importer", "open=[" + inputPath + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");

            // Save as TIFF
            saveAs("Tiff", outputPath);

            close("*"); // Close all images to free memory
        }
    }

    setBatchMode(false);
    print("Done.");
}
