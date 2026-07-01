/* batch_transform.ijm
 * Author: Noah B.C. Piper
 * SPDX-License-Identifier: BSD-3-Clause
 * ==============================================================================
 * Applies a geometric transformation (rotate or flip) to all .tif images in a folder.
 * Overwrites existing image files.
 * ==============================================================================
 * Input Requirements: Directory with .tif files.
 * Outputs: Transformed images saved in place (overwrites originals).
 * ==============================================================================
 */

macro "Batch Transform" {
	inputDir = getDirectory("Choose a Directory");
	list = getFileList(inputDir);
	options = newArray("Rotate 90 Degrees Right", "Rotate 90 Degrees Left", "Flip Horizontally", "Flip Vertically");

	Dialog.create("Batch Transformation");
	Dialog.addChoice("Select transformation:", options, "Rotate 90 Degrees Right");
	Dialog.show();

	method = Dialog.getChoice();

	for (i = 0; i < list.length; i++) {
		if (endsWith(list[i], ".tif")) {
			open(inputDir + list[i]);
			run(method);
			run("Save");
			close("*");
		}
	}
}
