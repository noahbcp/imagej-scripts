/* rotate_align.ijm
 * Author: Noah B.C. Piper
 * SPDX-License-Identifier: BSD-3-Clause
 * ==============================================================================
 * Rotates an image to the angle of a straight line selection.
 * Useful for aligning slightly skewed images (e.g. western blots).
 * ==============================================================================
 * Input Requirements: N/A
 * Outputs: Rotated image.
 * ==============================================================================
 */

macro "Rotate and Align" {

    getLine(x1, y1, x2, y2, width);

    // Horizontal & vertical displacement
    dy = y1 - y2;
    dx = x2 - x1;

    // Length & midpoint of original line
    length = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    cx = (x1 + x2) / 2;
    cy = (y1 + y2) / 2;

    delta = atan2(dy, dx) * 180.0 / PI;

    run("Arbitrarily...", "angle=" + delta + " interpolate");

    // Redraw line
    makeLine(cx - length / 2, cy, cx + length / 2, cy);
}
