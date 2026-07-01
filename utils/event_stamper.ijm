/* event_stamper.ijm
 * Author: Noah B.C. Piper
 * SPDX-License-Identifier: BSD-3-Clause
 * ==============================================================================
 * Adds a text label to a range of frames in a hyperstack. Useful for showing
 * when a condition begins and ends (e.g. a drug addition).
 * Does not 'burn in' the overlay.
 * ==============================================================================
 * Input Requirements: Image stack with time dimension.
 * Outputs: Image stack with overlay over specified frames.
 * ==============================================================================
 */

macro "Event Stamper" {
    // Ensure there's an image open
    if (nImages == 0) {
        showMessage("No image open", "Please open a hyperstack/stack with time frames first.");
        exit();
    }

    // Capture dimensions and confirm frames exist
    getDimensions(w, h, nC, nZ, nT);
    isHyper = Stack.isHyperstack;
    if (nT < 1) {
        showMessage("No time frames", "This image has no time dimension (nT=0).");
        exit();
    }

    // Defaults
    defaultStart = 1;
    defaultEnd = nT;
    defaultText = "Sample label";
    defaultSize = 18;
    defaultMargin = 30;
    defaultColor = "white"; // can be white/black/red/green/blue/yellow/cyan/magenta
    defaultBg = "none"; // none or black (simple backdrop)

    // Dialog
    Dialog.create("Add text on a frame range");
    Dialog.addNumber("Start frame (T):", defaultStart);
    Dialog.addNumber("End frame (T):", defaultEnd);
    Dialog.addString("Text:", defaultText, 30);
    Dialog.addNumber("Font size (px):", defaultSize);
    Dialog.addNumber("Margin (px):", defaultMargin);
    Dialog.addChoice("Text color:", newArray("white", "black", "red", "green", "blue", "yellow", "cyan", "magenta"), defaultColor);
    Dialog.addChoice("Background:", newArray("none", "black"), defaultBg);
    Dialog.show();

    tStart = Dialog.getNumber();
    tEnd = Dialog.getNumber();
    label = Dialog.getString();
    fSize = Dialog.getNumber();
    margin = Dialog.getNumber();
    tColor = Dialog.getChoice();
    bg = Dialog.getChoice();

    // Sanity checks and clamping
    if (tStart > tEnd) {
        tmp = tStart;
        tStart = tEnd;
        tEnd = tmp;
    }
    if (tStart < 1) tStart = 1;
    if (tEnd > nT) tEnd = nT;

    if (lengthOf(label) == 0) {
        showMessage("Empty text", "Please enter a non-empty text string.");
        exit();
    }

    // Prepare overlay (preserve any existing)
    run("Show Overlay");
    // We'll add one text ROI per frame so it only shows at that frame

    // Set font parameters for the TextRoi
    setFont("SansSerif", fSize, "antialiased");

    // Compute anchor position (lower-left with margin)
    x = margin;
    // In ImageJ, the y given to makeText() is the baseline of the text.
    // We push it up by ~margin from the bottom; small extra offset for safety.
    y = h - margin;

    // Remember current positions to restore later
    Stack.getPosition(channel, slice, frame);
    origC = channel;
    origZ = slice;
    origT = frame;

    // Loop through frames and add text at each frame's overlay position
    for (t = tStart; t <= tEnd; t++) {
        Stack.setFrame(t);
        // If hyperstack, keep current C/Z (or set to 1 if you want only C1/Z1)
        // Here we keep the current viewer indices:
        if (isHyper) Stack.setPosition(origC, origZ, t);

        // Create text selection and push to overlay
        makeText(label, x, y);

        // Apply simple visibility styling
        // Text color
        setColor(tColor);

        // Optional simple background: draw a small opaque box behind text by duplicating the ROI as a filled rectangle.
        // Because macro language doesn't expose font metrics, we use a quick-and-safe backdrop width.
        if (bg == "black") {
            // Duplicate the text ROI bounds to create a rectangle backdrop
            // Get bounding box of the current selection
            getSelectionBounds(bx, by, bw, bh);
            // Slight padding
            pad = round(fSize * 0.35);
            // Save current color then set fill to black
            saveCol = getValue("overlay.foreground"); // store but won't reapply directly
            setColor("black");
            // Make and add filled rectangle behind text for contrast (assigned to this frame)
            makeRectangle(maxOf(0, bx - pad), maxOf(0, by - pad), bw + 2 * pad, bh + 2 * pad);
            run("Add Selection..."); // adds rectangle at this frame
            // Recreate the text selection (the rectangle replaced it)
            makeText(label, x, y);
            setColor(tColor);
        }

        // Add the text selection to the overlay at this frame
        run("Add Selection...");
        // Clear active selection so next loop is clean
        run("Select None");
    }

    // Restore original position
    if (isHyper)
        Stack.setPosition(origC, origZ, origT);
    else
        Stack.setFrame(origT);

    // Ensure overlay is visible
    run("Show Overlay");
    // Optional: keep selections editable in overlay
    // i.e. run("Overlay Options...", "stroke=1 width=1 draw");
}

// Small helpers for safety
function maxOf(a, b) {
    if (a > b) return a;
    else return b;
}
