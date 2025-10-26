// ====================================================
// Advanced Multi-Channel Zeiss CZI Image Processing Macro (Multi-File)
// Author: Clifton Lewis (modified)
// Date modified: 04 Mar 2025 (modified)
// Version: 1.6
//
// Description: This macro processes multiple z‑section CZI images 
// from the same embryo. The user first specifies the number of z‑sections.
// For each file, the macro opens the file and interactively prompts for 
// channel renaming, color assignment, and brightness/contrast adjustments.
// For the first (reference) file, the user then selects a reference channel
// and landmarks, from which transformation parameters are computed and stored.
// These global transformation parameters are then automatically applied 
// to all subsequent files (across all open channel windows). Finally, each
// file's individual channels are saved as PNG and TIFF and a composite TIFF 
// image is created and saved before moving on to the next file.
// =============================================================

// ===== PART A: INITIALISATION AND UTILITY FUNCTIONS =====

// Print start message to the ImageJ Log window.
print(">>> Multi-channel Zeiss Image Processing Macro started:");

// Constants for point coordinates.
X = 0;
Y = 1;

// UI parameters.
padX = 100;       // Horizontal padding for cropping (pixels)
padY = 100;       // Vertical padding for cropping (pixels)
uiScale = 1.0;    // Scale factor (unused for now)

// Global transformation parameters (computed from reference file).
globalTx = 0;
globalTy = 0;
globalAngleDeg = 0;
globalCropX = 0;
globalCropY = 0;
globalCropW = 0;
globalCropH = 0;
globalDoFlip = false;
globalDoFlipHorizontally = false;  // New global flag for horizontal flipping

// Global variable to store the reference channel index.
refIndex = -1;

// Global variable for the output directory.
outputDir = "";

// --- Utility Functions ---

// Function to extract resolution information from image metadata
function getResolutionInfo() {
    image_info = getImageInfo();
    lines = split(image_info, '\n');
    for (i = 0; i < lengthOf(lines); i += 1) {
        if (startsWith(lines[i], "Resolution:")) {
            return lines[i];
        }
    }
    return "Resolution: 1 pixels per unit"; // Default fallback if no resolution is found
}

// Print a point (for debugging).
function printPoint(name, point) {
    print(name + " (" + point[X] + ", " + point[Y] + ")");
}

// Get a point via mouse click (with timeout).
function getPoint(name, point) {
    leftButton = 16;
    getCursorLoc(x, y, z, flags);
    buttonDown = ((flags & leftButton) != 0);
    counter = 0;
    maxCounter = 3000; // ~30 seconds timeout
    while (isOpen("Log") && counter < maxCounter) {
        getCursorLoc(x, y, z, flags);
        if ((flags & leftButton) != 0) {
            if (!buttonDown) {
                buttonDown = true;
                print(name + ": (" + x + ", " + y + ")");
                point[X] = x;
                point[Y] = y;
                return;
            }
        } else if (buttonDown) {
            buttonDown = false;
        }
        wait(10);
        counter += 1;
    }
    if (counter >= maxCounter) {
        exit("Error: Timeout waiting for point selection for " + name);
    }
}

// Update a point's coordinates after translation and rotation.
function update(point, angle, tx, ty) {
    newPoint = newArray(2);
    newPoint[X] = point[X] + tx;
    newPoint[Y] = point[Y] + ty;
    newPoint[X] -= getWidth() / 2;
    newPoint[Y] -= getHeight() / 2;
    x_val = newPoint[X];
    y_val = newPoint[Y];
    // Reverse angle for rotation correction.
    angle = -angle;
    newPoint[X] = x_val * Math.cos(angle) - y_val * Math.sin(angle);
    newPoint[Y] = x_val * Math.sin(angle) + y_val * Math.cos(angle);
    newPoint[X] += getWidth() / 2;
    newPoint[Y] += getHeight() / 2;
    point[X] = newPoint[X];
    point[Y] = newPoint[Y];
}

// (Optional) Add a temporary overlay for a point.
function addPointOverlay(point, color) {
    radius = 10 * uiScale;
    makeOval(point[X] - radius, point[Y] - radius, 2 * radius, 2 * radius);
    setColor(color);
    run("Add Selection...");
}

// --- Increase Canvas Size ---
// Increases the canvas size by adding half the diagonal to the dimensions.
function increaseCanvasSize() {
    w = getWidth();
    h = getHeight();
    diag = Math.sqrt(w * w + h * h);
    newWidth = w + diag/2;
    newHeight = h + diag/2;
    // Increase the canvas size, centering the original image.
    run("Canvas Size...", "width=" + newWidth + " height=" + newHeight + " position=Center");
    print("Canvas size increased: new width = " + newWidth + ", new height = " + newHeight);
}

// --- mergeChannels ---
// This function dynamically builds the merge parameters based on available channel windows,
// then runs the "Merge Channels..." command.
function mergeChannels() {
    mergeParams = "";
    titles = getList("image.titles");
    for (c = 1; c <= channels; c++) {
        found = false;
        for (j = 0; j < lengthOf(titles); j++) {
            if (startsWith(titles[j], "C" + c)) {
                mergeParams += " c" + c + "=[" + titles[j] + "]";
                found = true;
                break;
            }
        }
        if (!found) {
            print("Warning: Channel " + c + " window not found.");
        }
    }
    if (mergeParams != "") {
        run("Merge Channels...", mergeParams + " create");
    } else {
        exit("Error: No channel windows found to merge.");
    }
}

// --- getBaseName ---
// Extracts the base file name (without path and extension) from the given file path.
function getBaseName(path) {
    slashIndex = lastIndexOf(path, "/");
    backslashIndex = lastIndexOf(path, "\\");
    if (backslashIndex > slashIndex) {
        slashIndex = backslashIndex;
    }
    base = substring(path, slashIndex + 1);
    dotIndex = lastIndexOf(base, ".");
    if (dotIndex > 0) {
        base = substring(base, 0, dotIndex);
    }
    return base;
}

// --- Saving Functions ---
// These functions save each channel using the assigned channel name and the file base name,
// and save the merged image. The output directory is prepended to each file name.
function savePNGChannels(baseName, channelNames) {
    for (c = 0; c < channels; c++) {
        // Set the active channel in the current hyperstack.
        Stack.setChannel(c + 1);
        saveAs("PNG", outputDir + baseName + "_" + channelNames[c] + ".png");
    }
}

function saveTIFFChannels(baseName, channelNames) {
    for (c = 0; c < channels; c++) {
        Stack.setChannel(c + 1);
        saveAs("TIFF", outputDir + baseName + "_" + channelNames[c] + ".tif");
    }
}

function savePNGMerge(baseName) {
    // Attempt to select the composite window (adjust the name if needed)
    if (isOpen("Composite")) {
        selectWindow("Composite");
    } else {
        print("Warning: 'Composite' window not found; saving active window instead.");
    }
    saveAs("PNG", outputDir + baseName + "_Merged.png");
}

function saveTIFFMerge(baseName) {
    if (isOpen("Composite")) {
        selectWindow("Composite");
    } else {
        print("Warning: 'Composite' window not found; saving active window instead.");
    }
    saveAs("TIFF", outputDir + baseName + "_Merged.tif");
}

// --- closeSplitChannels ---
// Closes all open image windows whose title starts with "C".
function closeSplitChannels() {
    titles = getList("image.titles");
    for (j = 0; j < lengthOf(titles); j++) {
        if (startsWith(titles[j], "C")) {
            close(titles[j]);
        }
    }
}

// --- Transformation Functions ---

// applyTransform(): For the reference file only.
// Computes slopes, center, rotation angle, and translation offsets based on selected landmarks,
// updates landmarks, defines the crop region, and stores the transformation parameters globally.
// Then, applies the same transformation to all open channel windows.
// At the end, the scale bar (line only, no text) is drawn and the image is ensured to have the head on the left.
function applyTransform() {
    if (tail[X] == head[X]) {
        exit("Error: Head and Tail landmarks have identical X coordinates.");
    }
    m1 = (tail[Y] - head[Y]) / (tail[X] - head[X]);
    
    if (dors[X] == vent[X]) {
        exit("Error: Dorsal and Ventral landmarks have identical X coordinates.");
    }
    m2 = (dors[Y] - vent[Y]) / (dors[X] - vent[X]);

    c1 = head[Y] - m1 * head[X];
    c2 = vent[Y] - m2 * vent[X];

    printPoint("m1, c1", newArray(m1, c1));
    printPoint("m2, c2", newArray(m2, c2));

    xCenter = (c1 - c2) / (m2 - m1);
    yCenter = m1 * xCenter + c1;
    print("Center: (" + xCenter + ", " + yCenter + ")");

    angle = PI + Math.atan2(head[Y] - yCenter, head[X] - xCenter);
    angle_deg = angle * 180 / PI;
    print("Angle: " + angle_deg);

    tx = (getWidth() / 2) - xCenter;
    ty = (getHeight() / 2) - yCenter;

    // Store transformation parameters globally.
    globalTx = tx;
    globalTy = ty;
    globalAngleDeg = angle_deg;

    // Update landmark positions.
    update(head, angle, tx, ty);
    update(tail, angle, tx, ty);
    update(dors, angle, tx, ty);
    update(vent, angle, tx, ty);

    top = Math.max(dors[Y], vent[Y]);
    bottom = Math.min(dors[Y], vent[Y]);
    right = Math.max(head[X], tail[X]);
    left = Math.min(head[X], tail[X]);
    rectX = left - padX;
    rectY = bottom - padY;
    rectW = 2 * padX + (right - left);
    rectH = 2 * padY + (top - bottom);

    globalCropX = rectX;
    globalCropY = rectY;
    globalCropW = rectW;
    globalCropH = rectH;

    // Apply the computed transformation to all channel windows.
    titles = getList("image.titles");
    for (k = 0; k < lengthOf(titles); k++) {
        if (startsWith(titles[k], "C")) {  // Assume channel windows start with "C"
            selectWindow(titles[k]);
            // Apply translation.
            makeRectangle(0, 0, getWidth(), getHeight());
            run("Translate...", "x=" + tx + " y=" + ty + " interpolation=None");
            // Apply rotation.
            run("Rotate...", "angle=" + (-angle_deg) + " grid=1 interpolation=Bilinear");
            // Apply cropping.
            makeRectangle(rectX, rectY, rectW, rectH);
            setColor(255, 0, 0);
            setLineWidth(20);
            run("Crop");
            // Apply vertical flip if needed.
            if (vent[Y] > dors[Y]) {
                run("Flip Vertically");
            }
            // Ensure head is on the left side.
            // If head's X coordinate is greater than tail's, then head is on the right.
            if (head[X] > tail[X]) {
                run("Flip Horizontally");
                globalDoFlipHorizontally = true;
            } else {
                globalDoFlipHorizontally = false;
            }
            // Draw scale bar as a line (no text).
            //run("Scale Bar...", "width=100 height=100 thickness=70 font=0 bold overlay");
        }
    }
}

// applyStoredTransform(): For subsequent files.
// Applies the stored transformation parameters to all open channel windows.
// At the end, the scale bar (line only, no text) is drawn and the stored horizontal flip is applied.
function applyStoredTransform() {
    titles = getList("image.titles");
    for (k = 0; k < lengthOf(titles); k++) {
        if (startsWith(titles[k], "C")) {
            selectWindow(titles[k]);
            makeRectangle(0, 0, getWidth(), getHeight());
            run("Translate...", "x=" + globalTx + " y=" + globalTy + " interpolation=None");
            run("Rotate...", "angle=" + (-globalAngleDeg) + " grid=1 interpolation=Bilinear");
            makeRectangle(globalCropX, globalCropY, globalCropW, globalCropH);
            setColor(255, 0, 0);
            setLineWidth(20);
            run("Crop");
            if (globalDoFlip) {
                run("Flip Vertically");
            }
            // Apply stored horizontal flip if needed.
            if (globalDoFlipHorizontally) {
                run("Flip Horizontally");
            }
            // Draw scale bar as a line (no text).
            //run("Scale Bar...", "width=100 height=100 thickness=70 font=0 bold overlay");
        }
    }
}

// ===== PART B: FILE SELECTION AND PROCESSING =====

// Ask the user for the number of input z‑section files.
Dialog.create("Input Z‑Sections");
Dialog.addNumber("Enter the number of z‑section files to process:", 1);
Dialog.show();
numZ = Dialog.getNumber();
print("Number of input z‑sections: " + numZ);

// Ask the user to choose an output directory for saving images.
outputDir = getDirectory("Choose output directory for saving images:");
// Ensure outputDir ends with a slash.
if (endsWith(outputDir, "/") == false && endsWith(outputDir, "\\") == false) {
    outputDir = outputDir + "/";
}
print("Output directory: " + outputDir);

// Process each file one by one.
for (i = 0; i < numZ; i++) {
    titleMsg = "Select CZI file for z‑section " + (i + 1);
    filePath = File.openDialog(titleMsg);
    if (filePath == "") {
        exit("File selection canceled. Macro terminated.");
    }
    print("\nProcessing file: " + filePath);
    
    // Get the base file name (without path or extension)
    baseName = getBaseName(filePath);

    // Open the file.
    run("Bio-Formats Importer", "open=[" + filePath + "] autoscale=false color_mode=Default view=Hyperstack stack_order=XYCZT");
    getDimensions(w, h, channels, slices, frames);

    // Extract and set scale information
    resolution_metadata = getResolutionInfo();
    resolution_tokens = split(resolution_metadata, " ");
    scale_value = parseFloat(resolution_tokens[1]);
    print("Scale value: " + scale_value + " pixels/um");
    run("Set Scale...", "distance=" + toString(scale_value) + " known=1 unit=um");

    // --- CHANNEL NAMING & COLOR ASSIGNMENT (Interactive for EACH file) ---
    defaultNames = newArray();
    for (c = 1; c <= channels; c++) {
        defaultNames[c - 1] = "Channel " + c;
    }
    listMsg = "Detected Channels:\n";
    for (c = 0; c < channels; c++) {
        listMsg += "- " + defaultNames[c] + "\n";
    }
    showMessage("Channel Names", listMsg);

    Dialog.create("Rename Channels (File " + (i + 1) + ")");
    for (c = 0; c < channels; c++) {
        Dialog.addString("New name for Channel " + (c + 1) + " (" + defaultNames[c] + "):", defaultNames[c]);
    }
    Dialog.show();
    localChannelNames = newArray(channels);
    for (c = 0; c < channels; c++) {
        nameInput = Dialog.getString();
        if (trim(nameInput) == "") {
            localChannelNames[c] = defaultNames[c];
        } else {
            localChannelNames[c] = nameInput;
        }
    }

    Dialog.create("Assign Channel Colors (File " + (i + 1) + ")");
    for (c = 0; c < channels; c++) {
        Dialog.addString("Color for \"" + localChannelNames[c] + "\" (hex code or LUT name):", "");
    }
    Dialog.show();
    localChannelColors = newArray(channels);
    for (c = 0; c < channels; c++) {
        localChannelColors[c] = trim(Dialog.getString());
    }

    // --- BRIGHTNESS/CONTRAST ADJUSTMENT (Interactive for EACH file) ---
    Dialog.create("Channel Brightness/Contrast Adjustment (File " + (i + 1) + ")");
    for (c = 0; c < channels; c++) {
        // Ensure the correct channel is active.
        Stack.setChannel(c + 1);
        run("Brightness/Contrast...");
        waitForUser("Adjust Channel " + (c + 1), "1. Adjust brightness/contrast for " + localChannelNames[c] +
                    ".\n2. Click 'Apply' in B&C window.\n3. Click 'OK' here.");
        selectWindow("B&C");
        run("Close");
    }

    // --- REFERENCE FILE PROCESSING ---
    if (i == 0) {
        Dialog.create("Select Reference Channel (Reference File)");
        Dialog.addChoice("Reference Channel:", localChannelNames, localChannelNames[0]);
        Dialog.show();
        refChoice = Dialog.getChoice();
        refIndex = -1;
        for (j = 0; j < lengthOf(localChannelNames); j++) {
            if (localChannelNames[j] == refChoice) {
                refIndex = j;
                break;
            }
        }
        if (refIndex == -1) {
            exit("Error: Reference channel not found.");
        }
        print("Selected reference channel: " + refChoice + " (index " + refIndex + ")");

        // Increase canvas size BEFORE splitting channels.
        increaseCanvasSize();

        // Split channels.
        run("Split Channels");
        expectedPrefix = "C" + (refIndex + 1);
        titles = getList("image.titles");
        refWin = "";
        for (j = 0; j < lengthOf(titles); j++) {
            if (startsWith(titles[j], expectedPrefix)) {
                refWin = titles[j];
                break;
            }
        }
        if (refWin == "") {
            print("Available windows:");
            for (j = 0; j < lengthOf(titles); j++) {
                print(titles[j]);
            }
            exit("Error: Reference channel window not found. Expected prefix: " + expectedPrefix);
        }
        selectWindow(refWin);

        print("\nSelect points in this order: HEAD, VENTRAL, TAIL, DORSAL");
        waitForUser("Point Selection", "Please ensure all channels are visible and proceed to select points:\n" +
                    "1. Click HEAD position\n" +
                    "2. Click VENTRAL (belly) position\n" +
                    "3. Click TAIL position\n" +
                    "4. Click DORS (back) position\n\nClick OK to start selection");

        // Landmark selection for transformation.
        head = newArray(2);   getPoint(" - HEAD", head);
        vent = newArray(2);   getPoint(" - VENTRAL", vent);
        tail = newArray(2);   getPoint(" - TAIL", tail);
        dors = newArray(2);   getPoint(" - DORS", dors);

        // Compute and apply transformation.
        applyTransform();
        // Wait for user to perform any additional modifications.
        waitForUser("Additional Modifications", "Perform any additional modifications now. Click OK to continue and save the images.");
        
        // Close all split channel windows.
        closeSplitChannels();
    } else {
        // --- SUBSEQUENT FILE PROCESSING ---
        Dialog.create("Apply Transformation?");
        Dialog.addChoice("Apply stored transformation to this file?", newArray("Yes", "No"), "Yes");
        Dialog.show();
        applyChoice = Dialog.getChoice();
        if (applyChoice == "Yes") {
            run("Split Channels");
            expectedPrefix = "C" + (refIndex + 1);
            titles = getList("image.titles");
            refWin = "";
            for (j = 0; j < lengthOf(titles); j++) {
                if (startsWith(titles[j], expectedPrefix)) {
                    refWin = titles[j];
                    break;
                }
            }
            if (refWin == "") {
                print("Reference channel window not found in file: " + filePath);
            } else {
                selectWindow(refWin);
                applyStoredTransform();
            }
            // Wait for user to perform any additional modifications.
            waitForUser("Additional Modifications", "Perform any additional modifications now. Click OK to continue and save the images.");
            
            // Close all split channel windows.
            closeSplitChannels();
        } else {
            print("User opted not to apply transformation for file: " + filePath);
        }
    }
    
    // Print message for successful processing of current file.
    print("Successfully completed processing file " + (i + 1) + " of " + numZ + ".");
}

// End of macro.
print("End of Macro: All files processed successfully.");
