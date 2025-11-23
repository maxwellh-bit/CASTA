macro "Spatial Coherence" {
	//Must have OrientationJ pulugin installed
	
	//This macro will generate a heatmap of orientation from an image of fibers
	
	//Choose image to be processed and open
	Dialog.create("Choose image");
		Dialog.addFile("Select SHG collagen image to process", "C:/Users/hamilms2/OneDrive - VUMC/Documents/Labs/Vanderbilt/Coffey/Data/Imaging/tumor images/2024.01.05 MC38 tumors 2HG/z-projections/Tumor 887-2.tif");
	Dialog.show();
	path = Dialog.getString();
	open(path);
	name = getTitle();
	run("Smooth");
	run("Smooth");
	run("Smooth");
	
	//Remove background
	//If background is not set to 0, OrientationJ will interpret the edges of the ROI as real edges and introduce edge effects.
	//Make sure to return the maximum before hitting apply.
	run("Brightness/Contrast...");
	waitForUser("Adjust thresholds to remove background, click Apply, then OK.");
	
	//Select ROI
	setTool("polygon");
	waitForUser("Select ROI to analyze and click OK.");
	
	//Sets batch mode to prevent images from displaying and speeds up macro
	setBatchMode(true);
	
	//Make mask to ensure correct area is analyzed
	//removes parts of the image not to be analyzed
	run("Create Mask");
	imageCalculator("AND create", "Mask", name);
	selectWindow(name);
	close();
	selectWindow("Result of Mask");
	
	//Choose location to save heatmap
	Dialog.create("Choose save folder");
		Dialog.addDirectory("Choose folder to save heatmap output", "C:/Users/hamilms2/OneDrive - VUMC/Documents/Labs/Vanderbilt/Coffey/Data/Imaging/tumor images/2024.01.05 MC38 tumors 2HG/z-projections");
	Dialog.show();
	save_path = Dialog.getString();
	
		//Ask the user what they want the resolution of the heatmap to be in pixels.
	Dialog.create("Heatmap resolution"); 
		Dialog.addNumber("Enter the resolution of the heatmap in pixels. Lower resolution (larger pixel size) increases the speed of the program.", 25);
		Dialog.show();
	heat_res = Dialog.getNumber();
	
	//Ask user the diameter of the ROI which will be used to calculate coherence
	//The diameter is the size of the area around a pixel that will be averaged to determine that pixel's value
	Dialog.create("Select box size"); 
		Dialog.addNumber("Enter the diameter of OrientationJ window in pixels. Must be larger than or equal to heatmap resolution.", 200);
		Dialog.show();
	win_size = Dialog.getNumber();

	
	//Get the dimentions of the current image
	height = getHeight();
	width = getWidth();
	
	//Generate a new image to store the heatmap and make it black
	newImage("Coherence Output", "tiff", width, height, "8-bit black");
	setForegroundColor(0, 0, 0);
	makeRectangle(0, 0, height, height);
	run("Fill", "slice");

	//Initialize a variable to store the progress and display progress bar
	prog = 0;
	showProgress(prog);
	
	//Interate through every region in the fiber image
	for (i = 0; i <= height - win_size; i = i + heat_res){
		for (j = 0; j <= width - win_size; j = j + heat_res){
			
			//First, test if the ROI is inside the masked area
			//If the region is outside of the mask, the value of the pixels will be 0
			selectWindow("Mask");
			x_mask = j + (win_size / 2);
			y_mask = i + (win_size / 2);
			mask_boolean = getPixel(x_mask, y_mask);
			
			if (mask_boolean > 0) {
			
				//Draw circle ROI with a diameter of the window size 
				selectWindow("Result of Mask");
				makeOval(j, i, win_size, win_size);
			
				//run orientationJ on the ROI and collect coherence value
				run("OrientationJ Dominant Direction");
				coh = getResult("Coherency [%]", 0);
				coh_float = parseFloat(coh);
				close("Log");
			
				//write the output to a new image encoding coherence
				pix_val = coh_float * 255;
				x_pix = j + (win_size/2) + (heat_res/2);
				y_pix = i + (win_size/2) - (heat_res/2);
				selectWindow("Coherence Output");
				setForegroundColor(pix_val, pix_val, pix_val);
				makeRectangle(x_pix, y_pix, heat_res, heat_res);
				run("Fill", "slice");
			
			}

		}

		//update progress bar
		//This is going to be off if the ROI is much smaller than the image
		prog = i / height;
		showProgress(prog);
		
	}
	
//make the heatmap image an actual heatmap and save
selectWindow("Coherence Output");
run("Fire");
save(save_path + "/" + name + "_coherence map.tif");

}