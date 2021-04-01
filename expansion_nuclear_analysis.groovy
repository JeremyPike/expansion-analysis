import java.util.ArrayList
import java.awt.Rectangle


import ij.IJ
import ij.plugin.Duplicator
import ij.gui.GenericDialog
import ij.gui.WaitForUserDialog
import ij.WindowManager
import ij.measure.ResultsTable

import fiji.plugin.trackmate.Model
import fiji.plugin.trackmate.Settings
import fiji.plugin.trackmate.detection.LogDetectorFactory
import fiji.plugin.trackmate.features.FeatureFilter
import fiji.plugin.trackmate.TrackMate
import fiji.plugin.trackmate.SelectionModel
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer
import fiji.plugin.trackmate.features.spot.SpotIntensityAnalyzerFactory
import fiji.plugin.trackmate.SpotCollection
import fiji.plugin.trackmate.Spot
import fiji.plugin.trackmate.tracking.sparselap.SimpleSparseLAPTrackerFactory
import fiji.plugin.trackmate.SelectionModel
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer
import fiji.plugin.trackmate.FeatureModel;
import fiji.plugin.trackmate.Spot
import fiji.plugin.trackmate.visualization.SpotColorGenerator

import loci.plugins.BF


/* Create user interface requesting necessary user intputs */

#@File(label = "Input directory", style = "directory") srcFile
#@String(label = "File extension", value=".ims") ext
#@Integer(label="Satellite spot channel", value=1) satChannel
#@Integer(label="Central spot channel", value=2) centralChannel
#@Double(label="Specify spot radius for satellite channel (microns)", value=0.75) spotRadiusC1
#@Double(label="Specify spot radius for central channel (microns)", value=0.75) spotRadiusC2
#@Double(label="Specify quality threshold for spots (satellite)", value=6.) spotThreshQualC1
#@Double(label="Specify quality threshold for spots (central)", value=7.5) spotThreshQualC2
#@Double(label="Specify rolling ball background subtration radius (satellite) (pixels)", value=40) ballSize1
#@Double(label="Specify rolling ball background subtration radius (central) (pixels)", value=14) ballSize2
#@Double(label="XY Voxel calibration (microns)", value=0.159) xyCal
#@Double(label="Z Voxel calibration (microns)", value=0.500) zCal
#@Double(label="Linking distance", value=0) linkingDist
#@Integer(label="Box crop width (microns)", value=10) cropWidth
#@Integer(label="Number of classes", value=5) numClasses
#@Boolean(label="Do you want tp save cropped site images", value=false) saveCrops


// close all open window
IJ.run("Close All", "")
// create results table to hold summary statistics
output = new ResultsTable()

// create empty TrackMate Settings object
settings = new Settings()
// use a Laplacian of Gaussian spot detector
settings.detectorFactory = new LogDetectorFactory()
// specify that intensity based features should be calculated for spots
settings.addSpotAnalyzerFactory(new SpotIntensityAnalyzerFactory<>())
// set detector settings as default for now
settings.detectorSettings = settings.detectorFactory.getDefaultSettings()

// process all data in folder with specified file extension
processFolder()
// show the output table
output.show("Results")
// save results table to csv file
outPath =  srcFile.getAbsolutePath() + File.separator + "outputStats.csv"
output.save(outPath)

def processFolder() {
	srcFile.eachFileRecurse {
		name = it.getName()
		if (name.endsWith(ext)) 
			processFile(it)
	}
}

def processFile(file) {
	
	println "Processing $file"

	// Opening the file using bioformats
	imp = BF.openImagePlus(file.getAbsolutePath())[0]
	imp.setDisplayMode(IJ.COMPOSITE)
	
	// Overwite any existing calibration with user specified values
	cali = imp.getCalibration()
	cali.pixelWidth = xyCal
	cali.pixelHeight = xyCal
	cali.pixelDepth = zCal
	cali.setUnit('micron')

	// Duplicate the satellite (C1) and central spot channels (C2)
	impC1 = new Duplicator().run(imp, satChannel, satChannel, 1, imp.nSlices, 1, 1)
	impC2 = new Duplicator().run(imp, centralChannel, centralChannel, 1, imp.nSlices, 1, 1)

	// Perform rolling ball background subtraction
	IJ.run(impC1, "Subtract Background...", "rolling=" + ballSize1 + " stack")
	IJ.run(impC2, "Subtract Background...", "rolling=" + ballSize2 + " stack")

	// TrackMate spot detections settings for satellite channel
	settings.setFrom(impC1)
	settings.detectorSettings['RADIUS'] = spotRadiusC1
	settings.detectorSettings['THRESHOLD'] = spotThreshQualC1
	settings.detectorSettings['TARGET_CHANNEL'] = 1
	// Detect spots using Laplacian of Gaussian filter
	modelC1 = detectSpots(settings)
	spotsC1 = modelC1.getSpots()

	// TrackMate spot detections settings for central channel
	settings.setFrom(impC2)
	settings.detectorSettings['RADIUS'] = spotRadiusC2
	settings.detectorSettings['THRESHOLD'] = spotThreshQualC2
	settings.detectorSettings['TARGET_CHANNEL'] = 1
	// Detect spots using Laplacian of Gaussian filter
	modelC2 = detectSpots(settings)
	spotsC2 = modelC2.getSpots()

	// Create tables of spot statstics for both channels and save as csv files
	spotTableC1 = calcAllSpotFeatures(modelC1)
	outSpotPathC1 =  file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf("."))  + "_spotStatsC1.csv"
	spotTableC1.save(outSpotPathC1)
	spotTableC2 = calcAllSpotFeatures(modelC2)
	outSpotPathC2 =  file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf("."))  + "_spotStatsC2.csv"
	spotTableC2.save(outSpotPathC2)

	// Display the detected spots over the image stack
	displaySpots(spotsC1, spotsC2, imp)
	
	// Ask the user if they want to do class assigments
	ok = askForBoolean("Do you want to do class assignments for this dataset?")
	if(ok) {

		// For convenience
		IJ.run("Brightness/Contrast...")
		
		// Construct graph by linking all "central" spots within linkingDist
		graph = constructGraph(spotsC2, linkingDist)
		
		// Find the connected components of the graph
		ccs = connectedComponents(graph)

		// Create table to hold site statistics
		siteTable = new ResultsTable()

		// Iterate through sites (connected components of the graph)
		count = 0
		for (cc in ccs) {
	
			// Find centre of mass for the site (siteCoords)
			siteCoords = new double[3]
			ccNumSpots = cc.getValue().size()
			for (spot in cc.getValue()) {
				for (d = 0; d < 3; d++) 
					siteCoords[d] += spot.getDoublePosition(d)		
			}
			for (d = 0; d < 3; d++) 
				siteCoords[d] = siteCoords[d] / ccNumSpots

			// Crop the stack around around the site centre of mass
			// Width of crop is a user specified parameter
			xMin = (int) Math.floor((siteCoords[0] - cropWidth / 2) / xyCal)
			yMin = (int) Math.floor((siteCoords[1] - cropWidth / 2) / xyCal)
			zMin = (int) Math.floor((siteCoords[2] - cropWidth / 2) / zCal) + 1
			zMax = (int) Math.floor((siteCoords[2] + cropWidth / 2) / zCal) + 1
			xyCropWidth = (int) Math.floor(cropWidth / xyCal
			imp.setRoi(xMin, yMin, xyCropWidth, xyCropWidth)
			impCrop = new Duplicator().run(imp, 1, 2, zMin, zMax, 1, 1)

			// Display cropped site stack
			impCrop.setDisplayMode(IJ.COMPOSITE)
			// Add spot representing the site
			displaySpot = new Spot(cropWidth / 2, cropWidth / 2, cropWidth / 2, 2, 0)
			siteSpot = new SpotCollection()
			siteSpot.add(displaySpot, 0);
			displaySpots(siteSpot, impCrop)

			// Enhance constrast in relvent channels, resize windows etc
			// For convenience
			impCrop.setPosition(1, (int) Math.ceil(impCrop.nSlices / 2), 1)
			impCrop.setC(satChannel)
			IJ.run(impCrop, "Enhance Contrast", "saturated=0.35")
			impCrop.setC(centralChannel)
			IJ.run(impCrop, "Enhance Contrast", "saturated=0.35")
			wind = WindowManager.getCurrentWindow()
			wind.setLocationAndSize((int) wind.LEFT_ALIGNMENT, (int) wind.TOP_ALIGNMENT, wind.width * 5, wind.height * 5)

			// Ask user to define a class for the current site
			cat = askForClass(numClasses)
			
			// Add site ID, COM and class to the site stats output table
			siteTable.incrementCounter()
			siteTable.addValue("ID", count)
			siteTable.addValue("POSITION_X", siteCoords[0])
			siteTable.addValue("POSITION_Y", siteCoords[1])
			siteTable.addValue("POSITION_Z", siteCoords[2])
			siteTable.addValue("Class", cat)

			// If requested save the cropped site stack
			if (saveCrops) {
				outCropPath = file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf("."))  + "_class" + cat + "_siteID" + count + ".tif"
				IJ.saveAs(impCrop, "Tiff", outCropPath);
			}
	
			// Close cropped image		
			impCrop.close()
			count++
		}
	
		// Save the site class table
		outSitePath =  file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf("."))  + "_siteStats.csv"
		siteTable.save(outSitePath)
	
		// Add summary measurments to overall output table
		output.incrementCounter()
		output.addValue("File", file.getAbsolutePath())
		output.addValue("Num spots C1", spotsC1.getNSpots(true))
		output.addValue("Num spots C2", spotsC2.getNSpots(true))
		output.addValue("Num sites", ccs.size())
	}	

	// Close image
	imp.close()
	
}

/**
 * Detect spots using TrackMate
 * 
 * @param settings
 *            the Settings for TrackMate
 * 
 * @return the TrackMate model
 * 
 */
def detectSpots(settings) {
	
	// Create empty TrackMate Model
	model = new Model()
	// Create TrackMate object to peform spot detection
    trackmate = new TrackMate(model, settings)
  
    ok = trackmate.checkInput()
    if (ok == false) 
    	println trackmate.getErrorMessage()
    	
    // Detect spots and execute filtering
    trackmate.execDetection()
	trackmate.computeSpotFeatures(true)
	trackmate.execSpotFiltering(true)

	return model
}	

/**
 * Displays TrackMate SpotCollection
 * 
 * @param spots
 *            TrackMate SpotCollection
 * 
 * @param img
 *            ImagePlus to display spots on 
 */
def displaySpots(spots, imp) {

	model = new Model();
	model.setSpots(spots, false);
	selectionModel = new SelectionModel(model)
	view =  new HyperStackDisplayer(model, selectionModel, imp)
	view.render()
	view.refresh()
}

/**
 * Displays two TrackMate SpotCollection over an image
 * 
 * @param spots1
 *            1st TrackMate SpotCollection
 * @param spots2
 *            2nd TrackMate SpotCollection
 * @param img
 *            ImagePlus to display spots on 
 */
def displaySpots(spots1, spots2, imp) {
	allSpots = new SpotCollection()
	for (spot in spots1.iterable(true)) {
		spot.putFeature("CHANNEL", 2);
		allSpots.add(spot, 0);
	}
	for (spot in spots2.iterable(true)) {
		spot.putFeature("CHANNEL", 1);
		allSpots.add(spot, 0);
	}
	modelAll = new Model();
	modelAll.setSpots(allSpots, false);
	selectionModel = new SelectionModel(modelAll)
	view =  new HyperStackDisplayer(modelAll, selectionModel, imp)
	spotColorGen = new SpotColorGenerator(modelAll) 
	spotColorGen.setFeature("CHANNEL")
	spotColorGen.setMinMax(0, 2)
	view.setDisplaySettings("SpotColoring", spotColorGen)
	view.render()
}


/**
 * Constructs graph linking all spots within specified distance
 * 
 * @param spots
 *            TrackMate SpotCollection
 * @param linkingDist
 *            Distance threshold for spot linking
 * @return graph
 *            HashMap, containing a list of spots linked to each spot in the collection
 */
def constructGraph(spots, linkingDist) {
	
	graph = new HashMap<Spot, ArrayList<Spot>>()
	
	// Loop through spots
	for (sp1 in spots.iterable(true)) {
		// Add spot to graph, along with empty ArrayList to contain linked spots
		graph.put(sp1, new ArrayList<Spot>())
		// Get spot position
		spotCoords1 = [sp1.getDoublePosition(0), sp1.getDoublePosition(1), sp1.getDoublePosition(2)] 
		// Loop through spots again
		for (sp2 in spots.iterable(true)) {
			// Get distance to sp1
			spotCoords2 = [sp2.getDoublePosition(0), sp2.getDoublePosition(1), sp2.getDoublePosition(2)] 
			dist = coordDistance(spotCoords1, spotCoords2)
			// Add to list if less than specified distance (and not 0)
			if (dist <= linkingDist && dist > 0) 
				graph.get(sp1).add(sp2)	
		}
	}
	return graph
}


def coordDistance(coords1, coords2) {
	sum = 0
	for (i = 0; i < coords1.size(); i++) 
		sum += (coords1[i] - coords2[i]) * (coords1[i] - coords2[i])
	return Math.sqrt(sum)
	
}


/**
 * Finds connected components in graph
 * 
 * @param graph
 *            HashMap, containing a list of spots linked to each spot in the collection
 * @return ccs
 *            HashMap containing pairs of IDs and Spots for each connected component of the graph
 */
def connectedComponents(graph) {
	count = 0
    marked = new HashMap<Spot, Boolean>()
	for(key in graph.keySet()) 
		marked.put(key, false)
	ccs = new HashMap<Integer, ArrayList<Spot>>()
    for (v in graph.keySet()) {
        if (!marked.get(v)) {
            ccs.put(count, new ArrayList<Spot>())
            dfs(graph, v)
            count++
        }
    }
    return ccs
}

def dfs(graph, v) {
    marked.put(v, true)
	ccs.get(count).add(v)
    for (w in graph.get(v)) {
        if (!marked.get(w)) 
            dfs(graph, w)    
    }
}

/**
 * Creates a user dialog asking user to specify a class
 * 
 * @param numClasses
 *            Number of allowed classes
 * @return val
 *            Class specifed by user
 */
def askForClass(numClasses) {
	new WaitForUserDialog("ok?").show()
	gd = new GenericDialog("Define class")
	gd.addNumericField("Define class, maxium: " + numClasses + ", 0: no class", 0, 0)
	gd.showDialog()
	while(true) {
		val = gd.getNextNumber()
		if (gd.invalidNumber() == false && val >= 0 && val <= numClasses) {
			break
		} else {
			println "please return a valid integer" 
			gd.showDialog()
		}
	}
	return val
}

/**
 * Creates a user dialog asking user for a Boolean value
 * 
 * @param message
 *            Display message
 * @return val
 *            Boolean value from user
 */
def askForBoolean(message) {
	gd = new GenericDialog("Quesion")
	gd.addCheckbox(message, true);
	gd.showDialog()
	val = gd.getCheckbox()
	return val
}


/**
 * Construct a ResultsTable containing spot statstics for all visible spots
 * 
 * @param model
 *            the Model containg the SpotCollection
 * 
 * @return the ResultsTable containg the Spot statistics
 * 
 */
def calcAllSpotFeatures(model) {
	
	// Create table
	spotTable = new ResultsTable()
	fm = model.getFeatureModel()
	
	// Get list of all spot features
	spotFeatures = model.getFeatureModel().getSpotFeatures()

	// Iterate though all visible spots
	spots = model.getSpots()
	for (spot in spots.iterable(true)) {

			// add spot name and ID to table
			spotTable.incrementCounter()
			spotTable.addLabel(spot.getName())
			spotTable.addValue( "ID", spot.ID())

			// add remaining spot features to table
			for (feature in spotFeatures ) {
				val = spot.getFeature(feature)
				spotTable.addValue(feature, val.doubleValue())
	
			}
		}
	return spotTable
}

