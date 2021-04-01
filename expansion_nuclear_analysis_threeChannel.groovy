import java.util.ArrayList
import java.awt.Rectangle
import java.util.Random

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



#@Integer(label="Satallite spot channel 1", value=1) satChannel1
#@Integer(label="Satallite spot channel 2", value=3) satChannel2
#@Integer(label="Central spot channel", value=2) centralChannel

#@File(label = "Input directory", style = "directory") srcFile
#@String(label = "File extension", value=".ims") ext

#@Double(label="Specify spot radius (microns) for sat 1", value=0.75) spotRadiusS1
#@Double(label="Specify spot radius (microns) for sat 2", value=0.75) spotRadiusS2
#@Double(label="Specify spot radius (microns) for central channel", value=0.75) spotRadiusC

#@Double(label="Specify quality threshold for spots for sat 1", value=6.) spotThreshQualS1
#@Double(label="Specify quality threshold for spots for sat 1", value=6.) spotThreshQualS2
#@Double(label="Specify quality threshold for spots for central channel", value=7.5) spotThreshQualC

#@Double(label="Specify ball radius C1 (pixels) for spots for sat 1", value=40) ballSizeS1
#@Double(label="Specify ball radius C2 (pixels) for spots for sat 1", value=40) ballSizeS2
#@Double(label="Specify ball radius C2 (pixels) for central channel", value=14) ballSizeC


#@Double(label="XY Voxel calibration (microns)", value=0.159) xyCal
#@Double(label="Z Voxel calibration (microns)", value=0.500) zCal




#@Double(label="Linking distance", value=0) linkingDist
#@Integer(label="Box crop width (microns)", value=10) cropWidth
#@Integer(label="Max display channel", value=2) maxDisplayChannel

#@Integer(label="Number of classes", value=5) numClasses


#@Boolean(label="Save cropped site images", value=false) saveCrops

IJ.run("Close All", "")
// output table to hold summary stats
output = new ResultsTable()

// create empty TrackMate Settings
settings = new Settings()
// use a Laplacian of Gaussian detector
settings.detectorFactory = new LogDetectorFactory()
// specify that inensity based features should be calculated for spots
settings.addSpotAnalyzerFactory(new SpotIntensityAnalyzerFactory<>())
// specify the quality threshold for spots
//settings.addSpotFilter(new FeatureFilter('QUALITY', spotThreshQual, true))

// set detector settings
settings.detectorSettings = settings.detectorFactory.getDefaultSettings()



println settings.detectorSettings

rand = new Random()

processFolder()
// show the output table
output.show("Results")
// save to csv file
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
	
	
	cali = imp.getCalibration()
	cali.pixelWidth = xyCal
	cali.pixelHeight = xyCal
	cali.pixelDepth = zCal
	cali.setUnit('micron')
	imp.setDisplayMode(IJ.COMPOSITE)

	impS1 = new Duplicator().run(imp, satChannel1, satChannel1, 1, imp.nSlices, 1, 1)
	IJ.run(impS1, "Subtract Background...", "rolling=" + ballSizeS1 + " stack")
	impS2 = new Duplicator().run(imp, satChannel2, satChannel2, 1, imp.nSlices, 1, 1)
	IJ.run(impS1, "Subtract Background...", "rolling=" + ballSizeS1 + " stack")
	impC = new Duplicator().run(imp, centralChannel, centralChannel, 1, imp.nSlices, 1, 1)
	IJ.run(impC, "Subtract Background...", "rolling=" + ballSizeC + " stack")


	// point TrackMate at image
	settings.setFrom(impS1)
	settings.detectorSettings['RADIUS'] = spotRadiusS1
	settings.detectorSettings['THRESHOLD'] = spotThreshQualS1
	settings.detectorSettings['TARGET_CHANNEL'] = 1
	
	modelS1 = detectSpots(settings)
	spotsS1 = modelS1.getSpots()

	settings.setFrom(impS2)
	settings.detectorSettings['RADIUS'] = spotRadiusS2
	settings.detectorSettings['THRESHOLD'] = spotThreshQualS2
	settings.detectorSettings['TARGET_CHANNEL'] = 1
	
	modelS2 = detectSpots(settings)
	spotsS2 = modelS2.getSpots()

	settings.setFrom(impC)
	settings.detectorSettings['RADIUS'] = spotRadiusC
	settings.detectorSettings['THRESHOLD'] = spotThreshQualC
	settings.detectorSettings['TARGET_CHANNEL'] = 1
	
	modelC = detectSpots(settings)
	spotsC = modelC.getSpots()


	spotTableS1 = calcAllSpotFeatures(modelS1)
	outSpotPathS1 =  file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf("."))  + "_spotStatsS1.csv"
	spotTableS1.save(outSpotPathS1)
	spotTableS2 = calcAllSpotFeatures(modelS2)
	outSpotPathS2 =  file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf("."))  + "_spotStatsS2.csv"
	spotTableS2.save(outSpotPathS2)
	spotTableC = calcAllSpotFeatures(modelC)
	outSpotPathC =  file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf("."))  + "_spotStatsC.csv"
	spotTableC.save(outSpotPathC)



	displaySpots(spotsS1, spotsS2, spotsC, imp)
	
	

	
	
	
	graph = constructGraph(spotsC, linkingDist)
	
	// find the connected components of the graph
	ccs = connectedComponents(graph)
	IJ.run("Brightness/Contrast...")

	// Create table
	siteTable = new ResultsTable()
	count = 0
	for (cc in ccs) {

		// find centre of mass for each connected component
		siteCoords = new double[3]
		ccNumSpots = cc.getValue().size()
		siteSpots = new SpotCollection()
		
		for (spot in cc.getValue()) {

			for (d = 0; d < 3; d++) 
				siteCoords[d] += spot.getDoublePosition(d)		
				
			
			
		}
		for (d = 0; d < 3; d++) 
			siteCoords[d] = siteCoords[d] / ccNumSpots


		displaySpot = new Spot(cropWidth / 2, cropWidth / 2, cropWidth / 2, 2, 0)
		siteSpots.add(displaySpot, 0);
		
		xMin = (int) Math.floor((siteCoords[0] - cropWidth / 2) / xyCal)
		yMin = (int) Math.floor((siteCoords[1] - cropWidth / 2) / xyCal)
		zMin = (int) Math.floor((siteCoords[2] - cropWidth / 2) / zCal) + 1
		zMax = (int) Math.floor((siteCoords[2] + cropWidth / 2) / zCal) + 1
		xyCropWidth = (int) Math.floor(cropWidth / xyCal)

	
		imp.setRoi(xMin, yMin, xyCropWidth, xyCropWidth)

		
		
		impCrop = new Duplicator().run(imp, 1, maxDisplayChannel, zMin, zMax, 1, 1)
		impCrop.setDisplayMode(IJ.COMPOSITE)
		
		displaySpots(siteSpots, impCrop)
		
		//impCrop.show()
		impCrop.setPosition(1, (int) Math.ceil(impCrop.nSlices / 2), 1)

		for (d = 1; d <= maxDisplayChannel; d++) {
			impCrop.setC(d)
			IJ.run(impCrop, "Enhance Contrast", "saturated=0.35")
		}

		wind = WindowManager.getCurrentWindow()
		wind.setLocationAndSize((int) wind.LEFT_ALIGNMENT, (int) wind.TOP_ALIGNMENT, wind.width * 5, wind.height * 5)
		
		cat = askForClass(numClasses)
		IJ.log("cat: " + cat)
		//cat = rand.nextInt(numClasses + 1)
		
		
		siteTable.incrementCounter()
		siteTable.addValue("ID", count)
		siteTable.addValue("POSITION_X", siteCoords[0])
		siteTable.addValue("POSITION_Y", siteCoords[1])
		siteTable.addValue("POSITION_Z", siteCoords[2])
		siteTable.addValue("Class", cat)

		if (saveCrops) {

			outCropPath = file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf("."))  + "_class" + cat + "_siteID" + count + ".tif"
			IJ.saveAs(impCrop, "Tiff", outCropPath);
		}


		
		impCrop.close()
		count++
		
	
	
	
		
		outSitePath =  file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf("."))  + "_siteStats.csv"
		siteTable.save(outSitePath)
	
	
		// add stats to output table
		output.incrementCounter()
		output.addValue("File", file.getAbsolutePath())
	
		output.addValue("Num spots S1", spotsS1.getNSpots(true))
		output.addValue("Num spots S2", spotsS2.getNSpots(true))
		output.addValue("Num spots C", spotsC.getNSpots(true))
		output.addValue("Num sites", count)
	}	
	
	imp.close()
	
	IJ.run("Close All", "")
}

def detectSpots(settings) {
	// create empty TrackMate Model
	model = new Model()
	// create TrackMate object to peform spot detection and tracking
    trackmate = new TrackMate(model, settings)
    
    ok = trackmate.checkInput()
    if (ok == false) 
    	println trackmate.getErrorMessage()
    // Find spots
    trackmate.execDetection()
	trackmate.computeSpotFeatures(true)
	trackmate.execSpotFiltering(true)

	return model
}	

def displaySpots(spots, imp) {

	model = new Model();
	model.setSpots(spots, false);
	
	selectionModel = new SelectionModel(model)
	view =  new HyperStackDisplayer(model, selectionModel, imp)
	view.render()
	view.refresh()
}

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


def displaySpots(spots1, spots2, spots3, imp) {
	allSpots = new SpotCollection()
	for (spot in spots1.iterable(true)) {
		spot.putFeature("CHANNEL", 2);
		allSpots.add(spot, 0);
	}
		for (spot in spots2.iterable(true)) {
		spot.putFeature("CHANNEL", 3);
		allSpots.add(spot, 0);
	}
	for (spot in spots3.iterable(true)) {
		spot.putFeature("CHANNEL", 1);
		allSpots.add(spot, 0);
	}
	modelAll = new Model();
	modelAll.setSpots(allSpots, false);
	
	selectionModel = new SelectionModel(modelAll)
	view =  new HyperStackDisplayer(modelAll, selectionModel, imp)
	spotColorGen = new SpotColorGenerator(modelAll) 
	spotColorGen.setFeature("CHANNEL")
	spotColorGen.setMinMax(0, 3)
	view.setDisplaySettings("SpotColoring", spotColorGen)
	view.render()
}



def coordDistance(coords1, coords2) {

	sum = 0
	for (i = 0; i < coords1.size(); i++) 
		sum += (coords1[i] - coords2[i]) * (coords1[i] - coords2[i])
	return Math.sqrt(sum)
	
}

def constructGraph(spots, linkingDist) {
	
	graph = new HashMap<Spot, ArrayList<Spot>>()
	for (sp1 in spots.iterable(true)) {
		graph.put(sp1, new ArrayList<Spot>())
		spotCoords1 = [sp1.getDoublePosition(0), sp1.getDoublePosition(1), sp1.getDoublePosition(2)] 
		for (sp2 in spots.iterable(true)) {
			spotCoords2 = [sp2.getDoublePosition(0), sp2.getDoublePosition(1), sp2.getDoublePosition(2)] 
			dist = coordDistance(spotCoords1, spotCoords2)
			if (dist <= linkingDist && dist > 0) 
				graph.get(sp1).add(sp2)	
		}
	}
	return graph
}

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

