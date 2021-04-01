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



#@Integer(label="Satallite spot channel", value=1) satChannel
#@Integer(label="Central spot channel", value=2) centralChannel

#@File(label = "Input directory", style = "directory") srcFile
#@String(label = "File extension", value=".ims") ext
#@Double(label="Specify spot radius (microns)", value=0.75) spotRadiusC1
#@Double(label="Specify spot radius (microns)", value=0.75) spotRadiusC2
#@Double(label="Specify quality threshold for spots", value=6.) spotThreshQualC1
#@Double(label="Specify quality threshold for spots", value=7.5) spotThreshQualC2
#@Double(label="Specify ball radius C1 (pixels)", value=40) ballSize1
#@Double(label="Specify ball radius C2 (pixels)", value=14) ballSize2

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

	impC1 = new Duplicator().run(imp, satChannel, satChannel, 1, imp.nSlices, 1, 1)
	IJ.run(impC1, "Subtract Background...", "rolling=" + ballSize1 + " stack")
	impC2 = new Duplicator().run(imp, centralChannel, centralChannel, 1, imp.nSlices, 1, 1)
	IJ.run(impC2, "Subtract Background...", "rolling=" + ballSize2 + " stack")

	///////////////////////////////////////
	//impC1.show()
	//impC2.show()
	//new WaitForUserDialog("ok?").show()
	/////////////////////////////////////

	// point TrackMate at image
	settings.setFrom(impC1)
	settings.detectorSettings['RADIUS'] = spotRadiusC1
	settings.detectorSettings['THRESHOLD'] = spotThreshQualC1
	settings.detectorSettings['TARGET_CHANNEL'] = 1
	
	modelC1 = detectSpots(settings)
	spotsC1 = modelC1.getSpots()
	
	settings.setFrom(impC2)
	settings.detectorSettings['RADIUS'] = spotRadiusC2
	settings.detectorSettings['THRESHOLD'] = spotThreshQualC2
	settings.detectorSettings['TARGET_CHANNEL'] = 1
	
	modelC2 = detectSpots(settings)
	spotsC2 = modelC2.getSpots()


	spotTableC1 = calcAllSpotFeatures(modelC1)
	outSpotPathC1 =  file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf("."))  + "_spotStatsC1.csv"
	spotTableC1.save(outSpotPathC1)
	spotTableC2 = calcAllSpotFeatures(modelC2)
	outSpotPathC2 =  file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf("."))  + "_spotStatsC2.csv"
	spotTableC2.save(outSpotPathC2)

	
	

	

	displaySpots(spotsC1, spotsC2, imp)
	
	
	ok = askForClass(2)

	if(ok == 1) {
	
	
		graph = constructGraph(spotsC2, linkingDist)
		
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
		}
	
	
	
		
		outSitePath =  file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf("."))  + "_siteStats.csv"
		siteTable.save(outSitePath)
	
	
		// add stats to output table
		output.incrementCounter()
		output.addValue("File", file.getAbsolutePath())
	
		output.addValue("Num spots C1", spotsC1.getNSpots(true))
		output.addValue("Num spots C2", spotsC2.getNSpots(true))
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

