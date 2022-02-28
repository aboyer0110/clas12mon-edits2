package org.clas.detectors;

import org.clas.viewer.DetectorMonitor;
import org.jlab.geom.prim.Point3D;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.Axis;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.List;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.graphics.IDataSetPlotter;



public class RICHmonitor  extends DetectorMonitor {

    private static final int NPMT        = 391; //number of PMTs in the RICH
    private static final int NANODE      = 64; //number of channels in PMT
    private static final int NTILE       = 138; //number of front end tiles
    private static final double MAXMAP   = 5; //maximum something?
    private static final double MAXPMT   = 64*4;//maximum something? Related to channels
    private static final double MINPMT   = MAXPMT/1000;//Minimum something? Again related to channels
    private static final double MAXPIXEL = 20;//Maximum pixel? Not sure the meaning
    private static final double MINPIXEL = MAXPIXEL/1000;//Minimum pixel? Not sure the meaning
    private static final double MAXTIME  = 6;//Max time interval between counts maybe?
    private static final int NPMTROWS        = 23;//number of tile rows
    private static final int NPIXELROWS      = 8;//number of pixel rows in PMT
    private static final int NPIXELCOLUMNS   = 8;//number of pixel columns
    private static final double PIXELSIZE    = 1;//size of a channel in the GUI (1 pixel?)
    private static final double PMTCLEARANCE = 1;//Not sure what this is
  
    private static final double LE  = 100;//Data line for TDC Leading edge
    private static final double TOT = 60;//Upper dataline threshold for TDC hit time
    private static final double XT  = 25;//Lower dataline Threshold for TDC hit time
    private final int[] CHAN2PIX = {60, 58, 59, 57, 52, 50, 51, 49, 44, 42, 43, 41, 36, 34, 35, 33, 28, 26, 27, 25, 20, 18, 19, 17, 12, 10, 
                                    11, 9, 4, 2, 3, 1, 5, 7, 6, 8, 13, 15, 14, 16, 21, 23, 22, 24, 29, 31, 30, 32, 37, 39, 38, 40, 45, 47, 
                                    46, 48, 53, 55, 54, 56, 61, 63, 62, 64};//array converts MAROC ID to pixel number
    private final Integer[] TWOTILERS = {3, 5, 7, 12, 15, 19, 24, 28, 33, 39, 44, 50, 57, 63, 70, 78, 85, 93, 102, 110, 119, 129, 138};//Tiles which contain Two PMTs
    private final int[] FIRSTPMTS = {1, 7, 14, 22, 31, 41, 52, 64, 77, 91, 106, 122, 139, 157, 176, 196, 217, 239, 262, 286, 311, 337, 364};//First PMTs in row   
            
    private IndexedList<Integer> tileToPMT = null;//empty tileToPMT list before beginning
    
    public RICHmonitor(String name) {
        super(name);

        this.setDetectorTabNames("Occupancy and time","Occupancy Map","TDC"); //Initialize sub-tabs in the RICH tab
        this.tileToPMT = this.setTiletoPMTMap();//initialize a tileToPMT list?
        this.init(false);
    }
    
    private IndexedList setTiletoPMTMap() {
        IndexedList<Integer> list = new IndexedList<>(2);
        List<Integer> twoTilers = Arrays.asList(TWOTILERS);
        int pmt=0;
        for(int i=0; i<NTILE; i++) {//Start identifier at zero. If it is less than 138, increment
            int tile = i+1;//tile # = 1 + loop index
            for(int j=0; j<3; j++) {//Start identifier at zero. If it is less than 3, increment. This tracks the # of PMT in the tile
                if(j==1 && twoTilers.contains(tile)) continue; //if j==1 and the tile matches one of the twoTilers, skip.
                pmt++;
                list.add(pmt, tile, j);//build the TiletoPMTMap
            }
        }
        return list;
    }


    @Override
    public void createHistos() {
        H2F hi_pmt_leading_edge = new H2F("hi_pmt_leading_edge", "TDC Hit Leading Edge Time", NPMT, +0.5, NPMT+0.5, 100, 0, 300);//TDC Leading edge plot
        hi_pmt_leading_edge.setTitleX("PMT");//X-axis label for TDC leading edge
        hi_pmt_leading_edge.setTitleY("Time (ns)");//Y-axis label for TDC leading edge
        H2F hi_pmt_duration     = new H2F("hi_pmt_duration", "TDC Hit Time Over Threshold",   NPMT, +0.5, NPMT+0.5, 100, 0,  100);//TDC Hit time over threshold
        hi_pmt_duration.setTitleX("PMT");//X-axis label for Hit time plot
        hi_pmt_duration.setTitleY("Time (ns)");//Y-axis label for Hit time plot
        
        H1F delta_tdc = new H1F("delta_tdc", "Time (ns)", "Counts", 100, 0, 100);
        delta_tdc.setTitle("RICH Delta TDC");
        //delta_tdc.getXaxis().;
        delta_tdc.setLineWidth(2);
        delta_tdc.setOptStat(1111);
        
        H1F hi_pmt_occupancy    = new H1F("hi_pmt_occupancy", "PMT",   "Counts", NPMT, +0.5, NPMT+0.5);//PMT occupancy plot
        hi_pmt_occupancy.setTitle("PMT Hit Occupancy");//Title PMT occupancy plot
        hi_pmt_occupancy.setFillColor(25);//define a color to use
        hi_pmt_occupancy.setOptStat("10");//displays the # of entries in top right
        H1F hi_pmt_max          = new H1F("hi_pmt_max", " ", 1, 0.5, NPMT+0.5);//Sets a maximum occupancy?
        hi_pmt_max.setLineWidth(2);//Draws a max line on plot
        hi_pmt_max.setLineColor(2);//colors the max line
        H1F hi_pix_occupancy    = new H1F("hi_pix_occupancy", "Pixel", "Counts", NPMT*NANODE, -0.5, NPMT*NANODE-0.5);//pixel occupancy
        hi_pix_occupancy.setTitle("Pixel Hit Occupancy");//title pixel occupancy plot
        hi_pix_occupancy.setFillColor(25);//define a color to use. same of PMT occupancy
        hi_pix_occupancy.setOptStat("10");//displays the # of entries in the top right
        H1F hi_pix_max          = new H1F("hi_pix_max", " ", 1, 0.5, NPMT*NANODE+0.5);//define a max pixel occupancy        
        hi_pix_max.setLineWidth(2);//Draws a max line on plot
        hi_pix_max.setLineColor(2);//colors the max line
        H2F hi_scaler           = new H2F("hi_scaler", "TDC Hit Map",260, -130, 130, 207, 0, 207); //plots the TDC time over threshold (length of signal)
        H1F hi_summary          = new H1F("summary", "PMT",   "Counts", NPMT, +0.5, NPMT+0.5);//Counts plot on summary tab
        hi_summary.setTitle("RICH");//title summary plot
        hi_summary.setFillColor(25);//fill color
        DataGroup dg = new DataGroup(1,5); //creating a data group "dg", 1x5
        dg.addDataSet(hi_pmt_leading_edge, 0);//(ID of data set, order)
        dg.addDataSet(hi_pmt_duration,     1);
        
        dg.addDataSet(delta_tdc, 1);
        
        dg.addDataSet(hi_pmt_occupancy,    2);
        dg.addDataSet(hi_pmt_max,          2);
        dg.addDataSet(hi_pix_occupancy,    3);
        dg.addDataSet(hi_pix_max,          3);
        dg.addDataSet(hi_scaler,           4);
        this.getDataGroup().add(dg,0,0,0);//adds something to dg
        DataGroup sum = new DataGroup(1,1);//creating new data group "sum", 1x1
        sum.addDataSet(hi_summary, 0);//(ID of data set, order)
        this.setDetectorSummary(sum);//Sets the data group for the summary tab?
    }


    @Override
    public void plotHistos() {
        this.getDetectorCanvas().getCanvas("Occupancy and time").divide(2, 2);//divide canvas into 4 plots
        for(String tab: this.getDetectorTabNames()) {//pulls detector tab names
            this.getDetectorCanvas().getCanvas(tab).setGridX(false);//no X grid
            this.getDetectorCanvas().getCanvas(tab).setGridY(false);//no Y grid
        }
        
        DataLine YlineTOT = new DataLine(TOT,0,TOT,30e5);
        DataLine YlineXT = new DataLine(XT,0,XT,30e5);
        this.getDetectorCanvas().getCanvas("TDC").draw(this.getDataGroup().getItem(0,0,0).getH1F("delta_tdc"));
        this.getDetectorCanvas().getCanvas("TDC").draw(YlineXT);
        this.getDetectorCanvas().getCanvas("TDC").draw(YlineTOT);
        this.getDetectorCanvas().getCanvas("TDC").update();
        
        Axis pmtAxis = this.getDataGroup().getItem(0,0,0).getH2F("hi_pmt_duration").getXAxis();//Determines X-axis of the TDC hit duration plot
        DataLine lineLE  = new DataLine(pmtAxis.min(),LE,  pmtAxis.max(), LE);//Creates a data line for I think the leading edge plot
        DataLine lineTOT = new DataLine(pmtAxis.min(),TOT, pmtAxis.max(), TOT);//Creates a data line (top) for the TDC hit duration
        DataLine lineXT  = new DataLine(pmtAxis.min(),XT,  pmtAxis.max(), XT);//Creates a data line (bottom) for the TDC hit duration
        this.getDetectorCanvas().getCanvas("Occupancy and time").cd(0);//Get the top left canvas
        this.getDetectorCanvas().getCanvas("Occupancy and time").getPad(0).getAxisY().setLog(true);//set log scale
        this.getDetectorCanvas().getCanvas("Occupancy and time").draw(this.getDataGroup().getItem(0,0,0).getH1F("hi_pmt_occupancy"));//Draw the PMT occupancy plot
        this.getDetectorCanvas().getCanvas("Occupancy and time").draw(this.getDataGroup().getItem(0,0,0).getH1F("hi_pmt_max"),"same");//Draws the max line?
        this.getDetectorCanvas().getCanvas("Occupancy and time").cd(1);//Gets the top right canvas
        this.getDetectorCanvas().getCanvas("Occupancy and time").getPad(1).getAxisY().setLog(true);//set log scale
        this.getDetectorCanvas().getCanvas("Occupancy and time").draw(this.getDataGroup().getItem(0,0,0).getH1F("hi_pix_occupancy"));//Draws the Pixel occupancy
        this.getDetectorCanvas().getCanvas("Occupancy and time").draw(this.getDataGroup().getItem(0,0,0).getH1F("hi_pix_max"),"same");//Draws the max line?
        this.getDetectorCanvas().getCanvas("Occupancy and time").cd(2);//Get the bottom left canvas
        this.getDetectorCanvas().getCanvas("Occupancy and time").getPad(2).getAxisZ().setLog(getLogZ());//set log scale
        this.getDetectorCanvas().getCanvas("Occupancy and time").draw(this.getDataGroup().getItem(0,0,0).getH2F("hi_pmt_leading_edge"));//Draws the leading edge plot
        this.getDetectorCanvas().getCanvas("Occupancy and time").draw(lineLE);//Draws the data line on the leading edge plot
        this.getDetectorCanvas().getCanvas("Occupancy and time").cd(3);//grabs the bottom right canvas
        this.getDetectorCanvas().getCanvas("Occupancy and time").getPad(3).getAxisZ().setLog(getLogZ());//set log scale
        this.getDetectorCanvas().getCanvas("Occupancy and time").draw(this.getDataGroup().getItem(0,0,0).getH2F("hi_pmt_duration"));//Draws the TDC hit duration
        this.getDetectorCanvas().getCanvas("Occupancy and time").draw(lineXT);//draws the lower data line
        this.getDetectorCanvas().getCanvas("Occupancy and time").draw(lineTOT);//draws the upper data line
        this.getDetectorCanvas().getCanvas("Occupancy and time").update();//update the GUI canvas
        this.getDetectorCanvas().getCanvas("Occupancy Map").getPad(0).setPalette("kRainBow");//set the color map for the PMT window
        this.getDetectorCanvas().getCanvas("Occupancy Map").draw(this.getDataGroup().getItem(0,0,0).getH2F("hi_scaler"));//draw the hi_scaler plot on pad(0)
        this.getDetectorCanvas().getCanvas("Occupancy Map").getPad(0).getAxisZ().setLog(true);//set log scale
        this.getDetectorCanvas().getCanvas("Occupancy Map").update();//update
        this.DrawTiles();//draw tile outlines
    }


    @Override
    public void processEvent(DataEvent event) {
        
        // process event info and save into data group
        if(event.hasBank("RICH::tdc")) {//if the event has RICH info then...
            Map<Integer, ArrayList<TDCHit>> tdcMap = new HashMap<>();//create a HashMap that maps data to integers
 
            DataBank  bank = event.getBank("RICH::tdc");//Store the data bank from the event
            int rows = bank.rows();//initialize a # of rows of data based on the bank
            
            for(int i = 0; i < rows; i++) {//For all rows of data...
                int sector = bank.getByte("sector",i);     //4 by default (only 1 RICH at the time of writing)
                int  layer = bank.getByte("layer",i);      //byte variable, ranges from -127 to 127
                int   tile = layer & 0xFF;                 //conversion of byte to int variable, ranges from 1 to 138 (Tile Number)
                int   comp = bank.getShort("component",i); //short variable, comp is the MAROC ID shifted by 1 (ranges 1-192)
                int    tdc = bank.getInt("TDC",i);         //TDC value
                int  order = bank.getByte("order",i);      //order specifies leading or trailing edge. 1=leading, 0=trailing

                int anode = CHAN2PIX[(comp-1) % NANODE];//from 1 to 64
                int asic  = (comp-1) / NANODE;//ASIC value assigns which PMT in a tile it is, ranges 0, 1, or 2
                int pmt   = this.tileToPMT.getItem(tile, asic);//Set the PMT based off the ASIC info determined from the data bank
                int pixel = (pmt-1)*NANODE + anode-1;//from 0 to 25023

                if(tdc>0) {
                    
                    if(!tdcMap.containsKey(pixel) && order==1) {//if the tdc map contains a pixel with order 1
                        tdcMap.put(pixel, new ArrayList<>());//adds mapping to the tdcmap
                        tdcMap.get(pixel).add(new TDCHit(tdc));//retrive pixel and add tdc info to it
                    }
                    else {//if tdc<=0
                        if(order==1) {//if its leading
                            tdcMap.get(pixel).add(new TDCHit(tdc));//retrieve pixel and add the tdc info to it
                        }
                        else {//if order==0 (trailing)
                            TDCHit last = tdcMap.get(pixel).get(tdcMap.get(pixel).size()-1);//Not sure what this means, it stores the trailing info somewhere
                            if(last.getDuration()==0) last.setTrailingEdge(tdc);//again, not sure?
                        }
                    }
                }

            }
            for(int pixel : tdcMap.keySet()) {//keySet maps the elements of a hashmap to a set
                for(TDCHit hit : tdcMap.get(pixel)) {//retrieve the pixel
                    if(hit.getDuration()>0) {//if the hit duration is not zero, aka, if there's a detection
                        int pmt   = pixel/NANODE + 1;//determine the PMT in question
                        int anode = pixel%NANODE + 1;//determine the anode in question
                        this.getDataGroup().getItem(0,0,0).getH2F("hi_pmt_leading_edge").fill(pmt,hit.getLedingEdge());//fills the leading edge plot
                        this.getDataGroup().getItem(0,0,0).getH2F("hi_pmt_duration").fill(pmt,hit.getDuration());//fills the TDC hit duration plot
                    
                        this.getDataGroup().getItem(0,0,0).getH1F("delta_tdc").fill(hit.getDuration());
                        
                        Point3D pxy = this.getCoordinates(pmt, anode);//retrieve the X,Y,Z coords for the given PMT and Pixel
                        this.getDataGroup().getItem(0,0,0).getH2F("hi_scaler").fill(pxy.x(), pxy.y());//fill the (X,Y) location of the PMT and Pixel with the occupancy data
                    
                        this.getDetectorSummary().getH1F("summary").fill(pmt);//fill the RICH summary plot
                        this.getDataGroup().getItem(0,0,0).getH1F("hi_pmt_occupancy").fill(pmt);//fill the PMT occupancy plot
                        this.getDataGroup().getItem(0,0,0).getH1F("hi_pix_occupancy").fill(pixel);//fill the pixel occupancy plot
                    }
                }
            }
        }
    }





    @Override
    public void analysisUpdate() {
        double nentries = this.getDataGroup().getItem(0,0,0).getH2F("hi_scaler").getEntries();//Define the number of entries as a double variable by retrieving them from the occupancy plot
        double average  = nentries/NPMT/NANODE;//calculate the average of entries
        
        this.getDataGroup().getItem(0,0,0).getH1F("delta_tdc").setBinContent(0, average*MAXTIME);
        
        this.getDataGroup().getItem(0,0,0).getH1F("hi_pmt_max").setBinContent(0, average*MAXPMT);//replace the existing PMT occupancy histogram content, ie update it
        this.getDataGroup().getItem(0,0,0).getH1F("hi_pix_max").setBinContent(0, average*MAXPIXEL);//replace the existing Pixel occupancy histogram content
        this.setYAxisMin(this.getDetectorCanvas().getCanvas("Occupancy and time").getPad(0),average*MINPMT);////automatically set the axis bound on the top left?
        this.setYAxisMin(this.getDetectorCanvas().getCanvas("Occupancy and time").getPad(1),average*MINPIXEL);//automatically set the axis bound on the top right?
        this.getDetectorCanvas().getCanvas("Occupancy and time").getPad(2).getAxisZ().setRange(0, average*MAXTIME);//set the Z-axis range on the bottom left. What is the Z-axis? The scale on the right?
        this.getDetectorCanvas().getCanvas("Occupancy and time").getPad(3).getAxisZ().setRange(0, average*MAXTIME);//set the Z-axis range on the bottom right.
        this.getDetectorCanvas().getCanvas("Occupancy Map").getPad().getAxisZ().setRange(0, average*MAXMAP);//Set the Z-axis range. Is the Z-axis the scale?
    }
    
    private void setYAxisMin(EmbeddedPad pad, double min) {//Function to determine the Y-axis min (see above)
        double max = 0;//clear max at the beginning
        for(IDataSetPlotter hh : pad.getDatasetPlotters()) {//DatasetPlotters sends info to the client for processing
            max = Math.max(max, 1.2*hh.getDataSet().getMax());//Calculate the maximum
        }
        min = Math.max(min,0.1);//Calculate the minimum
        if(max<=min) max = min*10;//if the max is less than min, then change the max
        pad.getAxisY().setRange(min, max);//Use these calculations to set the range
    }

    private class TDCHit {//Class involving TDC assignments
        private int leadingEdge;//initialize leading edge
        private int trailingEdge;//initialize trailing edge
        
        TDCHit(int tdc) {//Function called in processEvent
            this.leadingEdge = tdc;//set leading edge equal to tdc
        }
        
        public void setTrailingEdge(int tdc) {//function called in processEvent
            this.trailingEdge = tdc;//set trailing edge equal to tdc
        }
        
        public int getLedingEdge() {//return the leading edge
            return this.leadingEdge;
        }
        
        public int getDuration() {//return the duration of the PMT signal
            if(this.trailingEdge>this.leadingEdge)
                return this.trailingEdge-this.leadingEdge;//duration = difference b/w trailing and leading
            else//if the leading edge is larger, return nothing
                return 0;
        }
    }
    
    private void DrawTiles() {//function that draws the Tile GUI on the Occupancy window
        for(int i=0; i<NTILE; i++) {//Initialize loop index; provided its less than the total tiles, increment and then...
            int tile = i+1;//assign tile number
            int pmt0 = this.tileToPMT.getItem(tile,2);//right-most PMT I think, not reflective of the mapping
            int pmt2 = this.tileToPMT.getItem(tile,0);//left-most PMT I think, not reflective of the mapping
            Point3D p1 = this.getCoordinates(pmt0, 1);//define a 3D point. pmt0 is the index, 1 is the coordinate for the object. (bottom left of 1st tile)
            Point3D p2 = this.getCoordinates(pmt0, 57);//define a 3D point. pmt0 is the index, 57 is the coordinate for the object. (top left of 1st tile)
            Point3D p3 = this.getCoordinates(pmt2, 64);//define a 3D point. pmt2 is the index, 64 is the coordinate for the object. (top right of 1st tile)
            Point3D p4 = this.getCoordinates(pmt2, 8);//define a 3D point. pmt2 is the index, 8 is the coordinate for the object. (bottom right of 1st tile)
            DataLine line1 = new DataLine(p1.x()-0.5, p1.y()+1.5, p2.x()-0.5, p2.y()-0.5);
            DataLine line2 = new DataLine(p2.x()-0.5, p2.y()-0.5, p3.x()+1.5, p3.y()-0.5);
            DataLine line3 = new DataLine(p3.x()+1.5, p3.y()-0.5, p4.x()+1.5, p4.y()+1.5);
            DataLine line4 = new DataLine(p4.x()+1.5, p4.y()+1.5, p1.x()-0.5, p1.y()+1.5);//create data lines based on these points
            this.getDetectorCanvas().getCanvas("Occupancy Map").draw(line1);
            this.getDetectorCanvas().getCanvas("Occupancy Map").draw(line2);
            this.getDetectorCanvas().getCanvas("Occupancy Map").draw(line3);
            this.getDetectorCanvas().getCanvas("Occupancy Map").draw(line4);//draw the lines. These lines are the outlines of the tiles
        }
    }
        
    private Point3D getCoordinates(int pmt, int anode) {//function to retrieve coordinates

        /* Finding the row and column of the pmt */
        int row = this.getPMTRow(pmt);//calls function to determine PMT row
        int col = this.getPMTColumn(pmt, row);//calls function to determine PMT column

        /* Finding the position of the first anode of the row */
        Point3D anode1 = getRowAnode1(row);//see function below
        /* finding the local coordinates of the anode in the pmt */
        Point3D local = getLocalCoordinates(anode);

        double x = local.x() + anode1.x() + (col - 1) * (NPIXELCOLUMNS * PIXELSIZE + PMTCLEARANCE);
        double y = local.y() + anode1.y();
        //if (anode == 1) cout << "  localx=" << localx << "   localy=" << localy << "   x=" << *x << "  y=" << *y << endl;

        return new Point3D(x,y,0);
    }

    private Point3D getRowAnode1(int row) {
        /* Return the position of anode 1 of the leftmost pmt of the row 
        Anode 1 is top left pixel of the PMT 
        */
        double anode1x = - ((row-1)*0.5+3) * (NPIXELCOLUMNS * PIXELSIZE + PMTCLEARANCE);
        anode1x = Math.floor(anode1x);
        double anode1y =    row * (NPIXELCOLUMNS * PIXELSIZE + PMTCLEARANCE);

        return new Point3D(anode1x, anode1y, 0);
    }

    private Point3D getLocalCoordinates(int anode){
        /* anode start from 1 
        column is from 1 to 8
        row is from -1 to -8
        anode 1 is in column 1 and row -1
        anode 64 is in column 8 and row -8
        */

        /* row and column of the anode */
        int col =   1 + (anode - 1) % NPIXELROWS;
        int row = -(1 + (anode - 1) / NPIXELROWS);

        //cout << "anode=" << anode << "  row=" << row << "  col=" << column << endl;
        double localx = (col - 0) * PIXELSIZE;
        double localy = (row + 0) * PIXELSIZE;

        return new Point3D(localx, localy, 0);
    }

    private int getPMTRow(int pmt) {
        int row = NPMTROWS;
        for (int r = 0; r < NPMTROWS; r++) {
            if (pmt < this.FIRSTPMTS[r]) {
                row = r;
                break;
            }
        }
        return row;
    }

    private int getPMTColumn(int pmt, int row) {

        int col = 1 + pmt - this.FIRSTPMTS[row - 1];
        int nCols = getNColumns(row);
        col = 1 + nCols - col;

        return col;
    }

    private int getNColumns(int row) {
        if (row == NPMTROWS) {
            return NPMT - FIRSTPMTS[row - 1] + 1;
        }
        else {
            return FIRSTPMTS[row] - FIRSTPMTS[row - 1];
        }
    }

    
    
}
