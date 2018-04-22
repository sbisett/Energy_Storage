/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pkg2salt;

import javafx.application.Application;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.Scene;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Button;
import javafx.scene.layout.StackPane;
import javafx.stage.Stage;

/**
 *
 * @author erich
 */
public abstract class XYPlot extends Application {
	
    public static void MakePlot(Stage plot_stage, double[] xrange, double[] yrange, 
    		double[] xrange2, double[] yrange2, 
    		double[] xrange3, double[] yrange3, 
    		double[] xrange4, double[] yrange4, 
    		double[] xrange5, double[] yrange5, 
    		double[] xrange6, double[] yrange6,
    		double[] xrange7, double[] yrange7, 
    		double[] xrange8, double[] yrange8,
    		double[] xrange9, double[] yrange9,
    		String xlabel, String ylabel, String title, double[] yaxis) 
    {
	    //this loop determines a y max range
    		int l,m,n;
	    double max_x, max_y=0.;
	    n=yrange.length;
	     for (l=0; (l<(n-1)); l++) 
	     {
	        if(yrange[l]>max_y) 
	        {
	            max_y = yrange[l] + 100; //change if not in the range when plotted
	        }//if 
	     }// for l
	     
	     //creates x and y axis labels
	     //y axis is based on the passed values from Salt_Tanks_Main.java
	     max_x=xrange[xrange.length-1];
	     final NumberAxis x_axis = new NumberAxis(0,max_x, max_x/10.);
	     final NumberAxis y_axis = new NumberAxis (yaxis[0], yaxis[1], yaxis[1]/10.); //change if not in the range when plotted
	     final ScatterChart<Number,Number> sc = new ScatterChart<Number,Number>(x_axis,y_axis);
	     x_axis.setLabel(xlabel);
	     y_axis.setLabel(ylabel);
	     
	     //creates a data series for each sensitivity combination
	     XYChart.Series s1 = new XYChart.Series();
	     XYChart.Series s2 = new XYChart.Series();
	     XYChart.Series s3 = new XYChart.Series();
	     XYChart.Series s4 = new XYChart.Series();
	     XYChart.Series s5 = new XYChart.Series();
	     XYChart.Series s6 = new XYChart.Series();
	     XYChart.Series s7 = new XYChart.Series();
	     XYChart.Series s8 = new XYChart.Series();
	     XYChart.Series s9 = new XYChart.Series();
	     
	     //creates a key (names all the data series)
	     sc.setTitle(title);
	     s1.setName("550K, 350K, 20W/m2/K");
	     s2.setName("525K, 375K, 20W/m2/K");
	     s3.setName("500K, 400K, 20W/m2/K");
	     s4.setName("550K, 350K, 35W/m2/K");
	     s5.setName("525K, 375K, 35W/m2/K");
	     s6.setName("500K, 400K, 35W/m2/K");
	     s7.setName("550K, 350K, 50W/m2/K");
	     s8.setName("525K, 375K, 50W/m2/K");
	     s9.setName("500K, 400K, 50W/m2/K");
	     
	     //assigns arrays passed from Salt_Tanks_Main.java to the created data series
	     for (m=0; m<xrange.length-1; m++) 
	     {
	       s1.getData().add(new XYChart.Data(xrange[m],yrange[m]));
	       s2.getData().add(new XYChart.Data(xrange2[m],yrange2[m]));
	       s3.getData().add(new XYChart.Data(xrange3[m],yrange3[m]));
	       s4.getData().add(new XYChart.Data(xrange4[m],yrange4[m]));
	       s5.getData().add(new XYChart.Data(xrange5[m],yrange5[m]));
	       s6.getData().add(new XYChart.Data(xrange6[m],yrange6[m]));
	       s7.getData().add(new XYChart.Data(xrange7[m],yrange7[m]));
	       s8.getData().add(new XYChart.Data(xrange8[m],yrange8[m]));
	       s9.getData().add(new XYChart.Data(xrange9[m],yrange9[m]));  
	     }// for m
	     
	     //plots the data
	     sc.getData().addAll(s1, s2, s3,s4,s5,s6,s7,s8,s9);
	     Scene scene  = new Scene(sc, 500, 400);
	     plot_stage.setScene(scene);
	     
    }//public
    
    public static void main(String[] args) 
    {
        launch(args);
    }
}