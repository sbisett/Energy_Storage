/* This calculates the temperature/volume of the hot and cold tanks for 3 temperature ranges and 3 heat transfer 
 * coefficients. 
 */

package pkg2salt;

import javafx.application.Application;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.layout.StackPane;
import javafx.stage.Stage;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.io.FileNotFoundException;


/**
 *
 * @author erich
 */
public class Salt_Tanks_Main extends Application 
{
    
    @Override
    public void start(Stage primaryStage) throws IOException 
    {

            Stage fig2=new Stage();
            Stage fig3=new Stage();
            Stage fig4=new Stage();
            
            // These parameters are the inputs
            
            //may need to be changed according to the tank size
            //all sensitivity graphs in my box files use a tank height of 6
            final double TANK_RADIUS=6.;  // m - both tanks have same dimensions
            final double TANK_HEIGHT=6;  // m  // in the results on the Milestone 3 report a tank height of 12 was used
            
            final double STARTING_HOT_TANK_VOL_FRACT=0.75;  // total fluid volume equals the volume of one tank.  This is the fraction of the total fluid that's initially in the hot tank

            final double P_IN = 50.E6;  // J/s - rate at which THERMAL power is transferred in when charging
            final double P_OUT = 30.E6; // J/s - rate at which THERMAL power is transferred out when discharging
            
            //makes hot and cold tank target temperature arrays 
            double T_HL[] = new double[3];
            double T_CL[] = new double[3];
            
            //*** 9 combinations of variables are used for the sensitivity analysis (based on the temperature ranges, hx_cf, and type of salt)
            
            //ranges for hot and cold tank target temperatures are set (K)
            //change for sensitivity purposes
            T_HL[0] = 550;
            T_HL[1] = 525;
            T_HL[2] = 500;
            T_CL[0] = 350;
            T_CL[1] = 375;
            T_CL[2] = 400;
            
            double T_inf = 300.;  // K - the ambient temperature outside the tank
            double hx_cf[] = new double[] {20, 35, 50}; // heat transfer coeff for natural circulation convective heat xfer to environment - W/m2/K
            		//these are supposed to be a range of values to test out for sensitivity 

            //input which salt properties to use -- 3 options for sensitivity analysis
            //based on the salt entered, the c_p and rho will be assigned and used throughout the rest of the program 
            String salt = "Solar Salt"; //Solar Salt, Hitec XL or Glauber
            double c_p = 0;  // J/kg-K  value initiated with 0 
            double rho = 0; // kg/m3 value initiated with 0
            if (salt.equals("Solar Salt"))
            {
            		c_p = 1496;
            		rho = 1899;
            }
            else if(salt.equals("Hitec XL"))
            {
	            	c_p = 1447;
	        		rho = 1992;
            }
            else if(salt.equals("Glauber"))
            {
	            	c_p = 3349;
	        		rho = 1460;
            }
            else 
            {
            		System.out.println("Error in salt input");
            }
   
            		//sets the number of events (discharge, nothing, charge)
	            int N_EVENTS = 8;  // remember to change this; also remember java indexes arrays starting at 0
	            int[] event = new int[N_EVENTS];
	            double[] endt_event = new double[N_EVENTS];
	            int[] tsteps_per_event = new int[N_EVENTS];
	            
	            //the event array is initiated with 0, 1 or 2 to correspond to a charge, constant, or discharge cycle
	            event[0]=0;  // discharge
	            event[1]=1;  // do nothing for a while
	            event[2]=2;  // charge
	            event[3]=1;
	            event[4]=0;
	            event[5]=2;
	            event[6]=1;
	            event[7]=0;
	            
	            // can string together as many events 0, 1 and 2 as you want, but code will give garbage results if either tank becomes totally empty
	            // or if heat loss to environment is too high
	            // this can be fixed with some more coding but for now just keep an eye on the plots
	            
	            //sets the end time of each event
	            endt_event[0]=1800.;  // event 1 ends after 1800 s (30 minutes)
	            endt_event[1]=3600.;  // event 2 ends at 3600 s (hence lasts 1800 s)
	            endt_event[2]=5400.;
	            endt_event[3]=12000.;
	            endt_event[4]=15000.;
	            endt_event[5]=17600.;
	            endt_event[6]=36000.;
	            endt_event[7]=38400.;
	
	            double dt = 100.;  // time step in seconds - I haven't tested convergence, could possibly jack this way up
	
	      // end of input section 
	            
	            int i,j,k;
	            //calculates tank volume and area of the sides 
	            final double TANK_TOTVOL = Math.PI*TANK_RADIUS*TANK_RADIUS*TANK_HEIGHT;
	            final double TANK_SIDESA = 2*Math.PI*TANK_RADIUS*TANK_HEIGHT;
	            
	            int nt = 0;
	            
	            //calculates number of time steps in the fist event 
	            tsteps_per_event[0]=(int)Math.round(endt_event[0]/dt);
	            endt_event[0]=dt*(double)tsteps_per_event[0];
	            
	            //calculates total number of time steps for all events
	            nt = tsteps_per_event[0];
		            for(i=1; i<N_EVENTS; i++) 
		            {
		                tsteps_per_event[i]=(int)Math.round((endt_event[i]-endt_event[i-1])/dt);
		                endt_event[i]=endt_event[i-1]+dt*(double)tsteps_per_event[i];
		                nt=nt+tsteps_per_event[i];
		            }//for int i 
	            
		           
	            double[] t = new double[nt];
	            t[0]=0.;
		            for (i=1; i<nt; i++) 
		            {
		                t[i]=t[i-1]+dt;
		            }//for int i
		        
		        //creates arrays of length nt for temperatures, volumes, and masses
		        double[] Th = new double[nt];
	            double[] Tc = new double[nt];
	            double[] Vh = new double[nt];
	            double[] Vc = new double[nt];
	            double[] mh = new double[nt];
	            double[] mc = new double[nt];
	            //double[] mdotchg = new double[nt];
	            //double[] mdotdis = new double[nt];
		           
	            //creates temporary 2D arrays for Th, Tc, Vh, and Vc
	            //these will be used to store the results for Th, Tc, Vh, and Vc for all nine sensitivity trials  
		        double [][] Th_temp = new double[9][nt];
		        double [][] Tc_temp = new double[9][nt]; 
		        double [][] Vh_temp = new double[9][nt];
		        double [][] Vc_temp = new double[9][nt]; 
	       
		        int count = 0;
		        
		        //these two loops will iterate through the different sensitivity inputs
		        //h iterates through the heat transfer coefficients
		        //g iterates though the temperature ranges 
		    for(int h=0; h<3; h++)
	        {
		        for(int g=0; g<3; g++) 
		        {   
		        		//inital assignments/calculations for Th, Tc, Vc, and Vh at time step zero
		            Th[0]=T_HL[g];
		            Tc[0]=T_CL[g];
		            mh[0]=TANK_TOTVOL*rho*STARTING_HOT_TANK_VOL_FRACT;
		            mc[0]=TANK_TOTVOL*rho*(1.-STARTING_HOT_TANK_VOL_FRACT);
		            Vh[0]=mh[0]/rho;
		            Vc[0]=mc[0]/rho;
		           //mdotchg[0]=0.;
		           //mdotdis[0]=0.;
		            int ix = 1;
		            double mt;
		                    
		            //calculated temperatures and volumes at the time steps for each case
		            //the formulas used for Th, Tc, Vh, and Vc can be found in the report -- derived by Schneider 
		            for(i=0; i<N_EVENTS; i++)  
		            {
			            switch (event[i]) 
			            {
			            case 0:  // discharge
			                     for(j=0; j<tsteps_per_event[i]; j++) 
			                     {
			                         Th[ix]=Th[ix-1]-dt*(hx_cf[h]*(Th[ix-1]-T_inf)*TANK_SIDESA/(c_p*rho*Vh[ix-1]));
			                         Tc[ix]=Tc[ix-1]-dt*(hx_cf[h]*(Tc[ix-1]-T_inf)*TANK_SIDESA/(c_p*rho*Vc[ix-1]));
			                         Vh[ix]=Vh[ix-1]-dt*(P_OUT/(rho*c_p*(Th[ix-1]-T_CL[g]))); 
			                         Vc[ix]=Vc[ix-1]+dt*(P_OUT/(rho*c_p*(Th[ix-1]-T_CL[g])));
			                         mt=dt*(P_OUT/(rho*c_p*(Th[ix-1]-T_CL[g])));
			                         ix++;
			                     }//for int j
			                     ix--;
			                break;
			
			            case 1:  // just heat loss
			                    for(j=0; j<tsteps_per_event[i]; j++) 
			                    {
			                         Th[ix]=Th[ix-1]-dt*(hx_cf[h]*(Th[ix-1]-T_inf)*TANK_SIDESA/(c_p*rho*Vh[ix-1]));
			                         Tc[ix]=Tc[ix-1]-dt*(hx_cf[h]*(Tc[ix-1]-T_inf)*TANK_SIDESA/(c_p*rho*Vc[ix-1]));
			                         Vh[ix]=Vh[ix-1];
			                         Vc[ix]=Vc[ix-1];
			                         ix++;
			                     }//for int j
			                    ix--;
			                     break;
			            case 2:  // charge
			                     for(j=0; j<tsteps_per_event[i]; j++) 
			                     {
			                         Th[ix]=Th[ix-1]+dt*(P_IN*(T_HL[g]-Th[ix-1])/(Vh[ix-1]*c_p*rho*(T_HL[g]-Tc[ix-1]))-hx_cf[h]*(Th[ix-1]-T_inf)*TANK_SIDESA/(c_p*rho*Vh[ix-1]));//**
			                         Tc[ix]=Tc[ix-1]-dt*(P_IN*(T_HL[g]-Th[ix-1])/(Vh[ix-1]*c_p*rho*(T_HL[g]-Tc[ix-1]))+hx_cf[h]*(Tc[ix-1]-T_inf)*TANK_SIDESA/(c_p*rho*Vc[ix-1]));//**
			                         Vh[ix]=Vh[ix-1]+dt*(P_IN/(rho*c_p*(T_HL[g]-Tc[ix-1])));
			                         Vc[ix]=Vc[ix-1]-dt*(P_IN/(rho*c_p*(T_HL[g]-Tc[ix-1]))); 
			                         ix++;
			                     }//for int j
			                     break;
			            default: // shouldn't be here
			                     break;
			            }//switch 
		            }//for int i
		            
		            		//once all temperatures and volumes are calculated for a given sensitivity trial (1-9), the results are stored
		            		//in the 2D temporary arrays
		            		for(int r=0; r<nt; r++)
			            {
			            	Th_temp[count][r] = Th[r];
			            	Tc_temp[count][r] = Tc[r];
			            	Vh_temp[count][r] = Vh[r];
			            	Vc_temp[count][r] = Vc[r];
			            }//for int r
		            		
		            		count ++;
		            		
	        		}//for int g
            	}//for int h
	        
		    //these 1D arrays are created to hold the temperature and volume results of each sensitivity trial
		    //9 possible combos for Th, Tc, Vh, and Vc
	        double [] Th1 = new double[nt];
	        double [] Tc1 = new double[nt];
	        double [] Vh1 = new double[nt]; 
	        double [] Vc1 = new double[nt]; 
	        
	        double [] Th2 = new double[nt];
	        double [] Tc2 = new double[nt]; 
	        double [] Vh2 = new double[nt];
	        double [] Vc2 = new double[nt];
	        
	        double [] Th3 = new double[nt];
	        double [] Tc3 = new double[nt];
	        double [] Vh3 = new double[nt]; 
	        double [] Vc3 = new double[nt]; 
	        
	        double [] Th4 = new double[nt];
	        double [] Tc4 = new double[nt];
	        double [] Vh4 = new double[nt]; 
	        double [] Vc4 = new double[nt]; 
	        
	        double [] Th5 = new double[nt];
	        double [] Tc5 = new double[nt]; 
	        double [] Vh5 = new double[nt];
	        double [] Vc5 = new double[nt];
	        
	        double [] Th6 = new double[nt];
	        double [] Tc6 = new double[nt];
	        double [] Vh6 = new double[nt]; 
	        double [] Vc6 = new double[nt]; 
	        
	        double [] Th7 = new double[nt];
	        double [] Tc7 = new double[nt];
	        double [] Vh7 = new double[nt]; 
	        double [] Vc7 = new double[nt]; 
	        
	        double [] Th8 = new double[nt];
	        double [] Tc8 = new double[nt]; 
	        double [] Vh8 = new double[nt];
	        double [] Vc8 = new double[nt];
	        
	        double [] Th9 = new double[nt];
	        double [] Tc9 = new double[nt];
	        double [] Vh9 = new double[nt]; 
	        double [] Vc9 = new double[nt]; 
	        
	        //assignments from the 2D temporary array are made to the individual 1D arrays defined above 
	        for(int y=0; y<nt; y++)
	        {
	        		Th1[y]=Th_temp[0][y];
	        		Tc1[y]=Tc_temp[0][y];
	        		Vh1[y]=Vh_temp[0][y];
	        		Vc1[y]=Vc_temp[0][y];
	        		
	        		Th2[y]=Th_temp[1][y];
	        		Tc2[y]=Tc_temp[1][y];
	        		Vh2[y]=Vh_temp[1][y];
	        		Vc2[y]=Vc_temp[1][y];
	        		
	        		Th3[y]=Th_temp[2][y];
	        		Tc3[y]=Tc_temp[2][y];
	        		Vh3[y]=Vh_temp[2][y];
	        		Vc3[y]=Vc_temp[2][y];
	        		
	        		Th4[y]=Th_temp[3][y];
	        		Tc4[y]=Tc_temp[3][y];
	        		Vh4[y]=Vh_temp[3][y];
	        		Vc4[y]=Vc_temp[3][y];
	        		
	        		Th5[y]=Th_temp[4][y];
	        		Tc5[y]=Tc_temp[4][y];
	        		Vh5[y]=Vh_temp[4][y];
	        		Vc5[y]=Vc_temp[4][y];
	        		
	        		Th6[y]=Th_temp[5][y];
	        		Tc6[y]=Tc_temp[5][y];
	        		Vh6[y]=Vh_temp[5][y];
	        		Vc6[y]=Vc_temp[5][y];
	        		
	        		Th7[y]=Th_temp[6][y];
	        		Tc7[y]=Tc_temp[6][y];
	        		Vh7[y]=Vh_temp[6][y];
	        		Vc7[y]=Vc_temp[6][y];
	        		
	        		Th8[y]=Th_temp[7][y];
	        		Tc8[y]=Tc_temp[7][y];
	        		Vh8[y]=Vh_temp[7][y];
	        		Vc8[y]=Vc_temp[7][y];
	        		
	        		Th9[y]=Th_temp[8][y];
	        		Tc9[y]=Tc_temp[8][y];
	        		Vh9[y]=Vh_temp[8][y];
	        		Vc9[y]=Vc_temp[8][y];
	        }// for y
	        
	       //y axis bounds are defined for the temperature and volume plots 
	       //***these values will change for each salt type
	       //kind of a guess and check situation (if nothing is showing up on the plot, then change the range)
	       //1  --> Th
	       //2  --> Tc
	       //3  --> Vh
	       //4  --> Vc
	       double [] yaxis1 = new double [] {400, 600};
	       double [] yaxis2 = new double [] {300, 400};
	       double [] yaxis3 = new double [] {100, 700};
	       double [] yaxis4 = new double [] {0, 600};
	             
	       //1D arrays are passed to the XYPlot file
	       //axis titles and yaxis bounds are also passed
            XYPlot.MakePlot(primaryStage,t,Th1,t,Th2,t,Th3,t,Th4,t,Th5,t,Th6,t,Th7,t,Th8,t,Th9,"time[s]","Hot tank temperature [K]", "Hot Tank Temperature for "+ salt, yaxis1); 
            primaryStage.show();
            
            XYPlot.MakePlot(fig2,t,Tc1,t,Tc2,t,Tc3,t,Tc4,t,Tc5,t,Tc6,t,Tc7,t,Tc8,t,Tc9,"time[s]","Cold tank temperature [K]", "Cold Tank Temperature for " + salt, yaxis2); 
            fig2.show();
            
            XYPlot.MakePlot(fig3,t,Vh1,t,Vh2,t,Vh3, t,Vh4, t,Vh5, t,Vh6,t,Vh7,t,Vh8,t,Vh9,"time[s]","Hot tank fluid volume [m3]","Hot Tank Volume for " + salt, yaxis3); 
            fig3.show();
            
            XYPlot.MakePlot(fig4,t,Vc1,t,Vc2,t,Vc3,t,Vc4,t,Vc5,t,Vc6,t,Vc7,t,Vc8,t,Vc9,"time[s]","Cold tank fluid volume [m3]","Cold Tank Volume for " + salt, yaxis4); 
            fig4.show();
           
           /*
            //this section creates a csv file of the data for the high and low sensitivity trials for Th, Tc, Vh, and Vc
            //this section should be commented out unless csv data is wanted
         	FileWriter writer = new FileWriter("/Users/kaylakelley/Desktop/JavaTempVoulumes_Glauber.csv");
	        for (int g = 0; g < nt; g++)
	        {
	        	    writer.append(String.valueOf(Tc3[g]));
	      		writer.append(",");
	      		writer.append(String.valueOf(Tc7[g]));
	      		writer.append(",");
	      		writer.append(",");
	      		writer.append(String.valueOf(Th1[g]));
		        writer.append(",");
		        writer.append(String.valueOf(Th9[g]));
		        writer.append(",");
	      		writer.append(",");
	      		writer.append(String.valueOf(Vc3[g]));
		        writer.append(",");
		        writer.append(String.valueOf(Vc7[g]));
		        writer.append(",");
	      		writer.append(",");
	      		writer.append(String.valueOf(Vh3[g]));
		        writer.append(",");
		        writer.append(String.valueOf(Vh7[g]));
	      		writer.append("\n");
	        } 		       
            writer.close();
             */
          
    }//public void start

    
    /**
     * @param args the command line arguments
     */
	    public static void main(String[] args) 
	    {
	        launch(args);
	    }
	    
}// public class main
